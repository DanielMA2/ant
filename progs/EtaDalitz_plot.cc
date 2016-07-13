/**
  * @file EtaDalitz_plot.cc
  * @brief Plotting tool for TTrees from the ant::analysis::physics::calibration::PID_Energy_etaDalitz physics class
  */

#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/physics/calibration/PID_Energy_etaDalitz.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"
#include "base/vec/vec2.h"
#include "base/ProgressCounter.h"

#include "expconfig/ExpConfig.h"

#include "tree/TAntHeader.h"

#include <list>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile sig_atomic_t interrupt = false;
static double binScale = 1.;


template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    const decltype(physics::PID_Energy_etaDalitz::makeChannels()) channels;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
        channels(physics::PID_Energy_etaDalitz::makeChannels())
    {
        using histstyle::Mod_t;

        this->GetHist(0, "Data",   Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Signal", Mod_t::MakeLine(kRed, 2));
        // mctrue is never >= 3 (and < 9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(2, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(3, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = unsigned(f.Tree.channel);

        using histstyle::Mod_t;

        auto get_bkg_name = [] (const unsigned mctrue) {
            const auto entry = physics::PID_Energy_etaDalitz::reaction_channels.channels.find(int(mctrue));

            if (entry != physics::PID_Energy_etaDalitz::reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        using histstyle::Mod_t;
        const Hist_t& hist = mctrue < 10 ? this->GetHist(mctrue)
                                         : this->GetHist(mctrue,
                                                         get_bkg_name(mctrue),
                                                         Mod_t::MakeLine(histstyle::color_t::Get(mctrue-10), 1, kGray+1)
                                                         );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if (mctrue > 0) {
            this->GetHist(2).Fill(f);
            if (mctrue >= 10)
                this->GetHist(3).Fill(f);
        }
    }
};

bool Contains(const interval<double>& i, const std::vector<double>& d)
{
    for (const auto& v : d)
        if (i.Contains(v))
            return true;

    return false;
}

double max(const std::vector<double>& data)
{
    return *max_element(data.cbegin(), data.cend());
}


// define the structs containing the histograms and the cuts
struct Hist_t {

    using Tree_t = physics::PID_Energy_etaDalitz::Tree_t;

    struct Fill_t {
        const Tree_t& Tree;

        Fill_t(const Tree_t& t) : Tree(t) {}

        double TaggW() const {
            return Tree.TaggW;
        }
    };

    template <typename Hist>
    using fillfunc_t = std::function<void(Hist*, const Fill_t&)>;

    template <typename Hist>
    struct HistFiller_t {
        fillfunc_t<Hist> func;
        Hist* h;
        HistFiller_t(Hist* hist, fillfunc_t<Hist> f): func(f), h(hist) {}
        void Fill(const Fill_t& data) const {
            func(h, data);
        }
    };

    template <typename Hist>
    struct HistMgr : std::list<HistFiller_t<Hist>> {

        using list<HistFiller_t<Hist>>::list;

        void Fill(const Fill_t& data) const {
            for (auto& h : *this)
                h.Fill(data);
        }
    };

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings pid_channels = BinSettings(
                ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID)->GetNChannels());

    const BinSettings Ebins    = Bins(1200, 0, 1200);

    const BinSettings Chi2Bins = BinSettings(250, 0, 25);
    const BinSettings probbins = BinSettings(250, 0,  1);

    const BinSettings IMbins   = Bins(1200,   0, 1200);
    const BinSettings MMbins   = Bins(1200, 400, 1600);

    const BinSettings pThetaBins = Bins( 125,  0,   50);
    const BinSettings pEbins     = Bins( 250,  0, 1000);
    const BinSettings PSAABins   = Bins(  60, 20,   60);
    const BinSettings PSARBins   = Bins( 100,  0,  450);
    const BinSettings TaggChBins = BinSettings(47);

    const BinSettings TaggTime   = BinSettings(240, -30, 30);
    const BinSettings CoplBins   = Bins(300, 0, 30);

    const BinSettings vetoEbins  = Bins(200, 0, 10);

    HistogramFactory HistFac;

    void AddTH1(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &bins, const string &name, fillfunc_t<TH1D> f)
    {
        h1.emplace_back(HistFiller_t<TH1D>(
                            HistFac.makeTH1D(title, xlabel, ylabel, bins, name), f));
    }

    void AddTH2(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &xbins, const BinSettings& ybins,
                const string &name, fillfunc_t<TH2D> f)
    {
        h2.emplace_back(HistFiller_t<TH2D>(
                            HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name), f));
    }

    Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf) {

        AddTH1("KinFitChi2", "#chi^{2}", "#", Chi2Bins, "KinFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.chi2, f.TaggW());
        });

        AddTH1("KinFitProb", "probability", "#", probbins, "KinFitProb",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.probability, f.TaggW());
        });

        AddTH1("3 photon IM", "3#gamma IM [MeV]", "#", IMbins, "etaIM",
               [] (TH1D* h, const Fill_t& f) {
//            TLorentzVector eta(0,0,0,0);
//            for(const auto& g : f.Tree.photons())
//                eta += g;
//            h->Fill(eta.M(), f.TaggW());
            h->Fill(f.Tree.eta().M(), f.TaggW());
        });

        AddTH1("3 photon IM fitted",  "3#gamma IM fit [MeV]", "#", IMbins, "etaIM_fitted",
               [] (TH1D* h, const Fill_t& f) {
//            TLorentzVector eta(0,0,0,0);
//            for(const auto& g : f.Tree.photons_fitted())
//                eta += g;
//            h->Fill(eta.M(), f.TaggW());
            h->Fill(f.Tree.eta_fit().M(), f.TaggW());
        });

        AddTH1("Missing Mass", "MM [MeV]", "", MMbins, "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });

        AddTH2("PID 2 charged 1 neutral", "PID Energy [MeV]", "PID Channel", vetoEbins, pid_channels, "eegPID",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_vetoE().at(i), f.Tree.photons_vetoChannel().at(i), f.TaggW());
        });

        AddTH2("Cluster Size vs. Energy", "Energy [MeV]", "Cluster Size", Ebins, BinSettings(50), "clusterSize_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_clusterSize().at(i), f.TaggW());
        });

        AddTH1("Tagger Time - CB Average Time", "t [ns]", "#", TaggTime, "TaggTime",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
        });

        AddTH2("dEvEproton", "E [MeV]", "dE [MeV]", Ebins, vetoEbins, "dEvE",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.p_fitted().Energy() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p_vetoE);
        });

        AddTH1("nCands", "# Candidates", "#", BinSettings(4, 3, 7), "nCands",
                [] (TH1D* h, const Fill_t& f) {

            h->Fill(f.Tree.nCands, f.TaggW());
        });

    }

    void Fill(const Fill_t& f) const
    {
        h1.Fill(f);
        h2.Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        vector<TH1*> v;
        v.reserve(h1.size()+h2.size());
        for (auto& e : h1)
            v.emplace_back(e.h);
        for (auto& e: h2)
            v.emplace_back(e.h);
        return v;
    }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts()
    {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                                 //{"Prob>0.02+mm", [] (const Fill_t& f) { return f.Tree.probability > 0.02 && f.Tree.mm().M()<1100 && f.Tree.mm().M() > 780; } }
                              {"Prob>0.02", [] (const Fill_t& f) { return f.Tree.probability > .02; }},
                              {"Prob>0.05", [] (const Fill_t& f) { return f.Tree.probability > .05; }}
                             });

        auto antiPi0Cut = [] (const Fill_t& f, const double low = 102., const double high = 170.) {
            const interval<double> pion_cut(low, high);
            TLorentzVector pi0;
            const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};

            for (const auto pi0_comb : pi0_combs) {
                pi0 = TLorentzVector(0., 0., 0., 0.);

                for (const auto idx : pi0_comb)
                    pi0 += TParticle(ParticleTypeDatabase::Photon, f.Tree.photons().at(idx));

                // check anti pi^0 cut
                if (pion_cut.Contains(pi0.M()))
                    return false;
            }

            return true;
        };

        auto IM2d_lin_cut = [] (const Fill_t& f) {
            const auto photons = f.Tree.photons();
            const auto vetos = f.Tree.photons_vetoE();

            const unsigned size = photons.size();
            unsigned idx[size];
            for (unsigned i = 0; i < size; i++)
                idx[i] = i;

            // sort idx according to the photons' veto energies
            sort(idx, idx+size, [vetos] (unsigned i, unsigned j) { return vetos[i] > vetos[j]; });

            // IM of the two clusters with the highest veto energies (e+ and e-)
            const double eeIM = (photons.at(idx[0]) + photons.at(idx[1])).M();

            return eeIM < (1.15*f.Tree.eta().M() - 170);
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"anti pi0", antiPi0Cut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"lin cut",  IM2d_lin_cut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"MM < 1030",  [] (const Fill_t& f) { return f.Tree.mm().M() < 1030; }}
                          });

        auto cleanEvent = [] (const Fill_t& f) {
            return f.Tree.nCands == 4;
        };

        auto dontcareclean = [] (const Fill_t&) {
            return true;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"clean",         cleanEvent},
                              {"dontcare",   dontcareclean}
                          });

        return cuts;
    }

};


int main(int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i", "input", "Input file", true, "", "input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch", "Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m", "maxevents", "Process only max events", false, "maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("", "tree", "Tree name", false, "T1", "treename");
    auto cmd_pres = cmd.add<TCLAP::SwitchArg>("p", "", "Presentation Mode", false);
    auto cmd_binscale = cmd.add<TCLAP::ValueArg<double>>("B", "bin-scale", "Bin Scaling", false, 1.0, "bins");

    cmd.parse(argc, argv);


    // open the TRint app as early as possible to prevent ROOT to create a new one automatically
    // which will cause problems because of bad ROOT internal data / pointer handling, might cause segfaults
    argc = 0;  // prevent TRint to parse any cmdline
    TRint app("EtaDalitz_plot", &argc, argv, nullptr, 0, true);


    // set signal handler after starting TRint, otherwise it will be overwritten with ROOT handlers
    signal(SIGINT, [] (int) {
        LOG(WARNING) << "Processing interrupted";
        interrupt = true;
    });


    // general styling settings for transparence
    gStyle->SetFillColor(0);
    gStyle->SetFillStyle(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);

    if (cmd_pres->isSet()) {
        gStyle->SetLabelSize(.05f, "XYZ");
        gStyle->SetTextSize(.05f);
        gStyle->SetCanvasBorderSize(0);
    }

    if (cmd_binscale->isSet()) {
        binScale = cmd_binscale->getValue();
    }


    WrapTFileInput input;
    string setup_name;
    try {

        input.OpenFile(cmd_input->getValue());

        auto header = input.GetSharedClone<TAntHeader>("AntHeader");

        if (!header) {
            LOG(WARNING) << "No TAntHeader found in " << cmd_input->getValue();
            return 1;
        }

        setup_name = header->SetupName;

    } catch (const std::runtime_error& e) {
        LOG(ERROR) << "Can't open " << cmd_input->getValue() << " " << e.what();
    }

    if (setup_name.empty())
        setup_name = "Setup_2014_07_EPT_Prod";

    ExpConfig::Setup::Get(setup_name);


    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if (!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree " + treename + " in input file");
        if (expected_entries >= 0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree " + treename + " does not have entries == " + to_string(expected_entries));
        if (wraptree->Matches(t, false)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };


    Hist_t::Tree_t tree;

    if (!link_branches("PID_Energy_etaDalitz/tree", addressof(tree), -1)) {
        LOG(ERROR) << "Cannot link branches of tree";
        //return 1;
    }

    const auto entries = tree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if (cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           WrapTFileOutput::mode_t::recreate,
                                                           true); // cd into masterFile upon creation
    }


    HistogramFactory HistFac("EtaDalitz");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");

    auto signal_hists = cuttree::Make<MCTrue_Splitter<Hist_t>>(HistFac,
                                              sanitized_treename,
                                              Hist_t::GetCuts()
                                              );

    LOG(INFO) << "Tree entries = " << entries;

    auto max_entries = entries;
    if (cmd_maxevents->isSet() && cmd_maxevents->getValue().back() < entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    long long entry = 0;
    double last_percent = 0;
    ProgressCounter::Interval = 3;
    ProgressCounter progress(
                [&entry, entries, &last_percent] (std::chrono::duration<double> elapsed) {
        const double percent = 100.*entry/entries;
        const double speed = (percent - last_percent)/elapsed.count();
        LOG(INFO) << setw(2) << setprecision(4) << "Processed " << percent << " %, ETA: " << ProgressCounter::TimeToStr((100-percent)/speed);
        last_percent = percent;
    });

    while (entry < max_entries) {
        if (interrupt)
            break;

        tree.Tree->GetEntry(entry);
        cuttree::Fill<MCTrue_Splitter<Hist_t>>(signal_hists, {tree});

        entry++;

        ProgressCounter::Tick();
    }

    LOG(INFO) << "Analyzed " << entry << " events, speed "
              << entry/progress.GetTotalSecs() << " event/s";

    if (!cmd_batchmode->isSet()) {
        if (!std_ext::system::isInteractive())
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        else {
            if (masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if (masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
