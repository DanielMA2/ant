#include "EtapDalitzMC.h"

#include "utils/Combinatorics.h"
#include "utils/Matcher.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "base/Logger.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::physics;


EtapDalitzMC::PerChannel_t::PerChannel_t(const std::string& Name, const string& Title, HistogramFactory& hf):
    title(Title),
    name(Name)
{
    // only filled if not just the reference channel is analysed
    if (Settings_t::get().reference_only())
        return;

    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);
    const BinSettings pid_channels(detector->GetNChannels());
    const BinSettings veto_energy(1000, 0, 10);
    const BinSettings energy(1200);

    steps = hf.makeTH1D(title + " Accepted Events", "step", "#", BinSettings(10), name + " steps");

    etapIM = hf.makeTH1D(title + " IM #eta' all comb", "IM [MeV]", "#", energy, name + " etapIM");
    etapIM_kinfit = hf.makeTH1D(title + " IM #eta' kinfitted", "IM [MeV]", "#", energy, name + " etapIM_kinfit");
    MM = hf.makeTH1D(title + " Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), name + " MM");
    hCopl = hf.makeTH1D(title + " Coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + " hCopl");
    trueZVertex = hf.makeTH1D(title + " true Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " trueZ");
    true_rec_ZVertex = hf.makeTH2D(title + " true vs. reconstructed Z Vertex", "rec. vz [cm]", "true vz [cm]",
                                   BinSettings(100, -10, 10), BinSettings(100, -10, 10), name + " true_vs_rec_vz");

    treefitChi2 = hf.makeTH1D(title + " treefitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " treefitChi2");
    treefitProb = hf.makeTH1D(title + " treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefitProb");
    treefitIter = hf.makeTH1D(title + " treefitted # Iterations", "#iterations", "#", BinSettings(20), name + " treefitIter");
    treefit_freeZ_chi2 = hf.makeTH1D(title + " free Z treefitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " treefit_freeZ_chi2");
    treefit_freeZ_prob = hf.makeTH1D(title + " free Z treefitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " treefit_freeZ_prob");
    treefit_freeZ_iter = hf.makeTH1D(title + " free Z treefitted # Iterations", "#iterations", "#", BinSettings(20), name + " treefit_freeZ_iter");
    kinfitChi2 = hf.makeTH1D(title + " kinfitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " kinfitChi2");
    kinfitProb = hf.makeTH1D(title + " kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfitProb");
    kinfitIter = hf.makeTH1D(title + " kinfitted # Iterations", "#iterations", "#", BinSettings(20), name + " kinfitIter");
    kinfit_freeZ_chi2 = hf.makeTH1D(title + " free Z kinfitted #chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), name + " kinfit_freeZ_chi2");
    kinfit_freeZ_prob = hf.makeTH1D(title + " free Z kinfitted Probability", "probability", "#", BinSettings(500, 0, 1), name + " kinfit_freeZ_prob");
    kinfit_freeZ_iter = hf.makeTH1D(title + " free Z kinfitted # Iterations", "#iterations", "#", BinSettings(20), name + " kinfit_freeZ_iter");
    kinfit_ZVertex = hf.makeTH1D(title + " kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_ZVertex");
    kinfit_freeZ_ZVertex = hf.makeTH1D(title + " free Z kinfitted Z Vertex", "z [cm]", "#", BinSettings(100, -10, 10), name + " kinfit_freeZ_ZVertex");
    treefit_ZVertex = hf.makeTH1D(title + " treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_ZVertex");
    treefit_freeZ_ZVertex = hf.makeTH1D(title + " free Z treefitted Z Vertex", "z [cm]", "#", BinSettings(300, -30, 30), name + " treefit_freeZ_ZVertex");
    antiPionProb = hf.makeTH1D(title + " Probability anti-#pi Fit", "probability", "#", BinSettings(500, 0, 1), name + " antiPionProb");

    if (Settings_t::get().less_plots())
        return;

    steps_vs_IMee = hf.makeTH2D("Accepted Events vs. Dilepton Mass", "IM(e^{+}e^{-}) [MeV]", "Steps",
                                BinSettings(20, 0, 1000), BinSettings(10), name + " steps_vs_IMee");

    etapIM_kinfit_freeZ = hf.makeTH1D(title + " IM #eta' free Z kinfitted", "IM [MeV]", "#", energy, name + " etapIM_kinfit_freeZ");
    etapIM_treefit = hf.makeTH1D(title + " IM #eta' treefitted", "IM [MeV]", "#", energy, name + " etapIM_treefit");
    etapIM_treefit_freeZ = hf.makeTH1D(title + " IM #eta' free Z treefitted", "IM [MeV]", "#", energy, name + " etapIM_treefit_freeZ");

    etapIM_cand = hf.makeTH1D(title + " IM #eta' candidates", "IM [MeV]", "#", energy, name + " etapIM_cand");
    etapIM_final = hf.makeTH1D(title + " IM #eta' final", "IM [MeV]", "#", energy, name + " etapIM_final");
    IM2d = hf.makeTH2D(title + " IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), name + " IM2d");
    hCopl_final = hf.makeTH1D(title + " Coplanarity #eta - proton final", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), name + " hCopl_final");

    effect_rad = hf.makeTH1D(title + " Effective Radius", "R", "#", BinSettings(500, 0, 50), name + " effect_rad");
    effect_rad_E = hf.makeTH2D(title + " Effective Radius vs. Cluster Energy", "E [MeV]", "R", energy, BinSettings(500, 0, 50), name + " effect_rad_E");
    cluster_size = hf.makeTH1D(title + " Cluster Size", "N", "#", BinSettings(50), name + " cluster_size");
    cluster_size_E = hf.makeTH2D(title + " Cluster Size vs. Cluster Energy", "E [MeV]", "N", energy, BinSettings(50), name + " cluster_size_E");
    lat_moment = hf.makeTH1D(title + " Lateral Moment", "L", "#", BinSettings(200, 0, 1), name + " lat_moment");
    lat_moment_E = hf.makeTH2D(title + " Lateral Moment vs. Cluster Energy", "E [MeV]", "L", energy, BinSettings(200, 0, 1), name + " lat_moment_E");

    proton_E_theta = hf.makeTH2D(title + " proton", "E [MeV]", "#vartheta [#circ]", energy, BinSettings(360, 0, 180), name + " e_theta");
}

void EtapDalitzMC::PerChannel_t::Show()
{
    //canvas("Per Channel: " + title) << drawoption("colz") << IM2d << endc;
    canvas("Per Channel: " + title) << steps
                                    << etapIM_kinfit
                                    //<< etapIM_treefit
                                    //<< etapIM_final
                                    << kinfitChi2
                                    << kinfitProb
                                    << kinfit_freeZ_chi2
                                    << kinfit_freeZ_prob
                                    << treefitChi2
                                    << treefitProb
                                    << endc;
}

void EtapDalitzMC::PerChannel_t::Fill(const TEventData& d)
{
    if (Settings_t::get().less_plots())
        return;

    auto particles = d.ParticleTree ?
                         utils::ParticleTypeList::Make(d.ParticleTree) :
                         utils::ParticleTypeList::Make(d.Candidates);
    const auto& protons = particles.Get(ParticleTypeDatabase::Proton);
    if (!protons.empty()) {
        const auto& p = protons.at(0);
        proton_E_theta->Fill(p->Ek(), p->Theta()*TMath::RadToDeg());
    }
}

mev_t EtapDalitzMC::calcEnergySum(const TParticleList& particles) const
{
    double esum = 0.;

    for (const TParticlePtr& p : particles)
        if (geo.DetectorFromAngles(p->Theta(), p->Phi()) == Detector_t::Type_t::CB)
            esum += p->Ek();

    for (const auto& p : getGeoAcceptedDetector(particles, Detector_t::Type_t::CB))
        esum += p->Ek();

    return esum;
}

TParticleList EtapDalitzMC::getGeoAccepted(const TParticleList &particles) const
{
    TParticleList list;

    for (auto& p : particles)
        if (geo.DetectorFromAngles(p->Theta(), p->Phi()) != Detector_t::Any_t::None)
            list.emplace_back(p);

    return list;
}

TParticleList EtapDalitzMC::getGeoAcceptedDetector(const TParticleList &particles,
                                                   const Detector_t::Type_t d) const
{
    TParticleList list;

    for (auto& p : particles)
        if (geo.DetectorFromAngles(p->Theta(), p->Phi()) == d)
            list.emplace_back(p);

    return list;
}

template <typename T, typename>
T EtapDalitzMC::geoAccepted(const TParticleList& particles) const
{
    auto n = count_if(particles.begin(), particles.end(), [this] (const TParticlePtr p) {
        return geo.DetectorFromAngles(p->Theta(), p->Phi()) != Detector_t::Any_t::None;
    });

    return static_cast<T>(n);
}

template <typename T, typename>
T EtapDalitzMC::geoAccepted(const TCandidateList& cands) const
{
    auto n = count_if(cands.begin(), cands.end(), [this] (const TCandidate& c) {
        return geo.DetectorFromAngles(c.Theta, c.Phi) != Detector_t::Any_t::None;
    });

    return static_cast<T>(n);
}

template <typename T, typename>
T EtapDalitzMC::geoAcceptedDetector(const TParticleList& particles, const Detector_t::Type_t d) const
{
    auto n = count_if(particles.begin(), particles.end(), [this, d] (const TParticlePtr p) {
        return geo.DetectorFromAngles(p->Theta(), p->Phi()) == d;
    });

    return static_cast<T>(n);
}

template <typename T, typename>
T EtapDalitzMC::geoAcceptedDetector(const TCandidateList& cands, const Detector_t::Type_t d) const
{
    auto n = count_if(cands.begin(), cands.end(), [this, d] (const TCandidate& c) {
        return geo.DetectorFromAngles(c.Theta, c.Phi) == d;
    });

    return static_cast<T>(n);
}

EtapDalitzMC::EtapDalitzMC(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    model_MC(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 make_shared<utils::UncertaintyModels::FitterSergey>()
                 )),
    kinfit(model_MC, opts->HasOption("SigmaZ"), EtapDalitz::MakeFitSettings(20)),
    kinfit_freeZ(model_MC, true,                EtapDalitz::MakeFitSettings(20)),
    treefitter_etap(etap_3g(), model_MC,
                    opts->HasOption("SigmaZ"), {}, EtapDalitz::MakeFitSettings(20)
                    ),
    treefitter_etap_freeZ(etap_3g(), model_MC,
                          true, {}, EtapDalitz::MakeFitSettings(20)
                          )
{
    promptrandom.AddPromptRange({-3, 2});
    promptrandom.AddRandomRange({-35, -10});
    promptrandom.AddRandomRange({10, 35});

    // initialize settings
    settings.init(opts->Get<bool>("reference", 0),
                  opts->Get<bool>("reference_only", 0),
                  opts->Get<bool>("less_plots", 0));

    cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    sig.CreateBranches(HistFac.makeTTree("signal"));
    if (settings.reference()) {
        if (settings.reference_only())
            LOG(INFO) << "Only Reference channel will be analysed";
        else
            LOG(INFO) << "Reference channel included in analysis";
        ref.CreateBranches(HistFac.makeTTree("ref"));
        etap2g = new Etap2gMC("Etap2gMC", opts);
        etap2g->linkTree(ref);
        etap2g->setPromptRandom(promptrandom);
    }

    // skip all the remaining instantiations of the ctor if just the reference channel should be analysed
    if (settings.reference_only())
        return;

    mc.CreateBranches(HistFac.makeTTree("MC"));

    if (settings.less_plots())
        LOG(INFO) << "Less histograms will be created and stored";

    h_counts = HistFac.makeTH1D("Events per Channel", "channel", "#", BinSettings(20), "h_counts");
    h_nCands = HistFac.makeTH1D("Number of Candidates", "#Candidates", "#", BinSettings(30), "h_nCands");
    missed_channels = HistFac.makeTH1D("Unlisted Channels", "", "Total Events seen", BinSettings(20), "missed_channels");
    found_channels  = HistFac.makeTH1D("Listed Channels", "", "Total Events seen", BinSettings(20), "found_channels");

    const auto IMee_bins = BinSettings(20, 0, 1000);
    const auto energy_2d = BinSettings(300, 600, 1200);
    const string IMee_label = "IM(e^{+}e^{-}) [MeV]";
    h_etapIM_vs_IMee = HistFac.makeTH2D("#eta' IM vs. Dilepton Mass", IMee_label, "IM [MeV]",
                                        IMee_bins, energy_2d, "h_etapIM_vs_IMee");
    h_MM_vs_IMee = HistFac.makeTH2D("MM vs. Dilepton Mass", IMee_label, "MM [MeV]",
                                    IMee_bins, energy_2d, "h_MM_vs_IMee");
    h_etapIM_fitted_vs_IMee = HistFac.makeTH2D("Fitted #eta' vs. Dilepton Mass", IMee_label, "IM [MeV",
                                               IMee_bins, energy_2d, "h_etapIM_fitted_vs_IMee");

    h_IMee_true = HistFac.makeTH1D("True IM(e+e-)", "IM(e+e-) [MeV]", "#", BinSettings(100, 0, 1000), "h_IMee_true");

    const BinSettings energybins(1000, 0, 10);

    if (!settings.less_plots()) {
        h_cluster_CB = HistFac.makeTH1D("#Cluster CB", "#Cluster", "#", BinSettings(20), "h_cluster_CB");
        h_cluster_TAPS = HistFac.makeTH1D("#Cluster TAPS", "#Cluster", "#", BinSettings(20), "h_cluster_TAPS");

        h_pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "h_pTheta");
        h_protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", energybins, "h_protonVeto");
        h_etapIM_final = HistFac.makeTH1D("IM #eta' final", "IM [MeV]", "#", BinSettings(1200), "h_etapIM_final");
        h_IM2d = HistFac.makeTH2D("IM(e+e-) vs IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(1200), BinSettings(1200), "h_IM2d");
        h_etap = HistFac.makeTH2D("Kinematics #eta'", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(360, 0, 180), "h_etap");
        h_proton = HistFac.makeTH2D("Kinematics p", "Energy [MeV]", "#vartheta [#circ]", BinSettings(1200), BinSettings(160, 0, 80), "h_proton");
        h_subIM_2g = HistFac.makeTH1D("#pi^{0} Candidate sub IM 2#gamma", "IM [MeV]", "#", BinSettings(1600, 0, 400), "h_subIM_2g");
        h_subIM_2g_fit = HistFac.makeTH1D("#pi^{0} Candidate sub IM 2#gamma after KinFit", "IM [MeV]", "#", BinSettings(1600, 0, 400), "h_subIM_2g_fit");

        h_energy_deviation = HistFac.makeTH1D("Energy Deviation #eta' Pluto - Geant", "#DeltaE [MeV]", "#", BinSettings(600, -300, 300), "h_MCdeltaE");
        h_fsClE_vs_pluto_geant_dE = HistFac.makeTH2D("#eta' FS Energy vs. Energy Difference #eta' Pluto - Geant", "#DeltaE [MeV]", "E [MeV]",
                                                     BinSettings(400, -200, 200), BinSettings(700), "h_fsClE_MCdE");
        h_theta_vs_vz = HistFac.makeTH2D("#vartheta #eta' vs. Kinfit v_{z}", "v_{z} [cm]", "#vartheta [#circ]",
                                         BinSettings(100, -10, 10), BinSettings(90), "h_vz_theta");
        h_theta_vs_pluto_geant_dE = HistFac.makeTH2D("#vartheta Dependence of Energy Difference #eta' Pluto - Geant", "dE [MeV]", "#vartheta [#circ]",
                                                     BinSettings(400, -200, 200), BinSettings(45, 0, 90), "h_theta_MCdE");
        h_vz_vs_pluto_geant_dE = HistFac.makeTH2D("Kinfit v_{z} Dependence of Energy Difference #eta' Pluto - Geant", "dE [MeV]", "v_{z} [cm]",
                                                  BinSettings(400, -200, 200), BinSettings(100, -10, 10), "h_vz_MCdE");
        h_delta_vz_vs_pluto_geant_dE = HistFac.makeTH2D("v_{z} Difference vs Energy Difference #eta' Pluto - Geant", "#DeltaE [MeV]", "#Deltav_{z} [cm]",
                                                        BinSettings(400, -200, 200), BinSettings(100, -5, 5), "h_dvz_MCdE");

        // histograms to check devitation between true and reconstructed MC events (resolution)
        // energy-related
        const auto e_sigma = BinSettings(600, -150, 150);
        const auto e_sigma_2d = BinSettings(150, -150, 150);
        const auto theta = BinSettings(160);
        const auto e_true = BinSettings(200, 0, 1000);
        h_energy_resolution_g = HistFac.makeTH1D("Energy Deviation #gamma Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_g");
        h_energy_resolution_em = HistFac.makeTH1D("Energy Deviation e^{-} Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_em");
        h_energy_resolution_ep = HistFac.makeTH1D("Energy Deviation e^{+} Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_ep");
        h_energy_resolution_g_fit = HistFac.makeTH1D("Energy Deviation fitted #gamma Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_g_fit");
        h_energy_resolution_em_fit = HistFac.makeTH1D("Energy Deviation fitted e^{-} Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_em_fit");
        h_energy_resolution_ep_fit = HistFac.makeTH1D("Energy Deviation fitted e^{+} Pluto - Geant", "#sigmaE [MeV]", "#", e_sigma, "h_MCsigmaE_ep_fit");
        h_energy_resolution_vs_theta_g = HistFac.makeTH2D("Energy Deviation vs. #vartheta #gamma Pluto - Geant", "#vartheta [#circ]", "#sigmaE [MeV]",
                                                          theta, e_sigma_2d, "h_MCsigmaE_theta_g");
        h_energy_resolution_vs_theta_em = HistFac.makeTH2D("Energy Deviation vs. #vartheta e^{-} Pluto - Geant", "#vartheta [#circ]", "#sigmaE [MeV]",
                                                           theta, e_sigma_2d, "h_MCsigmaE_theta_em");
        h_energy_resolution_vs_theta_ep = HistFac.makeTH2D("Energy Deviation vs. #vartheta e^{+} Pluto - Geant", "#vartheta [#circ]", "#sigmaE [MeV]",
                                                           theta, e_sigma_2d, "h_MCsigmaE_theta_ep");
        h_energy_resolution_vs_trueE_g = HistFac.makeTH2D("Energy Deviation vs. E_{true} #gamma Pluto - Geant", "E_{true} [MeV]", "#sigmaE [MeV]",
                                                          e_true, e_sigma_2d, "h_MCsigmaE_trueE_g");
        h_energy_resolution_vs_trueE_em = HistFac.makeTH2D("Energy Deviation vs. E_{true} e^{-} Pluto - Geant", "E_{true} [MeV]", "#sigmaE [MeV]",
                                                           e_true, e_sigma_2d, "h_MCsigmaE_trueE_em");
        h_energy_resolution_vs_trueE_ep = HistFac.makeTH2D("Energy Deviation vs. E_{true} e^{+} Pluto - Geant", "E_{true} [MeV]", "#sigmaE [MeV]",
                                                           e_true, e_sigma_2d, "h_MCsigmaE_trueE_ep");
        h_energy_resolution_vs_theta_g_fit = HistFac.makeTH2D("Energy Deviation vs. #vartheta fitted #gamma Pluto - Geant",
                                                              "#vartheta [#circ]", "#sigmaE [MeV]", theta, e_sigma_2d, "h_MCsigmaE_theta_g_fit");
        h_energy_resolution_vs_theta_em_fit = HistFac.makeTH2D("Energy Deviation vs. #vartheta fitted e^{-} Pluto - Geant",
                                                               "#vartheta [#circ]", "#sigmaE [MeV]", theta, e_sigma_2d, "h_MCsigmaE_theta_em_fit");
        h_energy_resolution_vs_theta_ep_fit = HistFac.makeTH2D("Energy Deviation vs. #vartheta fitted e^{+} Pluto - Geant",
                                                               "#vartheta [#circ]", "#sigmaE [MeV]", theta, e_sigma_2d, "h_MCsigmaE_theta_ep_fit");
        h_energy_resolution_vs_trueE_g_fit = HistFac.makeTH2D("Energy Deviation vs. E_{true} fitted #gamma Pluto - Geant",
                                                              "E_{true} [MeV]", "#sigmaE [MeV]", e_true, e_sigma_2d, "h_MCsigmaE_trueE_g_fit");
        h_energy_resolution_vs_trueE_em_fit = HistFac.makeTH2D("Energy Deviation vs. E_{true} fitted e^{-} Pluto - Geant",
                                                               "E_{true} [MeV]", "#sigmaE [MeV]", e_true, e_sigma_2d, "h_MCsigmaE_trueE_em_fit");
        h_energy_resolution_vs_trueE_ep_fit = HistFac.makeTH2D("Energy Deviation vs. E_{true} fitted e^{+} Pluto - Geant",
                                                               "E_{true} [MeV]", "#sigmaE [MeV]", e_true, e_sigma_2d, "h_MCsigmaE_trueE_ep_fit");
        // theta-related
        const auto theta_sigma = BinSettings(400, -50, 50);
        const auto theta_sigma_2d = BinSettings(100, -50, 50);
        const auto energy = BinSettings(200, 0, 1000);
        const auto theta_true = BinSettings(160);
        h_theta_resolution_g = HistFac.makeTH1D("Theta Deviation #gamma Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_g");
        h_theta_resolution_em = HistFac.makeTH1D("Theta Deviation e^{-} Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_em");
        h_theta_resolution_ep = HistFac.makeTH1D("Theta Deviation e^{+} Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_ep");
        h_theta_resolution_g_fit = HistFac.makeTH1D("Theta Deviation fitted #gamma Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_g_fit");
        h_theta_resolution_em_fit = HistFac.makeTH1D("Theta Deviation fitted e^{-} Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_em_fit");
        h_theta_resolution_ep_fit = HistFac.makeTH1D("Theta Deviation fitted e^{+} Pluto - Geant", "#sigma#vartheta [#circ]", "#", theta_sigma, "h_MCsigmaTheta_ep_fit");
        h_theta_resolution_vs_energy_g = HistFac.makeTH2D("Theta Deviation vs. Energy #gamma Pluto - Geant", "E [MeV]", "#sigma#vartheta [#circ]",
                                                          energy, theta_sigma_2d, "h_MCsigmaTheta_E_g");
        h_theta_resolution_vs_energy_em = HistFac.makeTH2D("Theta Deviation vs. Energy e^{-} Pluto - Geant", "E [MeV]", "#sigma#vartheta [#circ]",
                                                           energy, theta_sigma_2d, "h_MCsigmaTheta_E_em");
        h_theta_resolution_vs_energy_ep = HistFac.makeTH2D("Theta Deviation vs. Energy e^{+} Pluto - Geant", "E [MeV]", "#sigma#vartheta [#circ]",
                                                           energy, theta_sigma_2d, "h_MCsigmaTheta_E_ep");
        h_theta_resolution_vs_trueTheta_g = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} #gamma Pluto - Geant", "#vartheta_{true} [#circ]",
                                                             "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_g");
        h_theta_resolution_vs_trueTheta_em = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} e^{-} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                              "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_em");
        h_theta_resolution_vs_trueTheta_ep = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} e^{+} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                              "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_ep");
        h_theta_resolution_vs_energy_g_fit = HistFac.makeTH2D("Theta Deviation vs. Energy fitted #gamma Pluto - Geant",
                                                              "E [MeV]", "#sigma#vartheta [#circ]", energy, theta_sigma_2d, "h_MCsigmaTheta_E_g_fit");
        h_theta_resolution_vs_energy_em_fit = HistFac.makeTH2D("Theta Deviation vs. Energy fitted e^{-} Pluto - Geant",
                                                               "E [MeV]", "#sigma#vartheta [#circ]", energy, theta_sigma_2d, "h_MCsigmaTheta_E_em_fit");
        h_theta_resolution_vs_energy_ep_fit = HistFac.makeTH2D("Theta Deviation vs. Energy fitted e^{+} Pluto - Geant",
                                                               "E [MeV]", "#sigma#vartheta [#circ]", energy, theta_sigma_2d, "h_MCsigmaTheta_E_ep_fit");
        h_theta_resolution_vs_trueTheta_g_fit = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} fitted #gamma Pluto - Geant", "#vartheta_{true} [#circ]",
                                                                 "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_g_fit");
        h_theta_resolution_vs_trueTheta_em_fit = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} fitted e^{-} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                                  "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_em_fit");
        h_theta_resolution_vs_trueTheta_ep_fit = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} fitted e^{+} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                                  "#sigma#vartheta [#circ]", theta_true, theta_sigma_2d, "h_MCsigmaTheta_trueTheta_ep_fit");

        // histograms to check the dilepton mass dependence of certain kinematics
        const auto energy_bins = BinSettings(120, 0, 1200);
        h_IMee = HistFac.makeTH1D("Dilepton Mass", IMee_label, "#", BinSettings(1000), "h_IMee");
        h_IMee_fraction3CB1TAPS_total = HistFac.makeTH1D("Total Fraction of 3CB&1TAPS Cluster vs Dilepton Mass",
                                                         IMee_label, "Total Efficiency", IMee_bins,
                                                         "h_IMee_fraction3CB1TAPS_total", true);  // use Sumw2
        h_IMee_fraction3CB1TAPS_acceptance = HistFac.makeTH1D("Relative Fraction of 3CB&1TAPS Cluster to geo. Acceptance vs Dilepton Mass",
                                                              IMee_label, "Relative Efficiency", IMee_bins,
                                                              "h_IMee_fraction3CB1TAPS_acceptance", true);  // use Sumw2
        h_IMee_fraction3CB1TAPS_Trigger4Cl = HistFac.makeTH1D("Relative Fraction of 3CB&1TAPS Cluster vs Dilepton Mass",
                                                              IMee_label, "Relative Efficiency", IMee_bins,
                                                              "h_IMee_fraction3CB1TAPS_Trigger4Cl", true);  // use Sumw2
        h_nCands_vs_IMee = HistFac.makeTH2D("Number of Candidates vs Dilepton Mass", IMee_label, "#Candidates",
                                            IMee_bins, BinSettings(30), "h_nCands_vs_IMee");
        h_openingAngle_vs_IMee = HistFac.makeTH2D("Dilepton Opening Angle vs Dilepton Mass", IMee_label, "Opening Angle [#circ]",
                                                  BinSettings(100, 0, 1000), BinSettings(360, 0, 180), "h_openingAngle_vs_IMee");
        h_CBEsum = HistFac.makeTH1D("CB E_{sum}", "E_{sum} [MeV]", "#", BinSettings(1600), "h_CBEsum");
        h_CBEsum_true = HistFac.makeTH1D("CB E_{sum} true (based on CB geo. Acceptance)", "E_{sum} [MeV]", "#", BinSettings(1600), "h_CBEsum_true");
        h_CBEsum_vs_IMee = HistFac.makeTH2D("CB E_{sum} vs Dilepton Mass", IMee_label, "E_{sum} [MeV]",
                                            IMee_bins, BinSettings(800, 0, 1600), "h_CBEsum_vs_IMee");
        h_E_vs_IMee_eCharged_true = HistFac.makeTH2D("True e^{#pm} Energy vs Dilepton Mass", IMee_label, "E_{true} [MeV]",
                                                     IMee_bins, energy_bins, "h_E_vs_IMee_eCharged_true");
        h_E_vs_IMee_photon_true = HistFac.makeTH2D("True #gamma Energy vs Dilepton Mass", IMee_label, "E_{true} [MeV]",
                                                   IMee_bins, energy_bins, "h_E_vs_IMee_photon_true");
        h_E_vs_IMee_proton_true = HistFac.makeTH2D("True p Energy vs Dilepton Mass", IMee_label, "E_{true} [MeV]",
                                                   IMee_bins, BinSettings(70, 0, 700), "h_E_vs_IMee_proton_true");
        h_E_vs_IMee_eCharged_rec = HistFac.makeTH2D("Reconstructed e^{#pm} Energy vs Dilepton Mass",
                                                    IMee_label, "E_{rec} [MeV]", IMee_bins, energy_bins, "h_E_vs_IMee_eCharged_rec");
        h_E_vs_IMee_photon_rec = HistFac.makeTH2D("Reconstructed #gamma Energy vs Dilepton Mass",
                                                  IMee_label, "E_{rec} [MeV]", IMee_bins, energy_bins, "h_E_vs_IMee_photon_rec");
        h_E_vs_IMee_proton_rec = HistFac.makeTH2D("Reconstructed p Energy vs Dilepton Mass",
                                                  IMee_label, "E_{rec} [MeV]", IMee_bins, BinSettings(70, 0, 700), "h_E_vs_IMee_proton_rec");
        // IM(e+e-) count rate for total and relative efficiencies
        h_IMee_total = HistFac.makeTH1D("Dilepton Mass", IMee_label, "#", IMee_bins, "h_IMee_total");
        h_IMee_acceptance = HistFac.makeTH1D("Acceptance Dilepton Mass", IMee_label, "Acceptance", IMee_bins, "h_IMee_acceptance", true);
        h_IMee_Trigger4Cl = HistFac.makeTH1D("Dilepton Mass", IMee_label, "#", IMee_bins, "h_IMee_Trigger4Cl");

        // test histograms for some checks related to radiative corrections
        h_radCorr_x = HistFac.makeTH1D("Radiative Corrections: x", "x", "#", BinSettings(1000, 0, 1), "h_radCorr_x");
        h_radCorr_y = HistFac.makeTH1D("Radiative Corrections: y", "y", "#", BinSettings(1000, 0, 2), "h_radCorr_y");
        h_radCorr_checkBoundaries = HistFac.makeTH2D("Rad. Corrections: x and y divided by max boundary",
                                                     "x/x_{max}", "y/y_{max}", BinSettings(100, 0, 2), BinSettings(200, 0, 10), "h_radCorr_checkBoundaries");
        h_radCorr_y_vs_IMee = HistFac.makeTH2D("Radiative Corrections: y vs Dilepton Mass", IMee_label, "y",
                                               BinSettings(1000), BinSettings(20, 0, 1), "h_radCorr_y_vs_IMee");
        h_radCorr_y_vs_x = HistFac.makeTH2D("Radiative Corrections: y vs x", "x", "y", BinSettings(1000, 0, 1), BinSettings(20, 0, 1), "h_radCorr_y_vs_x");

        // test for cone prediction of proton candidate
        const auto theta_diff = BinSettings(160, -20, 20);
        h_angle_miss_res = HistFac.makeTH1D("Opening Angle Missing Momentum and TAPS Cluster", "#alpha [#circ]", "#", BinSettings(100, 0, 50), "h_angle_miss_res");
        h_theta_miss_res = HistFac.makeTH1D("Resolution #vartheta_{miss} in respect to TAPS Cluster", "#vartheta [#circ]", "#", theta_diff, "h_theta_miss_res");
        h_theta_vs_phi_miss_res = HistFac.makeTH2D("Resolution #vartheta_{miss} vs. #varphi_{miss} in respect to TAPS Cluster",
                                                   "#varphi [#circ]", "#vartheta [#circ]", BinSettings(360,-180,180), theta_diff, "h_theta_vs_phi_miss_res");
    }

    // get target information
    const auto target = ExpConfig::Setup::Get().GetTargetProperties();

    // set sigma to 0 for unmeasured --> free z vertex
    kinfit_freeZ.SetZVertexSigma(0);
    kinfit_freeZ.SetTarget(target.length);
    treefitter_etap_freeZ.SetZVertexSigma(0);
    treefitter_etap_freeZ.SetTarget(target.length);

    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        LOG(INFO) << "Fit Z vertex enabled with sigma = " << sigma_z;
        kinfit.SetZVertexSigma(sigma_z);
        treefitter_etap.SetZVertexSigma(sigma_z);
    }
}

void EtapDalitzMC::ProcessEvent(const TEvent& event, manager_t&)
{
    // only process MC
    if (!(sig.MCtrue = event.Reconstructed().ID.isSet(TID::Flags_t::MC)))
        return;

    const auto& data = event.Reconstructed();
    const auto& cands = data.Candidates;
    const auto& particletree = event.MCTrue().ParticleTree;

    if (!particletree) {
        LOG(ERROR) << "No particle tree found, only Geant or Pluto file provided, not both";
        return;
    }

    // check if the current event is the signal
    const bool signalMC = particletree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_eeg),
                                                utils::ParticleTools::MatchByParticleName);

    imee = std_ext::NaN;

    // signalMC histograms, matching; skip if reference_only is specified
    if (signalMC && !settings.reference_only()) {
        // first handle the leptons
        TParticleList mctrue(utils::ParticleTools::FindParticles(ParticleTypeDatabase::eCharged, particletree));
        assert(mctrue.size() == 2);

        imee = (*mctrue.front() + *mctrue.back()).M();
        mc.imee = imee;
        mc.opening = std_ext::radian_to_degree(TParticle::CalcAngle(mctrue.front(), mctrue.back()));
        h_IMee_true->Fill(imee);

        // now get the other particles
        mctrue.emplace_back(utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon, particletree));
        mctrue.emplace_back(utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree));

        // do some matching
        const auto matched = utils::match1to1(mctrue, cands.get_ptr_list(),
                                              [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
            return p1->Angle(*p2); },
                                              IntervalD(0., std_ext::degree_to_radian(15.)));

        // get the true particles
        const auto em = utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus, particletree);
        const auto ep = utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus, particletree);
        const auto g = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon, particletree);
        const auto p = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree);

        // store true information in the tree
        for (const auto& tp : {em, ep, g, p}) {
            mc.names().emplace_back(tp->Type().Name());
            mc.energies_true().push_back(tp->Ek());
            mc.thetas_true().push_back(std_ext::radian_to_degree(tp->Theta()));
            mc.phis_true().push_back(std_ext::radian_to_degree(tp->Phi()));
        }

        if (matched.size() == mctrue.size()) {
            const auto matched_em = utils::FindMatched(matched, em);
            const auto matched_ep = utils::FindMatched(matched, ep);
            const auto matched_g = utils::FindMatched(matched, g);
            const auto matched_p = utils::FindMatched(matched, p);

            // store matched information in the tree
            // order of the list used is the same, hence the information order should be the same
            for (const auto& p : {matched_em, matched_ep, matched_g, matched_p}) {
                mc.energies().push_back(p->CaloEnergy);
                mc.thetas().push_back(std_ext::radian_to_degree(p->Theta));
                mc.phis().push_back(std_ext::radian_to_degree(p->Phi));
            }

            if (!settings.less_plots()) {
                h_E_vs_IMee_eCharged_rec->Fill(mc.imee, matched_em->CaloEnergy);
                h_E_vs_IMee_eCharged_rec->Fill(mc.imee, matched_ep->CaloEnergy);
                h_E_vs_IMee_photon_rec->Fill(mc.imee, matched_g->CaloEnergy);
                h_E_vs_IMee_proton_rec->Fill(mc.imee, matched_p->CaloEnergy);
            }
        }

        const bool allAccepted = (geoAccepted({em, ep, g, p}) == mctrue.size());
        const bool etapFSinCB = (geoAcceptedDetector({em, ep, g}, Detector_t::Type_t::CB) == 3);
        const bool etapCBprotonTAPS = (etapFSinCB && geoAcceptedDetector({p}, Detector_t::Type_t::TAPS) == 1);

        if (!settings.less_plots()) {
            h_IMee->Fill(mc.imee);
            h_IMee_total->Fill(mc.imee);
            h_openingAngle_vs_IMee->Fill(mc.imee, mc.opening);

            if (allAccepted) {
                h_IMee_acceptance->Fill(mc.imee);
                if (etapCBprotonTAPS)
                    h_IMee_fraction3CB1TAPS_acceptance->Fill(mc.imee);
            }

            h_E_vs_IMee_eCharged_true->Fill(mc.imee, em->Ek());
            h_E_vs_IMee_eCharged_true->Fill(mc.imee, ep->Ek());
            h_E_vs_IMee_photon_true->Fill(mc.imee, g->Ek());
            h_E_vs_IMee_proton_true->Fill(mc.imee, p->Ek());

            // radiative corrections, check x and y variables
            const auto etap = utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree);
            double q2 = (*ep + *em).M2();
            double im2 = etap->M2();
            double x = q2/im2;
            //double y = meson->Vect4()*(lp-lm);
            double y = (*etap).Dot(*ep - *em);
            double nu2 = 4*(*em).M2()/im2;
            double beta = sqrt(1-nu2/x);
            //cout << "[DEBUG]   nu = " << sqrt(nu2) << endl;
            //cout << "[DEBUG]   beta = " << beta << endl;
            // use absolute value for y since correction values are just provided for positive y values
            // y should be symmetric so under this assumption everything should be fine
            y = 2*abs(y)/etap->M2()/(1-x);
            //cout << "[DEBUG]   x elem [" << nu2 << " , 1]" << endl;
            //cout << "[DEBUG]   y elem [0 , " << beta << "]" << endl;

            /*if ((x < nu2) || (x > 1.))
                cerr << "x value outside of kinematical bounds: x = " << x << " not in [" << nu2 << " , 1]" << endl;
            if ((y < 0.) || (y > beta))
                cerr << "y value outside of kinematical bounds: y = " << y << " not in [0 , " << beta << "]" << endl;*/
            h_radCorr_x->Fill(x);
            h_radCorr_y->Fill(y);
            h_radCorr_checkBoundaries->Fill(x, y/beta);
            h_radCorr_y_vs_IMee->Fill(mc.imee, y);
            h_radCorr_y_vs_x->Fill(x, y);
        }
    }

    // CB Esum
    if (!settings.less_plots() && !settings.reference_only()) {
        h_CBEsum_true->Fill(event.MCTrue().Trigger.CBEnergySum);

        auto CBsum = [] (double sum, const TCandidate& c) {
            if (c.Detector & Detector_t::Type_t::CB)
                sum += c.CaloEnergy;
            return sum;
        };

        double CBEsum = accumulate(cands.begin(), cands.end(), 0., CBsum);
        h_CBEsum->Fill(CBEsum);
        h_CBEsum_vs_IMee->Fill(mc.imee, CBEsum);
    }

    triggersimu.ProcessEvent(event);
    sig.init();

    sig.channel = reaction_channels.identify(event.MCTrue().ParticleTree);
    if (!sig.channel)  // assign other_index in case of an empty or unknown particle tree for MC (tagged as data otherwise)
        sig.channel = reaction_channels.other_index;
    sig.trueZVertex = event.MCTrue().Target.Vertex.z;  // NaN in case of data

    sig.nCands = cands.size();
    sig.CBSumE = triggersimu.GetCBEnergySum();
    sig.CBAvgTime = triggersimu.GetRefTiming();

    // up to this point no histograms are filled or tree entries are written, important if reference_only is set

    if (settings.reference()) {
        ref.init();
        ref.MCtrue = sig.MCtrue;
        ref.channel = sig.channel;
        ref.trueZVertex = sig.trueZVertex;
        ref.nCands = sig.nCands;
        ref.CBSumE = sig.CBSumE;
        ref.CBAvgTime = sig.CBAvgTime;

        etap2g->Process(event);

        if (settings.reference_only())
            return;
    }

    // starting the signal channel processing, this point is only reached if reference_only is not specified

    if (sig.channel == ReactionChannelList_t::other_index)
        missed_channels->Fill(utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree).c_str(), 1);
    else
        found_channels->Fill(sig.channel);

    // identify the currently processed channel
    channel_id(event, chan_id);

    // manage histogram structure for different channels, get histograms for current channel
    auto h = manage_channel_histograms_get_current(sig.MCtrue, event);
    h.trueZVertex->Fill(sig.trueZVertex);

    mc.multiplicity = cands.size();
    h_nCands->Fill(sig.nCands);
    if (!settings.less_plots())
        h_nCands_vs_IMee->Fill(mc.imee, sig.nCands);

    const auto fill_steps = [&] (const string& step) {
        h.steps->Fill(step.c_str(), 1);
        if (!settings.less_plots())
            h.steps_vs_IMee->Fill(imee, step.c_str(), 1);
    };

    fill_steps("seen");

    // histogram amount of CB and TAPS clusters
    size_t nCB = 0, nTAPS = 0;
    count_clusters(cands, nCB, nTAPS);
    mc.nCB = nCB;
    mc.nTAPS = nTAPS;
    if (!settings.less_plots()) {
        h_cluster_CB->Fill(nCB);
        h_cluster_TAPS->Fill(nTAPS);
        if (nCB == 3 && nTAPS == 1) {
            h_IMee_fraction3CB1TAPS_total->Fill(mc.imee);
            h_IMee_fraction3CB1TAPS_Trigger4Cl->Fill(mc.imee);
        }
        if (sig.nCands == 4 && triggersimu.HasTriggered())
            h_IMee_Trigger4Cl->Fill(mc.imee);
    }

    // up to this point are no selection or cut criteria applied
    // save the MC tree here to have all MC information available
    mc.fillAndReset();

    if (!triggersimu.HasTriggered())
        return;
    fill_steps("triggered");

    if (!isfinite(sig.CBAvgTime))
        return;
    fill_steps("CBAvgTime OK");


//    if (cands.size() != Cuts_t::N_FINAL_STATE)
//        return;
//    h.steps->Fill("#cands", 1);

    // q2 preselection on MC data
    if (Cuts_t::Q2_PRESELECTION) {
        if (!q2_preselection(event.MCTrue(), Cuts_t::Q2_MIN_VALUE))
            return;
        stringstream ss;
        ss << "MC q2 > " << Cuts_t::Q2_MIN_VALUE;
        fill_steps(ss.str());
    }

    //const auto mass_etap = ParticleTypeDatabase::EtaPrime.Mass();
    //const interval<double> etap_im({mass_etap-Cuts_t::ETAP_SIGMA, mass_etap+Cuts_t::ETAP_SIGMA});

    utils::ProtonPhotonCombs proton_photons(cands);


    double best_prob_fit = -std_ext::inf;
    // loop over all tagger hits
    for (const TTaggerHit& taggerhit : data.TaggerHits) {
        sig.reset();

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        fill_steps("time window");

        sig.TaggW = promptrandom.FillWeight();
        sig.TaggE = taggerhit.PhotonEnergy;
        sig.TaggT = taggerhit.Time;
        sig.TaggCh = taggerhit.Channel;


        // check the resolution between the missing momentum of the proton compared to the reconstructed TAPS cluster
        if (!settings.less_plots()
                && (nCB == 3 && nTAPS == 1)) {
            LorentzVec photon_sum;
            TParticlePtr proton;
            for (const auto& c : cands.get_iter())
                if (c->Detector & Detector_t::Type_t::CB)
                    photon_sum += TParticle(ParticleTypeDatabase::Photon, c);
                else if (c->Detector & Detector_t::Type_t::TAPS)
                    proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, c);

            const auto beam_target = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass());
            const auto miss_momentum = beam_target - photon_sum;

            h_angle_miss_res->Fill(std_ext::radian_to_degree(proton->Angle(miss_momentum)));
            h_theta_miss_res->Fill(std_ext::radian_to_degree(proton->Theta() - miss_momentum.Theta()));
            h_theta_vs_phi_miss_res->Fill(std_ext::radian_to_degree(proton->Phi() - miss_momentum.Phi()),
                                          std_ext::radian_to_degree(proton->Theta() - miss_momentum.Theta()));
        }

//        particle_combs_t selection = proton_photons()
//                .Observe([h] (const std::string& s) { h.steps->Fill(s.c_str(), 1.); }, "[S] ")
//                // require 3 photons and allow discarded energy of 100 MeV
//                .FilterMult(settings.n_final_state_etap, settings.max_discarded_energy)
//                .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(settings.mm_window_size).Round())  // MM cut on proton mass
//                .FilterCustom([=] (const particle_comb_t& p) {
//                    // ensure the possible proton candidate is kinematically allowed
//                    if (std_ext::radian_to_degree(p.Proton->Theta()) > settings.max_proton_theta)
//                        return true;
//                    return false;
//                }, "proton #vartheta")
//                .FilterCustom([] (const particle_comb_t& p) {
//                    // require 2 PID entries for the eta' candidate
//                    if (std::count_if(p.Photons.begin(), p.Photons.end(),
//                                      [](TParticlePtr g){ return g->Candidate->VetoEnergy; }) < 2)
//                        return true;
//                    return false;
//                }, "2 PIDs");

//        if (selection.empty())
//            continue;
//        fill_steps("Selection");

        //begin of test to only use one combination with 3CB and 1TAPS cluster
        size_t nCB = 0, nTAPS = 0;
        count_clusters(cands, nCB, nTAPS);
        if (nCB != 3 || nTAPS != 1)
            return;
        fill_steps("3CB&1TAPS");

        fake_comb_t comb;
        comb.reset();

        for (auto c : cands.get_iter())
            if (c->Detector & Detector_t::Type_t::CB)
                comb.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, c));
            else if (c->Detector & Detector_t::Type_t::TAPS)
                comb.Proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, c);

        comb.calc_values(taggerhit);

        // prefilter events
        // check if MM within defined window
        if (!ParticleTypeDatabase::Proton.GetWindow(settings.mm_window_size).Round().Contains(comb.MissingMass))
            continue;
        fill_steps("MM prefilter");
        // check if there are at least 2 PID entries
        if (std::count_if(comb.Photons.begin(), comb.Photons.end(),
                          [](TParticlePtr g){ return g->Candidate->VetoEnergy; }) < 2)
            continue;
        fill_steps("2 PIDs prefilter");
        // check proton cone prediction
        const double theta_sigma = 1.9;
        const double phi_sigma = 10.2;
        const auto miss_momentum = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass()) - comb.PhotonSum;
        //if (std_ext::radian_to_degree(comb.Proton->Angle(miss_momentum)) > 5)
        if (std_ext::radian_to_degree(std_ext::abs_diff(comb.Proton->Theta(), miss_momentum.Theta())) > 2*theta_sigma
                || std_ext::radian_to_degree(std_ext::abs_diff(comb.Proton->Phi(), miss_momentum.Phi())) > 2*phi_sigma)
            continue;
        fill_steps("Proton Cone");
        //test end

        // find best combination for each Tagger hit
        best_prob_fit = -std_ext::inf;
        // #combinations: binomial coefficient (n\\k)
        vector<double> IM_2g(3, std_ext::NaN);
        vector<double> IM_2g_fit(3, std_ext::NaN);

        //for (const auto& cand : selection) {
        for (const auto& cand : {comb}) {
            // do the fitting and check if the combination is better than the previous best
            if (!doFit_checkProb(taggerhit, cand, h, sig, best_prob_fit))
                continue;

            sig.DiscardedEk = cand.DiscardedEk;

            // run a kinematic fit with lepton candidates treated as charged pions and set the probability in the to-be-written tree
            sig.prob_antiPionFit = anti_pion_fit(taggerhit, cand);
            h.antiPionProb->Fill(sig.prob_antiPionFit);

            if (settings.less_plots())
                continue;
            utils::ParticleTools::FillIMCombinations(IM_2g.begin(), 2, cand.Photons);
            utils::ParticleTools::FillIMCombinations(IM_2g_fit.begin(), 2, sig.photons_kinfitted());
        }

        // only fill tree if a valid combination for the current Tagger hit was found
        if (!isfinite(best_prob_fit))
            continue;

        if (!settings.less_plots()) {
            for (const auto& im : IM_2g)
                h_subIM_2g->Fill(im, sig.TaggW);
            for (const auto& im : IM_2g_fit)
                h_subIM_2g_fit->Fill(im, sig.TaggW);
        }

        sig.Tree->Fill();
        fill_steps("Tree filled");
        h.true_rec_ZVertex->Fill(sig.kinfit_ZVertex, sig.trueZVertex);
    }

    if (!isfinite(best_prob_fit))
        return;
    fill_steps("best comb");

    h_counts->Fill(chan_id.decaystring.c_str(), 1);

    if (settings.less_plots())
        return;

    // histograms to investigate deviations between Pluto and Geant as well as MC and data
    {
        const auto etapMC = utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree);
        if (etapMC) {
            const double dE = etapMC->E - sig.etap().E();
            const double theta = std_ext::radian_to_degree(etapMC->Theta());
            h_energy_deviation->Fill(dE);
            h_fsClE_vs_pluto_geant_dE->Fill(dE, etapMC->E - etapMC->M());
            h_theta_vs_vz->Fill(sig.kinfit_ZVertex, theta);
            h_theta_vs_pluto_geant_dE->Fill(dE, theta);
            h_vz_vs_pluto_geant_dE->Fill(dE, sig.kinfit_ZVertex);
            h_delta_vz_vs_pluto_geant_dE->Fill(dE, sig.trueZVertex - sig.kinfit_ZVertex);
        } else
            LOG_N_TIMES(1, WARNING) << "(MC debug hists) No eta' found, only eta' decays used for MC investigation";

        // do some matching if we have the right particle tree
        if (signalMC) {
            TParticleList mctrue(utils::ParticleTools::FindParticles(ParticleTypeDatabase::eCharged, particletree));
            mctrue.emplace_back(utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon, particletree));
            // first for the reconstructed particles
            const auto matched = utils::match1to1(mctrue, cands.get_ptr_list(),
                                                  [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
                return p1->Angle(*p2); },
                                                  IntervalD(0., std_ext::degree_to_radian(15.)));
            // then for the fitted ones
            TParticleList fitted_photons;
            for (const auto& p : sig.photons_kinfitted())
                fitted_photons.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, p));
            const auto matched_fit = utils::match1to1(mctrue, fitted_photons,
                                                      TParticle::CalcAngle,
                                                      IntervalD(0., std_ext::degree_to_radian(15.)));

            const auto em = utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus, particletree);
            const auto ep = utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus, particletree);
            const auto g = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon, particletree);

            if (matched.size() == mctrue.size()) {
                const auto matched_g = utils::FindMatched(matched, g);
                const auto matched_em = utils::FindMatched(matched, em);
                const auto matched_ep = utils::FindMatched(matched, ep);

                // energy resolution
                h_energy_resolution_g->Fill(g->Ek() - matched_g->CaloEnergy);
                h_energy_resolution_em->Fill(em->Ek() - matched_em->CaloEnergy);
                h_energy_resolution_ep->Fill(ep->Ek() - matched_ep->CaloEnergy);
                h_energy_resolution_vs_theta_g->Fill(std_ext::radian_to_degree(matched_g->Theta),
                                                     g->Ek() - matched_g->CaloEnergy);
                h_energy_resolution_vs_theta_em->Fill(std_ext::radian_to_degree(matched_em->Theta),
                                                      em->Ek() - matched_em->CaloEnergy);
                h_energy_resolution_vs_theta_ep->Fill(std_ext::radian_to_degree(matched_ep->Theta),
                                                      ep->Ek() - matched_ep->CaloEnergy);
                h_energy_resolution_vs_trueE_g->Fill(g->Ek(), g->Ek() - matched_g->CaloEnergy);
                h_energy_resolution_vs_trueE_em->Fill(em->Ek(), em->Ek() - matched_em->CaloEnergy);
                h_energy_resolution_vs_trueE_ep->Fill(ep->Ek(), ep->Ek() - matched_ep->CaloEnergy);
                // theta resolution
                h_theta_resolution_g->Fill(std_ext::radian_to_degree(g->Theta() - matched_g->Theta));
                h_theta_resolution_em->Fill(std_ext::radian_to_degree(em->Theta() - matched_em->Theta));
                h_theta_resolution_ep->Fill(std_ext::radian_to_degree(ep->Theta() - matched_ep->Theta));
                h_theta_resolution_vs_energy_g->Fill(matched_g->CaloEnergy,
                                                     std_ext::radian_to_degree(g->Theta() - matched_g->Theta));
                h_theta_resolution_vs_energy_em->Fill(matched_em->CaloEnergy,
                                                      std_ext::radian_to_degree(em->Theta() - matched_em->Theta));
                h_theta_resolution_vs_energy_ep->Fill(matched_ep->CaloEnergy,
                                                      std_ext::radian_to_degree(ep->Theta() - matched_ep->Theta));
                h_theta_resolution_vs_trueTheta_g->Fill(std_ext::radian_to_degree(g->Theta()),
                                                        std_ext::radian_to_degree(g->Theta() - matched_g->Theta));
                h_theta_resolution_vs_trueTheta_em->Fill(std_ext::radian_to_degree(em->Theta()),
                                                         std_ext::radian_to_degree(g->Theta() - matched_g->Theta));
                h_theta_resolution_vs_trueTheta_ep->Fill(std_ext::radian_to_degree(ep->Theta()),
                                                         std_ext::radian_to_degree(g->Theta() - matched_g->Theta));
            } else
                LOG_N_TIMES(10, WARNING) << "(MC debug hists) Couldn't match all reconstructed FS particles";

            if (matched_fit.size() == mctrue.size()) {
                const auto matched_g = utils::FindMatched(matched_fit, g);
                const auto matched_em = utils::FindMatched(matched_fit, em);
                const auto matched_ep = utils::FindMatched(matched_fit, ep);

                // energy resolution
                h_energy_resolution_g_fit->Fill(g->Ek() - matched_g->Ek());
                h_energy_resolution_em_fit->Fill(em->Ek() - matched_em->Ek());
                h_energy_resolution_ep_fit->Fill(ep->Ek() - matched_ep->Ek());
                h_energy_resolution_vs_theta_g_fit->Fill(std_ext::radian_to_degree(matched_g->Theta()),
                                                         g->Ek() - matched_g->Ek());
                h_energy_resolution_vs_theta_em_fit->Fill(std_ext::radian_to_degree(matched_em->Theta()),
                                                          em->Ek() - matched_em->Ek());
                h_energy_resolution_vs_theta_ep_fit->Fill(std_ext::radian_to_degree(matched_ep->Theta()),
                                                          ep->Ek() - matched_ep->Ek());
                h_energy_resolution_vs_trueE_g_fit->Fill(g->Ek(), g->Ek() - matched_g->Ek());
                h_energy_resolution_vs_trueE_em_fit->Fill(em->Ek(), em->Ek() - matched_em->Ek());
                h_energy_resolution_vs_trueE_ep_fit->Fill(ep->Ek(), ep->Ek() - matched_ep->Ek());
                // theta resolution
                h_theta_resolution_g_fit->Fill(std_ext::radian_to_degree(g->Theta() - matched_g->Theta()));
                h_theta_resolution_em_fit->Fill(std_ext::radian_to_degree(em->Theta() - matched_em->Theta()));
                h_theta_resolution_ep_fit->Fill(std_ext::radian_to_degree(ep->Theta() - matched_ep->Theta()));
                h_theta_resolution_vs_energy_g_fit->Fill(matched_g->Ek(),
                                                         std_ext::radian_to_degree(g->Theta() - matched_g->Theta()));
                h_theta_resolution_vs_energy_em_fit->Fill(matched_em->Ek(),
                                                          std_ext::radian_to_degree(em->Theta() - matched_em->Theta()));
                h_theta_resolution_vs_energy_ep_fit->Fill(matched_ep->Ek(),
                                                          std_ext::radian_to_degree(ep->Theta() - matched_ep->Theta()));
                h_theta_resolution_vs_trueTheta_g_fit->Fill(std_ext::radian_to_degree(g->Theta()),
                                                            std_ext::radian_to_degree(g->Theta() - matched_g->Theta()));
                h_theta_resolution_vs_trueTheta_em_fit->Fill(std_ext::radian_to_degree(em->Theta()),
                                                             std_ext::radian_to_degree(g->Theta() - matched_g->Theta()));
                h_theta_resolution_vs_trueTheta_ep_fit->Fill(std_ext::radian_to_degree(ep->Theta()),
                                                             std_ext::radian_to_degree(g->Theta() - matched_g->Theta()));
            } else
                LOG_N_TIMES(10, WARNING) << "(MC debug hists) Couldn't match all fitted FS particles";
        } else
            LOG_N_TIMES(1, INFO) << "(MC debug hists) No eta' Dalitz decay, skip particle matching";
    }

    auto get_veto_energies = [] (vector<TSimpleParticle> particles)
    {
        vector<double> veto_energies;
        for (const auto& p : particles)
            veto_energies.emplace_back(p.VetoE);

        return veto_energies;
    };

    const auto etap_fs = sig.photons();
    const auto veto_energies = get_veto_energies(etap_fs);

    // get sorted indices of the eta' final state according to their Veto energies
    const auto sorted_idx = std_ext::get_sorted_indices_desc(veto_energies);

    // do an anti pi0 cut on the combinations e+g and e-g
    // (assuming the photon deposited the least energy in the PIDs)
    if (Cuts_t::ANTI_PI0_CUT) {
        const interval<double> pion_cut(Cuts_t::ANTI_PI0_LOW, Cuts_t::ANTI_PI0_HIGH);
        LorentzVec pi0;
        const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};
        for (const auto pi0_comb : pi0_combs) {
            pi0 = LorentzVec({0,0,0,0});
            for (const auto idx : pi0_comb)
                pi0 += TParticle(ParticleTypeDatabase::Photon, etap_fs.at(sorted_idx[idx]));
            // apply an anti pion cut
            if (pion_cut.Contains(pi0.M()))
                return;
        }
        h.steps->Fill("anti #pi^{0} cut", 1);
    }

    const TSimpleParticle proton = sig.p;
    LorentzVec etap({0,0,0,0});
    for (const auto& g : etap_fs)
        etap += TParticle(ParticleTypeDatabase::Photon, g);
    h.etapIM_cand->Fill(etap.M());
    h_protonVeto->Fill(proton.VetoE);
    h_pTheta->Fill(std_ext::radian_to_degree(proton.Theta()));

    ///\todo: still needed here? check final selection via ProtonPhotonCombs
    const auto idx1 = sorted_idx[0];
    const auto idx2 = sorted_idx[1];
    const auto l1 = etap_fs.at(idx1);
    const auto l2 = etap_fs.at(idx2);
    // suppress conversion decays
    if (sig.photons_vetoChannel().at(idx1) == sig.photons_vetoChannel().at(idx2))
        return;
    fill_steps("distinct PID");
    const double eeIM = (TParticle(ParticleTypeDatabase::eMinus, l1)
                         + TParticle(ParticleTypeDatabase::eMinus, l2)).M();
    h_IM2d->Fill(etap.M(), eeIM);
    h.IM2d->Fill(etap.M(), eeIM);

    // test effective cluster radius to distinguish between leptons and charged pions
    double effective_radius = sig.photons_effect_radius().at(idx1);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l1.Energy(), effective_radius);
    }
    effective_radius = sig.photons_effect_radius().at(idx2);
    if (isfinite(effective_radius)) {
        h.effect_rad->Fill(effective_radius);
        h.effect_rad_E->Fill(l2.Energy(), effective_radius);
    }
    double lateral_moment = sig.photons_lat_moment().at(idx1);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l1.Energy(), lateral_moment);
    }
    lateral_moment = sig.photons_lat_moment().at(idx2);
    if (isfinite(lateral_moment)) {
        h.lat_moment->Fill(lateral_moment);
        h.lat_moment_E->Fill(l2.Energy(), lateral_moment);
    }

    // test cluster size compared to energy
    h.cluster_size->Fill(l1.ClusterSize);
    h.cluster_size->Fill(l2.ClusterSize);
    h.cluster_size_E->Fill(l1.Energy(), l1.ClusterSize);
    h.cluster_size_E->Fill(l2.Energy(), l2.ClusterSize);

    h.etapIM_final->Fill(etap.M());
    h_etapIM_final->Fill(etap.M());
    h.hCopl_final->Fill(std_ext::radian_to_degree(abs(etap.Phi() - proton.Phi())) - 180.);
}

void EtapDalitzMC::Finish()
{
    if (settings.less_plots() || settings.reference_only())
        return;

    //TH1D* h_copy = HistFac.clone(h_IMee, "h_IMee_copy");
    //h_IMee_fraction3CB1TAPS->Sumw2();
    h_IMee_fraction3CB1TAPS_acceptance->Divide(h_IMee_acceptance);
    h_IMee_acceptance->Divide(h_IMee_total);
    h_IMee_fraction3CB1TAPS_total->Divide(h_IMee_total);
    h_IMee_fraction3CB1TAPS_Trigger4Cl->Divide(h_IMee_Trigger4Cl);
}

void EtapDalitzMC::ShowResult()
{
    if (settings.reference_only())
        return;

    for (auto& entry : channels)
        entry.second.Show();

    if (settings.less_plots())
        return;

    canvas(GetName()) << drawoption("colz") << h_IM2d << endc;

    canvas(GetName() + ": Efficiency 3CB1TAPS")
            << h_IMee_fraction3CB1TAPS_total
            << h_IMee_fraction3CB1TAPS_Trigger4Cl << endc;

    canvas(GetName() + ": Geometrical Acceptance")
            << h_IMee_acceptance
            << h_IMee_fraction3CB1TAPS_acceptance << endc;

    canvas(GetName() + ": Energy Resolution eta' FS")
            << h_energy_resolution_g
            << h_energy_resolution_em
            << h_energy_resolution_ep
            << drawoption("colz")
            << h_energy_resolution_vs_theta_g
            << h_energy_resolution_vs_theta_em
            << h_energy_resolution_vs_theta_ep
            << h_energy_resolution_vs_trueE_g
            << h_energy_resolution_vs_trueE_em
            << h_energy_resolution_vs_trueE_ep << endc;
    canvas(GetName() + ": Energy Resolution Fitted Particles")
            << h_energy_resolution_g_fit
            << h_energy_resolution_em_fit
            << h_energy_resolution_ep_fit
            << drawoption("colz")
            << h_energy_resolution_vs_theta_g_fit
            << h_energy_resolution_vs_theta_em_fit
            << h_energy_resolution_vs_theta_ep_fit
            << h_energy_resolution_vs_trueE_g_fit
            << h_energy_resolution_vs_trueE_em_fit
            << h_energy_resolution_vs_trueE_ep_fit << endc;

    canvas(GetName() + ": Theta Resolution eta' FS")
            << h_theta_resolution_g
            << h_theta_resolution_em
            << h_theta_resolution_ep
            << drawoption("colz")
            << h_theta_resolution_vs_energy_g
            << h_theta_resolution_vs_energy_em
            << h_theta_resolution_vs_energy_ep
            << h_theta_resolution_vs_trueTheta_g
            << h_theta_resolution_vs_trueTheta_em
            << h_theta_resolution_vs_trueTheta_ep << endc;
    canvas(GetName() + ": Theta Resolution Fitted Particles")
            << h_theta_resolution_g_fit
            << h_theta_resolution_em_fit
            << h_theta_resolution_ep_fit
            << drawoption("colz")
            << h_theta_resolution_vs_energy_g_fit
            << h_theta_resolution_vs_energy_em_fit
            << h_theta_resolution_vs_energy_ep_fit
            << h_theta_resolution_vs_trueTheta_g_fit
            << h_theta_resolution_vs_trueTheta_em_fit
            << h_theta_resolution_vs_trueTheta_ep_fit << endc;

//    list<TH1*> hists;
//    for (auto& entry : channels) {
//        hists.push_back(entry.second.proton_E_theta);
//    }

//    hists.sort([](const TH1* a, const TH1* b) {return a->GetEntries() > b->GetEntries();});

//    int i=0;
//    for (auto& h : hists) {
//        c << h;
//        i++;
//        if (i>=9)
//            break;
//    }

//    c << endc;
}

bool EtapDalitzMC::doFit_checkProb(const TTaggerHit& taggerhit,
                                   const particle_comb_t& comb,
                                   PerChannel_t& h,
                                   SigTree_t& t,
                                   double& best_prob_fit)
{
    LorentzVec etap_kinfit({0,0,0,0});
    LorentzVec etap_treefit({0,0,0,0});
    LorentzVec etap_kinfit_freeZ({0,0,0,0});
    LorentzVec etap_treefit_freeZ({0,0,0,0});

    LorentzVec etap = comb.PhotonSum;
    h.etapIM->Fill(etap.M(), t.TaggW);
    h_etapIM_vs_IMee->Fill(imee, etap.M(), t.TaggW);
    h.MM->Fill(comb.MissingMass, t.TaggW);
    h_MM_vs_IMee->Fill(imee, comb.MissingMass, t.TaggW);

    const auto fill_steps = [&] (const string& step) {
        h.steps->Fill(step.c_str(), 1);
        if (!settings.less_plots())
            h.steps_vs_IMee->Fill(imee, step.c_str(), 1);
    };


    /* start with the kinematic fitting */

    // treefit
    APLCON::Result_t treefit_result;

    treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    // works this way because only one combination needs to be fitted
    treefitter_etap.NextFit(treefit_result);

    if (settings.use_treefit) {
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        fill_steps("treefit");
    }

    // treefit free Z vertex
    APLCON::Result_t treefit_freeZ_result;

    treefitter_etap_freeZ.PrepareFits(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    treefitter_etap_freeZ.NextFit(treefit_freeZ_result);


    // kinfit
    auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

    if (!settings.use_treefit) {
        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            return false;
        fill_steps("kinfit");
    }

    // kinfit free Z vertex
    auto kinfit_freeZ_result = kinfit_freeZ.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);


    const double kinfit_prob = kinfit_result.Probability;
    const double treefit_prob = treefit_result.Probability;

    h.treefitChi2->Fill(treefit_result.ChiSquare);
    h.treefitProb->Fill(treefit_prob);
    h.treefitIter->Fill(treefit_result.NIterations);
    h.treefit_freeZ_chi2->Fill(treefit_freeZ_result.ChiSquare);
    h.treefit_freeZ_prob->Fill(treefit_freeZ_result.Probability);
    h.treefit_freeZ_iter->Fill(treefit_freeZ_result.NIterations);
    h.kinfitChi2->Fill(kinfit_result.ChiSquare);
    h.kinfitProb->Fill(kinfit_prob);
    h.kinfitIter->Fill(kinfit_result.NIterations);
    h.kinfit_freeZ_chi2->Fill(kinfit_freeZ_result.ChiSquare);
    h.kinfit_freeZ_prob->Fill(kinfit_freeZ_result.Probability);
    h.kinfit_freeZ_iter->Fill(kinfit_freeZ_result.NIterations);

    // determine which probability should be used to find the best candidate combination
    const double prob = settings.use_treefit ? treefit_prob : kinfit_prob;

    if (Cuts_t::PROBABILITY_CUT) {
        if (prob < Cuts_t::PROBABILITY)
            return false;
        fill_steps("probability");
    }

    if (!settings.less_plots()) {
        h_etap->Fill(etap.E - etap.M(), std_ext::radian_to_degree(etap.Theta()), t.TaggW);
        h_proton->Fill(comb.Proton->E - comb.Proton->M(), std_ext::radian_to_degree(comb.Proton->Theta()), t.TaggW);
    }

    // check if a better probability has been found
    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;


    /* Gather relevant fitter information and update branches */

    t.set_proton_information(comb.Proton);
    t.set_photon_information(comb.Photons, true);  // store lateral moment and effective cluster radius
    t.set_additional_photon_information(comb.Photons);

    t.p_effect_radius = tools.effective_radius(comb.Proton->Candidate);
    t.p_lat_moment    = tools.lat_moment(comb.Proton->Candidate);

    t.etap = etap;
    t.mm   = comb.MissingMass;
    t.copl = std_ext::radian_to_degree(abs(etap.Phi() - comb.Proton->Phi())) - 180.;

    // now handle the different fitted particle information separately

    // kinfit
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : kinfit.GetFittedPhotons())
            etap_kinfit += *g;
        h.etapIM_kinfit->Fill(etap_kinfit.M(), t.TaggW);
        h_etapIM_fitted_vs_IMee->Fill(imee, etap_kinfit.M(), t.TaggW);

        h.kinfit_ZVertex->Fill(kinfit.GetFittedZVertex());

        // update tree branches
        t.set_kinfit_information(kinfit, kinfit_result);
        t.etap_kinfit = etap_kinfit;
    }

    // treefit
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : treefitter_etap.GetFittedPhotons())
            etap_treefit += *g;
        if (!settings.less_plots())
            h.etapIM_treefit->Fill(etap_treefit.M(), t.TaggW);

        h.treefit_ZVertex->Fill(treefitter_etap.GetFittedZVertex());

        // update tree branches
        t.set_treefit_information(treefitter_etap, treefit_result);
        t.etap_treefit = etap_treefit;
    }

    // handling free Z vertex fits starting here
    if (kinfit_freeZ_result.Status != APLCON::Result_Status_t::Success
            && treefit_freeZ_result.Status != APLCON::Result_Status_t::Success)
        return true;

    // update tree branches
    t.set_fit_freeZ_results(kinfit_freeZ, treefitter_etap_freeZ,
                            kinfit_freeZ_result, treefit_freeZ_result);

    // kinfit with free Z vertex
    if (kinfit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit_freeZ.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : kinfit_freeZ.GetFittedPhotons())
            etap_kinfit_freeZ += *g;
        if (!settings.less_plots())
            h.etapIM_kinfit_freeZ->Fill(etap_kinfit_freeZ.M(), t.TaggW);

        h.kinfit_freeZ_ZVertex->Fill(kinfit_freeZ.GetFittedZVertex());

        t.etap_kinfit_freeZ = etap_kinfit_freeZ;
    }

    // treefit with free Z vertex
    if (treefit_freeZ_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap_freeZ.GetFitParticles().size() == settings.n_final_state);

        for (const auto& g : treefitter_etap_freeZ.GetFittedPhotons())
            etap_treefit_freeZ += *g;
        if (!settings.less_plots())
            h.etapIM_treefit_freeZ->Fill(etap_treefit_freeZ.M(), t.TaggW);

        h.treefit_freeZ_ZVertex->Fill(treefitter_etap_freeZ.GetFittedZVertex());

        t.etap_treefit_freeZ = etap_treefit_freeZ;
    }

    return true;
}

double EtapDalitzMC::anti_pion_fit(const TTaggerHit& taggerhit, const particle_comb_t& comb)
{
    fake_comb_t cand;
    cand.reset();
    cand.Proton = comb.Proton;

    auto leptons = get_sorted_indices_vetoE(comb.Photons);

    assert(leptons.size() == comb.Photons.size());

    // test the hypothesis of the two possible photon clusters with the highest veto energy to be charged pions
    cand.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::PiCharged, comb.Photons.at(leptons[0])->Candidate));
    cand.Photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::PiCharged, comb.Photons.at(leptons[1])->Candidate));
    cand.Photons.emplace_back(comb.Photons.at(leptons[2]));

    cand.calc_values(taggerhit);

    auto anti_fit_result = kinfit.DoFit(taggerhit.PhotonEnergy, cand.Proton, cand.Photons);

    if (anti_fit_result.Status != APLCON::Result_Status_t::Success)
        return -std_ext::inf;

    return anti_fit_result.Probability;
}

EtapDalitzMC::PerChannel_t EtapDalitzMC::manage_channel_histograms_get_current(const bool MC, const TEvent& event)
{
    // check if the current production is known already, create new HistogramFactory otherwise
    auto prod = productions.find(chan_id.production);
    if (prod == productions.end()) {
        auto hf = new HistogramFactory(chan_id.production, HistFac, "");
        productions.insert({chan_id.production, *hf});
    }
    prod = productions.find(chan_id.production);
    auto hf = prod->second;

    // check if the decay channel is known already, if not insert it
    auto c = channels.find(chan_id.decaystring);
    if (c == channels.end())
        channels.insert({chan_id.decaystring, PerChannel_t(chan_id.decay_name, chan_id.decaystring, hf)});

    c = channels.find(chan_id.decaystring);
    if (MC && !Settings_t::get().less_plots())
        c->second.Fill(event.MCTrue());

    // return the histogram struct for the current channel
    return c->second;
}

unsigned EtapDalitzMC::ReactionChannelList_t::identify(const ant::TParticleTree_t& tree) const
{
    if (!tree)
        return 0;

    for (const auto& c : channels) {

        if (!c.second.tree)
            continue;

        if (tree->IsEqual(c.second.tree, utils::ParticleTools::MatchByParticleName))
            return c.first;
    }

    return other_index;
}

const EtapDalitzMC::ReactionChannelList_t EtapDalitzMC::reaction_channels = EtapDalitz::makeChannels();

const unsigned EtapDalitzMC::ReactionChannelList_t::other_index = EtapDalitz::ReactionChannelList_t::other_index;


/* Reference channel analysis */
Etap2gMC::Etap2gMC(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    less_plots(opts->Get<bool>("less_plots", 0)),
    model_MC(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 make_shared<utils::UncertaintyModels::FitterSergey>()
                 )),
    kinfit(model_MC,
           opts->HasOption("SigmaZ"), EtapDalitz::MakeFitSettings(20)
           ),
    treefitter_etap(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g),
                    model_MC, opts->HasOption("SigmaZ"), {}, EtapDalitz::MakeFitSettings(20)
                    )
{
    ept = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    mc.CreateBranches(HistFac.makeTTree("MC"));

    h_taggChannel_vs_trueIM = HistFac.makeTH2D("EPT Channel vs. true IM", "IM(#gamma#gamma) [MeV]", "EPT Channel",
                                               BinSettings(50, ParticleTypeDatabase::EtaPrime.GetWindow(10)),  // only small IM range needed, true MC
                                               BinSettings(ept->GetNChannels()), "h_taggCh_vs_trueIM");

    if (!less_plots) {
        const auto energy = BinSettings(240, 0, 1200);
        const auto theta = BinSettings(160);
        const auto e_sigma = BinSettings(300, -150, 150);
        const auto theta_sigma = BinSettings(200, -50, 50);
        // histograms to check devitation between true and reconstructed MC events (resolution)
        // energy-related
        h_energy_resolution_vs_theta_g1 = HistFac.makeTH2D("Energy Deviation vs. #vartheta #gamma_{1} Pluto - Geant", "#vartheta [#circ]", "#sigmaE [MeV]",
                                                          theta, e_sigma, "h_MCsigmaE_theta_g1");
        h_energy_resolution_vs_theta_g2 = HistFac.makeTH2D("Energy Deviation vs. #vartheta #gamma_{2} Pluto - Geant", "#vartheta [#circ]", "#sigmaE [MeV]",
                                                           theta, e_sigma, "h_MCsigmaE_theta_g2");
        h_energy_resolution_vs_trueE_g1 = HistFac.makeTH2D("Energy Deviation vs. E_{true} #gamma_{1} Pluto - Geant", "E_{true} [MeV]", "#sigmaE [MeV]",
                                                           energy, e_sigma, "h_MCsigmaE_trueE_g1");
        h_energy_resolution_vs_trueE_g2 = HistFac.makeTH2D("Energy Deviation vs. E_{true} #gamma_{2} Pluto - Geant", "E_{true} [MeV]", "#sigmaE [MeV]",
                                                           energy, e_sigma, "h_MCsigmaE_trueE_g2");
        h_energy_resolution_vs_theta_g1_fit = HistFac.makeTH2D("Energy Deviation vs. #vartheta fitted #gamma_{1} Pluto - Geant",
                                                               "#vartheta [#circ]", "#sigmaE [MeV]", theta, e_sigma, "h_MCsigmaE_theta_g1_fit");
        h_energy_resolution_vs_theta_g2_fit = HistFac.makeTH2D("Energy Deviation vs. #vartheta fitted #gamma_{2} Pluto - Geant",
                                                               "#vartheta [#circ]", "#sigmaE [MeV]", theta, e_sigma, "h_MCsigmaE_theta_g2_fit");
        h_energy_resolution_vs_trueE_g1_fit = HistFac.makeTH2D("Energy Deviation vs. E_{true} fitted #gamma_{1} Pluto - Geant",
                                                               "E_{true} [MeV]", "#sigmaE [MeV]", energy, e_sigma, "h_MCsigmaE_trueE_g1_fit");
        h_energy_resolution_vs_trueE_g2_fit = HistFac.makeTH2D("Energy Deviation vs. E_{true} fitted #gamma_{2} Pluto - Geant",
                                                               "E_{true} [MeV]", "#sigmaE [MeV]", energy, e_sigma, "h_MCsigmaE_trueE_g2_fit");

        // theta-related
        h_theta_resolution_vs_energy_g1 = HistFac.makeTH2D("Theta Deviation vs. Energy #gamma_{1} Pluto - Geant", "E [MeV]", "#sigma#vartheta [#circ]",
                                                           energy, theta_sigma, "h_MCsigmaTheta_E_g1");
        h_theta_resolution_vs_energy_g2 = HistFac.makeTH2D("Theta Deviation vs. Energy #gamma_{2} Pluto - Geant", "E [MeV]", "#sigma#vartheta [#circ]",
                                                           energy, theta_sigma, "h_MCsigmaTheta_E_g2");
        h_theta_resolution_vs_trueTheta_g1 = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} #gamma_{1} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                              "#sigma#vartheta [#circ]", theta, theta_sigma, "h_MCsigmaTheta_trueTheta_g1");
        h_theta_resolution_vs_trueTheta_g2 = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} #gamma_{2} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                             "#sigma#vartheta [#circ]", theta, theta_sigma, "h_MCsigmaTheta_trueTheta_g2");
        h_theta_resolution_vs_energy_g1_fit = HistFac.makeTH2D("Theta Deviation vs. Energy fitted #gamma_{1} Pluto - Geant",
                                                               "E [MeV]", "#sigma#vartheta [#circ]", energy, theta_sigma, "h_MCsigmaTheta_E_g1_fit");
        h_theta_resolution_vs_energy_g2_fit = HistFac.makeTH2D("Theta Deviation vs. Energy fitted #gamma_{2} Pluto - Geant",
                                                               "E [MeV]", "#sigma#vartheta [#circ]", energy, theta_sigma, "h_MCsigmaTheta_E_g2_fit");
        h_theta_resolution_vs_trueTheta_g1_fit = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} fitted #gamma_{1} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                                  "#sigma#vartheta [#circ]", theta, theta_sigma, "h_MCsigmaTheta_trueTheta_g1_fit");
        h_theta_resolution_vs_trueTheta_g2_fit = HistFac.makeTH2D("Theta Deviation vs. #vartheta_{true} fitted #gamma_{2} Pluto - Geant", "#vartheta_{true} [#circ]",
                                                                  "#sigma#vartheta [#circ]", theta, theta_sigma, "h_MCsigmaTheta_trueTheta_g2_fit");
    } else
        LOG(INFO) << "Less histograms will be created and stored";

    if (opts->HasOption("SigmaZ")) {
        double sigma_z = 0.;
        std_ext::copy_if_greater(sigma_z, opts->Get<double>("SigmaZ", 0.));
        LOG(INFO) << "Fit Z vertex enabled with sigma = " << sigma_z;
        kinfit.SetZVertexSigma(sigma_z);
        treefitter_etap.SetZVertexSigma(sigma_z);
    }
}

void Etap2gMC::ProcessEvent(const TEvent& event, manager_t&)
{
    Process(event);
}

void Etap2gMC::Process(const TEvent& event)
{
    // only process MC
    if (!t->MCtrue)
        return;

    const auto& cands = event.Reconstructed().Candidates;
    const auto& particletree = event.MCTrue().ParticleTree;

    // check if the current event is the reference
    const bool refMC = particletree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g),
                                             utils::ParticleTools::MatchByParticleName);

    TParticleList etap_photons;  // used to store true photons, make sure the used order is preserved everywhere
    // work only with reference MC
    if (refMC) {
        // obtain eta' final state photons from tree using general approach
        particletree->Map_nodes([&etap_photons] (const TParticleTree_t& t) {
            const auto& parent = t->GetParent();
            if (!parent)
                return;
            if (parent->Get()->Type() == ParticleTypeDatabase::EtaPrime) {
                etap_photons.push_back(t->Get());
            }
        });

        // get all photons (should be two for reference)
        TParticleList mctrue(utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, particletree));
        assert(mctrue.size() == 2);

        h_taggChannel_vs_trueIM->Fill((*mctrue.front() + *mctrue.back()).M(),
                                      // max. one true Tagger hit in case of MC, but could be empty if outside of covered energy range
                                      // --> fill underflow bin in this case
                                      event.MCTrue().TaggerHits.empty() ? -1 : event.MCTrue().TaggerHits.front().Channel);

        mc.opening = std_ext::radian_to_degree(TParticle::CalcAngle(mctrue.front(), mctrue.back()));

        // fetch the recoil proton
        mctrue.emplace_back(utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree));

        // store true information in the tree
        for (const auto& tp : mctrue) {
            mc.names().emplace_back(tp->Type().Name());
            mc.energies_true().push_back(tp->Ek());
            mc.thetas_true().push_back(std_ext::radian_to_degree(tp->Theta()));
            mc.phis_true().push_back(std_ext::radian_to_degree(tp->Phi()));
        }

        // do some matching
        const auto matched = utils::match1to1(mctrue, cands.get_ptr_list(),
                                              [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
            return p1->Angle(*p2); },
                                              IntervalD(0., std_ext::degree_to_radian(15.)));

        if (matched.size() == mctrue.size()) {
            const auto matched_g1 = utils::FindMatched(matched, etap_photons.front());
            const auto matched_g2 = utils::FindMatched(matched, etap_photons.back());
            const auto matched_p = utils::FindMatched(matched, mctrue.back());

            // store matched information in the tree
            for (const auto& p : {matched_g1, matched_g2, matched_p}) {
                mc.energies().push_back(p->CaloEnergy);
                mc.thetas().push_back(std_ext::radian_to_degree(p->Theta));
                mc.phis().push_back(std_ext::radian_to_degree(p->Phi));
            }
        }
    }

    mc.multiplicity = cands.size();

    size_t nCB = 0, nTAPS = 0;
    count_clusters(cands, nCB, nTAPS);
    mc.nCB = nCB;
    mc.nTAPS = nTAPS;

    // up to this point are no selection or cut criteria applied
    // save the MC tree here to have all MC information available
    mc.fillAndReset();

    triggersimu.ProcessEvent(event);

    if (!triggersimu.HasTriggered())
        return;

    if (t->nCands != Etap2g::Cuts_t::N_FINAL_STATE)
        return;

    /* test combinations to determine best proton candidate */
    if (Etap2g::TEST_COMBS) {
        utils::ProtonPhotonCombs proton_photons(cands);
        double best_prob_fit = -std_ext::inf;

        // loop over all tagger hits
        for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {
            promptrandom->SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
            if (promptrandom->State() == PromptRandom::Case::Outside)
                continue;

            t->TaggW = promptrandom->FillWeight();
            t->TaggE = taggerhit.PhotonEnergy;
            t->TaggT = taggerhit.Time;
            t->TaggCh = taggerhit.Channel;

            // find best combination for each Tagger hit
            best_prob_fit = -std_ext::inf;

            particle_combs_t selection = proton_photons()
                    .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(Etap2g::Cuts_t::MM_WINDOW_SIZE).Round())  // MM cut on proton mass
                    .FilterCustom([=] (const particle_comb_t& p) {
                        // ensure the possible proton candidate is kinematically allowed
                        if (std_ext::radian_to_degree(p.Proton->Theta()) > Etap2g::Cuts_t::MAX_PROTON_THETA)
                            return true;
                        return false;
                    }, "proton #vartheta");
//                    .FilterCustom([] (const particle_comb_t& p) {
//                        // require that PID entries do not have more than .3 MeV
//                        if (std::count_if(p.Photons.begin(), p.Photons.end(),
//                                          [](TParticlePtr g){ return g->Candidate->VetoEnergy < .3; }) < 2)
//                            return true;
//                        return false;
//                    }, "PID energy");

            if (selection.empty())
                continue;

            for (const auto& cand : selection) {
                // do the fitting and check if the combination is better than the previous best
                if (!doFit_checkProb(taggerhit, cand.Proton, cand.Photons, best_prob_fit))
                    continue;
            }

            // only fill tree if a valid combination for the current Tagger hit was found
            if (!isfinite(best_prob_fit))
                continue;

            t->Tree->Fill();
        }

        return;
    }


    /* use simple selection requiring proton in TAPS and 2 photons in CB */
    TParticlePtr proton;
    TParticleList photons;

    if (!simple2CB1TAPS(cands, proton, photons))
        return;

    for (const TTaggerHit& taggerhit : event.Reconstructed().TaggerHits) {  // loop over all tagger hits
        t->reset();
        promptrandom->SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if (promptrandom->State() == PromptRandom::Case::Outside)
            continue;

        t->TaggW = promptrandom->FillWeight();
        t->TaggE = taggerhit.PhotonEnergy;
        t->TaggT = taggerhit.Time;
        t->TaggCh = taggerhit.Channel;

        /* kinematic fitting */
        // treefit
        APLCON::Result_t treefit_result;

        treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

        treefitter_etap.NextFit(treefit_result);  // no loop, just one possible combination

        if (Etap2g::USE_TREEFIT && treefit_result.Status != APLCON::Result_Status_t::Success)
            return;

        // kinfit

        auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

        if (Etap2g::USE_KINFIT && kinfit_result.Status != APLCON::Result_Status_t::Success)
            return;

        // do some matching with the fitted particles if it's the reference particle tree
        if (refMC && !less_plots) {
            // get all photons
            TParticleList mctrue(utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, particletree));
            assert(mctrue.size() == 2);

            // first for the reconstructed particles
            const auto matched = utils::match1to1(mctrue, photons,
                                                  TParticle::CalcAngle,
                                                  IntervalD(0., std_ext::degree_to_radian(15.)));
            // then for the fitted ones
            const auto matched_fit = utils::match1to1(mctrue, t->photons_kinfitted(),
                                                      [] (const TParticlePtr& p1, const TLorentzVector& p2) {
                    return p1->Angle(p2); },
                                                      IntervalD(0., std_ext::degree_to_radian(15.)));

            TParticleList fitted_photons;
            transform(etap_photons.begin(), etap_photons.end(), back_inserter(fitted_photons),
                      [matched_fit] (const TParticlePtr& p) {
                return make_shared<TParticle>(ParticleTypeDatabase::Photon, utils::FindMatched(matched_fit, p));
            });

            const auto& g1 = etap_photons.front();
            const auto& g2 = etap_photons.back();

            if (matched.size() == mctrue.size()) {
                const auto matched_g1 = utils::FindMatched(matched, g1);
                const auto matched_g2 = utils::FindMatched(matched, g2);

                // energy resolution
                h_energy_resolution_vs_theta_g1->Fill(std_ext::radian_to_degree(matched_g1->Candidate->Theta),
                                                      g1->Ek() - matched_g1->Candidate->CaloEnergy);
                h_energy_resolution_vs_theta_g2->Fill(std_ext::radian_to_degree(matched_g2->Candidate->Theta),
                                                      g2->Ek() - matched_g2->Candidate->CaloEnergy);
                h_energy_resolution_vs_trueE_g1->Fill(g1->Ek(), g1->Ek() - matched_g1->Candidate->CaloEnergy);
                h_energy_resolution_vs_trueE_g2->Fill(g2->Ek(), g2->Ek() - matched_g2->Candidate->CaloEnergy);
                // theta resolution
                h_theta_resolution_vs_energy_g1->Fill(matched_g1->Candidate->CaloEnergy,
                                                      std_ext::radian_to_degree(g1->Theta() - matched_g1->Candidate->Theta));
                h_theta_resolution_vs_energy_g2->Fill(matched_g2->Candidate->CaloEnergy,
                                                      std_ext::radian_to_degree(g2->Theta() - matched_g2->Candidate->Theta));
                h_theta_resolution_vs_trueTheta_g1->Fill(std_ext::radian_to_degree(g1->Theta()),
                                                         std_ext::radian_to_degree(g1->Theta() - matched_g1->Candidate->Theta));
                h_theta_resolution_vs_trueTheta_g2->Fill(std_ext::radian_to_degree(g2->Theta()),
                                                         std_ext::radian_to_degree(g2->Theta() - matched_g2->Candidate->Theta));
            } else
                LOG_N_TIMES(10, WARNING) << "(MC ref debug hists) Couldn't match all reconstructed FS particles";

            if (matched_fit.size() == mctrue.size()) {
                const auto matched_g1 = fitted_photons.front();
                const auto matched_g2 = fitted_photons.back();

                // energy resolution
                h_energy_resolution_vs_theta_g1_fit->Fill(std_ext::radian_to_degree(matched_g1->Theta()),
                                                          g1->Ek() - matched_g1->Ek());
                h_energy_resolution_vs_theta_g2_fit->Fill(std_ext::radian_to_degree(matched_g2->Theta()),
                                                          g2->Ek() - matched_g2->Ek());
                h_energy_resolution_vs_trueE_g1_fit->Fill(g1->Ek(), g1->Ek() - matched_g1->Ek());
                h_energy_resolution_vs_trueE_g2_fit->Fill(g2->Ek(), g2->Ek() - matched_g2->Ek());
                // theta resolution
                h_theta_resolution_vs_energy_g1_fit->Fill(matched_g1->Ek(),
                                                          std_ext::radian_to_degree(g1->Theta() - matched_g1->Theta()));
                h_theta_resolution_vs_energy_g2_fit->Fill(matched_g2->Ek(),
                                                          std_ext::radian_to_degree(g2->Theta() - matched_g2->Theta()));
                h_theta_resolution_vs_trueTheta_g1_fit->Fill(std_ext::radian_to_degree(g1->Theta()),
                                                             std_ext::radian_to_degree(g1->Theta() - matched_g1->Theta()));
                h_theta_resolution_vs_trueTheta_g2_fit->Fill(std_ext::radian_to_degree(g2->Theta()),
                                                             std_ext::radian_to_degree(g2->Theta() - matched_g2->Theta()));
            } else
                LOG_N_TIMES(10, WARNING) << "(MC ref debug hists) Couldn't match all fitted FS particles";
        }

        // fill the tree with the fitted values
        fill_tree(treefit_result, kinfit_result, proton, photons);
        t->Tree->Fill();
    }
}

void Etap2gMC::fill_tree(const APLCON::Result_t& treefit_result,
                         const APLCON::Result_t& kinfit_result,
                         const TParticlePtr proton,
                         const TParticleList& photons)
{
    LorentzVec etap;
    LorentzVec etap_kinfit;
    LorentzVec etap_treefit;

    /* check if the fits converged and fill the trees accordingly */

    // update branches with general particle and fitter information
    etap = sumlv(photons.begin(), photons.end());

    t->set_proton_information(proton);
    t->set_photon_information(photons, true);  // store lateral moment and effective cluster radius

    t->etap = etap;

    // kinfit
    if (kinfit_result.Status == APLCON::Result_Status_t::Success) {
        assert(kinfit.GetFitParticles().size() == Etap2g::Cuts_t::N_FINAL_STATE);

        auto kinfit_photons = kinfit.GetFittedPhotons();

        etap_kinfit = sumlv(kinfit_photons.begin(), kinfit_photons.end());

        // update tree branches
        t->set_kinfit_information(kinfit, kinfit_result);
        t->etap_kinfit = etap_kinfit;
    }

    // treefit
    if (treefit_result.Status == APLCON::Result_Status_t::Success) {
        assert(treefitter_etap.GetFitParticles().size() == Etap2g::Cuts_t::N_FINAL_STATE);

        auto treefit_photons = treefitter_etap.GetFittedPhotons();

        etap_treefit = sumlv(treefit_photons.begin(), treefit_photons.end());

        // update tree branches
        t->set_treefit_information(treefitter_etap, treefit_result);
        t->etap_treefit = etap_treefit;
    }
}

bool Etap2gMC::simple2CB1TAPS(const TCandidateList& cands,
                              TParticlePtr& proton,
                              TParticleList& photons)
{
    size_t nCB = 0, nTAPS = 0;
    photons.clear();

    for (auto p : cands.get_iter())
        if (p->Detector & Detector_t::Type_t::CB && nCB++ < 2)
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, p));
        else if (p->Detector & Detector_t::Type_t::TAPS && nTAPS++ < 1)
            proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, p);
        else if (nCB > 2 || nTAPS > 1)
            return false;

    return true;
}

bool Etap2gMC::doFit_checkProb(const TTaggerHit& taggerhit,
                               const TParticlePtr proton,
                               const TParticleList& photons,
                               double& best_prob_fit)
{
    LorentzVec etap({0,0,0,0});

    for (const auto& g : photons)
        etap += *g;

    /* kinematical checks to reduce computing time */
    const interval<double> mm = ParticleTypeDatabase::Proton.GetWindow(Etap2g::Cuts_t::MM_WINDOW_SIZE);

    LorentzVec missing = taggerhit.GetPhotonBeam() + LorentzVec::AtRest(ParticleTypeDatabase::Proton.Mass());
    missing -= etap;
    if (!mm.Contains(missing.M()))
        return false;


    /* now start with the kinematic fitting */
    // treefit
    APLCON::Result_t treefit_result;

    treefitter_etap.PrepareFits(taggerhit.PhotonEnergy, proton, photons);

    // works this way because only one combination needs to be fitted
    while (treefitter_etap.NextFit(treefit_result))
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            continue;

    if (Etap2g::USE_TREEFIT)
        if (treefit_result.Status != APLCON::Result_Status_t::Success)
            return false;

    // kinfit

    auto kinfit_result = kinfit.DoFit(taggerhit.PhotonEnergy, proton, photons);

    if (!Etap2g::USE_TREEFIT)
        if (kinfit_result.Status != APLCON::Result_Status_t::Success)
            return false;


    // determine which probability should be used to find the best candidate combination
    const double prob = Etap2g::USE_TREEFIT ? treefit_result.Probability : kinfit_result.Probability;

    if (!std_ext::copy_if_greater(best_prob_fit, prob))
        return false;

    fill_tree(treefit_result, kinfit_result, proton, photons);

    return true;
}

void Etap2gMC::setPromptRandom(PromptRandom::Switch& prs)
{
    promptrandom = &prs;
}

void Etap2gMC::linkTree(RefTree_t& ref)
{
    t = &ref;
}


AUTO_REGISTER_PHYSICS(EtapDalitzMC)
