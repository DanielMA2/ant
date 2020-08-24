#include "hstack.h"

#include "base/std_ext/string.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/container.h"
#include "base/std_ext/system.h"
#include "base/Logger.h"

#include "tree/stream_TBuffer.h"

#include "TBrowser.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "TDirectory.h"
#include "TPaveText.h"
#include "TVirtualX.h"

#include "cereal/access.hpp"

#include <list>
#include <iomanip>
#include <algorithm>

using namespace ant;
using namespace std;

hstack::options_t hstack::GlobalOptions;

hstack::hstack(const string& name, const std::string& title, bool simple_) :
    THStack(name.c_str(), title.c_str()),
    origtitle(title),
    simple(simple_)
{
    // simulate behaviour of TH1
    if(!simple)
        gDirectory->Append(this);

    // figure out cutnames now, as later we have
    // no idea where gDirectory points to (in TBrowser)
    // we traverse the directory structure upwards
    // to find the cutnames. A directory name is a cutname if
    // it contains an "h" directory on the same level (see cuttree::StackedHists)
    auto currDir = gDirectory;
    while(currDir) {
        auto parentDir = currDir->GetMotherDir();
        if(parentDir && parentDir->GetDirectory("h")) {
            cutnames.push_front(currDir->GetName());
        }
        currDir = parentDir;
    }
}

hstack::hstack() : THStack() {}

hstack::~hstack() {}

TH1* hstack::hist_t::GetPtr(const string& path)
{
    auto ptr = dynamic_cast<TH1*>(gDirectory->Get(path.c_str()));
    if(ptr == nullptr)
        VLOG(5) << "Could not Get() TH1* for path " << path;
    return ptr;
}

string hstack::hist_t::GetPath(const TH1* ptr)
{
    const std::string path = std_ext::formatter()
                             << ptr->GetDirectory()->GetPath()
                             << "/" << ptr->GetName();
    // strip the file of the full path
    const auto n = path.find_first_of(':');
    if(n == path.npos)
        return path;
    return path.substr(n+1);
}


hstack& hstack::operator<<(TH1* hist)
{
    if(simple)
        Add(hist, current_option.DrawOption.c_str());
    hists.emplace_back(hist, current_option);
    current_option.Z = 0; // reset z after each add
    return *this;
}

hstack& hstack::operator<<(const drawoption& c)
{
    current_option.DrawOption = c.Option();
    return *this;
}

hstack& hstack::operator<<(const ModOption_t& option)
{
    current_option = option;
    return *this;
}

namespace ant {

ostream& operator<<(ostream& s, const hstack& o)
{
    s << "ant::hstack: "  << o.GetName() << ": " << o.GetTitle() << '\n';
    s << "Global MC Scaling = " << o.GlobalOptions.MCScale << '\n';
    /// \todo print current options...?
    for(const auto& h : o.hists) {
        s << "  " << h.Path << '\n';
    }
    s << endl;
    return s;
}

} // namespace ant

void hstack::Print(const char*) const
{
    cout << *this;
}

void hstack::Print() const
{
    Print("");
}

void hstack::checkHists()
{
    auto it_hist = hists.begin();
    while(it_hist != hists.end()) {
        hist_t hist = *it_hist;
        if(hist.Ptr == nullptr) {
            // try to get it now, might be that after merge
            // the Ptr is not yet initialized
            hist.Ptr = hist_t::GetPtr(hist.Path);
            // erase if still unsuccessful
            if(hist.Ptr == nullptr) {
                it_hist = hists.erase(it_hist);
                continue;
            }
        }
        it_hist++;
    }
}



bool hstack::IsCompatible(const hstack& other) const
{
    // we do not check hists here, since that's what could be merged
    return string(GetName()) == string(other.GetName()) &&
            string(GetTitle()) == string(other.GetTitle());
}

bool hstack::hist_operation(const hstack& other, function<bool(TH1*, const TH1*)> operation)
{
    // If this function is called directly from a ROOT shell, the logger is not initialized.
    // This leads to a segfault if LOG(LEVEL) is called. Within Ant itself the logger is
    // properly set up. To prevent segfaults in case of issuing LOG(ERROR), do a simple
    // check if the session type is interactive and call SetupLogger() to init it.
    // Calling it multiple times does no harm.
    ///\todo think of a better way to identify if a ROOT shell is started
    if (std_ext::system::isInteractive())
        SetupLogger();

    if (!IsCompatible(other)) {
        LOG(ERROR) << "Provided hstack objects do not seem to be compatible";
        return false;
    }

    if (hists.size() != other.hists.size()) {
        LOG(ERROR) << "Amount of histograms contained in hstack objecs does not match";
        return false;
    }

    bool success = true;
    for (auto& this_it : hists) {
        const auto& other_it = find_if(other.hists.begin(), other.hists.end(),
                                       [&this_it] (const hist_t& h) {
            return strcmp(this_it.Ptr->GetTitle(), h.Ptr->GetTitle()) == 0;
        });
        if (other_it != other.hists.end())
            success = success && operation(this_it.Ptr, other_it->Ptr);
        else {
            LOG(ERROR) << "Histogram " << this_it.Ptr->GetTitle() << " not found in hstack";
            success = false;
        }
    }

    return success;
}

bool hstack::Add(const hstack& other, const double scaling)
{
    auto add = [scaling](TH1* h1, const TH1* h2){
        return h1->Add(h2, scaling);
    };

    return hist_operation(other, add);
}

bool hstack::Divide(const hstack& other)
{
    auto divide = [](TH1* h1, const TH1* h2){
        return h1->Divide(h2);
    };

    return hist_operation(other, divide);
}

bool hstack::Multiply(const hstack& other)
{
    auto multiply = [](TH1* h1, const TH1* h2){
        return h1->Multiply(h2);
    };

    return hist_operation(other, multiply);
}

void hstack::Rebin(const int group) {
    for (auto& it : hists)
        it.Ptr = it.Ptr->Rebin(group);
}

template<typename Att>
void AttCopy(const Att* src, Att* dest) {
    src->Copy(*dest);
}

bool hstack::hist_t::isDataHist() const
{
    /// \todo find better way to detect which histograms
    /// should NOT be scaled or drawn with error bar in legend.
    /// For now, we assume that data is always drawn
    /// with error bars.
    return std_ext::contains(Option.DrawOption, "E");
}

string hstack::hist_t::getTitleKey() const
{
    // it is expected that this method returns some unique short string
    // for each histogram. When stuff is build with the cuttree, this is indeed the case
    const string path = Ptr->GetDirectory()->GetPath();
    const auto n_slash = path.find_last_of('/');
    if(n_slash == string::npos)
        return Ptr->GetTitle();
    return path.substr(n_slash+1);
}

void hstack::buildIntelliTitle()
{
    // remove cutnames from title (we need origtitle for that)
    {
        auto title_parts = std_ext::tokenize_string(origtitle, ": ");
        for(auto it = title_parts.begin(); it != title_parts.end(); ) {
            if(std_ext::contains(cutnames, *it)) {
                it = title_parts.erase(it);
            }
            else {
                ++it;
            }
        }
        SetTitle(std_ext::concatenate_string(title_parts, ": ").c_str());
    }

    const auto height = 0.03*cutnames.size();

    // ensure intellititle is there
    if(!intellititle) {
        const auto& p = GlobalOptions.LegendPosition;
        intellititle = std_ext::make_unique<TPaveText>(
                           p.Start().Start(),
                           p.Start().Stop()-height,
                           p.Stop().Start(),
                           p.Start().Stop(),
                           "NDC"
                           );
        intellititle->SetFillColor(kWhite);
        intellititle->SetBorderSize(1);
        for(auto& cutname : cutnames)
            intellititle->AddText(cutname.c_str());
    }

    if(intellilegend) {
        const auto x1 = intellilegend->GetX1NDC();
        const auto y1 = intellilegend->GetY1NDC()-height;
        const auto x2 = intellilegend->GetX2NDC();
        const auto y2 = intellilegend->GetY1NDC();
        intellititle->SetX1NDC(x1);
        intellititle->SetY1NDC(y1);
        intellititle->SetX2NDC(x2);
        intellititle->SetY2NDC(y2);
    }
}

void hstack::updateIntelliLegend(TLegend& legend, std::list<wraphist_t> wraphists)
{

    if(GlobalOptions.FixLegendPosition) {
        const auto& p = GlobalOptions.LegendPosition;
        legend.SetX1NDC(p.Start().Start());
        legend.SetY1NDC(p.Start().Stop());
        legend.SetX2NDC(p.Stop().Start());
        legend.SetY2NDC(p.Stop().Stop());
    }

    legend.Clear();

    for(const auto& wraphist : wraphists) {
        const hist_t& hist = wraphist.Hist;
        auto title = hist.getTitleKey();
        if(hstack::GlobalOptions.ShowEntriesInLegend)
            title += std_ext::formatter() << " (" << setprecision(3) << wraphist.Entries << ")";
        auto entry = legend.AddEntry((TObject*)0, title.c_str(), hist.isDataHist() ? "lpfe" : "lpf");
        AttCopy<TAttLine>(hist.Ptr, entry);
        AttCopy<TAttFill>(hist.Ptr, entry);
        AttCopy<TAttMarker>(hist.Ptr, entry);
        if(hist.Option.BkgColor>=0) {
            // during legend building, modify fill color
            entry->SetFillColor(hist.Option.BkgColor);
            entry->SetFillStyle(1001);
            entry->SetLineWidth(2);
        }
    }
}

void hstack::Paint(const char* chopt)
{
    if(hists.empty())
        return;

    // ensure the histograms are there
    checkHists();

    if(simple) {
        THStack::Paint(chopt);
        return;
    }

    UpdateMCScaling();

    if(GlobalOptions.UseIntelliTitle)
        SetTitle("");
    else
        SetTitle(origtitle.c_str());

    if(GlobalOptions.YAxisRange.IsSane()) {
        SetMinimum(GlobalOptions.YAxisRange.Start());
        SetMaximum(GlobalOptions.YAxisRange.Stop());
    }
    else {
        // resets it to default -1111, well, THStack knows what to do then
        SetMinimum();
        SetMaximum();
    }

    // make list of wrapped hists, in order to kick out elements and
    // count entries...
    list<wraphist_t> tmp_hists(hists.begin(), hists.end());

    // remove all empty/hidden hists if requested
    const auto basehist = fHistogram ? fHistogram : hists.front().Ptr;
    const auto xaxis = basehist->GetXaxis();
    auto it_hist = tmp_hists.begin();
    while(it_hist != tmp_hists.end()) {
        const hist_t hist = it_hist->Hist;
        auto h = hist.Ptr;
        it_hist->Entries = h->GetDimension() > 1 ?
                               h->GetEntries() : h->Integral(xaxis->GetFirst(), xaxis->GetLast());
        const auto empty = GlobalOptions.IgnoreEmptyHist && it_hist->Entries < 1 && it_hist->Entries > -1;
        const auto hidden = GlobalOptions.tryGetHist(hist.getTitleKey()).Hidden;
        if(empty || hidden)
            it_hist = tmp_hists.erase(it_hist);
        else
            ++it_hist;
    }

    // if all were removed, use all again
    // as the THStack really looks weird if nothing is drawn
    if(tmp_hists.empty()) {
        tmp_hists.insert(tmp_hists.end(), hists.begin(), hists.end());
    }


    // build the legend before re-ordering and after removal
    if(GlobalOptions.UseIntelliLegend) {
        if(!intellilegend) {
            const auto& p = GlobalOptions.LegendPosition;
            intellilegend = std_ext::make_unique<TLegend>(
                                p.Start().Start(),
                                p.Start().Stop(),
                                p.Stop().Start(),
                                p.Stop().Stop()
                                );
        }
        // sort by name if it's Bkg_
        // little hack to get legends stable...
        tmp_hists.sort([] (const wraphist_t& a, const wraphist_t& b) {
            const hist_t& hist_a = a.Hist;
            const hist_t& hist_b = b.Hist;
            if(hist_a.Option.Z == hist_b.Option.Z) {
                const string title_a = hist_a.Ptr->GetTitle();
                const string title_b = hist_b.Ptr->GetTitle();
                if(std_ext::contains(title_a,"Bkg_") && std_ext::contains(title_b,"Bkg_"))
                    return title_a < title_b;
                else
                    return false;
            }
            return false;
        });
        updateIntelliLegend(*intellilegend, tmp_hists);
    }
    else {
        intellilegend = nullptr;
    }

    if(GlobalOptions.UseIntelliTitle) {
        buildIntelliTitle();
    }
    else {
        intellititle = nullptr;
    }

    if(GlobalOptions.UseYAxisEntriesPerBin) {
        string xaxis_title(xaxis->GetTitle());
        string xaxis_unit;
        {
            auto pos_Slash = xaxis_title.find_last_of("/");
            if(pos_Slash != string::npos) {
                xaxis_unit = std_ext::string_sanitize(xaxis_title.substr(pos_Slash+1));
            }
        }
        auto bin_width = (xaxis->GetXmax() - xaxis->GetXmin())/xaxis->GetNbins();
        bin_width = std::round(bin_width*100)/100;
        string yaxis_title = std_ext::formatter() << "Entries / " << bin_width << " " << xaxis_unit;
        const auto yaxis = basehist->GetYaxis();
        yaxis->SetTitle(yaxis_title.c_str());
    }

    // reverse first
    std::reverse(tmp_hists.begin(), tmp_hists.end());

    // then sort according to Z
    tmp_hists.sort([] (const wraphist_t& a, const wraphist_t& b) {
        return a.Hist.get().Option.Z < b.Hist.get().Option.Z;
    });

    if(fHists)
        fHists->Clear("nodelete");
    for(const auto& o : tmp_hists) {
        const hist_t& hist = o.Hist;
        Add(hist.Ptr, hist.Option.DrawOption.c_str());
    }

    // let's hope that this does not throw an exception
    // before we unreverse the fHists
    string chopt_str(chopt);
    chopt_str += "nostack"; // always draw with nostack
    THStack::Paint(chopt_str.c_str());

    if(!tmp_hists.empty() && basehist != fHistogram) {
        const auto yaxis = basehist->GetYaxis();
        GetXaxis()->SetTitle(xaxis->GetTitle());
        GetYaxis()->SetTitle(yaxis->GetTitle());
        THStack::Paint(chopt_str.c_str());
    }

    if(intellilegend && !gPad->FindObject(intellilegend.get()))
        intellilegend->Draw();

    if(intellititle && !gPad->FindObject(intellititle.get()))
        intellititle->Draw();
}



void hstack::SetGlobalYAxisRange(double low, double high)
{
    GlobalOptions.YAxisRange = {low, high};
    gPad->Modified();
    gPad->Update();
}

void hstack::FixLegendPosition(bool flag)
{
    TIter next(gPad->GetListOfPrimitives());
    while(TObject* obj = next()) {
        if(obj->IsA() == TLegend::Class()) {
            TLegend* legend = dynamic_cast<TLegend*>(obj);
            GlobalOptions.LegendPosition = {
                {legend->GetX1NDC(), legend->GetY1NDC()},
                {legend->GetX2NDC(), legend->GetY2NDC()}
            };
        }
    }
    GlobalOptions.FixLegendPosition = flag;
}
bool hstack::GetFixLegendPosition() const { return GlobalOptions.FixLegendPosition; }

void hstack::UseIntelliLegend(bool flag) { GlobalOptions.UseIntelliLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliLegend() const { return GlobalOptions.UseIntelliLegend; }
void hstack::UseIntelliTitle(bool flag) { GlobalOptions.UseIntelliTitle = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseIntelliTitle() const { return GlobalOptions.UseIntelliTitle; }
void hstack::IgnoreEmptyHist(bool flag) { GlobalOptions.IgnoreEmptyHist = flag; gPad->Modified(); gPad->Update();}
bool hstack::GetIgnoreEmptyHist() const { return GlobalOptions.IgnoreEmptyHist; }
void hstack::ShowEntriesInLegend(bool flag) { GlobalOptions.ShowEntriesInLegend = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetShowEntriesInLegend() const { return GlobalOptions.ShowEntriesInLegend; }
void hstack::UseYAxisEntriesPerBin(bool flag) { GlobalOptions.UseYAxisEntriesPerBin = flag; gPad->Modified(); gPad->Update(); }
bool hstack::GetUseYAxisEntriesPerBin() const { return GlobalOptions.UseYAxisEntriesPerBin; }

void hstack::UpdateMCScaling()
{
    // gather summation targets already while scaling...
    // the value of the map represents the old histogram
    map<string, std::vector<double>> items;

    // do the scaling perHist and globally
    for(hist_t& hist : hists) {
        if(hist.isDataHist())
            continue;
        // have a look if this hist was scaled already
        const auto&  titlekey = hist.getTitleKey();
        const auto&  perHist = GlobalOptions.tryGetHist(titlekey);
        if(!perHist.AddTo.empty())
            items[perHist.AddTo] = {};
        const double new_scale = GlobalOptions.MCScale*perHist.Scale;

        // prevent loss of information...
        if(new_scale <= 0)
            continue;
        const double eff_scale = new_scale/GlobalOptions.AppliedScale[hist.Ptr];
        hist.Ptr->Scale(eff_scale);
        GlobalOptions.AppliedScale[hist.Ptr] = new_scale;
    }

    // clear summation targets, and build full items
    // that's ok because title keys are hopefully unique
    for(const hist_t& hist : hists) {
        const auto&  titlekey = hist.getTitleKey();
        auto it_item = items.find(titlekey);
        if(it_item != items.end()) {
            auto axis = hist.Ptr->GetXaxis();
            for(int bin = 0; bin<=axis->GetNbins()+1; bin++) {
                hist.Ptr->SetBinContent(bin, 0);
            }
            it_item->second.resize(axis->GetNbins()+2, 0);
        }
        else {
            auto& v = items[titlekey]; // always adds to items, never modifies
            auto axis = hist.Ptr->GetXaxis();
            v.resize(axis->GetNbins()+2, 0);
            for(int bin = 0; bin<=axis->GetNbins()+1; bin++) {
                v[bin] = hist.Ptr->GetBinContent(bin);
            }
        }
    }

    if(hists.size() != items.size()) {
        LOG(ERROR) << "Strange, probably weird histogram names";
        return;
    }

    // sum all way up, but use old values from items to avoid double summation
    // this is pretty slow but it works...too lazy to think about trees and stuff
    // and it might never stop if cycles are present
    for(const hist_t& hist : hists) {
        const auto&  titlekey = hist.getTitleKey();
        const auto&  perHist = GlobalOptions.tryGetHist(titlekey);
        auto addTo = perHist.AddTo;
        while(!addTo.empty()) {
            // find the hist where this should be added to
            auto it_hist = std::find_if(hists.begin(), hists.end(), [addTo] (const hist_t& h) {
                return h.getTitleKey() == addTo;
            });
            if(it_hist == hists.end()) {
                LOG(ERROR) << "Strange, could not find " << addTo;
                break;
            }
            auto axis = it_hist->Ptr->GetXaxis();
            const auto& this_hist = items.at(titlekey);
            for(int bin = 0; bin<=axis->GetNbins()+1; bin++) {
                it_hist->Ptr->SetBinContent(bin, this_hist.at(bin) + it_hist->Ptr->GetBinContent(bin));
            }
            // go to next add to
            addTo = GlobalOptions.tryGetHist(addTo).AddTo;
        }
    }


}

Long64_t hstack::Merge(TCollection* li, TFileMergeInfo*)
{
    if(!li)
        return 0;
    TIter next(li);

    while(auto h = dynamic_cast<hstack*>(next())) {
        if(!IsCompatible(*h)) {
            LOG(ERROR) << "Skipping incompatible hstack:\n "
                       << *this << "\n"
                       << *h;
            continue;
        }
        // we add non-existing paths
        auto have_path = [] (const hists_t& hists, const string& path) {
            for(const auto& hist : hists)
                if(hist.Path == path)
                    return true;
            return false;
        };
        for(const auto& hist : h->hists) {
            if(!have_path(hists, hist.Path))
                hists.emplace_back(hist);
        }
    }

    return 0;
}

#include "TGFrame.h"
#include "TGTableLayout.h"
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGComboBox.h"
#include "TExec.h"

namespace ant {

// wrapper copied/modified from calibration/gui/ManagerWindow.cc
struct LambdaExec : TExec {
    LambdaExec(function<void()> action_) : action(action_) {}
    virtual void Exec(const char*) override {
        action();
    }
    static void Connect(TQObject* w, const vector<string>& mysignals, function<void()> action_) {
        auto a = new LambdaExec(action_);
        for(auto& signal : mysignals)
            w->Connect(signal.c_str(), "TExec", a , "Exec(=\"\")");
    }

private:
    function<void()> action;
};

struct hstack_Menu : TGTransientFrame {


    hstack_Menu(const hstack& s) :
        TGTransientFrame(gClient->GetRoot())
    {
        SetCleanup(kDeepCleanup);
        // Set a name to the main frame
        SetWindowName("ant::hstack Menu");

        auto frame = new TGVerticalFrame(this);

        AddFrame(frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

        auto make_scaleNumEntry = [] (const TGWindow* p, double initval) {
            auto e = new TGNumberEntry(p, initval, 5, -1, TGNumberEntry::kNESReal, TGNumberEntry::kNEAPositive);
            string formatted = std_ext::formatter() << std::fixed << setprecision(2) << initval;
            e->SetText(formatted.c_str());
            return e;
        };

        // add some global MC scaling factor
        auto frame_1 = new TGHorizontalFrame(this);
        frame_1->AddFrame(new TGLabel(frame_1,"Global MC Scaling: "));
        auto entry_globalMCScaling = make_scaleNumEntry(frame_1, hstack::GlobalOptions.MCScale);
        LambdaExec::Connect(entry_globalMCScaling, {"ValueSet(Long_t)","ValueChanged(Long_t)"}, [entry_globalMCScaling] () {
            hstack::GlobalOptions.MCScale = entry_globalMCScaling->GetNumber();
            gPad->Modified();
            gPad->Update();
        });

        frame_1->AddFrame(entry_globalMCScaling);
        frame->AddFrame(frame_1);

        auto frame_table = new TGCompositeFrame(frame);
        // make three columns, one row for each hist plus title row
        frame_table->SetLayoutManager(new TGTableLayout(frame_table, s.hists.size()+1, 4));

        // title row
        frame_table->AddFrame(new TGLabel(frame_table, "Name"),  new TGTableLayoutHints(0,1,0,1));
        frame_table->AddFrame(new TGLabel(frame_table, "Hide"),  new TGTableLayoutHints(1,2,0,1));
        frame_table->AddFrame(new TGLabel(frame_table, "Scale"), new TGTableLayoutHints(2,3,0,1));
        frame_table->AddFrame(new TGLabel(frame_table, "AddTo"), new TGTableLayoutHints(3,4,0,1));

        // check if anything in AddTo was set already set (is not-empty)
        const auto addToModified = std::find_if(hstack::GlobalOptions.PerHist.begin(),
                                           hstack::GlobalOptions.PerHist.end(),
                                           [] (const pair<string, hstack::options_t::hist_t>& h) {
            return !h.second.AddTo.empty();
        }) != hstack::GlobalOptions.PerHist.end();


        for(int row = 1; row <= int(s.hists.size()); row++) {
            const auto& h = s.hists.at(row-1);

            const auto& titlekey = h.getTitleKey();

            // Name / TitleKey
            frame_table->AddFrame(new TGLabel(frame_table, titlekey.c_str()), new TGTableLayoutHints(0,1,row,row+1));

            // Hide button
            {
                auto btn = new TGCheckButton(frame_table);
                btn->SetState(hstack::GlobalOptions.PerHist[titlekey].Hidden ? kButtonDown : kButtonUp);
                LambdaExec::Connect(btn, {"Clicked()"}, [titlekey,btn] () {
                    hstack::GlobalOptions.PerHist[titlekey].Hidden = btn->IsOn();
                    gPad->Modified();
                    gPad->Update();
                });
                frame_table->AddFrame(btn, new TGTableLayoutHints(1,2,row,row+1));
            }

            if(h.isDataHist())
                continue;

            // Scale numentry
            {
                auto numentry = make_scaleNumEntry(frame_table, hstack::GlobalOptions.tryGetHist(titlekey).Scale);
                LambdaExec::Connect(numentry, {"ValueSet(Long_t)","ValueChanged(Long_t)"}, [titlekey,numentry] () {
                    hstack::GlobalOptions.PerHist[titlekey].Scale = numentry->GetNumber();
                    gPad->Modified();
                    gPad->Update();
                });
                frame_table->AddFrame(numentry, new TGTableLayoutHints(2,3,row,row+1));
            }

            // AddTo combobox
            {
                auto combo = new TGComboBox(frame_table);
                combo->NewEntry(""); // empty

                int i_Sum_MC = -1;
                int i_Bkg_MC = -1;
                for(auto i = 0u; i<s.hists.size(); ++i)  {
                    const auto& titlekey = s.hists[i].getTitleKey();
                    combo->NewEntry(titlekey.c_str());
                    if(titlekey == "Sum_MC")
                        i_Sum_MC = i;
                    if(titlekey == "Bkg_MC")
                        i_Bkg_MC = i;
                }
                combo->Resize(100, 20); // combobox is too stupid to resize automagically *sigh*

                if(!addToModified) {
                    // select reasonable defaults for comboboxes
                    // we assume that histgrams Sum_MC and Bkg_MC are the i-th and (i+1)-th element
                    // and after Bkg_MC only non-signal/reference channels are present
                    // (this is the default and recommended ordering of the CutTree)
                    if(i_Bkg_MC>=0 && i_Sum_MC>=0
                       && i_Bkg_MC-1 == i_Sum_MC
                       && !h.isDataHist()
                       && titlekey != "Sum_MC")
                    {

                        if(row-1<i_Sum_MC || titlekey == "Bkg_MC") {
                            combo->Select(i_Sum_MC+2);
                            hstack::GlobalOptions.PerHist[titlekey].AddTo = "Sum_MC";
                        }
                        else {
                            combo->Select(i_Bkg_MC+2);
                            hstack::GlobalOptions.PerHist[titlekey].AddTo = "Bkg_MC";
                        }
                    }
                }
                else {
                    const auto& addTo = hstack::GlobalOptions.tryGetHist(titlekey).AddTo;
                    for(auto i=0u;i<s.hists.size();i++) {
                        if(s.hists.at(i).getTitleKey() == addTo) {
                            combo->Select(i+2);
                            break;
                        }
                    }
                }

                // add the callback
                LambdaExec::Connect(combo, {"Selected(Int_t)"}, [titlekey,combo] () {
                    auto e = dynamic_cast<TGTextLBEntry*>(combo->GetSelectedEntry());
                    hstack::GlobalOptions.PerHist[titlekey].AddTo = e->GetText()->GetString();
                    gPad->Modified();
                    gPad->Update();
                });

                frame_table->AddFrame(combo, new TGTableLayoutHints(3,4,row,row+1));
            }
        }

        frame_table->Resize();
        frame->AddFrame(frame_table, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

        // set focus
        gVirtualX->SetInputFocus(GetId());

        // Map all subwindows of main frame
        MapSubwindows();
        Resize(GetDefaultSize()); // this is used here to init layout algorithm
        MapWindow();
    }

};

} // namespace ant

void hstack::OpenStackMenu()
{
    // the most ugliest way to make the window appear...
    new hstack_Menu(*this);
}



// cereal business

template<class Archive>
void save(Archive& archive, const TNamed& m)
{
    archive(string(m.GetName()), string(m.GetTitle()));
}

template<class Archive>
void load(Archive& archive, TNamed& m)
{
    string name;
    string title;
    archive(name, title);
    m.SetName(name.c_str());
    m.SetTitle(title.c_str());
}

namespace cereal
{
template <class Archive>
struct specialize<Archive, hstack, cereal::specialization::member_serialize> {};
}

void hstack::Streamer(TBuffer& R__b)
{
    stream_TBuffer::DoBinary(R__b, *this);
}






