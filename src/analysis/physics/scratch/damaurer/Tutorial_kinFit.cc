#include "Tutorial_kinFit.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


scratch_damaurer_Tutorial_kinFit::scratch_damaurer_Tutorial_kinFit(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    const BinSettings bins_nClusters(20);
    const BinSettings tagger_time_bins(2000, -200, 200);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // our binnings, may write directly BinSettings(10) here
                                   "h_nClusters"         // ROOT object name, auto-generated if omitted
                                   );

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("t"));
}

void scratch_damaurer_Tutorial_kinFit::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);
    h_nClusters->Fill(event.Reconstructed().Clusters.size());

    for(auto& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        const double TaggWeight = promptrandom.FillWeight();
        h_TaggerTime->Fill(taggerhit.Time, TaggWeight);

        t.TaggW = TaggWeight;
        t.nClusters = event.Reconstructed().Clusters.size();
        t.Tree->Fill();
    } 
}

void scratch_damaurer_Tutorial_kinFit::ShowResult()
{
    // ShowResult is called after processing of events has finished,
    // and interactive mode (aka non-batchmode) is chosen

    // ant::canvas nice wrapper around TCanvas
    ant::canvas(GetName()+": Basic plots")
            << h_nClusters
            << h_TaggerTime
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

void scratch_damaurer_Tutorial_kinFit::Finish()
{
    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;
}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)
AUTO_REGISTER_PHYSICS(scratch_damaurer_Tutorial_kinFit)
