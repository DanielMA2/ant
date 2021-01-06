#include "omega_Dalitz.h"
#include "base/WrapTFile.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "expconfig/detectors/TAPS.h"
#include "utils/uncertainties/Interpolated.h"

#include <iostream>
#include <memory>

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

APLCON::Fit_Settings_t scratch_damaurer_omega_Dalitz::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

scratch_damaurer_omega_Dalitz::scratch_damaurer_omega_Dalitz(const string& name, OptionsPtr opts) :
    Physics(name, opts),

    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                             utils::UncertaintyModels::Interpolated::Type_t::MC,
                             make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true))
{

    fitter.SetZVertexSigma(3.0);

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    const BinSettings statistic_bins((nrCuts_total)*10,0,nrCuts_total);
    const BinSettings bins_tagger_time(2000, -200, 200);
    const BinSettings bins_nClusters(20);
    const BinSettings bins_nCand(10);
    const BinSettings bins_nNeuCand(6);
    const BinSettings bins_nChaCand(6);
    const BinSettings bins_Veto_Energy(500, 0, 10);
    const BinSettings bins_Calo_Energy(500, 0, 1200);
    const BinSettings BeamE_bins(100,0, 1600);
    const BinSettings Im_pi0_bins(200, 0, 600);

    auto hf_RecData_CandStat = new HistogramFactory("hf_RecData_CandStat", HistFac, "");
    auto hf_Tagger = new HistogramFactory("hf_Tagger", HistFac, "");
    auto hf_CaloEvsVetoE = new HistogramFactory("hf_CaloEvsVetoE", HistFac, "");
    auto hf_Candidates = new HistogramFactory("hf_Candidates", HistFac, "");
    auto hf_VetoE = new HistogramFactory("hf_VetoE", HistFac, "");
    auto hf_BeamE = new HistogramFactory("hf_BeamE", HistFac, "");
    auto hf_2g_IM = new HistogramFactory("hf_2g_IM", HistFac, "");

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = hf_Tagger->makeTH1D("Number of clusters", "nClusters", "#", bins_nClusters, "h_nClusters", true);

    h_TaggerTime = hf_Tagger->makeTH1D("Tagger Time",     // title
                                    "t [ns]", "#",     // xlabel, ylabel
                                    bins_tagger_time,  // our binnings
                                    "h_TaggerTime", true     // ROOT object name, auto-generated if omitted
                                    );

    h_nCandidates = hf_Candidates->makeTH1D("Number of candidates", "nCandidates", "#", bins_nCand, "h_nCandidates", true);
    h_nNeuChaCandidates = hf_Candidates->makeTH2D("Number of charged vs neutral candidates", "nChaCand", "nNeuCand", bins_nChaCand, bins_nNeuCand, "h_nNeuChaCandidates", true);

    h_RecData_Stat = hf_RecData_CandStat->makeTH1D("Amount of events after cuts:",     // title
                                     "Cuts", "#",     // xlabel, ylabel
                                     statistic_bins,  // our binnings
                                     "h_RecData_Stat", true    // ROOT object name, auto-generated if omitted
                                     );

    for (unsigned int i=0; i<nrCuts_total; i++){

    h_AllCaloEvsVetoE_CB[i] = hf_CaloEvsVetoE->makeTH2D("Deposited calo vs veto energies in CB "+cuts[i],     // title
                                     "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                     bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                     "h_AllCaloEvsVetoE_CB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                     );

    h_AllCaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE->makeTH2D("Deposited calo vs veto energies in TAPS "+cuts[i],     // title
                                     "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                     bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                     "h_AllCaloEvsVetoE_TAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                     );

    h_AllVetoE_CB[i] = hf_VetoE->makeTH1D("Deposited veto energies in CB "+cuts[i],     // title
                                     "Veto E [MeV]", "#",  // xlabel, ylabel
                                     bins_Veto_Energy,    // our binnings
                                     "h_AllVetoE_CB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                     );

    h_AllVetoE_TAPS[i] = hf_VetoE->makeTH1D("Deposited veto energies in TAPS "+cuts[i],     // title
                                     "Veto E [MeV]", "#",  // xlabel, ylabel
                                     bins_Veto_Energy,    // our binnings
                                     "h_AllVetoE_TAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                     );

    h_beamE[i] = hf_BeamE->makeTH1D("Initial photon beam energy "+cuts[i],     // title
                                                         "E_{photonbeam} [MeV]", "#",     // xlabel, ylabel
                                                         BeamE_bins,  // our binnings
                                                         "h_beamE_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                         );

    }

    for (unsigned int i=0; i<(nrCuts_total-1); i++){

    h_2g_IM[i] = hf_2g_IM->makeTH1D("IM(2g) "+cuts[i+1],     // title
                                         "IM(2g) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2g_IM_"+cuts[i+1], true     // ROOT object name, auto-generated if omitted
                                         );
    }

    //hist = HistFac.makeTH1D(" Accepted Events", "step", "#", BinSettings(10), "steps");

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("t"));  

}

void scratch_damaurer_omega_Dalitz::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    // set fitter uncertainty models
    fitter.SetUncertaintyModel(fit_model);

    const auto& data = event.Reconstructed();
    const auto& clusters = data.Clusters;
    const auto& candidates = data.Candidates;

    //Set up particle combinations for the kinFit
    //utils::ProtonPhotonCombs proton_photons(candidates);
    //particle_combs_t protphotcombs = proton_photons();

    //Fill e.g. the polar angle distribution into a histogram
    TCandidatePtrList all;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    //TCandidatePtrList protonCand;

    //TParticleList photons;
    //TParticleList protons;

    vector<bool> allCanInCB;
    vector<bool> allCanInTAPS;
    vector<bool> chaCanInCB;
    vector<bool> chaCanInTAPS;
    vector<bool> neuCanInCB;
    vector<bool> neuCanInTAPS;

    vector<double> allCanCluSize;
    vector<double> allCanCaloE;
    vector<double> allCanVetoE;
    vector<double> allCanTime;
    vector<double> allThe;
    vector<double> allPhi;
    vector<double> neuCanCluSize;
    vector<double> neuCanCaloE;
    vector<double> neuCanVetoE;
    vector<double> neuCanTime;
    vector<double> neuThe;
    vector<double> neuPhi;
    vector<double> chaCanCluSize;
    vector<double> chaCanCaloE;
    vector<double> chaCanVetoE;
    vector<double> chaCanTime;
    vector<double> chaThe;
    vector<double> chaPhi;

    for(const auto& cand : candidates.get_iter()) {

        //h_VetoEnergies->Fill(cand->VetoEnergy);

        all.emplace_back(cand);
        allThe.push_back(cand->Theta*radtodeg);
        allPhi.push_back(cand->Phi*radtodeg);
        allCanCluSize.push_back(cand->ClusterSize);
        allCanCaloE.push_back(cand->CaloEnergy);
        allCanVetoE.push_back(cand->VetoEnergy);
        allCanTime.push_back(cand->Time);

        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            allCanInCB.push_back(true);
        }
        else{
            allCanInCB.push_back(false);
        }
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            allCanInTAPS.push_back(true);
        }
        else{
            allCanInTAPS.push_back(false);
        }

        if(cand->VetoEnergy <= vetoEthreshold) {
            //photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, cand));
            neutral.emplace_back(cand);
            neuThe.push_back(cand->Theta*radtodeg);
            neuPhi.push_back(cand->Phi*radtodeg);
            neuCanCluSize.push_back(cand->ClusterSize);
            neuCanCaloE.push_back(cand->CaloEnergy);
            neuCanVetoE.push_back(cand->VetoEnergy);
            neuCanTime.push_back(cand->Time);

            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                neuCanInCB.push_back(true);
            }
            else{
                neuCanInCB.push_back(false);
            }
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                neuCanInTAPS.push_back(true);
            }
            else{
                neuCanInTAPS.push_back(false);
            }
        }
        else{

            charged.emplace_back(cand);
            chaThe.push_back(cand->Theta*radtodeg);
            chaPhi.push_back(cand->Phi*radtodeg);
            chaCanCluSize.push_back(cand->ClusterSize);
            chaCanCaloE.push_back(cand->CaloEnergy);
            chaCanVetoE.push_back(cand->VetoEnergy);
            chaCanTime.push_back(cand->Time);

            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                chaCanInCB.push_back(true);
            }
            else{
                chaCanInCB.push_back(false);
            }
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                chaCanInTAPS.push_back(true);
            }
            else{
                chaCanInTAPS.push_back(false);
            }
            /*
            if(cand->Theta*radtodeg < 50){
                protonCand.emplace_back(cand);
                protons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, cand));
            }
            */
        }

    }

    TParticle g[neu_nrSel];
    TLorentzVector L2g;

    if(neutral.size() == neu_nrSel && charged.size() == cha_nrSel){

        //proton_toFit = protons.at(0);

        for(int i = 0; i<neu_nrSel; i++){
            g[i] = TParticle(ParticleTypeDatabase::Photon, neutral[i]);
            L2g += (TLorentzVector)(g[i]);
        }

    }

    h_nClusters->Fill(clusters.size());

    /*
    double best_probability = std_ext::NaN;
    double fit_z_vert = -50;
    double fitbeamE = -50;
    */

    //test

    int cut_ind = 0;
    for(auto& taggerhit : data.TaggerHits) {

        //Get corrected Taggertime
        double cortagtime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        promptrandom.SetTaggerTime(cortagtime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        const double TaggWeight = promptrandom.FillWeight();

        stat[cut_ind]+=TaggWeight;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);

        h_TaggerTime->Fill(taggerhit.Time, TaggWeight);
        h_nCandidates->Fill(candidates.size(), TaggWeight);

        h_nNeuChaCandidates->Fill(charged.size(), neutral.size(), TaggWeight);

        for (unsigned int i=0; i<all.size(); i++){
            if(allCanInCB[i]){
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_CB[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
            if(allCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_TAPS[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
        }

        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        if(!(neutral.size() == neu_nrSel && charged.size() == cha_nrSel))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;

        for (unsigned int i=0; i<all.size(); i++){
            if(allCanInCB[i]){
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_CB[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
            if(allCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_TAPS[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
        }

        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-1]->Fill(L2g.M(),TaggWeight);

        if(!(InitialPhotonVec.E()>=Omega_Ethreshold))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;

        for (unsigned int i=0; i<all.size(); i++){
            if(allCanInCB[i]){
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_CB[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
            if(allCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_TAPS[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
        }

        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-1]->Fill(L2g.M(),TaggWeight);

        if(!(L2g.M()>(mpi0-0.4*mpi0) && L2g.M()<(mpi0+0.4*mpi0)))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;

        for (unsigned int i=0; i<all.size(); i++){
            if(allCanInCB[i]){
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_CB[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
            if(allCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(allCanCaloE[i],allCanVetoE[i],TaggWeight);
                h_AllVetoE_TAPS[cut_ind]->Fill(allCanVetoE[i],TaggWeight);
            }
        }

        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-1]->Fill(L2g.M(),TaggWeight);

        t.TaggW = TaggWeight;
        t.nClusters = data.Clusters.size();
        t.Tree->Fill();

        /*
        //Particles for the later kinFit
        TParticlePtr proton;
        TParticleList photons;

        //Photon-proton comb:

        //Performing the kinFit
        APLCON::Result_t fitresult = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

        // check if the fit converged
        if (fitresult.Status != APLCON::Result_Status_t::Success)
            continue;

        // check if we found a better probability for this fit and copy it if true, continue otherwise
        if (!std_ext::copy_if_greater(best_probability, fitresult.Probability))
            continue;

        // retrieve the fitted photon and proton information as well as the number of iterations
        auto fitted_proton = fitter.GetFittedProton();
        auto fitted_photons = fitter.GetFittedPhotons();
        auto iterations = fitresult.NIterations;
        */


    } 

}

void scratch_damaurer_omega_Dalitz::ShowResult()
{

    // ShowResult is called after processing of events has finished,
    // and interactive mode (aka non-batchmode) is chosen

    // ant::canvas nice wrapper around TCanvas

    int lower_edge = 0;
    int upper_edge = nrCuts_total;
    int number_of_bins = nrCuts_total*10;
    int steps = (int)(number_of_bins/(upper_edge-lower_edge));

    for (unsigned int i=0; i<nrCuts_total; i++){
        cout << "Amount events after " << i << " applied cuts: " << stat[i] << endl;
        h_RecData_Stat->SetBinContent(i*steps+1, stat[i]);
        h_RecData_Stat->GetXaxis()->SetBinLabel(i*steps+1,cuts[i].c_str());
    }

    ant::canvas(GetName()+": Event statistics after applied cuts:")
            << h_RecData_Stat
            << endc; // actually draws the canvas

    ant::canvas c_Tagger_hists(GetName()+": Basic plots");
            c_Tagger_hists << h_TaggerTime;
            c_Tagger_hists << h_nClusters;
            c_Tagger_hists << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW");
            c_Tagger_hists << endc; // actually draws the canvas;

    ant::canvas c_nCandidates(GetName()+": Number of candidates");
            c_nCandidates << h_nCandidates;
            c_nCandidates << drawoption("pcolz");
            c_nCandidates << h_nNeuChaCandidates;
            c_nCandidates << endc; // actually draws the canvas;

    ant::canvas c_CaloVsVetoE(GetName()+": Calo vs Veto energies");
            c_CaloVsVetoE << drawoption("pcolz");
            for (unsigned int i=0; i<nrCuts_total; i++){
            c_CaloVsVetoE << h_AllCaloEvsVetoE_CB[i];
            c_CaloVsVetoE << h_AllCaloEvsVetoE_TAPS[i];
            }
            c_CaloVsVetoE << endc; // actually draws the canvas

    ant::canvas c_VetoE(GetName()+": Veto energies");
            for (unsigned int i=0; i<nrCuts_total; i++){
            c_VetoE << h_AllVetoE_CB[i];
            c_VetoE << h_AllVetoE_TAPS[i];
            }
            c_VetoE << endc; // actually draws the canvas

    ant::canvas c_BeamE(GetName()+": Initial photonbeam energy");
            c_BeamE << drawoption("HIST");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c_BeamE << h_beamE[i];
    }
            c_BeamE << endc; // actually draws the canvas

    ant::canvas c_2g_IM(GetName()+": IM(2neu) combinations");
    for (unsigned int i=0; i<(nrCuts_total-1); i++){
            c_2g_IM << h_2g_IM[i];
    }
            c_2g_IM << endc; // actually draws the canvas

}

void scratch_damaurer_omega_Dalitz::Finish()
{

    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;

}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)

AUTO_REGISTER_PHYSICS(scratch_damaurer_omega_Dalitz)
