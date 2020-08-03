#include "MC_pw_ppi0g_p3g.h"
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

#include <iostream>
#include <memory>

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_damaurer_MC_pw_ppi0g_p3g::scratch_damaurer_MC_pw_ppi0g_p3g(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    taps = make_shared<expconfig::detector::TAPS_2013_11>(false, false, false);

    const BinSettings tagger_time_bins(1000, -200, 200);
    const BinSettings Veto_Energy_bins(1000, 0, 10);
    const BinSettings Calo_Energy_bins(1000, 0, 1200);
    const BinSettings initialBeamE_bins(1000,1,1.6);
    const BinSettings cos_bins(100,-1,1);
    //const BinSettings cos_bins_BackToBack(1000,-1,1);
    const BinSettings sqrt_S_bins(100,s_square_min,s_square_max);
    const BinSettings phi_bins(100, -180, 180);

    const BinSettings theta_bins(180, 0, 180);
    const BinSettings proton_theta_bins(180, 0, 90);
    const BinSettings theta_bins_CB(140, 19, 160);
    const BinSettings theta_bins_TAPS(25, 0, 25);

    const BinSettings Im_proton_bins(100, 0, 1800);
    const BinSettings Im_omega_bins(100, 0, 1800);
    const BinSettings Im_pi0_bins(100, 0, 1000);
    const BinSettings wpi0_Q_bins(100, 0, 1500);

    const BinSettings E_photon_bins(100, 0, 1600);
    const BinSettings E_proton_bins(100, 900, 2000);
    const BinSettings E_omega_bins(100, 400, 1800);
    const BinSettings E_pi0_bins(100, 0, 1600);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    auto hfTagger = new HistogramFactory("Tagger", HistFac, "");
    auto hfOnlyEnergy = new HistogramFactory("All_the_candidate_energies", HistFac, "");
    auto hfOnlyAngles = new HistogramFactory("All_the_candidate_angles", HistFac, "");
    auto hfCrossCheck = new HistogramFactory("Cross_section_check", HistFac, "");
    auto hf_wpi0g_EvsTheta = new HistogramFactory("Energy_vs_theta_dists", HistFac, "");
    auto hf_wpi0g_IM = new HistogramFactory("Invariant_mass_dists", HistFac, "");
    //auto hf_extras = new HistogramFactory("Some_additional_extra_hists", HistFac, "");

    for (unsigned int i=0; i<nrCuts; i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        string myString2[nrCuts] = {"CUT#1_Sel3Neu1Cha", "CUT#2_OmegaEthreshold", "CUT#3_ImMissingParticle_+-20%mp", "CUT#4_SelMinM(2neu-mpi0)_+-20%mpi0"};

        h_missingProton_Im[i] = hf_wpi0g_IM->makeTH1D("Im(missingP) after"+myString2[i],     // title
                                                             "IM(missing_particle) [MeV]", "#",     // xlabel, ylabel
                                                             Im_proton_bins,  // our binnings
                                                             "h_missingProton_Im"+myString2[i], false     // ROOT object name, auto-generated if omitted
                                                             );

        h_wOnly3g_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. Omega invariant mass",     // title
                                            "IM(3g) [MeV]", "#",     // xlabel, ylabel
                                            Im_omega_bins,  // our binnings
                                            "h_wOnly3g_Im"+myString, false     // ROOT object name, auto-generated if omitted
                                            );

        h_2gMassComb[i] = hf_wpi0g_IM->makeTH1D("3-Body Omega decay mass-combinations",     // title
                                             "IM(2g comb) [MeV]", "#",     // xlabel, ylabel
                                             Im_pi0_bins,  // our binnings
                                             "h_2gMassComb"+myString, false     // ROOT object name, auto-generated if omitted
                                             );

        h_ChaCaloVSVetoEnergies_TAPS[i] = hfTagger->makeTH2D("Charged candidates (EVeto > 0): Deposited Calo VS Veto energy in TAPS",     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_ChaCaloVSVetoEnergies_TAPS"+myString, true    // ROOT object name, auto-generated if omitted
                                         );

        h_ChaCaloVSVetoEnergies_CB[i] = hfTagger->makeTH2D("Charged candidates (EVeto > 0): Calo VS Veto energy in CB",     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_ChaCaloVSVetoEnergies_CB"+myString, true    // ROOT object name, auto-generated if omitted
                                         );
    }

    for (unsigned int i=0; i<nrCuts_IM_pi0; i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        h_Pi0Only2g_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. Pi0 invariant mass",     // title
                                         "IM(pi0) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_Pi0Only2g_Im"+myString, false     // ROOT object name, auto-generated if omitted
                                         );

    }


    h_AllCaloVSVetoEnergies_TAPS = hfTagger->makeTH2D("All candidates: Deposited Calo VS Veto energy in TAPS",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_AllCaloVSVetoEnergies_TAPS", true    // ROOT object name, auto-generated if omitted
                                     );

    h_AllCaloVSVetoEnergies_CB = hfTagger->makeTH2D("All candidates: Deposited Calo VS Veto energy in CB",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_AllCaloVSVetoEnergies_CB", true    // ROOT object name, auto-generated if omitted
                                     );


    h_NeuCaloVSVetoEnergies_TAPS = hfTagger->makeTH2D("Uncharged candidates: Deposited Calo VS Veto energy in TAPS",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_NeuCaloVSVetoEnergies_TAPS", true    // ROOT object name, auto-generated if omitted
                                     );

    h_NeuCaloVSVetoEnergies_CB = hfTagger->makeTH2D("Uncharged candidates: Deposited Calo VS Veto energy in CB",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_NeuCaloVSVetoEnergies_CB", true    // ROOT object name, auto-generated if omitted
                                     );

    h_TaggerTime = hfTagger->makeTH1D("Tagger Time",     // title
                                    "t [ns]", "#",     // xlabel, ylabel
                                    tagger_time_bins,  // our binnings
                                    "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                    );
    
    h_InitialBeamE = hfOnlyEnergy->makeTH1D("Photon beam - Energy distribution",     // title
                                   "E [GeV]", "#",     // xlabel, ylabel
                                   initialBeamE_bins,  // our binnings
                                   "h_InitialBeamE"     // ROOT object name, auto-generated if omitted
                                   );

    h_nClusters = hfTagger->makeTH1D("Number of Clusters", "nClusters", "#", BinSettings(15), "h_nClusters");

    h_VetoEnergies = hfTagger->makeTH1D("Veto Energies", "E [MeV]", "#", Veto_Energy_bins, "h_VetoEnergies");


    h_NeuPolarAngles = hfOnlyAngles->makeTH1D("Neutral candidate polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_NeuPolarAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_NeuAzimuthAngles = hfOnlyAngles->makeTH1D("Neutral candidate azimuthangles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_NeuAzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_ChaPolarAngles = hfOnlyAngles->makeTH1D("Charged particles polarangles",     // title
                                                    "#Theta [deg]", "#",     // xlabel, ylabel
                                                    theta_bins,  // our binnings
                                                    "h_ChaPolarAngles", false     // ROOT object name, auto-generated if omitted
                                                    );

    h_ChaAzimuthAngles = hfOnlyAngles->makeTH1D("Charged particles azimuthangles",     // title
                                                    "#Phi [deg]", "#",     // xlabel, ylabel
                                                    phi_bins,  // our binnings
                                                    "h_ChaAzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                                    );

    h_3gPolarAngles = hfOnlyAngles->makeTH1D("Rec. Photons polarangles",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_3gPolarAngles", false     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesCB = hfOnlyAngles->makeTH1D("Rec. Photons polarangles only in CB",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins_CB,  // our binnings
                                             "h_3gPolarAnglesCB", false     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesTAPS = hfOnlyAngles->makeTH1D("Rec. Photons polarangles only in TAPS",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins_TAPS,  // our binnings
                                             "h_3gPolarAnglesTAPS", false     // ROOT object name, auto-generated if omitted
                                             );

    h_3gAzimuthAngles = hfOnlyAngles->makeTH1D("Rec. Photons azimuthangles",     // title
                                               "#Phi [deg]", "#",     // xlabel, ylabel
                                               phi_bins,  // our binnings
                                               "h_3gAzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                               );

    h_3g_EvTheta_CB = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in CB", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_CB, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_CB", true    // ROOT object name, auto-generated if omitted
                                       );

    h_3g_EvTheta_TAPS = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in TAPS", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_TAPS, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_TAPS", true    // ROOT object name, auto-generated if omitted
                                       );

    h_p_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. proton energy vs theta", //title
                                                 "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                 theta_bins, E_proton_bins,    // our binnings
                                                 "h_p_EvTheta", true    // ROOT object name, auto-generated if omitted
                                                 );

    h_w_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. omega meson energy vs theta", //title
                                                 "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                 theta_bins, E_omega_bins,    // our binnings
                                                 "h_w_EvTheta", true    // ROOT object name, auto-generated if omitted
                                                 );

    h_wg_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. omega decay-photon energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_photon_bins,    // our binnings
                                       "h_wg_EvTheta", true    // ROOT object name, auto-generated if omitted
                                       );

    h_wpi0_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. Pi0 energy vs theta", //title
                                        "Theta [deg]","E [MeV]", // xlabel, ylabel
                                        theta_bins, E_pi0_bins,    // our binnings
                                        "h_wpi0_EvTheta", true    // ROOT object name, auto-generated if omitted
                                        );

    h_wpi02g_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. Pi0 decay-photons energy vs theta", //title
                                        "Theta [deg]","E [MeV]", // xlabel, ylabel
                                        theta_bins, E_pi0_bins,    // our binnings
                                        "h_wpi02g_EvTheta", true    // ROOT object name, auto-generated if omitted
                                        );

    h_doubly_wp_DCS_reconstructed_lab = hfCrossCheck->makeTH2D("Cross section check in lab-frame", //title
                                                                        "cos_Theta","sqrt_S [MeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_wp_DCS_reconstructed_lab", true    // ROOT object name, auto-generated if omitted
                                                                        );
    
    h_doubly_wp_DCS_reconstructed_cmFrame = hfCrossCheck->makeTH2D("Cross section check in cm-frame", //title
                                                                        "cos_Theta*","sqrt_S [MeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_wp_DCS_reconstructed_cmFrame", true    // ROOT object name, auto-generated if omitted
                                                                        );

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("tree"));

}

void scratch_damaurer_MC_pw_ppi0g_p3g::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event); 

    //Fill e.g. the polar angle distribution into a histogram
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList all;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    vector<double> allCanCluSize;
    vector<double> allCanCaloE;
    vector<double> allCanVetoE;
    vector<double> allThe;
    vector<double> allPhi;
    vector<double> neuCanCluSize;
    vector<double> neuCanCaloE;
    vector<double> neuCanVetoE;
    vector<double> neuThe;
    vector<double> neuPhi;
    vector<double> chaCanCluSize;
    vector<double> chaCanCaloE;
    vector<double> chaCanVetoE;
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
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
            h_AllCaloVSVetoEnergies_CB->Fill(cand->CaloEnergy,cand->VetoEnergy);
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
            h_AllCaloVSVetoEnergies_TAPS->Fill(cand->CaloEnergy,cand->VetoEnergy);
        if(cand->VetoEnergy <= vetoEthreshold) {
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                h_NeuCaloVSVetoEnergies_CB->Fill(cand->CaloEnergy,cand->VetoEnergy);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                h_NeuCaloVSVetoEnergies_TAPS->Fill(cand->CaloEnergy,cand->VetoEnergy);
            neutral.emplace_back(cand);
            neuThe.push_back(cand->Theta*radtodeg);
            neuPhi.push_back(cand->Phi*radtodeg);
            neuCanCluSize.push_back(cand->ClusterSize);
            neuCanCaloE.push_back(cand->CaloEnergy);
            neuCanVetoE.push_back(cand->VetoEnergy);
        }       
        else{
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                h_ChaCaloVSVetoEnergies_CB[0]->Fill(cand->CaloEnergy,cand->VetoEnergy);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                h_ChaCaloVSVetoEnergies_TAPS[0]->Fill(cand->CaloEnergy,cand->VetoEnergy);
            charged.emplace_back(cand);
            chaThe.push_back(cand->Theta*radtodeg);
            chaPhi.push_back(cand->Phi*radtodeg);
            chaCanCluSize.push_back(cand->ClusterSize);
            chaCanCaloE.push_back(cand->CaloEnergy);
            chaCanVetoE.push_back(cand->VetoEnergy);
        }

    }

    //getting access to the pi0 decay photons

    double wpi0gOnly3g_IM=-100;
    TParticle g[neu_nrSel];
    TParticle proton_tmp;

    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;
    vector<double> wpi0g_times;

    if(neuCanCaloE.size() == neu_nrSel && chaCanCaloE.size() == cha_nrSel){

        for(int i = 0; i<neu_nrSel; i++){
            g[i] = TParticle(ParticleTypeDatabase::Photon, neutral[i]);
        }

        proton_tmp = TParticle(ParticleTypeDatabase::Proton, charged[0]);

        TParticle wpi0gOnly3g;
        double wpi0gOnly3g_time = 0;
        for (const auto& photon : neutral) {
            wpi0gOnly3g += TParticle(ParticleTypeDatabase::Photon, photon);
            wpi0gOnly3g_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                PhotonTheCB.push_back(photon->Theta*radtodeg);
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                PhotonTheTAPS.push_back(photon->Theta*radtodeg);
            PhotonThe.push_back(photon->Theta*radtodeg);
            PhotonPhi.push_back(photon->Phi*radtodeg);
            PhotonE.push_back(photon->CaloEnergy);
        }
        wpi0gOnly3g_IM = wpi0gOnly3g.M();
        wpi0gOnly3g_time = wpi0gOnly3g_time/2.;
        wpi0g_times.push_back(wpi0gOnly3g_time);

    }

    //Looping over the taggerhits

    TParticle proton;
    TParticle omega_tmp;

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) { // Event loop

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);

        //Adding selections & filling histograms

        const double weight = promptrandom.FillWeight();

        h_TaggerTime->Fill(taggerhit.Time, weight);
        h_InitialBeamE->Fill(InitialPhotonVec.E()/1000,weight);

        if(!(neuCanCaloE.size() == neu_nrSel && chaCanCaloE.size() == cha_nrSel))
            continue;

        TLorentzVector L3g;
        for(int i = 0; i<neu_nrSel; i++){
            L3g += (TLorentzVector)(g[i]);
        }

        omega_tmp = TParticle(ParticleTypeDatabase::Omega,(TLorentzVector)(L3g));

        TLorentzVector Lw_tmp = (TLorentzVector)(omega_tmp);
        TLorentzVector Lp_tmp = (TLorentzVector)(proton_tmp);

        TLorentzVector LmissingProton = Linitial-Lw_tmp;     

        h_wOnly3g_Im[0]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[0]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[0]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[0]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[0]->Fill((g[1]+g[2]).M(),weight);

        if(!(InitialPhotonVec.E()>=Omega_Ethreshold))
            continue;

        for (unsigned int i=0; i<chaThe.size(); i++){
            if(chaThe[i]>=20 && chaThe[i]<=160)
                h_ChaCaloVSVetoEnergies_CB[1]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
            if(chaThe[i]<20 && chaThe[i]>=2)
                h_ChaCaloVSVetoEnergies_TAPS[1]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
        }

        h_wOnly3g_Im[1]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[1]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[1]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[1]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[1]->Fill((g[1]+g[2]).M(),weight);

        for (unsigned int i=0; i<neuThe.size(); i++){
            h_NeuPolarAngles->Fill(neuThe[i],weight);
        }

        for (unsigned int i=0; i<neuPhi.size(); i++){
            h_NeuAzimuthAngles->Fill(neuPhi[i],weight);
        }

        double cut_shift1 = mp*0.2;
        //t.ProtonIm = proton_missing_particle.M();
        //double cut_shift2 = mw*0.2;

        if(!(LmissingProton.M() > (mp-cut_shift1) && LmissingProton.M() < (mp+cut_shift1)))
            continue;

        h_wOnly3g_Im[2]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[2]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[2]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[2]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[2]->Fill((g[1]+g[2]).M(),weight);

        TLorentzVector Lp = Lp_tmp;
        TLorentzVector Lw = Lw_tmp;
        TLorentzVector Lw_boosted = Lw;
        Lw_boosted.Boost(-Linitial.BoostVector());

        h_doubly_wp_DCS_reconstructed_lab->Fill(cos(Lw.Theta()),Linitial.M(),weight);
        h_doubly_wp_DCS_reconstructed_cmFrame->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        //LorentzBoost-Stuff with opening angle selections:

        /*
        TParticle wpi0_tmp1 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1+g2));
        TParticle wpi0_tmp2 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1+g3));
        TParticle wpi0_tmp3 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g2+g3));

        TParticle wg_tmp1 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g3));
        TParticle wg_tmp2 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g2));
        TParticle wg_tmp3 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1));

        TLorentzVector wpi0_tmp1_boosted = (TLorentzVector)wpi0_tmp1;
        TLorentzVector wpi0_tmp2_boosted = (TLorentzVector)wpi0_tmp2;
        TLorentzVector wpi0_tmp3_boosted = (TLorentzVector)wpi0_tmp3;

        TLorentzVector wg_tmp1_boosted = (TLorentzVector)wg_tmp1;
        TLorentzVector wg_tmp2_boosted = (TLorentzVector)wg_tmp2;
        TLorentzVector wg_tmp3_boosted = (TLorentzVector)wg_tmp3;

        wpi0_tmp1_boosted.Boost(-Lw.BoostVector());
        wpi0_tmp2_boosted.Boost(-Lw.BoostVector());
        wpi0_tmp3_boosted.Boost(-Lw.BoostVector());

        wg_tmp1_boosted.Boost(-Lw.BoostVector());
        wg_tmp2_boosted.Boost(-Lw.BoostVector());
        wg_tmp3_boosted.Boost(-Lw.BoostVector());

        double cmAngleComb1 = cos(wpi0_tmp1_boosted.Angle(wg_tmp1_boosted.Vect()));
        double cmAngleComb2 = cos(wpi0_tmp2_boosted.Angle(wg_tmp2_boosted.Vect()));
        double cmAngleComb3 = cos(wpi0_tmp3_boosted.Angle(wg_tmp3_boosted.Vect()));

        double array[] = {1.2,1.1,1.4,1.5,0.9};
        double min_i;
        double min = 10;
        int arraySize = sizeof(array)/sizeof(array[0]);

        for (int i=0; i<arraySize; i++){
            //array[i] = i;
            if(array[i] <= min){min = array[i]; min_i = i;}
            cout << "Arrayelement = " << array[i] << endl;
        }

        */

        TLorentzVector w_decayComb[neu_nrSel] = {(g[0] + g[1]),(g[0] + g[2]),(g[1] + g[2])};
        double inv_mass_comb[neu_nrSel] = {sqrt(Lw.M2()-w_decayComb[1].M2()-w_decayComb[2].M2()),sqrt(Lw.M2()-w_decayComb[0].M2()-w_decayComb[2].M2()),sqrt(Lw.M2()-w_decayComb[0].M2()-w_decayComb[1].M2())};
        double pi0_absMassDiff[neu_nrSel];

        TLorentzVector wpi0g[nr_pi0g];
        TLorentzVector wpi0;
        TLorentzVector wg;

        TParticle wpi0_tmp[neu_nrSel];
        TParticle wg_tmp[neu_nrSel];

        for(int i=0; i<neu_nrSel;i++){
            wpi0_tmp[i] = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(w_decayComb[i]));
            wg_tmp[i] = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g[neu_nrSel-1-i]));
            pi0_absMassDiff[i] = abs(inv_mass_comb[i]-mpi0);
        }

        if(pi0_absMassDiff[0] <= pi0_absMassDiff[1] && pi0_absMassDiff[0] <= pi0_absMassDiff[2]){wpi0 = wpi0_tmp[0]; wg = wg_tmp[0]; wpi0g[0] = g[0]; wpi0g[1] = g[1];}
        else if(pi0_absMassDiff[1] <= pi0_absMassDiff[0] && pi0_absMassDiff[1] <= pi0_absMassDiff[2]){wpi0 = wpi0_tmp[1]; wg = wg_tmp[1]; wpi0g[0] = g[0]; wpi0g[1] = g[2];}
        else{wpi0 = wpi0_tmp[2]; wg = wg_tmp[2]; wpi0g[0] = g[1]; wpi0g[1] = g[2];}

        h_Pi0Only2g_Im[0]->Fill(wpi0.M(),weight);

        double cut_shift3 = mpi0*0.2;
        if(!(wpi0.M()>(mpi0-cut_shift3) && wpi0.M()<(mpi0+cut_shift3)))
            continue;

        h_Pi0Only2g_Im[1]->Fill(wpi0.M(),weight);

        h_wOnly3g_Im[3]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[3]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[3]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[3]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[3]->Fill((g[1]+g[2]).M(),weight);

        for (unsigned int i=0; i<chaThe.size(); i++){
            if(chaThe[i]>=20 && chaThe[i]<=160)
                h_ChaCaloVSVetoEnergies_CB[2]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
            if(chaThe[i]<20 && chaThe[i]>=2)
                h_ChaCaloVSVetoEnergies_TAPS[2]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
        }

        /*
        if(!(charged[0]->VetoEnergy>0.9))
                continue;

        if(!(proton_tmp.Theta()*radtodeg <= 80))
            continue;

        for (unsigned int i=0; i<chaThe.size(); i++){
            if(chaThe[i]>=20 && chaThe[i]<=160)
                h_ChaCaloVSVetoEnergies_CB[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
            if(chaThe[i]<20 && chaThe[i]>=2)
                h_ChaCaloVSVetoEnergies_TAPS[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
        }

        proton = TParticle(ParticleTypeDatabase::Proton,(TLorentzVector)(proton_tmp));

        h_wOnly3g_Im[4]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[4]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[4]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[4]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[4]->Fill((g[1]+g[2]).M(),weight);

        */

        h_p_EvTheta->Fill(Lp.Theta()*radtodeg,Lp.E(),weight);
        h_w_EvTheta->Fill(Lw.Theta()*radtodeg,Lw.E(),weight);
        h_wg_EvTheta->Fill(wg.Theta()*radtodeg,wg.E(),weight);
        h_wpi0_EvTheta->Fill(wpi0.Theta()*radtodeg,wpi0.E(),weight);

        for(int i=0; i<nr_pi0g;i++){
            h_wpi02g_EvTheta->Fill(wpi0g[i].Theta()*radtodeg,wpi0g[i].E(),weight);
        }

        //weight_res += weight;

        //Filling the tree for further analysis

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = event.Reconstructed().Clusters.size();
        t.PhotonAzimuthAngles = PhotonPhi;
        t.PhotonPolarAngles = PhotonThe;
        t.PhotonPolarAnglesCB = PhotonTheCB;
        t.PhotonPolarAnglesTAPS = PhotonTheTAPS;
        t.Tree->Fill();     
    }
    h_nClusters->Fill(event.Reconstructed().Clusters.size());

}

void scratch_damaurer_MC_pw_ppi0g_p3g::ShowResult()
{    

   ant::canvas(GetName()+": All candidates deposited Calo vs Veto energy")
           << drawoption("pcolz")
           << h_AllCaloVSVetoEnergies_CB
           << h_AllCaloVSVetoEnergies_TAPS
           << endc; // actually draws the canvas

   ant::canvas(GetName()+": Deposited Calo vs Veto energy of expected uncharged candidates.")
           << drawoption("pcolz")
           << h_NeuCaloVSVetoEnergies_CB
           << h_NeuCaloVSVetoEnergies_TAPS
           << endc; // actually draws the canvas

   for (unsigned int i=0; i<nrCuts_IM; i++){

       stringstream ss;
       ss << i;
       string myString = ss.str();
       ant::canvas(GetName()+": Deposited Calo vs Veto energy of expected charged candidates after "+myString+ " applied cuts")
               << drawoption("pcolz")
               << h_ChaCaloVSVetoEnergies_CB[i]
               << h_ChaCaloVSVetoEnergies_TAPS[i]
               << endc; // actually draws the canvas
   }


    for (unsigned int i=0; i<nrCuts_IM; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. Omega & missing particle: Invariant masses after "+myString+ " applied cuts")
                << h_missingProton_Im[i]
                << h_wOnly3g_Im[i]
                << h_2gMassComb[i]
                << endc; // actually draws the canvas

    }

    for (unsigned int i=0; i<nrCuts_IM_pi0; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. pi0: Invariant mass after "+myString+ " applied cuts")
                << h_Pi0Only2g_Im[i]
                << endc; // actually draws the canvas

    }

    ant::canvas(GetName()+": Rec. Omega & protons: EvsTheta distribution")
            << drawoption("pcolz")
            << h_w_EvTheta
            << h_p_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Rec. wpi0g decay: EvsTheta distribution")
            << drawoption("pcolz")
            << h_wg_EvTheta
            << h_wpi0_EvTheta
            << h_wpi02g_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Cross section checks")
            << drawoption("Surf")
            << h_doubly_wp_DCS_reconstructed_lab
            << h_doubly_wp_DCS_reconstructed_cmFrame
            << endc; // actually draws the canvas

}

void scratch_damaurer_MC_pw_ppi0g_p3g::Finish()
{
    cout << "please work!" << endl;
/*
    int max_particles_detected = h_ALL_PolarAngles->GetEntries();
    int Photons_max_detected = h_3gPolarAngles->GetEntries();
    int MissingCHA_max_detected = h_missingChaPolarAngles->GetEntries();

    PhotonPolarAngles_TAPS_frac = (h_3gPolarAnglesTAPS->GetEntries())/Photons_max_detected;
    PhotonPolarAngles_CB_frac = (h_3gPolarAnglesCB->GetEntries())/Photons_max_detected;
    //ProtonPolarAngles_TAPS_frac = ProtonPolarAngles_TAPS/MissingCHA_max_detected;
    //ProtonPolarAngles_CB_frac = ProtonPolarAngles_CB/MissingCHA_max_detected;

    PhotonPolarAngles_TAPS_err = sqrt((PhotonPolarAngles_TAPS_frac*(1-PhotonPolarAngles_TAPS_frac))/Photons_max_detected);
    PhotonPolarAngles_CB_err = sqrt((PhotonPolarAngles_CB_frac*(1-PhotonPolarAngles_CB_frac))/Photons_max_detected);
    //ProtonPolarAngles_TAPS_err = sqrt((ProtonPolarAngles_TAPS_frac*(1-ProtonPolarAngles_TAPS_frac))/MissingCHA_max_detected);
    //ProtonPolarAngles_CB_err = sqrt((ProtonPolarAngles_CB_frac*(1-ProtonPolarAngles_CB_frac))/MissingCHA_max_detected);

    cout << "" << endl;
    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;
    cout << "Detection efficiency: " << max_particles_detected/max_particles << endl;
    cout << "Max number of detected photons: " << Photons_max_detected << endl;
*/
    /*
    int max_particles_detected = h_PolarAngles->GetEntries();
    int Photons_max_detected = h_PhotonPolarAngles->GetEntries();
    int Protons_max_detected = h_ProtonPolarAngles->GetEntries();

    PhotonPolarAngles_TAPS_frac = (h_PhotonPolarAngles_TAPS->GetEntries())/Photons_max_detected;
    PhotonPolarAngles_CB_frac = (h_PhotonPolarAngles_CB->GetEntries())/Photons_max_detected;
    ProtonPolarAngles_TAPS_frac = ProtonPolarAngles_TAPS/Protons_max_detected;
    ProtonPolarAngles_CB_frac = ProtonPolarAngles_CB/Protons_max_detected;

    PhotonPolarAngles_TAPS_err = sqrt((PhotonPolarAngles_TAPS_frac*(1-PhotonPolarAngles_TAPS_frac))/Photons_max_detected);
    PhotonPolarAngles_CB_err = sqrt((PhotonPolarAngles_CB_frac*(1-PhotonPolarAngles_CB_frac))/Photons_max_detected);
    ProtonPolarAngles_TAPS_err = sqrt((ProtonPolarAngles_TAPS_frac*(1-ProtonPolarAngles_TAPS_frac))/Protons_max_detected);
    ProtonPolarAngles_CB_err = sqrt((ProtonPolarAngles_CB_frac*(1-ProtonPolarAngles_CB_frac))/Protons_max_detected);

    cout << "" << endl;
    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;
    cout << "Detection efficiency: " << max_particles_detected/max_particles << endl;
    cout << "Max number of detected photons: " << Photons_max_detected << endl;
    cout << "Max number of detected protons: " << Protons_max_detected << endl;
    cout << "Photon amount falling into TAPS: " << 100*PhotonPolarAngles_TAPS_frac << "% with error: " << 100*PhotonPolarAngles_TAPS_err << "%" << endl;
    cout << "Photon amount falling into CB: " << 100*PhotonPolarAngles_CB_frac << "% with error: " << 100*PhotonPolarAngles_CB_err << "%" << endl;
    cout << "Proton amount falling into TAPS: " << 100*ProtonPolarAngles_TAPS_frac << "% with error: " << 100*ProtonPolarAngles_TAPS_err << "%" << endl;
    cout << "Proton amount falling into CB: " << 100*ProtonPolarAngles_CB_frac << "% with error: " << 100*ProtonPolarAngles_CB_err << "%" << endl;
    cout << "" << endl;

    cout << "1) Both CB: " << (float)100*((float)count_first/(float)weight_res) << "%" << endl;
    cout << "2) CB & BaF2: " << (float)100*((float)count_second/(float)weight_res) << "%" << endl;
    cout << "3) CB & PbWO: " << (float)100*((float)count_third/(float)weight_res) << "%" << endl;
    cout << "4) Both BaF2: " << (float)100*((float)count_fourth/(float)weight_res) << "%" << endl;
    cout << "5) BaF2 & PbWO: " << (float)100*((float)count_fifth/(float)weight_res) << "%" << endl;
    cout << "6) Both PbWO: " << (float)100*((float)count_sixth/(float)weight_res) << "%" << endl;
    //cout << "7) PbWO & not detected: " << (float)100*((float)count_seventh/(float)weight_res) << "%" << endl;
    //cout << "8) BaF2 & not detected: " << (float)100*((float)count_eigth/(float)weight_res) << "%" << endl;
    //cout << "9) CB & not detected: " << (float)100*((float)count_nineth/(float)weight_res) << "%" << endl;
    //cout << "10) Both not detected: " << (float)100*((float)count_tenth/(float)weight_res) << "%" << endl;
    */

}

AUTO_REGISTER_PHYSICS(scratch_damaurer_MC_pw_ppi0g_p3g)
