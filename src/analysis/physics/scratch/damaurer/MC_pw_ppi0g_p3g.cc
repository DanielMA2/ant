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
    const BinSettings Im_omega_bins(100, 0, 1600);
    const BinSettings Im_pi0_bins(100, 0, 1000);

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

    for (unsigned int i=0; i<nrCuts_IM; i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        h_proton_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. proton invariant mass",     // title
                                                             "IM(proton) [MeV]", "#",     // xlabel, ylabel
                                                             Im_proton_bins,  // our binnings
                                                             "h_proton_Im"+myString, false     // ROOT object name, auto-generated if omitted
                                                             );
        h_wOnly3g_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. Omega invariant mass",     // title
                                            "IM(w) [MeV]", "#",     // xlabel, ylabel
                                            Im_omega_bins,  // our binnings
                                            "h_wOnly3g_Im"+myString, false     // ROOT object name, auto-generated if omitted
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

    h_ChaCaloVSVetoEnergies_TAPS = hfTagger->makeTH2D("Charged candidates: Deposited Calo VS Veto energy in TAPS",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_ChaCaloVSVetoEnergies_TAPS", true    // ROOT object name, auto-generated if omitted
                                     );

    h_ChaCaloVSVetoEnergies_CB = hfTagger->makeTH2D("Charged candidates: Calo VS Veto energy in CB",     // title
                                     "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                     Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                     "h_ChaCaloVSVetoEnergies_CB", true    // ROOT object name, auto-generated if omitted
                                     );



    /*
    for (unsigned int i=0; i<nrCuts_BackToBack; i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        h_pi0g_BackToBack[i] = hfOnlyAngles->makeTH1D("w_pi0g angles in cm of rec. omega meson",     // title
                                                     "cos_Theta*", "#",     // xlabel, ylabel
                                                     cos_bins_BackToBack,  // our binnings
                                                     "h_pi0g_BackToBack"+myString, false     // ROOT object name, auto-generated if omitted
                                                     );

    }
    */


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

    h_2gMassComb = hf_wpi0g_IM->makeTH1D("3-Body Omega decay mass-combinations",     // title
                                         "IM(2g) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2gMassComb", false     // ROOT object name, auto-generated if omitted
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
                                                                        "cos_Theta","sqrt_S [GeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_wp_DCS_reconstructed_lab", true    // ROOT object name, auto-generated if omitted
                                                                        );
    
    h_doubly_wp_DCS_reconstructed_cmFrame = hfCrossCheck->makeTH2D("Cross section check in cm-frame", //title
                                                                        "cos_Theta*","sqrt_S [GeV]", // xlabel, ylabel
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
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    vector<double> NeuCanCluSize;
    vector<double> NeuCanCaloE;
    vector<double> ChaCanCluSize;
    vector<double> ChaCanCaloE;
    vector<double> neuThe;
    vector<double> neuPhi;
    vector<double> chaThe;
    vector<double> chaPhi;
    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;
    vector<double> wpi0g_times;

    for(const auto& cand : candidates.get_iter()) {
        h_VetoEnergies->Fill(cand->VetoEnergy);
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
            h_AllCaloVSVetoEnergies_CB->Fill(cand->CaloEnergy,cand->VetoEnergy);
        if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
            h_AllCaloVSVetoEnergies_TAPS->Fill(cand->CaloEnergy,cand->VetoEnergy);
        neuThe.push_back(cand->Theta*radtodeg);
        neuPhi.push_back(cand->Phi*radtodeg);
        if(cand->VetoEnergy < vetoEthreshold) {
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                h_NeuCaloVSVetoEnergies_CB->Fill(cand->CaloEnergy,cand->VetoEnergy);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                h_NeuCaloVSVetoEnergies_TAPS->Fill(cand->CaloEnergy,cand->VetoEnergy);
            neutral.emplace_back(cand);
            NeuCanCluSize.push_back(cand->ClusterSize);
            NeuCanCaloE.push_back(cand->CaloEnergy);
        }       
        else{
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                h_ChaCaloVSVetoEnergies_CB->Fill(cand->CaloEnergy,cand->VetoEnergy);
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                h_ChaCaloVSVetoEnergies_TAPS->Fill(cand->CaloEnergy,cand->VetoEnergy);
            charged.emplace_back(cand);
            ChaCanCluSize.push_back(cand->ClusterSize);
            ChaCanCaloE.push_back(cand->CaloEnergy);
            chaThe.push_back(cand->Theta*radtodeg);
            chaPhi.push_back(cand->Phi*radtodeg);
        }

    }

    //getting access to the pi0 decay photons

    double wpi0gOnly3g_IM=-100;
    TParticle g1;
    TParticle g2;
    TParticle g3;
    TParticle proton;

    if(NeuCanCaloE.size() == neu_nrSel && ChaCanCaloE.size() == cha_nrSel){

        g1 = TParticle(ParticleTypeDatabase::Photon, neutral[0]);
        g2 = TParticle(ParticleTypeDatabase::Photon, neutral[1]);
        g3 = TParticle(ParticleTypeDatabase::Photon, neutral[2]);
        proton = TParticle(ParticleTypeDatabase::Proton, charged[0]);

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

    TParticle proton_tmp;
    TParticle omega_tmp;

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) { // Event loop

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);

        omega_tmp = TParticle(ParticleTypeDatabase::Omega,(TLorentzVector)(g1+g2+g3));
        proton_tmp = TParticle(ParticleTypeDatabase::Proton,(TLorentzVector)(proton));

        TLorentzVector Lw_tmp = (TLorentzVector)(omega_tmp);
        TLorentzVector Lp_tmp = (TLorentzVector)(proton_tmp);

        //Adding selections & filling histograms

        const double weight = promptrandom.FillWeight();

        h_TaggerTime->Fill(taggerhit.Time, weight);
        h_InitialBeamE->Fill(InitialPhotonVec.E()/1000,weight);

        if(NeuCanCaloE.size() == neu_nrSel && ChaCanCaloE.size() == cha_nrSel ){

            h_wOnly3g_Im[0]->Fill(omega_tmp.M(),weight);
            h_proton_Im[0]->Fill(proton_tmp.M(),weight);

            for (unsigned int i=0; i<neuThe.size(); i++){
                h_NeuPolarAngles->Fill(neuThe[i],weight);
            }

            for (unsigned int i=0; i<neuPhi.size(); i++){
                h_NeuAzimuthAngles->Fill(neuPhi[i],weight);
            }

            //t.ProtonIm = proton_missing_particle.M();

            double cut_shift1 = mp*0.2;
            double cut_shift2 = mw*0.2;
            if(proton_tmp.M()>(mp-cut_shift1) && proton_tmp.M()<(mp+cut_shift1) && omega_tmp.M()>(mw-cut_shift2) && omega_tmp.M()<(mw+cut_shift2)){

                h_proton_Im[1]->Fill(proton_tmp.M(),weight);
                h_wOnly3g_Im[1]->Fill(omega_tmp.M(),weight);

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

                TLorentzVector w_decayComb1 = (TLorentzVector)(g1 + g2);
                TLorentzVector w_decayComb2 = (TLorentzVector)(g1 + g3);
                TLorentzVector w_decayComb3 = (TLorentzVector)(g2 + g3);

                TParticle wpi0_tmp1 = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(g1+g2));
                TParticle wpi0_tmp2 = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(g1+g3));
                TParticle wpi0_tmp3 = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(g2+g3));

                h_Pi0Only2g_Im[0]->Fill((wpi0_tmp1).M(),weight);
                h_Pi0Only2g_Im[0]->Fill((wpi0_tmp2).M(),weight);
                h_Pi0Only2g_Im[0]->Fill((wpi0_tmp3).M(),weight);

                TParticle wg_tmp1 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g3));
                TParticle wg_tmp2 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g2));
                TParticle wg_tmp3 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1));

                TLorentzVector wpi0;
                TLorentzVector wpi0g1;
                TLorentzVector wpi0g2;

                TLorentzVector wg;

                double inv_mass_comb1 = sqrt(Lw.M2()-w_decayComb2.M2()-w_decayComb3.M2());
                double inv_mass_comb2 = sqrt(Lw.M2()-w_decayComb1.M2()-w_decayComb3.M2());
                double inv_mass_comb3 = sqrt(Lw.M2()-w_decayComb1.M2()-w_decayComb2.M2());

                h_2gMassComb->Fill(inv_mass_comb1,weight);
                h_2gMassComb->Fill(inv_mass_comb2,weight);
                h_2gMassComb->Fill(inv_mass_comb3,weight);

                double pi0_absMassDiff1 = abs(inv_mass_comb1-mpi0);
                double pi0_absMassDiff2 = abs(inv_mass_comb2-mpi0);
                double pi0_absMassDiff3 = abs(inv_mass_comb3-mpi0);

                if(pi0_absMassDiff1 <= pi0_absMassDiff2 && pi0_absMassDiff1 <= pi0_absMassDiff3){wpi0 = wpi0_tmp1; wg = wg_tmp1; wpi0g1 = g1; wpi0g2 = g2;}
                else if(pi0_absMassDiff2 <= pi0_absMassDiff1 && pi0_absMassDiff2 <= pi0_absMassDiff3){wpi0 = wpi0_tmp2; wg = wg_tmp2; wpi0g1 = g1; wpi0g2 = g3;}
                else{wpi0 = wpi0_tmp3; wg = wg_tmp3; wpi0g1 = g2; wpi0g2 = g3;}

                double cut_shift3 = mpi0*0.2;

                h_Pi0Only2g_Im[1]->Fill(wpi0.M(),weight);

                if(wpi0.M()>(mpi0-cut_shift3) && wpi0.M()<(mpi0+cut_shift3)){
                    h_Pi0Only2g_Im[2]->Fill(wpi0.M(),weight);

                    h_p_EvTheta->Fill(Lp.Theta()*radtodeg,Lp.E(),weight);
                    h_w_EvTheta->Fill(Lw.Theta()*radtodeg,Lw.E(),weight);
                    h_wg_EvTheta->Fill(wg.Theta()*radtodeg,wg.E(),weight);
                    h_wpi0_EvTheta->Fill(wpi0.Theta()*radtodeg,wpi0.E(),weight);
                    h_wpi02g_EvTheta->Fill(wpi0g1.Theta()*radtodeg,wpi0g1.E(),weight);
                    h_wpi02g_EvTheta->Fill(wpi0g2.Theta()*radtodeg,wpi0g2.E(),weight);
                }

                weight_res += weight;
            }
        }

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
/*
    ant::canvas(GetName()+": tagger stuff")
            << h_TaggerTime
            << h_nClusters
            << h_VetoEnergies
            << h_InitialBeamE
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Polarangles in whole setup")
            << h_ALL_PolarAngles
            << h_3gPolarAngles
            << h_missingChaPolarAngles
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Azimuthangles in whole setup")
            << h_ALL_AzimuthAngles
            << h_3gAzimuthAngles
            << h_missingChaAzimuthAngles
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": 3 Photons Polarangles in CB/TAPS")
            << h_3gPolarAnglesCB
            << h_3gPolarAnglesTAPS
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Particle Kinematics: Energy vs angles")
            << drawoption("pcolz")
            << h_w_EvTheta
            << h_wg_EvTheta
            << h_wpi0_EvTheta
            << h_missingP_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Photon energy vs angles CB/TAPS")
            << drawoption("pcolz")
            << h_3g_EvTheta_CB
            << h_3g_EvTheta_TAPS
            << endc; // actually draws the canvas

    for (unsigned int i=0; i<nrCuts_BackToBack; i++){

            stringstream ss;
            ss << i;
            string myString = ss.str();

            ant::canvas(GetName()+": w_pi0g angles in cm of rec. omega meson after "+myString+ " applied cuts")
                << h_pi0g_BackToBack[i]
                << endc; // actually draws the canvas

    }
*/


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

   ant::canvas(GetName()+": Deposited Calo vs Veto energy of expected charged candidates.")
           << drawoption("pcolz")
           << h_ChaCaloVSVetoEnergies_CB
           << h_ChaCaloVSVetoEnergies_TAPS
           << endc; // actually draws the canvas


    for (unsigned int i=0; i<nrCuts_IM; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. particles invariant masses after "+myString+ " applied cuts")
                << h_proton_Im[i]
                << h_wOnly3g_Im[i]
                << endc; // actually draws the canvas

    }

    for (unsigned int i=0; i<nrCuts_IM_pi0; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. pi0 invariant mass after "+myString+ " applied cuts")
                << h_Pi0Only2g_Im[i]
                << endc; // actually draws the canvas

    }

    ant::canvas(GetName()+": Rec. Protons & Omega mesons EvsTheta distribution")
            << drawoption("pcolz")
            << h_w_EvTheta
            << h_p_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Rec. wpi0g decay EvsTheta distribution")
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
