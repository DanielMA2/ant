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
    const BinSettings Veto_Energy_bins(1000, 0, 40);
    const BinSettings initialBeamE_bins(1000,1,1.6);
    const BinSettings cos_bins(100,-1,1);
    const BinSettings sqrt_S_bins(100,s_square_min,s_square_max);
    const BinSettings phi_bins(100, -180, 180);

    const BinSettings theta_bins(180, 0, 180);
    const BinSettings theta_bins_CB(140, 19, 160);
    const BinSettings theta_bins_TAPS(25, 0, 25);

    const BinSettings Im_missing_proton_bins(100, 0, 1600);
    const BinSettings Im_omega_bins(100, 0, 1600);
    const BinSettings Im_pi0_bins(100, 0, 1600);

    const BinSettings E_photon_bins(100, 0, 1600);
    const BinSettings E_missingP_bins(100, 0, 1600);
    const BinSettings E_omega_bins(100, 0, 1600);
    const BinSettings E_pi0_bins(100, 0, 1600);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    auto hfTagger = new HistogramFactory("Tagger", HistFac, "");
    auto hfOnlyEnergy = new HistogramFactory("All the candidate energies", HistFac, "");
    auto hfOnlyAngles = new HistogramFactory("All the candidate angles", HistFac, "");
    auto hfCrossCheck = new HistogramFactory("Cross section check", HistFac, "");
    auto hf_wpi0g_EvsTheta = new HistogramFactory("Energy vs theta dists", HistFac, "");
    auto hf_wpi0g_IM = new HistogramFactory("Invariant mass dists", HistFac, "");
    auto hf_extras = new HistogramFactory("Some additional extra-hists", HistFac, "");

    for (unsigned int i=0; i<nrCuts; i++){

        h_missing_proton_Im[i] = hf_wpi0g_IM->makeTH1D("Missing particle invariant mass",     // title
                                                             "M [MeV]", "#",     // xlabel, ylabel
                                                             Im_missing_proton_bins,  // our binnings
                                                             "h_missing_proton_Im", false     // ROOT object name, auto-generated if omitted
                                                             );
        h_wOnly3g_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. Omega invariant mass",     // title
                                            "M [MeV]", "#",     // xlabel, ylabel
                                            Im_omega_bins,  // our binnings
                                            "h_wOnly3g_Im", false     // ROOT object name, auto-generated if omitted
                                            );
        h_Pi0Only2g_Im[i] = hf_wpi0g_IM->makeTH1D("Rec. Pi0 invariant mass",     // title
                                             "M [MeV]", "#",     // xlabel, ylabel
                                             Im_pi0_bins,  // our binnings
                                             "h_Pi0Only2g_Im", false     // ROOT object name, auto-generated if omitted
                                             );
        h_pi0g_BackToBack[i] = hfOnlyAngles->makeTH1D("w_pi0g angles in cm of rec. omega meson",     // title
                                                         "cos_Theta* [deg]", "#",     // xlabel, ylabel
                                                         cos_bins,  // our binnings
                                                         "h_pi0g_BackToBack", false     // ROOT object name, auto-generated if omitted
                                                         );

    }

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


    h_ALL_PolarAngles = hfOnlyAngles->makeTH1D("All candidate polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_ALL_PolarAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_ALL_AzimuthAngles = hfOnlyAngles->makeTH1D("All candidate azimuthangles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_ALL_AzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_missingChaPolarAngles = hfOnlyAngles->makeTH1D("Missing particle polarangles",     // title
                                                    "#Theta [deg]", "#",     // xlabel, ylabel
                                                    theta_bins,  // our binnings
                                                    "h_missingChaPolarAngles", false     // ROOT object name, auto-generated if omitted
                                                    );

    h_missingChaAzimuthAngles = hfOnlyAngles->makeTH1D("Missing particle azimuthangles",     // title
                                                    "#Phi [deg]", "#",     // xlabel, ylabel
                                                    phi_bins,  // our binnings
                                                    "h_missingChaAzimuthAngles", false     // ROOT object name, auto-generated if omitted
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

    h_missingP_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Missing particle energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_missingP_bins,    // our binnings
                                       "h_missingP_EvTheta", true    // ROOT object name, auto-generated if omitted
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

    h_doubly_wp_DCS_reconstructed_lab = hfCrossCheck->makeTH2D("Cross section check in lab-frame", //title
                                                                        "cos_Theta [deg]","sqrt_S [GeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_wp_DCS_reconstructed_lab", true    // ROOT object name, auto-generated if omitted
                                                                        );
    
    h_doubly_wp_DCS_reconstructed_cmFrame = hfCrossCheck->makeTH2D("Cross section check in cm-frame", //title
                                                                        "cos_Theta* [deg]","sqrt_S [GeV]", // xlabel, ylabel
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
    vector<double> The;
    vector<double> Phi;
    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;
    vector<double> wpi0g_times;

    for(const auto& cand : candidates.get_iter()) {
        h_VetoEnergies->Fill(cand->VetoEnergy);
        The.push_back(cand->Theta*radtodeg);
        Phi.push_back(cand->Phi*radtodeg);
        if(cand->VetoEnergy == 0) {
            neutral.emplace_back(cand);
            NeuCanCluSize.push_back(cand->ClusterSize);
            NeuCanCaloE.push_back(cand->CaloEnergy);
        }       
        else{
            charged.emplace_back(cand);
            ChaCanCluSize.push_back(cand->ClusterSize);
            ChaCanCaloE.push_back(cand->CaloEnergy);
        }

    }

    //getting access to the pi0 decay photons

    double wpi0gOnly3g_IM=-100;
    TParticle w;
    TParticle pi0;
    TParticle g1;
    TParticle g2;
    TParticle g3;

    if(NeuCanCaloE.size() == 3){

        g1 = TParticle(ParticleTypeDatabase::Photon, neutral[0]);
        g2 = TParticle(ParticleTypeDatabase::Photon, neutral[1]);
        g3 = TParticle(ParticleTypeDatabase::Photon, neutral[2]);

        /*
        auto g1channel = neutral[0]->FindCaloCluster()->CentralElement;
        auto g2channel = neutral[1]->FindCaloCluster()->CentralElement;
        auto g3channel = neutral[2]->FindCaloCluster()->CentralElement;       
        if(neutral[0]->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
        {isg1CB = true;}
        if(neutral[0]->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            if(taps->IsPbWO4(g1channel)){isg1PbWO = true;}else{isg1BaF2 = true;}
        }
        //if(!isg1CB && !isg1BaF2 && !isg1PbWO)
        //{isg1undetected = true;}

        if(neutral[1]->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
        {isg2CB = true;}
        if(neutral[1]->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            if(taps->IsPbWO4(g2channel)){isg2PbWO = true;}else{isg2BaF2 = true;}
        }
        //if(!isg2CB && !isg2BaF2 && !isg2PbWO)
        //{isg2undetected = true;}
        */

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

    TParticle proton_missing_particle;
    TParticle omega;
    TParticle pion;
    for (const auto& taggerhit : event.Reconstructed().TaggerHits) { // Event loop

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());

        proton_missing_particle = TParticle(ParticleTypeDatabase::Proton,(TLorentzVector)(InitialPhotonVec+InitialProtonVec-(g1+g2+g3)));
        omega = TParticle(ParticleTypeDatabase::Omega,(TLorentzVector)(g1+g2+g3));

        TLorentzVector Lw = (TLorentzVector)(omega);
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);
        TLorentzVector LwBoosted = Lw;
        LwBoosted.Boost(-Linitial.BoostVector());

        //Adding selections & filling histograms

        const double weight = promptrandom.FillWeight();

        h_TaggerTime->Fill(taggerhit.Time, weight);
        h_InitialBeamE->Fill(InitialPhotonVec.E()/1000,weight);

        /*
        if(NeuCanCaloE.size() == 2){

            h_ProtonIm->Fill(proton_missing_particle.M(),weight);
            t.ProtonIm = proton_missing_particle.M();

            if(proton_missing_particle.M()<(mp+sigma) && proton_missing_particle.M()>(mp-sigma)){

                h_PhotonIm->Fill(wpi0g_3g_IM,weight);
                t.PhotonIm = wpi0g_3g_IM;

                for (unsigned int i=0; i<The.size(); i++){
                    h_PolarAngles->Fill(The[i],weight);
                }

                for (unsigned int i=0; i<Phi.size(); i++){
                    h_AzimuthAngles->Fill(Phi[i],weight);
                }

                for (unsigned int i=0; i<PhotonThe.size(); i++){
                    h_PhotonPolarAngles->Fill(PhotonThe[i],weight);
                    h_PhotonETheta->Fill(PhotonThe[i],PhotonE[i],weight);
                    //if (PhotonThe[i] < 20){PhotonPolarAngles_TAPS += 1.0;}
                    //if (PhotonThe[i] >= 20){PhotonPolarAngles_CB += 1.0;}
                    h_doubly_ap_DCS_reconstructed_lab->Fill(cos(Lw.Theta()),Linitial.M(),weight);
                    h_doubly_ap_DCS_reconstructed_cmFrame->Fill(cos(Linitial.Angle(LwBoosted.Vect())),Linitial.M(),weight);
                }

                for (unsigned int i=0; i<PhotonTheCB.size(); i++){
                    h_PhotonPolarAngles_CB->Fill(PhotonTheCB[i],weight);
                    h_PhotonETheta_CB->Fill(PhotonTheCB[i],PhotonE[i],weight);
                }

                for (unsigned int i=0; i<PhotonTheTAPS.size(); i++){
                    h_PhotonPolarAngles_TAPS->Fill(PhotonTheTAPS[i],weight);
                    h_PhotonETheta_TAPS->Fill(PhotonTheTAPS[i],PhotonE[i],weight);
                }

                h_ProtonPolarAngles->Fill(proton_missing_particle.Theta()*radtodeg,weight);
                h_ProtonETheta->Fill(proton_missing_particle.Theta()*radtodeg,proton_missing_particle.E,weight);
                if (proton_missing_particle.Theta()*radtodeg < 20){ProtonPolarAngles_TAPS += 1.0;}
                if (proton_missing_particle.Theta()*radtodeg >= 20 && proton_missing_particle.Theta()*radtodeg < 160){ProtonPolarAngles_CB += 1.0;}

                for (unsigned int i=0; i<PhotonPhi.size(); i++){
                    h_PhotonAzimuthAngles->Fill(PhotonPhi[i],weight);
                    h_PhotonEPhi->Fill(PhotonPhi[i],PhotonE[i],weight);
                }

                h_ProtonAzimuthAngles->Fill(proton_missing_particle.Phi()*radtodeg,weight);
                h_ProtonEPhi->Fill(proton_missing_particle.Phi()*radtodeg,proton_missing_particle.E,weight);

                h_ProtonPeak->Fill(proton_missing_particle.M(),weight);

                //Filling the statistical histogram

                if (isg1CB && isg2CB){count_first+=promptrandom.FillWeight();}
                if ((isg1CB && isg2BaF2)||(isg2CB && isg1BaF2)){count_second+=promptrandom.FillWeight();}
                if ((isg1CB && isg2PbWO)||(isg2CB && isg1PbWO)){count_third+=promptrandom.FillWeight();}
                if (isg1BaF2 && isg2BaF2){count_fourth+=promptrandom.FillWeight();}
                if ((isg1BaF2 && isg2PbWO)||(isg2BaF2 && isg1PbWO)){count_fifth+=promptrandom.FillWeight();}
                if (isg1PbWO && isg2PbWO){count_sixth+=promptrandom.FillWeight();}
                //if ((isg1PbWO && isg2undetected)||(isg2PbWO && isg1undetected)){count_seventh++;}
                //if ((isg1BaF2 && isg2undetected)||(isg2BaF2 && isg1undetected)){count_eigth++;}
                //if ((isg1CB && isg2undetected)||(isg2CB && isg1undetected)){count_nineth++;}
                //if (isg1undetected && isg2undetected){count_tenth++;}

                weight_res += weight;
            }
        }
        */

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

    for (unsigned int i=0; i<nrCuts; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. particles invariant masses after "+myString+ " applied cuts")
                << h_missing_proton_Im[i]
                << h_wOnly3g_Im[i]
                << h_Pi0Only2g_Im[i]
                << endc; // actually draws the canvas

        ant::canvas(GetName()+": w_pi0g angles in cm of rec. omega meson after "+myString+ " applied cuts")
                << h_pi0g_BackToBack[i]
                << endc; // actually draws the canvas

    }

    ant::canvas(GetName()+": Cross section checks")
            << drawoption("Surf")
            << h_doubly_wp_DCS_reconstructed_lab
            << h_doubly_wp_DCS_reconstructed_cmFrame
            << endc; // actually draws the canvas

}

void scratch_damaurer_MC_pw_ppi0g_p3g::Finish()
{

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
