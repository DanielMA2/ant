#include "pw_ppi0g_p3g_kinFit.h"
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

APLCON::Fit_Settings_t scratch_damaurer_pw_ppi0g_p3g_kinFit::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

scratch_damaurer_pw_ppi0g_p3g_kinFit::scratch_damaurer_pw_ppi0g_p3g_kinFit(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    /*
    fitter(utils::UncertaintyModels::Interpolated::makeAndLoad(),
    true  // enable fit z vertex
    )
    */

    /*
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                             utils::UncertaintyModels::Interpolated::Type_t::MC,
                             make_shared<utils::UncertaintyModels::FitterSergey>())),
    */

    fit_model_data(utils::UncertaintyModels::Interpolated::makeAndLoad(
                   utils::UncertaintyModels::Interpolated::Type_t::Data,
                   // use Sergey as starting point
                   make_shared<utils::UncertaintyModels::FitterSergey>()
                  )
          ),

    fit_model_mc(utils::UncertaintyModels::Interpolated::makeAndLoad(
                 utils::UncertaintyModels::Interpolated::Type_t::MC,
                 // use Sergey as starting point
                 make_shared<utils::UncertaintyModels::FitterSergey>()
                )
          ),

    fitter(nullptr, opts->Get<bool>("FitZVertex", true))

{

    fitter.SetZVertexSigma(3.0);
    taps = make_shared<expconfig::detector::TAPS_2013_11>(false, false, false);

    //const BinSettings time_bins(200, 0,10);
    //const BinSettings timeDiffCorTaggCB_bins(200, -5,10);
    //const BinSettings timeDiffCorTaggTAPS_bins(200, -5,20);

    const BinSettings bins_tagger_time(2000, -200, 200);

    const BinSettings Veto_Energy_bins(500, 0, 12);
    const BinSettings Calo_Energy_bins(500, 0, 1200);
    const BinSettings initialBeamE_bins(200,0, 1600);
    const BinSettings cos_bins(100,-1,1);
    const BinSettings statistic_bins(number_of_bins,0,nrCuts_total);
    const BinSettings sqrt_S_bins(100,s_square_min,s_square_max);

    const BinSettings phi_bins(720, -180, 180);
    const BinSettings theta_bins(360, 0, 180);
    const BinSettings PIDtime_bins(482,-60.5,60.5);

    const BinSettings Im_proton_bins(200, 0, 1500);
    const BinSettings Im_omega_bins(200, 0, 1200);
    const BinSettings Im_pi0_bins(200, 0, 1000);

    //const BinSettings Ekin_neu_bins(200, 0, 1600);
    //const BinSettings Ekin_cha_bins(200, 0, 1000);

    const BinSettings E_photon_bins(200, 0, 1600);
    const BinSettings E_proton_bins(200, 900, 2000);
    const BinSettings E_omega_bins(200, 400, 1800);
    const BinSettings E_pi0_bins(200, 0, 1600);
    const BinSettings CB_Esum_bins(250,0.,2000.);
    //const BinSettings cos_bins_BackToBack(1000,-1,1);
    //const BinSettings wpi0_Q_bins(100, 0, 1500);
    //const BinSettings proton_theta_bins(180, 0, 90);

    //KinFit binSettings:
    const BinSettings defEk_bins = BinSettings(201,-200,200);
    const BinSettings defPhi_bins = BinSettings(101,-100,100);
    const BinSettings defTheta_bins = BinSettings(101,-100,100);
    vector<BinSettings> partTypeEkBins;
    partTypeEkBins.emplace_back(E_proton_bins);
    partTypeEkBins.emplace_back(E_photon_bins);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    auto hf_missingP_IM = new HistogramFactory("hf_missingP_IM", HistFac, "");
    auto hf_3g_IM = new HistogramFactory("hf_3g_IM", HistFac, "");
    auto hf_2gComb_IM = new HistogramFactory("hf_2gComb_IM", HistFac, "");
    auto hf_2g_OpeningAngles = new HistogramFactory("hf_2g_OpeningAngles", HistFac, "");
    auto hf_2gPi0_IM = new HistogramFactory("hf_2gPi0_IM", HistFac, "");
    auto hf_BackToBack = new HistogramFactory("hf_BackToBack", HistFac, "");
    auto hf_EvsTheta_CB = new HistogramFactory("hf_EvsTheta_CB", HistFac, "");
    auto hf_EvsTheta_TAPS = new HistogramFactory("hf_EvsTheta_TAPS", HistFac, "");
    auto hf_EvsPhi = new HistogramFactory("hf_EvsPhi", HistFac, "");
    auto hf_CrossSection = new HistogramFactory("hf_CrossSection", HistFac, "");
    auto hf_CaloEvsVetoE = new HistogramFactory("hf_CaloEvsVetoE", HistFac, "");
    auto hf_RecData_CandStat = new HistogramFactory("hf_RecData_CandStat", HistFac, "");
    auto hf_EvsT_CB = new HistogramFactory("hf_EvsT_CB", HistFac, "");
    auto hf_EvsT_TAPS = new HistogramFactory("hf_EvsT_TAPS", HistFac, "");
    auto hf_Extras = new HistogramFactory("hf_Extras", HistFac, "");

    auto hf_Tagger = new HistogramFactory("hf_Tagger", HistFac, "");
    auto hf_Candidates = new HistogramFactory("hf_Candidates", HistFac, "");
    auto hf_wpi0g_EvsTheta = new HistogramFactory("hf_wpi0g_EvsTheta", HistFac, "");

    auto hf_VetoE = new HistogramFactory("hf_VetoE", HistFac, "");
    auto hf_CBEsum = new HistogramFactory("hf_CB_Esum", HistFac, "");

    //auto hf_OnlyEnergy = new HistogramFactory("All_the_candidate_energies", HistFac, "");
    //auto hf_OnlyAngles = new HistogramFactory("All_the_candidate_angles", HistFac, "");
    //auto hf_extras = new HistogramFactory("Some_additional_extra_hists", HistFac, "");

    //KinFit histfactories:
    auto hf_OverviewKF = new HistogramFactory("hf_OverviewKF", HistFac, "");
    auto hf_CombCBKF = new HistogramFactory("hf_CombCBKF", HistFac, "");
    auto hf_CombTAPSKF = new HistogramFactory("hf_CombTAPSKF", HistFac, "");
    auto hf_invmassKF = new HistogramFactory("hf_invmassKF", HistFac, "");
    auto hfPullsCBKF = new HistogramFactory("hf_PullsCBKF", HistFac, "");
    auto hfPullsTAPSKF = new HistogramFactory("hf_PullsTAPSKF", HistFac, "");
    auto hf_CaloEvsVetoE_CB_SideCheck = new HistogramFactory("hf_CaloEvsVetoE_CB_SideCheck", HistFac, "");
    auto hf_CaloEvsVetoE_TAPS_SideCheck = new HistogramFactory("hf_CaloEvsVetoE_TAPS_SideCheck", HistFac, "");
    auto hf_PIDEvsTime = new HistogramFactory("hf_PIDEvsTime", HistFac, "");

    h_RecData_Stat = hf_RecData_CandStat->makeTH1D("Amount of events after cuts:",     // title
                                     "Cuts", "#",     // xlabel, ylabel
                                     statistic_bins,  // our binnings
                                     "h_RecData_Stat", true    // ROOT object name, auto-generated if omitted
                                     );

    h_RecData_relStat = hf_RecData_CandStat->makeTH1D("rel. amount of events after cuts:",     // title
                                      "Cuts", "rel. #",     // xlabel, ylabel
                                      statistic_bins,  // our binnings
                                      "h_RecData_relStat", true    // ROOT object name, auto-generated if omitted
                                      );

    h_nClusters = hf_Tagger->makeTH1D("Number of Clusters", "nClusters", "#", BinSettings(20), "h_nClusters");


    /*
    h_InitialBeamE = hf_OnlyEnergy->makeTH1D("Photon beam - Energy distribution",     // title
                                   "E [GeV]", "#",     // xlabel, ylabel
                                   initialBeamE_bins,  // our binnings
                                   "h_InitialBeamE"     // ROOT object name, auto-generated if omitted
                                   );

    h_VetoEnergies = hf_Tagger->makeTH1D("Veto Energies", "E [MeV]", "#", Veto_Energy_bins, "h_VetoEnergies");


    h_NeuPolarAngles = hf_OnlyAngles->makeTH1D("Neutral candidate polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_NeuPolarAngles", true     // ROOT object name, auto-generated if omitted
                                       );

    h_NeuAzimuthAngles = hf_OnlyAngles->makeTH1D("Neutral candidate azimuthangles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_NeuAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                       );

    h_ChaPolarAngles = hf_OnlyAngles->makeTH1D("Charged particles polarangles",     // title
                                                    "#Theta [deg]", "#",     // xlabel, ylabel
                                                    theta_bins,  // our binnings
                                                    "h_ChaPolarAngles", true     // ROOT object name, auto-generated if omitted
                                                    );

    h_ChaAzimuthAngles = hf_OnlyAngles->makeTH1D("Charged particles azimuthangles",     // title
                                                    "#Phi [deg]", "#",     // xlabel, ylabel
                                                    phi_bins,  // our binnings
                                                    "h_ChaAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                                    );

    h_3gPolarAngles = hf_OnlyAngles->makeTH1D("Rec. Photons polarangles",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_3gPolarAngles", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesCB = hf_OnlyAngles->makeTH1D("Rec. Photons polarangles only in CB",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_3gPolarAnglesCB", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesTAPS = hf_OnlyAngles->makeTH1D("Rec. Photons polarangles only in TAPS",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_3gPolarAnglesTAPS", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gAzimuthAngles = hf_OnlyAngles->makeTH1D("Rec. Photons azimuthangles",     // title
                                               "#Phi [deg]", "#",     // xlabel, ylabel
                                               phi_bins,  // our binnings
                                               "h_3gAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                               );

    h_3g_EvTheta_CB = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in CB", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_CB", true    // ROOT object name, auto-generated if omitted
                                       );

    h_3g_EvTheta_TAPS = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in TAPS", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_TAPS", true    // ROOT object name, auto-generated if omitted
                                       );
    */

    for (unsigned int i=0; i<nrCuts_total; i++){
        /*
        stringstream ss;
        ss << i;
        string myString = ss.str();
        */

        h_nCandidates[i] = hf_Candidates->makeTH1D("Number of Candidates", "nCandidates "+cuts[i], "#", BinSettings(10), "h_nCandidates_"+cuts[i]);
        h_nNeuChaCandidates[i] = hf_Candidates->makeTH2D("Number of charged vs neutral candidates "+cuts[i], "nChaCand", "nNeuCand", BinSettings(6), BinSettings(8), "h_nNeuChaCandidates_"+cuts[i], true);

        h_TaggerTime[i] = hf_Tagger->makeTH1D("Tagger Time "+cuts[i],     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        bins_tagger_time,  // our binnings
                                        "h_TaggerTime_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                        );

        h_AllCaloEvsVetoE_CB[i] = hf_CaloEvsVetoE->makeTH2D("Deposited Calo- VS Veto-energies in CB "+cuts[i],     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_AllCaloEvsVetoE_CB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_AllCaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE->makeTH2D("Deposited Calo- VS Veto-energies in TAPS "+cuts[i],     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_AllCaloEvsVetoE_TAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_neuEkinVSThetaCB[i] = hf_EvsTheta_CB->makeTH2D("Photons Ekin vs Theta in CB "+cuts[i], //title
                                                "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                theta_bins, Calo_Energy_bins,    // our binnings
                                                "h_neuEkinVSThetaCB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                );

        h_neuEkinVSThetaTAPS[i] = hf_EvsTheta_TAPS->makeTH2D("Photons Ekin vs Theta in TAPS "+cuts[i], //title
                                                "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                theta_bins, Calo_Energy_bins,    // our binnings
                                                "h_neuEkinVSThetaTAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                );

        h_neuEkinVSPhi[i] = hf_EvsPhi->makeTH2D("Photons Ekin vs Phi "+cuts[i], //title
                                            "#Phi [deg]","E_kin [MeV]",  // xlabel, ylabel
                                            phi_bins, Calo_Energy_bins,    // our binnings
                                            "h_neuEkinVSPhi_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                            );

        h_chaEkinVSThetaCB[i] = hf_EvsTheta_CB->makeTH2D("Charged Ekin vs Theta in CB "+cuts[i], //title
                                                  "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                  theta_bins, Calo_Energy_bins,    // our binnings
                                                  "h_chaEkinVSThetaCB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                  );

        h_chaEkinVSThetaTAPS[i] = hf_EvsTheta_TAPS->makeTH2D("Charged Ekin vs Theta in TAPS "+cuts[i], //title
                                                  "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                  theta_bins, Calo_Energy_bins,    // our binnings
                                                  "h_chaEkinVSThetaTAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                  );

        h_chaEkinVSPhi[i] = hf_EvsPhi->makeTH2D("Charged Ekin vs Phi "+cuts[i], //title
                                              "#Phi [deg]","E_kin [MeV]",  // xlabel, ylabel
                                              phi_bins, Calo_Energy_bins,    // our binnings
                                              "h_chaEkinVSPhi_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                              );

        h_neuCaloE[i] = hf_Extras->makeTH1D("Neutral Calorimeter energy "+cuts[i], //title
                                            "Calo E [MeV]", "#",  // xlabel, ylabel
                                            Calo_Energy_bins,    // our binnings
                                            "h_neuCaloE_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                            );

        h_chaCaloE[i] = hf_Extras->makeTH1D("Charged Calorimeter energy "+cuts[i], //title
                                            "Calo E [MeV]", "#",  // xlabel, ylabel
                                            Calo_Energy_bins,    // our binnings
                                            "h_chaCaloE_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                            );

        h_beamE[i] = hf_Extras->makeTH1D("Initial photon beam energy "+cuts[i],     // title
                                                             "E_{photonbeam} [MeV]", "#",     // xlabel, ylabel
                                                             initialBeamE_bins,  // our binnings
                                                             "h_beamE_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_neuEvsT_CB[i] = hf_EvsT_CB->makeTH2D("Neutrals: Calo energy vs time in CB "+cuts[i],     // title
                                                             "Calo E [MeV]","t_{CB} [ns]",     // xlabel, ylabel
                                                             Calo_Energy_bins,PIDtime_bins,  // our binnings
                                                             "h_neuEvsT_CB_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_neuEvsT_TAPS[i] = hf_EvsT_TAPS->makeTH2D("Neutrals: Calo energy vs time in TAPS "+cuts[i],     // title
                                                             "Calo E [MeV]","t_{TAPS} [ns]",     // xlabel, ylabel
                                                             Calo_Energy_bins,PIDtime_bins,  // our binnings
                                                             "h_neuEvsT_TAPS_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_chaEvsT_CB[i] = hf_EvsT_CB->makeTH2D("Charged: Calo energy vs time in CB "+cuts[i],     // title
                                                             "Calo E [MeV]","t_{CB} [ns]",     // xlabel, ylabel
                                                             Calo_Energy_bins,PIDtime_bins,  // our binnings
                                                             "h_chaEvsT_CB_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_chaEvsT_TAPS[i] = hf_EvsT_TAPS->makeTH2D("Charged: Calo energy vs time in TAPS "+cuts[i],     // title
                                                             "Calo E [MeV]","t_{TAPS} [ns]",     // xlabel, ylabel
                                                             Calo_Energy_bins,PIDtime_bins,  // our binnings
                                                             "h_chaEvsT_TAPS_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_AllVetoE_CB[i] = hf_VetoE->makeTH1D("Deposited veto energies in CB "+cuts[i],     // title
                                         "Veto E [MeV]", "#",  // xlabel, ylabel
                                         Veto_Energy_bins,    // our binnings
                                         "h_AllVetoE_CB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_AllVetoE_TAPS[i] = hf_VetoE->makeTH1D("Deposited veto energies in TAPS "+cuts[i],     // title
                                         "Veto E [MeV]", "#",  // xlabel, ylabel
                                         Veto_Energy_bins,    // our binnings
                                         "h_AllVetoE_TAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_CBEsum[i] = hf_CBEsum->makeTH1D("CB Esum "+cuts[i],"CB Esum [MeV]","#",CB_Esum_bins,"h_CBEsum_"+cuts[i],true);

    }

    for (unsigned int i=0; i<nrCuts_Sel; i++){
        /*
        stringstream ss;
        ss << i;
        string myString = ss.str();
        */

        h_neuLowestCaloE[i] = hf_Extras->makeTH1D("Lowest neutral energy clusters "+cuts[i+nrCuts_beforeSel], //title
                                                  "Calo E [MeV]", "#",  // xlabel, ylabel
                                                  Calo_Energy_bins,    // our binnings
                                                  "h_neuLowestCaloE_"+cuts[i+nrCuts_beforeSel], true    // ROOT object name, auto-generated if omitted
                                                  );

        h_missingP_IM[i] = hf_missingP_IM->makeTH1D("Im(missingP) "+cuts[i+nrCuts_beforeSel],     // title
                                                             "IM(missingP) [MeV]", "#",     // xlabel, ylabel
                                                             Im_proton_bins,  // our binnings
                                                             "h_missingP_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_3g_IM[i] = hf_3g_IM->makeTH1D("IM(3g) "+cuts[i+nrCuts_beforeSel],     // title
                                            "IM(3g) [MeV]", "#",     // xlabel, ylabel
                                            Im_omega_bins,  // our binnings
                                            "h_3g_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                            );

        h_2gComb_IM[i] = hf_2gComb_IM->makeTH1D("IM(2g) combinations "+cuts[i+nrCuts_beforeSel],     // title
                                             "IM(2g) [MeV]", "#",     // xlabel, ylabel
                                             Im_pi0_bins,  // our binnings
                                             "h_2gComb_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                             );

        h_2gComb_OpeningAngles[i] = hf_2g_OpeningAngles->makeTH1D("2g opening angle combinations "+cuts[i+nrCuts_beforeSel],     // title
                                             "2g comb opening Angles [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_2gComb_OpeningAngles_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                             );

        h_2gLowestClusterE_OpeningAngles[i] = hf_2g_OpeningAngles->makeTH1D("2g opening angles of low energetic clusters "+cuts[i+nrCuts_beforeSel],     // title
                                             "2g opening Angles [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_2gLowestClusterE_OpeningAngles_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                             );

        h_wp_BackToBack[i] = hf_BackToBack->makeTH1D("Omega Proton Back-To-Back check "+cuts[i+nrCuts_beforeSel],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_wp_BackToBack_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                         );

        h_doublyDCScm_gp_wp[i] = hf_CrossSection->makeTH2D("gp_wp cross section in cm-frame "+cuts[i+nrCuts_beforeSel], //title
                                                                            "cos_Theta*","sqrt_S [MeV]", // xlabel, ylabel
                                                                            cos_bins, sqrt_S_bins,    // our binnings
                                                                            "h_doublyDCScm_gp_wp_"+cuts[i+nrCuts_beforeSel], true    // ROOT object name, auto-generated if omitted
                                                                            );

        h_Proton_PIDEvsTime[i] = hf_PIDEvsTime->makeTH2D("Protons: PID energy vs time "+cuts[i+nrCuts_beforeSel],     // title
                                                                          "PID time [ns]", "PID E [MeV]",  // xlabel, ylabel
                                                                          PIDtime_bins,BinSettings(200,0.,20),    // our binnings
                                                                          "h_Proton_PIDEvsTime_"+cuts[i+nrCuts_beforeSel], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Proton_CaloEvsVetoE_CB[i] = hf_CaloEvsVetoE_CB_SideCheck->makeTH2D("Protons: Deposited calo vs veto energies in CB "+cuts[i+nrCuts_beforeSel],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_Proton_CaloEvsVetoE_CB_"+cuts[i+nrCuts_beforeSel], true    // ROOT object name, auto-generated if omitted
                                         );

        h_Proton_CaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE_TAPS_SideCheck->makeTH2D("Protons: Deposited calo vs veto energies in TAPS "+cuts[i+nrCuts_beforeSel],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_Proton_CaloEvsVetoE_TAPS_"+cuts[i+nrCuts_beforeSel], true    // ROOT object name, auto-generated if omitted
                                         );

        h_2gPi0_IM[i] = hf_2gPi0_IM->makeTH1D("IM(2gPi0) "+cuts[i+nrCuts_beforeSel],     // title
                                         "IM(2gPi0) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2gPi0_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                         );

        h_EpivsIM3g[i] = hf_Extras->makeTH2D("E_{#Pi⁰} vs IM(3g)"+cuts[i+nrCuts_beforeSel],     // title
                                        "IM(3g) [MeV]", "E_{#Pi}",     // xlabel, ylabel
                                        Im_omega_bins,E_pi0_bins,  // our binnings
                                        "h_EpivsIM3g_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                        );

        h_wpi0g_BackToBack[i] = hf_BackToBack->makeTH1D("Pi0 Gamma Back-To-Back check "+cuts[i+nrCuts_beforeSel],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_wpi0g_BackToBack_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                         );

        h_pi0gg_BackToBack[i] = hf_BackToBack->makeTH1D("2Gamma Back-To-Back check "+cuts[i+nrCuts_beforeSel],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_pi0gg_BackToBack_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                         );



    }

/*
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforePi0); i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        h_2gPi0_IM[i] = hf_2gPi0_IM->makeTH1D("IM(2gPi0) "+cuts[i+nrCuts_beforePi0],     // title
                                         "IM(2gPi0) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2gPi0_IM_"+cuts[i+nrCuts_beforePi0], true     // ROOT object name, auto-generated if omitted
                                         );

        h_EpivsIM3g[i] = hf_Extras->makeTH2D("E_{#Pi⁰} vs IM(3g)"+cuts[i+nrCuts_beforePi0],     // title
                                        "IM(3g) [MeV]", "E_{#Pi}",     // xlabel, ylabel
                                        Im_omega_bins,E_pi0_bins,  // our binnings
                                        "h_EpivsIM3g_"+cuts[i+nrCuts_beforePi0], true     // ROOT object name, auto-generated if omitted
                                        );

        h_wpi0g_BackToBack[i] = hf_BackToBack->makeTH1D("Pi0 Gamma Back-To-Back check "+cuts[i+nrCuts_beforePi0],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_wpi0g_BackToBack_"+cuts[i+nrCuts_beforePi0], true     // ROOT object name, auto-generated if omitted
                                         );

        h_pi0gg_BackToBack[i] = hf_BackToBack->makeTH1D("2Gamma Back-To-Back check "+cuts[i+nrCuts_beforePi0],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_pi0gg_BackToBack_"+cuts[i+nrCuts_beforePi0], true     // ROOT object name, auto-generated if omitted
                                         );

    }
    */

    //KinFit-overview:
    //h_Steps = hf_OverviewKF->makeTH1D("Steps","","",BinSettings(10),"h_Steps");
    for (unsigned int i=0; i<nrCuts_KF; i++){
        h_Probability[i] = hf_OverviewKF->makeTH1D(Form("Fit probability %s",cuts[i+nrCuts_beforeKF].c_str()),"P(#chi^{2})","#",BinSettings(1000,0.,1.),Form("h_Probability_%s",cuts[i+nrCuts_beforeKF].c_str()),true);
        h_Fit_zvert[i] = hf_OverviewKF->makeTH1D(Form("Fitted z-vertex %s",cuts[i+nrCuts_beforeKF].c_str()),"z [cm]","#",BinSettings(50,-15,15),Form("h_Fit_zvert_%s",cuts[i+nrCuts_beforeKF].c_str()),true);
        h_fitEbeam[i] = hf_OverviewKF->makeTH1D(Form("Fitted beam energy %s",cuts[i+nrCuts_beforeKF].c_str()),"E [MeV]","#",initialBeamE_bins,Form("h_fitEbeam_%s",cuts[i+nrCuts_beforeKF].c_str()),true);
        h_IM3g_Fit[i] = hf_invmassKF->makeTH1D(Form("Fitted photons invariant mass %s",cuts[i+nrCuts_beforeKF].c_str()),"Im(3g_fit) [MeV]","#",Im_omega_bins,Form("h_IM3g_Fit_%s",cuts[i+nrCuts_beforeKF].c_str()),true);
        h_IM2gPi0_Fit[i] = hf_invmassKF->makeTH1D(Form("IM(2gPi0) for fitted photons %s",cuts[i+nrCuts_beforeKF].c_str()),"IM(2gPi0) [MeV]", "#",Im_pi0_bins,Form("h_IM2gPi0_Fit_%s",cuts[i+nrCuts_beforeKF].c_str()), true);

        h_p_EvTheta[i] = hf_wpi0g_EvsTheta->makeTH2D(Form("Rec. proton energy vs theta after %s",cuts[i+nrCuts_beforeKF].c_str()), //title
                                                     "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                     theta_bins, E_proton_bins,    // our binnings
                                                     Form("h_p_EvTheta_%s",cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                     );

        h_w_EvTheta[i] = hf_wpi0g_EvsTheta->makeTH2D(Form("Rec. omega meson energy vs theta after %s",cuts[i+nrCuts_beforeKF].c_str()), //title
                                                     "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                     theta_bins, E_omega_bins,    // our binnings
                                                     Form("h_w_EvTheta_%s",cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                     );

        h_wg_EvTheta[i] = hf_wpi0g_EvsTheta->makeTH2D(Form("Rec. omega decay-photon energy vs theta after %s",cuts[i+nrCuts_beforeKF].c_str()), //title
                                           "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                           theta_bins, E_photon_bins,    // our binnings
                                           Form("h_wg_EvTheta_%s",cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                           );

        h_wpi0_EvTheta[i] = hf_wpi0g_EvsTheta->makeTH2D(Form("Rec. Pi0 energy vs theta after %s",cuts[i+nrCuts_beforeKF].c_str()), //title
                                            "Theta [deg]","E [MeV]", // xlabel, ylabel
                                            theta_bins, E_pi0_bins,    // our binnings
                                            Form("h_wpi0_EvTheta_%s",cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                            );

        h_wpi02g_EvTheta[i] = hf_wpi0g_EvsTheta->makeTH2D(Form("Rec. Pi0 decay-photons energy vs theta after %s",cuts[i+nrCuts_beforeKF].c_str()), //title
                                            "Theta [deg]","E [MeV]", // xlabel, ylabel
                                            theta_bins, E_pi0_bins,    // our binnings
                                            Form("h_wpi02g_EvTheta_%s",cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                            );

        for (unsigned int j=0; j<nrPartType; j++){ //0 refers to proton, 1 to photons

            h_Ek_dev_CB[i][j] = hf_CombCBKF->makeTH2D(Form("CB: Energy deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "E_{rec} [MeV]","E_{fit}-E_{rec} [MeV]", // xlabel, ylabel
                                                partTypeEkBins.at(j), defEk_bins,    // our binnings
                                                Form("h_Ek_dev_CB_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            h_Theta_dev_CB[i][j] = hf_CombCBKF->makeTH2D(Form("CB: Theta deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "#Theta_{rec} [deg]","#Theta_{fit}-#Theta_{rec} [deg]", // xlabel, ylabel
                                                theta_bins, defTheta_bins,    // our binnings
                                                Form("h_Theta_dev_CB_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            h_Phi_dev_CB[i][j] = hf_CombCBKF->makeTH2D(Form("CB: Phi deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "#Phi_{rec} [deg]","#Phi_{fit}-#Phi_{rec} [deg]", // xlabel, ylabel
                                                phi_bins, defPhi_bins,    // our binnings
                                                Form("h_Phi_dev_CB_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            h_Ek_dev_TAPS[i][j] = hf_CombTAPSKF->makeTH2D(Form("TAPS: Energy deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "E_{rec} [MeV]","E_{fit}-E_{rec} [MeV]", // xlabel, ylabel
                                                partTypeEkBins.at(j), defEk_bins,    // our binnings
                                                Form("h_Ek_dev_TAPS_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            h_Theta_dev_TAPS[i][j] = hf_CombTAPSKF->makeTH2D(Form("TAPS: Theta deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "#Theta_{rec} [deg]","#Theta_{fit}-#Theta_{rec} [deg]", // xlabel, ylabel
                                                theta_bins, defTheta_bins,    // our binnings
                                                Form("h_Theta_dev_TAPS_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            h_Phi_dev_TAPS[i][j] = hf_CombTAPSKF->makeTH2D(Form("TAPS: Theta deviation of rec. vs fitted %s after %s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                "#Phi_{rec} [deg]","#Phi_{fit}-#Phi_{rec} [deg]", // xlabel, ylabel
                                                phi_bins, defPhi_bins,    // our binnings
                                                Form("h_Phi_dev_TAPS_%s_%s",fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                );

            for (unsigned int k=0; k<nrFitVars; k++){

                  h_PartPulls_CB[i][j][k] = hfPullsCBKF->makeTH1D(Form("CB: %s pulls of fitted %s after %s",fitvarnameCB[k].c_str(),fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                               Form("%s pull",fitvarnameCB[k].c_str()),"#", // xlabel, ylabel
                                                               BinSettings(200,-10.,10.),   // our binnings
                                                               Form("h_Pulls_CB_%s_%s_%s",fitPartName[j].c_str(),fitvarnameCB[k].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                               );
                  h_PartPulls_TAPS[i][j][k] = hfPullsTAPSKF->makeTH1D(Form("TAPS: %s pulls of fitted %s after %s",fitvarnameTA[k].c_str(),fitPartName[j].c_str(),cuts[i+nrCuts_beforeKF].c_str()), //title
                                                                   Form("%s pull",fitvarnameTA[k].c_str()),"#", // xlabel, ylabel
                                                                   BinSettings(200,-10.,10.),   // our binnings
                                                                   Form("h_Pulls_TAPS_%s_%s_%s",fitPartName[j].c_str(),fitvarnameTA[k].c_str(),cuts[i+nrCuts_beforeKF].c_str()), true    // ROOT object name, auto-generated if omitted
                                                                   );

            }

        }

    }


    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("tree"));

}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event); 

    const auto& data = event.Reconstructed();
    const auto& candidates = data.Candidates;

    //set fitter uncertainty model:
    //fitter.SetUncertaintyModel(fit_model);

    // choose uncertainty depending on Data/MC input
    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);
    fitter.SetUncertaintyModel(is_MC ? fit_model_mc : fit_model_data);

    h_nClusters->Fill(data.Clusters.size());

    //Fill e.g. the polar angle distribution into a histogram
    TCandidatePtrList all;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    TCandidatePtrList protonCand;

    TParticleList photons;
    TParticleList protons;

    vector<bool> chaCanInCB;
    vector<bool> chaCanInTAPS;
    vector<bool> chaCanInPID;
    vector<bool> neuCanInCB;
    vector<bool> neuCanInTAPS;
    vector<bool> protonCanInCB;
    vector<bool> protonCanInTAPS;
    vector<bool> protonCanInPID;

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
    vector<double> protonCanCluSize;
    vector<double> protonCanCaloE;
    vector<double> protonCanVetoE;
    vector<double> protonCanTime;
    vector<double> protonThe;
    vector<double> protonPhi;

    for(const auto& cand : candidates.get_iter()) {
        //h_VetoEnergies->Fill(cand->VetoEnergy);
        all.emplace_back(cand);

        if(cand->VetoEnergy <= vetoEthreshold) {
            photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, cand));
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

            if(cand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                chaCanInPID.push_back(true);
            }
            else{
                chaCanInPID.push_back(false);
            }

            //if(cand->Theta*radtodeg < 50){

            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                protonCanInCB.push_back(true);
            }
            else{
                protonCanInCB.push_back(false);
            }
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                protonCanInTAPS.push_back(true);
            }
            else{
                protonCanInTAPS.push_back(false);
            }

            if(cand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                protonCanInPID.push_back(true);
            }
            else{
                protonCanInPID.push_back(false);
            }

            protonCand.emplace_back(cand);
            protons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, cand));
            protonThe.push_back(cand->Theta*radtodeg);
            protonPhi.push_back(cand->Phi*radtodeg);
            protonCanCluSize.push_back(cand->ClusterSize);
            protonCanCaloE.push_back(cand->CaloEnergy);
            protonCanVetoE.push_back(cand->VetoEnergy);
            protonCanTime.push_back(cand->Time);
            //}
        }

    }

    //getting access to the pi0 decay photons

    double wpi0gOnly3g_IM=-100;
    TParticle g[neu_nrSel];
    TParticle proton_tmp;
    TParticlePtr proton_toFit;

    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;
    vector<double> wpi0g_times;

    if(photons.size() == neu_nrSel && protons.size() == cha_nrSel){

        proton_toFit = protons.at(0);

        for(int i = 0; i<neu_nrSel; i++){
            g[i] = TParticle(ParticleTypeDatabase::Photon, neutral[i]);
        }

        proton_tmp = TParticle(ParticleTypeDatabase::Proton, protonCand[0]);

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
        wpi0gOnly3g_time = wpi0gOnly3g_time/3.;
        wpi0g_times.push_back(wpi0gOnly3g_time);

    }

    //Looping over the taggerhits

    TParticle omega_tmp;
    double corTaggTime;

    double best_probability = std_ext::NaN;
    double fit_z_vert = -50;
    double fitbeamE = -50;

    int cut_ind;

    for (const auto& taggerhit : data.TaggerHits) { // Event loop

        corTaggTime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        promptrandom.SetTaggerTime(corTaggTime);

        const double weight = promptrandom.FillWeight();

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);

        //Adding selections & filling histograms

        cut_ind=0;

        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        cut_ind++;

        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(!triggersimu.HasTriggered())
            continue;

        cut_ind++;

        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(!(photons.size() == neu_nrSel && protons.size() == cha_nrSel))
            continue;

        TLorentzVector L3g;

        TLorentzVector g0 = (TLorentzVector)g[0];
        TLorentzVector g1 = (TLorentzVector)g[1];
        TLorentzVector g2 = (TLorentzVector)g[2];

        for(int i = 0; i<neu_nrSel; i++){
            L3g += (TLorentzVector)(g[i]);
        }

        omega_tmp = TParticle(ParticleTypeDatabase::Omega,(TLorentzVector)(L3g));

        TLorentzVector Lw_tmp = (TLorentzVector)(omega_tmp);
        TLorentzVector Lp_tmp = (TLorentzVector)(proton_tmp);

        TLorentzVector LmissingProton = Linitial-Lw_tmp;

        TLorentzVector Lw = Lw_tmp;
        TLorentzVector Lp = Lp_tmp;
        TLorentzVector Lw_boosted = Lw;
        TLorentzVector Lp_boosted = Lp;
        Lw_boosted.Boost(-Linitial.BoostVector());
        Lp_boosted.Boost(-Linitial.BoostVector());
        double wp_angle = cos(Lw_boosted.Angle(Lp_boosted.Vect()));

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

        TLorentzVector missing_pi0 = Linitial-Lp_tmp-wg;
        TLorentzVector missing_Omega = Linitial-Lp_tmp;
        TLorentzVector wpi0_boosted = wpi0;
        wpi0_boosted.Boost(-missing_Omega.BoostVector());
        TLorentzVector wg_boosted = wg;
        wg_boosted.Boost(-missing_Omega.BoostVector());
        TLorentzVector pi0g1_boosted = wpi0g[0];
        pi0g1_boosted.Boost(-missing_pi0.BoostVector());
        TLorentzVector pi0g2_boosted = wpi0g[1];
        pi0g2_boosted.Boost(-missing_pi0.BoostVector());

        double pi0gg_angle = cos(pi0g1_boosted.Angle(pi0g2_boosted.Vect()));
        double wpi0g_angle = cos(wpi0_boosted.Angle(wg_boosted.Vect()));

        TLorentzVector clusterElowestg;
        double clusterElowestgEnergy;
        double lowestAngleToSplitoff;

        if(neuCanCaloE[0] <= neuCanCaloE[1] && neuCanCaloE[0] <= neuCanCaloE[2]){
            clusterElowestg = (TLorentzVector)g0;
            clusterElowestgEnergy = neuCanCaloE[0];
            if(clusterElowestg.Angle(g1.Vect()) <= clusterElowestg.Angle(g2.Vect())){
                lowestAngleToSplitoff = clusterElowestg.Angle(g1.Vect());
            }
            else{
                lowestAngleToSplitoff = clusterElowestg.Angle(g2.Vect());
            }
        }
        else if(neuCanCaloE[1] <= neuCanCaloE[0] && neuCanCaloE[1] <= neuCanCaloE[2]){
            clusterElowestg = (TLorentzVector)g1;
            clusterElowestgEnergy = neuCanCaloE[1];
            if(clusterElowestg.Angle(g0.Vect()) <= clusterElowestg.Angle(g2.Vect())){
                lowestAngleToSplitoff = clusterElowestg.Angle(g0.Vect());
            }
            else{
                lowestAngleToSplitoff = clusterElowestg.Angle(g2.Vect());
            }
        }
        else{
            clusterElowestg = (TLorentzVector)g2;
            clusterElowestgEnergy = neuCanCaloE[2];
            if(clusterElowestg.Angle(g0.Vect()) <= clusterElowestg.Angle(g1.Vect())){
                lowestAngleToSplitoff = clusterElowestg.Angle(g0.Vect());
            }
            else{
                lowestAngleToSplitoff = clusterElowestg.Angle(g1.Vect());
            }
        }

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

        //if(!(LmissingProton.M() > (mpFit-1.5*sigmaMissingP) && LmissingProton.M() < (mpFit+1.5*sigmaMissingP)))
        if(!(LmissingProton.M() > (mp-0.2*mp) && LmissingProton.M() < (mp+0.2*mp)))
            continue;

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

        if(!(InitialPhotonVec.E()>=Omega_Ethreshold))
            continue;

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

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

        //----------------------------------------------------------------------------------------------------------------------------------

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

        h_2gPi0_IM[cut_ind-nrCuts_beforePi0]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforePi0]->Fill(omega_tmp.M(),wpi0.E(),weight);

        TLorentzVector missing_pi0 = Linitial-Lp_tmp-wg;
        TLorentzVector missing_Omega = Linitial-Lp_tmp;
        TLorentzVector wpi0_boosted = wpi0;
        wpi0_boosted.Boost(-missing_Omega.BoostVector());
        TLorentzVector wg_boosted = wg;
        wg_boosted.Boost(-missing_Omega.BoostVector());
        TLorentzVector pi0g1_boosted = wpi0g[0];
        pi0g1_boosted.Boost(-missing_pi0.BoostVector());
        TLorentzVector pi0g2_boosted = wpi0g[1];
        pi0g2_boosted.Boost(-missing_pi0.BoostVector());

        double pi0gg_angle = cos(pi0g1_boosted.Angle(pi0g2_boosted.Vect()));
        double wpi0g_angle = cos(wpi0_boosted.Angle(wg_boosted.Vect()));

        h_wpi0g_BackToBack[cut_ind-nrCuts_beforePi0]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforePi0]->Fill(pi0gg_angle,weight);

        */

        /*
        if(!(wpi0.M()>(mpi0Fit-3*sigmaPi0IM) && wpi0.M()<(mpi0Fit+3*sigmaPi0IM)))
            continue;
        */

        if(!(wpi0.M()>(mpi0-0.4*mpi0) && wpi0.M()<(mpi0+0.4*mpi0)))
            continue;

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

        if(!(neuCanCaloE[0] > neuCaloEthreshhold && neuCanCaloE[1] > neuCaloEthreshhold && neuCanCaloE[2] > neuCaloEthreshhold))
            continue;

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
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

        //weight_res += weight;


        //Setting up a kinFit

        std::vector<utils::Fitter::FitParticle> fitParticles;
        vector<LorentzVec> fitPhotons;

        APLCON::Result_t fitresult = fitter.DoFit(taggerhit.PhotonEnergy, proton_toFit, photons);

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

        fitPhotons.clear();
        for(const auto& fph : fitter.GetFittedPhotons())
            fitPhotons.push_back(*fph);

        fitParticles.clear();
        fitParticles = fitter.GetFitParticles();

        fit_z_vert = fitter.GetFittedZVertex();
        vec3 vertshift{0,0,fit_z_vert};
        fitbeamE = fitter.GetFittedBeamE();

        long double fitphotcomb1 = (*fitted_photons.at(0) + *fitted_photons.at(1)).M();
        long double fitphotcomb2 = (*fitted_photons.at(0) + *fitted_photons.at(2)).M();
        long double fitphotcomb3 = (*fitted_photons.at(1) + *fitted_photons.at(2)).M();
        long double fitphotcombdiff[neu_nrSel] = {abs(fitphotcomb1-mpi0),abs(fitphotcomb2-mpi0),abs(fitphotcomb3-mpi0)};

        long double pi02g_fitM;

        if(fitphotcombdiff[0]<fitphotcombdiff[1] && fitphotcombdiff[0]<fitphotcombdiff[2]){
            pi02g_fitM = fitphotcomb1;
        }
        else if(fitphotcombdiff[1]<fitphotcombdiff[2] && fitphotcombdiff[1]<fitphotcombdiff[0]){
            pi02g_fitM = fitphotcomb2;
        }
        else{
            pi02g_fitM = fitphotcomb3;
        }

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,weight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(fit_z_vert,weight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(fitbeamE,weight);

        h_IM3g_Fit[cut_ind-nrCuts_beforeKF]->Fill((*fitted_photons.at(0) + *fitted_photons.at(1) + *fitted_photons.at(2)).M(),weight);
        h_IM2gPi0_Fit[cut_ind-nrCuts_beforeKF]->Fill(pi02g_fitM,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

        h_p_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(Lp.Theta()*radtodeg,Lp.E(),weight);
        h_w_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(Lw.Theta()*radtodeg,Lw.E(),weight);
        h_wg_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wg.Theta()*radtodeg,wg.E(),weight);
        h_wpi0_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wpi0.Theta()*radtodeg,wpi0.E(),weight);

        for(int i=0; i<nr_pi0g;i++){
            h_wpi02g_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wpi0g[i].Theta()*radtodeg,wpi0g[i].E(),weight);
        }

        //Pulls:

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){

            h_Ek_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Ek(),fitted_proton->Ek()-proton_toFit->Ek(),weight);
            h_Theta_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Theta()*radtodeg,(fitted_proton->p + vertshift).Theta()*radtodeg - proton_toFit->Theta()*radtodeg,weight);
            h_Phi_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Phi()*radtodeg,fitted_proton->Phi()*radtodeg-proton_toFit->Phi()*radtodeg,weight);
            for (unsigned int j=0; j<nrFitVars; j++){
                h_PartPulls_CB[cut_ind-nrCuts_beforeKF][0][j]->Fill(fitParticles.at(0).GetPulls().at(j),weight);
            }

        }
        else{

            h_Ek_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Ek(),fitted_proton->Ek()-proton_toFit->Ek(),weight);
            h_Theta_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Theta()*radtodeg,(fitted_proton->p + vertshift).Theta()*radtodeg - proton_toFit->Theta()*radtodeg,weight);
            h_Phi_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Phi()*radtodeg,fitted_proton->Phi()*radtodeg-proton_toFit->Phi()*radtodeg,weight);
            for (unsigned int j=0; j<nrFitVars; j++){
                h_PartPulls_TAPS[cut_ind-nrCuts_beforeKF][0][j]->Fill(fitParticles.at(0).GetPulls().at(j),weight);
            }
        }

        for (unsigned int i=0; i<photons.size(); i++){

            if(photons.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){

                h_Ek_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Ek(),fitted_photons.at(i)->Ek()-photons.at(i)->Ek(),weight);
                h_Theta_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Theta()*radtodeg,(fitted_photons.at(i)->p + vertshift).Theta()*radtodeg - photons.at(i)->Theta()*radtodeg,weight);
                h_Phi_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Phi()*radtodeg,fitted_photons.at(i)->Phi()*radtodeg-photons.at(i)->Phi()*radtodeg,weight);
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_CB[cut_ind-nrCuts_beforeKF][1][j]->Fill(fitParticles.at(i+1).GetPulls().at(j),weight);
                }

            }
            else{

                h_Ek_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Ek(),fitted_photons.at(i)->Ek()-photons.at(i)->Ek(),weight);
                h_Theta_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Theta()*radtodeg,(fitted_photons.at(i)->p + vertshift).Theta()*radtodeg - photons.at(i)->Theta()*radtodeg,weight);
                h_Phi_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Phi()*radtodeg,fitted_photons.at(i)->Phi()*radtodeg-photons.at(i)->Phi()*radtodeg,weight);
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_TAPS[cut_ind-nrCuts_beforeKF][1][j]->Fill(fitParticles.at(i+1).GetPulls().at(j),weight);
                }
            }

        }

        if(!(best_probability>0.01))
            continue;

        cut_ind++;
        stat[cut_ind]+=weight;
        h_RecData_Stat->Fill(cut_ind, weight);

        h_TaggerTime[cut_ind]->Fill(taggerhit.Time, weight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wp_angle,weight);
        h_wpi0g_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(wpi0g_angle,weight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeSel]->Fill(pi0gg_angle,weight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,weight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(fit_z_vert,weight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(fitbeamE,weight);

        h_IM3g_Fit[cut_ind-nrCuts_beforeKF]->Fill((*fitted_photons.at(0) + *fitted_photons.at(1) + *fitted_photons.at(2)).M(),weight);
        h_IM2gPi0_Fit[cut_ind-nrCuts_beforeKF]->Fill(pi02g_fitM,weight);

        h_2gPi0_IM[cut_ind-nrCuts_beforeSel]->Fill(wpi0.M(),weight);
        h_EpivsIM3g[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),wpi0.E(),weight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(LmissingProton.M(),weight);
        h_3g_IM[cut_ind-nrCuts_beforeSel]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[cut_ind-nrCuts_beforeSel]->Fill((g[1]+g[2]).M(),weight);

        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g1.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g0.Angle(g2.Vect())*radtodeg,weight);
        h_2gComb_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(g1.Angle(g2.Vect())*radtodeg,weight);

        h_neuLowestCaloE[cut_ind-nrCuts_beforeSel]->Fill(clusterElowestgEnergy,weight);
        h_2gLowestClusterE_OpeningAngles[cut_ind-nrCuts_beforeSel]->Fill(lowestAngleToSplitoff*radtodeg,weight);

        h_doublyDCScm_gp_wp[cut_ind-nrCuts_beforeSel]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);

        h_nCandidates[cut_ind]->Fill(candidates.size(), weight);
        h_nNeuChaCandidates[cut_ind]->Fill(charged.size(), neutral.size(), weight);

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),weight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[cut_ind]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuCaloE[cut_ind]->Fill(neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaCB[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_CB[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
            if(neuCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(neuCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuEkinVSThetaTAPS[cut_ind]->Fill(neuThe[i],neuCanCaloE[i],weight);
                h_neuEvsT_TAPS[cut_ind]->Fill(neuCanCaloE[i],neuCanTime[i],weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[cut_ind]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaCaloE[cut_ind]->Fill(chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllVetoE_CB[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaCB[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_CB[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
            if(chaCanInTAPS[i]){
                h_AllVetoE_TAPS[cut_ind]->Fill(chaCanVetoE[i],weight);
                h_AllCaloEvsVetoE_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaEkinVSThetaTAPS[cut_ind]->Fill(chaThe[i],chaCanCaloE[i],weight);
                h_chaEvsT_TAPS[cut_ind]->Fill(chaCanCaloE[i],chaCanTime[i],weight);
            }
        }

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->CaloEnergy,proton_toFit->Candidate->VetoEnergy,weight);
        }
        if(proton_toFit->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeSel]->Fill(proton_toFit->Candidate->FindVetoCluster()->Time,proton_toFit->Candidate->FindVetoCluster()->Energy,weight);
        }

        h_p_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(Lp.Theta()*radtodeg,Lp.E(),weight);
        h_w_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(Lw.Theta()*radtodeg,Lw.E(),weight);
        h_wg_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wg.Theta()*radtodeg,wg.E(),weight);
        h_wpi0_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wpi0.Theta()*radtodeg,wpi0.E(),weight);

        for(int i=0; i<nr_pi0g;i++){
            h_wpi02g_EvTheta[cut_ind-nrCuts_beforeKF]->Fill(wpi0g[i].Theta()*radtodeg,wpi0g[i].E(),weight);
        }

        //Pulls:

        if(proton_toFit->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){

            h_Ek_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Ek(),fitted_proton->Ek()-proton_toFit->Ek(),weight);
            h_Theta_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Theta()*radtodeg,(fitted_proton->p + vertshift).Theta()*radtodeg - proton_toFit->Theta()*radtodeg,weight);
            h_Phi_dev_CB[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Phi()*radtodeg,fitted_proton->Phi()*radtodeg-proton_toFit->Phi()*radtodeg,weight);
            for (unsigned int j=0; j<nrFitVars; j++){
                h_PartPulls_CB[cut_ind-nrCuts_beforeKF][0][j]->Fill(fitParticles.at(0).GetPulls().at(j),weight);
            }

        }
        else{

            h_Ek_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Ek(),fitted_proton->Ek()-proton_toFit->Ek(),weight);
            h_Theta_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Theta()*radtodeg,(fitted_proton->p + vertshift).Theta()*radtodeg - proton_toFit->Theta()*radtodeg,weight);
            h_Phi_dev_TAPS[cut_ind-nrCuts_beforeKF][0]->Fill(proton_toFit->Phi()*radtodeg,fitted_proton->Phi()*radtodeg-proton_toFit->Phi()*radtodeg,weight);
            for (unsigned int j=0; j<nrFitVars; j++){
                h_PartPulls_TAPS[cut_ind-nrCuts_beforeKF][0][j]->Fill(fitParticles.at(0).GetPulls().at(j),weight);
            }
        }

        for (unsigned int i=0; i<photons.size(); i++){

            if(photons.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){

                h_Ek_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Ek(),fitted_photons.at(i)->Ek()-photons.at(i)->Ek(),weight);
                h_Theta_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Theta()*radtodeg,(fitted_photons.at(i)->p + vertshift).Theta()*radtodeg - photons.at(i)->Theta()*radtodeg,weight);
                h_Phi_dev_CB[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Phi()*radtodeg,fitted_photons.at(i)->Phi()*radtodeg-photons.at(i)->Phi()*radtodeg,weight);
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_CB[cut_ind-nrCuts_beforeKF][1][j]->Fill(fitParticles.at(i+1).GetPulls().at(j),weight);
                }

            }
            else{

                h_Ek_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Ek(),fitted_photons.at(i)->Ek()-photons.at(i)->Ek(),weight);
                h_Theta_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Theta()*radtodeg,(fitted_photons.at(i)->p + vertshift).Theta()*radtodeg - photons.at(i)->Theta()*radtodeg,weight);
                h_Phi_dev_TAPS[cut_ind-nrCuts_beforeKF][1]->Fill(photons.at(i)->Phi()*radtodeg,fitted_photons.at(i)->Phi()*radtodeg-photons.at(i)->Phi()*radtodeg,weight);
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_TAPS[cut_ind-nrCuts_beforeKF][1][j]->Fill(fitParticles.at(i+1).GetPulls().at(j),weight);
                }
            }

        }

        //Filling the tree for further analysis

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = data.Clusters.size();
        t.PhotonAzimuthAngles = PhotonPhi;
        t.PhotonPolarAngles = PhotonThe;
        t.PhotonPolarAnglesCB = PhotonTheCB;
        t.PhotonPolarAnglesTAPS = PhotonTheTAPS;
        t.Tree->Fill();     
    }

}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::ShowResult()
{    

    cout << "" << endl;

    for (unsigned int i=0; i<nrCuts_total; i++){
        cout << "Amount events after " << i << " applied cuts: " << h_RecData_Stat->GetBinContent(i*steps+1) << endl;
        //h_RecData_Stat->SetBinContent(i*steps+1, stat[i]);
        //h_RecData_Stat->SetBinError(i*steps+1, sqrt(stat[i]));
        h_RecData_Stat->GetXaxis()->SetBinLabel(i*steps+1,cuts[i].c_str());
    }

    cout << "" << endl;

    for (unsigned int i=0; i<(nrCuts_total-1); i++){
        cout << "Relative amount of events after " << (i+1) << " applied cuts: " << 100*(h_RecData_Stat->GetBinContent((i+1)*steps+1)/h_RecData_Stat->GetBinContent(1)) << " %" << endl;
        //h_RecData_relStat->SetBinContent((i+1)*steps+1, stat[i+1]/stat[0]);
        h_RecData_relStat->SetBinContent((i+1)*steps+1, (h_RecData_Stat->GetBinContent((i+1)*steps+1)/h_RecData_Stat->GetBinContent(1)));
        h_RecData_relStat->SetBinError((i+1)*steps+1, sqrt(pow((sqrt(h_RecData_Stat->GetBinError((i+1)*steps+1))/h_RecData_Stat->GetBinContent(1)),2)+pow((sqrt(h_RecData_Stat->GetBinError(1))*h_RecData_Stat->GetBinContent((i+1)*steps+1))/(pow(h_RecData_Stat->GetBinContent(1),2)),2)));
        h_RecData_relStat->GetXaxis()->SetBinLabel((i+1)*steps+1,cuts[i+1].c_str());
    }

    cout << "" << endl;

    ant::canvas(GetName()+": Event-Statistics after applied cuts:")
            << drawoption("HIST")
            << h_RecData_Stat
            << h_RecData_relStat
            << endc; // actually draws the canvas

    /*

    ant::canvas c1(GetName()+": IM(missingP)");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c1 << h_missingP_IM[i];
    }
            c1 << endc; // actually draws the canvas

    ant::canvas c2(GetName()+": IM(3neu)");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c2 << h_3g_IM[i];
    }
            c2 << endc; // actually draws the canvas


    ant::canvas c3(GetName()+": IM(2neu) combinations");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c3 << h_2gComb_IM[i];
    }
            c3 << endc; // actually draws the canvas

    ant::canvas c4(GetName()+": gp_wp cross section in cm-frame");
            c4 << drawoption("Surf");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c4 << h_doublyDCScm_gp_wp[i];
    }
            c4 << endc; // actually draws the canvas

    ant::canvas c5(GetName()+": IM(2neuPi0)");
    for (unsigned int i=0; i<nrCuts_pi0; i++){
            c5 << h_2gPi0_IM[i];
    }
            c5 << endc; // actually draws the canvas

    ant::canvas c6(GetName()+": Deposited Calo- VS Veto-energies in CB");
            c6 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c6 << h_AllCaloEvsVetoE_CB[i];
    }
            c6 << endc; // actually draws the canvas

    ant::canvas c7(GetName()+": Deposited Calo- VS Veto-energies in TAPS");
            c7 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c7 << h_AllCaloEvsVetoE_TAPS[i];
    }
            c7 << endc; // actually draws the canvas

    ant::canvas c8(GetName()+": Initial photonbeam energy");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c8 << h_beamE[i];
    }
            c8 << endc; // actually draws the canvas

*/

/*

    ant::canvas c9(GetName()+": Photons Ekin vs Theta");
            c9 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c9 << h_neuEkinVSTheta[i];
    }
            c9 << endc; // actually draws the canvas

    ant::canvas c10(GetName()+": charged Ekin vs Theta");
            c10 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c10 << h_chaEkinVSTheta[i];
    }
            c10 << endc; // actually draws the canvas

    ant::canvas c11(GetName()+": Photons Ekin vs Phi");
            c11 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c11 << h_neuEkinVSPhi[i];
    }
            c11 << endc; // actually draws the canvas

    ant::canvas c12(GetName()+": charged Ekin vs Phi");
            c12 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c12 << h_chaEkinVSPhi[i];
    }
            c12 << endc; // actually draws the canvas



    ant::canvas c17(GetName()+": charged time-diff to cor. Tagg time in CB");
            c17 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c17 << h_chaTimeDiffCorTaggCB[i];
    }
            c17 << endc; // actually draws the canvas

    ant::canvas c18(GetName()+": charged time-diff to cor. Tagg time in TAPS");
            c18 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c18 << h_chaTimeDiffCorTaggTAPS[i];
    }
            c18 << endc; // actually draws the canvas

    ant::canvas c19(GetName()+": neutrals time-diff to cor. Tagg time in CB");
            c19 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c19 << h_neuTimeDiffCorTaggCB[i];
    }
            c19 << endc; // actually draws the canvas

    ant::canvas c20(GetName()+": neutrals time-diff to cor. Tagg time in TAPS");
            c20 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c20 << h_neuTimeDiffCorTaggTAPS[i];
    }
            c20 << endc; // actually draws the canvas


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

    for (unsigned int i=0; i<nrCuts_pi0; i++){

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

    ant::canvas c_KinFit_overview(GetName()+": Kinfit overview");
            for (unsigned int i=0; i<nrCuts_KF; i++){
            c_KinFit_overview << h_Fit_zvert[i];
            c_KinFit_overview << h_Probability[i];
            c_KinFit_overview << h_fitEbeam[i];
            }
            c_KinFit_overview << endc;

    ant::canvas c_KinFit_invmass(GetName()+": KinFit invmasses");
            for (unsigned int i=0; i<nrCuts_KF; i++){
            c_KinFit_invmass << h_IM3g_Fit[i];
            c_KinFit_invmass << h_IM2gPi0_Fit[i];
            }
            c_KinFit_invmass << endc; // actually draws the canvas


    ant::canvas c_KinFit_CB(GetName()+"Rec. vs fitted particles in CB");
            c_KinFit_CB << drawoption("pcolz");
            for (unsigned int i=0; i<nrPartType; i++){
                c_KinFit_CB << h_Ek_dev_CB[i];
                c_KinFit_CB << h_Theta_dev_CB[i];
                c_KinFit_CB << h_Phi_dev_CB[i];
            }
            c_KinFit_CB << endc; // actually draws the canvas

    ant::canvas c_KinFit_TAPS(GetName()+"Rec. vs fitted particles in TAPS");
            c_KinFit_TAPS << drawoption("pcolz");
            for (unsigned int i=0; i<nrPartType; i++){
                c_KinFit_TAPS << h_Ek_dev_TAPS[i];
                c_KinFit_TAPS << h_Theta_dev_TAPS[i];
                c_KinFit_TAPS << h_Phi_dev_TAPS[i];
            }
            c_KinFit_TAPS << endc; // actually draws the canvas

    ant::canvas c_kfPulls_CB(GetName()+"KinFit CB pulls");
            for (unsigned int i=0; i<nrPartType; i++){
                for (unsigned int j=0; j<nrFitVars; j++){
                    c_kfPulls_CB << h_PartPulls_CB[i][j];
                }
            }
            c_kfPulls_CB << endc; // actually draws the canvas

    ant::canvas c_kfPulls_TAPS(GetName()+"KinFit TAPS pulls");
            for (unsigned int i=0; i<nrPartType; i++){
                for (unsigned int j=0; j<nrFitVars; j++){
                    c_kfPulls_TAPS << h_PartPulls_TAPS[i][j];
                }
            }
            c_kfPulls_TAPS << endc; // actually draws the canvas

    ant::canvas c_CBEsum(GetName()+": CB Esum");
            //c_CBEsum << drawoption("HIST");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c_CBEsum << h_CBEsum[i];
    }
            c_CBEsum << endc; // actually draws the canvas

    */

}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::Finish()
{
    cout << "please work!" << endl;

    cout << "" << endl;

    LOG(INFO) << "Fit Model Statistics Data:\n" << *fit_model_data;
    LOG(INFO) << "Fit Model Statistics MC:\n" << *fit_model_mc;

    cout << "" << endl;

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

AUTO_REGISTER_PHYSICS(scratch_damaurer_pw_ppi0g_p3g_kinFit)
