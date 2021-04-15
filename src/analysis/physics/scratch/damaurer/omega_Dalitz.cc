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
    fitter(nullptr, opts->Get<bool>("FitZVertex", true)),
    fitter_freeZ(nullptr, opts->Get<bool>("FitZVertex", true))
{

    // set standard fitter with constrained Zvertex
    fitter.SetZVertexSigma(3.0);

    // get target information
    const auto target = ExpConfig::Setup::Get().GetTargetProperties();

    // set sigma to 0 for unmeasured --> free z vertex
    fitter_freeZ.SetZVertexSigma(0);
    fitter_freeZ.SetTarget(target.length);  // double expected, target length in cm

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    const BinSettings statistic_bins(number_of_bins,0,nrCuts_total);
    const BinSettings pid_bins((eePID)*20,0,eePID);
    const BinSettings bins_tagger_time(2000, -200, 200);
    const BinSettings bins_nClusters(20);
    const BinSettings bins_nCand(10);
    const BinSettings bins_nNeuCand(6);
    const BinSettings bins_nChaCand(6);
    const BinSettings bins_kfFails(5);
    const BinSettings bins_Veto_Energy(500, 0, 12);
    const BinSettings bins_Calo_Energy(500, 0, 1200);
    const BinSettings BeamE_bins(100,0, 1600);
    const BinSettings zVert_bins(50,-15,15);
    const BinSettings zVertFree_bins(50,-20,20);
    const BinSettings Im_pi0_bins(200, 0, 1000);
    const BinSettings Im_proton_bins(200, 0, 1800);
    const BinSettings Im_omega_bins(200, 0, 1200);
    const BinSettings IM_ee_bins(200, 0, 800);
    const BinSettings kf_prob_bins(1000,0.,1.);
    const BinSettings CB_Esum_bins(250,0.,2000.);
    const BinSettings timeDiffCorTaggCB_bins(50, 0,10);
    const BinSettings timeDiffCorTaggTAPS_bins(100, 0,25);
    const BinSettings PIDtime_bins(241,-60.5,60.5);
    const BinSettings effR_bins(100, 0,25);
    const BinSettings nCrystal_bins(30);
    const BinSettings cos_bins(100,-1,1);
    const BinSettings theta_bins(360, 0, 180);
    const BinSettings phi_bins(720, -180, 180);
    const BinSettings PIDelement_bins(30);

    auto hf_RecData_CandStat = new HistogramFactory("hf_RecData_CandStat", HistFac, "");
    auto hf_Tagger = new HistogramFactory("hf_Tagger", HistFac, "");
    auto hf_CaloEvsVetoE = new HistogramFactory("hf_CaloEvsVetoE", HistFac, "");
    auto hf_Candidates = new HistogramFactory("hf_Candidates", HistFac, "");
    auto hf_VetoE = new HistogramFactory("hf_VetoE", HistFac, "");
    auto hf_BeamE = new HistogramFactory("hf_BeamE", HistFac, "");
    auto hf_2g_IM = new HistogramFactory("hf_2g_IM", HistFac, "");
    auto hf_2gee_IM = new HistogramFactory("hf_2gee_IM", HistFac, "");
    auto hf_missingP_IM = new HistogramFactory("hf_missingP_IM", HistFac, "");
    auto hf_KFfails = new HistogramFactory("hf_KFfails", HistFac, "");
    auto hf_OverviewKF = new HistogramFactory("hf_OverviewKF", HistFac, "");
    auto hf_invmassKF = new HistogramFactory("hf_invmassKF", HistFac, "");
    auto hf_CBEsum = new HistogramFactory("hf_CB_Esum", HistFac, "");

    auto hfPullsCBKF = new HistogramFactory("KinFit_CBpulls", HistFac, "");
    auto hfPullsTAPSKF = new HistogramFactory("KinFit_TAPSpulls", HistFac, "");

    auto hf_OverviewKF_freeZ = new HistogramFactory("hf_OverviewKF_freeZ", HistFac, "");
    auto hf_cluster_effRvsCaloE = new HistogramFactory("hf_cluster_effRvsCaloE", HistFac, "");
    auto hf_cluster_nCrystalsvsCaloE = new HistogramFactory("hf_cluster_nCrystalsvsCaloE", HistFac, "");

    //Side-check hf
    auto hf_CaloEvsVetoE_CB_SideCheck = new HistogramFactory("hf_CaloEvsVetoE_CB_SideCheck", HistFac, "");
    auto hf_CaloEvsVetoE_TAPS_SideCheck = new HistogramFactory("hf_CaloEvsVetoE_TAPS_SideCheck", HistFac, "");
    auto hf_TimeDiffCorTagg_CB_SideCheck = new HistogramFactory(" hf_TimeDiffCorTagg_CB_SideCheck", HistFac, "");
    auto hf_TimeDiffCorTagg_TAPS_SideCheck = new HistogramFactory(" hf_TimeDiffCorTagg_TAPS_SideCheck", HistFac, "");

    auto hf_BackToBack = new HistogramFactory("hf_BackToBack", HistFac, "");

    auto hf_eeChecks = new HistogramFactory("hf_eeChecks", HistFac, "");
    auto hf_PIDEvsTime = new HistogramFactory("hf_PIDEvsTime", HistFac, "");
    auto hf_EvsThetaCB = new HistogramFactory("hf_EvsThetaCB", HistFac, "");
    auto hf_EvsThetaTAPS = new HistogramFactory("hf_EvsThetaTAPS", HistFac, "");
    auto hf_EvsPhi = new HistogramFactory("hf_EvsPhi", HistFac, "");
    auto hf_PID_stuff = new HistogramFactory("hf_PID_stuff", HistFac, "");

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = hf_Tagger->makeTH1D("Number of clusters", "nClusters", "#", bins_nClusters, "h_nClusters", true);
    h_mmpFails = hf_KFfails->makeTH1D("Fit mm(p) fails","nrFails","#",bins_kfFails,"h_mmpFails",true);
    h_kfFails = hf_KFfails->makeTH1D("Fit status fails","nrFails","#",bins_kfFails,"h_kfFails",true);
    h_totalFails = hf_KFfails->makeTH1D("Total fit fails","nrFails","#",bins_kfFails,"h_totalFails",true);

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

    h_RecData_relStat = hf_RecData_CandStat->makeTH1D("rel. amount of events after cuts:",     // title
                                      "Cuts", "rel. #",     // xlabel, ylabel
                                      statistic_bins,  // our binnings
                                      "h_RecData_relStat", true    // ROOT object name, auto-generated if omitted
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

    h_CBEsum[i] = hf_CBEsum->makeTH1D("CB Esum "+cuts[i],"CB Esum [MeV]","#",CB_Esum_bins,"h_CBEsum_"+cuts[i],true);

    h_beamE[i] = hf_BeamE->makeTH1D("Initial photon beam energy "+cuts[i],     // title
                                                         "E_{photonbeam} [MeV]", "#",     // xlabel, ylabel
                                                         BeamE_bins,  // our binnings
                                                         "h_beamE_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                         );

    }

    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){

    h_2g_IM[i] = hf_2g_IM->makeTH1D("IM(2g) "+cuts[i+nrCuts_beforeSel],     // title
                                         "IM(2g) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2g_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                         );

    h_missingP_IM[i] = hf_missingP_IM->makeTH1D("Im(missingP) "+cuts[i+nrCuts_beforeSel],     // title
                                                         "IM(missingP) [MeV]", "#",     // xlabel, ylabel
                                                         Im_proton_bins,  // our binnings
                                                         "h_missingP_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                                         );

    h_2gee_IM[i] = hf_2gee_IM->makeTH1D("IM(2gee) "+cuts[i+nrCuts_beforeSel],     // title
                                        "IM(2gee) [MeV]", "#",     // xlabel, ylabel
                                        Im_omega_bins,  // our binnings
                                        "h_2gee_IM_"+cuts[i+nrCuts_beforeSel], true     // ROOT object name, auto-generated if omitted
                                        );

    }

    for (unsigned int i=0; i<nrCuts_KF; i++){
        h_Probability[i] = hf_OverviewKF->makeTH1D(Form("Fit probability %s",cuts_KF[i].c_str()),"P(#chi^{2})","#",kf_prob_bins,Form("h_Probability_%s",cuts_KF[i].c_str()),true);
        h_Fit_zvert[i] = hf_OverviewKF->makeTH1D(Form("Fitted z-vertex %s",cuts_KF[i].c_str()),"z [cm]","#",zVert_bins,Form("h_Fit_zvert_%s",cuts_KF[i].c_str()),true);
        h_fitEbeam[i] = hf_OverviewKF->makeTH1D(Form("Fitted beam energy %s",cuts_KF[i].c_str()),"E [MeV]","#",BeamE_bins,Form("h_fitEbeam_%s",cuts_KF[i].c_str()),true);
        h_IMmissingP_Fit[i] = hf_invmassKF->makeTH1D(Form("Fitted mm(p) %s",cuts_KF[i].c_str()),"mm(fitted_p) [MeV]","#",Im_proton_bins,Form("h_IMmissingP_Fit_%s",cuts_KF[i].c_str()),true);
        h_IM2gee_Fit[i] = hf_invmassKF->makeTH1D(Form("Fitted 2gee invariant mass %s",cuts_KF[i].c_str()),"Im(fitted_2gee) [MeV]","#",Im_omega_bins,Form("h_IM2gee_Fit_%s",cuts_KF[i].c_str()),true);
        h_IM2g_Fit[i] = hf_invmassKF->makeTH1D(Form("Fitted 2g invariant mass %s",cuts_KF[i].c_str()),"Im(fitted_2g) [MeV]", "#",Im_pi0_bins,Form("h_IM2g_Fit_%s",cuts_KF[i].c_str()), true);
        h_Probability_freeZ[i] = hf_OverviewKF_freeZ->makeTH1D(Form("Free Z Kinfitter probability %s",cuts_KF[i].c_str()),"P(#chi^{2})","#",kf_prob_bins,Form("h_Probability_freeZ_%s",cuts_KF[i].c_str()),true);
        h_Fit_zvert_freeZ[i] = hf_OverviewKF_freeZ->makeTH1D(Form("Fitted unconstrained z-vertex %s",cuts_KF[i].c_str()),"z [cm]","#",zVertFree_bins,Form("h_Fit_zvert_unconstrained_%s",cuts_KF[i].c_str()), true);

        //-------------------PID stuff--------------------------------------:

        h_IMmissingP_Fit_samePID[i] = hf_PID_stuff->makeTH1D(Form("Fitted mm(p) with ee same PID %s",cuts_KF[i].c_str()),"mm(fitted_p) [MeV]","#",Im_proton_bins,Form("h_IMmissingP_Fit_samePID_%s",cuts_KF[i].c_str()),true);
        h_IM2gee_Fit_samePID[i] = hf_PID_stuff->makeTH1D(Form("Fitted 2gee invariant mass with ee same PID %s",cuts_KF[i].c_str()),"Im(fitted_2gee) [MeV]","#",Im_omega_bins,Form("h_IM2gee_Fit_samePID_%s",cuts_KF[i].c_str()),true);
        h_IM2g_Fit_samePID[i] = hf_PID_stuff->makeTH1D(Form("Fitted 2g invariant mass with ee same PID %s",cuts_KF[i].c_str()),"Im(fitted_2g) [MeV]", "#",Im_pi0_bins,Form("h_IM2g_Fit_samePID_%s",cuts_KF[i].c_str()), true);
        h_IMmissingP_Fit_diffPID[i] = hf_PID_stuff->makeTH1D(Form("Fitted mm(p) with ee diff PID %s",cuts_KF[i].c_str()),"mm(fitted_p) [MeV]","#",Im_proton_bins,Form("h_IMmissingP_Fit_diffPID_%s",cuts_KF[i].c_str()),true);
        h_IM2gee_Fit_diffPID[i] = hf_PID_stuff->makeTH1D(Form("Fitted 2gee invariant mass with ee diff PID %s",cuts_KF[i].c_str()),"Im(fitted_2gee) [MeV]","#",Im_omega_bins,Form("h_IM2gee_Fit_diffPID_%s",cuts_KF[i].c_str()),true);
        h_IM2g_Fit_diffPID[i] = hf_PID_stuff->makeTH1D(Form("Fitted 2g invariant mass with ee diff PID %s",cuts_KF[i].c_str()),"Im(fitted_2g) [MeV]", "#",Im_pi0_bins,Form("h_IM2g_Fit_diffPID_%s",cuts_KF[i].c_str()), true);

        //-------------------------------------------------------------------

        h_cluster_effRvsCaloE[i] = hf_cluster_effRvsCaloE->makeTH2D("No Protons: Cluster energy vs effective radius "+cuts_KF[i],     // title
                                         "Cluster E [MeV]", "effR",  // xlabel, ylabel
                                         bins_Calo_Energy, effR_bins,    // our binnings
                                         "h_cluster_effRvsCaloE_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_cluster_nCrystalsvsCaloE[i] = hf_cluster_nCrystalsvsCaloE->makeTH2D("No Protons: : Cluster energy vs number of crystals "+cuts_KF[i],     // title
                                         "Cluster E [MeV]", "nCrystals",  // xlabel, ylabel
                                         bins_Calo_Energy, nCrystal_bins,    // our binnings
                                         "h_cluster_nCrystalsvsCaloE_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        //Some side check hists:

        h_NoProton_CaloEvsVetoE_CB[i] = hf_CaloEvsVetoE_CB_SideCheck->makeTH2D("No protons: Deposited calo vs veto energies in CB "+cuts_KF[i],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                         "h_NoProton_CaloEvsVetoE_CB_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_NoProton_CaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE_TAPS_SideCheck->makeTH2D("No protons: Deposited calo vs veto energies in TAPS "+cuts_KF[i],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                         "h_NoProton_CaloEvsVetoE_TAPS_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_Proton_CaloEvsVetoE_CB[i] = hf_CaloEvsVetoE_CB_SideCheck->makeTH2D("Protons: Deposited calo vs veto energies in CB "+cuts_KF[i],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                         "h_Proton_CaloEvsVetoE_CB_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_Proton_CaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE_TAPS_SideCheck->makeTH2D("Protons: Deposited calo vs veto energies in TAPS "+cuts_KF[i],     // title
                                         "Calo E [MeV]", "Veto E [MeV]",  // xlabel, ylabel
                                         bins_Calo_Energy, bins_Veto_Energy,    // our binnings
                                         "h_Proton_CaloEvsVetoE_TAPS_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_Proton_PIDEvsTime[i] = hf_PIDEvsTime->makeTH2D("Protons: PID energy vs time "+cuts_KF[i],     // title
                                                                          "PID time [ns]", "PID E [MeV]",  // xlabel, ylabel
                                                                          PIDtime_bins,BinSettings(200,0.,20),    // our binnings
                                                                          "h_Proton_PIDEvsTime_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_NoProton_PIDEvsTime[i] = hf_PIDEvsTime->makeTH2D("No Protons: PID energy vs time "+cuts_KF[i],     // title
                                                                          "PID time [ns]", "PID E [MeV]",  // xlabel, ylabel
                                                                          PIDtime_bins,BinSettings(200,0.,20),    // our binnings
                                                                          "h_NoProton_PIDEvsTime_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Photon_EvsThetaCB[i] = hf_EvsThetaCB->makeTH2D("Photons in CB: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Photon_EvsThetaCB_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_NoProton_EvsThetaCB[i] = hf_EvsThetaCB->makeTH2D("No Protons in CB: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_NoProton_EvsThetaCB_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Proton_EvsThetaCB[i] = hf_EvsThetaCB->makeTH2D("Protons in CB: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Proton_EvsThetaCB_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Photon_EvsThetaTAPS[i] = hf_EvsThetaTAPS->makeTH2D("Photons in TAPS: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Photon_EvsThetaTAPS_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_NoProton_EvsThetaTAPS[i] = hf_EvsThetaTAPS->makeTH2D("No Protons in TAPS: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_NoProton_EvsThetaTAPS_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Proton_EvsThetaTAPS[i] = hf_EvsThetaTAPS->makeTH2D("Protons in TAPS: Energy vs theta "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          theta_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Proton_EvsThetaTAPS_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );
        h_Photon_EvsPhi[i] = hf_EvsPhi->makeTH2D("Photons: Energy vs phi "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          phi_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Photon_EvsPhi_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_NoProton_EvsPhi[i] = hf_EvsPhi->makeTH2D("No Protons: Energy vs phi "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          phi_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_NoProton_EvsPhi_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );

        h_Proton_EvsPhi[i] = hf_EvsPhi->makeTH2D("Protons: Energy vs phi "+cuts_KF[i],     // title
                                                                          "#Theta [deg]", "Calo E [MeV]",  // xlabel, ylabel
                                                                          phi_bins,bins_Calo_Energy,    // our binnings
                                                                          "h_Proton_EvsPhi_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                                                          );


        h_Proton_TimeDiffCorTaggCB[i] = hf_TimeDiffCorTagg_CB_SideCheck->makeTH1D("Protons: Time diff to cor tagg-time in CB "+cuts_KF[i],     // title
                                                             "(t_{CB}-t_{Tagger}) [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggCB_bins,  // our binnings
                                                             "h_Proton_TimeDiffCorTaggCB_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_Proton_TimeDiffCorTaggTAPS[i] = hf_TimeDiffCorTagg_TAPS_SideCheck->makeTH1D("Protons: Time diff to cor tagg-time in TAPS "+cuts_KF[i],     // title
                                                             "t_{TAPS}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggTAPS_bins,  // our binnings
                                                             "h_Proton_TimeDiffCorTaggTAPS_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_NoProton_TimeDiffCorTaggCB[i] = hf_TimeDiffCorTagg_CB_SideCheck->makeTH1D("No protons: Time diff to cor tagg-time in CB "+cuts_KF[i],     // title
                                                             "(t_{CB}-t_{Tagger}) [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggCB_bins,  // our binnings
                                                             "h_NoProton_TimeDiffCorTaggCB_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_NoProton_TimeDiffCorTaggTAPS[i] = hf_TimeDiffCorTagg_TAPS_SideCheck->makeTH1D("No protons: Time diff to cor tagg-time in TAPS "+cuts_KF[i],     // title
                                                             "t_{TAPS}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggTAPS_bins,  // our binnings
                                                             "h_NoProton_TimeDiffCorTaggTAPS_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_Photon_TimeDiffCorTaggCB[i] = hf_TimeDiffCorTagg_CB_SideCheck->makeTH1D("Photons: Time diff to cor tagg-time in CB "+cuts_KF[i],     // title
                                                             "(t_{CB}-t_{Tagger}) [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggCB_bins,  // our binnings
                                                             "h_Photon_TimeDiffCorTaggCB_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_Photon_TimeDiffCorTaggTAPS[i] = hf_TimeDiffCorTagg_TAPS_SideCheck->makeTH1D("Photons: Time diff to cor tagg-time in TAPS "+cuts_KF[i],     // title
                                                             "t_{TAPS}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggTAPS_bins,  // our binnings
                                                             "h_Photon_TimeDiffCorTaggTAPS_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_IMee[i] = hf_eeChecks->makeTH1D("IM(ee) "+cuts_KF[i],     // title
                                                             "IM(ee) [MeV]", "#",     // xlabel, ylabel
                                                             IM_ee_bins,  // our binnings
                                                             "h_IMee_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_IMeeFit[i] = hf_eeChecks->makeTH1D("IM(fitted_ee) "+cuts_KF[i],     // title
                                                             "IM(fitted_ee) [MeV]", "#",     // xlabel, ylabel
                                                             IM_ee_bins,  // our binnings
                                                             "h_IMeeFit_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_IMeeFit_samePID[i] = hf_PID_stuff->makeTH1D("IM(fitted_ee) with ee same PID "+cuts_KF[i],     // title
                                                             "IM(fitted_ee) [MeV]", "#",     // xlabel, ylabel
                                                             IM_ee_bins,  // our binnings
                                                             "h_IMeeFit_samePID_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_IMeeFit_diffPID[i] = hf_PID_stuff->makeTH1D("IM(fitted_ee) with ee diff PID "+cuts_KF[i],     // title
                                                             "IM(fitted_ee) [MeV]", "#",     // xlabel, ylabel
                                                             IM_ee_bins,  // our binnings
                                                             "h_IMeeFit_diffPID_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_wp_BackToBack[i] = hf_BackToBack->makeTH1D("Omega Proton Back-To-Back check "+cuts_KF[i],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_wp_BackToBack_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                         );

        h_dilee_BackToBack[i] = hf_BackToBack->makeTH1D("Dilepton el pos Back-To-Back check "+cuts_KF[i],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_dilee_BackToBack_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                         );

        h_pi0gg_BackToBack[i] = hf_BackToBack->makeTH1D("Pi0 gg Back-To-Back check "+cuts_KF[i],     // title
                                         "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                         cos_bins,  // our binnings
                                         "h_pi0gg_BackToBack_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                         );

        h_wpi0dil_BackToBack[i] = hf_BackToBack->makeTH1D("Pi0 dilepton Back-To-Back check "+cuts_KF[i],     // title
                                                          "cos(#Theta_{cms})", "#",     // xlabel, ylabel
                                                          cos_bins,  // our binnings
                                                          "h_wpi0dil_BackToBack_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                          );


        h_eePID[i] = hf_PID_stuff->makeTH1D("el/pos PID statistics "+cuts_KF[i],     // title
                                                          "PID element", "#",     // xlabel, ylabel
                                                          pid_bins,  // our binnings
                                                          "h_eePID_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                          );

        h_eePID[i]->GetXaxis()->SetBinLabel(1+1*steps_PID,PIDstat[0].c_str());
        h_eePID[i]->GetXaxis()->SetBinLabel(1+2*steps_PID,PIDstat[1].c_str());

        h_eePIDelementNumbers[i] = hf_PID_stuff->makeTH2D("PID element numbers of leptons "+cuts_KF[i],     // title
                                         "l1 PID nr", "l2 PID nr",  // xlabel, ylabel
                                         PIDelement_bins, PIDelement_bins,    // our binnings
                                         "h_eePIDelementNumbers_"+cuts_KF[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_eeOpeningAngles[i] = hf_eeChecks->makeTH1D("el/pos opening angles "+cuts_KF[i],     // title
                                                          "#Theta [deg]", "#",     // xlabel, ylabel
                                                          theta_bins,  // our binnings
                                                          "h_eeOpeningAngles_"+cuts_KF[i], true     // ROOT object name, auto-generated if omitted
                                                          );

    }

    for (unsigned int i=0; i<nrPartType; i++){ //0 refers to proton, 1 to photons

        for (unsigned int j=0; j<nrFitVars; j++){

              h_PartPulls_CB[i][j] = hfPullsCBKF->makeTH1D(Form("CB: %s pulls of fitted %s",fitvarnameCB[j].c_str(),fitPartName[i].c_str()), //title
                                                           Form("%s pull",fitvarnameCB[j].c_str()),"#", // xlabel, ylabel
                                                           BinSettings(200,-10.,10.),   // our binnings
                                                           Form("h_Pulls_CB_%s_%s",fitPartName[i].c_str(),fitvarnameCB[j].c_str()), true    // ROOT object name, auto-generated if omitted
                                                           );
              h_PartPulls_TAPS[i][j] = hfPullsTAPSKF->makeTH1D(Form("TAPS: %s pulls of fitted %s",fitvarnameTA[j].c_str(),fitPartName[i].c_str()), //title
                                                               Form("%s pull",fitvarnameTA[j].c_str()),"#", // xlabel, ylabel
                                                               BinSettings(200,-10.,10.),   // our binnings
                                                               Form("h_Pulls_TAPS_%s_%s",fitPartName[i].c_str(),fitvarnameCB[j].c_str()), true    // ROOT object name, auto-generated if omitted
                                                               );

        }
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

    const auto& data = event.Reconstructed();
    const auto& clusters = data.Clusters;
    const auto& candidates = data.Candidates;

    // set fitter uncertainty models
    //fitter.SetUncertaintyModel(fit_model);
    //fitter_freeZ.SetUncertaintyModel(fit_model);

    // choose uncertainty depending on Data/MC input
    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);
    fitter.SetUncertaintyModel(is_MC ? fit_model_mc : fit_model_data);
    fitter_freeZ.SetUncertaintyModel(is_MC ? fit_model_mc : fit_model_data);

    //Set up particle combinations for the kinFit
    //utils::ProtonPhotonCombs proton_photons(candidates);
    //particle_combs_t protphotcombs = proton_photons();

    //Fill e.g. the polar angle distribution into a histogram
    TCandidatePtrList all;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    //TCandidatePtrList protonCand;

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

    TParticleList photons;
    TParticleList photonCombs[nrCombs];
    TParticleList protonCombs[nrCombs];
    TCandidatePtrList photonCandForFit[nrCombs];
    TCandidatePtrList protonCandForFit[nrCombs];

    if(neutral.size() == neu_nrSel && charged.size() == cha_nrSel){

        //proton_toFit = protons.at(0);

        for(int i = 0; i<neu_nrSel; i++){
            g[i] = TParticle(ParticleTypeDatabase::Photon, neutral[i]);
            L2g += (TLorentzVector)(g[i]);
            photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, neutral.at(i)));

            for(int j = 0; j<nrCombs; j++){
                photonCombs[j].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, neutral.at(i)));
                photonCandForFit[j].emplace_back(neutral.at(i));
            }

        }

        for(int i = 0; i<nrCombs; i++){
            protonCombs[i].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, charged.at(i)));
            protonCandForFit[i].emplace_back(charged.at(i));
        }

        //fill in the rest combinations:

        //photonCombs[0].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::eMinus, charged.at(1)));
        //photonCombs[0].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::ePlus, charged.at(2)));
        photonCombs[0].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(1)));
        photonCombs[0].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(2)));
        photonCandForFit[0].emplace_back(charged.at(1));
        photonCandForFit[0].emplace_back(charged.at(2));

        //photonCombs[1].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::eMinus, charged.at(0)));
        //photonCombs[1].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::ePlus, charged.at(2)));
        photonCombs[1].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(0)));
        photonCombs[1].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(2)));
        photonCandForFit[1].emplace_back(charged.at(0));
        photonCandForFit[1].emplace_back(charged.at(2));

        //photonCombs[2].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::eMinus, charged.at(0)));
        //photonCombs[2].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::ePlus, charged.at(1)));
        photonCombs[2].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(0)));
        photonCombs[2].emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, charged.at(1)));
        photonCandForFit[2].emplace_back(charged.at(0));
        photonCandForFit[2].emplace_back(charged.at(1));

    }

    h_nClusters->Fill(clusters.size());

    double best_probability;
    double bestFitted_Zvert;
    double bestFitted_BeamE;  
    double bestFitted_mm_proton;
    double bestFitted_mm_2gee;
    double bestFitted_mm_2g;
    double bestFitted_mm_2e;
    double bestFitted_Iterations;
    double best_mm_proton;
    double best_mm_2gee;

    int nr_mmpFails;
    int nr_kfFails;

    double best_probability_freeZ;
    double bestFitted_Zvert_freeZ;
    double bestFitted_Iterations_freeZ;

    /*
    double bestFitted_BeamE_freeZ;
    double bestFitted_mm_proton_freeZ;
    double bestFitted_mm_2gee_freeZ;
    double bestFitted_mm_2g_freeZ;
    double bestFitted_mm_2e_freeZ;
    double best_mm_proton_freeZ;
    double best_mm_2gee_freeZ;
    */

    int nr_mmpFails_freeZ;
    int nr_kfFails_freeZ;

    double mm_protonCombs[nrCombs];
    double mm_2geeCombs[nrCombs];
    double mm_proton;

    TLorentzVector L4g;
    TLorentzVector L4gCombs[nrCombs];

    TParticlePtr bestFitted_proton;
    TParticleList bestFitted_photons;

    TParticlePtr bestFitted_proton_freeZ;
    TParticleList bestFitted_photons_freeZ;

    int cut_ind;

    for(auto& taggerhit : data.TaggerHits) {

        //Get corrected Taggertime
        double cortagtime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        promptrandom.SetTaggerTime(cortagtime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;    

        const double TaggWeight = promptrandom.FillWeight();

        cut_ind=0;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        if(!triggersimu.HasTriggered())
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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


        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        if(!(neutral.size() == neu_nrSel && charged.size() == cha_nrSel))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);

        for(int i = 0; i<nrCombs; i++){
            L4gCombs[i] = (TLorentzVector)(*photonCombs[i].at(0)+*photonCombs[i].at(1)+*photonCombs[i].at(2)+*photonCombs[i].at(3));
            mm_2geeCombs[i] = L4gCombs[i].M();
            mm_protonCombs[i] = (Linitial - L4gCombs[i]).M();
            h_2gee_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_2geeCombs[i],TaggWeight);
            h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_protonCombs[i],TaggWeight);
        }

        if(!(InitialPhotonVec.E()>=Omega_Ethreshold))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);

        for(int i = 0; i<nrCombs; i++){
            h_2gee_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_2geeCombs[i],TaggWeight);
            h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_protonCombs[i],TaggWeight);
        }

        if(!(L2g.M()>(mpi0-0.4*mpi0) && L2g.M()<(mpi0+0.4*mpi0)))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        //just a test:

        /*
        TLorentzVector L2g;
        for (unsigned int i=0; i<neu_nrSel; i++){
            L2g += (TLorentzVector)(*photons.at(i));
            //L2g += (TLorentzVector)(*photonCombs[0].at(i));
        }
        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);
        */

        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);
        for(int i = 0; i<nrCombs; i++){
            h_2gee_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_2geeCombs[i],TaggWeight);
            h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(mm_protonCombs[i],TaggWeight);
        }

        //FreeZ kinfitter for experimenting:

        nr_mmpFails_freeZ = 0;
        nr_kfFails_freeZ = 0;

        best_probability_freeZ = std_ext::NaN;
        bestFitted_Zvert_freeZ = std_ext::NaN;
        bestFitted_Iterations_freeZ = std_ext::NaN;

        /*
        bestFitted_BeamE_freeZ = std_ext::NaN;
        bestFitted_mm_proton_freeZ = std_ext::NaN;
        bestFitted_mm_2gee_freeZ = std_ext::NaN;
        bestFitted_mm_2g_freeZ = std_ext::NaN;
        best_mm_proton_freeZ = std_ext::NaN;
        best_mm_2gee_freeZ = std_ext::NaN;
        */

        std::vector<utils::Fitter::FitParticle> bestFitParticles_freeZ;
        bestFitParticles_freeZ.clear();
        bestKFindex_freeZ = 0;

        for(int i = 0; i<nrCombs; i++){

            L4g = (TLorentzVector)(*photonCombs[i].at(0)+*photonCombs[i].at(1)+*photonCombs[i].at(2)+*photonCombs[i].at(3));
            mm_proton = (Linitial-L4g).M();

            if(!(mm_proton > (mp-0.2*mp) && mm_proton < (mp+0.2*mp))){
                nr_mmpFails_freeZ += 1;
                continue;
            }

            //Performing the kinFit
            APLCON::Result_t fitresult_freeZ = fitter_freeZ.DoFit(taggerhit.PhotonEnergy, protonCombs[i].at(0), photonCombs[i]);

            // check if the fit converged
            if (fitresult_freeZ.Status != APLCON::Result_Status_t::Success){
                nr_kfFails_freeZ += 1;
                continue;
            }

            // check if we found a better probability for this fit and copy it if true, continue otherwise
            if (!std_ext::copy_if_greater(best_probability_freeZ, fitresult_freeZ.Probability))
                continue;

            bestKFindex_freeZ = i;

            // retrieve the fitted photon and proton information as well as the number of iterations
            bestFitted_proton_freeZ = fitter_freeZ.GetFittedProton();
            bestFitted_photons_freeZ = fitter_freeZ.GetFittedPhotons();
            bestFitted_Zvert_freeZ = fitter_freeZ.GetFittedZVertex();
            bestFitted_Iterations_freeZ = fitresult_freeZ.NIterations;

            bestFitParticles_freeZ.clear();
            bestFitParticles_freeZ = fitter_freeZ.GetFitParticles();

            /*
            bestFitted_BeamE_freeZ = fitter_freeZ.GetFittedBeamE();
            best_mm_proton_freeZ = mm_proton;
            best_mm_2gee_freeZ = L4g.M();

            bestFitted_mm_2gee_freeZ = ((TLorentzVector)(*bestFitted_photons_freeZ.at(0) + *bestFitted_photons_freeZ.at(1) + *bestFitted_photons_freeZ.at(2) + *bestFitted_photons_freeZ.at(3))).M();
            bestFitted_mm_2g_freeZ = ((TLorentzVector)(*bestFitted_photons_freeZ.at(0) + *bestFitted_photons_freeZ.at(1))).M();
            bestFitted_mm_proton_freeZ = ((TLorentzVector)(LorentzVec({0, 0, bestFitted_BeamE_freeZ}, bestFitted_BeamE_freeZ) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass()))-(TLorentzVector)(*bestFitted_photons_freeZ.at(0) + *bestFitted_photons_freeZ.at(1) + *bestFitted_photons_freeZ.at(2) + *bestFitted_photons_freeZ.at(3))).M();
            */

        }

        /*
        if(((nr_mmpFails_freeZ + nr_kfFails_freeZ) != nrCombs) && best_probability_freeZ > 0.01){
            h_Probability_freeZ->Fill(best_probability_freeZ,TaggWeight);
            h_Fit_zvert_freeZ->Fill(bestFitted_Zvert_freeZ,TaggWeight);
        }
        */

        //Standard & main kinfitter for this analysis:

        nr_mmpFails = 0;
        nr_kfFails = 0;

        best_probability = std_ext::NaN;
        bestFitted_Zvert = std_ext::NaN;
        bestFitted_BeamE = std_ext::NaN;
        bestFitted_mm_proton = std_ext::NaN;
        bestFitted_mm_2gee = std_ext::NaN;
        bestFitted_mm_2g = std_ext::NaN;
        bestFitted_mm_2e = std_ext::NaN;
        bestFitted_Iterations = std_ext::NaN;
        best_mm_proton = std_ext::NaN;
        best_mm_2gee = std_ext::NaN;

        std::vector<utils::Fitter::FitParticle> bestFitParticles;
        bestFitParticles.clear();
        bestKFindex = 0;

        for(int i = 0; i<nrCombs; i++){

            L4g = (TLorentzVector)(*photonCombs[i].at(0)+*photonCombs[i].at(1)+*photonCombs[i].at(2)+*photonCombs[i].at(3));
            mm_proton = (Linitial-L4g).M();

            if(!(mm_proton > (mp-0.2*mp) && mm_proton < (mp+0.2*mp))){
                nr_mmpFails += 1;
                continue;
            }

            //Performing the kinFit
            APLCON::Result_t fitresult = fitter.DoFit(taggerhit.PhotonEnergy, protonCombs[i].at(0), photonCombs[i]);

            // check if the fit converged
            if (fitresult.Status != APLCON::Result_Status_t::Success){
                nr_kfFails += 1;
                continue;
            }

            // check if we found a better probability for this fit and copy it if true, continue otherwise
            if (!std_ext::copy_if_greater(best_probability, fitresult.Probability))
                continue;

            bestKFindex = i;

            // retrieve the fitted photon and proton information as well as the number of iterations
            bestFitted_proton = fitter.GetFittedProton();
            bestFitted_photons = fitter.GetFittedPhotons();
            bestFitted_Zvert = fitter.GetFittedZVertex();
            bestFitted_BeamE = fitter.GetFittedBeamE();
            bestFitted_Iterations = fitresult.NIterations;

            best_mm_proton = mm_proton;
            best_mm_2gee = L4g.M();

            bestFitted_mm_2gee = ((TLorentzVector)(*bestFitted_photons.at(0) + *bestFitted_photons.at(1) + *bestFitted_photons.at(2) + *bestFitted_photons.at(3))).M();
            bestFitted_mm_2g = ((TLorentzVector)(*bestFitted_photons.at(0) + *bestFitted_photons.at(1))).M();
            bestFitted_mm_proton = ((TLorentzVector)(LorentzVec({0, 0, bestFitted_BeamE}, bestFitted_BeamE) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass()))-(TLorentzVector)(*bestFitted_photons.at(0) + *bestFitted_photons.at(1) + *bestFitted_photons.at(2) + *bestFitted_photons.at(3))).M();
            bestFitted_mm_2e = ((TLorentzVector)(*bestFitted_photons.at(2) + *bestFitted_photons.at(3))).M();

            bestFitParticles.clear();
            bestFitParticles = fitter.GetFitParticles();
        }

        h_mmpFails->Fill(nr_mmpFails,TaggWeight);
        h_kfFails->Fill(nr_kfFails,TaggWeight);
        h_totalFails->Fill(nr_mmpFails+nr_kfFails,TaggWeight);

        if(!((nr_mmpFails + nr_kfFails) != nrCombs && (nr_mmpFails_freeZ + nr_kfFails_freeZ) != nrCombs))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

        TLorentzVector Lproton_tmp = (TLorentzVector)(*protonCombs[bestKFindex].at(0));
        TLorentzVector Lomega_tmp = (TLorentzVector)(*photonCombs[bestKFindex].at(0)+*photonCombs[bestKFindex].at(1)+*photonCombs[bestKFindex].at(2)+*photonCombs[bestKFindex].at(3));
        TLorentzVector Lpion_tmp = (TLorentzVector)(*photonCombs[bestKFindex].at(0)+*photonCombs[bestKFindex].at(1));
        TLorentzVector Ldil_tmp = (TLorentzVector)(*photonCombs[bestKFindex].at(2)+*photonCombs[bestKFindex].at(3));

        TLorentzVector missing_omega = Linitial-Lproton_tmp;
        TLorentzVector missing_pi0 = Linitial-Lproton_tmp-Ldil_tmp;
        TLorentzVector pi0_g1 = (TLorentzVector)(*photonCombs[bestKFindex].at(0));
        TLorentzVector pi0_g2 = (TLorentzVector)(*photonCombs[bestKFindex].at(1));

        TLorentzVector missing_dilepton = Linitial-Lproton_tmp-Lpion_tmp;
        TLorentzVector dil_l1 = (TLorentzVector)(*photonCombs[bestKFindex].at(2));
        TLorentzVector dil_l2 = (TLorentzVector)(*photonCombs[bestKFindex].at(3));

        TLorentzVector Lp_boosted = Lproton_tmp;
        TLorentzVector Lw_boosted = Lomega_tmp;
        Lp_boosted.Boost(-Linitial.BoostVector());
        Lw_boosted.Boost(-Linitial.BoostVector());
        double wp_angle = cos(Lp_boosted.Angle(Lw_boosted.Vect()));

        TLorentzVector Lpion_boosted = Lpion_tmp;
        TLorentzVector Ldil_boosted = Ldil_tmp;
        Lpion_boosted.Boost(-missing_omega.BoostVector());
        Ldil_boosted.Boost(-missing_omega.BoostVector());
        double piondil_angle = cos(Lpion_boosted.Angle(Ldil_boosted.Vect()));

        TLorentzVector Lpi0g1_boosted = pi0_g1;
        TLorentzVector Lpi0g2_boosted = pi0_g2;
        Lpi0g1_boosted.Boost(-missing_pi0.BoostVector());
        Lpi0g2_boosted.Boost(-missing_pi0.BoostVector());
        double piongg_angle = cos(Lpi0g1_boosted.Angle(Lpi0g2_boosted.Vect()));

        TLorentzVector Ldil_l1_boosted = dil_l1;
        TLorentzVector Ldil_l2_boosted = dil_l2;
        Ldil_l1_boosted.Boost(-missing_dilepton.BoostVector());
        Ldil_l2_boosted.Boost(-missing_dilepton.BoostVector());
        double dilee_angle = cos(Ldil_l1_boosted.Angle(Ldil_l2_boosted.Vect()));

        double ee_openingAngle = dil_l1.Angle(dil_l2.Vect())*radtodeg;

        bool eeInPID = false;
        bool eeSamePID = false;
        int l1eenr = 0;
        int l2eenr = 0;

        if(photonCombs[bestKFindex].at(2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID && photonCombs[bestKFindex].at(3)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            eeInPID = true;
            l1eenr = photonCombs[bestKFindex].at(2)->Candidate->FindVetoCluster()->CentralElement;
            l2eenr = photonCombs[bestKFindex].at(3)->Candidate->FindVetoCluster()->CentralElement;
            if(l1eenr == l2eenr){
                eeSamePID = true;
            }
        }

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-nrCuts_beforeSel]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }


        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //Pulls:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            for (unsigned int i=0; i<nrFitVars; i++){
                h_PartPulls_CB[0][i]->Fill(bestFitParticles.at(0).GetPulls().at(i),TaggWeight);
            }
        }
        else{
            for (unsigned int i=0; i<nrFitVars; i++){
                h_PartPulls_TAPS[0][i]->Fill(bestFitParticles.at(0).GetPulls().at(i),TaggWeight);
            }
        }

        for (unsigned int i=0; i<nrPhotons; i++){

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_CB[1][j]->Fill(bestFitParticles.at(i+1).GetPulls().at(j),TaggWeight);
                }
            }
            else{
                for (unsigned int j=0; j<nrFitVars; j++){
                    h_PartPulls_TAPS[1][j]->Fill(bestFitParticles.at(i+1).GetPulls().at(j),TaggWeight);
                }
            }

        }

        //effR calculation for charged cand (not proton!):
        // CaloCluster in the candidate is the CB, get all hits in the cluster --> crystals
        // this method doesn't make sense for clusters with only one or two contributing crystals
        double effR_l1 = std_ext::NaN;
        double e_l1 = std_ext::NaN;
        double effR_l2 = std_ext::NaN;
        double e_l2 = std_ext::NaN;
        double nCrystals_l1 = std_ext::NaN;
        double nCrystals_l2 = std_ext::NaN;

        TClusterHitList crystals_l1 = photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Hits;
        TClusterHitList crystals_l2 = photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Hits;
        // get a (x,y,z) vector to the CB position of the central crystal of this candidate
        vec3 central_l1 = cb_detector->GetPosition(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->CentralElement);
        vec3 central_l2 = cb_detector->GetPosition(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->CentralElement);

        nCrystals_l1 = crystals_l1.size();
        nCrystals_l2 = crystals_l2.size();
        //double nCrystals_l1 = photonCombs[bestKFindex].at(2)->Candidate->Clusters.size();
        //double nCrystals_l2 = photonCombs[bestKFindex].at(3)->Candidate->Clusters.size();

        if (nCrystals_l1 > 2){
            effR_l1 = 0;
            e_l1 = 0;
            // loop over all crystals
            for (TClusterHit crystal_l1 : crystals_l1) {
                // r is the angle (distance) between the current crystal and the center crystal of the cluster
                const double r = std_ext::radian_to_degree(central_l1.Angle(cb_detector->GetPosition(crystal_l1.Channel)));
                effR_l1 += r*r*crystal_l1.Energy;
                e_l1 += crystal_l1.Energy;
            }
        }

        effR_l1 = sqrt(effR_l1/e_l1);

        if (nCrystals_l2 > 2){
            effR_l2 = 0;
            e_l2 = 0;
            // loop over all crystals
            for (TClusterHit crystal_l2 : crystals_l2) {
                // r is the angle (distance) between the current crystal and the center crystal of the cluster
                const double r = std_ext::radian_to_degree(central_l2.Angle(cb_detector->GetPosition(crystal_l2.Channel)));
                effR_l2 += r*r*crystal_l2.Energy;
                e_l2 += crystal_l2.Energy;
            }
        }

        effR_l2 = sqrt(effR_l2/e_l2);  // normalize energy weighted distance of the crystals to the total energy

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);

        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        //----------------------------------------------------------------------------------------

        if(!(best_probability > bestprob_cutval && best_probability_freeZ > bestprob_cutval))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-nrCuts_beforeSel]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-nrCuts_beforeSel]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-nrCuts_beforeSel]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }

        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);

        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        //-----------------------------------------------------------------------------------

        if(!(bestFitted_Zvert_freeZ > -10 && bestFitted_Zvert_freeZ < 5))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-2]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-2]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-2]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }


        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);

        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        //-----------------------------------------------------------------------------------

        Double_t clusterEl1 = photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy;
        Double_t clusterEl2 = photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy;

        Double_t x1_effR = 75;
        Double_t y1_effR = 10;
        Double_t x2_effR = 650;
        Double_t y2_effR = 6.7;
        Double_t shift_effR = 0;

        Double_t m_effR = (y2_effR-y1_effR)/(x2_effR-x1_effR);
        Double_t b_effR = y2_effR-(y2_effR-y1_effR)/(x2_effR-x1_effR)*x2_effR;

        if((isfinite(effR_l1) && clusterEl1 > 75 && clusterEl1 < 950 && effR_l1 >= (m_effR*clusterEl1+b_effR+shift_effR)) || (isfinite(effR_l2) && clusterEl2 > 75 && clusterEl2 < 950 && effR_l2 >= (m_effR*clusterEl2+b_effR+shift_effR)))
            continue;

        Double_t x1_effR2 = 75;
        Double_t y1_effR2 = 1;
        Double_t x2_effR2 = 950;
        Double_t y2_effR2 = 3.5;
        Double_t shift_effR2 = 0.5;

        Double_t m_effR2 = (y2_effR2-y1_effR2)/(x2_effR2-x1_effR2);
        Double_t b_effR2 = y2_effR2-(y2_effR2-y1_effR2)/(x2_effR2-x1_effR2)*x2_effR2;

        if((isfinite(effR_l1) && clusterEl1 > 75 && clusterEl1 < 950 && effR_l1 <= (m_effR2*clusterEl1+b_effR2+shift_effR2)) || (isfinite(effR_l2) && clusterEl2 > 75 && clusterEl2 < 950 && effR_l2 <= (m_effR2*clusterEl2+b_effR2+shift_effR2)))
            continue;

        //box cut for testing stuff:
        /*
        if((isfinite(effR_l1) && clusterEl1 > 80 && clusterEl1 < 700 && effR_l1 > 7) || (isfinite(effR_l2) && clusterEl2 > 80 && clusterEl2 < 700 && effR_l2 > 7))
            continue;
        */

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-2]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-2]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-2]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }

        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);


        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        //-----------------------------------------------------------------------------------

        Double_t x1_nCryst_first = 80;
        Double_t y1_nCryst_first = 1;
        Double_t x2_nCryst_first = 700;
        Double_t y2_nCryst_first = 10.5;
        Double_t shift_nCryst_first = 0;

        Double_t m_nCryst_first = (y2_nCryst_first-y1_nCryst_first)/(x2_nCryst_first-x1_nCryst_first);
        Double_t b_nCryst_first = y2_nCryst_first-(y2_nCryst_first-y1_nCryst_first)/(x2_nCryst_first-x1_nCryst_first)*x2_nCryst_first;

        if((isfinite(nCrystals_l1) && clusterEl1 < 900 && nCrystals_l1 <= (m_nCryst_first*clusterEl1+b_nCryst_first+shift_nCryst_first)) || (isfinite(nCrystals_l2) && clusterEl2 < 900 && nCrystals_l2 <= (m_nCryst_first*clusterEl2+b_nCryst_first+shift_nCryst_first)))
            continue;

        Double_t x1_nCryst_second = 0;
        Double_t y1_nCryst_second = 9;
        Double_t x2_nCryst_second = 600;
        Double_t y2_nCryst_second = 19;
        Double_t shift_nCryst_second = 0;

        Double_t m_nCryst_second = (y2_nCryst_second-y1_nCryst_second)/(x2_nCryst_second-x1_nCryst_second);
        Double_t b_nCryst_second = y2_nCryst_second-(y2_nCryst_second-y1_nCryst_second)/(x2_nCryst_second-x1_nCryst_second)*x2_nCryst_second;

        if((isfinite(nCrystals_l1) && clusterEl1 < 600 && nCrystals_l1 >= (m_nCryst_second*clusterEl1+b_nCryst_second+shift_nCryst_second)) || (isfinite(nCrystals_l2) && clusterEl2 < 600 && nCrystals_l2 >= (m_nCryst_second*clusterEl2+b_nCryst_second+shift_nCryst_second)))
            continue;

        if((isfinite(nCrystals_l1) && clusterEl1 >= 600 && nCrystals_l1 >= 19) || (isfinite(nCrystals_l2) && clusterEl2 >= 600 && nCrystals_l2 >= 19))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-2]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-2]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-2]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }

        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);


        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        //-----------------------------------------------------------------------------------

        if(!(photonCombs[bestKFindex].at(2)->Candidate->VetoEnergy > lepton_VetoE_LOWERthreshold && photonCombs[bestKFindex].at(2)->Candidate->VetoEnergy < lepton_VetoE_UPPERthreshold && photonCombs[bestKFindex].at(3)->Candidate->VetoEnergy > lepton_VetoE_LOWERthreshold && photonCombs[bestKFindex].at(3)->Candidate->VetoEnergy < lepton_VetoE_UPPERthreshold))
            continue;

        cut_ind++;
        stat[cut_ind]+=TaggWeight;
        h_RecData_Stat->Fill(cut_ind, TaggWeight);

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

        h_CBEsum[cut_ind]->Fill(triggersimu.GetCBEnergySum(),TaggWeight);
        h_beamE[cut_ind]->Fill(InitialPhotonVec.E(),TaggWeight);

        h_2g_IM[cut_ind-2]->Fill(L2g.M(),TaggWeight);

        h_missingP_IM[cut_ind-2]->Fill(best_mm_proton,TaggWeight);
        h_2gee_IM[cut_ind-2]->Fill(best_mm_2gee,TaggWeight);

        h_Probability[cut_ind-nrCuts_beforeKF]->Fill(best_probability,TaggWeight);
        h_Fit_zvert[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert,TaggWeight);
        h_fitEbeam[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_BeamE,TaggWeight);
        h_IM2gee_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
        h_IM2g_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
        h_IMmissingP_Fit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
        h_Probability_freeZ[cut_ind-nrCuts_beforeKF]->Fill(best_probability_freeZ,TaggWeight);
        h_Fit_zvert_freeZ[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_Zvert_freeZ,TaggWeight);

        h_wp_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(wp_angle,TaggWeight);
        h_wpi0dil_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piondil_angle,TaggWeight);
        h_dilee_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(dilee_angle,TaggWeight);
        h_pi0gg_BackToBack[cut_ind-nrCuts_beforeKF]->Fill(piongg_angle,TaggWeight);

        h_eeOpeningAngles[cut_ind-nrCuts_beforeKF]->Fill(ee_openingAngle,TaggWeight);

        //PID stuff:

        if(eeInPID){
            h_eePIDelementNumbers[cut_ind-nrCuts_beforeKF]->Fill(l1eenr,l2eenr,TaggWeight);
        }

        if(eeSamePID){
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(2,TaggWeight);
            h_IM2gee_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_samePID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        }
        else{
            h_eePID[cut_ind-nrCuts_beforeKF]->Fill(1,TaggWeight);
            h_IM2gee_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2gee,TaggWeight);
            h_IM2g_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2g,TaggWeight);
            h_IMmissingP_Fit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_proton,TaggWeight);
            h_IMeeFit_diffPID[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);
        }

        //PID E vs time:

        if(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
            h_Proton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Time,protonCombs[bestKFindex].at(0)->Candidate->FindVetoCluster()->Energy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_NoProton_PIDEvsTime[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Time,photonCombs[bestKFindex].at(i+2)->Candidate->FindVetoCluster()->Energy,TaggWeight);
            }

        }

        //-----------------Filling IM(ee) hists----------------------------------------------------------

        h_IMee[cut_ind-nrCuts_beforeKF]->Fill((dil_l1+dil_l2).M(),TaggWeight);
        h_IMeeFit[cut_ind-nrCuts_beforeKF]->Fill(bestFitted_mm_2e,TaggWeight);

        //-----------------Filling cluster-hists of charged particles (not protons!)---------------------

        if (isfinite(effR_l1)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,effR_l1,TaggWeight);
        }
        if (isfinite(effR_l2)){
            h_cluster_effRvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,effR_l2,TaggWeight);
        }
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(2)->Candidate->FindCaloCluster()->Energy,nCrystals_l1,TaggWeight);
        h_cluster_nCrystalsvsCaloE[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(3)->Candidate->FindCaloCluster()->Energy,nCrystals_l2,TaggWeight);


        //-----------------side-checks & more hists----------------------------------------------------

        h_Proton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Phi)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_Proton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }
        if(protonCombs[bestKFindex].at(0)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            h_Proton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,protonCombs[bestKFindex].at(0)->Candidate->VetoEnergy,TaggWeight);
            h_Proton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Time)-cortagtime,TaggWeight);
            h_Proton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((protonCombs[bestKFindex].at(0)->Candidate->Theta)*radtodeg,protonCombs[bestKFindex].at(0)->Candidate->CaloEnergy,TaggWeight);
        }

        for (unsigned int i=0; i<neu_nrSel; i++){

            h_Photon_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_Photon_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_Photon_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Time)-cortagtime,TaggWeight);
                h_Photon_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i)->Candidate->CaloEnergy,TaggWeight);
            }

            h_NoProton_EvsPhi[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Phi)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);

            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                h_NoProton_CaloEvsVetoE_CB[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaCB[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
            if(photonCombs[bestKFindex].at(i+2)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                h_NoProton_CaloEvsVetoE_TAPS[cut_ind-nrCuts_beforeKF]->Fill(photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,photonCombs[bestKFindex].at(i+2)->Candidate->VetoEnergy,TaggWeight);
                h_NoProton_TimeDiffCorTaggTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Time)-cortagtime,TaggWeight);
                h_NoProton_EvsThetaTAPS[cut_ind-nrCuts_beforeKF]->Fill((photonCombs[bestKFindex].at(i+2)->Candidate->Theta)*radtodeg,photonCombs[bestKFindex].at(i+2)->Candidate->CaloEnergy,TaggWeight);
            }
        }

        t.TaggW = TaggWeight;
        t.nClusters = data.Clusters.size();
        t.Tree->Fill();

    } 

}

void scratch_damaurer_omega_Dalitz::ShowResult()
{

    // ShowResult is called after processing of events has finished,
    // and interactive mode (aka non-batchmode) is chosen

    // ant::canvas nice wrapper around TCanvas

    for (unsigned int i=0; i<nrCuts_total; i++){
        cout << "Amount events after " << i << " applied cuts: " << h_RecData_Stat->GetBinContent(i*steps+1) << endl;
        //h_RecData_Stat->SetBinContent(i*steps+1, stat[i]);
        //h_RecData_Stat->SetBinError(i*steps+1, sqrt(stat[i]));
        h_RecData_Stat->GetXaxis()->SetBinLabel(i*steps+1,cuts[i].c_str());
    }

    for (unsigned int i=0; i<(nrCuts_total-1); i++){
        cout << "Relative amount of events after " << (i+1) << " applied cuts: " << 100*(h_RecData_Stat->GetBinContent((i+1)*steps+1)/h_RecData_Stat->GetBinContent(1)) << " %" << endl;
        //h_RecData_relStat->SetBinContent((i+1)*steps+1, stat[i+1]/stat[0]);
        h_RecData_relStat->SetBinContent((i+1)*steps+1, (h_RecData_Stat->GetBinContent((i+1)*steps+1)/h_RecData_Stat->GetBinContent(1)));
        h_RecData_relStat->SetBinError((i+1)*steps+1, sqrt(pow((sqrt(h_RecData_Stat->GetBinError((i+1)*steps+1))/h_RecData_Stat->GetBinContent(1)),2)+pow((sqrt(h_RecData_Stat->GetBinError(1))*h_RecData_Stat->GetBinContent((i+1)*steps+1))/(pow(h_RecData_Stat->GetBinContent(1),2)),2)));
        h_RecData_relStat->GetXaxis()->SetBinLabel((i+1)*steps+1,cuts[i+1].c_str());
    }

    ant::canvas(GetName()+": Event statistics after applied cuts:")
            << drawoption("HIST")
            << h_RecData_Stat
            << h_RecData_relStat
            << endc; // actually draws the canvas

    /*
    ant::canvas c_Tagger_hists(GetName()+": Some basic hists");
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

    ant::canvas c_CBEsum(GetName()+": CB Esum");
            //c_CBEsum << drawoption("HIST");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c_CBEsum << h_CBEsum[i];
    }
            c_CBEsum << endc; // actually draws the canvas



    ant::canvas c_2g_IM(GetName()+": IM(2g)");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c_2g_IM << h_2g_IM[i];
    }
            c_2g_IM << endc; // actually draws the canvas

    ant::canvas c_mmp(GetName()+": mm(p)");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c_mmp << h_missingP_IM[i];
    }
            c_mmp << endc; // actually draws the canvas

    ant::canvas c_IM_2gee(GetName()+": IM(2gee)");
    for (unsigned int i=0; i<(nrCuts_total-nrCuts_beforeSel); i++){
            c_IM_2gee << h_2gee_IM[i];
    }
            c_IM_2gee << endc; // actually draws the canvas

    ant::canvas c_KF_hists_Fails(GetName()+": Kinfit fail-hists");
            c_KF_hists_Fails << h_mmpFails;
            c_KF_hists_Fails << h_kfFails;
            c_KF_hists_Fails << h_totalFails;
            c_KF_hists_Fails << endc; // actually draws the canvas;

    ant::canvas c_KF_Overview(GetName()+": KinFit overview");
            //c_KF_Overview << drawoption("pcolz");
            for (unsigned int i=0; i<nrCuts_KF; i++){
            c_KF_Overview << h_Probability[i];
            c_KF_Overview << h_Fit_zvert[i];
            c_KF_Overview << h_fitEbeam[i];
            }
            c_KF_Overview << endc; // actually draws the canvas

    ant::canvas c_KF_invMasses(GetName()+": KinFit invariant masses");
            //c_KF_Overview << drawoption("pcolz");
            for (unsigned int i=0; i<nrCuts_KF; i++){
            c_KF_invMasses << h_IM2g_Fit[i];
            c_KF_invMasses << h_IMmissingP_Fit[i];
            c_KF_invMasses << h_IM2gee_Fit[i];
            }
            c_KF_invMasses << endc; // actually draws the canvas

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

    ant::canvas c_KF_Overview_freeZ(GetName()+": KinFit_freeZ overview");
            //c_KF_Overview << drawoption("pcolz");
            for (unsigned int i=0; i<nrCuts_KF; i++){
            c_KF_Overview_freeZ << h_Probability_freeZ[i];
            c_KF_Overview_freeZ << h_Fit_zvert_freeZ[i];
            }
            c_KF_Overview_freeZ << endc; // actually draws the canvas

    */

}

void scratch_damaurer_omega_Dalitz::Finish()
{

    cout << "please work!" << endl;

    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;

    LOG(INFO) << "Fit Model Statistics Data:\n" << *fit_model_data;
    LOG(INFO) << "Fit Model Statistics MC:\n" << *fit_model_mc;

}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)

AUTO_REGISTER_PHYSICS(scratch_damaurer_omega_Dalitz)
