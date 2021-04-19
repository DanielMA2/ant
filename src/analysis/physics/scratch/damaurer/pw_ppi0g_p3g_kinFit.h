// better than using the standard ifndef include guards...
#pragma once

// physics classes need to derive from this interface
#include "physics/Physics.h"
#include "TLorentzVector.h"

// Ant provides many utility classes,
// such as support for prompt-random handling,
// or trigger simulation
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

//Include Kinematic fitter
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/ProtonPhotonCombs.h"

// handle ROOT's TTree
#include "base/WrapTTree.h"

//include more stuff, if necessary:
#include <iostream>

// the physics classes reside in this nested namespace
namespace ant {

namespace expconfig {
namespace detector {
    struct TAPS;
}}
namespace analysis {
namespace physics {

class scratch_damaurer_pw_ppi0g_p3g_kinFit : public Physics {
public:
    scratch_damaurer_pw_ppi0g_p3g_kinFit(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

    struct tree_t : WrapTTree {
        /*
        ADD_BRANCH_T(double, InitialBeamE)
        ADD_BRANCH_T(double, ProtonIm)
        ADD_BRANCH_T(double, PhotonIm)
        ADD_BRANCH_T(double, ProtonE)
        ADD_BRANCH_T(double, ProtonPolarAngles)
        ADD_BRANCH_T(double, missing_proton_Im)
        */
        ADD_BRANCH_T(unsigned, nClusters)
        ADD_BRANCH_T(double, TaggW)
        ADD_BRANCH_T(std::vector<double>, PhotonAzimuthAngles)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAngles)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAnglesCB)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAnglesTAPS)

    };

protected:

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    //Distinguish between uncertainty models between MC and exp. data:
    using model_t = std::shared_ptr<const utils::UncertaintyModels::Interpolated>;
    model_t fit_model_data;
    model_t fit_model_mc;

    //utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;

    //-- Histograms

    static const int nrCuts_total = 8;
    static const int nrCuts_beforeSel = 2;
    static const int nrCuts_beforePi0 = 4;
    static const int nrCuts_beforeKF = 6;

    static const int nrCuts_pi0 = nrCuts_total-nrCuts_beforePi0;
    static const int nrCuts_Sel = nrCuts_total-nrCuts_beforeSel;
    static const int nrCutsKF = nrCuts_total-nrCuts_beforeKF;

    static const int neu_nrSel = 3;
    static const int cha_nrSel = 1;

    static const int nrPartType = 2;
    static const int nrFitVars = 4;

    const int nr_pi0g = 2;

private:

    tree_t t;

    TH1D* h_beamE[nrCuts_total];
    TH1D* h_missingP_IM[nrCuts_Sel];
    TH1D* h_3g_IM[nrCuts_Sel];
    TH1D* h_2gComb_IM[nrCuts_Sel];
    TH2D* h_doublyDCScm_gp_wp[nrCuts_Sel];
    TH2D* h_doublyDCSlab_gp_wp[nrCuts_Sel];

    TH1D* h_wp_BackToBack[nrCuts_Sel];
    TH2D* h_Proton_PIDEvsTime[nrCuts_Sel];

    TH2D* h_Proton_CaloEvsVetoE_CB[nrCuts_Sel];
    TH2D* h_Proton_CaloEvsVetoE_TAPS[nrCuts_Sel];

    TH2D* h_AllCaloEvsVetoE_CB[nrCuts_total];
    TH2D* h_AllCaloEvsVetoE_TAPS[nrCuts_total];

    TH1D* h_neuCaloE[nrCuts_total];
    TH1D* h_chaCaloE[nrCuts_total];
    TH1D* h_AllVetoE_CB[nrCuts_total];
    TH1D* h_AllVetoE_TAPS[nrCuts_total];
    TH2D* h_neuEkinVSThetaCB[nrCuts_total];
    TH2D* h_neuEkinVSThetaTAPS[nrCuts_total];
    TH2D* h_chaEkinVSThetaCB[nrCuts_total];
    TH2D* h_chaEkinVSThetaTAPS[nrCuts_total];
    TH2D* h_neuEkinVSPhi[nrCuts_total];
    TH2D* h_chaEkinVSPhi[nrCuts_total];

    TH1D* h_neuTimeDiffCorTaggCB[nrCuts_total];
    TH1D* h_neuTimeDiffCorTaggTAPS[nrCuts_total];
    TH1D* h_chaTimeDiffCorTaggCB[nrCuts_total];
    TH1D* h_chaTimeDiffCorTaggTAPS[nrCuts_total];

    TH1D* h_2gPi0_IM[nrCuts_pi0];
    TH1D* h_wpi0g_BackToBack[nrCuts_pi0];
    TH1D* h_pi0gg_BackToBack[nrCuts_pi0];
    TH2D* h_EpivsIM3g[nrCuts_pi0];

    //TH1D* h_pi0g_BackToBack[nrCuts_BackToBack];

    TH1D* h_TaggerTime;
    TH1D* h_InitialBeamE;
    TH1D* h_VetoEnergies;
    TH1D* h_nClusters;
    TH1D* h_nCandidates;
    TH1D* h_nClusters_pr;

    TH1D* h_NeuAzimuthAngles;
    TH1D* h_NeuPolarAngles;
    TH1D* h_ChaPolarAngles;
    TH1D* h_ChaAzimuthAngles;
    TH1D* h_3gPolarAngles;
    TH1D* h_3gPolarAnglesCB;
    TH1D* h_3gPolarAnglesTAPS;
    TH1D* h_3gAzimuthAngles;

    TH1D* h_RecData_relStat;
    TH1D* h_RecData_Stat;

    TH2D* h_3g_EvTheta_CB;
    TH2D* h_3g_EvTheta_TAPS;
    TH2D* h_p_EvTheta;
    TH2D* h_w_EvTheta;
    TH2D* h_wg_EvTheta;
    TH2D* h_wpi0_EvTheta;
    TH2D* h_wpi02g_EvTheta;
    //TH1D* h_Reconstructed_Data_Statistics;

    //KinFit-overview:
    //TH1D* h_Steps;
    TH1D* h_Probability[nrCutsKF];
    TH1D* h_Fit_zvert[nrCutsKF];
    TH1D* h_fitEbeam[nrCutsKF];
    TH1D* h_IM3g_Fit[nrCutsKF];
    TH1D* h_IM2gPi0_Fit[nrCutsKF];

    TH2D* h_Ek_dev_CB[nrPartType];
    TH2D* h_Theta_dev_CB[nrPartType];
    TH2D* h_Phi_dev_CB[nrPartType];

    TH2D* h_Ek_dev_TAPS[nrPartType];
    TH2D* h_Theta_dev_TAPS[nrPartType];
    TH2D* h_Phi_dev_TAPS[nrPartType];

    TH1D *h_PartPulls_CB[nrPartType][nrFitVars];
    TH1D *h_PartPulls_TAPS[nrPartType][nrFitVars];

    TH1D* h_CBEsum[nrCuts_total];

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    std::string cuts[nrCuts_total] = {"CUT#0_NoCuts", "CUT#1_CBEsum", "CUT#2_Sel3Neu1Cha", "CUT#3_ImMissingParticle_+-2sigma_mp", "CUT#4_OmegaEthreshold", "CUT#5_SelMinM(2neu-mpi0)_+-2sigma_mpi0", "CUT#6_neuCaloEthreshold", "CUT#7_kinFit_prob_CL1%"};
    std::string cutsKF[nrCutsKF] = {"after_kinFit","after_1%CL_CUT"};
    std::string fitPartName[nrPartType] ={"protons" , "photons"};
    std::string fitvarnameCB[nrFitVars] = {"invEk","theta","phi","R"};
    std::string fitvarnameTA[nrFitVars] = {"invEk","Rxy","phi","L"};

    double max_particles = 1000000;
    double vetoEthreshold = 0.2;
    double neuCaloEthreshhold = 60;
    long double mpi0 = 134.9766;
    long double mp = 938.2720813;
    long double mw = 782.65;

    int lower_edge = 0;
    int upper_edge = nrCuts_total;
    int number_of_bins = nrCuts_total*15;
    int steps = (int)(number_of_bins/(upper_edge-lower_edge));

    double_t energyGamma_min = 1105;
    double_t energyGamma_max = 1500;

    double_t s_square_min = sqrt(2*mp*energyGamma_min+mp*mp);
    double_t s_square_max = sqrt(2*mp*energyGamma_max+mp*mp);

    long double Omega_Ethreshold = (mw*mw+2*mw*mp)/(2*mp);

    double PhotonPolarAngles_TAPS_frac;
    double ProtonPolarAngles_TAPS_frac;
    double PhotonPolarAngles_CB_frac;
    double ProtonPolarAngles_CB_frac;

    double PhotonPolarAngles_TAPS_err;
    double PhotonPolarAngles_CB_err;
    double ProtonPolarAngles_TAPS_err;
    double ProtonPolarAngles_CB_err;

    long double stat[nrCuts_total] = {0};

    double weight_res = 0;

    double mpFit = 926.144;
    double mpi0Fit = 136.905;

    double sigmaMissingP = 61.1838;
    double sigmaPi0IM = 16.4186;

    //double sigmaMissingP = 50.138;
    //double sigmaPi0IM = 12.494;

    std::shared_ptr<expconfig::detector::TAPS> taps;

};

}}}
