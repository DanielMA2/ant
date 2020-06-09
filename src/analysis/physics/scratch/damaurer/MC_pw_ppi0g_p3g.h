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

// handle ROOT's TTree
#include "base/WrapTTree.h"

// the physics classes reside in this nested namespace
namespace ant {

namespace expconfig {
namespace detector {
    struct TAPS;
}}
namespace analysis {
namespace physics {

class scratch_damaurer_MC_pw_ppi0g_p3g : public Physics {
public:
    scratch_damaurer_MC_pw_ppi0g_p3g(const std::string& name, OptionsPtr opts);
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
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

//-- Histograms
static const int nrCuts = 2;

private:

    tree_t t;

    TH1D *h_missing_proton_Im[nrCuts];
    TH1D *h_wOnly3g_Im[nrCuts];
    TH1D *h_Pi0Only2g_Im[nrCuts];
    TH1D *h_pi0g_BackToBack[nrCuts];

    TH1D* h_TaggerTime;
    TH1D* h_InitialBeamE;
    TH1D* h_VetoEnergies;
    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;

    TH1D* h_ALL_AzimuthAngles;
    TH1D* h_ALL_PolarAngles;
    TH1D* h_missingChaPolarAngles;
    TH1D* h_missingChaAzimuthAngles;
    TH1D* h_3gPolarAngles;
    TH1D* h_3gPolarAnglesCB;
    TH1D* h_3gPolarAnglesTAPS;
    TH1D* h_3gAzimuthAngles;

    TH2D* h_3g_EvTheta_CB;
    TH2D* h_3g_EvTheta_TAPS;
    TH2D* h_missingP_EvTheta;
    TH2D* h_w_EvTheta;
    TH2D* h_wg_EvTheta;
    TH2D* h_wpi0_EvTheta;

    TH2D* h_doubly_wp_DCS_reconstructed_lab;
    TH2D* h_doubly_wp_DCS_reconstructed_cmFrame;
    //TH1D* h_Reconstructed_Data_Statistics;

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    double max_particles = 1000000;

    long double m = 134.9766;
    long double mp = 938.2720813;
    long double mw = 762.65;
    long double energyGamma_min = 1105;
    long double energyGamma_max = 1500;

    double PhotonPolarAngles_TAPS_frac;
    double ProtonPolarAngles_TAPS_frac;
    double PhotonPolarAngles_CB_frac;
    double ProtonPolarAngles_CB_frac;

    double PhotonPolarAngles_TAPS_err;
    double PhotonPolarAngles_CB_err;
    double ProtonPolarAngles_TAPS_err;
    double ProtonPolarAngles_CB_err;

    double_t s_square_min = sqrt(2*mp*energyGamma_min+mp*mp);
    double_t s_square_max = sqrt(2*mp*energyGamma_max+mp*mp);

    std::shared_ptr<expconfig::detector::TAPS> taps;

};

}}}
