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
namespace analysis {
namespace physics {

class scratch_damaurer_MC_ppi0_pgg : public Physics {
public:
    scratch_damaurer_MC_ppi0_pgg(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

    struct tree_t : WrapTTree {
        ADD_BRANCH_T(unsigned, nClusters)
        ADD_BRANCH_T(double, TaggW)
        ADD_BRANCH_T(double, InitialBeamE)
        ADD_BRANCH_T(double, ProtonIm)
        ADD_BRANCH_T(double, PhotonIm)
        ADD_BRANCH_T(double, ProtonE)
        ADD_BRANCH_T(double, ProtonPolarAngles)
        ADD_BRANCH_T(double, ProtonAzimuthAngles)
        ADD_BRANCH_T(std::vector<double>, PolarAngles)
        ADD_BRANCH_T(std::vector<double>, AzimuthAngles)
        ADD_BRANCH_T(std::vector<double>, PhotonE)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAngles)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAnglesCB)
        ADD_BRANCH_T(std::vector<double>, PhotonPolarAnglesTAPS)
        ADD_BRANCH_T(std::vector<double>, PhotonAzimuthAngles)
    };

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

private:
    tree_t t;
    TH1D* h_TaggerTime;
    TH1D* h_InitialBeamE;
    TH1D* h_VetoEnergies;
    TH1D* h_PolarAngles;
    TH1D* h_ProtonPolarAngles;
    TH1D* h_PhotonPolarAngles;
    TH1D* h_PhotonPolarAngles_CB;
    TH1D* h_PhotonPolarAngles_TAPS;
    TH1D* h_AzimuthAngles;
    TH1D* h_ProtonAzimuthAngles;
    TH1D* h_PhotonAzimuthAngles;
    TH1D* h_nClusters;
    TH1D* h_nClusters_pr;
    TH1D* h_PhotonIm;
    TH1D* h_ProtonIm;
    TH1D* h_ProtonImAbove80;
    TH2D* h_ProtonETheta;
    TH2D* h_PhotonETheta;
    TH2D* h_PhotonETheta_CB;
    TH2D* h_PhotonETheta_TAPS;
    TH2D* h_ProtonEPhi;
    TH2D* h_PhotonEPhi;
    TH1D* h_alp_Photons_CB_TAPS;
    TH1D* h_Reconstructed_Data_Statistics;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    double PhotonPolarAngles_TAPS = 0.0;
    double ProtonPolarAngles_TAPS = 0.0;
    double PhotonPolarAngles_CB = 0.0;
    double ProtonPolarAngles_CB = 0.0;

    double count_first = 0;
    double count_second = 0;
    double count_third = 0;
    double count_fourth = 0;
    double count_fifth = 0;
    double count_sixth = 0;
    double count_seventh = 0;
    double count_eigth = 0;
    double count_nineth = 0;
    double count_tenth = 0;

    double max_particles = 300000.0;
    double Photons_max = 200000.0;
    double Protons_max = 100000.0;

    double PhotonPolarAngles_TAPS_frac;
    double ProtonPolarAngles_TAPS_frac;
    double PhotonPolarAngles_CB_frac;
    double ProtonPolarAngles_CB_frac;

    double PhotonPolarAngles_TAPS_err;
    double PhotonPolarAngles_CB_err;
    double ProtonPolarAngles_TAPS_err;
    double ProtonPolarAngles_CB_err;

    Int_t number_of_bins = 110;
    Int_t lower_edge = 0;
    Int_t upper_edge = 11;
    int steps = (int)(number_of_bins/(upper_edge-lower_edge));
};

}}}
