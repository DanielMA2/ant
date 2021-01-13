// better than using the standard ifndef include guards...
#pragma once

//test

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

//Including detecor properties
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "expconfig/detectors/CB.h"

//Include Kinematic fitter
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "analysis/utils/ProtonPhotonCombs.h"

// the physics classes reside in this nested namespace
namespace ant {
namespace analysis {
namespace physics {

// choose a nice name for your class
class scratch_damaurer_omega_Dalitz : public Physics {
public:

    // use "struct" to avoid all those public keywords...
    struct tree_t : WrapTTree {
        // define two branches with ADD_BRANCH_T (again, this is black macro magic)
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(unsigned, nClusters)
        // you may even use more complex types with non-default ctor, for example
//        ADD_BRANCH_T(std::vector<double>, Numbers, 3) // vector Numbers is three items large now
    };

protected:

    // use this instance to handle prompt-random weighting
    PromptRandom::Switch promptrandom;
    // this encapsulates trigger timings and other related things (CBESum)
    utils::TriggerSimulation triggersimu;

    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<expconfig::detector::TAPSVeto> veto_detector;

    static const int nrCuts_total = 6;
    std::string cuts[nrCuts_total] = {"CUT#0_NoCuts","CUT#1_Sel2Neu3Cha","CUT#2_OmegaEthreshold","CUT#3_IM(2g)_mpi0+-0.4mpi0","CUT#4_mm(p)ANDkf","CUT#5_kinFit_prob_CL1%"};

    double max_particles = 1000000;
    double vetoEthreshold = 0.2;
    long double mpi0 = 134.9766;
    long double mp = 938.2720813;
    long double mw = 782.65;

    static const int neu_nrSel = 2;
    static const int cha_nrSel = 3;

    long double stat[nrCuts_total] = {0};

    long double Omega_Ethreshold = (mw*mw+2*mw*mp)/(2*mp);
    static const int nrCombs = 3;

private:

    TH1D* h_RecData_Stat;

    TH1D* h_TaggerTime;
    TH1D* h_nClusters;   
    TH1D* h_nCandidates;
    TH1D* h_Probability;
    TH2D* h_nNeuChaCandidates;
    TH1D* h_mmpFails;
    TH1D* h_kfFails;
    TH1D* h_totalFails;

    TH1D* h_AllVetoE_CB[nrCuts_total];
    TH1D* h_AllVetoE_TAPS[nrCuts_total];
    TH1D* h_beamE[nrCuts_total];
    TH2D* h_AllCaloEvsVetoE_CB[nrCuts_total];
    TH2D* h_AllCaloEvsVetoE_TAPS[nrCuts_total];
    //TH1D* hist;

    TH1D* h_2g_IM[nrCuts_total-1];
    TH1D* h_missingP_IM[nrCuts_total-1];
    TH1D* h_2gee_IM[nrCuts_total-1];

    tree_t t;

public:
    // physics need to implement this public constructor...
    scratch_damaurer_omega_Dalitz(const std::string& name, OptionsPtr opts);

    // ...and the following method (use override keyword)
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    // optional method (called in interactive mode to show histograms)
    virtual void ShowResult() override;
    virtual void Finish() override;

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);
};


}}} // end of ant::analysis::physics
