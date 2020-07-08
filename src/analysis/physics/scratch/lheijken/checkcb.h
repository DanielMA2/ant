#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/CB.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A playground class for checking cb stuff
 *
 */
class scratch_lheijken_checkcb: public Physics {
protected:

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<expconfig::detector::PID> pid_detector;

    calibration::converter::GeSiCa_SADC adc_converter;
    calibration::converter::CATCH_TDC   tdc_converter;

    TH1D* hTrigRefTiming;
    TH1D *hDRHCBesum_1Cl, *hDRHCBesum_2Cl, *hDRHCBesum_3Cl, *hDRHCBesum_4Cl, *hDRHCBesum_5Cl, *hDRHCBesum;
    TH2D* hDRHUncalTimeAll;
    TH2D* hDRHCalTimeAll;
    TH2D* hDRHUncalTimeFirst;
    TH2D* hDRHCalTimeFirst;
    TH2D* hDRHUncalEnergy;
    TH2D* hDRHCalEnergy;
    TH2D* hDRHUncalTimeMult;
    TH2D* hDRHUncalEnergyMult;
    TH2D* hDRHUnCalEn_wTiming;
    TH2D* hDRHUnCalEn_woTiming;
    TH2D* hDRHCalEn_wTiming;
    TH2D* hDRHCalEn_woTiming;
    TH2D* hCHEnergy;
    TH2D* hCHTime;
//    TH3D* hCHTimeRawE;
    TH2D* hCBPhiPIDPhi;
    TH2D* hClTime;
    TH2D* hClEnergy;
    TH2D* hClTimeVsEnergy;
    TH2D* hCandClTime;
    TH2D* hCandClEnergy;
protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();

public:
    scratch_lheijken_checkcb(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_checkcb() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
