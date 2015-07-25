#pragma once

#include "Detector_t.h"
#include "unpacker/UnpackerAcqu.h"
#include <stdexcept>

namespace ant {
namespace expconfig {
namespace detector {
struct Trigger :
        Detector_t,
        UnpackerAcquConfig
{

    Trigger() : Detector_t(Detector_t::Type_t::Trigger) {}

    virtual TVector3 GetPosition(unsigned) const override {
        // when you ask the trigger detector for positions,
        // this is certainly a bug :)
        throw std::runtime_error("The trigger detector knows nothing about positions.");
    }

    const LogicalChannel_t Reference_CATCH_TaggerCrate = {Type, Channel_t::Type_t::Timing, 1000};
    const LogicalChannel_t Reference_CATCH_CBCrate = {Type, Channel_t::Type_t::Timing, 1001};

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

};


struct Trigger_2014 : Trigger {
    const unsigned Scaler_Exptrigger_1MHz = 10;
    const unsigned Scaler_Beampolmon_1MHz = 20;

    virtual bool Matches(const THeaderInfo& headerInfo) const override;
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

};

}}} // namespace ant::expconfig::detector
