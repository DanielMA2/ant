#pragma once

#include "calibration/Calibration.h"
#include "base/Detector_t.h"
#include "base/OptionsList.h"
#include "base/Detector_t.h"

#include "tree/TID.h" // for TKeyValue, TID

#include <memory>


namespace ant {

namespace calibration {

class DataManager;

class ClusterSmearing :
        public Calibration::Module, // this makes this module abstract
        public ReconstructHook::Clusters
{

public:
    // ReconstructHook
    virtual void ApplyTo(clusters_t& clusters) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

    ClusterSmearing(std::shared_ptr<ClusterDetector_t> det,
           std::shared_ptr<DataManager> calmgr);
    virtual ~ClusterSmearing();

protected:


    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >&, const ant::OptionsPtr) override;

    const Detector_t::Type_t DetectorType;
    const unsigned nelements;

    std::shared_ptr<DataManager> calibrationManager;

    struct SigmaInterpolator;
    std::unique_ptr<SigmaInterpolator> interpolator;

};

}}  // namespace ant::calibration