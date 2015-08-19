#pragma once

#include "Calibration.h"

class TGraph;

namespace ant{
namespace calibration{

}
}

namespace ant {


namespace expconfig {
namespace detector {
class CB;
}}

namespace calibration {

class DataManager;

namespace gui {
class FitGaus;
}


class CB_TimeWalk :
        public Calibration::Module,
        public ReconstructHook::Clusters
{
public:
    CB_TimeWalk(
            const std::shared_ptr<expconfig::detector::CB>& cb,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~CB_TimeWalk();

    class ThePhysics : public analysis::Physics {
    protected:
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        TH3D* h_timewalk;

    public:
        ThePhysics(const std::string& name, const std::shared_ptr<expconfig::detector::CB>& cb);

        virtual void ProcessEvent(const analysis::data::Event& event) override;
        virtual void Finish() override ;
        virtual void ShowResult() override;
    }; // ThePhysics

    class TheGUI : public gui::Manager_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::CB> cb_detector;
        std::shared_ptr<gui::FitGaus> func;

//        gui::CalCanvas* c_singlechannel;
//        gui::CalCanvas* c_result;

//        TH1*  h_projection = nullptr;
//        TGraph* h_result;

        std::map< unsigned, std::vector<double> > fitParameters;


    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::CB>& cb
               );

        virtual std::string GetHistogramName() const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI();
        virtual std::list<gui::CalCanvas*> GetCanvases() const;

        virtual void StartRange(const interval<TID>& range);
        virtual DoFitReturn_t DoFit(TH1* hist, unsigned channel);
        virtual void DisplayFit();
        virtual void StoreFit(unsigned channel);
        virtual bool FinishRange();
        virtual void StoreFinishRange(const interval<TID>& range);
    }; // TheGUI

    virtual std::unique_ptr<analysis::Physics> GetPhysicsModule() override;
    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::Manager_traits> >& guis) override;

    virtual void ApplyTo(clusters_t& sorted_clusters) override;

    // Updateable_traits interface
    virtual std::vector<std::list<TID> > GetChangePoints() const override;
    virtual void Update(std::size_t index, const TID& id) override;

protected:
    struct timewalk_t {
        double offset;
        double scale;
        double shift;
        double exponent;
        double calc(const double E) const {
            return offset + scale/std::pow(E + shift, exponent);
        }
        TKeyValue<std::vector<double>> get() const {
            return {1,{offset,scale,shift,exponent}};
        }
        timewalk_t(const TKeyValue<std::vector<double>>& kv) :
            offset(kv.Value[0]),
            scale(kv.Value[1]),
            shift(kv.Value[2]),
            exponent(kv.Value[3])
        {}
        timewalk_t() :
            offset(0),
            scale(1),
            shift(0),
            exponent(-1)
        {}
    };

    std::vector<timewalk_t> timewalks;

    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<DataManager> calibrationManager;





};

}}
