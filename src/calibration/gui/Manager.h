#pragma once

#include "Manager_traits.h"
#include "AvgBuffer.h"

#include "base/interval.h"
#include "tree/TDataRecord.h"

#include <memory>
#include <list>
#include <vector>
#include <string>

class TH1D;
class TH2D;
class TFile;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvasMode;

class Manager {

protected:

    using myBuffer_t = AvgBuffer<TH2D, interval<TID>>;

    std::shared_ptr<Manager_traits> module;
    myBuffer_t buffer;

    struct input_file_t {

        input_file_t(const std::string& FileName, const interval<TID>& R):
            filename(FileName),
            range(R) {}

        bool operator< (const input_file_t& other) const;

        const std::string filename;
        const ant::interval<ant::TID> range;
    };

    std::list<input_file_t> input_files;

    struct state_t {
        state_t() :
              is_init(false),
              breakpoint_fit(false),
              breakpoint_finish(false)
        {}

        std::list<input_file_t>::iterator it_file;
        myBuffer_t::const_iterator it_buffer;
        int channel;

        bool is_init;
        bool breakpoint_fit;
        bool breakpoint_finish;
    };
    state_t state;


    std::unique_ptr<CalCanvasMode> mode;

    void BuildInputFiles(const std::vector<std::string>& filenames);

    void FillWorklistFromFiles();

    struct SignalConnection_t {

        std::string receiver_class;
        void* receiver = nullptr;
        std::string slot;
    };
    SignalConnection_t signalConnection;

    int maxChannels;
    bool DoInit();

public:
    std::string SetupName;

    Manager(const std::vector<std::string>& inputfiles, unsigned avglength);

    virtual void SetModule(const std::shared_ptr<Manager_traits>& module_) {
        module = move(module_);
    }

    virtual void ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot);

    virtual bool Run();

    virtual ~Manager();

};

}
}
}
