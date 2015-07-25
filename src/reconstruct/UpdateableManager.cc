#include "UpdateableManager.h"

#include "Reconstruct_traits.h"

#include "tree/TDataRecord.h" // for TID

#include "base/Logger.h"

#include <list>
#include <memory>
#include <map>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


UpdateableManager::UpdateableManager(
        const TID& headerInfo_ID,
        const std::list<std::shared_ptr<Updateable_traits> >& updateables
        )
{

    // ask each updateable for its update points and
    // build the list of changePoints

    map<TID, shared_ptr_list<Updateable_traits> > sorted_updateables;
    for(const shared_ptr<Updateable_traits>& updateable : updateables) {
        list<TID> all_changePoints = updateable->GetChangePoints();

        // no change points means the updateable is actually constant
        if(all_changePoints.empty())
            continue;

        // the following extraction relies on the changepoints being sorted in time
        // we don't require the updateables to provide a sorted list
        all_changePoints.sort();

        // scan the change points and build interesting_changePoints
        const TID* lastBeforeHeaderInfo = nullptr;
        list<TID> interesting_changePoints;
        for(const TID& changePoint : all_changePoints) {
            // ignore too early change points, but remember
            // the last one seen as a pointer
            if(changePoint < headerInfo_ID) {
                lastBeforeHeaderInfo = addressof(changePoint);
                continue;
            }
            interesting_changePoints.emplace_back(move(changePoint));
        }

        // check if we should add a timepoint before the headerInfo
        // in order to init the updateable correctly at the first detectorRead
        // in DoReconstruct()
        if(lastBeforeHeaderInfo != nullptr) {
            // prepend this last changepoint before headerInfo timepoint
            // only if the interesting changePoints are empty so far,
            // or the first interesting point is not equal the header timepoint
            if(interesting_changePoints.empty()
               || headerInfo_ID != interesting_changePoints.front())
                interesting_changePoints.emplace_front(move(*lastBeforeHeaderInfo));
        }

        // now the interesting points are built, add the updateable to the map
        for(const TID& changePoint : interesting_changePoints) {
            sorted_updateables[changePoint].emplace_back(move(updateable));
        }
    }

    // we rely on the fact that the map is sorted,
    // and simply convert it to a list for easier handling
    // in UpdateParameters()
    for(auto it_updateable : sorted_updateables) {
        // move the whole map pair into the list
        changePoints.emplace_back(move(it_updateable));
        }
}

void UpdateableManager::UpdateParameters(const TID& currentPoint)
{
    // it might be that the current point lies far in the future
    // so calling Update() more than once is necessary
    while(!changePoints.empty() && changePoints.front().first <= currentPoint) {

        unsigned nUpdateables = 0;
        for(const auto& updateable : changePoints.front().second) {
            updateable->Update(changePoints.front().first);
            nUpdateables++;
        }

        VLOG(7) << "Updated parameters for " << nUpdateables
                << " calibrations and detectors at ID=" << currentPoint;

        // go to next change point (if any)
        changePoints.pop_front();
    }
}
