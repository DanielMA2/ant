#include "TSimpleParticle.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;


TSimpleParticle::TSimpleParticle(const TParticle& particle):
    TLorentzVector(particle),
    Time(std_ext::NaN),
    VetoE(std_ext::NaN),
    TrackerE(std_ext::NaN),
    ShortE(std_ext::NaN),
    ClusterSize(0)
{
    Mass = particle.Type().Mass();
    if (particle.Candidate)
    {
        Time        = particle.Candidate->Time;
        VetoE       = particle.Candidate->VetoEnergy;
        TrackerE    = particle.Candidate->TrackerEnergy;
        ClusterSize = particle.Candidate->ClusterSize;

        const auto caloCluster = particle.Candidate->FindCaloCluster();
        if (caloCluster) {
            ShortE  = caloCluster->ShortEnergy;
            TouchesHole = caloCluster->HasFlag(TCluster::Flags_t::TouchesHoleCentral);
        }
    }
}


