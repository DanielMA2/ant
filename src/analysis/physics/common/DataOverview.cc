#include "physics/common/DataOverview.h"
#include "plot/root_draw.h"
#include "expconfig/ExpConfig.h"
#include <cassert>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

TaggerOverview::TaggerOverview(const string &name, PhysOptPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_hits(50);
    const BinSettings bins_energy(200, 1400, 1600);
    const BinSettings bins_channels(47);
    const BinSettings bins_times(1000,-50,50);

//    const auto setup = ExpConfig::Setup::GetLastFound();
//    const auto tagger = setup ? setup->GetDetector<TaggerDetector_t>() : nullptr;

//    if(!tagger) {
//        throw std::runtime_error("No Tagger in Setup!");
//    }

//    const BinSettings bins_channels(tagger->GetNChannels());

//    const auto Emin   = tagger->GetPhotonEnergy(0);
//    const auto Emax   = tagger->GetPhotonEnergy(tagger->GetNChannels()-1);
//    const auto width  = (Emax - Emin)/tagger->GetNChannels();

//    const BinSettings bins_energy(tagger->GetNChannels(), Emin-width/2, Emax+width/2);

    nHitsEvent = HistFac.makeTH1D("Tagger Hits / Enent ", "# Hits/Event",   "",     bins_hits,    "TaggerHitsPerEvent");
    nHitsEvent->SetFillColor(kRed);

    Channels   = HistFac.makeTH1D("Tagger Channels hit ", "Channel Number", "Hits", bins_channels," TaggerChannels");
    Channels->SetFillColor(kYellow);

    Energies   = HistFac.makeTH1D("Tagged Photon Energies ", "E_{#gamma,tag} [MeV]", "", bins_energy, "TaggedEnergies");
    Energies->SetFillColor(kBlue);

    Times      = HistFac.makeTH1D("Tagger Hit Times ", "Time [ns]", "", bins_times, "TaggedTimes");
    Times->SetFillColor(kGreen);

    channel_correlation = HistFac.makeTH2D("Tagger Channel Correlation","Channel", "Channel", bins_channels, bins_channels, "ChannelCorreleation");
}

TaggerOverview::~TaggerOverview()
{
}

void TaggerOverview::ProcessEvent(const Event &event)
{
    const auto taggerhits = (mode == Mode::Reconstructed) ? event.Reconstructed().TaggerHits() : event.MCTrue().TaggerHits();

    for(auto hit=taggerhits.cbegin(); hit!=taggerhits.cend(); ++hit) {
        Channels->Fill((*hit)->Channel());
        Energies->Fill((*hit)->PhotonEnergy());
        Times->Fill((*hit)->Time());

        for(auto hit2 = next(hit); hit2!=taggerhits.cend(); ++hit2) {
            channel_correlation->Fill((*hit)->Channel(), (*hit2)->Channel());
            channel_correlation->Fill((*hit2)->Channel(), (*hit)->Channel());
        }
    }
    nHitsEvent->Fill(taggerhits.size());
}

void TaggerOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << nHitsEvent
            << Channels
            << Energies
            << Times
            << drawoption("colz") << channel_correlation
            << endc;

}


DataOverviewBase::DataOverviewBase(const string &name, PhysOptPtr opts):
    Physics(name, opts)
{
    if(opts->Get<string>("Mode") == "Reconstructed")
        mode = Mode::Reconstructed;
    else if(opts->Get<string>("Mode") == "MCTrue")
        mode = Mode::MCTrue;

    HistFac.SetTitlePrefix(GetMode());
}

DataOverviewBase::~DataOverviewBase()
{}

string DataOverviewBase::GetMode() const
{
    if(mode == Mode::Reconstructed) {
        return "Reconstructed";
    } else
        return "MCTrue";
}

const Event::Data &DataOverviewBase::GetBranch(const Event &event) const
{
   return (mode == Mode::Reconstructed) ? event.Reconstructed() : event.MCTrue();
}


TriggerOverview::TriggerOverview(const string &name, PhysOptPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_errors(100);
    const BinSettings bins_multiplicity(10);
    const BinSettings bins_energy(1600);

    CBESum       = HistFac.makeTH1D("CB Energy Sum",  "CB Energy Sum [MeV]", "", bins_energy,      "CBESum");
    Multiplicity = HistFac.makeTH1D("Multiplicity",   "# Hits",              "", bins_multiplicity,"Multiplicity");
    nErrorsEvent = HistFac.makeTH1D("Errors / Event", "# errors",            "", bins_errors,      "nErrrorsEvent");
}

TriggerOverview::~TriggerOverview()
{}

void TriggerOverview::ProcessEvent(const Event &event)
{
    const auto TriggerInfo = GetBranch(event).TriggerInfos();

    CBESum->Fill(TriggerInfo.CBEenergySum());
    Multiplicity->Fill(TriggerInfo.Multiplicity());
    nErrorsEvent->Fill(TriggerInfo.Errors().size());
}

void TriggerOverview::ShowResult()
{
    canvas(this->GetName()+" "+GetMode())
            << CBESum
            << Multiplicity
            << nErrorsEvent
            << endc;
}



void ParticleOverview::SetBinLabels(TH1D *hist, const ParticleTypeDatabase::TypeList_t &types)
{
    assert(int(types.size()) <= hist->GetNbinsX());
    int i=1;
    for(const auto& type : types) {
        hist->GetXaxis()->SetBinLabel(i++, type->PrintName().c_str());
    }
}

ParticleOverview::ParticleOverview(const string &name, PhysOptPtr opts):
    DataOverviewBase(name, opts)
{
    const BinSettings bins_particles(15);

    nParticles    = HistFac.makeTH1D("Particles / Event", "# Particles", "", bins_particles, "nParticles");
    particleTypes = HistFac.makeTH1D("Particle Types",    "Type",        "", BinSettings(unsigned(ParticleTypeDatabase::DetectableTypes().size())), "ParticleTypes");
    SetBinLabels(particleTypes, ParticleTypeDatabase::DetectableTypes());

    for(const ParticleTypeDatabase::Type* type : ParticleTypeDatabase::DetectableTypes()) {
        nType[type] = HistFac.makeTH1D(type->PrintName() + " / Event", "# Particles", "", bins_particles, type->PrintName());
    }
}

ParticleOverview::~ParticleOverview()
{

}

void ParticleOverview::ProcessEvent(const Event &event)
{
    const auto& particles = GetBranch(event).Particles().GetAll();
    nParticles->Fill(particles.size());

    for(const auto& p : particles) {
        particleTypes->Fill(p->Type().PrintName().c_str(), 1.0);
    }

    for(const ParticleTypeDatabase::Type* type : ParticleTypeDatabase::DetectableTypes()) {
        nType[type]->Fill(GetBranch(event).Particles().Get(*type).size());
    }
}

void ParticleOverview::ShowResult()
{
    canvas c(GetName());

    c << nParticles
      << particleTypes;

    for(const auto& entry : nType) {
        c<< entry.second;
    }

    c  << endc;
}

AUTO_REGISTER_PHYSICS(ParticleOverview)
AUTO_REGISTER_PHYSICS(TaggerOverview)
AUTO_REGISTER_PHYSICS(TriggerOverview)

