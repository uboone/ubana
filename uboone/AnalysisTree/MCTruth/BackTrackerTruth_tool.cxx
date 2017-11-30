
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <cmath>
#include <algorithm>

namespace truth
{
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BackTrackerTruth::BackTrackerTruth(fhicl::ParameterSet const & pset)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("BackTrackerTruth") << "BackTrackerTruth configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BackTrackerTruth::~BackTrackerTruth()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BackTrackerTruth::reconfigure(fhicl::ParameterSet const & pset)
{
    fTrackProducerLabel = pset.get<art::InputTag>("TrackProducerLabel", "");
}
    
//----------------------------------------------------------------------
const sim::ParticleList& BackTrackerTruth::ParticleList() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->ParticleList();
}

    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToParticle(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToParticle(id);
}
    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToMotherParticle(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToMotherParticle(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::TrackIDToMCTruth(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToMCTruth(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::ParticleToMCTruth(const simb::MCParticle* p) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->ParticleToMCTruth(p);
}
    
//----------------------------------------------------------------------
std::vector<const simb::MCParticle*> BackTrackerTruth::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->MCTruthToParticles(mct);
}
    
//----------------------------------------------------------------------
const std::vector< art::Ptr<simb::MCTruth> >&  BackTrackerTruth::MCTruthVector() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->MCTruthVector();
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(recob::Hit const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToTrackID(hit);
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToTrackID(*hit);
}

//----------------------------------------------------------------------
const std::vector<std::vector<art::Ptr<recob::Hit>>> BackTrackerTruth::TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                 std::vector<int> const& tkIDs) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDsToHits(allhits,tkIDs);
}
    
//----------------------------------------------------------------------
// plist is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> BackTrackerTruth::HitToEveID(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToEveID(hit);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfEveIDs();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfTrackIDs();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfEveIDs(hits);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfTrackIDs(hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitCollectionPurity(trackIDs, hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitChargeCollectionPurity(trackIDs, hits);
}
    
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionEfficiency(std::set<int> trackIDs,
                                            std::vector< art::Ptr<recob::Hit> > const& hits,
                                            std::vector< art::Ptr<recob::Hit> > const& allhits,
                                            geo::View_t const& view) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                  std::vector< art::Ptr<recob::Hit> > const& hits,
                                                  std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                  geo::View_t                         const& view) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitChargeCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToXYZ(hit);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointToXYZ(art::Ptr<recob::SpacePoint>         const& spt,
                                                 art::Event                          const& evt,
                                                 std::string                         const& label) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->SpacePointToXYZ(spt, evt, label);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->SpacePointHitsToXYZ(hits);
}

//----------------------------------------------------------------------------
    
//DEFINE_ART_CLASS_TOOL(BackTrackerTruth)
}
