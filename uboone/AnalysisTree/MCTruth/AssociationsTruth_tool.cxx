
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"

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
AssociationsTruth::AssociationsTruth(fhicl::ParameterSet const & pset) :
    fMCTruthAssociations(pset.get<fhicl::ParameterSet>("MCTruthAssociations"))
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("AssociationsTruth") << "AssociationsTruth configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
AssociationsTruth::~AssociationsTruth()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void AssociationsTruth::reconfigure(fhicl::ParameterSet const & pset)
{
    fAssnsProducerLabels = pset.get<std::vector<art::InputTag>>("AssnsProducerLabels");
    fG4ProducerLabel     = pset.get<art::InputTag>             ("G4ProducerLabel");
}

//----------------------------------------------------------------------------
/// Rebuild method -> rebuild the basic maps to get truth information
///
/// Arguments:
///
/// event - the art event used to extract all information
///
void AssociationsTruth::Rebuild(const art::Event& evt)
{
    // Create a container for testing
    HitParticleAssociationsVec partHitAssnsVec;

    // Get a handle for the associations...
    art::Handle<art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>> partHitAssnsHandle;
    
    for(const auto& assnsProducerLabel : fAssnsProducerLabels)
    {
        evt.getByLabel(assnsProducerLabel, partHitAssnsHandle);
    
        if (!partHitAssnsHandle.isValid())
        {
            throw cet::exception("AssociationsTruth") << "===>> NO MCParticle <--> Hit associations found for run/subrun/event: " << evt.run() << "/" << evt.subRun() << "/" << evt.id().event() << std::endl;
        }
    
        partHitAssnsVec.emplace_back(&*partHitAssnsHandle);
    }
    
    // Recover the associations between MCTruth and MCParticles
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    evt.getByLabel(fG4ProducerLabel, mcParticleHandle);
    
    std::vector<art::Ptr<simb::MCParticle>> mcParticlePtrVec;
    art::fill_ptr_vector(mcParticlePtrVec, mcParticleHandle);
    
    MCTruthAssns mcTruthAssns(mcParticleHandle, evt, fG4ProducerLabel);
    
    // Pass this to the truth associations code
    fMCTruthAssociations.setup(partHitAssnsVec, mcParticlePtrVec, mcTruthAssns, *fGeometry, *fDetectorProperties);
    
    // Ugliness to follow! Basically, we need to build the "particle list" and the current implementation of
    // that code requires a copy...
    const MCTruthParticleList& locParticleList = fMCTruthAssociations.getParticleList();
    
    // Clear the current container
    fParticleList.clear();
    
    // Now we add particles back in one at a time...
    for(const auto& element : locParticleList)
    {
        fParticleList.Add(new simb::MCParticle(*(element.second)));
    }
    
    // Just to be consistent with the backtracker...
    fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
}

//----------------------------------------------------------------------
const sim::ParticleList& AssociationsTruth::ParticleList() const
{
    // Unfortunately, this requires special handling at the moment...
    return fParticleList;
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& AssociationsTruth::TrackIDToMCTruth(int const& id) const
{
    return fMCTruthAssociations.TrackIDToMCTruth(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& AssociationsTruth::ParticleToMCTruth(const simb::MCParticle* p) const
{
    return fMCTruthAssociations.ParticleToMCTruth(p);
}
    
//----------------------------------------------------------------------
std::vector<const simb::MCParticle*> AssociationsTruth::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    return fMCTruthAssociations.MCTruthToParticles(mct);
}
    
//----------------------------------------------------------------------
const std::vector< art::Ptr<simb::MCTruth> >&  AssociationsTruth::MCTruthVector() const
{
    return fMCTruthAssociations.MCTruthVector();
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> AssociationsTruth::HitToTrackID(recob::Hit const& hit) const
{
    std::vector<truth::TrackIDE> locTrackIDEVec = fMCTruthAssociations.HitToTrackID(&hit);
    std::vector<sim::TrackIDE>   outputVec;
    
    outputVec.reserve(locTrackIDEVec.size());

    for(const auto& trackIDE : locTrackIDEVec) outputVec.emplace_back(trackIDE.trackID,trackIDE.energyFrac,trackIDE.energy,trackIDE.numElectrons);
    
    return outputVec;
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> AssociationsTruth::HitToTrackID(art::Ptr<recob::Hit> const& hit) const
{
    return HitToTrackID(*hit.get());
}

//----------------------------------------------------------------------
const std::vector<std::vector<art::Ptr<recob::Hit>>> AssociationsTruth::TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                       std::vector<int> const& tkIDs) const
{
    return fMCTruthAssociations.TrackIDsToHits(allhits,tkIDs);
}
    
//----------------------------------------------------------------------
// plist is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> AssociationsTruth::HitToEveID(art::Ptr<recob::Hit> const& hit) const
{
    std::vector<truth::TrackIDE> locTrackIDEVec = fMCTruthAssociations.HitToEveID(hit);
    std::vector<sim::TrackIDE>   outputVec;
    
    outputVec.reserve(locTrackIDEVec.size());
    
    for(const auto& trackIDE : locTrackIDEVec) outputVec.emplace_back(trackIDE.trackID,trackIDE.energyFrac,trackIDE.energy,trackIDE.numElectrons);
    
    return outputVec;
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfEveIDs() const
{
    return fMCTruthAssociations.GetSetOfEveIDs();
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfTrackIDs() const
{
    return fMCTruthAssociations.GetSetOfTrackIDs();
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.GetSetOfEveIDs(hits);
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.GetSetOfTrackIDs(hits);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.HitCollectionPurity(trackIDs, hits);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitChargeCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.HitChargeCollectionPurity(trackIDs, hits);
}
    
    
//----------------------------------------------------------------------
double AssociationsTruth::HitCollectionEfficiency(std::set<int>                              trackIDs,
                                                  std::vector< art::Ptr<recob::Hit> > const& hits,
                                                  std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                  geo::View_t const&                         view) const
{
    return fMCTruthAssociations.HitCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                        std::vector< art::Ptr<recob::Hit> > const& hits,
                                                        std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                        geo::View_t                         const& view) const
{
    return fMCTruthAssociations.HitChargeCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
{
    return fMCTruthAssociations.HitToXYZ(hit);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::SpacePointToXYZ(art::Ptr<recob::SpacePoint> const& spt,
                                                       art::Event                  const& evt,
                                                       std::string                 const& label) const
{
    // Get hits that make up this space point.
    art::PtrVector<recob::SpacePoint> spv;
    spv.push_back(spt);
    art::FindManyP<recob::Hit> fmh(spv, evt, label);
    std::vector< art::Ptr<recob::Hit> > hitv = fmh.at(0);
    
    // make a PtrVector
    art::PtrVector<recob::Hit> hits;
    for(size_t h = 0; h < hitv.size(); ++h) hits.push_back(hitv[h]);
    
    return this->SpacePointHitsToXYZ(hits);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const
{
    return fMCTruthAssociations.SpacePointHitsToXYZ(hits);
}

//----------------------------------------------------------------------------
    
//DEFINE_ART_CLASS_TOOL(AssociationsTruth)
}
