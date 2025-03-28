/**
 * @file    MCTruthAssociations.cpp
 * @brief   Does something with the tracks (implementation file).
 * @author  Tracy Usher (usher@slac.stanford.edu
 * @date    October 24, 2017
 * @see     MCTruthAssociations.h
 * 
 */

#include "MCTruthAssociations.h"

// LArSoft libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools
#include "nug4/ParticleNavigation/EmEveIdCalculator.h"

// canvas libraries
#include "canvas/Persistency/Common/FindMany.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <algorithm> // std::count_if()

namespace truth
{

MCTruthAssociations::MCTruthAssociations(const fhicl::ParameterSet& config)
{
    fMinHitEnergyFraction = config.get<float>("MinHitEnergyFraction", 0.);
}

void MCTruthAssociations::setup(const HitParticleAssociationsVec&  partToHitAssnsVec,
                                const MCTruthParticleAssociations& truthToPartAssns,
                                const geo::GeometryCore&           geometry,
                                const detinfo::DetectorPropertiesData& detectorProperties)
{
    // Keep track of input services
    fGeometry           = &geometry;
    fDetectorProperties = &detectorProperties;
    
    fHitPartAssnsVec.clear();
    
    // Loop through the input vector of associations
    for(const auto* partToHitAssns : partToHitAssnsVec)
    {
        // Note that one element of our struct is non copiable which means
        // we need to construct in place before filling
        fHitPartAssnsVec.emplace_back();
    
        HitPartAssnsStruct& hitPartAssns = fHitPartAssnsVec.back();
    
        // Clear the maps in case they were previously filled
        hitPartAssns.fHitToPartVecMap.clear();
        hitPartAssns.fPartToHitVecMap.clear();
    
        // Build out the maps between hits/particles
        for(HitParticleAssociations::const_iterator partHitItr = partToHitAssns->begin(); partHitItr != partToHitAssns->end(); partHitItr++)
        {
            const art::Ptr<simb::MCParticle>&       mcParticle = partHitItr->first;
            const art::Ptr<recob::Hit>&             recoHit    = partHitItr->second;
            const anab::BackTrackerHitMatchingData* data       = partHitItr->data;
        
            hitPartAssns.fHitToPartVecMap[recoHit.get()].insert(PartMatchDataPair(mcParticle.get(),data));
            hitPartAssns.fPartToHitVecMap[mcParticle.get()].insert(HitMatchDataPair(recoHit.get(),data));
        }
        
        mf::LogDebug("MCTruthAssociations") << "Built maps with " << hitPartAssns.fHitToPartVecMap.size() << " hits, " << hitPartAssns.fPartToHitVecMap.size() << "\n";
    }
    
    // Note that there is only one instance of MCTruth <--> MCParticle associations so we do this external to the above loop
    fParticleList.clear();
    fMCTruthVec.clear();
    fTrackIDToMCTruthIndex.clear();
    
    for(const auto& truthPartAssn : truthToPartAssns)
    {
        const art::Ptr<simb::MCTruth>&    mcTruth    = truthPartAssn.first;
        const art::Ptr<simb::MCParticle>& mcParticle = truthPartAssn.second;
        
        fParticleList.Add(mcParticle.get());

        if (fTrackIDToMCTruthIndex.find(mcParticle->TrackId()) == fTrackIDToMCTruthIndex.end())
        {
            fMCTruthVec.emplace_back(mcTruth);
            fTrackIDToMCTruthIndex[mcParticle->TrackId()] = mcTruth;
        }
    }
    
    return;
}
    
const MCTruthParticleList& MCTruthAssociations::getParticleList() const
{
    return fParticleList;
}

// Return a pointer to the simb::MCParticle object corresponding to the given TrackID
const simb::MCParticle* MCTruthAssociations::TrackIDToParticle(int const& id) const
{
    // Pointer to return
    const simb::MCParticle* mcParticle(0);
    
    MCTruthParticleList::const_iterator partItr = fParticleList.find(id);
    
    if (partItr != fParticleList.end()) mcParticle = partItr->second;
    
    if(!mcParticle)
    {
        mf::LogWarning("MCTruthAssociations") << "can't find particle with track id "
        << id << " in MCTruthParticleList"
        << " returning null pointer";
    }
    
    return mcParticle;
}
    
const simb::MCParticle* MCTruthAssociations::TrackIDToMotherParticle(int const& id) const
{
    // get the mother id from the particle navigator
    // the EveId was adopted in the Rebuild method
    return this->TrackIDToParticle(fParticleList.EveId(abs(id)));
}

// Get art::Ptr<> to simb::MCTruth and related information
const art::Ptr<simb::MCTruth>& MCTruthAssociations::TrackIDToMCTruth(int const& id) const
{
    // find the entry in the MCTruth collection for this track id
    MCTruthTrackIDMap::const_iterator trackTruthItr = fTrackIDToMCTruthIndex.find(abs(id));
    
    if (trackTruthItr == fTrackIDToMCTruthIndex.end())
        throw cet::exception("MCTruthAssociations") << "attempting to find MCTruth index for "
        << "out of range value: " << id
        << "/" << fMCTruthVec.size() << "\n";
    
    return trackTruthItr->second;
}
    
const art::Ptr<simb::MCTruth>& MCTruthAssociations::ParticleToMCTruth(const simb::MCParticle* p) const
{
    return this->TrackIDToMCTruth(p->TrackId());
}
    
std::vector<const simb::MCParticle*> MCTruthAssociations::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    std::vector<const simb::MCParticle*> ret;
    
    // I'm slightly uncertain what is going on here since I naively would think a track ID is
    // unique to an MCParticle... but this is what was done previously so copied here...
    // sim::ParticleList::value_type is a pair (track ID, particle pointer)
    for (const MCTruthParticleList::value_type& TrackIDpair : fParticleList)
    {
        if (TrackIDToMCTruth(TrackIDpair.first) == mct) ret.push_back(TrackIDpair.second);
    }
    
    return ret;
}
    
const MCTruthTruthVec& MCTruthAssociations::MCTruthVector() const
{
    return fMCTruthVec;
}
    
// this method will return the Geant4 track IDs of
// the particles contributing ionization electrons to the identified hit
std::vector<sim::TrackIDE> MCTruthAssociations::HitToTrackID(const recob::Hit* hit) const
{
    std::vector<sim::TrackIDE> trackIDEs;
    
    // Apply brute force loop through all possible collections
    for(const auto& hitPartAssns : fHitPartAssnsVec)
    {
        HitToPartVecMap::const_iterator hitMatchPairItr = hitPartAssns.fHitToPartVecMap.find(hit);
    
        if (hitMatchPairItr != hitPartAssns.fHitToPartVecMap.end())
        {
            // float totalE(0.); // unused
        
            for (const auto& matchPair : hitMatchPairItr->second)
            {
                const simb::MCParticle*                 part = matchPair.first;
                const anab::BackTrackerHitMatchingData* data = matchPair.second;
                
                sim::TrackIDE info;
                info.trackID      = part->TrackId();
                info.energyFrac   = data->ideFraction;
                info.energy       = data->energy;
                info.numElectrons = data->numElectrons;
            
                // totalE += data->energy; // unused
            
                trackIDEs.push_back(info);
            }
            
            break;
        }
    }
    
    return trackIDEs;
}
    
std::vector<sim::TrackIDE> MCTruthAssociations::HitToTrackID(const art::Ptr<recob::Hit>& hit) const
{
    return HitToTrackID(hit.get());
}

//----------------------------------------------------------------------
const std::vector<std::vector<art::Ptr<recob::Hit>>> MCTruthAssociations::TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allHits,
                                                                                         std::vector<int> const&                  tkIDs) const
{
    // returns a subset of the hits in the allhits collection that are matched
    // to MC particles listed in tkIDs
    
    // temporary vector of TrackIDs and Ptrs to hits so only one
    // loop through the (possibly large) allhits collection is needed
    std::unordered_map<int, std::vector<art::Ptr<recob::Hit>>> trackIDHitVecMap;
    
    for(const auto& hit : allHits)
    {
        std::vector<sim::TrackIDE> trackIDEVec = this->HitToTrackID(hit.get());
        
        for(const auto trackIDE : trackIDEVec)
        {
            for(auto trackID : tkIDs)
            {
                if(trackID == trackIDE.trackID && trackIDE.energyFrac > fMinHitEnergyFraction)
                {
                    trackIDHitVecMap[trackID].push_back(hit);
                }
            }
        }
    }
    
    // now build the truHits vector that will be returned to the caller
    std::vector<std::vector<art::Ptr<recob::Hit>>> truHits;
    
    for(const auto& trackHitPair : trackIDHitVecMap) truHits.emplace_back(trackHitPair.second);
    
    return truHits;
}

//----------------------------------------------------------------------
// the particle list is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> MCTruthAssociations::HitToEveID(const recob::Hit* hit) const
{
    std::vector<sim::TrackIDE> trackIDEVec = this->HitToTrackID(hit);
    
    // Need the particle list, want the "biggest" one to make sure we have all the particles
    const MCTruthParticleList& particleList = getParticleList();

    // Create a map which will go between eve ids and TrackIDEs
    std::unordered_map<int, sim::TrackIDE> idToTrackIDEMap;
    
    for(const auto& trackID : trackIDEVec)
    {
        int eveTrackID  = particleList.EveId(trackID.trackID);
        
        std::unordered_map<int, sim::TrackIDE>::iterator eveElemItr = idToTrackIDEMap.find(eveTrackID);
        
        // If this id is already in the map we are incrementing values...
        if (eveElemItr != idToTrackIDEMap.end())
        {
            eveElemItr->second.energyFrac   += trackID.energyFrac;
            eveElemItr->second.energy       += trackID.energy;
            eveElemItr->second.numElectrons += trackID.numElectrons;
        }
        // Otherwise we are adding for the first time
        else
        {
            // First insertion, initialize to the values in the current trackID, then make sure we have eve track id
            idToTrackIDEMap[eveTrackID]         = trackID;
            idToTrackIDEMap[eveTrackID].trackID = eveTrackID;
        }
    }
    
    // Now create an output container and move info to it
    std::vector<sim::TrackIDE> eveTrackIDEVec;
    
    for(const auto& eveIDEPair : idToTrackIDEMap) eveTrackIDEVec.emplace_back(eveIDEPair.second);
    
    return eveTrackIDEVec;
}

// method to return the XYZ position of the weighted average energy deposition for a given hit
std::vector<double> MCTruthAssociations::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
{
    mf::LogWarning("MCTruthAssociations") << " ** HitToXYZ currently not implemented in MCTruthAssociations, return a zero point";
    return std::vector<double>() = {0.,0.,0.};
}

// method to return the XYZ position of a space point (unweighted average XYZ of component hits).
std::vector<double> MCTruthAssociations::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const
{
    mf::LogWarning("MCTruthAssociations") << " ** SpacePointHitsToXYZ currently not implemented in MCTruthAssociations, return a zero point";
    return std::vector<double>() = {0.,0.,0.};
}

// method to return the fraction of hits in a collection that come from the specified Geant4 track ids
double MCTruthAssociations::HitCollectionPurity(std::set<int>                              trackIDs,
                                                const std::vector< art::Ptr<recob::Hit> >& hitVec) const
{
    // get the list of EveIDs that correspond to the hits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = hitVec.size();;
    float desired = 0.;
    
    // don't have to check the plane in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct plane by definition then.
    for(const auto& hit : hitVec)
    {
        std::vector<sim::TrackIDE> hitTrackIDEVec = this->HitToTrackID(hit.get());
        
        for(const auto& trackIDE : hitTrackIDEVec)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end())
            {
                desired += 1.;
                break;
            }
        }
    }
    
    double purity = 0.;
    if(total > 0) purity = desired/total;
    
    return purity;
}

// method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
// represented in a collection of hits
double MCTruthAssociations::HitCollectionEfficiency(std::set<int>                              trackIDs,
                                                    const std::vector< art::Ptr<recob::Hit> >& hitVec,
                                                    const std::vector< art::Ptr<recob::Hit> >& allHitVec,
                                                    const geo::View_t&                         view) const
{
    // get the list of EveIDs that correspond to the hits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;
    
    // don't have to check the view in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(const auto& hit : hitVec)
    {
        std::vector<sim::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        // don't worry about hits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(const auto& trackIDE : hitTrackIDs)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end() &&
                trackIDE.energyFrac             >= fMinHitEnergyFraction)
            {
                desired += 1.;
                break;
            }
        }
    }
    
    // now figure out how many hits in the whole collection are associated with this id
    for(const auto& hit : allHitVec)
    {
        
        // check that we are looking at the appropriate view here
        // in the case of 3D objects we take all hits
        if(hit->View() != view && view != geo::k3D ) continue;
        
        std::vector<sim::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);

        for (const auto& trackIDE : hitTrackIDs)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end() &&
                trackIDE.energyFrac             >= fMinHitEnergyFraction)
            {
                total += 1.;
                break;
            }
        }
    }
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;
    
    return efficiency;
}

// method to return the fraction of charge in a collection that come from the specified Geant4 track ids
double MCTruthAssociations::HitChargeCollectionPurity(std::set<int>                              trackIDs,
                                                      const std::vector< art::Ptr<recob::Hit> >& hitVec) const
{
    // get the list of EveIDs that correspond to the hits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = 0;
    float desired = 0.;
    
    // don't have to check the view in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(const auto& hit : hitVec)
    {
        std::vector<sim::TrackIDE> hitTrackIDEVec = this->HitToTrackID(hit);
        
        total += hit->Integral();
        
        for(const auto& trackIDE : hitTrackIDEVec)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end())
            {
                desired += hit->Integral();
                break;
            }
        }

    }
    
    double purity = 0.;
    if(total > 0) purity = desired/total;
    
    return purity;
}

// method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are
// represented in a collection of hits
double MCTruthAssociations::HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                          const std::vector< art::Ptr<recob::Hit> >& hitVec,
                                                          const std::vector< art::Ptr<recob::Hit> >& allHitVec,
                                                          const geo::View_t&                         view) const
{
    // get the list of EveIDs that correspond to the hits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;
    
    // don't have to check the view in the hits collection because
    // those are assumed to be from the object we are testing and will
    // the correct view by definition then.
    for(const auto& hit : hitVec)
    {
        std::vector<sim::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        // don't worry about hits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(const auto& trackIDE : hitTrackIDs)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end() &&
                trackIDE.energyFrac             >= fMinHitEnergyFraction)
            {
                desired += hit->Integral();
                break;
            }
        }
    }
    
    // now figure out how many hits in the whole collection are associated with this id
    for(const auto& hit : allHitVec)
    {
        
        // check that we are looking at the appropriate view here
        // in the case of 3D objects we take all hits
        if(hit->View() != view && view != geo::k3D ) continue;
        
        std::vector<sim::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        for (const auto& trackIDE : hitTrackIDs)
        {
            if (trackIDs.find(trackIDE.trackID) != trackIDs.end() &&
                trackIDE.energyFrac             >= fMinHitEnergyFraction)
            {
                total += hit->Integral();
                break;
            }
        }
    }
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;
    
    return efficiency;
}

// method to return all EveIDs corresponding to the current sim::ParticleList
std::set<int> MCTruthAssociations::GetSetOfEveIDs() const
{
    std::set<int>              eveIDs;
    const MCTruthParticleList& particleList = getParticleList();
    
    for(const auto& pl : particleList) eveIDs.insert(particleList.EveId(pl.first));
    
    return eveIDs;
}

// method to return all TrackIDs corresponding to the current sim::ParticleList
std::set<int> MCTruthAssociations::GetSetOfTrackIDs() const
{
    // fParticleList::value_type is a pair (track, particle pointer)
    std::set<int>              trackIDs;
    const MCTruthParticleList& particleList = getParticleList();
    
    for (const auto& pl: particleList) trackIDs.insert(pl.first);
    
    return trackIDs;
}

// method to return all EveIDs corresponding to the given list of hits
std::set<int> MCTruthAssociations::GetSetOfEveIDs(const std::vector< art::Ptr<recob::Hit> >& hitVec) const
{
    std::set<int> eveIDs;
    
    for(const auto& hit : hitVec)
    {
        const std::vector<sim::TrackIDE> ideVec = HitToEveID(hit.get());
        
        // loop over the ides and extract the track ids
        for(const auto& trackIDE : ideVec) eveIDs.insert(trackIDE.trackID);
    }
    
    return eveIDs;
}

// method to return all TrackIDs corresponding to the given list of hits
std::set<int> MCTruthAssociations::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hitVec) const
{
    std::set<int> trackIDs;
    
    for(const auto& hit : hitVec)
    {
        std::vector<sim::TrackIDE> trackIDEVec = this->HitToTrackID(hit);
        
        for(const auto& trackIDE : trackIDEVec) trackIDs.insert(trackIDE.trackID);
    }
    
    return trackIDs;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double MCTruthAssociations::length(const recob::Track* track) const
{
  return track->Length();
}

// Length of MC particle.
//----------------------------------------------------------------------------
double MCTruthAssociations::length(const simb::MCParticle& part, double dx,
                                   TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                                   unsigned int tpc, unsigned int cstat) const
{
    // Need the readout window size...
    double readOutWindowSize = fDetectorProperties->ReadOutWindowSize();
    
    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;
    
    // Loop over the complete collection of trajectory points
    for(int i = 0; i < n; ++i)
    {
        TVector3 pos = part.Position(i).Vect();

        // Need to make sure this position is in an active region of the TPC
        // If the particle is not in the cryostat then we skip
        geo::TPCGeo const* tpcGeo = fGeometry->PositionToTPCptr(geo::vect::toPoint(pos));
        if (tpcGeo == nullptr) {
          continue;
        }
            
        auto const activePos = geo::vect::toPoint(pos) - tpcGeo->GetActiveVolumeCenter();
            
        if (std::fabs(activePos.X()) > tpcGeo->ActiveHalfWidth() or
            std::fabs(activePos.Y()) > tpcGeo->ActiveHalfHeight() or
            std::fabs(activePos.Z()) > 0.5 * tpcGeo->ActiveLength()) continue;
            
        double const xpos_in_tpc = activePos.X() + tpcGeo->ActiveHalfWidth();
        
        // Make fiducial cuts.
        // There are two sets here:
        // 1) We check the original x,y,z position of the trajectory points and require they be
        //    within the confines of the physical TPC
        // 2) We then check the timing of the presumed hit and make sure it lies within the
        //    readout window for this simulation
        pos[0] += dx;
        double const ticks = fDetectorProperties->ConvertXToTicks(xpos_in_tpc,  geo::PlaneID{tpcGeo->ID(), 0});

        // Currently it appears that the detector properties are not getting initialized properly and the returned
        // number of ticks is garbage so we need to skip this for now. Will be important when we get CR's
        if(ticks >= 0. && ticks < readOutWindowSize)
        {
            if(first)
            {
                start = pos;
                startmom = part.Momentum(i).Vect();
            }
            else
            {
                disp -= pos;
                result += disp.Mag();
            }
            first = false;
            disp = pos;
            end = pos;
            endmom = part.Momentum(i).Vect();
        }
    }
    
    return result;
}
    
} // end of namespace
