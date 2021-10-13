#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "RecoTruthMCParticle.h"

#ifndef RECOTRUTHMCPARTICLE_CXX
#define RECOTRUTHMCPARTICLE_CXX

void RecoTruth::TrackToMCParticle(art::Event const & evt, const art::Handle<std::vector<recob::Hit> >& hit_handle, std::string track_producer){
//void RecoTruth::TrackToMCParticle(art::Event const & evt, art::Ptr<recob::Track track, simb::MCParticle MCparticle){
  art::InputTag MCParticleModuleLabel { "largeant" };
  art::InputTag BackTrackerLabel { "gaushitTruthMatch" };
  art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,evt,BackTrackerLabel);

  std::vector<art::Ptr<simb::MCParticle>> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  //particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);

  // Get hit associated to the corresponding tracks by lar_pandora
  lar_pandora::TrackVector trackVec;
  lar_pandora::TracksToHits hitsIntrack;
  lar_pandora::LArPandoraHelper::CollectTracks(evt, track_producer, trackVec, hitsIntrack);
  std::cout<<"track size: "<<trackVec.size()<<" / " <<hitsIntrack.size()<<std::endl;
  for(int n = 0; n < (int) trackVec.size(); n++){
    std::cout<<"hit in track size: "<<hitsIntrack.at(trackVec[n]).size()<<std::endl;
    //particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
  }
}


#endif
