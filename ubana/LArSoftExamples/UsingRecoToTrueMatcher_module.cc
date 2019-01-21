////////////////////////////////////////////////////////////////////////
// Class:       UsingRecoToTrueMatcher
// Plugin Type: analyzer (art v3_00_00)
// File:        UsingRecoToTrueMatcher_module.cc
//
// Generated at Mon Jan 21 15:13:20 2019 by Adam Lister using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"           
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TTree.h"

class UsingRecoToTrueMatcher;


class UsingRecoToTrueMatcher : public art::EDAnalyzer {
public:
  explicit UsingRecoToTrueMatcher(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UsingRecoToTrueMatcher(UsingRecoToTrueMatcher const&) = delete;
  UsingRecoToTrueMatcher(UsingRecoToTrueMatcher&&) = delete;
  UsingRecoToTrueMatcher& operator=(UsingRecoToTrueMatcher const&) = delete;
  UsingRecoToTrueMatcher& operator=(UsingRecoToTrueMatcher&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  int run;
  int subRun;
  int event;

  std::string fTrackLabel;
  std::string fHitLabel;
  std::string fHitTrackAssnLabel;
  std::string fHitTruthAssnLabel;
};


UsingRecoToTrueMatcher::UsingRecoToTrueMatcher(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  
  fTrackLabel = p.get< std::string >("TrackLabel", "pandoraNu::McRecoStage2");
  fHitLabel   = p.get< std::string >("HitLabel", "pandoraCosmicHitRemoval::McRecoStage2");
  
  fHitTrackAssnLabel = p.get< std::string >("HitTrackAssnLabel", "pandoraNu::McRecoStage2");
  fHitTruthAssnLabel = p.get< std::string >("HitTruthAssnLabel", "pandoraCosmicHitRemoval::McRecoStage2");

}

void UsingRecoToTrueMatcher::analyze(art::Event const& e)
{

    run = e.run();
    subRun = e.subRun();
    event = e.event();

    std::cout << "[UsingRecoToTrueMatcher] Processing event " << 
        run << "." <<
        subRun << "." <<
        event << std::endl;

    // get handles to necessary information
    art::Handle< std::vector< recob::Track > > trackHandle;
    e.getByLabel(fTrackLabel, trackHandle);
    std::vector< art::Ptr< recob::Track > > trackPtrVector;
    art::fill_ptr_vector(trackPtrVector, trackHandle);

    art::Handle< std::vector< recob::Hit > > hitHandle;
    e.getByLabel(fHitLabel, hitHandle);

    art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fHitTrackAssnLabel); 
    art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit(hitHandle, e, fHitTruthAssnLabel);

    // loop tracks in event
    for (size_t i_track = 0; i_track < trackPtrVector.size(); i_track++){

        art::Ptr< recob::Track > thisTrack = trackPtrVector.at(i_track);

        // find associated hits
        std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(thisTrack.key());
    
        // 
        // Truth matching happens here!
        //

        std::unordered_map< int, double> trkide;
        double maxe = -1;
        double tote = 0;
        
        simb::MCParticle const* matchedParticle = NULL;
       
        std::vector< simb::MCParticle const*> particleVec;
        std::vector< anab::BackTrackerHitMatchingData const* > matchVec;

        for (size_t i_hit = 0; i_hit < hits.size(); ++i_hit){

            // clear vectors
            particleVec.clear();
            matchVec.clear();
            particlesPerHit.get(hits[i_hit].key(), particleVec, matchVec);

            // loop over particles which deposit energy in this hit
            for (size_t i_particle = 0; i_particle < particleVec.size(); ++i_particle){

                // store energy per track id
                trkide[ particleVec[i_particle]->TrackId() ] += matchVec[i_particle]->energy;

                // store total energy deposited
                tote += matchVec[i_particle]->energy;

                if ( trkide[ particleVec[i_particle]->TrackId() ] > maxe ){
                    // keep track of maximum
                    maxe = trkide[ particleVec[i_particle]->TrackId() ];
                    matchedParticle = particleVec[i_particle];
                }

            }

        }

        float purity = maxe/tote;

        std::cout << "[UsingRecoToTrueMatcher] Particle Matched to track with id " << thisTrack->ID() << " is MCParticle with PDG code " << matchedParticle->PdgCode() << std::endl; 
        std::cout << "[UsingRecoToTrueMatcher] Match purity/completeness: " << purity << std::endl;

    }

}

void UsingRecoToTrueMatcher::beginJob()
{
  // Implementation of optional member function here.
  
}

DEFINE_ART_MODULE(UsingRecoToTrueMatcher)
