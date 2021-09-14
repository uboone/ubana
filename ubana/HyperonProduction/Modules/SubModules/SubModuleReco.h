#ifndef _SubModuleReco_h_
#define _SubModuleReco_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Headers/LLR_PID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"
#include "ubana/HyperonProduction/Alg/MeandEdXCalculator.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"


#include "TVector3.h"

using std::string;

namespace hyperon {

struct RecoData {

   TVector3 RecoPrimaryVertex = TVector3(-1000,-1000,-1000);

   int NPrimaryDaughters; 
   int NPrimaryTrackDaughters;
   int NPrimaryShowerDaughters;

   std::vector<RecoParticle> TrackPrimaryDaughters;
   std::vector<RecoParticle> ShowerPrimaryDaughters;

   std::vector<TVector3> TrackStarts;

   size_t TrueMuonIndex = -1;
   size_t TrueDecayProtonIndex = -1;
   size_t TrueDecayPionIndex = -1;

   bool GoodReco = false;

};


class SubModuleReco {

   public:

      SubModuleReco();
      SubModuleReco(art::Event const& e,bool isdata,string pfparticlelabel,string tracklabel,
                        string showerlabel,string vertexlabel,string pidlabel,string calolabel,string hitlabel,
                        string hittruthassnlabel,string trackhitassnlabel,string metadatalabel,string g4label);

      SubModuleReco(art::Event const& e,bool isdata,fhicl::ParameterSet pset);

      void PrepareInfo(); 
      TVector3 GetPrimaryVertex();
      void SetIndices(bool IsSignal=false);

      RecoData GetInfo();

   private:

      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;

      art::Handle<std::vector<recob::Track>> Handle_Track;
      std::vector<art::Ptr<recob::Track>> Vect_Track;

      art::Handle<std::vector<recob::Shower>> Handle_Shower;
      std::vector<art::Ptr<recob::Shower>> Vect_Shower;

      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;

      RecoParticle MakeRecoParticle(art::Ptr<recob::PFParticle> pfp);

      art::FindManyP<recob::Vertex>* Assoc_PFParticleVertex;
      art::FindManyP<recob::Track>* Assoc_PFParticleTrack;
      art::FindManyP<recob::Shower>* Assoc_PFParticleShower;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* Assoc_PFParticleMetadata;
      art::FindManyP<recob::Hit>* Assoc_TrackHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* Assoc_MCParticleBacktracker;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID;

      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
      SubModuleG4Truth* G4T = nullptr;

      MeandEdXCalculator dEdXCalc;

      RecoData theData;
      size_t neutrinoID = 99999;

      void GetPFPMetadata(art::Ptr<recob::PFParticle> pfp,RecoParticle &P);
      void GetTrackData(art::Ptr<recob::PFParticle> pfp,RecoParticle &P);
      void TruthMatch(art::Ptr<recob::Track> trk,RecoParticle &P);
      void GetPIDs(art::Ptr<recob::Track> trk,RecoParticle &P);
      void GetVertexData(art::Ptr<recob::PFParticle> pfp,RecoParticle &P);

      bool IsData;

};

}

#endif


