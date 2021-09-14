#ifndef _Helpers_h_
#define _Helpers_h_

//local includes
#include "SimParticle.h"
#include "RecoParticle.h"

//larsoft objects
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

//uboonecode objects
#include "nusimdata/SimulationBase/MCParticle.h"

//root objects
#include "TVector3.h"

// Helper functions used to transform larsoft objects into SimParticle and RecoParticle

inline SimParticle MakeSimParticle(simb::MCParticle Part){

   SimParticle S;

   S.SetKinematics(Part.Momentum(),Part.Mass());
   S.PDG = Part.PdgCode();

   S.SetPositions(Part.Position(),Part.EndPosition());

   return S;

}


// Helper function for setting track variables in Reco Particle
inline void SetTrackVariables(RecoParticle &P , art::Ptr<recob::Track> trk){

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   P.TrackLength = trk->Length();

   P.TrackDirectionX = trk->StartDirection().X();
   P.TrackDirectionY = trk->StartDirection().Y();
   P.TrackDirectionZ = trk->StartDirection().Z();

   trkf::TrackMomentumCalculator trkm{0};

   P.ProtonMomentum = trkm.GetTrackMomentum(trk->Length(),2212);	
   P.MuonMomentum = trkm.GetTrackMomentum(trk->Length(),13);

   geo::Point_t point = { trk->Start().X() ,  trk->Start().Y() ,  trk->Start().Z() };
   geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

   // Offset track start

   //sce - forward
   TVector3 Start( trk->Start().X() + sce_corr.X() , trk->Start().Y() - sce_corr.Y() , trk->Start().Z() - sce_corr.Z() );

   point = { trk->End().X() ,  trk->End().Y() ,  trk->End().Z() };
   sce_corr = SCE->GetPosOffsets(point);

   //sce - forward
   TVector3 End( trk->End().X() + sce_corr.X() , trk->End().Y() - sce_corr.Y() , trk->End().Z() - sce_corr.Z() );

   P.SetTrackPositions( Start , End );


}

#endif
