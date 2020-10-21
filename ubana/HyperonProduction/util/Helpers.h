#include "SimParticle.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "TVector3.h"
#include "nusimdata/SimulationBase/MCParticle.h"


// Helper function for constructing SimParticle

SimParticle MakeSimParticle(art::Ptr<simb::MCParticle> Part){

SimParticle S;

S.SetKinematics(Part->Momentum(),Part->Mass());
S.PDG = Part->PdgCode();

S.SetPositions(Part->Position(),Part->EndPosition());

return S;



}

