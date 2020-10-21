
#include "ubana/HyperonProduction/util/RecoParticle.h"

//choose muon candidate from list of reco'd particles
//return -1 if no suitable partilce found

int Muon_ID(std::vector<RecoParticle> RecoParticles , double Three_Plane_PID_Cut=0.6){

int i_longest=-1;
int length=-1;

for(size_t i=0;i<RecoParticles.size();i++){

//if proton like or too far from reco PV - skip!
if( RecoParticles.at(i).TrackPID > Three_Plane_PID_Cut || RecoParticles.at(i).Displacement > 10) continue;

if(RecoParticles.at(i).TrackLength > length) {
i_longest = i;
length = RecoParticles.at(i).TrackLength;

}


}

return i_longest;

}


