#ifndef _Muon_ID_h_
#define _Muon_ID_h_

#include "ubana/HyperonProduction/util/RecoParticle.h"

//choose muon candidate from list of reco'd particles
//return -1 if no suitable partilce found

int Muon_ID(std::vector<RecoParticle> RecoParticles , double LLR_PID_Cut=-0.35, double min_length=18){

int i_longest=-1;
int length=-1;



for(size_t i=0;i<RecoParticles.size();i++){

//if proton like or too far from reco PV or too short - skip!
if( RecoParticles.at(i).Track_LLR_PID < LLR_PID_Cut || RecoParticles.at(i).Displacement > 10 || RecoParticles.at(i).TrackLength < min_length) continue;

if(RecoParticles.at(i).TrackLength > length) {
i_longest = i;
length = RecoParticles.at(i).TrackLength;

}


}

return i_longest;

}

#endif
