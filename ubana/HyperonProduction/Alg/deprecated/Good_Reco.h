#ifndef _Good_Reco_h_
#define _Good_Reco_h_

#include "ubana/HyperonProduction/util/RecoParticle.h"

//goes through list of reconstructed particles, checks if there is a proton and pion
//from the decay reconstructed

bool Good_Reco(std::vector<RecoParticle> Particles){

bool got_proton=false;
bool got_pion=false;

for(size_t i=0;i<Particles.size();i++){

if(Particles.at(i).TrackTruePDG == 2212 && Particles.at(i).TrackTrueOrigin == 2) got_proton = true;
if(Particles.at(i).TrackTruePDG == -211 && Particles.at(i).TrackTrueOrigin == 2) got_pion = true;

}

return (got_proton && got_pion);

}


#endif
