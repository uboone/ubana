#ifndef _Track_Length_Cut_h_
#define _Track_Length_Cut_h_

#include "ubana/HyperonProduction/util/RecoParticle.h"

//Reject events with tracks exceeded a certain length
//(excluding muon candidate)


bool Track_Length_Cut(std::vector<RecoParticle> RecoParticles , int i_muon , double cut_leading=65 , double cut_subleading=35){

std::vector<double>lengths;

for(size_t i=0;i<RecoParticles.size();i++){

if(i == (size_t)i_muon) continue;

lengths.push_back(RecoParticles.at(i).TrackLength);

}

//place list of track lengths in ascending order
std::sort(lengths.begin(),lengths.end());

if(lengths.at(lengths.size()-1) > cut_leading) return false;
if(lengths.at(lengths.size()-2) > cut_subleading) return false;

return true;


}

#endif
