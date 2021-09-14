#ifndef _Gap_Cut_h_
#define _Gap_Cut_h_

#include "ubana/HyperonProduction/util/RecoParticle.h"

//reject events with tracks exceeded a certain length
//(excluding muon candidate)

bool Gap_Cut(std::vector<RecoParticle> RecoParticles , int i_muon, double min_gap=0, double max_gap=80){


//get position of muon candidate
double muon_x = RecoParticles.at(i_muon).X;
double muon_y = RecoParticles.at(i_muon).Y;
double muon_z = RecoParticles.at(i_muon).Z;

std::vector<int> accepted;

for(size_t i=0;i<RecoParticles.size();i++){

if(i == (size_t)i_muon) continue;

double d = sqrt( (RecoParticles.at(i).X - muon_x)*(RecoParticles.at(i).X - muon_x) 
	  + (RecoParticles.at(i).Y - muon_y)*(RecoParticles.at(i).Y - muon_y)
	  + (RecoParticles.at(i).Z - muon_z)*(RecoParticles.at(i).Z - muon_z) );

if(d < max_gap && d > min_gap) accepted.push_back(i); 

}


if(accepted.size() != 2) return false;

//check relative positions of accepted tracks

RecoParticle Accepted1 = RecoParticles.at(accepted.at(0));
RecoParticle Accepted2 = RecoParticles.at(accepted.at(1));

double d = sqrt( (Accepted1.X - Accepted2.X)*(Accepted1.X - Accepted2.X)
	       + (Accepted1.Y - Accepted2.Y)*(Accepted1.Y - Accepted2.Y)
	       + (Accepted1.Z - Accepted2.Z)*(Accepted1.Z - Accepted2.Z) );
 
//std::cout << d << std::endl;

if(d > 1) return false;

return true;


}

#endif
