#ifndef _MeandEdX_h_
#define _MeandEdX_h_

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "TVector3.h"

//calculate average dEdX  for a track
std::vector<std::pair<int,double>> MeandEdX( std::vector<art::Ptr<anab::Calorimetry>> Calo ){

std::vector<std::pair<int,double>> dEdX;

for(size_t i_plane = 0;i_plane < Calo.size();i_plane++){

art::Ptr<anab::Calorimetry> thisPlaneCalo = Calo.at(i_plane);

//go through each trajectory point, calculate the length of the step and the E dep in it

double totalE=0;
double totalX=0;

//sometimes this vector is empty, causes crash below, skip plane if it is
if( thisPlaneCalo->XYZ().size() < 2 ){
dEdX.push_back( std::make_pair(thisPlaneCalo->PlaneID().Plane,-1) );
continue;
}


for(size_t i_point = 0;i_point < thisPlaneCalo->XYZ().size()-1;i_point++){

//this point
anab::Point_t thisPos = thisPlaneCalo->XYZ().at(i_point);

//next point
anab::Point_t nextPos = thisPlaneCalo->XYZ().at(i_point+1);

//step vector
TVector3 D( thisPos.X() - nextPos.X() , thisPos.X() - nextPos.X() , thisPos.X() - nextPos.X() );

totalX += D.Mag();
totalE += thisPlaneCalo->dEdx().at(i_point)*D.Mag();

}

dEdX.push_back( std::make_pair(thisPlaneCalo->PlaneID().Plane,totalE/totalX) );

}

return dEdX;
}



#endif
