#ifndef _MeandEdXCalculator_cxx_
#define _MeandEdXCalculator_cxx_

#include "MeandEdXCalculator.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

MeandEdXCalculator::MeandEdXCalculator(){}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MeandEdXCalculator::SetTophatThresh(double thresh){

   TophatThresh = thresh;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Single plane Mean dE/dX

double MeandEdXCalculator::GetMeandEdX(art::Ptr<anab::Calorimetry> calo){

   double totalE=0;
   double totalX=0;

   // Sometimes this vector is empty, causes crash below, skip plane if it is
   if(calo->XYZ().size() < 2) return -1;

   for(size_t i_point = 0;i_point < calo->XYZ().size()-1;i_point++){

      anab::Point_t thisPos = calo->XYZ().at(i_point);
      anab::Point_t nextPos = calo->XYZ().at(i_point+1);

      // Step vector
      TVector3 D(thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X());

      totalX += D.Mag();
      totalE += calo->dEdx().at(i_point)*D.Mag();

   }

   return totalE/totalX;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MeandEdXCalculator::PlaneWeight(art::Ptr<recob::Track> track,int i_pl){

   // Find track angle in plane of the wires
   double y = track->End().y() - track->Start().y();
   double z = track->End().z() - track->Start().z();

   TVector3 trackvec(0, y, z);
   trackvec = trackvec.Unit();
   TVector3 zaxis(0, 0, 1);
   double costhetayz = trackvec.Dot(zaxis);
   double thetayz = TMath::ACos(costhetayz);
   if ((y < 0) && (thetayz > 0)) thetayz *= -1;

   double theta_towires = 0;
   if (i_pl == 0) theta_towires = std::min(std::abs(plane0_wireangle - thetayz), std::abs((-1*(6.28-plane0_wireangle) - thetayz)));
   if (i_pl == 1) theta_towires = std::min(std::abs(plane1_wireangle - thetayz), std::abs((-1*(6.28-plane1_wireangle) - thetayz)));
   if (i_pl == 2) theta_towires = std::min(std::abs(plane2_wireangle - thetayz), std::abs((-1*(6.28-plane2_wireangle) - thetayz)));

   double angle_planeweight = sin(theta_towires)*sin(theta_towires);
   if (angle_planeweight < TophatThresh) angle_planeweight = 0;
   if (angle_planeweight != 0) angle_planeweight = 1;

   return angle_planeweight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

dEdXStore MeandEdXCalculator::ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v){

   dEdXStore theStore;

   double TotaldEdX=0;
   double TotalWeight=0;

   for(size_t i_pl=0;i_pl<calo_v.size();i_pl++){

      int plane = calo_v.at(i_pl)->PlaneID().Plane;
        
      if(plane != 0 && plane != 1 && plane != 2) continue;        

      double dEdX = GetMeandEdX(calo_v.at(i_pl));

      // Catch default fills
      if(dEdX < 0) continue;

      double thisPlaneWeight = PlaneWeight(track,plane);

      if(plane == 0){
         theStore.Weight_Plane0 = thisPlaneWeight;       
         theStore.Plane0 = dEdX;
      }
      if(plane == 1){
         theStore.Weight_Plane1 = thisPlaneWeight;       
         theStore.Plane1 = dEdX;
      }
      if(plane == 2){
         theStore.Weight_Plane2 = thisPlaneWeight;       
         theStore.Plane2 = dEdX;
      }

      TotaldEdX += dEdX*thisPlaneWeight;
      TotalWeight += thisPlaneWeight;
   }

   if(TotalWeight > 0) theStore.ThreePlaneAverage = TotaldEdX/TotalWeight;

return theStore;

}


#endif
