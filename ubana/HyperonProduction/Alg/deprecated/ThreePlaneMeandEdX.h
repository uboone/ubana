#ifndef _ThreePlaneMeandEdX_h_
#define _ThreePlaneMeandEdX_h_

#include "lardataobj/RecoBase/Track.h"

//get weight for plane
//based on code by Pip Hamilton in ubana/ParticleID/Modules/ParticleId_module.cc
double PlaneWeight( art::Ptr<recob::Track> track , int i_pl ){

   //Define wire plane angles
   double plane0_wireangle = 30*6.28/360.0;
   double plane1_wireangle = -30*6.28/360.0;
   double plane2_wireangle = 90*6.28/360.0;

   //Define weighting threshold
   double tophat_thresh = 0.175;


   //Find track angle in plane of the wires
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
   if (angle_planeweight < tophat_thresh) angle_planeweight = 0;
   if (angle_planeweight != 0) angle_planeweight = 1;


   return angle_planeweight;


}

/*
std::vector<std::pair<int,double>> MeandEdX(std::vector<art::Ptr<anab::Calorimetry>> Calo){

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
*/


//combine dEdX scores to make 3 plane dEdX

double ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<std::pair<int,double>> dEdXs){

   double dEdX=0;
   double Weight=0;

   for(size_t i_pl=0;i_pl<dEdXs.size();i_pl++){

      //catch default fills
      if( dEdXs.at(i_pl).first < 0 || dEdXs.at(i_pl).second < 0 ) continue;

      double thisPlaneWeight=PlaneWeight(track,dEdXs.at(i_pl).first);

      dEdX += dEdXs.at(i_pl).second*thisPlaneWeight;
      Weight += thisPlaneWeight;


   }

   if( Weight == 0 ) return -1;	

   return dEdX/Weight;
}


/*
//combine Chi2 scores to form three plane Chi2

double ThreePlaneChi2( art::Ptr<recob::Track> track, std::vector<std::pair<int,double>> Chi2s ){

   double Chi2=0;
   double Weight=0;


   for(size_t i_pl=0;i_pl<Chi2s.size();i_pl++){

      //catch default fills
      if( Chi2s.at(i_pl).first < 0 || Chi2s.at(i_pl).second < 0 ) continue;

      double thisPlaneWeight=PlaneWeight(track,Chi2s.at(i_pl).first);

      Chi2 += Chi2s.at(i_pl).second*thisPlaneWeight;
      Weight += thisPlaneWeight;


   }
   return Chi2/Weight;

}
*/

#endif
