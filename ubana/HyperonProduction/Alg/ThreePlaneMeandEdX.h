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


//combine dEdX scores to make 3 plane dEdX

double ThreePlaneMeandEdX( art::Ptr<recob::Track> track, std::vector<std::pair<int,double>> dEdXs ){


	double dEdX=0;
	double Weight=0;

	for(size_t i_pl=0;i_pl<dEdXs.size();i_pl++){

		//catch default fills
		if( dEdXs.at(i_pl).first < 0 || dEdXs.at(i_pl).second < 0 ) continue;

		double thisPlaneWeight=PlaneWeight(track,dEdXs.at(i_pl).first);

		dEdX += dEdXs.at(i_pl).second*thisPlaneWeight;
		Weight += thisPlaneWeight;


	}

	return dEdX/Weight;
}



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

#endif
