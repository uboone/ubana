#ifndef _Position_Match_h_
#define _Position_Match_h_

#include "TVector3.h"
#include "lardataobj/RecoBase/Track.h"

//calculates the distance of closest approach between point Position and track trk
//option to exclude the first b cm of the track from calculation

std::pair<double,double> Track_Closest_Approach(art::Ptr<recob::Track> trk, TVector3 Position ){
	
//distance of closest approach
double d_min=1e10;

//distance along track the approach ocurs
double x_min=-1;

//go through the track trajectory points,check track vertex separation
for(size_t i_p=0;i_p<trk->NumberTrajectoryPoints();i_p++){

//skip default fills
if(trk->LocationAtPoint(i_p).X() == -999 || trk->LocationAtPoint(i_p).Y() == -999 || trk->LocationAtPoint(i_p).Z() == -999) continue;

//calculate distance from start
double x = trk->Length() - trk->Length(i_p);

//check proximity of Position
TVector3 TrackPosition(trk->LocationAtPoint(i_p).X(),trk->LocationAtPoint(i_p).Y(),trk->LocationAtPoint(i_p).Z());

//just check distance between position and each traj point - should be ok so long
//as distance between traj points << dmax
double d = sqrt( (Position - TrackPosition)*(Position - TrackPosition) );

if(d<d_min) { d_min = d; x_min = x;}

}

return std::make_pair(x_min,d_min);
}

#endif
