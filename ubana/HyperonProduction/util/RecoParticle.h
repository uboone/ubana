#ifndef _RecoParticle_h_
#define _RecoParticle_h_

#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class RecoParticle : public TObject{
#else
class RecoParticle {
#endif

public:

RecoParticle(){}
~RecoParticle(){}

int Index;

//general reco info
int PDG; //Pandora PDG code (11 or 13)
double TrackShowerScore;
double X,Y,Z;
double Displacement; //distance from RECO PV

//track info
double TrackLength;
double TrackDirectionX,TrackDirectionY,TrackDirectionZ;
double TrackStartX,TrackStartY,TrackStartZ;
double TrackEndX,TrackEndY,TrackEndZ;
double TrackMuonClosestApproachPosition,TrackMuonClosestApproachDistance;


//track PID info
//3 plane PID score
double TrackPID; 
//Mean dE/dX scores
double MeandEdX_Plane0,MeandEdX_Plane1,MeandEdX_Plane2,MeandEdX_ThreePlane;
//Nicolo's PID
double Track_LLR_PID;

//track kinematics
double ProtonMomentum,MuonMomentum;


//truth info
bool HasTruth; //false if reco particle has no corresponding MC particle (eg if its a cosmic overlay!)
int TrackTruePDG;
double TrackTrueE,TrackTruePx,TrackTruePy,TrackTruePz;
double TrackTrueModMomentum;
double TrackTrueKE;
double TrackTrueLength;
int TrackTrueOrigin; // 1 - primary , 2 - hyperon decay, 3 - other
double TrackTrueTotalEdep;
double TrackEdepPurity;



void SetVertex(TVector3 V);
void SetTrackPositions(TVector3 Start,TVector3 End);
void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(RecoParticle,1);
#endif


};


void RecoParticle::SetVertex(TVector3 V){

X = V.X();
Y = V.Y();
Z = V.Z();

}

void RecoParticle::SetTrackPositions(TVector3 Start,TVector3 End){

TrackStartX = Start.X();
TrackStartY = Start.Y();
TrackStartZ = Start.Z();

TrackEndX = End.X();
TrackEndY = End.Y();
TrackEndZ = End.Z();

}



void RecoParticle::Print(){

std::cout << "Reco Info:" << std::endl;
std::cout << "PDG Code: " << PDG << "  Track/Shower score: " << TrackShowerScore << std::endl;
std::cout << "Track length: " << TrackLength << "  PID score: " << TrackPID <<  std::endl;
std::cout << "Truth Info:" << std::endl;
std::cout << "PDG: " << TrackTruePDG << "  Origin: " << TrackTrueOrigin << std::endl;
std::cout << "Length: " << TrackTrueLength << "  KE: " << TrackTrueKE << std::endl;


}




#endif
