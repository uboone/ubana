#ifndef _RecoParticle_h_
#define _RecoParticle_h_

#include <iostream>

#include "TLorentzVector.h"
#include "TVector3.h"

#include "ubana/HyperonProduction/Headers/FV.h"

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

// General reco info
int PDG; // Pandora PDG code (11 or 13)
double TrackShowerScore;
double X,Y,Z;
double Displacement; // Distance from PV

// Track variables
double TrackLength=0;
double TrackDirectionX=0,TrackDirectionY=0,TrackDirectionZ=0;
double TrackStartX=0,TrackStartY=0,TrackStartZ=0;
double TrackEndX=0,TrackEndY=0,TrackEndZ=0;
double TrackContained;
double TrackPID; // 3 plane PID score
double MeandEdX_Plane0,MeandEdX_Plane1,MeandEdX_Plane2,MeandEdX_ThreePlane; // Mean dE/dX scores
double Track_LLR_PID; // LLR PID
double ProtonMomentum,MuonMomentum; // Track kinematics

// Truth info
bool HasTruth; // False if reco particle has no corresponding MC particle
int TrackTruePDG;
double TrackTrueE,TrackTruePx,TrackTruePy,TrackTruePz;
double TrackTrueModMomentum;
double TrackTrueKE;
double TrackTrueLength;
int TrackTrueOrigin; // 1 - primary , 2 - hyperon decay, 3 - other, 4 - kaon decay, 5 - Sigma0 decay
double TrackTruthPurity;

inline void SetVertex(TVector3 V);
inline void SetTrackPositions(TVector3 Start,TVector3 End);
inline void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(RecoParticle,1);
#endif

};


inline void RecoParticle::SetVertex(TVector3 V){

X = V.X();
Y = V.Y();
Z = V.Z();

}

inline void RecoParticle::SetTrackPositions(TVector3 Start,TVector3 End){

TrackStartX = Start.X();
TrackStartY = Start.Y();
TrackStartZ = Start.Z();

TrackEndX = End.X();
TrackEndY = End.Y();
TrackEndZ = End.Z();


// Set containment
TrackContained = inActiveTPC(TVector3(TrackEndX,TrackEndY,TrackEndZ));

}



inline void RecoParticle::Print(){

std::cout << "Reco Info:" << std::endl;
std::cout << "PDG Code: " << PDG << "  Track/Shower score: " << TrackShowerScore << std::endl;
std::cout << "Track length: " << TrackLength << "  PID score: " << TrackPID <<  std::endl;
std::cout << "Truth Info:" << std::endl;
std::cout << "PDG: " << TrackTruePDG << "  Origin: " << TrackTrueOrigin << std::endl;
std::cout << "Length: " << TrackTrueLength << "  KE: " << TrackTrueKE << std::endl;

}




#endif
