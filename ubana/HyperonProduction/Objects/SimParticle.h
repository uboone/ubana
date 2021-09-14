#ifndef _SimParticle_h_
#define _SimParticle_h_

#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class SimParticle : public TObject{
#else
class SimParticle {
#endif

public:

SimParticle() {}
~SimParticle() {}

int PDG=0;
double E=0,Px=0,Py=0,Pz=0;
double ModMomentum=0;
double KE=0;
double StartX=0,StartY=0,StartZ=0;
double EndX=0,EndY=0,EndZ=0;
double Travel=0; // Distance travelled
double Theta=0;
double Phi=0;
int Origin=0; //0 = neutrino , 1 = primary , 2 = hyperon decay, 3 = other, 4 = Kaon decay, 5 = Sigma0 decay

inline void SetKinematics(TLorentzVector P, double Mass);
inline void SetPositions(TLorentzVector Start, TLorentzVector End);
inline void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(SimParticle,1);
#endif

};

inline void SimParticle::SetKinematics(TLorentzVector P, double Mass){

//Momentum = P;
E = P.E();
Px = P.X();
Py = P.Y();
Pz = P.Z();
ModMomentum = sqrt(Px*Px + Py*Py + Pz*Pz);
KE = E - Mass;

Theta = (180/3.1416)*TMath::ACos(Pz/ModMomentum);
Phi = (180/3.1416)*TMath::ASin(Py/sqrt(Px*Px+Py*Py));

}

inline void SimParticle::SetPositions(TLorentzVector Start, TLorentzVector End){ 

//Vertex = {Start.X() , Start.Y() , Start.Z()};	
StartX = Start.X();
StartY = Start.Y();
StartZ = Start.Z();

EndX = End.X();
EndY = End.Y();
EndZ = End.Z();

Travel = sqrt( (StartX-EndX)*(StartX-EndX) + (StartY-EndY)*(StartY-EndY) + (StartZ-EndZ)*(StartZ-EndZ) );


}

inline void SimParticle::Print(){

std::cout << "PDG: " << PDG << "  Origin: " << Origin << std::endl;
std::cout << "Length: " << Travel << "  KE: " << KE << std::endl;

}

#endif
