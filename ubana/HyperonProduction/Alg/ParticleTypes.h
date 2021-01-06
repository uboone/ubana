#ifndef _ParticleTypes_h_
#define _ParticleTypes_h_

bool isHyperon(int pdg){ 

	if(pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 3112) return true;

	else 
		return false;
}


bool isPion(int pdg){

	if(pdg == 111 || abs(pdg) == 211) return true;

	else 
		return false;
}


bool isNucleon(int pdg){

	if(pdg == 2112 || pdg == 2212) return true;

	else return false;
}


bool isLepton(int pdg){

	if(abs(pdg) == 11 || abs(pdg) == 13) return true;

	else return false;


}

bool isNeutrino(int pdg){

	if(abs(pdg) == 12 || abs(pdg) == 14) return true;

	else return false;


}

bool isKaon(int pdg){

if(abs(pdg) == 321 || pdg == 311) return true;

else return false;


}


#endif
