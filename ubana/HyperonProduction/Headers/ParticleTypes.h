#ifndef _ParticleTypes_h_
#define _ParticleTypes_h_

inline bool isHyperon(int pdg){ return pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 3112; }

inline bool isPion(int pdg){ return pdg == 111 || abs(pdg) == 211; }

inline bool isNucleon(int pdg){ return pdg == 2112 || pdg == 2212; }

inline bool isLepton(int pdg){ return abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15; }

inline bool isNeutrino(int pdg){ return abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16; }

inline bool isKaon(int pdg){ return abs(pdg) == 321 || pdg == 311; }

#endif
