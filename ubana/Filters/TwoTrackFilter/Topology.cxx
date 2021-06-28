// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef TOPOLOGY_CXX
#define TOPOLOGY_CXX

#include "Topology.h"
// _________________________________________________________________________________________________________________________________________________________________________________________________

int Topology::TopologyLabel(int Nmuons, int Nelectrons, int Npiplus_abTH, int Npiplus_blTH, int Npiminus_abTH, int Npiminus_blTH, int Npi0, int Nprotons_abTH, int Nprotons_blTH, int PDG, int CCNC, bool ifbeam, bool ifFV){
  // In the current version, Proton Momentum threshold is set to be 300MeV, which should be reviewed soon
  // 1. NuMuCC0pi0p in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 0 && Npiminus_abTH == 0 && Nprotons_abTH == 0 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 1;
  // 2. NuMuCC0pi1p in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 0 && Npiminus_abTH == 0 && Nprotons_abTH == 1 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 2;
  // 3. NuMuCC0pi2p in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 0 && Npiminus_abTH == 0 && Nprotons_abTH == 2 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 3;
  // 4. NuMuCC0piNp in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 0 && Npiminus_abTH == 0 && Nprotons_abTH > 2 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 4;
  // 5. NuMuCC1pi+Xp in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 1 && Npiminus_abTH == 0 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 5;
  // 6. NuMuCC1pi-Xp in FV
  if (Nmuons > 0 && Npi0 == 0 && Npiplus_abTH == 0 && Npiminus_abTH == 1 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 6;
  // 7. NuMuCC1pi0Xp in FV
  if (Nmuons > 0 && Npi0 == 1 && Npiplus_abTH == 0 && Npiminus_abTH == 0 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 7;
  // 8. NuMuCCNpiXp in FV
  if (Nmuons > 0 && (Npi0 + Npiplus_abTH + Npiminus_abTH) > 1 && ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == 14) return 8;
  // 9. Anti NuMu CC in FV
  if (ifbeam == 1 && ifFV == 1 && CCNC == 0 && PDG == -14) return 9;
  // 10. Nue / Anti-Nue CC in FV
  if (ifbeam == 1 && ifFV == 1 && CCNC == 0 && abs(PDG) == 12) return 10;
  // 11. NC
  if (ifbeam == 1 && ifFV == 1 && CCNC == 1) return 11;
  // 12. true neutrino vertex out of FV
  if (ifbeam == 1 && ifFV == 0) return 12;
  // 13. cosmic
  if (ifbeam == 0) return 13;
  // 14. other
  else return 14;
}


#endif
