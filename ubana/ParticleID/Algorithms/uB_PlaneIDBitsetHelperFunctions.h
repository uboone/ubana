// MicroBooNE-specific helper functions to translate bitset planeIDs used in anab::ParticleID

#ifndef _UB_PLANEIDBITSETHELPER_H_
#define _UB_PLANEIDBITSETHELPER_H_

#include <bitset>

namespace UBPID{

  /// Returns planeID from bitset, if exactly 1 plane is set in the bitset.
  /// Plane 2 = collection (y) plane, Plane 1 = induction (v) plane, Plane 0 = induction (u) plane
  inline int uB_getSinglePlane(std::bitset<5> planeid){
    if ((planeid[0]+planeid[1]+planeid[2]+planeid[3]+planeid[4]>1) || (planeid[0]+planeid[1]+planeid[2]+planeid[3]+planeid[4]==0) || planeid[0]==1 || planeid[1]==1){
      std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a single MicroBooNE plane for bitset " << planeid << ". Returning -1 (invalid planeID)." << std::endl;
      return -1;
    }
    else if (planeid[4]==1) return 2;
    else if (planeid[3]==1) return 1;
    else if (planeid[2]==1) return 0;

    // Default: invalid return
    return -1;
  }

  /// Returns the bitset corresponding to a given (single) uBooNE plane.
  inline std::bitset<5> uB_SinglePlaneGetBitset(int planeid){
    if (planeid<0 || planeid>2){
      std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a bitset for MicroBooNE planeid " << planeid << ". Returning 00000 (invalid planeID)." << std::endl;
      return std::bitset<5>("00000");
    }
    else if (planeid==2) return std::bitset<5>("10000");
    else if (planeid==1) return std::bitset<5>("01000");
    else if (planeid==0) return std::bitset<5>("00100");

    // Default: invalid return
    return std::bitset<5>("00000");
  }

}

#endif
