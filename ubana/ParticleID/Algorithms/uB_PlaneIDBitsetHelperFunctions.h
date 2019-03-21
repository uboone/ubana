// MicroBooNE-specific helper functions to translate bitset planeIDs used in anab::ParticleID

#ifndef _UB_PLANEIDBITSETHELPER_H_
#define _UB_PLANEIDBITSETHELPER_H_

#include <bitset>

namespace UBPID{

  /// Returns planeID from bitset, if exactly 1 plane is set in the bitset.
  /// Plane 2 = collection (y) plane, Plane 1 = induction (v) plane, Plane 0 = induction (u) plane
  inline int uB_getSinglePlane(std::bitset<8> planeid){
    if (planeid.none() || planeid.count() > 1 || (planeid.count()==1 && !(planeid.test(0) || planeid.test(1) || planeid.test(2)))){
      std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a single MicroBooNE plane for bitset " << planeid << ". Returning -1 (invalid planeID)." << std::endl;
      return -1;
    }
    else if (planeid.test(0)) return 2;
    else if (planeid.test(1)) return 1;
    else if (planeid.test(2)) return 0;

    // Default: invalid return
    return -1;
  }

  /// Returns the bitset corresponding to a given (single) uBooNE plane.
  inline std::bitset<8> uB_SinglePlaneGetBitset(int planeid){
    std::bitset<8> bs;
    if (planeid<0 || planeid>2){
      std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a bitset for MicroBooNE planeid " << planeid << ". Returning invalid planeID." << std::endl;
      return bs;
    }
    else bs.set(2-planeid); // minus 2 because that is the collection plane

    // Default: invalid return
    return bs;
  }

}

#endif
