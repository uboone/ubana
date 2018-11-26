#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"

namespace fidvol{

  class fiducialVolume{

    public:
      
      std::vector<float> setFiducialVolume(std::vector<float> fv, fhicl::ParameterSet const & p);

      void printFiducialVolume(std::vector<float> fv);

      bool isInFiducialVolume(TVector3 xyz, std::vector<float> fv);

      std::vector<float> getTpcDimensions();

  };

}

#endif

