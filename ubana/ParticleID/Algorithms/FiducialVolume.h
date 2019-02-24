#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"

namespace fidvol{

  class fiducialVolume{

    public:
      
      std::vector<double> setFiducialVolume(std::vector<double> fv, fhicl::ParameterSet const & p);

      void printFiducialVolume(std::vector<double> fv);

      bool isInFiducialVolume(TVector3 xyz, std::vector<double> fv);

      std::vector<double> getTpcDimensions();

  };

}

#endif

