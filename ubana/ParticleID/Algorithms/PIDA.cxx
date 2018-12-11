/**
 * Re-implementation of PIDA algorithm originally developed by B. Baller.
 */

#include "PIDA.h"

namespace particleid{


  double PIDA::getPida(std::vector<float> dEdx, std::vector<float> resRange, std::string method){


    kde::KernelDensityEstimator kde;
    double pida = -1.0;

    std::vector<double> pidaValues;
    if (resRange.size() == 0 || resRange.size() == 1) return -1;
    for (size_t i = 1; i < resRange.size()-1; i++){
      // Don't take resRange values higher than 30cm
      if (resRange.at(i)>30) continue;
      
      pidaValues.push_back(dEdx.at(i)*std::pow(resRange.at(i),0.42));

    }

    if (method == "mean")
      pida = TMath::Mean(pidaValues.begin(), pidaValues.end());

    else if (method == "median")
      pida = TMath::Median(pidaValues.size(), &pidaValues[0]);

    else if (method == "kde")
      pida = kde.getKernelDensityMpv(pidaValues);

    return pida;

  }

}
