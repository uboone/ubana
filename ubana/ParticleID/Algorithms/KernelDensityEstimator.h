#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "TF1.h"
#include "TFormula.h"
#include "TString.h"

#include <vector>
#include <string>

namespace kde{

  class KernelDensityEstimator{

    public:

      static Double_t epoKernelFunction(Double_t *x, Double_t *p);

      static Double_t gausKernelFunction(Double_t *x, Double_t *p);

      double getLocalDensity(std::vector<TF1*> kernels, int testpoint, double pilotBandwith);

      TF1* getKernel(double normalisation, double kernelMean, double bandwith, std::string kernelType);

      double getKernelDensityMpv(std::vector<double> pidaVals);

  };

}


#endif
