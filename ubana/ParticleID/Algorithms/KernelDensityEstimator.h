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

      float getLocalDensity(std::vector<TF1*> kernels, int testpoint, float pilotBandwith);

      TF1* getKernel(float normalisation, float kernelMean, float bandwith, std::string kernelType);

      float getKernelDensityMpv(std::vector<float> pidaVals);

  };

}


#endif
