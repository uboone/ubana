/**
 * This uses a Gaussian Kernel implementation of the kernel density estimator.
 * KDEs are used to predict the underlying distribution of sparse data sets by
 * replacing each data point with a gaussian and summing the individual
 * Gaussians to produce some distribution.
 *
 * This implementation is a slightly more advanced implementation of this, in
 * particular it scales the Gaussian width in an anti-correlated way to the
 * local density of data points.
 */

#include "KernelDensityEstimator.h"

namespace kde{

  Double_t KernelDensityEstimator::epoKernelFunction(Double_t* x, Double_t* p){

    x[0] = ((x[0] - p[1])/p[2]) ;
    if (std::abs(x[0]) < 1)
      return p[0] * ((3./4.) * (1 - std::pow(x[0],2)));
    else
      return 0;

  }

  Double_t KernelDensityEstimator::gausKernelFunction(Double_t* x, Double_t* p){

    //x = ((x - testPoint)/bandwith);
    return p[0] * (1./(std::sqrt(2*3.1415*std::pow(p[2],2))))*std::exp(-std::pow((x[0]-p[1]),2)/(2*std::pow(p[2],2)));

  }

  float KernelDensityEstimator::getLocalDensity(std::vector<TF1*> kernels, int testpoint, float pilotBandwith){

    TF1* testKernel = kernels.at(testpoint);

    float W = 0;
    for (size_t i = 0; i < kernels.size(); i++){
      W += testKernel->Eval(kernels.at(i)->GetParameter(1));
    }

    float P = W/(kernels.size()*pilotBandwith);

    return P;


  }

  TF1* KernelDensityEstimator::getKernel(float kernelNormalisation, float kernelMean, float kernelBandwith, std::string kernelType){

    TF1* kernel;

    if (kernelType == "epo"){
      kernel = new TF1("kernel", kde::KernelDensityEstimator::epoKernelFunction, -100, 100, 3);
      kernel->SetParameters((double)kernelNormalisation,(double)kernelMean,(double)kernelBandwith);
    }
    else {
      //kernel = new TF1("kernel", kde::KernelDensityEstimator::gausKernelFunction, -100, 100, 3);
      //kernel = new TF1("kernel", "[0] * (1./(std::sqrt(2*3.1415*std::pow([2],2))))*std::exp(-std::pow((x-[1]),2)/(2*std::pow([2],2)))", -100, 100);
      kernel = new TF1("kernel", "gaus", -100, 100);
      kernel->SetParameters((double)kernelNormalisation,(double)kernelMean,(double)kernelBandwith);
    }

    return kernel;


  }

  float KernelDensityEstimator::getKernelDensityMpv(std::vector<float> pidaVals){

    // to be fhiclised
    std::string kernelType = "gaus";
    bool doAdaptiveKde = true;
    float pilotBandwith = 1.0; // for gaussian kernel

    kde::KernelDensityEstimator kde;

    TF1* kernel;
    std::vector<TF1*> kernels;
    float kernelNormalisation = 1.0;
    float kernelMean = -1;

    // KDE is built up of many gaussians.
    // Easiest way to do this is to append the formulas
    // to a string and then use that to construct
    // final KDE
    std::string kdeFormula = "0.";

    //std::cout << "[KDE]  Using KDE method to evaluate PIDA. " << std::endl;
    //std::cout << "[KDE]  Using Kernel Type: " << kernelType << std::endl;

    for (size_t i = 0; i < pidaVals.size(); i++){

      kernelMean = pidaVals.at(i);

      kernel = kde.getKernel(kernelNormalisation, kernelMean, pilotBandwith, kernelType);

      kernels.push_back(kernel);

      if (!doAdaptiveKde){
        //std::cout << "[KDE]  Using static baseline." << std::endl;
        kdeFormula.append("+"+kernel->GetFormula()->GetExpFormula("P"));
      }
    }

    if (doAdaptiveKde){

      //std::cout << "[KDE]  Using adaptive baseline." << std::endl;

      float adaptiveBandwith = -1;

      for (size_t i = 0; i < kernels.size(); i++){

        kernel = kernels.at(i);

        if (kernelType == "epo") adaptiveBandwith = pilotBandwith*0.025/(std::sqrt(kernels.size())*std::log((kde.getLocalDensity(kernels, i, pilotBandwith)+1)));
        else adaptiveBandwith = pilotBandwith*2.5/(std::sqrt(kernels.size())*std::log((kde.getLocalDensity(kernels, i, pilotBandwith)+1)));

        float kernelMean = kernel->GetParameter(1);

        kernel = (TF1*)getKernel(kernelNormalisation, kernelMean, adaptiveBandwith, kernelType);

        float kernelIntegral = kernel->Integral(-100, 100);

        if (!std::isnormal(kernelIntegral)){
          return -1;
        }

        kernelNormalisation = 1./kernelIntegral;

        //std::cout << "[KDE]  >> kernelIntegral:      "<< kernelIntegral << std::endl;
        //std::cout << "[KDE]  >> kernelNormalisation: "<< kernelNormalisation << std::endl;
        //std::cout << "[KDE]  >> kernelMean:          "<< kernelMean << std::endl;
        //std::cout << "[KDE]  >> kernelBandwith:      "<< adaptiveBandwith << std::endl;

        // re-do to get correct normalisation
        kernel = (TF1*)getKernel(kernelNormalisation, kernelMean, adaptiveBandwith, kernelType);
        kdeFormula.append("+"+kernel->GetFormula()->GetExpFormula("P"));

      }

    }

    const char* kdeFormulaChar = kdeFormula.c_str();
    TF1* kernelDensityEstimation = new TF1("kernelDensityEstimation", kdeFormulaChar, -100, 100);

    float kdeMpv = kernelDensityEstimation->GetMaximumX();

    return kdeMpv;
  }

}
