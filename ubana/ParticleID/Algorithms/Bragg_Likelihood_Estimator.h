// Implementation file for algorithm that calculates likelihood from comparison of data
// dEdx vs residual range to predicted Bragg peak
//
// Takes in dE/dx and residual range vectors, and a particle species, and returns
// the likelihood for the data to have been produced for that
// given particle species
//
// Likelihood is calculated by comparing the measured dE/dx at a given residual range
// to the expected dE/dx at that residual range. Assumes that dE/dx is Landau-Gaussian distributed
// with a mean given by the theoretical prediction in the Theory_dEdx_resrange class
// (calculated by Bruce Baller), and with Landau/Gaussian widths which are user configurable.
//
// Can be run for the following particle species: muon, pion, kaon, proton
//
// Tracks are fit in both directions (forwards and backwards) to account for reconstruction
// getting the track direction wrong, and the highest likelihood is returned

#ifndef BRAGG_LIKELIHOOD_H
#define BRAGG_LIKELIHOOD_H

// cpp
#include <vector>
#include <iostream>

// ROOT
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

// art
#include "fhiclcpp/ParameterSet.h"

// local
#include "Theory_dEdx_resrange.h"
#include "LandauGaussian.h"

namespace particleid{

  class Bragg_Likelihood_Estimator{

  public:
    void configure(fhicl::ParameterSet const &p);
    void printConfiguration();
    float getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum);
    float getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum, float &shift);

  //private:
    std::vector<float> gausWidth_mu;
    std::vector<float> gausWidth_pi;
    std::vector<float> gausWidth_k;
    std::vector<float> gausWidth_p;
    std::vector<float> gausWidth_mip;
    std::vector<float> landauWidth_mu;
    std::vector<float> landauWidth_pi;
    std::vector<float> landauWidth_k;
    std::vector<float> landauWidth_p;
    std::vector<float> landauWidth_mip;
    float offset_p;
    float offset_mu;
    float offset_pi;
    float offset_k;
    float offset_mip;
    int nHitsToDrop;
    float endPointFloatShort;
    float endPointFloatLong;
    float endPointFloatStepSize;

    bool checkRange;
  };

}

#endif
