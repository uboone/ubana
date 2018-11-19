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
    double getLikelihood(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward, int planenum);
    double getLikelihood(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward, int planenum, double &shift);

  //private:
    std::vector<double> gausWidth_mu;
    std::vector<double> gausWidth_pi;
    std::vector<double> gausWidth_k;
    std::vector<double> gausWidth_p;
    std::vector<double> gausWidth_mip;
    std::vector<double> landauWidth_mu;
    std::vector<double> landauWidth_pi;
    std::vector<double> landauWidth_k;
    std::vector<double> landauWidth_p;
    std::vector<double> landauWidth_mip;
    double offset_p;
    double offset_mu;
    double offset_pi;
    double offset_k;
    double offset_mip;
    int nHitsToDrop;
    double endPointFloatShort;
    double endPointFloatLong;
    double endPointFloatStepSize;

    bool checkRange;
  };

}

#endif
