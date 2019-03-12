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
#include <string>

// ROOT
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"

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

    double getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum);
    double getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum, double &shift);


  //private:
    int nHitsToDrop;
    double endPointFloatShort;
    double endPointFloatLong;
    double endPointFloatStepSize;

    bool checkRange;

    std::string LikelihoodMapsFileName;
    TFile *LikelihoodMapsFile;
    TH2F *h_lmap[5][3];
    // first index mu=0, pi=1, k=2, p=3, mip=4
    // second index plane number
  };

}

#endif
