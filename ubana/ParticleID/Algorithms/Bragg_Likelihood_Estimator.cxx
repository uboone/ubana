// Implementation file for algorithm that calculates likelihood from comparison of data
// dEdx vs residual range to predicted Bragg peak
//
// Takes in dE/dx and residual range vectors, and a particle species, and returns
// the likelihood, L, for the data to have been produced for that
// given particle species
//
// Likelihood is calculated by comparing the measured dE/dx at a given residual range
// to the expected dE/dx at that residual range. Assumes that dE/dx is Landau distributed
// with a mean given by the theoreticl prediction in the Theory_dEdx_resrange class
// (calculated by Bruce Baller), and a width that is user-configurable. The default widths
// contained in the fcl files have been measured from v06_26_01_13 data and simulation.
//
// Can be run for the following particle species: muon, pion, Kaon, proton
//
// Tracks are fit in both directions (forwards and backwards) to account for reconstruction
// getting the track direction wrong, and the highest likelihood is returned

#ifndef BRAGG_LIKELIHOOD_CXX
#define BRAGG_LIKELIHOOD_CXX

#include "Bragg_Likelihood_Estimator.h"

namespace particleid{

  void Bragg_Likelihood_Estimator::configure(fhicl::ParameterSet const &p){
    nHitsToDrop    = p.get<int>("NHitsToDrop", 1);
    endPointFloatShort    = p.get<double>("EndPointFloatShort", -1.0);
    endPointFloatLong     = p.get<double>("EndPointFloatLong" , 1.0);
    endPointFloatStepSize = p.get<double>("EndPointFloatStepSize", 0.05);

    checkRange = p.get<bool>("CheckRange", true);

    LikelihoodMapsFileName = p.get<std::string>("LikelihoodMapsFile");
    LikelihoodMapsFile = new TFile(LikelihoodMapsFileName.c_str(),"read");
    // first index mu=0, pi=1, k=2, p=3, mip=4
    // second index plane number
    h_lmap[0][0] = (TH2F*)LikelihoodMapsFile->Get("h_mu_bragglikelihoodmap_plane0");
    h_lmap[1][0] = (TH2F*)LikelihoodMapsFile->Get("h_pi_bragglikelihoodmap_plane0");
    h_lmap[2][0] = (TH2F*)LikelihoodMapsFile->Get("h_k_bragglikelihoodmap_plane0");
    h_lmap[3][0] = (TH2F*)LikelihoodMapsFile->Get("h_p_bragglikelihoodmap_plane0");
    h_lmap[4][0] = (TH2F*)LikelihoodMapsFile->Get("h_mip_bragglikelihoodmap_plane0");
    h_lmap[0][1] = (TH2F*)LikelihoodMapsFile->Get("h_mu_bragglikelihoodmap_plane1");
    h_lmap[1][1] = (TH2F*)LikelihoodMapsFile->Get("h_pi_bragglikelihoodmap_plane1");
    h_lmap[2][1] = (TH2F*)LikelihoodMapsFile->Get("h_k_bragglikelihoodmap_plane1");
    h_lmap[3][1] = (TH2F*)LikelihoodMapsFile->Get("h_p_bragglikelihoodmap_plane1");
    h_lmap[4][1] = (TH2F*)LikelihoodMapsFile->Get("h_mip_bragglikelihoodmap_plane1");
    h_lmap[0][2] = (TH2F*)LikelihoodMapsFile->Get("h_mu_bragglikelihoodmap_plane2");
    h_lmap[1][2] = (TH2F*)LikelihoodMapsFile->Get("h_pi_bragglikelihoodmap_plane2");
    h_lmap[2][2] = (TH2F*)LikelihoodMapsFile->Get("h_k_bragglikelihoodmap_plane2");
    h_lmap[3][2] = (TH2F*)LikelihoodMapsFile->Get("h_p_bragglikelihoodmap_plane2");
    h_lmap[4][2] = (TH2F*)LikelihoodMapsFile->Get("h_mip_bragglikelihoodmap_plane2");
  }

  void Bragg_Likelihood_Estimator::printConfiguration(){

    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] PRINTING CONFIGURATION: " << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Number of Hits to Drop : " << nHitsToDrop << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point float long   : " << endPointFloatLong  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point float short  : " << endPointFloatShort  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point step size    : " << endPointFloatStepSize  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Reading likelihood maps from    : " << LikelihoodMapsFileName  << std::endl;

  }

  // Here is a dummy method for the case where you don't want to save the shift. It just calls the second method below but with a dummy variable for "shift"
  double Bragg_Likelihood_Estimator::getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum)

  {
    double dummy;
    return Bragg_Likelihood_Estimator::getLikelihood( dEdx, resRange, particlehypothesis, forward, planenum, dummy);
  }


  double Bragg_Likelihood_Estimator::getLikelihood(std::vector<float> dEdx, std::vector<float> resRange, int particlehypothesis, bool forward, int planenum, double &shift)

  {

    /**
     * Get likelihood map histogram for given particle hypothesis
     * This is only available for muon, proton, kaon, pion, and MIP
     * Return an error if user tries to request a different particle type
     */


    int absph = TMath::Abs(particlehypothesis);
    int i_particle; // mu=0, pi=1, k=2, p=3, mip=4

    switch(absph){
      case 13: // muon
        i_particle = 0;
        break;
      case 211: // pion
        i_particle = 1;
        break;
      case 321: // kaon
        i_particle = 2;
        break;
      case 2212: // proton
        i_particle = 3;
        break;
      case 0: // special case: fit to MIP region of muon prediction with no Bragg peak
        i_particle = 4;
        break;
      default:
        std::cout << "[ParticleID::Bragg_Likelihood_Estimator] ERROR: cannot calculate theoretical prediction for given particle hypothesis: " << particlehypothesis << ". Theoretical predictions are only available for charged muons (+/-13), pions (+/-211), kaons (+/-321), protons (2212), and non-Bragg MIP region (0)" << std::endl;
        std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Exiting." << std::endl;
        throw;
    } // switch


    /**
     * Now loop through hits (entries in dEdx and resRange vectors), and get likelihood
     * rr_shift allows us to shift the residual range so that we
     * can account for end point resolution
     */

    double likelihoodNdf = 0;
    int n_hits_used_total = 0;
    double rr_shift_final = 0.;

    for (double rr_shift = endPointFloatShort; rr_shift < endPointFloatLong; rr_shift = rr_shift+endPointFloatStepSize){

      double likelihood = 0.;
      int n_hits_used = 0;

      for (size_t i_hit=0; i_hit < resRange.size(); i_hit++){

        size_t rr_index;
        // Fit tracks "forward" (i.e. in the direction they already have)
        if (forward){
          rr_index = i_hit;
          if ((int)rr_index >= (int)resRange.size() - nHitsToDrop){
            continue;
          }
        }
        // Fit tracks "backward"
        else{
          rr_index = (resRange.size()-1)-i_hit;
          if ((int)i_hit < nHitsToDrop){
            continue;
          }
        }

        double resrg_i = resRange.at(rr_index)+rr_shift;
        double dEdx_i  = dEdx.at(i_hit);

        /**
         * Theory values are only defined up to 30 cm residual Range so we
         * can't compare beyond that
         Note: check range only if particle hypothesis is not 0 (i.e. don't bother checking range for MIP hypothesis because it doesn't look at the Bragg peak)
         */
        if (checkRange && (resrg_i > 30.0 || resrg_i < 0.0)) continue;

        // Evaluate likelihood
        int bin = h_lmap[i_particle][planenum]->FindBin(resrg_i,dEdx_i);
        double likelihood_i = h_lmap[i_particle][planenum]->GetBinContent(bin);
        if (likelihood_i == 0){
          continue;
        }
        else{
          n_hits_used++;
        }

        likelihood += likelihood_i;

      } // Loop over i_hit

      if (likelihood/n_hits_used > likelihoodNdf){
        likelihoodNdf = likelihood/n_hits_used;
        n_hits_used_total = n_hits_used;
        rr_shift_final = rr_shift;
      }

    } // residual range shift

    // Set final residual range shift
    if (shift){
      shift = rr_shift_final;
    }

    if (n_hits_used_total == 0)
      return -9999999;
    else return likelihoodNdf;

  } // getLikelihood


} // particleid namespace

#endif
