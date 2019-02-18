/*
  This code is intended to generate 2D probability maps for the Bragg Likelihood PID algorithm (in which we treat charge deposition as modelled by a Landau-Gaussian function that peaks at the value of dE/dx expected for a given distance from a Bragg peak).

  Run with:
    cd Scripts
    root -b
    root [0] .L Make2DLikelihoodMaps.C
    root [1] Make2DLikelihoodMaps()
    root [2] .q
    cp BraggLikelihoodMaps.root

  The code takes about 3 minutes to run.

  The variables defined at the top of this code (gausWidth and landauWidth for each particle species) should be determined from calibrated MC and data. For details, see DocDB .
  *** NOTE: THESE VALUES WILL NEED TO BE RECALCULATED EVERY TIME THERE IS A NEW CALIBRATION OR PRODUCTION ***

  The output file should be copied to the folder uboonedata/ParticleID/ and to stashcache: /pnfs/uboone/persistent/stash/ParticleID_LLH_Maps/
*/

// cpp and ROOT
#include <vector>
#include <iostream>

// local includes
#include "../Algorithms/Theory_dEdx_resrange.cxx"
#include "../Algorithms/LandauGaussian.h"

void Make2DLikelihoodMaps(){

  // Note: vector is {plane0, plane1, plane2}

  // MC with dE/dx smearing applied, and data (MCC 8)
  std::vector<std::vector<double>> gausWidth = {
    {0.33, 0.53, 0.20}, // mu
    {0.33, 0.53, 0.20}, // pi
    {0.33, 0.53, 0.20}, // k
    {0.45, 0.60, 0.31}, // p
    {0.33, 0.53, 0.20} // MIP
  };
  std::vector<std::vector<double>> landauWidth = {
    {0.10, 0.12, 0.09}, // mu
    {0.10, 0.12, 0.09}, // pi
    {0.10, 0.12, 0.09}, // k
    {0.23, 0.19, 0.13}, // p
    {0.10, 0.12, 0.09} // MIP
  };

  particleid::Theory_dEdx_resrange theory;
  TGraph *theorypred[5] = {theory.g_ThdEdxRR_Muon, theory.g_ThdEdxRR_Pion, theory.g_ThdEdxRR_Kaon, theory.g_ThdEdxRR_Proton, theory.g_ThdEdxRR_MuonNoBragg};
  TF1 *langaus = new TF1("langaus", landauGaussian, 0, 100, 4);
  TF1 *landau = new TF1("landau", "TMath::Landau(x, [0], [1], [2])", 0, 100);

  TFile *fout = new TFile("BraggLikelihoodMaps.root","RECREATE");

  // These are the TH2Ds that we will use to store the likelihood maps
  TH2F *h_likelihoodmap[5][3];
  for (int plane=0; plane<3; plane++){
    TString name = Form("h_mu_bragglikelihoodmap_plane%d",plane);
    h_likelihoodmap[0][plane] = new TH2F(name.Data(),"Muon Bragg peak;dE/dx (MeV/cm);Residual range (cm)",107,0,32.1,2000,0,100);

    name = Form("h_pi_bragglikelihoodmap_plane%d",plane);
    h_likelihoodmap[1][plane] = new TH2F(name.Data(),"Pion Bragg peak;dE/dx (MeV/cm);Residual range (cm)",107,0,32.1,2000,0,100);

    name = Form("h_k_bragglikelihoodmap_plane%d",plane);
    h_likelihoodmap[2][plane] = new TH2F(name.Data(),"Kaon Bragg peak;dE/dx (MeV/cm);Residual range (cm)",107,0,32.1,2000,0,100);

    name = Form("h_p_bragglikelihoodmap_plane%d",plane);
    h_likelihoodmap[3][plane] = new TH2F(name.Data(),"Proton Bragg peak;dE/dx (MeV/cm);Residual range (cm)",107,0,32.1,2000,0,100);

    name = Form("h_mip_bragglikelihoodmap_plane%d",plane);
    h_likelihoodmap[4][plane] = new TH2F(name.Data(),"MIP (no Bragg peak);dE/dx (MeV/cm);Residual range (cm)",107,0,32.1,2000,0,100);
  }

  // --- Fill likelihood histograms --- //
  // Loop through particle types
  // In order: mu, pi, k, p, MIP
  for (int ipdg=0; ipdg<5; ipdg++){
    // Loop through planes
    for (int plane=0; plane<3; plane++){
      // Create Landau of correct width, and calculate offset between MPV and mean
      // N.b. we want to do this with the Landau not the Landau-Gaussian
      langaus->SetParameters(landauWidth.at(ipdg).at(plane), 10, 1, gausWidth.at(ipdg).at(plane));
      landau->SetParameters(10, landauWidth.at(ipdg).at(plane), 0);
      double landau_mean = landau->Mean(0, 100);
      double landau_mpv  = landau->GetMaximumX();
      double landau_mean_mpv_offset = landau_mean - landau_mpv;

      // Loop through bins in histograms and fill likelihoods
      for (int binx=1; binx<h_likelihoodmap[0][0]->GetXaxis()->GetNbins()+1; binx++){
        for (int biny=1; biny<h_likelihoodmap[0][0]->GetYaxis()->GetNbins()+1; biny++){

          double dEdx_i = h_likelihoodmap[0][0]->GetYaxis()->GetBinCenter(biny);
          double resrg_i = h_likelihoodmap[0][0]->GetXaxis()->GetBinCenter(binx);

          langaus->SetParameters(landauWidth.at(ipdg).at(plane),theorypred[ipdg]->Eval(resrg_i,0,"S")-landau_mean_mpv_offset,1, gausWidth.at(ipdg).at(plane));

          // Evaluate likelihood
          double likelihood_i = langaus->Eval(dEdx_i);
          h_likelihoodmap[ipdg][plane]->SetBinContent(binx,biny,likelihood_i);

        } // end loop over biny
      } // end loop over binx


      // Finally: save histograms
      fout->cd();
      h_likelihoodmap[ipdg][plane]->Write();
    } // end loop over planes
  } // end loop over particle types (ipdg)

  fout->Close();

}
