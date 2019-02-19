/**
 * example script for how to fit a Landau Gaussian distribution
 */ 

#include "../Algorithms/LandauGaussian.h"

void FitLandauGaus()
{
  TString plane("y");
  TF1 *langaus = new TF1("langaus", landauGaussian, 0, 10, 4);
  langaus->SetParNames("Landau width","Peak value","Normalisation","Gaussian width");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gROOT->ForceStyle();

  // All tracks
  TCanvas *c1 = new TCanvas();
  TTree *tree = (TTree*)_file0->Get("pidvalid/pidTree");

  TH1D* h = new TH1D("h", "", 100, 0, 10);
  tree->Draw("track_dEdx_perhit_"+plane+" >> h", "track_resrange_perhit_"+plane+" > 100 && track_resrange_perhit_"+plane+" < 150");

  langaus->SetParameters(0.2,1.7,1000,0.1);
  // Uncomment for data fit
  langaus->FixParameter(0,0.09);
  h->Fit(langaus,"","",1.0,3.);
  h->Draw();
  
  c1->SaveAs("simsmeared_dEdx_"+plane+"plane.png");
}
