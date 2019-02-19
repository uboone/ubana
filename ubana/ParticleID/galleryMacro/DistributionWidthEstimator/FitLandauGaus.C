#include "../../Algorithms/LandauGaussian.h"

void FitLandauGaus(std::string filename, bool fixGaus=false)
{
  TF1 *langaus = new TF1("langaus", landauGaussian, 0, 10, 4);
  langaus->SetParNames("Landau width","Peak value","Normalisation","Gaussian width");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gROOT->ForceStyle();

  TFile *f = new TFile(filename.c_str(),"open");
  TFile *fout = new TFile("output_LanGausFit.root","recreate");

  // All tracks
  TCanvas *c1 = new TCanvas();

  TH1D *hdEdx_all = (TH1D*)f->Get("hdEdx_all_plane2_BNB");
  hdEdx_all->GetXaxis()->SetRangeUser(0,10);
  langaus->SetParameters(0.2,1.7,1000,0.1);
  // Uncomment for data fit
  langaus->FixParameter(0,0.09);
  hdEdx_all->Fit(langaus,"","",1.,10.);
  hdEdx_all->Draw();
  hdEdx_all->Write("hdEdx_all");
  TF1 *fit_hdEdx_all = (TF1*)langaus->Clone("fit_hdEdx_all");
  fit_hdEdx_all->Write();
  c1->Write("canv_hdEdx_all");
  c1->Print("fitted_hdEdx_all.pdf");

  // Get Gaussian width from this fit and fix for other distributions
  double GausWidth = fit_hdEdx_all->GetParameter(3);


  langaus->SetParameters(0.1,1.7,1000,GausWidth);
  if (fixGaus) langaus->FixParameter(3, GausWidth);
  // True muons
  TH1D *hdEdx_Muon = (TH1D*)f->Get("hdEdx_Muon_plane2");
  if (hdEdx_Muon->Integral() > 0){
    hdEdx_Muon->GetXaxis()->SetRangeUser(0,10);
    TCanvas *c2 = new TCanvas();
    hdEdx_Muon->Fit(langaus,"","",1.,10.);
    hdEdx_Muon->Draw();
    hdEdx_Muon->Write("hdEdx_Muon");
    TF1 *fit_hdEdx_Muon = (TF1*)langaus->Clone("fit_hdEdx_Muon");
    fit_hdEdx_Muon->Write();
    c2->Write("canv_hdEdx_Muon");
    c2->Print("fitted_hdEdx_Muon.pdf");
  }


  langaus->SetParameters(0.1,2.2,100,GausWidth);
  if (fixGaus) langaus->FixParameter(3, GausWidth);
  // True protons
  TH1D *hdEdx_Proton = (TH1D*)f->Get("hdEdx_Proton_plane2");
  if (hdEdx_Proton->Integral() > 0){
    hdEdx_Proton->Rebin(2);
    hdEdx_Proton->GetXaxis()->SetRangeUser(0,10);
    TCanvas *c3 = new TCanvas();
    hdEdx_Proton->Fit(langaus,"","",1.5,10.);
    hdEdx_Proton->Draw();
    hdEdx_Proton->Write("hdEdx_Proton");
    TF1 *fit_hdEdx_Proton = (TF1*)langaus->Clone("fit_hdEdx_Proton");
    fit_hdEdx_Proton->Write();
    c3->Write("canv_hdEdx_Proton");
    c3->Print("fitted_hdEdx_Proton.pdf");
  }
}
