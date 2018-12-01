/**
 * ROOT script to measure the landau component of the dE/dx width with a fixed 
 * Gaussian width, as a function of residual range
 */

#include "../Algorithms/LandauGaussian.h"

void FitLanGausPerUnitLength(){

  /**
   * setup
   */
  const int n_dedx_bins = 50;
  const double max_dedx = 20;
  const int n_resrg_bins = 25;

  double landau_widths[n_resrg_bins];
  double chisq_ndfs[n_resrg_bins];
  double resrgs[n_resrg_bins];

  TFile *file = new TFile("/uboone/data/users/alister1/particleID/particleId.root", "read");
  TFile *outfile = new TFile("outfile.root", "RECREATE");
  TTree *tree = (TTree*)file->Get("pidvalid/pidTree");

  /**
   * Initialise lan-gaus func
   */
  TF1 *langaus = new TF1("langaus", landauGaussian, 0, 10, 4);
  langaus->SetParNames("Landau width","Peak value","Normalisation","Gaussian width");

  /**
   * Get 2D dE/dx - residual range histogram
   */
  TH2D* h_protons = new TH2D("h_protons", ";Residual Range (cm); dE/dx (MeV/cm)", n_resrg_bins, 0, 30, n_dedx_bins, 0, max_dedx);

  tree->Draw("track_dEdx_perhit[][2]:track_resrange_perhit[][2] >> h_protons", 
      "true_PDG == 2212 && true_start_x > 10 && true_start_x < 246 && true_end_x > 10 && true_end_x < 246 && true_start_y > -106.5 && true_start_y < 106.5 && true_end_y > -106.5 && true_end_y < 106.5 && true_start_z > 10 && true_start_z < 1030 && true_end_z > 10 && true_end_z < 1030 && true_end_momentum == 0");

  outfile->cd();
  h_protons->Draw("colz");
  h_protons->Write();

  /** 
   * For each x bin, fit the Landau Gaussian and extract the landau component
   * width
   */
  for (int i = 0; i < h_protons->GetNbinsX(); i++){

    TH1D* h = new TH1D(Form("h_%fcm", i*(double)30./n_resrg_bins), Form(";dE/dx at %d cm residual range;",i*30/n_resrg_bins), n_dedx_bins, 0, max_dedx);
    
    for(int j = 0; j < h_protons->GetNbinsY(); j++){

      h->SetBinContent(j, h_protons->GetBinContent(i,j));

    }
 
    h->GetXaxis()->SetRangeUser(2.5,max_dedx);
    langaus->SetParameters(0.1,h->GetBinLowEdge(h->GetMaximumBin()),h->GetMaximum(),0.1);
    langaus->FixParameter(3, 0.25);
    h->Fit(langaus,"","",2.5,max_dedx);

    if (h->Integral() != 0){

      landau_widths[i] = h->GetFunction("langaus")->GetParameter(0);
      chisq_ndfs[i] = h->GetFunction("langaus")->GetChisquare()/h->GetFunction("langaus")->GetNDF();

    }
    else{

      landau_widths[i] = 0;
      chisq_ndfs[i] = 0;

    }

    resrgs[i] = i*30.0/n_resrg_bins;

    h->Write();

  }

  /**
   * Draw width and chisq as a function of residual range to graphs, 
   * and save these to file
   */
  
  TCanvas* c2 = new TCanvas("c2", "", 500, 500);
  c2->cd();

  TGraph *gr_landau_widths = new TGraph(n_resrg_bins, resrgs, landau_widths);
  TGraph *gr_chisq_ndf = new TGraph(n_resrg_bins, resrgs, chisq_ndfs);

  gr_landau_widths->SetMarkerStyle(20);
  gr_landau_widths->SetLineStyle(0);
  gr_landau_widths->SetLineWidth(0);
  gr_landau_widths->SetMarkerColor(kGreen+1);
  gr_landau_widths->SetTitle(";Residual Range (cm); Landau Width (MeV/cm)");
  gr_chisq_ndf->SetMarkerStyle(20);
  gr_chisq_ndf->SetLineStyle(0);
  gr_chisq_ndf->SetLineWidth(0);
  gr_chisq_ndf->SetMarkerColor(kAzure+1);

  gr_chisq_ndf->Draw("AP");
  gr_landau_widths->Draw("sameP");

  TF1* expon = new TF1("expon", "[0]*TMath::Exp([1]*x)+[2]");
  TF1* xsqr  = new TF1("xsqr", "[0]*std::pow(x,2)+[1]");
  gr_landau_widths->Fit("expon", "", "", 3.0, 30);

  TLegend *leg = new TLegend(0.4, 0.75, 0.79, 0.85);
  leg->AddEntry(gr_landau_widths, "Landau Widths");
  leg->AddEntry(gr_chisq_ndf, "Chisq/NDF of fit");
  leg->Draw("same");

  c2->Write();

}


