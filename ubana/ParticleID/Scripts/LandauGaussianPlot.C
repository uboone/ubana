#include "LandauGaussianPlot.h"

/**
 * This uses the landau-gaussian in LandauGaussian.h to produce a 
 * "probability map" for protons, muons, and muons with no Bragg Peak. It
 * also produces a profile TGraphAsymmErrors where the point represents
 * the MPV of the distribution and the upper and lower limits represent
 * the 84th and 16th quantiles.
 * This is mostly for demonstration purposes.
 */

void LandauGaussianPlot(){

  const int nxbins = 500;
  const double xbinlow = 0;
  const double xbinhigh = 30;
  const int nybins = 500;
  const double ybinlow = 0;
  const double ybinhigh = 20;

  // default data:
  /**
    double muonlandauwidth = 0.09;
    double muongaussianwidth = 0.2;
    double protonlandauwidth = 0.19;
    double protongaussianwidth = 0.2;
    */
  // default MC:
  double muonlandauwidth = 0.09;
  double muongaussianwidth = 0.09;
  double pionlandauwidth = 0.09;
  double piongaussianwidth = 0.09;
  double kaonlandauwidth = 0.09;
  double kaongaussianwidth = 0.09;
  double protonlandauwidth = 0.19;
  double protongaussianwidth = 0.09;

  double muonMeanMPVCorrection = 0.60302250;
  double pionMeanMPVCorrection = 0.60302250;
  double kaonMeanMPVCorrection = 0.60302250;
  double protonMeanMPVCorrection = 1.1297;

  Theory_dEdx_resrange();

  TF1 *langaus = new TF1("langaus", landauGaussian, 0, 100, 4);
  langaus->SetNpx(10000);

  /** muon */

  /** variables for creating TGraphAsymmErrors*/
  double muon_resRg[nxbins];
  double muon_resRg_up[nxbins];
  double muon_resRg_down[nxbins];
  double muon_mpv[nxbins];
  double muon_errup[nxbins];
  double muon_errdown[nxbins];

  TH2D *h_muon = new TH2D("h_muon", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  double xbinwidth = h_muon->GetXaxis()->GetBinLowEdge(nxbins) / h_muon->GetNbinsX();
  double halfbinwidthx = xbinwidth/2.0;
  double ybinwidth = h_muon->GetYaxis()->GetBinLowEdge(nybins) / h_muon->GetNbinsY();
  double halfbinwidthy = ybinwidth/2.0; 

  for (int i = 0; i < nxbins; i++){

    langaus->SetParameters(muonlandauwidth, g_ThdEdxRR_Muon->Eval((i)*xbinhigh/(double)nxbins, 0, "S")-muonMeanMPVCorrection, 1.0, muongaussianwidth);

    for (int j = 0; j < nybins; j++){

      h_muon->SetBinContent(i,j,langaus->Eval(h_muon->GetYaxis()->GetBinCenter(j)));

    }

    /** 
     * project out chosen x bin in order to use
     * GetQuantiles method which is only available in TH1s, then find
     * the 16th and 84th quantiles
     */

    TH1D* h_muon_proj_y = (TH1D*)h_muon->ProjectionY("h_muon_proj_y", i, i+1);
    Double_t yq[2];
    Double_t xq[2] = {0.16, 0.84};
    h_muon_proj_y->GetQuantiles(1,&yq[0], &xq[0]);
    h_muon_proj_y->GetQuantiles(1,&yq[1], &xq[1]);

    muon_resRg[i] = i*xbinhigh/double(nxbins) - halfbinwidthx;
    muon_resRg_up[i] = 0;
    muon_resRg_down[i] = 0;
    muon_mpv[i] = h_muon_proj_y->GetBinCenter(h_muon_proj_y->GetMaximumBin());
    muon_errup[i] = yq[1] - muon_mpv[i];
    muon_errdown[i] = muon_mpv[i] - yq[0];

  }

  TGraphAsymmErrors* gr_muon_profile_x = new TGraphAsymmErrors(nxbins, muon_resRg, muon_mpv, muon_resRg_down, muon_resRg_up, muon_errdown, muon_errup);
  gr_muon_profile_x->SetName("gr_muon_profile_x");

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->cd();
  h_muon->SetContour(100);
  h_muon->Draw("col2");

  /** muonnobragg */

  /** variables for creating TGraphAsymmErrors*/
  double muonnobragg_resRg[nxbins];
  double muonnobragg_resRg_up[nxbins];
  double muonnobragg_resRg_down[nxbins];
  double muonnobragg_mpv[nxbins];
  double muonnobragg_errup[nxbins];
  double muonnobragg_errdown[nxbins];

  TH2D *h_muonnobragg = new TH2D("h_muonnobragg", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  for (int i = 0; i < nxbins; i++){

    langaus->SetParameters(muonlandauwidth, g_ThdEdxRR_MuonNoBragg->Eval((i)*xbinhigh/(double)nxbins, 0, "S"), 1.0, muongaussianwidth);

    for (int j = 0; j < nybins; j++){

      h_muonnobragg->SetBinContent(i,j,langaus->Eval(h_muonnobragg->GetYaxis()->GetBinCenter(j)));

    }

    /** 
     * project out chosen x bin in order to use
     * GetQuantiles method which is only available in TH1s, then find
     * the 16th and 84th quantiles
     */

    TH1D* h_muonnobragg_proj_y = (TH1D*)h_muonnobragg->ProjectionY("h_muonnobragg_proj_y", i, i+1);
    Double_t yq[2];
    Double_t xq[2] = {0.16, 0.84};
    h_muonnobragg_proj_y->GetQuantiles(1,&yq[0], &xq[0]);
    h_muonnobragg_proj_y->GetQuantiles(1,&yq[1], &xq[1]);

    muonnobragg_resRg[i] = i*xbinhigh/double(nxbins) - halfbinwidthx;
    muonnobragg_resRg_up[i] = 0;
    muonnobragg_resRg_down[i] = 0;
    muonnobragg_mpv[i] = h_muonnobragg_proj_y->GetBinCenter(h_muonnobragg_proj_y->GetMaximumBin());
    muonnobragg_errdown[i] = muonnobragg_mpv[i] - yq[0];
    muonnobragg_errup[i] = yq[1] - muonnobragg_mpv[i];

  }

  TGraphAsymmErrors* gr_muonnobragg_profile_x = new TGraphAsymmErrors(nxbins, muonnobragg_resRg, muonnobragg_mpv, muonnobragg_resRg_down, muonnobragg_resRg_up, muonnobragg_errdown, muonnobragg_errup);
  gr_muonnobragg_profile_x->SetName("gr_muonnobragg_profile_x");

  TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
  c3->cd();
  h_muonnobragg->SetContour(100);
  h_muonnobragg->Draw("col2");

  /** proton */

  /** variables for creating TGraphAsymmErrors*/
  double proton_resRg[nxbins];
  double proton_resRg_up[nxbins];
  double proton_resRg_down[nxbins];
  double proton_mpv[nxbins];
  double proton_errup[nxbins];
  double proton_errdown[nxbins];


  TH2D *h_proton = new TH2D("h_proton", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  for (int i = 0; i < nxbins; i++){

    langaus->SetParameters(protonlandauwidth, (g_ThdEdxRR_Proton->Eval((i)*xbinhigh/(double)nxbins, 0, "S"))-protonMeanMPVCorrection, 1.0, protongaussianwidth);

    for (int j = 0; j < nybins; j++){

      h_proton->SetBinContent(i,j,langaus->Eval(h_proton->GetYaxis()->GetBinCenter(j)));

    }

    /** 
     * project out chosen x bin in order to use
     * GetQuantiles method which is only available in TH1s, then find
     * the 16th and 84th quantiles
     */

    TH1D* h_proton_proj_y = (TH1D*)h_proton->ProjectionY("h_proton_proj_y", i, i+1);
    Double_t yq[2];
    Double_t xq[2] = {0.16, 0.84};
    h_proton_proj_y->GetQuantiles(1,&yq[0], &xq[0]);
    h_proton_proj_y->GetQuantiles(1,&yq[1], &xq[1]);

    proton_resRg[i] = i*xbinhigh/double(nxbins) - halfbinwidthx;
    proton_resRg_up[i] = 0;
    proton_resRg_down[i] = 0;
    proton_mpv[i] = h_proton_proj_y->GetBinCenter(h_proton_proj_y->GetMaximumBin());
    proton_errdown[i] = proton_mpv[i] - yq[0];
    proton_errup[i] = yq[1] - proton_mpv[i];

  }

  TGraphAsymmErrors* gr_proton_profile_x = new TGraphAsymmErrors(nxbins, proton_resRg, proton_mpv, proton_resRg_down, proton_resRg_up, proton_errdown, proton_errup);
  gr_proton_profile_x->SetName("gr_proton_profile_x");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();
  h_proton->SetContour(100);
  h_proton->Draw("col2");

  /** pion */

  /** variables for creating TGraphAsymmErrors*/
  double pion_resRg[nxbins];
  double pion_resRg_up[nxbins];
  double pion_resRg_down[nxbins];
  double pion_mpv[nxbins];
  double pion_errup[nxbins];
  double pion_errdown[nxbins];


  TH2D *h_pion = new TH2D("h_pion", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  for (int i = 0; i < nxbins; i++){

    langaus->SetParameters(pionlandauwidth, (g_ThdEdxRR_Pion->Eval((i)*xbinhigh/(double)nxbins, 0, "S"))-pionMeanMPVCorrection, 1.0, piongaussianwidth);

    for (int j = 0; j < nybins; j++){

      h_pion->SetBinContent(i,j,langaus->Eval(h_pion->GetYaxis()->GetBinCenter(j)));

    }

    /** 
     * project out chosen x bin in order to use
     * GetQuantiles method which is only available in TH1s, then find
     * the 16th and 84th quantiles
     */

    TH1D* h_pion_proj_y = (TH1D*)h_pion->ProjectionY("h_pion_proj_y", i, i+1);
    Double_t yq[2];
    Double_t xq[2] = {0.16, 0.84};
    h_pion_proj_y->GetQuantiles(1,&yq[0], &xq[0]);
    h_pion_proj_y->GetQuantiles(1,&yq[1], &xq[1]);

    pion_resRg[i] = i*xbinhigh/double(nxbins) - halfbinwidthx;
    pion_resRg_up[i] = 0;
    pion_resRg_down[i] = 0;
    pion_mpv[i] = h_pion_proj_y->GetBinCenter(h_pion_proj_y->GetMaximumBin());
    pion_errdown[i] = pion_mpv[i] - yq[0];
    pion_errup[i] = yq[1] - pion_mpv[i];

  }

  TGraphAsymmErrors* gr_pion_profile_x = new TGraphAsymmErrors(nxbins, pion_resRg, pion_mpv, pion_resRg_down, pion_resRg_up, pion_errdown, pion_errup);
  gr_pion_profile_x->SetName("gr_pion_profile_x");

  TCanvas *c4 = new TCanvas("c4", "c4", 500, 500);
  c4->cd();
  h_pion->SetContour(100);
  h_pion->Draw("col2");

  /** kaon */

  /** variables for creating TGraphAsymmErrors*/
  double kaon_resRg[nxbins];
  double kaon_resRg_up[nxbins];
  double kaon_resRg_down[nxbins];
  double kaon_mpv[nxbins];
  double kaon_errup[nxbins];
  double kaon_errdown[nxbins];


  TH2D *h_kaon = new TH2D("h_kaon", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  for (int i = 0; i < nxbins; i++){

    langaus->SetParameters(kaonlandauwidth, (g_ThdEdxRR_Kaon->Eval((i)*xbinhigh/(double)nxbins, 0, "S"))-kaonMeanMPVCorrection, 1.0, kaongaussianwidth);

    for (int j = 0; j < nybins; j++){

      h_kaon->SetBinContent(i,j,langaus->Eval(h_kaon->GetYaxis()->GetBinCenter(j)));

    }

    /** 
     * project out chosen x bin in order to use
     * GetQuantiles method which is only available in TH1s, then find
     * the 16th and 84th quantiles
     */

    TH1D* h_kaon_proj_y = (TH1D*)h_kaon->ProjectionY("h_kaon_proj_y", i, i+1);
    Double_t yq[2];
    Double_t xq[2] = {0.16, 0.84};
    h_kaon_proj_y->GetQuantiles(1,&yq[0], &xq[0]);
    h_kaon_proj_y->GetQuantiles(1,&yq[1], &xq[1]);

    kaon_resRg[i] = i*xbinhigh/double(nxbins) - halfbinwidthx;
    kaon_resRg_up[i] = 0;
    kaon_resRg_down[i] = 0;
    kaon_mpv[i] = h_kaon_proj_y->GetBinCenter(h_kaon_proj_y->GetMaximumBin());
    kaon_errdown[i] = kaon_mpv[i] - yq[0];
    kaon_errup[i] = yq[1] - kaon_mpv[i];

  }

  TGraphAsymmErrors* gr_kaon_profile_x = new TGraphAsymmErrors(nxbins, kaon_resRg, kaon_mpv, kaon_resRg_down, kaon_resRg_up, kaon_errdown, kaon_errup);
  gr_kaon_profile_x->SetName("gr_kaon_profile_x");

  TCanvas *c5 = new TCanvas("c5", "c5", 500, 500);
  c5->cd();
  h_kaon->SetContour(100);
  h_kaon->Draw("col2");


  /** combined */

  TH2D *h_combined = new TH2D("h_combined", ";Residual Range (cm); dE/dx (MeV/cm)", nxbins, xbinlow, xbinhigh, nybins, ybinlow, ybinhigh);

  for (int i = 0; i < nxbins*nybins; i++){

    h_combined->SetBinContent(i, std::max(std::max(h_proton->GetBinContent(i), h_muon->GetBinContent(i)), h_muonnobragg->GetBinContent(i)));

  }

  // Write histograms and graphs to file
  TFile* f = new TFile("likelihoodMaps.root", "RECREATE");
  f->cd();
  g_ThdEdxRR_Proton->SetName("gr_proton");
  g_ThdEdxRR_Proton->Write();
  g_ThdEdxRR_Muon->SetName("gr_muon");
  g_ThdEdxRR_Muon->Write();
  g_ThdEdxRR_MuonNoBragg->SetName("gr_muonnobragg");
  g_ThdEdxRR_MuonNoBragg->Write();
  g_ThdEdxRR_Pion->SetName("gr_pion");
  g_ThdEdxRR_Pion->Write();
  g_ThdEdxRR_Kaon->SetName("gr_kaon");
  g_ThdEdxRR_Kaon->Write();
  h_muon->Write();
  h_muonnobragg->Write();
  h_proton->Write();
  h_pion->Write();
  h_kaon->Write();
  h_combined->Write();
  gr_muon_profile_x->Write();
  gr_muonnobragg_profile_x->Write();
  gr_proton_profile_x->Write();


}

void Theory_dEdx_resrange(){
  // Constructor: set points in TGraphs and initialise TSpline3s

  // Muon, no Bragg
  g_ThdEdxRR_MuonNoBragg = new TGraph(107);
  g_ThdEdxRR_MuonNoBragg->SetPoint(0,   31.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(1,   31.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(2,   31.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(3,   31.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(4,   30.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(5,   30.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(6,   30.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(7,   29.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(8,   29.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(9,   29.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(10,  28.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(11,  28.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(12,  28.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(13,  28.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(14,  27.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(15,  27.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(16,  27.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(17,  26.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(18,  26.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(19,  26.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(20,  25.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(21,  25.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(22,  25.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(23,  25.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(24,  24.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(25,  24.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(26,  24.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(27,  23.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(28,  23.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(29,  23.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(30,  22.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(31,  22.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(32,  22.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(33,  22.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(34,  21.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(35,  21.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(36,  21.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(37,  20.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(38,  20.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(39,  20.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(40,  19.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(41,  19.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(42,  19.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(43,  19.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(44,  18.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(45,  18.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(46,  18.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(47,  17.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(48,  17.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(49,  17.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(50,  16.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(51,  16.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(52,  16.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(53,  16.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(54,  15.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(55,  15.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(56,  15.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(57,  14.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(58,  14.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(59,  14.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(60,  13.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(61,  13.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(62,  13.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(63,  13.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(64,  12.75, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(65,  12.45, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(66,  12.15, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(67,  11.85, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(68,  11.55, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(69,  11.25, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(70,  10.95, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(71,  10.65, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(72,  10.35, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(73,  10.05, 1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(74,  9.75,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(75,  9.45,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(76,  9.15,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(77,  8.85,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(78,  8.55,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(79,  8.25,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(80,  7.95,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(81,  7.65,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(82,  7.35,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(83,  7.05,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(84,  6.75,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(85,  6.45,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(86,  6.15,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(87,  5.85,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(88,  5.55,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(89,  5.25,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(90,  4.95,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(91,  4.65,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(92,  4.35,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(93,  4.05,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(94,  3.75,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(95,  3.45,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(96,  3.15,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(97,  2.85,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(98,  2.55,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(99,  2.25,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(100, 1.95,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(101, 1.65,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(102, 1.35,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(103, 1.05,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(104, 0.75,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(105, 0.45,  1.6);
  g_ThdEdxRR_MuonNoBragg->SetPoint(106, 0.15,  1.6);


  // Proton
  g_ThdEdxRR_Proton = new TGraph(107);
  g_ThdEdxRR_Proton->SetPoint(0,   31.95, 4.14);
  g_ThdEdxRR_Proton->SetPoint(1,   31.65, 4.16);
  g_ThdEdxRR_Proton->SetPoint(2,   31.35, 4.17);
  g_ThdEdxRR_Proton->SetPoint(3,   31.05, 4.18);
  g_ThdEdxRR_Proton->SetPoint(4,   30.75, 4.20);
  g_ThdEdxRR_Proton->SetPoint(5,   30.45, 4.21);
  g_ThdEdxRR_Proton->SetPoint(6,   30.15, 4.23);
  g_ThdEdxRR_Proton->SetPoint(7,   29.85, 4.25);
  g_ThdEdxRR_Proton->SetPoint(8,   29.55, 4.26);
  g_ThdEdxRR_Proton->SetPoint(9,   29.25, 4.28);
  g_ThdEdxRR_Proton->SetPoint(10,  28.95, 4.29);
  g_ThdEdxRR_Proton->SetPoint(11,  28.65, 4.31);
  g_ThdEdxRR_Proton->SetPoint(12,  28.35, 4.33);
  g_ThdEdxRR_Proton->SetPoint(13,  28.05, 4.34);
  g_ThdEdxRR_Proton->SetPoint(14,  27.75, 4.36);
  g_ThdEdxRR_Proton->SetPoint(15,  27.45, 4.38);
  g_ThdEdxRR_Proton->SetPoint(16,  27.15, 4.40);
  g_ThdEdxRR_Proton->SetPoint(17,  26.85, 4.42);
  g_ThdEdxRR_Proton->SetPoint(18,  26.55, 4.43);
  g_ThdEdxRR_Proton->SetPoint(19,  26.25, 4.45);
  g_ThdEdxRR_Proton->SetPoint(20,  25.95, 4.47);
  g_ThdEdxRR_Proton->SetPoint(21,  25.65, 4.49);
  g_ThdEdxRR_Proton->SetPoint(22,  25.35, 4.51);
  g_ThdEdxRR_Proton->SetPoint(23,  25.05, 4.53);
  g_ThdEdxRR_Proton->SetPoint(24,  24.75, 4.55);
  g_ThdEdxRR_Proton->SetPoint(25,  24.45, 4.57);
  g_ThdEdxRR_Proton->SetPoint(26,  24.15, 4.60);
  g_ThdEdxRR_Proton->SetPoint(27,  23.85, 4.62);
  g_ThdEdxRR_Proton->SetPoint(28,  23.55, 4.64);
  g_ThdEdxRR_Proton->SetPoint(29,  23.25, 4.66);
  g_ThdEdxRR_Proton->SetPoint(30,  22.95, 4.69);
  g_ThdEdxRR_Proton->SetPoint(31,  22.65, 4.71);
  g_ThdEdxRR_Proton->SetPoint(32,  22.35, 4.73);
  g_ThdEdxRR_Proton->SetPoint(33,  22.05, 4.76);
  g_ThdEdxRR_Proton->SetPoint(34,  21.75, 4.78);
  g_ThdEdxRR_Proton->SetPoint(35,  21.45, 4.81);
  g_ThdEdxRR_Proton->SetPoint(36,  21.15, 4.83);
  g_ThdEdxRR_Proton->SetPoint(37,  20.85, 4.86);
  g_ThdEdxRR_Proton->SetPoint(38,  20.55, 4.89);
  g_ThdEdxRR_Proton->SetPoint(39,  20.25, 4.92);
  g_ThdEdxRR_Proton->SetPoint(40,  19.95, 4.94);
  g_ThdEdxRR_Proton->SetPoint(41,  19.65, 4.97);
  g_ThdEdxRR_Proton->SetPoint(42,  19.35, 5.00);
  g_ThdEdxRR_Proton->SetPoint(43,  19.05, 5.03);
  g_ThdEdxRR_Proton->SetPoint(44,  18.75, 5.07);
  g_ThdEdxRR_Proton->SetPoint(45,  18.45, 5.10);
  g_ThdEdxRR_Proton->SetPoint(46,  18.15, 5.13);
  g_ThdEdxRR_Proton->SetPoint(47,  17.85, 5.16);
  g_ThdEdxRR_Proton->SetPoint(48,  17.55, 5.20);
  g_ThdEdxRR_Proton->SetPoint(49,  17.25, 5.23);
  g_ThdEdxRR_Proton->SetPoint(50,  16.95, 5.27);
  g_ThdEdxRR_Proton->SetPoint(51,  16.65, 5.31);
  g_ThdEdxRR_Proton->SetPoint(52,  16.35, 5.35);
  g_ThdEdxRR_Proton->SetPoint(53,  16.05, 5.39);
  g_ThdEdxRR_Proton->SetPoint(54,  15.75, 5.43);
  g_ThdEdxRR_Proton->SetPoint(55,  15.45, 5.47);
  g_ThdEdxRR_Proton->SetPoint(56,  15.15, 5.51);
  g_ThdEdxRR_Proton->SetPoint(57,  14.85, 5.56);
  g_ThdEdxRR_Proton->SetPoint(58,  14.55, 5.60);
  g_ThdEdxRR_Proton->SetPoint(59,  14.25, 5.65);
  g_ThdEdxRR_Proton->SetPoint(60,  13.95, 5.70);
  g_ThdEdxRR_Proton->SetPoint(61,  13.65, 5.75);
  g_ThdEdxRR_Proton->SetPoint(62,  13.35, 5.80);
  g_ThdEdxRR_Proton->SetPoint(63,  13.05, 5.85);
  g_ThdEdxRR_Proton->SetPoint(64,  12.75, 5.91);
  g_ThdEdxRR_Proton->SetPoint(65,  12.45, 5.97);
  g_ThdEdxRR_Proton->SetPoint(66,  12.15, 6.03);
  g_ThdEdxRR_Proton->SetPoint(67,  11.85, 6.09);
  g_ThdEdxRR_Proton->SetPoint(68,  11.55, 6.15);
  g_ThdEdxRR_Proton->SetPoint(69,  11.25, 6.22);
  g_ThdEdxRR_Proton->SetPoint(70,  10.95, 6.29);
  g_ThdEdxRR_Proton->SetPoint(71,  10.65, 6.36);
  g_ThdEdxRR_Proton->SetPoint(72,  10.35, 6.44);
  g_ThdEdxRR_Proton->SetPoint(73,  10.05, 6.52);
  g_ThdEdxRR_Proton->SetPoint(74,  9.75,  6.60);
  g_ThdEdxRR_Proton->SetPoint(75,  9.45,  6.68);
  g_ThdEdxRR_Proton->SetPoint(76,  9.15,  6.77);
  g_ThdEdxRR_Proton->SetPoint(77,  8.85,  6.87);
  g_ThdEdxRR_Proton->SetPoint(78,  8.55,  6.97);
  g_ThdEdxRR_Proton->SetPoint(79,  8.25,  7.08);
  g_ThdEdxRR_Proton->SetPoint(80,  7.95,  7.19);
  g_ThdEdxRR_Proton->SetPoint(81,  7.65,  7.30);
  g_ThdEdxRR_Proton->SetPoint(82,  7.35,  7.43);
  g_ThdEdxRR_Proton->SetPoint(83,  7.05,  7.56);
  g_ThdEdxRR_Proton->SetPoint(84,  6.75,  7.70);
  g_ThdEdxRR_Proton->SetPoint(85,  6.45,  7.85);
  g_ThdEdxRR_Proton->SetPoint(86,  6.15,  8.02);
  g_ThdEdxRR_Proton->SetPoint(87,  5.85,  8.19);
  g_ThdEdxRR_Proton->SetPoint(88,  5.55,  8.38);
  g_ThdEdxRR_Proton->SetPoint(89,  5.25,  8.58);
  g_ThdEdxRR_Proton->SetPoint(90,  4.95,  8.81);
  g_ThdEdxRR_Proton->SetPoint(91,  4.65,  9.05);
  g_ThdEdxRR_Proton->SetPoint(92,  4.35,  9.32);
  g_ThdEdxRR_Proton->SetPoint(93,  4.05,  9.61);
  g_ThdEdxRR_Proton->SetPoint(94,  3.75,  9.94);
  g_ThdEdxRR_Proton->SetPoint(95,  3.45,  10.32);
  g_ThdEdxRR_Proton->SetPoint(96,  3.15,  10.74);
  g_ThdEdxRR_Proton->SetPoint(97,  2.85,  11.23);
  g_ThdEdxRR_Proton->SetPoint(98,  2.55,  11.80);
  g_ThdEdxRR_Proton->SetPoint(99,  2.25,  12.48);
  g_ThdEdxRR_Proton->SetPoint(100, 1.95,  13.31);
  g_ThdEdxRR_Proton->SetPoint(101, 1.65,  14.35);
  g_ThdEdxRR_Proton->SetPoint(102, 1.35,  15.71);
  g_ThdEdxRR_Proton->SetPoint(103, 1.05,  17.59);
  g_ThdEdxRR_Proton->SetPoint(104, 0.75,  20.44);
  g_ThdEdxRR_Proton->SetPoint(105, 0.45,  25.48);
  g_ThdEdxRR_Proton->SetPoint(106, 0.15,  38.12);

  // Kaon
  g_ThdEdxRR_Kaon = new TGraph(107);
  g_ThdEdxRR_Kaon->SetPoint(0,   31.95, 3.3);
  g_ThdEdxRR_Kaon->SetPoint(1,   31.65, 3.3);
  g_ThdEdxRR_Kaon->SetPoint(2,   31.35, 3.3);
  g_ThdEdxRR_Kaon->SetPoint(3,   31.05, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(4,   30.75, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(5,   30.45, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(6,   30.15, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(7,   29.85, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(8,   29.55, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(9,   29.25, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(10,  28.95, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(11,  28.65, 3.4);
  g_ThdEdxRR_Kaon->SetPoint(12,  28.35, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(13,  28.05, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(14,  27.75, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(15,  27.45, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(16,  27.15, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(17,  26.85, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(18,  26.55, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(19,  26.25, 3.5);
  g_ThdEdxRR_Kaon->SetPoint(20,  25.95, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(21,  25.65, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(22,  25.35, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(23,  25.05, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(24,  24.75, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(25,  24.45, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(26,  24.15, 3.6);
  g_ThdEdxRR_Kaon->SetPoint(27,  23.85, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(28,  23.55, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(29,  23.25, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(30,  22.95, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(31,  22.65, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(32,  22.35, 3.7);
  g_ThdEdxRR_Kaon->SetPoint(33,  22.05, 3.8);
  g_ThdEdxRR_Kaon->SetPoint(34,  21.75, 3.8);
  g_ThdEdxRR_Kaon->SetPoint(35,  21.45, 3.8);
  g_ThdEdxRR_Kaon->SetPoint(36,  21.15, 3.8);
  g_ThdEdxRR_Kaon->SetPoint(37,  20.85, 3.8);
  g_ThdEdxRR_Kaon->SetPoint(38,  20.55, 3.9);
  g_ThdEdxRR_Kaon->SetPoint(39,  20.25, 3.9);
  g_ThdEdxRR_Kaon->SetPoint(40,  19.95, 3.9);
  g_ThdEdxRR_Kaon->SetPoint(41,  19.65, 3.9);
  g_ThdEdxRR_Kaon->SetPoint(42,  19.35, 3.9);
  g_ThdEdxRR_Kaon->SetPoint(43,  19.05, 4.0);
  g_ThdEdxRR_Kaon->SetPoint(44,  18.75, 4.0);
  g_ThdEdxRR_Kaon->SetPoint(45,  18.45, 4.0);
  g_ThdEdxRR_Kaon->SetPoint(46,  18.15, 4.0);
  g_ThdEdxRR_Kaon->SetPoint(47,  17.85, 4.1);
  g_ThdEdxRR_Kaon->SetPoint(48,  17.55, 4.1);
  g_ThdEdxRR_Kaon->SetPoint(49,  17.25, 4.1);
  g_ThdEdxRR_Kaon->SetPoint(50,  16.95, 4.1);
  g_ThdEdxRR_Kaon->SetPoint(51,  16.65, 4.2);
  g_ThdEdxRR_Kaon->SetPoint(52,  16.35, 4.2);
  g_ThdEdxRR_Kaon->SetPoint(53,  16.05, 4.2);
  g_ThdEdxRR_Kaon->SetPoint(54,  15.75, 4.3);
  g_ThdEdxRR_Kaon->SetPoint(55,  15.45, 4.3);
  g_ThdEdxRR_Kaon->SetPoint(56,  15.15, 4.3);
  g_ThdEdxRR_Kaon->SetPoint(57,  14.85, 4.3);
  g_ThdEdxRR_Kaon->SetPoint(58,  14.55, 4.4);
  g_ThdEdxRR_Kaon->SetPoint(59,  14.25, 4.4);
  g_ThdEdxRR_Kaon->SetPoint(60,  13.95, 4.4);
  g_ThdEdxRR_Kaon->SetPoint(61,  13.65, 4.5);
  g_ThdEdxRR_Kaon->SetPoint(62,  13.35, 4.5);
  g_ThdEdxRR_Kaon->SetPoint(63,  13.05, 4.6);
  g_ThdEdxRR_Kaon->SetPoint(64,  12.75, 4.6);
  g_ThdEdxRR_Kaon->SetPoint(65,  12.45, 4.6);
  g_ThdEdxRR_Kaon->SetPoint(66,  12.15, 4.7);
  g_ThdEdxRR_Kaon->SetPoint(67,  11.85, 4.7);
  g_ThdEdxRR_Kaon->SetPoint(68,  11.55, 4.8);
  g_ThdEdxRR_Kaon->SetPoint(69,  11.25, 4.8);
  g_ThdEdxRR_Kaon->SetPoint(70,  10.95, 4.9);
  g_ThdEdxRR_Kaon->SetPoint(71,  10.65, 4.9);
  g_ThdEdxRR_Kaon->SetPoint(72,  10.35, 5.0);
  g_ThdEdxRR_Kaon->SetPoint(73,  10.05, 5.0);
  g_ThdEdxRR_Kaon->SetPoint(74,  9.75,  5.1);
  g_ThdEdxRR_Kaon->SetPoint(75,  9.45,  5.2);
  g_ThdEdxRR_Kaon->SetPoint(76,  9.15,  5.2);
  g_ThdEdxRR_Kaon->SetPoint(77,  8.85,  5.3);
  g_ThdEdxRR_Kaon->SetPoint(78,  8.55,  5.4);
  g_ThdEdxRR_Kaon->SetPoint(79,  8.25,  5.5);
  g_ThdEdxRR_Kaon->SetPoint(80,  7.95,  5.5);
  g_ThdEdxRR_Kaon->SetPoint(81,  7.65,  5.6);
  g_ThdEdxRR_Kaon->SetPoint(82,  7.35,  5.7);
  g_ThdEdxRR_Kaon->SetPoint(83,  7.05,  5.8);
  g_ThdEdxRR_Kaon->SetPoint(84,  6.75,  5.9);
  g_ThdEdxRR_Kaon->SetPoint(85,  6.45,  6.0);
  g_ThdEdxRR_Kaon->SetPoint(86,  6.15,  6.2);
  g_ThdEdxRR_Kaon->SetPoint(87,  5.85,  6.3);
  g_ThdEdxRR_Kaon->SetPoint(88,  5.55,  6.4);
  g_ThdEdxRR_Kaon->SetPoint(89,  5.25,  6.6);
  g_ThdEdxRR_Kaon->SetPoint(90,  4.95,  6.7);
  g_ThdEdxRR_Kaon->SetPoint(91,  4.65,  6.9);
  g_ThdEdxRR_Kaon->SetPoint(92,  4.35,  7.1);
  g_ThdEdxRR_Kaon->SetPoint(93,  4.05,  7.3);
  g_ThdEdxRR_Kaon->SetPoint(94,  3.75,  7.6);
  g_ThdEdxRR_Kaon->SetPoint(95,  3.45,  7.9);
  g_ThdEdxRR_Kaon->SetPoint(96,  3.15,  8.2);
  g_ThdEdxRR_Kaon->SetPoint(97,  2.85,  8.5);
  g_ThdEdxRR_Kaon->SetPoint(98,  2.55,  9.0);
  g_ThdEdxRR_Kaon->SetPoint(99,  2.25,  9.5);
  g_ThdEdxRR_Kaon->SetPoint(100, 1.95,  10.1);
  g_ThdEdxRR_Kaon->SetPoint(101, 1.65,  10.9);
  g_ThdEdxRR_Kaon->SetPoint(102, 1.35,  11.9);
  g_ThdEdxRR_Kaon->SetPoint(103, 1.05,  13.3);
  g_ThdEdxRR_Kaon->SetPoint(104, 0.75,  15.4);
  g_ThdEdxRR_Kaon->SetPoint(105, 0.45,  19.2);
  g_ThdEdxRR_Kaon->SetPoint(106, 0.15,  28.8);


  // Pion
  g_ThdEdxRR_Pion = new TGraph(107);
  g_ThdEdxRR_Pion->SetPoint(0,   31.95, 2.4);
  g_ThdEdxRR_Pion->SetPoint(1,   31.65, 2.4);
  g_ThdEdxRR_Pion->SetPoint(2,   31.35, 2.4);
  g_ThdEdxRR_Pion->SetPoint(3,   31.05, 2.4);
  g_ThdEdxRR_Pion->SetPoint(4,   30.75, 2.4);
  g_ThdEdxRR_Pion->SetPoint(5,   30.45, 2.4);
  g_ThdEdxRR_Pion->SetPoint(6,   30.15, 2.4);
  g_ThdEdxRR_Pion->SetPoint(7,   29.85, 2.4);
  g_ThdEdxRR_Pion->SetPoint(8,   29.55, 2.4);
  g_ThdEdxRR_Pion->SetPoint(9,   29.25, 2.4);
  g_ThdEdxRR_Pion->SetPoint(10,  28.95, 2.4);
  g_ThdEdxRR_Pion->SetPoint(11,  28.65, 2.4);
  g_ThdEdxRR_Pion->SetPoint(12,  28.35, 2.4);
  g_ThdEdxRR_Pion->SetPoint(13,  28.05, 2.4);
  g_ThdEdxRR_Pion->SetPoint(14,  27.75, 2.5);
  g_ThdEdxRR_Pion->SetPoint(15,  27.45, 2.5);
  g_ThdEdxRR_Pion->SetPoint(16,  27.15, 2.5);
  g_ThdEdxRR_Pion->SetPoint(17,  26.85, 2.5);
  g_ThdEdxRR_Pion->SetPoint(18,  26.55, 2.5);
  g_ThdEdxRR_Pion->SetPoint(19,  26.25, 2.5);
  g_ThdEdxRR_Pion->SetPoint(20,  25.95, 2.5);
  g_ThdEdxRR_Pion->SetPoint(21,  25.65, 2.5);
  g_ThdEdxRR_Pion->SetPoint(22,  25.35, 2.5);
  g_ThdEdxRR_Pion->SetPoint(23,  25.05, 2.5);
  g_ThdEdxRR_Pion->SetPoint(24,  24.75, 2.5);
  g_ThdEdxRR_Pion->SetPoint(25,  24.45, 2.5);
  g_ThdEdxRR_Pion->SetPoint(26,  24.15, 2.5);
  g_ThdEdxRR_Pion->SetPoint(27,  23.85, 2.5);
  g_ThdEdxRR_Pion->SetPoint(28,  23.55, 2.5);
  g_ThdEdxRR_Pion->SetPoint(29,  23.25, 2.6);
  g_ThdEdxRR_Pion->SetPoint(30,  22.95, 2.6);
  g_ThdEdxRR_Pion->SetPoint(31,  22.65, 2.6);
  g_ThdEdxRR_Pion->SetPoint(32,  22.35, 2.6);
  g_ThdEdxRR_Pion->SetPoint(33,  22.05, 2.6);
  g_ThdEdxRR_Pion->SetPoint(34,  21.75, 2.6);
  g_ThdEdxRR_Pion->SetPoint(35,  21.45, 2.6);
  g_ThdEdxRR_Pion->SetPoint(36,  21.15, 2.6);
  g_ThdEdxRR_Pion->SetPoint(37,  20.85, 2.6);
  g_ThdEdxRR_Pion->SetPoint(38,  20.55, 2.6);
  g_ThdEdxRR_Pion->SetPoint(39,  20.25, 2.6);
  g_ThdEdxRR_Pion->SetPoint(40,  19.95, 2.6);
  g_ThdEdxRR_Pion->SetPoint(41,  19.65, 2.7);
  g_ThdEdxRR_Pion->SetPoint(42,  19.35, 2.7);
  g_ThdEdxRR_Pion->SetPoint(43,  19.05, 2.7);
  g_ThdEdxRR_Pion->SetPoint(44,  18.75, 2.7);
  g_ThdEdxRR_Pion->SetPoint(45,  18.45, 2.7);
  g_ThdEdxRR_Pion->SetPoint(46,  18.15, 2.7);
  g_ThdEdxRR_Pion->SetPoint(47,  17.85, 2.7);
  g_ThdEdxRR_Pion->SetPoint(48,  17.55, 2.7);
  g_ThdEdxRR_Pion->SetPoint(49,  17.25, 2.8);
  g_ThdEdxRR_Pion->SetPoint(50,  16.95, 2.8);
  g_ThdEdxRR_Pion->SetPoint(51,  16.65, 2.8);
  g_ThdEdxRR_Pion->SetPoint(52,  16.35, 2.8);
  g_ThdEdxRR_Pion->SetPoint(53,  16.05, 2.8);
  g_ThdEdxRR_Pion->SetPoint(54,  15.75, 2.8);
  g_ThdEdxRR_Pion->SetPoint(55,  15.45, 2.8);
  g_ThdEdxRR_Pion->SetPoint(56,  15.15, 2.9);
  g_ThdEdxRR_Pion->SetPoint(57,  14.85, 2.9);
  g_ThdEdxRR_Pion->SetPoint(58,  14.55, 2.9);
  g_ThdEdxRR_Pion->SetPoint(59,  14.25, 2.9);
  g_ThdEdxRR_Pion->SetPoint(60,  13.95, 2.9);
  g_ThdEdxRR_Pion->SetPoint(61,  13.65, 2.9);
  g_ThdEdxRR_Pion->SetPoint(62,  13.35, 3.0);
  g_ThdEdxRR_Pion->SetPoint(63,  13.05, 3.0);
  g_ThdEdxRR_Pion->SetPoint(64,  12.75, 3.0);
  g_ThdEdxRR_Pion->SetPoint(65,  12.45, 3.0);
  g_ThdEdxRR_Pion->SetPoint(66,  12.15, 3.0);
  g_ThdEdxRR_Pion->SetPoint(67,  11.85, 3.1);
  g_ThdEdxRR_Pion->SetPoint(68,  11.55, 3.1);
  g_ThdEdxRR_Pion->SetPoint(69,  11.25, 3.1);
  g_ThdEdxRR_Pion->SetPoint(70,  10.95, 3.1);
  g_ThdEdxRR_Pion->SetPoint(71,  10.65, 3.2);
  g_ThdEdxRR_Pion->SetPoint(72,  10.35, 3.2);
  g_ThdEdxRR_Pion->SetPoint(73,  10.05, 3.2);
  g_ThdEdxRR_Pion->SetPoint(74,  9.75,  3.3);
  g_ThdEdxRR_Pion->SetPoint(75,  9.45,  3.3);
  g_ThdEdxRR_Pion->SetPoint(76,  9.15,  3.3);
  g_ThdEdxRR_Pion->SetPoint(77,  8.85,  3.4);
  g_ThdEdxRR_Pion->SetPoint(78,  8.55,  3.4);
  g_ThdEdxRR_Pion->SetPoint(79,  8.25,  3.4);
  g_ThdEdxRR_Pion->SetPoint(80,  7.95,  3.5);
  g_ThdEdxRR_Pion->SetPoint(81,  7.65,  3.5);
  g_ThdEdxRR_Pion->SetPoint(82,  7.35,  3.6);
  g_ThdEdxRR_Pion->SetPoint(83,  7.05,  3.6);
  g_ThdEdxRR_Pion->SetPoint(84,  6.75,  3.7);
  g_ThdEdxRR_Pion->SetPoint(85,  6.45,  3.7);
  g_ThdEdxRR_Pion->SetPoint(86,  6.15,  3.8);
  g_ThdEdxRR_Pion->SetPoint(87,  5.85,  3.9);
  g_ThdEdxRR_Pion->SetPoint(88,  5.55,  3.9);
  g_ThdEdxRR_Pion->SetPoint(89,  5.25,  4.0);
  g_ThdEdxRR_Pion->SetPoint(90,  4.95,  4.1);
  g_ThdEdxRR_Pion->SetPoint(91,  4.65,  4.2);
  g_ThdEdxRR_Pion->SetPoint(92,  4.35,  4.3);
  g_ThdEdxRR_Pion->SetPoint(93,  4.05,  4.4);
  g_ThdEdxRR_Pion->SetPoint(94,  3.75,  4.6);
  g_ThdEdxRR_Pion->SetPoint(95,  3.45,  4.7);
  g_ThdEdxRR_Pion->SetPoint(96,  3.15,  4.9);
  g_ThdEdxRR_Pion->SetPoint(97,  2.85,  5.1);
  g_ThdEdxRR_Pion->SetPoint(98,  2.55,  5.3);
  g_ThdEdxRR_Pion->SetPoint(99,  2.25,  5.6);
  g_ThdEdxRR_Pion->SetPoint(100, 1.95,  5.9);
  g_ThdEdxRR_Pion->SetPoint(101, 1.65,  6.4);
  g_ThdEdxRR_Pion->SetPoint(102, 1.35,  6.9);
  g_ThdEdxRR_Pion->SetPoint(103, 1.05,  7.7);
  g_ThdEdxRR_Pion->SetPoint(104, 0.75,  8.9);
  g_ThdEdxRR_Pion->SetPoint(105, 0.45,  11.0);
  g_ThdEdxRR_Pion->SetPoint(106, 0.15,  16.5);


  // Muon
  g_ThdEdxRR_Muon = new TGraph(107);
  g_ThdEdxRR_Muon->SetPoint(0,   31.95, 2.3);
  g_ThdEdxRR_Muon->SetPoint(1,   31.65, 2.3);
  g_ThdEdxRR_Muon->SetPoint(2,   31.35, 2.3);
  g_ThdEdxRR_Muon->SetPoint(3,   31.05, 2.3);
  g_ThdEdxRR_Muon->SetPoint(4,   30.75, 2.3);
  g_ThdEdxRR_Muon->SetPoint(5,   30.45, 2.3);
  g_ThdEdxRR_Muon->SetPoint(6,   30.15, 2.3);
  g_ThdEdxRR_Muon->SetPoint(7,   29.85, 2.3);
  g_ThdEdxRR_Muon->SetPoint(8,   29.55, 2.3);
  g_ThdEdxRR_Muon->SetPoint(9,   29.25, 2.3);
  g_ThdEdxRR_Muon->SetPoint(10,  28.95, 2.3);
  g_ThdEdxRR_Muon->SetPoint(11,  28.65, 2.3);
  g_ThdEdxRR_Muon->SetPoint(12,  28.35, 2.3);
  g_ThdEdxRR_Muon->SetPoint(13,  28.05, 2.3);
  g_ThdEdxRR_Muon->SetPoint(14,  27.75, 2.3);
  g_ThdEdxRR_Muon->SetPoint(15,  27.45, 2.3);
  g_ThdEdxRR_Muon->SetPoint(16,  27.15, 2.3);
  g_ThdEdxRR_Muon->SetPoint(17,  26.85, 2.3);
  g_ThdEdxRR_Muon->SetPoint(18,  26.55, 2.3);
  g_ThdEdxRR_Muon->SetPoint(19,  26.25, 2.3);
  g_ThdEdxRR_Muon->SetPoint(20,  25.95, 2.3);
  g_ThdEdxRR_Muon->SetPoint(21,  25.65, 2.4);
  g_ThdEdxRR_Muon->SetPoint(22,  25.35, 2.4);
  g_ThdEdxRR_Muon->SetPoint(23,  25.05, 2.4);
  g_ThdEdxRR_Muon->SetPoint(24,  24.75, 2.4);
  g_ThdEdxRR_Muon->SetPoint(25,  24.45, 2.4);
  g_ThdEdxRR_Muon->SetPoint(26,  24.15, 2.4);
  g_ThdEdxRR_Muon->SetPoint(27,  23.85, 2.4);
  g_ThdEdxRR_Muon->SetPoint(28,  23.55, 2.4);
  g_ThdEdxRR_Muon->SetPoint(29,  23.25, 2.4);
  g_ThdEdxRR_Muon->SetPoint(30,  22.95, 2.4);
  g_ThdEdxRR_Muon->SetPoint(31,  22.65, 2.4);
  g_ThdEdxRR_Muon->SetPoint(32,  22.35, 2.4);
  g_ThdEdxRR_Muon->SetPoint(33,  22.05, 2.4);
  g_ThdEdxRR_Muon->SetPoint(34,  21.75, 2.4);
  g_ThdEdxRR_Muon->SetPoint(35,  21.45, 2.4);
  g_ThdEdxRR_Muon->SetPoint(36,  21.15, 2.4);
  g_ThdEdxRR_Muon->SetPoint(37,  20.85, 2.5);
  g_ThdEdxRR_Muon->SetPoint(38,  20.55, 2.5);
  g_ThdEdxRR_Muon->SetPoint(39,  20.25, 2.5);
  g_ThdEdxRR_Muon->SetPoint(40,  19.95, 2.5);
  g_ThdEdxRR_Muon->SetPoint(41,  19.65, 2.5);
  g_ThdEdxRR_Muon->SetPoint(42,  19.35, 2.5);
  g_ThdEdxRR_Muon->SetPoint(43,  19.05, 2.5);
  g_ThdEdxRR_Muon->SetPoint(44,  18.75, 2.5);
  g_ThdEdxRR_Muon->SetPoint(45,  18.45, 2.5);
  g_ThdEdxRR_Muon->SetPoint(46,  18.15, 2.5);
  g_ThdEdxRR_Muon->SetPoint(47,  17.85, 2.5);
  g_ThdEdxRR_Muon->SetPoint(48,  17.55, 2.6);
  g_ThdEdxRR_Muon->SetPoint(49,  17.25, 2.6);
  g_ThdEdxRR_Muon->SetPoint(50,  16.95, 2.6);
  g_ThdEdxRR_Muon->SetPoint(51,  16.65, 2.6);
  g_ThdEdxRR_Muon->SetPoint(52,  16.35, 2.6);
  g_ThdEdxRR_Muon->SetPoint(53,  16.05, 2.6);
  g_ThdEdxRR_Muon->SetPoint(54,  15.75, 2.6);
  g_ThdEdxRR_Muon->SetPoint(55,  15.45, 2.6);
  g_ThdEdxRR_Muon->SetPoint(56,  15.15, 2.6);
  g_ThdEdxRR_Muon->SetPoint(57,  14.85, 2.7);
  g_ThdEdxRR_Muon->SetPoint(58,  14.55, 2.7);
  g_ThdEdxRR_Muon->SetPoint(59,  14.25, 2.7);
  g_ThdEdxRR_Muon->SetPoint(60,  13.95, 2.7);
  g_ThdEdxRR_Muon->SetPoint(61,  13.65, 2.7);
  g_ThdEdxRR_Muon->SetPoint(62,  13.35, 2.7);
  g_ThdEdxRR_Muon->SetPoint(63,  13.05, 2.8);
  g_ThdEdxRR_Muon->SetPoint(64,  12.75, 2.8);
  g_ThdEdxRR_Muon->SetPoint(65,  12.45, 2.8);
  g_ThdEdxRR_Muon->SetPoint(66,  12.15, 2.8);
  g_ThdEdxRR_Muon->SetPoint(67,  11.85, 2.8);
  g_ThdEdxRR_Muon->SetPoint(68,  11.55, 2.8);
  g_ThdEdxRR_Muon->SetPoint(69,  11.25, 2.9);
  g_ThdEdxRR_Muon->SetPoint(70,  10.95, 2.9);
  g_ThdEdxRR_Muon->SetPoint(71,  10.65, 2.9);
  g_ThdEdxRR_Muon->SetPoint(72,  10.35, 2.9);
  g_ThdEdxRR_Muon->SetPoint(73,  10.05, 3.0);
  g_ThdEdxRR_Muon->SetPoint(74,  9.75,  3.0);
  g_ThdEdxRR_Muon->SetPoint(75,  9.45,  3.0);
  g_ThdEdxRR_Muon->SetPoint(76,  9.15,  3.0);
  g_ThdEdxRR_Muon->SetPoint(77,  8.85,  3.1);
  g_ThdEdxRR_Muon->SetPoint(78,  8.55,  3.1);
  g_ThdEdxRR_Muon->SetPoint(79,  8.25,  3.1);
  g_ThdEdxRR_Muon->SetPoint(80,  7.95,  3.2);
  g_ThdEdxRR_Muon->SetPoint(81,  7.65,  3.2);
  g_ThdEdxRR_Muon->SetPoint(82,  7.35,  3.3);
  g_ThdEdxRR_Muon->SetPoint(83,  7.05,  3.3);
  g_ThdEdxRR_Muon->SetPoint(84,  6.75,  3.4);
  g_ThdEdxRR_Muon->SetPoint(85,  6.45,  3.4);
  g_ThdEdxRR_Muon->SetPoint(86,  6.15,  3.5);
  g_ThdEdxRR_Muon->SetPoint(87,  5.85,  3.5);
  g_ThdEdxRR_Muon->SetPoint(88,  5.55,  3.6);
  g_ThdEdxRR_Muon->SetPoint(89,  5.25,  3.7);
  g_ThdEdxRR_Muon->SetPoint(90,  4.95,  3.7);
  g_ThdEdxRR_Muon->SetPoint(91,  4.65,  3.8);
  g_ThdEdxRR_Muon->SetPoint(92,  4.35,  3.9);
  g_ThdEdxRR_Muon->SetPoint(93,  4.05,  4.0);
  g_ThdEdxRR_Muon->SetPoint(94,  3.75,  4.1);
  g_ThdEdxRR_Muon->SetPoint(95,  3.45,  4.3);
  g_ThdEdxRR_Muon->SetPoint(96,  3.15,  4.4);
  g_ThdEdxRR_Muon->SetPoint(97,  2.85,  4.6);
  g_ThdEdxRR_Muon->SetPoint(98,  2.55,  4.8);
  g_ThdEdxRR_Muon->SetPoint(99,  2.25,  5.0);
  g_ThdEdxRR_Muon->SetPoint(100, 1.95,  5.3);
  g_ThdEdxRR_Muon->SetPoint(101, 1.65,  5.7);
  g_ThdEdxRR_Muon->SetPoint(102, 1.35,  6.2);
  g_ThdEdxRR_Muon->SetPoint(103, 1.05,  6.9);
  g_ThdEdxRR_Muon->SetPoint(104, 0.75,  7.9);
  g_ThdEdxRR_Muon->SetPoint(105, 0.45,  9.8);
  g_ThdEdxRR_Muon->SetPoint(106, 0.15,  14.6);

}
