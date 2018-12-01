/**
 * Plotting script to plot a profile TGraphAsymmErrors (produced by 
 * LandauGaussianPlot.C) on top of a dE/dx versus residual range histogram
 * from simulation or data
 */

void plotOverlay(){

  TCanvas *c1 = new TCanvas("c1", "", 500, 500);

  TH2D* h = (TH2D*)_file0->Get("hdEdx_rr_fullrange_all");
  h->GetYaxis()->SetRangeUser(0,20);
  h->SetContour(10000);

  TGraphAsymmErrors *gr_proton      = (TGraphAsymmErrors*)_file1->Get("gr_proton_profile_x");
  TGraphAsymmErrors *gr_muon        = (TGraphAsymmErrors*)_file1->Get("gr_muon_profile_x");
  TGraphAsymmErrors *gr_muonnobragg = (TGraphAsymmErrors*)_file1->Get("gr_muonnobragg_profile_x");

  gr_proton->SetLineColor(kRed);
  gr_proton->SetMarkerColor(kRed);
  gr_proton->SetMarkerStyle(20);
  gr_proton->SetMarkerSize(0.6);
  gr_muon->SetLineColor(kRed);
  gr_muon->SetMarkerColor(kRed);
  gr_muon->SetMarkerStyle(20);
  gr_muon->SetMarkerSize(0.6);
  gr_muonnobragg->SetLineColor(kRed);
  gr_muonnobragg->SetMarkerColor(kRed);
  gr_muonnobragg->SetMarkerStyle(20);
  gr_muonnobragg->SetMarkerSize(0.6);

  h->Draw("col");
  gr_proton->Draw("p");
  gr_muon->Draw("samep");
  gr_muonnobragg->Draw("samep");

  c1->SaveAs("overlay.png");

}
