void makeNewPlots(){

  TH1D* hmumu = (TH1D*)_file0->Get("h_trueMuon_muonLikelihoodZeroToOne");
  TH1D* hmup  = (TH1D*)_file0->Get("h_trueMuon_protonLikelihoodZeroToOne");

  TH1D* hpmu = (TH1D*)_file0->Get("h_trueProton_muonLikelihoodZeroToOne");
  TH1D* hpp  = (TH1D*)_file0->Get("h_trueProton_protonLikelihoodZeroToOne");

  TH1D* hmpull = (TH1D*)_file0->Get("h_trueMuon_MipConsistencyPull");
  TH1D* hppull = (TH1D*)_file0->Get("h_trueProton_MipConsistencyPull");

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->cd();
  hmumu->SetLineColor(kAzure+1);
  hmumu->SetFillColor(kAzure+1);
  hmumu->SetFillStyle(3345);
  hmup->SetLineColor(kRed-3);
  hmup->SetFillColor(kRed-3);
  hmup->SetFillStyle(3354);
 
  hmumu->SetTitle(";True Muon Likelihood;");
  
  hmumu->Draw();
  hmup->Draw("same");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();
  hpmu->SetLineColor(kAzure+1);
  hpmu->SetFillColor(kAzure+1);
  hpmu->SetFillStyle(3345);
  hpp->SetLineColor(kRed-3);
  hpp->SetFillColor(kRed-3);
  hpp->SetFillStyle(3354);
  
  hpmu->SetTitle(";True Proton Likelihood;");

  hpmu->Draw();
  hpp->Draw("same");

  TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
  c3->cd();

  hmpull->SetLineColor(kAzure+1);
  hmpull->SetFillColor(kAzure+1);
  hmpull->SetFillStyle(3345);
  hppull->SetLineColor(kRed-3);
  hppull->SetFillColor(kRed-3);
  hppull->SetFillStyle(3354);
  
  hmpull->SetTitle(";proton/muon pulls;");

  hmpull->Draw();
  hppull->Draw("same");



}
