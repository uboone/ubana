void makeChisqSeparatorPlots(){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->cd();

  TH2D* h1 = (TH2D*)_file0->Get("hTrueMuonChisqProtonChisqMuon");
  TH2D* h2 = (TH2D*)_file0->Get("hTrueProtonChisqProtonChisqMuon");

  h1->SetMarkerColor(kAzure+1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.4);
  h1->GetXaxis()->SetRangeUser(0, 15);
  h1->GetYaxis()->SetRangeUser(0, 15);

  h2->SetMarkerColor(kRed-3);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.4);

  h1->Draw();
  h2->Draw("same");

  TF1* f = new TF1("f", "x", -5, 30);
  f->SetLineColor(kGray);
  f->SetLineStyle(2);
  f->Draw("same");
  h1->SetLineWidth(0);
  h2->SetLineWidth(0);

  TLegend* leg = new TLegend(0.5, 0.75, 0.85, 0.85);
  leg->AddEntry(h1, "True Muon");
  leg->AddEntry(h2, "True Proton");
  leg->Draw("samw");


  c1->SaveAs("Chisq2Dplot.png");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();

  TH1D* h3 = (TH1D*)_file0->Get("hTrueProtonMuonChisqRatio");
  h3->GetXaxis()->SetRangeUser(0,1);
  h3->GetXaxis()->SetTitle("Proton-Like Score");
  h3->SetFillColor(kRed-3);
  h3->SetFillStyle(3345);
  h3->Draw();

  TH1D* h4 = (TH1D*)_file0->Get("hTrueMuonMuonChisqRatio");
  h4->SetFillColor(kAzure+1);
  h4->SetFillStyle(3354);
  h4->Draw("same");

  TLegend *leg2 = new TLegend(0.15, 0.75, 0.5, 0.85);
  leg2->AddEntry(h3, "True Proton");
  leg2->AddEntry(h4, "True Muon");
  leg2->Draw("same");

  c2->SaveAs("protonLikeScore.png");

  TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
  c3->cd();

  THStack *hs = new THStack("hs", "");
  hs->Add(h3);
  hs->Add(h4);
  hs->Draw();
  hs->GetXaxis()->SetRangeUser(0,1);
  hs->GetXaxis()->SetTitle("Proton-Like Score");
  hs->Draw();
  leg2->Draw("same");
  c3->SaveAs("protonLikeScoreStack.png");

  TCanvas *c4 = new TCanvas("c4", "c4", 500, 500);
  TH2D* h5 = (TH2D*)_file0->Get("hTrueMuonMuonChisqRatioRecoLength");
  TH2D* h6 = (TH2D*)_file0->Get("hTrueProtonMuonChisqRatioRecoLength");
  h5->SetMarkerColor(kAzure+1);
  h5->SetMarkerStyle(20);
  h5->SetMarkerSize(0.4);
  h6->SetMarkerColor(kRed-3);
  h6->SetMarkerStyle(20);
  h6->SetMarkerSize(0.4);
  h5->Draw();
  h6->Draw("same");

  c4->SaveAs("protonLikeScoreVsRecoLength.png");
}
