void makePlot(){

  TH1D* hMP = (TH1D*)_file0->Get("hTrueMuonChisqProton");
  TH1D* hMM = (TH1D*)_file0->Get("hTrueMuonChisqMuon");
  TH1D* hPP = (TH1D*)_file0->Get("hTrueProtonChisqProton");
  TH1D* hPM = (TH1D*)_file0->Get("hTrueProtonChisqMuon");

  hMP->Rebin(10);
  hMM->Rebin(10);
//  hPP->Rebin(10);
//  hPM->Rebin(10);

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->cd();
 
  hMM->SetTitle(";Chisq;");

  hMP->SetLineColor(kRed);
  hMM->Draw();
  hMP->Draw("same");

  TLegend *leg1 = new TLegend(0.15, 0.75, 0.5, 0.85);
  leg1->AddEntry(hMM, "Muon under Mu. assmp.");
  leg1->AddEntry(hMP, "Muon under Pr. assmp.");
  
  leg1->Draw("same");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();
 
  hPM->SetTitle(";Chisq;");
  hPM->GetXaxis()->SetRangeUser(-10, 10); 
  hPM->SetLineColor(kRed);
  hPM->Draw();
  hPP->Draw("same");

  TLegend *leg2 = new TLegend(0.15, 0.75, 0.5, 0.85);
  leg2->AddEntry(hPP, "Proton under Pr. assmp.");
  leg2->AddEntry(hPM, "Proton under Mu. assmp.");
  
  leg2->Draw("same");

  c1->SaveAs("trueMuons.png");
  c2->SaveAs("trueProtons.png");

}

