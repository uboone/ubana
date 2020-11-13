void makeTheoryCurvePlot(){

  TGraph* gr_proton = (TGraph*)_file0->Get("gr_proton");
  TGraph* gr_muon = (TGraph*)_file0->Get("gr_muon");
  TGraph* gr_pion = (TGraph*)_file0->Get("gr_pion");
  TGraph* gr_kaon = (TGraph*)_file0->Get("gr_kaon");

  gr_proton->SetLineColor(TColor::GetColor(215, 48, 39));
  gr_muon->SetLineColor(TColor::GetColor(8,64,129));
  gr_pion->SetLineColor(TColor::GetColor(166,217,106));
  gr_kaon->SetLineColor(TColor::GetColor(133,1,98));
  gr_proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  gr_muon->SetMarkerColor(TColor::GetColor(8,64,129));
  gr_pion->SetMarkerColor(TColor::GetColor(166,217,106));
  gr_kaon->SetMarkerColor(TColor::GetColor(133,1,98));
  gr_proton->SetLineWidth(3);
  gr_muon->SetLineWidth(3);
  gr_pion->SetLineWidth(3);
  gr_kaon->SetLineWidth(3);

  gr_proton->SetTitle(";Residual Range (cm);dE/dx (MeV/cm)");
  gr_proton->GetXaxis()->SetRangeUser(0,30);
  gr_proton->GetYaxis()->SetTitleOffset(0.8);

  TCanvas *c1 = new TCanvas("c1", "", 600, 400);

  c1->SetRightMargin (0.07);
  c1->SetLeftMargin(0.1);

  gr_proton->Draw();
  gr_muon->Draw("same");
  gr_pion->Draw("same");
  gr_kaon->Draw("same");

  TLegend *leg = new TLegend(0.75, 0.69, 0.90, 0.89);
  leg->AddEntry(gr_proton, "Proton");
  leg->AddEntry(gr_kaon, "Kaon");
  leg->AddEntry(gr_pion, "Pion");
  leg->AddEntry(gr_muon, "Muon");

  leg->Draw("same");

  c1->SaveAs("theoryCurves.png");

  TCanvas *c2 = new TCanvas("c2", "", 500, 500);
  c2->cd();

  gr_muon->GetXaxis()->SetRangeUser(0,30);
  gr_muon->SetTitle(";Residual Range (cm); dE/dx (MeV/cm)");
  
  TGraph *gr_muon_shiftl = new TGraph(107);
  double shift = 0;
  for (int i = 0; i < 107; i++){

    gr_muon_shiftl->SetPoint(i, i*(30./107.), gr_muon->Eval(i*30./107. - 2.0));

  }

  TGraph *gr_muon_shiftr = new TGraph(107);
  shift = 0;
  for (int i = 0; i < 107; i++){

    gr_muon_shiftr->SetPoint(i, i*(30./107.), gr_muon->Eval(i*30./107. + 2.0));

  }


  double x[107]; double y[107]; double exu[107]; double exd[107]; double eyu[107]; double eyd[107];
  for (int i = 0; i < 107; i++){

    x[i] = i*(30./107.);
    y[i] = gr_muon->Eval(i*30./107.);
    exu[i] = 0;
    exd[i] = 0;
    eyu[i] =  0 - gr_muon->Eval(i*30./107.) + gr_muon->Eval(i*30./107. - 2);
    eyd[i] =  gr_muon->Eval(i*30./107.) - gr_muon->Eval(i*30./107. + 2);

  }

  TGraphAsymmErrors *gr_shade = new TGraphAsymmErrors(107, x, y, exd, exu, eyd, eyu);

  gr_shade->SetFillColor(kGray);
  gr_shade->SetTitle(";Residual Range (cm);dE/dx (MeV/cm)");
  gr_shade->GetYaxis()->SetRangeUser(0,20);
  gr_shade->GetXaxis()->SetRangeUser(0,30);
  gr_shade->Draw("al3");

  gr_muon->SetLineColor(kBlack);
  gr_muon->SetLineWidth(1);

  gr_muon->Draw("samec");

  gr_muon_shiftl->SetLineColor(kGreen+1);
  gr_muon_shiftl->SetMarkerColor(kGreen+1);
  gr_muon_shiftl->SetFillStyle(0);
  gr_muon_shiftl->Draw("samec");

  gr_muon_shiftr->SetLineColor(kRed-3);

  gr_muon_shiftr->SetMarkerColor(kRed-3);
  gr_muon_shiftr->SetFillStyle(0);
  gr_muon_shiftr->Draw("samec");

  TLegend *leg2 = new TLegend(0.4, 0.74, 0.83, 0.89);
  leg2->AddEntry(gr_muon, "Nominal");
  leg2->AddEntry(gr_muon_shiftl, "Extended 2 cm", "l");
  leg2->AddEntry(gr_muon_shiftr, "Truncated 2 cm", "l");

  leg2->Draw("same");

  c2->SaveAs("EndPointResolution.png");
}
