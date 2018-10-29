void getInterval1D(TH1D *hist, TH1D &h68, TH1D &h90)
{
  std::cout << "getting interval" << std::endl;
  TH1D *hCopy = (TH1D*)hist->Clone("hCopy");
  TH1D *hCopy68 = (TH1D*)hist->Clone("hCopy68");
  TH1D *hCopy90 = (TH1D*)hist->Clone("hCopy90");

  hCopy68->Reset();
  hCopy90->Reset();

  double integral = hCopy->Integral();
  double tsum = 0;

  //std::cout << integral << std::endl;

  while((tsum / integral) < 0.9)
    {
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / integral < 0.68)
	{
	  hCopy->SetBinContent(bin, 0.0);
	  hCopy68->SetBinContent(bin, tmax);
	  hCopy90->SetBinContent(bin, tmax);
	}
      if ((tsum / integral < 0.9) && (tsum / integral > 0.68))
	{
	  hCopy->SetBinContent(bin, 0.0);
	  hCopy90->SetBinContent(bin, tmax);
	}
    }

  h90 = (*hCopy90);
  h68 = (*hCopy68);
}

void getInterval1D(TH1D *hist){

  TH1D *h68 = (TH1D*)hist->Clone("h68");
  TH1D *h90 = (TH1D*)hist->Clone("h90");

  getInterval1D(hist,*h68,*h90);

  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->SetFillStyle(0);
  //h68->SetLineColor(kBlack);
  h68->SetLineWidth(0);
  h68->SetFillColor(kGray);
  h68->SetFillStyle(1001);
  //h90->SetLineColor(kBlack);
  h90->SetLineWidth(0);
  h90->SetFillColor(kGray+2);
  h90->SetFillStyle(1001);

  gStyle->SetOptStat(0);
  h90->Draw("hist");
  h68->Draw("hist same");
  hist->Draw("hist e1 same");

  // Get interval limits and print to screen
  bool foundlowerlimit68 = false;
  bool foundlowerlimit90 = false;
  bool foundupperlimit68 = false;
  bool foundupperlimit90 = false;
  double ll68 = 0., ll90 = 0., ul68 = 0., ul90 = 0.;
  for (int i_bin=1; i_bin<h68->GetXaxis()->GetNbins()+1; i_bin++){
    if (!foundlowerlimit68){
      if (h68->GetBinContent(i_bin)>0) {
	ll68 = h68->GetXaxis()->GetBinLowEdge(i_bin);
	foundlowerlimit68 = true;
      }
    }
    else if (!foundupperlimit68){
      if (h68->GetBinContent(i_bin)==0){
	ul68 = h68->GetXaxis()->GetBinLowEdge(i_bin);
	foundupperlimit68 = true;
      }
    }
    
    if (!foundlowerlimit90){
      if (h90->GetBinContent(i_bin)>0){
	ll90 = h90->GetXaxis()->GetBinLowEdge(i_bin);
	foundlowerlimit90 = true;
      }
    }
    else if (!foundupperlimit90){
      if (h90->GetBinContent(i_bin)==0){
	ul90 = h90->GetXaxis()->GetBinLowEdge(i_bin);
	foundupperlimit90 = true;
      }
    }
  }

  TString entry68;
  entry68.Form("68%%: %.1f - %.1f",ll68,ul68);
  TString entry90;
  entry90.Form("90%%: %.1f - %.1f",ll90,ul90);
  
  TLegend *l = new TLegend(0.59,0.64,0.88,0.87);
  l->SetTextFont(132);
  l->SetFillStyle(0);
  l->SetLineColor(kBlack);
  l->AddEntry(hist,hist->GetTitle(),"l");
  l->AddEntry(h68,entry68.Data(),"f");
  l->AddEntry(h90,entry90.Data(),"f");
  l->Draw();

  std::cout << "To double-check: " << std::endl
	    << "   68% credible interval integral is " << h68->Integral()/hist->Integral() << std::endl
	    << "   90% credible interval integral is " << h90->Integral()/hist->Integral() << std::endl;
  
  std::cout << "Found 68% credible interval: " << ll68 << " - " << ul68 << std::endl;
  std::cout << "Found 90% credible interval: " << ll90 << " - " << ul90 << std::endl;
  
}
