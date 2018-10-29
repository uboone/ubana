void plotDataMcFromTree(){


  int nbins = 40;
  int binlow = 0;
  int binhigh = 20;
//  double protonScaling = 1.75;
//  double muonScaling = 1.48;

  double protonScaling = 1.563;
  double muonScaling = 1.563; 
  
  //double protonScaling = 1.0;
  //double muonScaling = 1.0;

  //TFile *f_bnbcos  = new TFile("/uboone/data/users/alister1/particleID/180423-ParticleId/pid_bnbcos.root", "read");
  TFile *f_bnbcos  = new TFile("pidtest.root", "read");

  //TFile *f_offbeam = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_offbeam.root", "read");
  //TFile *f_onbeam  = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_onbeam.root", "read");

  TTree *t_bnbcos  = (TTree*)f_bnbcos->Get("pidvalid/pidTree"); 
  //TTree *t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree"); 
  //TTree *t_onbeam  = (TTree*)f_onbeam->Get("pidvalid/pidTree"); 

  std::vector<TString> plotNames = {
    "track_neglogl_fwd_mu",
    "track_neglogl_fwd_p",
    "track_neglogl_fwd_pi",
    "track_neglogl_fwd_k",
    "track_neglogl_fwd_mip",
    "track_PIDA_mean",
    "track_PIDA_kde"};

  for (int i = 0; i < plotNames.size(); i++){
    /**
     * set branch addresses for t_bnbcos
     */

    int true_PDG;
    double variableOfInterest;
    TString plotName = plotNames.at(i);

    if (plotName == "track_PIDA_mean" 
        || plotName == "track_PIDA_median"
        || plotName == "track_PIDA_kde"){

      nbins =  30;
      binhigh = 30;
      protonScaling = 1.0;
      muonScaling = 1.0;

    }

    t_bnbcos->SetBranchAddress("true_PDG", &true_PDG);
    t_bnbcos->SetBranchAddress(plotName, &variableOfInterest);

    TH1D* h_bnbcos_p  = new TH1D("h_bnbcos_p", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_mu = new TH1D("h_bnbcos_mu", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_pi = new TH1D("h_bnbcos_pi", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_k  = new TH1D("h_bnbcos_k", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_other = new TH1D("h_bnbcos_other", ";-2NegLL_p;", nbins, binlow, binhigh);

    for (int i = 0; i < t_bnbcos->GetEntries(); i++){

      t_bnbcos->GetEntry(i);

      if (std::abs(true_PDG) == 2212)
        h_bnbcos_p->Fill(variableOfInterest*1./protonScaling);
      else if (std::abs(true_PDG) == 13)
        h_bnbcos_mu->Fill(variableOfInterest*1./muonScaling);
      else if (std::abs(true_PDG) == 211)
        h_bnbcos_pi->Fill(variableOfInterest*1./muonScaling);
      else if (std::abs(true_PDG) == 321)
        h_bnbcos_k->Fill(variableOfInterest*1./muonScaling);
      else h_bnbcos_other->Fill(variableOfInterest*1./muonScaling); 

    }

    h_bnbcos_p->SetLineWidth(0);
    h_bnbcos_mu->SetLineWidth(0);
    h_bnbcos_pi->SetLineWidth(0);
    h_bnbcos_k->SetLineWidth(0);
    h_bnbcos_other->SetLineWidth(0);
    h_bnbcos_p->SetFillColor(TColor::GetColor(215, 48, 39));
    h_bnbcos_mu->SetFillColor(TColor::GetColor(8,64,129));
    h_bnbcos_pi->SetFillColor(TColor::GetColor(166,217,106));
    h_bnbcos_k->SetFillColor(TColor::GetColor(133,1,98));
    h_bnbcos_other->SetFillColor(TColor::GetColor(197,197,197));
/*
    h_bnbcos_p->Sumw2();
    h_bnbcos_mu->Sumw2();
    h_bnbcos_pi->Sumw2();
    h_bnbcos_k->Sumw2();
    h_bnbcos_other->Sumw2();
*/

    TH1D* h_total = (TH1D*)h_bnbcos_p->Clone("h_total");
    h_total->Add(h_bnbcos_mu);
    h_total->Add(h_bnbcos_pi);
    h_total->Add(h_bnbcos_k);
    h_total->Add(h_bnbcos_other);  

    h_total->Sumw2();

    h_bnbcos_p->Scale(1./h_total->Integral());
    h_bnbcos_mu->Scale(1./h_total->Integral());
    h_bnbcos_pi->Scale(1./h_total->Integral());
    h_bnbcos_k->Scale(1./h_total->Integral());
    h_bnbcos_other->Scale(1./h_total->Integral());
    h_total->Scale(1./h_total->Integral());
    
    THStack *hs = new THStack("hs", "hs");
    hs->Add(h_bnbcos_p);
    hs->Add(h_bnbcos_mu);
    hs->Add(h_bnbcos_pi);
    hs->Add(h_bnbcos_k);
    hs->Add(h_bnbcos_other);

    hs->SetTitle(";"+plotName+";");
/*
    TH1D* h_offbeam = new TH1D("h_offbeam", ";;", nbins, binlow, binhigh);
    t_offbeam->Draw(plotName+" >> h_offbeam");
    h_offbeam->Scale(0.76);
    TH1D* h_onbeam = new TH1D("h_onbeam", ";;", nbins, binlow, binhigh);
    t_onbeam->Draw(plotName+" >> h_onbeam");

    TH1D* h_onminusoff = (TH1D*)h_onbeam->Clone("h_onminusoff");
    h_onminusoff->Add(h_offbeam, -1);
    h_onminusoff->Sumw2();
    h_onminusoff->Scale(1./h_onminusoff->Integral());
*/
    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

    c1->cd();
    //hs->SetMaximum(std::max(h_onminusoff->GetMaximum(), h_total->GetMaximum())*1.15);
    hs->Draw();
    h_total->SetFillStyle(3345);
    h_total->SetFillColor(kGray+2);
    h_total->Draw("sameE2");
  //  h_onminusoff->SetMarkerStyle(20);
  //  h_onminusoff->SetMarkerSize(0.6);
  //  h_onminusoff->Draw("samepE1");

    TLegend *leg = new TLegend(0.55, 0.64, 0.89, 0.89);
    if (plotName == "track_neglogl_p"){
      std::cout << "gotya" << std::endl;
      leg->SetX1(0.15);
      leg->SetX2(0.45);
    }
    leg->AddEntry(h_bnbcos_p, "True proton");
    leg->AddEntry(h_bnbcos_mu, "True Muon");
    leg->AddEntry(h_bnbcos_pi, "True Pion");
    leg->AddEntry(h_bnbcos_k, "True kaon");
    leg->AddEntry(h_bnbcos_other, "True Other");
    //leg->AddEntry(h_onminusoff, "Data");
    leg->SetLineWidth(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    c1->SaveAs(plotName+".png");

    h_bnbcos_p->Delete();
    h_bnbcos_mu->Delete();
    h_bnbcos_pi->Delete();
    h_bnbcos_k->Delete();
    h_bnbcos_other->Delete();
  //  h_offbeam->Delete();
  //  h_onbeam->Delete();

  }
}
