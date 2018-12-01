void plotDataMcFromTree_likelihoodcomparison(){

  int nbins_muminusp = 50;
  int binlow_muminusp = -250;
  int binhigh_muminusp = 150;
  int nbins_mipminusp = 50;
  int binlow_mipminusp = -450;
  int binhigh_mipminusp = 50;
  int nbins_minmumipminusp = 50;
  int binlow_minmumipminusp = -450;
  int binhigh_minmumipminusp = 50;

  //  double protonScaling = 1.0;
  //  double muonScaling = 1.18;

  //  double protonScaling = 1.75;
  //  double muonScaling = 1.9;

  double protonScaling = 1.563;
  double muonScaling = 1.563; 

  //double protonScaling = 1.0;
  //double muonScaling = 1.0;

  TFile *f_bnbcos  = new TFile("/uboone/data/users/alister1/particleID/180423-ParticleId/pid_bnbcos.root", "read");
  TFile *f_offbeam = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_offbeam.root", "read");
  TFile *f_onbeam  = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_onbeam.root", "read");

  TTree *t_bnbcos  = (TTree*)f_bnbcos->Get("pidvalid/pidTree"); 
  TTree *t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree"); 
  TTree *t_onbeam  = (TTree*)f_onbeam->Get("pidvalid/pidTree"); 

  /**
   * set branch addresses for t_bnbcos
   */

  int true_PDG;
  double bnbcos_track_neglogl_mu;
  double bnbcos_track_neglogl_mip;
  double bnbcos_track_neglogl_p;
  double offbeam_track_neglogl_mu;
  double offbeam_track_neglogl_mip;
  double offbeam_track_neglogl_p;
  double onbeam_track_neglogl_mu;
  double onbeam_track_neglogl_mip;
  double onbeam_track_neglogl_p;

  t_bnbcos->SetBranchAddress("true_PDG", &true_PDG);
  t_bnbcos->SetBranchAddress("track_neglogl_mip", &bnbcos_track_neglogl_mip);
  t_bnbcos->SetBranchAddress("track_neglogl_mu", &bnbcos_track_neglogl_mu);
  t_bnbcos->SetBranchAddress("track_neglogl_p", &bnbcos_track_neglogl_p);
  t_offbeam->SetBranchAddress("track_neglogl_mip", &offbeam_track_neglogl_mip);
  t_offbeam->SetBranchAddress("track_neglogl_mu", &offbeam_track_neglogl_mu);
  t_offbeam->SetBranchAddress("track_neglogl_p", &offbeam_track_neglogl_p);
  t_onbeam->SetBranchAddress("track_neglogl_mip", &onbeam_track_neglogl_mip);
  t_onbeam->SetBranchAddress("track_neglogl_mu", &onbeam_track_neglogl_mu);
  t_onbeam->SetBranchAddress("track_neglogl_p", &onbeam_track_neglogl_p);

  TH1D* h_bnbcos_neglogl_muminusp_p  = new TH1D("h_bnbcos_neglogl_muminusp_p", ";-2NegLL_p;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_bnbcos_neglogl_muminusp_mu = new TH1D("h_bnbcos_neglogl_muminusp_mu", ";-2NegLL_p;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_bnbcos_neglogl_muminusp_pi = new TH1D("h_bnbcos_neglogl_muminusp_pi", ";-2NegLL_p;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_bnbcos_neglogl_muminusp_k  = new TH1D("h_bnbcos_neglogl_muminusp_k", ";-2NegLL_p;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_bnbcos_neglogl_muminusp_other = new TH1D("h_bnbcos_neglogl_muminusp_other", ";-2NegLL_p;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_bnbcos_neglogl_mipminusp_p  = new TH1D("h_bnbcos_neglogl_mipminusp_p", ";-2NegLL_p;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_bnbcos_neglogl_mipminusp_mu = new TH1D("h_bnbcos_neglogl_mipminusp_mu", ";-2NegLL_p;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_bnbcos_neglogl_mipminusp_pi = new TH1D("h_bnbcos_neglogl_mipminusp_pi", ";-2NegLL_p;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_bnbcos_neglogl_mipminusp_k  = new TH1D("h_bnbcos_neglogl_mipminusp_k", ";-2NegLL_p;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_bnbcos_neglogl_mipminusp_other = new TH1D("h_bnbcos_neglogl_mipminusp_other", ";-2NegLL_p;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_bnbcos_neglogl_minmumipminusp_p  = new TH1D("h_bnbcos_neglogl_minmumipminusp_p", ";-2NegLL_p;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);
  TH1D* h_bnbcos_neglogl_minmumipminusp_mu = new TH1D("h_bnbcos_neglogl_minmumipminusp_mu", ";-2NegLL_p;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);
  TH1D* h_bnbcos_neglogl_minmumipminusp_pi = new TH1D("h_bnbcos_neglogl_minmumipminusp_pi", ";-2NegLL_p;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);
  TH1D* h_bnbcos_neglogl_minmumipminusp_k  = new TH1D("h_bnbcos_neglogl_minmumipminusp_k", ";-2NegLL_p;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);
  TH1D* h_bnbcos_neglogl_minmumipminusp_other = new TH1D("h_bnbcos_neglogl_minmumipminusp_other", ";-2NegLL_p;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);

  for (int i = 0; i < t_bnbcos->GetEntries(); i++){

    t_bnbcos->GetEntry(i);

    if (std::abs(true_PDG) == 2212){
      h_bnbcos_neglogl_muminusp_p->Fill(bnbcos_track_neglogl_mu*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_mipminusp_p->Fill(bnbcos_track_neglogl_mip*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_minmumipminusp_p->Fill(std::min(bnbcos_track_neglogl_mip, bnbcos_track_neglogl_mu)*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
    }
    else if (std::abs(true_PDG) == 13){
      h_bnbcos_neglogl_muminusp_mu->Fill(bnbcos_track_neglogl_mu*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_mipminusp_mu->Fill(bnbcos_track_neglogl_mip*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_minmumipminusp_mu->Fill(std::min(bnbcos_track_neglogl_mip, bnbcos_track_neglogl_mu)*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
    }
    else if (std::abs(true_PDG) == 211){
      h_bnbcos_neglogl_muminusp_pi->Fill(bnbcos_track_neglogl_mu*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_mipminusp_pi->Fill(bnbcos_track_neglogl_mip*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_minmumipminusp_pi->Fill(std::min(bnbcos_track_neglogl_mip, bnbcos_track_neglogl_mu)*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
    }
    else if (std::abs(true_PDG) == 321){
      h_bnbcos_neglogl_muminusp_k->Fill(bnbcos_track_neglogl_mu*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_mipminusp_k->Fill(bnbcos_track_neglogl_mip*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
      h_bnbcos_neglogl_minmumipminusp_k->Fill(std::min(bnbcos_track_neglogl_mip, bnbcos_track_neglogl_mu)*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling);
    }
    else{
      h_bnbcos_neglogl_muminusp_other->Fill(bnbcos_track_neglogl_mu*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling); 
      h_bnbcos_neglogl_mipminusp_other->Fill(bnbcos_track_neglogl_mip*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling); 
      h_bnbcos_neglogl_minmumipminusp_other->Fill(std::min(bnbcos_track_neglogl_mip, bnbcos_track_neglogl_mu)*1./muonScaling - bnbcos_track_neglogl_p*1./protonScaling); 
    }
  }

  TH1D* h_offbeam_neglogl_muminusp = new TH1D("h_offbeam_neglogl_muminusp", ";;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_onbeam_neglogl_muminusp = new TH1D("h_onbeam_neglogl_muminusp", ";;", nbins_muminusp, binlow_muminusp, binhigh_muminusp);
  TH1D* h_offbeam_neglogl_mipminusp = new TH1D("h_offbeam_neglogl_mipminusp", ";;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_onbeam_neglogl_mipminusp = new TH1D("h_onbeam_neglogl_mipminusp", ";;", nbins_mipminusp, binlow_mipminusp, binhigh_mipminusp);
  TH1D* h_offbeam_neglogl_minmumipminusp = new TH1D("h_offbeam_neglogl_minmumipminusp", ";;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);
  TH1D* h_onbeam_neglogl_minmumipminusp = new TH1D("h_onbeam_neglogl_minmumipminusp", ";;", nbins_minmumipminusp, binlow_minmumipminusp, binhigh_minmumipminusp);

  for (int i = 0; i < t_offbeam->GetEntries(); i++){

    t_offbeam->GetEntry(i);

    h_offbeam_neglogl_muminusp->Fill(offbeam_track_neglogl_mu - offbeam_track_neglogl_p);
    h_offbeam_neglogl_mipminusp->Fill(offbeam_track_neglogl_mip - offbeam_track_neglogl_p);
    h_offbeam_neglogl_minmumipminusp->Fill(std::min(offbeam_track_neglogl_mu, offbeam_track_neglogl_mip) - offbeam_track_neglogl_p);
  }

  for (int i = 0; i < t_onbeam->GetEntries(); i++){

    t_onbeam->GetEntry(i);

    h_onbeam_neglogl_muminusp->Fill(onbeam_track_neglogl_mu - onbeam_track_neglogl_p);
    h_onbeam_neglogl_mipminusp->Fill(onbeam_track_neglogl_mip - onbeam_track_neglogl_p);
    h_onbeam_neglogl_minmumipminusp->Fill(std::min(onbeam_track_neglogl_mu, onbeam_track_neglogl_mip) - onbeam_track_neglogl_p);

  }

  /**
   * muminusp
   */
  h_bnbcos_neglogl_muminusp_p->SetLineWidth(0);
  h_bnbcos_neglogl_muminusp_mu->SetLineWidth(0);
  h_bnbcos_neglogl_muminusp_pi->SetLineWidth(0);
  h_bnbcos_neglogl_muminusp_k->SetLineWidth(0);
  h_bnbcos_neglogl_muminusp_other->SetLineWidth(0);
  h_bnbcos_neglogl_muminusp_p->SetFillColor(TColor::GetColor(215, 48, 39));
  h_bnbcos_neglogl_muminusp_mu->SetFillColor(TColor::GetColor(8,64,129));
  h_bnbcos_neglogl_muminusp_pi->SetFillColor(TColor::GetColor(166,217,106));
  h_bnbcos_neglogl_muminusp_k->SetFillColor(TColor::GetColor(133,1,98));
  h_bnbcos_neglogl_muminusp_other->SetFillColor(TColor::GetColor(197,197,197));

  TH1D* h_total_neglogl_muminusp = (TH1D*)h_bnbcos_neglogl_muminusp_p->Clone("h_total_neglogl_muminusp");
  h_total_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_mu);
  h_total_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_pi);
  h_total_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_k);
  h_total_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_other);  
  h_total_neglogl_muminusp->Sumw2();

  h_bnbcos_neglogl_muminusp_p->Scale(1./h_total_neglogl_muminusp->Integral());
  h_bnbcos_neglogl_muminusp_mu->Scale(1./h_total_neglogl_muminusp->Integral());
  h_bnbcos_neglogl_muminusp_pi->Scale(1./h_total_neglogl_muminusp->Integral());
  h_bnbcos_neglogl_muminusp_k->Scale(1./h_total_neglogl_muminusp->Integral());
  h_bnbcos_neglogl_muminusp_other->Scale(1./h_total_neglogl_muminusp->Integral());
  h_total_neglogl_muminusp->Scale(1./h_total_neglogl_muminusp->Integral()); 

  THStack *hs_neglogl_muminusp = new THStack("s_neglogl_muminusp", ";s_neglogl_muminusp;");
  hs_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_p);
  hs_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_mu);
  hs_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_pi);
  hs_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_k);
  hs_neglogl_muminusp->Add(h_bnbcos_neglogl_muminusp_other);

  h_offbeam_neglogl_muminusp->Scale(0.76);

  TH1D* h_onminusoff_neglogl_muminusp = (TH1D*)h_onbeam_neglogl_muminusp->Clone("h_onminusoff_neglogl_muminusp");
  h_onminusoff_neglogl_muminusp->Add(h_offbeam_neglogl_muminusp, -1);
  h_onminusoff_neglogl_muminusp->Sumw2();

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  c1->cd();
  h_total_neglogl_muminusp->SetNdivisions(505);
  h_total_neglogl_muminusp->Draw();
  h_total_neglogl_muminusp->SetTitle(";2NegLL_p - 2NegLL_mu;"); 
  hs_neglogl_muminusp->Draw("same");
  h_total_neglogl_muminusp->SetFillColor(kGray+2);
  h_total_neglogl_muminusp->SetFillStyle(3345);
  h_total_neglogl_muminusp->Draw("E2same");
  h_onminusoff_neglogl_muminusp->SetMarkerStyle(20);
  h_onminusoff_neglogl_muminusp->SetMarkerSize(0.6);
  h_onminusoff_neglogl_muminusp->DrawNormalized("samepE0");

  TLegend *leg_muminusp = new TLegend(0.55, 0.64, 0.89, 0.89);
  leg_muminusp->AddEntry(h_bnbcos_neglogl_muminusp_p, "True proton");
  leg_muminusp->AddEntry(h_bnbcos_neglogl_muminusp_mu, "True Muon");
  leg_muminusp->AddEntry(h_bnbcos_neglogl_muminusp_pi, "True Pion");
  leg_muminusp->AddEntry(h_bnbcos_neglogl_muminusp_k, "True kaon");
  leg_muminusp->AddEntry(h_bnbcos_neglogl_muminusp_other, "True Other");
  leg_muminusp->AddEntry(h_onminusoff_neglogl_muminusp, "Data");
  leg_muminusp->SetLineWidth(0);
  leg_muminusp->SetFillStyle(0);
  leg_muminusp->Draw("same");

  c1->SaveAs("muminusp.png");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();
  h_bnbcos_neglogl_muminusp_mu->SetFillStyle(3354);
  h_bnbcos_neglogl_muminusp_mu->Draw();
  h_bnbcos_neglogl_muminusp_p->SetFillStyle(3345);
  h_bnbcos_neglogl_muminusp_p->Draw("same");

  /**
   * mipminusp
   */
  h_bnbcos_neglogl_mipminusp_p->SetLineWidth(0);
  h_bnbcos_neglogl_mipminusp_mu->SetLineWidth(0);
  h_bnbcos_neglogl_mipminusp_pi->SetLineWidth(0);
  h_bnbcos_neglogl_mipminusp_k->SetLineWidth(0);
  h_bnbcos_neglogl_mipminusp_other->SetLineWidth(0);
  h_bnbcos_neglogl_mipminusp_p->SetFillColor(TColor::GetColor(215, 48, 39));
  h_bnbcos_neglogl_mipminusp_mu->SetFillColor(TColor::GetColor(8,64,129));
  h_bnbcos_neglogl_mipminusp_pi->SetFillColor(TColor::GetColor(166,217,106));
  h_bnbcos_neglogl_mipminusp_k->SetFillColor(TColor::GetColor(133,1,98));
  h_bnbcos_neglogl_mipminusp_other->SetFillColor(TColor::GetColor(197,197,197));

  TH1D* h_total_neglogl_mipminusp = (TH1D*)h_bnbcos_neglogl_mipminusp_p->Clone("h_total_neglogl_mipminusp");
  h_total_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_mu);
  h_total_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_pi);
  h_total_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_k);
  h_total_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_other);  

  h_total_neglogl_mipminusp->Sumw2();

  h_bnbcos_neglogl_mipminusp_p->Scale(1./h_total_neglogl_mipminusp->Integral());
  h_bnbcos_neglogl_mipminusp_mu->Scale(1./h_total_neglogl_mipminusp->Integral());
  h_bnbcos_neglogl_mipminusp_pi->Scale(1./h_total_neglogl_mipminusp->Integral());
  h_bnbcos_neglogl_mipminusp_k->Scale(1./h_total_neglogl_mipminusp->Integral());
  h_bnbcos_neglogl_mipminusp_other->Scale(1./h_total_neglogl_mipminusp->Integral());
  h_total_neglogl_mipminusp->Scale(1./h_total_neglogl_mipminusp->Integral());

  THStack *hs_neglogl_mipminusp = new THStack("hs_neglogl_mipminusp", ";hs_neglogl_mipminusp;");
  hs_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_p);
  hs_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_mu);
  hs_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_pi);
  hs_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_k);
  hs_neglogl_mipminusp->Add(h_bnbcos_neglogl_mipminusp_other);

  h_offbeam_neglogl_mipminusp->Scale(0.76);

  TH1D* h_onminusoff_neglogl_mipminusp = (TH1D*)h_onbeam_neglogl_mipminusp->Clone("h_onminusoff_neglogl_mipminusp");
  h_onminusoff_neglogl_mipminusp->Add(h_offbeam_neglogl_mipminusp, -1);
  h_onminusoff_neglogl_mipminusp->Sumw2();

  TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);

  c3->cd();
  h_total_neglogl_mipminusp->SetNdivisions(505);
  h_total_neglogl_mipminusp->Draw();
  h_total_neglogl_mipminusp->SetTitle(";2NegLL_p - 2NegLL_mip;");
  hs_neglogl_mipminusp->Draw("same");
  h_total_neglogl_mipminusp->SetFillColor(kGray+2);
  h_total_neglogl_mipminusp->SetFillStyle(3345);
  h_total_neglogl_mipminusp->Draw("E2same");
  h_onminusoff_neglogl_mipminusp->SetMarkerStyle(20);
  h_onminusoff_neglogl_mipminusp->SetMarkerSize(0.6);
  h_onminusoff_neglogl_mipminusp->DrawNormalized("samepE0");

  TLegend *leg_mipminusp = new TLegend(0.15, 0.64, 0.40, 0.89);
  leg_mipminusp->AddEntry(h_bnbcos_neglogl_mipminusp_p, "True proton");
  leg_mipminusp->AddEntry(h_bnbcos_neglogl_mipminusp_mu, "True Muon");
  leg_mipminusp->AddEntry(h_bnbcos_neglogl_mipminusp_pi, "True Pion");
  leg_mipminusp->AddEntry(h_bnbcos_neglogl_mipminusp_k, "True kaon");
  leg_mipminusp->AddEntry(h_bnbcos_neglogl_mipminusp_other, "True Other");
  leg_mipminusp->AddEntry(h_onminusoff_neglogl_mipminusp, "Data");
  leg_mipminusp->SetLineWidth(0);
  leg_mipminusp->SetFillStyle(0);
  leg_mipminusp->Draw("same");

  c3->SaveAs("mipminusp.png");

  TCanvas *c4 = new TCanvas("c4", "c4", 500, 500);
  c4->cd();
  h_bnbcos_neglogl_mipminusp_mu->SetFillStyle(3354);
  h_bnbcos_neglogl_mipminusp_mu->Draw();
  h_bnbcos_neglogl_mipminusp_p->SetFillStyle(3345);
  h_bnbcos_neglogl_mipminusp_p->Draw("same");

  /**
   * minmipminusp
   */
  h_bnbcos_neglogl_minmumipminusp_p->SetLineWidth(0);
  h_bnbcos_neglogl_minmumipminusp_mu->SetLineWidth(0);
  h_bnbcos_neglogl_minmumipminusp_pi->SetLineWidth(0);
  h_bnbcos_neglogl_minmumipminusp_k->SetLineWidth(0);
  h_bnbcos_neglogl_minmumipminusp_other->SetLineWidth(0);
  h_bnbcos_neglogl_minmumipminusp_p->SetFillColor(TColor::GetColor(215, 48, 39));
  h_bnbcos_neglogl_minmumipminusp_mu->SetFillColor(TColor::GetColor(8,64,129));
  h_bnbcos_neglogl_minmumipminusp_pi->SetFillColor(TColor::GetColor(166,217,106));
  h_bnbcos_neglogl_minmumipminusp_k->SetFillColor(TColor::GetColor(133,1,98));
  h_bnbcos_neglogl_minmumipminusp_other->SetFillColor(TColor::GetColor(197,197,197));

  TH1D* h_total_neglogl_minmumipminusp = (TH1D*)h_bnbcos_neglogl_minmumipminusp_p->Clone("h_total_neglogl_minmumipminusp");
  h_total_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_mu);
  h_total_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_pi);
  h_total_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_k);
  h_total_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_other);  
  
  h_total_neglogl_minmumipminusp->Sumw2();
  h_bnbcos_neglogl_minmumipminusp_p->Scale(1./h_total_neglogl_minmumipminusp->Integral());
  h_bnbcos_neglogl_minmumipminusp_mu->Scale(1./h_total_neglogl_minmumipminusp->Integral());
  h_bnbcos_neglogl_minmumipminusp_pi->Scale(1./h_total_neglogl_minmumipminusp->Integral());
  h_bnbcos_neglogl_minmumipminusp_k->Scale(1./h_total_neglogl_minmumipminusp->Integral());
  h_bnbcos_neglogl_minmumipminusp_other->Scale(1./h_total_neglogl_minmumipminusp->Integral());
  h_total_neglogl_minmumipminusp->Scale(1./h_total_neglogl_minmumipminusp->Integral());

  THStack *hs_neglogl_minmumipminusp = new THStack("hs_neglogl_minmumipminusp", ";hs_neglogl_minmumipminusp;");
  hs_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_p);
  hs_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_mu);
  hs_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_pi);
  hs_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_k);
  hs_neglogl_minmumipminusp->Add(h_bnbcos_neglogl_minmumipminusp_other);

  h_offbeam_neglogl_minmumipminusp->Scale(0.76);

  TH1D* h_onminusoff_neglogl_minmumipminusp = (TH1D*)h_onbeam_neglogl_minmumipminusp->Clone("h_onminusoff_neglogl_minmumipminusp");
  h_onminusoff_neglogl_minmumipminusp->Add(h_offbeam_neglogl_minmumipminusp, -1);
  h_onminusoff_neglogl_minmumipminusp->Sumw2();

  TCanvas *c5 = new TCanvas("c5", "c5", 500, 500);

  c5->cd();
  h_total_neglogl_minmumipminusp->SetNdivisions(505);
  h_total_neglogl_minmumipminusp->Draw();
  h_total_neglogl_minmumipminusp->SetTitle(";2NegLL_p + min(-2NegLL_mip, -2NegLL_mu);");
  hs_neglogl_minmumipminusp->Draw("same");
  h_total_neglogl_minmumipminusp->SetFillColor(kGray+2);
  h_total_neglogl_minmumipminusp->SetFillStyle(3345);
  h_total_neglogl_minmumipminusp->Draw("sameE2");
  h_onminusoff_neglogl_minmumipminusp->SetMarkerStyle(20);
  h_onminusoff_neglogl_minmumipminusp->SetMarkerSize(0.6);
  h_onminusoff_neglogl_minmumipminusp->DrawNormalized("samepE0");

  int binSeparator = 42;

  double fullintegral_p = h_bnbcos_neglogl_minmumipminusp_p->Integral(-1, 51);
  double integralleft_p = h_bnbcos_neglogl_minmumipminusp_p->Integral(-1, binSeparator);
  double integralright_p = h_bnbcos_neglogl_minmumipminusp_p->Integral(binSeparator+1, 51);
  double fullintegral_mu = h_bnbcos_neglogl_minmumipminusp_mu->Integral(-1, 51);
  double integralleft_mu = h_bnbcos_neglogl_minmumipminusp_mu->Integral(-1, binSeparator);
  double integralright_mu = h_bnbcos_neglogl_minmumipminusp_mu->Integral(binSeparator+1, 51);
  double fullintegral_pi = h_bnbcos_neglogl_minmumipminusp_pi->Integral(-1, 51);
  double integralleft_pi = h_bnbcos_neglogl_minmumipminusp_pi->Integral(-1, binSeparator);
  double integralright_pi = h_bnbcos_neglogl_minmumipminusp_pi->Integral(binSeparator+1, 51);
  double fullintegral_total = h_total_neglogl_minmumipminusp->Integral(-1, 51);
  double integralleft_total = h_total_neglogl_minmumipminusp->Integral(-1, binSeparator);
  double integralright_total = h_total_neglogl_minmumipminusp->Integral(binSeparator+1, 51);

  std::cout << h_bnbcos_neglogl_minmumipminusp_p->GetBinLowEdge(42) << std::endl;
  std::cout << (integralleft_mu+integralright_mu)/fullintegral_mu << std::endl;
  std::cout << "fraction of mips ID'd by cut" << (integralleft_mu + integralleft_pi) /(fullintegral_mu+fullintegral_pi) << std::endl;
  std::cout << "contamination in MIP region" << 1 - (integralleft_mu+integralleft_pi)/(integralleft_total) << std::endl;
  std::cout << "fraction of protons ID'd by cut" << (integralright_p) /(fullintegral_p) << std::endl;
  std::cout << "contamination in proton region" << 1 - (integralright_p)/(integralright_total) << std::endl;
 

  TLegend *leg_minmumipminusp = new TLegend(0.15, 0.64, 0.45, 0.89);
  leg_minmumipminusp->AddEntry(h_bnbcos_neglogl_minmumipminusp_p, "True proton");
  leg_minmumipminusp->AddEntry(h_bnbcos_neglogl_minmumipminusp_mu, "True Muon");
  leg_minmumipminusp->AddEntry(h_bnbcos_neglogl_minmumipminusp_pi, "True Pion");
  leg_minmumipminusp->AddEntry(h_bnbcos_neglogl_minmumipminusp_k, "True kaon");
  leg_minmumipminusp->AddEntry(h_bnbcos_neglogl_minmumipminusp_other, "True Other");
  leg_minmumipminusp->AddEntry(h_onminusoff_neglogl_minmumipminusp, "Data");
  leg_minmumipminusp->SetLineWidth(0);
  leg_minmumipminusp->SetFillStyle(0);
  leg_minmumipminusp->Draw("same");

  c5->SaveAs("minmumipminusp.png");

  TCanvas *c6 = new TCanvas("c6", "c6", 500, 500);
  c6->cd();
  h_bnbcos_neglogl_minmumipminusp_mu->SetFillStyle(3354);
  h_bnbcos_neglogl_minmumipminusp_mu->Draw();
  h_bnbcos_neglogl_minmumipminusp_p->SetFillStyle(3345);
  h_bnbcos_neglogl_minmumipminusp_p->Draw("same");

}
