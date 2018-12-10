#include "plotFromTreeHeader.h"

void plotNewDataMcFromTree(std::string mcfile, std::string offbeamdatafile, std::string onbeamdatafile){


  int nbins = 40;
  int binlow = 0;
  int binhigh = 15;

  //double protonScaling = 1.563;
  //double muonScaling = 1.563;

  double protonScaling = 1.0;
  double muonScaling = 1.0;

  double offbeamScaling = 0.771; //0.78

  TFile *f_bnbcos  = new TFile(mcfile.c_str(), "read");
  TFile *f_offbeam  = new TFile(offbeamdatafile.c_str(), "read");
  TFile *f_onbeam  = new TFile(onbeamdatafile.c_str(), "read");

  TTree *t_bnbcos  = (TTree*)f_bnbcos->Get("pidvalid/pidTree");
  TTree *t_onbeam = (TTree*)f_onbeam->Get("pidvalid/pidTree");
  TTree *t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree");

  std::vector<TString> plotNames = {
    "track_neglogl_p",
    "track_neglogl_mu",
    "track_neglogl_pi",
    "track_neglogl_k",
    "track_neglogl_mip",
    "track_Chi2Muon",
    "track_Chi2Proton",
    "track_Chi2Kaon",
    "track_Chi2Pion"};

    treevars bnbcos_vars;
    settreevars(t_bnbcos, &bnbcos_vars);
    treevars offbeam_vars;
    settreevars(t_offbeam, &offbeam_vars);
    treevars onbeam_vars;
    settreevars(t_onbeam, &onbeam_vars);

    // what are these for?
    double track_PIDA_mean;
    double track_PIDA_kde;

  for (int j = 0; j < plotNames.size(); j++){

    if (j==6 || j==7){ // Chi2 proton and kaon plots
      nbins = 30;
      binlow = 0;
      binhigh = 300;
    }
    if (j==5 || j==8){ // Chi2 muon and pion plots
      nbins = 25;
      binlow = 0;
      binhigh = 125;
    }

    TH1D* h_bnbcos_p  = new TH1D("h_bnbcos_p", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_mu = new TH1D("h_bnbcos_mu", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_pi = new TH1D("h_bnbcos_pi", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_k  = new TH1D("h_bnbcos_k", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_other = new TH1D("h_bnbcos_other", ";-2NegLL_p;", nbins, binlow, binhigh);

    TH1D* h_offbeam = new TH1D("h_offbeam", ";;", nbins, binlow, binhigh);
    TH1D* h_onbeam = new TH1D("h_onbeam", ";;", nbins, binlow, binhigh);

    for (int i = 0; i < t_bnbcos->GetEntries(); i++){

      t_bnbcos->GetEntry(i);
      CalcPIDvars(&bnbcos_vars);

        std::vector<double> minLogLikelihoods = {
         bnbcos_vars.track_neglogl_p,
         bnbcos_vars.track_neglogl_mu,
         bnbcos_vars.track_neglogl_pi,
         bnbcos_vars.track_neglogl_k,
         bnbcos_vars.track_neglogl_mip,
         bnbcos_vars.track_chi2mu,
         bnbcos_vars.track_chi2p,
         bnbcos_vars.track_chi2k,
         bnbcos_vars.track_chi2pi
        };


      if (TMath::Abs(bnbcos_vars.true_PDG) == 2212)
        h_bnbcos_p->Fill(minLogLikelihoods.at(j)*1./protonScaling);
      else if (TMath::Abs(bnbcos_vars.true_PDG) == 13)
        h_bnbcos_mu->Fill(minLogLikelihoods.at(j)*1./muonScaling);
      else if (TMath::Abs(bnbcos_vars.true_PDG) == 211)
        h_bnbcos_pi->Fill(minLogLikelihoods.at(j)*1./muonScaling);
      else if (TMath::Abs(bnbcos_vars.true_PDG) == 321)
        h_bnbcos_k->Fill(minLogLikelihoods.at(j)*1./muonScaling);
      else
        h_bnbcos_other->Fill(minLogLikelihoods.at(j)*1./muonScaling);

    }

    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars);

      std::vector<double> minLogLikelihoods = {
         onbeam_vars.track_neglogl_p,
         onbeam_vars.track_neglogl_mu,
         onbeam_vars.track_neglogl_pi,
         onbeam_vars.track_neglogl_k,
         onbeam_vars.track_neglogl_mip,
         onbeam_vars.track_chi2mu,
         onbeam_vars.track_chi2p,
         onbeam_vars.track_chi2k,
         onbeam_vars.track_chi2pi
        };

      h_onbeam->Fill(minLogLikelihoods.at(j));

    }
    for (int i = 0; i < t_offbeam->GetEntries(); i++){

      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars);

      std::vector<double> minLogLikelihoods = {
         offbeam_vars.track_neglogl_p,
         offbeam_vars.track_neglogl_mu,
         offbeam_vars.track_neglogl_pi,
         offbeam_vars.track_neglogl_k,
         offbeam_vars.track_neglogl_mip,
         offbeam_vars.track_chi2mu,
         offbeam_vars.track_chi2p,
         offbeam_vars.track_chi2k,
         offbeam_vars.track_chi2pi
        };

      h_offbeam->Fill(minLogLikelihoods.at(j));

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
    h_bnbcos_p->SetMarkerColor(TColor::GetColor(215, 48, 39));
    h_bnbcos_mu->SetMarkerColor(TColor::GetColor(8,64,129));
    h_bnbcos_pi->SetMarkerColor(TColor::GetColor(166,217,106));
    h_bnbcos_k->SetMarkerColor(TColor::GetColor(133,1,98));
    h_bnbcos_other->SetMarkerColor(TColor::GetColor(197,197,197));

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

    hs->SetTitle(";"+plotNames.at(j)+";");


       h_offbeam->Scale(offbeamScaling);
       TH1D* h_onminusoff = (TH1D*)h_onbeam->Clone("h_onminusoff");

       h_onminusoff->Add(h_offbeam, -1);
       h_onminusoff->Sumw2();
       h_onminusoff->Scale(1./h_onminusoff->Integral());

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

    c1->cd();
    hs->SetMaximum(std::max(h_onminusoff->GetMaximum(), h_total->GetMaximum())*1.15);
    hs->Draw("hist");
    h_total->SetFillStyle(3345);
    h_total->SetFillColor(kGray+2);
    h_total->SetMarkerSize(0.); // bad hack because root keeps drawing markers and I can't make it stop
    h_total->Draw("sameE2");
    h_onminusoff->SetMarkerStyle(20);
    h_onminusoff->SetMarkerSize(0.6);
    h_onminusoff->Draw("samepE1");

    TLegend *leg = new TLegend(0.55, 0.64, 0.89, 0.89);
    if (plotNames.at(j) == "track_neglogl_p"){
      std::cout << "gotya" << std::endl;
      leg->SetX1(0.15);
      leg->SetX2(0.45);
    }
    leg->AddEntry(h_bnbcos_p, "True proton");
    leg->AddEntry(h_bnbcos_mu, "True muon");
    leg->AddEntry(h_bnbcos_pi, "True pion");
    leg->AddEntry(h_bnbcos_k, "True kaon");
    leg->AddEntry(h_bnbcos_other, "True other");
    leg->AddEntry(h_onminusoff, "Data (on-off)");
    leg->SetLineWidth(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    c1->SaveAs(plotNames.at(j)+".png");

    h_bnbcos_p->Delete();
    h_bnbcos_mu->Delete();
    h_bnbcos_pi->Delete();
    h_bnbcos_k->Delete();
    h_bnbcos_other->Delete();
    h_offbeam->Delete();
    h_onbeam->Delete();

  }

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->cd();

  TH1D* h_bnbcos_pidamean_p  = new TH1D("h_bnbcos_pidamean_p", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidamean_mu = new TH1D("h_bnbcos_pidamean_mu", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidamean_pi = new TH1D("h_bnbcos_pidamean_pi", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidamean_k  = new TH1D("h_bnbcos_pidamean_k", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidamean_other = new TH1D("h_bnbcos_pidamean_other", ";PIDa;", 40, 0, 30);
  h_bnbcos_pidamean_p->SetFillColor(TColor::GetColor(215, 48, 39));
  h_bnbcos_pidamean_mu->SetFillColor(TColor::GetColor(8,64,129));
  h_bnbcos_pidamean_pi->SetFillColor(TColor::GetColor(166,217,106));
  h_bnbcos_pidamean_k->SetFillColor(TColor::GetColor(133,1,98));
  h_bnbcos_pidamean_other->SetFillColor(TColor::GetColor(197,197,197));
  h_bnbcos_pidamean_p->SetLineWidth(0);
  h_bnbcos_pidamean_mu->SetLineWidth(0);
  h_bnbcos_pidamean_pi->SetLineWidth(0);
  h_bnbcos_pidamean_k->SetLineWidth(0);
  h_bnbcos_pidamean_other->SetLineWidth(0);

  t_bnbcos->Draw("track_PIDA_mean >> h_bnbcos_pidamean_p", "TMath::Abs(true_PDG) == 2212");
  t_bnbcos->Draw("track_PIDA_mean >> h_bnbcos_pidamean_mu", "TMath::Abs(true_PDG) == 13");
  t_bnbcos->Draw("track_PIDA_mean >> h_bnbcos_pidamean_pi", "TMath::Abs(true_PDG) == 211");
  t_bnbcos->Draw("track_PIDA_mean >> h_bnbcos_pidamean_k", "TMath::Abs(true_PDG) == 321");
  t_bnbcos->Draw("track_PIDA_mean >> h_bnbcos_pidamean_other", "TMath::Abs(true_PDG) != 2212 && TMath::Abs(true_PDG) != 13 && TMath::Abs(true_PDG) !=211 && TMath::Abs(true_PDG) !=321");

  TH1D* h_bnbcos_pidamean_total = new TH1D("h_bnbcos_pidamean_total", ";PIDa (mean method);", 40, 0, 30);
  h_bnbcos_pidamean_total->Add(h_bnbcos_pidamean_p);
  h_bnbcos_pidamean_total->Add(h_bnbcos_pidamean_mu);
  h_bnbcos_pidamean_total->Add(h_bnbcos_pidamean_pi);
  h_bnbcos_pidamean_total->Add(h_bnbcos_pidamean_k);
  h_bnbcos_pidamean_total->Add(h_bnbcos_pidamean_other);

  h_bnbcos_pidamean_total->SetFillColor(kBlack);
  h_bnbcos_pidamean_total->SetFillStyle(3345);
  h_bnbcos_pidamean_total->SetMarkerSize(0.);

  h_bnbcos_pidamean_total->Sumw2();

  h_bnbcos_pidamean_p->Scale(1./h_bnbcos_pidamean_total->Integral());
  h_bnbcos_pidamean_mu->Scale(1./h_bnbcos_pidamean_total->Integral());
  h_bnbcos_pidamean_pi->Scale(1./h_bnbcos_pidamean_total->Integral());
  h_bnbcos_pidamean_k->Scale(1./h_bnbcos_pidamean_total->Integral());
  h_bnbcos_pidamean_other->Scale(1./h_bnbcos_pidamean_total->Integral());
  h_bnbcos_pidamean_total->Scale(1./h_bnbcos_pidamean_total->Integral());

  THStack *hs_bnbcos_pidamean = new THStack("hs_bnbcos_pidamean", "");
  hs_bnbcos_pidamean->Add(h_bnbcos_pidamean_p);
  hs_bnbcos_pidamean->Add(h_bnbcos_pidamean_mu);
  hs_bnbcos_pidamean->Add(h_bnbcos_pidamean_pi);
  hs_bnbcos_pidamean->Add(h_bnbcos_pidamean_k);
  hs_bnbcos_pidamean->Add(h_bnbcos_pidamean_other);

  TH1D* h_onbeam_pidamean = new TH1D("h_onbeam_pidamean", "", 40, 0, 30);
  t_onbeam->Draw("track_PIDA_mean >> h_onbeam_pidamean");

  TH1D* h_offbeam_pidamean = new TH1D("h_offbeam_pidamean", "", 40, 0, 30);
  t_offbeam->Draw("track_PIDA_mean >> h_offbeam_pidamean");

  h_offbeam_pidamean->Scale(offbeamScaling);

  TH1D* h_onbeamminusoffbeam = (TH1D*)h_onbeam_pidamean->Clone("h_onbeamminusoffbeam");
  h_onbeamminusoffbeam->Add(h_offbeam_pidamean, -1);

  h_bnbcos_pidamean_total->Draw("E2");
  //h_onbeamminusoffbeam->DrawNormalized("p same");
  hs_bnbcos_pidamean->Draw("same hist");
  h_onbeamminusoffbeam->SetMarkerStyle(20);
  h_onbeamminusoffbeam->SetMarkerSize(0.6);
  h_onbeamminusoffbeam->DrawNormalized("samepE1");
  c2->SaveAs("pidamean.png");

  TCanvas *c3 = new TCanvas("c2", "c2", 500, 500);
  c3->cd();

  TH1D* h_bnbcos_pidakde_p  = new TH1D("h_bnbcos_pidakde_p", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidakde_mu = new TH1D("h_bnbcos_pidakde_mu", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidakde_pi = new TH1D("h_bnbcos_pidakde_pi", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidakde_k  = new TH1D("h_bnbcos_pidakde_k", ";PIDa;", 40, 0, 30);
  TH1D* h_bnbcos_pidakde_other = new TH1D("h_bnbcos_pidakde_other", ";PIDa;", 40, 0, 30);
  h_bnbcos_pidakde_p->SetFillColor(TColor::GetColor(215, 48, 39));
  h_bnbcos_pidakde_mu->SetFillColor(TColor::GetColor(8,64,129));
  h_bnbcos_pidakde_pi->SetFillColor(TColor::GetColor(166,217,106));
  h_bnbcos_pidakde_k->SetFillColor(TColor::GetColor(133,1,98));
  h_bnbcos_pidakde_other->SetFillColor(TColor::GetColor(197,197,197));
  h_bnbcos_pidakde_p->SetLineWidth(0);
  h_bnbcos_pidakde_mu->SetLineWidth(0);
  h_bnbcos_pidakde_pi->SetLineWidth(0);
  h_bnbcos_pidakde_k->SetLineWidth(0);
  h_bnbcos_pidakde_other->SetLineWidth(0);


  t_bnbcos->Draw("track_PIDA_kde >> h_bnbcos_pidakde_p", "TMath::Abs(true_PDG) == 2212");
  t_bnbcos->Draw("track_PIDA_kde >> h_bnbcos_pidakde_mu", "TMath::Abs(true_PDG) == 13");
  t_bnbcos->Draw("track_PIDA_kde >> h_bnbcos_pidakde_pi", "TMath::Abs(true_PDG) == 211");
  t_bnbcos->Draw("track_PIDA_kde >> h_bnbcos_pidakde_k", "TMath::Abs(true_PDG) == 321");
  t_bnbcos->Draw("track_PIDA_kde >> h_bnbcos_pidakde_other", "TMath::Abs(true_PDG) != 2212 && TMath::Abs(true_PDG) != 13 && TMath::Abs(true_PDG) !=211 && TMath::Abs(true_PDG) !=321");

  TH1D* h_bnbcos_pidakde_total = new TH1D("h_bnbcos_pidakde_total", ";PIDa (KDE method);", 40, 0, 30);
  h_bnbcos_pidakde_total->Add(h_bnbcos_pidakde_p);
  h_bnbcos_pidakde_total->Add(h_bnbcos_pidakde_mu);
  h_bnbcos_pidakde_total->Add(h_bnbcos_pidakde_pi);
  h_bnbcos_pidakde_total->Add(h_bnbcos_pidakde_k);
  h_bnbcos_pidakde_total->Add(h_bnbcos_pidakde_other);

  h_bnbcos_pidakde_total->Sumw2();

  h_bnbcos_pidakde_p->Scale(1./h_bnbcos_pidakde_total->Integral());
  h_bnbcos_pidakde_mu->Scale(1./h_bnbcos_pidakde_total->Integral());
  h_bnbcos_pidakde_pi->Scale(1./h_bnbcos_pidakde_total->Integral());
  h_bnbcos_pidakde_k->Scale(1./h_bnbcos_pidakde_total->Integral());
  h_bnbcos_pidakde_other->Scale(1./h_bnbcos_pidakde_total->Integral());
  h_bnbcos_pidakde_total->Scale(1./h_bnbcos_pidakde_total->Integral());

  THStack *hs_bnbcos_pidakde = new THStack("hs_bnbcos_pidakde", "");
  hs_bnbcos_pidakde->Add(h_bnbcos_pidakde_p);
  hs_bnbcos_pidakde->Add(h_bnbcos_pidakde_mu);
  hs_bnbcos_pidakde->Add(h_bnbcos_pidakde_pi);
  hs_bnbcos_pidakde->Add(h_bnbcos_pidakde_k);
  hs_bnbcos_pidakde->Add(h_bnbcos_pidakde_other);

  TH1D* h_onbeam_pidakde = new TH1D("h_onbeam_pidakde", "", 40, 0, 30);
  t_onbeam->Draw("track_PIDA_kde >> h_onbeam_pidakde");

  TH1D* h_offbeam_pidakde = new TH1D("h_offbeam_pidakde", "", 40, 0, 30);
  t_offbeam->Draw("track_PIDA_kde >> h_offbeam_pidakde");

  h_offbeam_pidakde->Scale(offbeamScaling);

  TH1D* h_onbeamminusoffbeamkde = (TH1D*)h_onbeam_pidakde->Clone("h_onbeamminusoffbeamkde");
  h_onbeamminusoffbeamkde->Add(h_offbeam_pidakde, -1);

  h_bnbcos_pidakde_total->SetFillColor(kBlack);
  h_bnbcos_pidakde_total->SetFillStyle(3345);
  h_bnbcos_pidakde_total->SetMarkerSize(0.);

  h_bnbcos_pidakde_total->Draw("E2");
  //h_onbeamminusoffbeamkde->DrawNormalized("psame");
  hs_bnbcos_pidakde->Draw("same hist");
  h_onbeamminusoffbeamkde->SetMarkerStyle(20);
  h_onbeamminusoffbeamkde->SetMarkerSize(0.6);
  h_onbeamminusoffbeamkde->DrawNormalized("samepE1");
  c3->SaveAs("pidakde.png");

}
