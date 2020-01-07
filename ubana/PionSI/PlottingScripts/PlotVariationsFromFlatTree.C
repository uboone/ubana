

void PlotVariationsFromFlatTree(std::string infilename){

  TFile *infile = new TFile(infilename.c_str(),"open");
  TTree *intree = (TTree*)infile->Get("output");

  // Variables to read from tree
  int event;
  double weight;
  double track_length;
  int PDG;
  double init_momentum;
  double final_momentum;
  double final_nonzero_momentum;
  std::string *varied_par = nullptr;
  double varied_par_val;
  std::string *final_proc = nullptr;

  intree->SetBranchAddress("event", &event);
  intree->SetBranchAddress("weight", &weight);
  intree->SetBranchAddress("track_length", &track_length);
  intree->SetBranchAddress("PDG", &PDG);
  intree->SetBranchAddress("init_momentum", &init_momentum);
  intree->SetBranchAddress("final_momentum", &final_momentum);
  intree->SetBranchAddress("final_nonzero_momentum", &final_nonzero_momentum);
  intree->SetBranchAddress("varied_par", &varied_par);
  intree->SetBranchAddress("varied_par_val", &varied_par_val);
  intree->SetBranchAddress("final_proc", &final_proc);

  // Parameters we will vary
  // std::vector< std::string > all_pars = {"fReacLow","fReacHigh","fAbs","fCex","fDCex","fPiProd","fInel"};
  std::vector< std::string > all_pars = {"fReacHigh","fReacLow","fAbs","fCex","fDCex","fPiProd","fInel"};
  int n_pars = all_pars.size();
  std::vector<double> par_upval(n_pars,1);
  std::vector<double> par_downval(n_pars,1);

  // Make histograms
  TH1D *h_nom = new TH1D("h_nom","Nominal;Initial particle momentum (MeV/c);Weighted no. events",20,0,2000);
  // TH1D *h_nom = new TH1D("h_nom","Nominal;Track length (cm);Weighted no. events",60,0,300);
  TH1D *h_var[n_pars][2]; // [varied parameter][up/down]
  TH1D *h_var_ratio[n_pars][2]; // [varied parameter][up/down]
  TH1D *h_nom_finalP = new TH1D("h_nom_finalP","Nominal;Final particle momentum (MeV/c);Weighted no. events",50,0,2000);
  TH1D *h_var_finalP[n_pars][2]; // [varied parameter][up/down]
  TH1D *h_var_ratio_finalP[n_pars][2]; // [varied parameter][up/down]
  TH1D *h_nom_xsec = new TH1D("h_nom_xsec","Nominal;Final particle momentum (MeV/c);Inferred Cross Section (arb. units)",50,0,2000);
  TH1D *h_var_xsec[n_pars][2]; // [varied parameter][up/down]
  TH1D *h_var_ratio_xsec[n_pars][2]; // [varied parameter][up/down]
  for (int i_par=0; i_par<n_pars; i_par++){
    for (int i_updown=0; i_updown<2; i_updown++){
      std::string name = "h_var_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      // std::string title = all_pars.at(i_par) + " variations;Initial particle momentum;Weighted no. events";
      h_var[i_par][i_updown] = (TH1D*)h_nom->Clone(name.c_str());

      name = "h_var_finalP_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      // std::string title = all_pars.at(i_par) + " variations;Initial particle momentum;Weighted no. events";
      h_var_finalP[i_par][i_updown] = (TH1D*)h_nom_finalP->Clone(name.c_str());

      name = "h_var_xsec_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      // std::string title = all_pars.at(i_par) + " variations;Initial particle momentum;Weighted no. events";
      h_var_xsec[i_par][i_updown] = (TH1D*)h_nom_xsec->Clone(name.c_str());
    } // i_updown
  } // i_var

  // Loop through tree and fill histograms
  for (int i_entry=0; i_entry<intree->GetEntries(); i_entry++){
    // std::cout << i_entry << std::endl;
    intree->GetEntry(i_entry);

    // Look at pi+s only
    if (PDG != 211) continue;

    // Look at inelastic scatters only
    // if (final_proc->compare(std::string("pi+Inelastic"))!=0) continue;

    if (varied_par->compare(std::string("nominal"))==0){
      h_nom->Fill(init_momentum);
      // h_nom->Fill(track_length);
      // h_nom_finalP->Fill(final_nonzero_momentum);
      h_nom_finalP->Fill(final_momentum);
      if (final_proc->compare(std::string("pi+Inelastic"))==0){
        h_nom_xsec->Fill(final_momentum);
      }
    } // varied_par == "nominal"
    else{
      for (int i_par=0; i_par<n_pars; i_par++){
        if (varied_par->compare(all_pars.at(i_par))!=0) continue;

        int i_updown = -999;
        if (varied_par_val > 1) {
          i_updown = 0;
          par_upval.at(i_par) = varied_par_val;
        }
        else {
          i_updown = 1;
          par_downval.at(i_par) = varied_par_val;
        }

        // Skip if weight is too high (because cascade weights can get crazy)
        if (TMath::Abs(weight) > 5 || std::isnan(weight)) continue;
        h_var[i_par][i_updown]->Fill(init_momentum,weight);
        // h_var[i_par][i_updown]->Fill(track_length,weight);
        // h_var_finalP[i_par][i_updown]->Fill(final_nonzero_momentum,weight);
        h_var_finalP[i_par][i_updown]->Fill(final_momentum,weight);
        if (final_proc->compare(std::string("pi+Inelastic"))==0){
          h_var_xsec[i_par][i_updown]->Fill(final_momentum,weight);
        }


      } // loop over i_par
    } // varied_par != "nominal"
  } // i_entry

  // Calculate xsec
  h_nom_xsec->Divide(h_nom_finalP);
  TH1D *h_nom_xsec_copy = (TH1D*)h_nom_xsec->Clone("h_nom_xsec_copy");
  for (int j=1; j<h_nom_xsec->GetNbinsX()+1; j++) {
    float p1 = h_nom_xsec_copy->GetBinContent(j);
    float p2 = h_nom_xsec_copy->GetBinContent(j-1);
    float v = 0;
    if (p1 > p2 && p1 < 1) {
      v = -1.0 * log((1.0 - p1) / (1.0 - p2));
    }
    h_nom_xsec->SetBinContent(j, v);
  }

  for (int i_par=0; i_par<n_pars; i_par++){
    for (int i_updown=0; i_updown<2; i_updown++){
      h_var_xsec[i_par][i_updown]->Divide(h_var_finalP[i_par][i_updown]);
      TH1D *h_var_xsec_copy = (TH1D*)h_var_xsec[i_par][i_updown]->Clone("h_var_xsec_copy");
      for (int j=1; j<h_var_xsec[i_par][i_updown]->GetNbinsX()+1; j++) {
        float p1 = h_var_xsec_copy->GetBinContent(j);
        float p2 = h_var_xsec_copy->GetBinContent(j-1);
        float v = 0;
        if (p1 > p2 && p1 < 1) {
          v = -1.0 * log((1.0 - p1) / (1.0 - p2));
        }
        h_var_xsec[i_par][i_updown]->SetBinContent(j, v);
      }
    }
  }



  // Style histograms
  h_nom->SetLineWidth(2);
  h_nom->SetLineColor(kBlue+3);
  h_nom_finalP->SetLineWidth(2);
  h_nom_finalP->SetLineColor(kBlue+3);
  h_nom_xsec->SetLineWidth(2);
  h_nom_xsec->SetLineColor(kBlue+3);
  for (int i_updown = 0; i_updown<2; i_updown++){
    h_var[0][i_updown]->SetLineWidth(2); // fReacLow
    h_var[0][i_updown]->SetLineColor(kBlack);
    h_var[0][0]->SetLineStyle(9);
    h_var[0][1]->SetLineStyle(2);
    h_var[1][i_updown]->SetLineWidth(2); // fReacHigh
    h_var[1][i_updown]->SetLineColor(kBlack);
    h_var[1][0]->SetLineStyle(9);
    h_var[1][1]->SetLineStyle(2);
    h_var[2][i_updown]->SetLineWidth(2); // fAbs
    h_var[2][i_updown]->SetLineColor(kGray+1);
    h_var[2][0]->SetLineStyle(9);
    h_var[2][1]->SetLineStyle(2);
    h_var[3][i_updown]->SetLineWidth(2); // fCex
    h_var[3][i_updown]->SetLineColor(kBlue);
    h_var[3][0]->SetLineStyle(9);
    h_var[3][1]->SetLineStyle(2);
    h_var[4][i_updown]->SetLineWidth(2); // fDCex
    h_var[4][i_updown]->SetLineColor(kGreen+2);
    h_var[4][0]->SetLineStyle(9);
    h_var[4][1]->SetLineStyle(2);
    h_var[5][i_updown]->SetLineWidth(2); // fPiProd
    h_var[5][i_updown]->SetLineColor(kOrange+2);
    h_var[5][0]->SetLineStyle(9);
    h_var[5][1]->SetLineStyle(2);
    h_var[6][i_updown]->SetLineWidth(2); // fInel
    h_var[6][i_updown]->SetLineColor(kRed);
    h_var[6][0]->SetLineStyle(9);
    h_var[6][1]->SetLineStyle(2);

    h_var_finalP[0][i_updown]->SetLineWidth(2); // fReacLow
    h_var_finalP[0][i_updown]->SetLineColor(kBlack);
    h_var_finalP[0][0]->SetLineStyle(9);
    h_var_finalP[0][1]->SetLineStyle(2);
    h_var_finalP[1][i_updown]->SetLineWidth(2); // fReacHigh
    h_var_finalP[1][i_updown]->SetLineColor(kBlack);
    h_var_finalP[1][0]->SetLineStyle(9);
    h_var_finalP[1][1]->SetLineStyle(2);
    h_var_finalP[2][i_updown]->SetLineWidth(2); // fAbs
    h_var_finalP[2][i_updown]->SetLineColor(kGray+1);
    h_var_finalP[2][0]->SetLineStyle(9);
    h_var_finalP[2][1]->SetLineStyle(2);
    h_var_finalP[3][i_updown]->SetLineWidth(2); // fCex
    h_var_finalP[3][i_updown]->SetLineColor(kBlue);
    h_var_finalP[3][0]->SetLineStyle(9);
    h_var_finalP[3][1]->SetLineStyle(2);
    h_var_finalP[4][i_updown]->SetLineWidth(2); // fDCex
    h_var_finalP[4][i_updown]->SetLineColor(kGreen+2);
    h_var_finalP[4][0]->SetLineStyle(9);
    h_var_finalP[4][1]->SetLineStyle(2);
    h_var_finalP[5][i_updown]->SetLineWidth(2); // fPiProd
    h_var_finalP[5][i_updown]->SetLineColor(kOrange+2);
    h_var_finalP[5][0]->SetLineStyle(9);
    h_var_finalP[5][1]->SetLineStyle(2);
    h_var_finalP[6][i_updown]->SetLineWidth(2); // fInel
    h_var_finalP[6][i_updown]->SetLineColor(kRed);
    h_var_finalP[6][0]->SetLineStyle(9);
    h_var_finalP[6][1]->SetLineStyle(2);

    h_var_xsec[0][i_updown]->SetLineWidth(2); // fReacLow
    h_var_xsec[0][i_updown]->SetLineColor(kBlack);
    h_var_xsec[0][0]->SetLineStyle(9);
    h_var_xsec[0][1]->SetLineStyle(2);
    h_var_xsec[1][i_updown]->SetLineWidth(2); // fReacHigh
    h_var_xsec[1][i_updown]->SetLineColor(kBlack);
    h_var_xsec[1][0]->SetLineStyle(9);
    h_var_xsec[1][1]->SetLineStyle(2);
    h_var_xsec[2][i_updown]->SetLineWidth(2); // fAbs
    h_var_xsec[2][i_updown]->SetLineColor(kGray+1);
    h_var_xsec[2][0]->SetLineStyle(9);
    h_var_xsec[2][1]->SetLineStyle(2);
    h_var_xsec[3][i_updown]->SetLineWidth(2); // fCex
    h_var_xsec[3][i_updown]->SetLineColor(kBlue);
    h_var_xsec[3][0]->SetLineStyle(9);
    h_var_xsec[3][1]->SetLineStyle(2);
    h_var_xsec[4][i_updown]->SetLineWidth(2); // fDCex
    h_var_xsec[4][i_updown]->SetLineColor(kGreen+2);
    h_var_xsec[4][0]->SetLineStyle(9);
    h_var_xsec[4][1]->SetLineStyle(2);
    h_var_xsec[5][i_updown]->SetLineWidth(2); // fPiProd
    h_var_xsec[5][i_updown]->SetLineColor(kOrange+2);
    h_var_xsec[5][0]->SetLineStyle(9);
    h_var_xsec[5][1]->SetLineStyle(2);
    h_var_xsec[6][i_updown]->SetLineWidth(2); // fInel
    h_var_xsec[6][i_updown]->SetLineColor(kRed);
    h_var_xsec[6][0]->SetLineStyle(9);
    h_var_xsec[6][1]->SetLineStyle(2);
  }

  // Make ratio histograms
  for (int i_par=0; i_par<n_pars; i_par++){
    for (int i_updown=0; i_updown<2; i_updown++){
      std::string name = "h_var_ratio_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      h_var_ratio[i_par][i_updown] = (TH1D*)h_var[i_par][i_updown]->Clone(name.c_str());
      // h_var_ratio[i_par][i_updown]->Add(h_nom,-1);
      h_var_ratio[i_par][i_updown]->Divide(h_nom);

      name = "h_var_ratio_finalP_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      h_var_ratio_finalP[i_par][i_updown] = (TH1D*)h_var_finalP[i_par][i_updown]->Clone(name.c_str());
      // h_var_ratio_finalP[i_par][i_updown]->Add(h_nom_finalP,-1);
      h_var_ratio_finalP[i_par][i_updown]->Divide(h_nom_finalP);

      name = "h_var_ratio_xsec_" + std::to_string(i_par) + "_" + std::to_string(i_updown);
      h_var_ratio_xsec[i_par][i_updown] = (TH1D*)h_var_xsec[i_par][i_updown]->Clone(name.c_str());
      // h_var_ratio_xsec[i_par][i_updown]->Add(h_nom_xsec,-1);
      h_var_ratio_xsec[i_par][i_updown]->Divide(h_nom_xsec);
    } // i_updown
  } // i_var

  // Make legend
  TH1D *dummy_nom = new TH1D();
  dummy_nom->SetLineColor(kBlue+3);
  dummy_nom->SetLineWidth(2);
  TH1D *dummy_upvar = new TH1D();
  dummy_upvar->SetLineColor(kBlack);
  dummy_upvar->SetLineWidth(2);
  dummy_upvar->SetLineStyle(9);
  TH1D *dummy_downvar = new TH1D();
  dummy_downvar->SetLineColor(kBlack);
  dummy_downvar->SetLineWidth(2);
  dummy_downvar->SetLineStyle(2);

  TH1D *h_nom_ratio = (TH1D*)h_nom->Clone("h_nom_ratio");
  h_nom_ratio->GetYaxis()->SetTitle("(Var-Nom)/Nom");
  h_nom_ratio->Divide(h_nom);
  TH1D *h_nom_ratio_finalP = (TH1D*)h_nom_finalP->Clone("h_nom_ratio_finalP");
  h_nom_ratio_finalP->GetYaxis()->SetTitle("(Var-Nom)/Nom");
  h_nom_ratio_finalP->Divide(h_nom_finalP);
  TH1D *h_nom_ratio_xsec = (TH1D*)h_nom_xsec->Clone("h_nom_ratio_xsec");
  h_nom_ratio_xsec->GetYaxis()->SetTitle("(Var-Nom)/Nom");
  h_nom_ratio_xsec->Divide(h_nom_xsec);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas();//"","",400,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.2);
  c1->cd();
  pad1->Draw();
  pad2->Draw();

  for (int i_par=0; i_par<n_pars; i_par++){

    // std::cout << h_nom->Integral() << std::endl;
    // std::cout << h_var[i_par][0]->Integral() << std::endl;
    // std::cout << h_var[i_par][1]->Integral() << std::endl << std::endl;
    c1->cd();
    pad1->cd();
    h_nom->SetLabelOffset(5);
    h_nom->SetTitleOffset(5);
    h_nom->Draw();
    h_var[i_par][1]->Draw("same");
    h_var[i_par][0]->Draw("same");

    TLegend *leg = new TLegend(0.5,0.6,0.88,0.88);
    leg->SetFillStyle(0);
    leg->SetLineColor(kWhite);
    std::string label = all_pars.at(i_par) + " = " + std::to_string(par_upval.at(i_par));
    leg->AddEntry(h_var[i_par][0],label.c_str(),"l");
    label = all_pars.at(i_par) + " = " + std::to_string(par_downval.at(i_par));
    leg->AddEntry(h_var[i_par][1],label.c_str(),"l");
    leg->AddEntry(h_nom,"Nominal","l");
    leg->Draw();

    pad2->cd();
    h_nom_ratio->Draw();
    h_nom_ratio->GetXaxis()->SetTitle("Initial particle momentum (MeV/c)");
    // h_nom_ratio->GetXaxis()->SetTitle("Track length (cm)");
    // h_nom_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
    h_nom_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    h_nom_ratio->GetYaxis()->SetNdivisions(010);
    h_nom_ratio->GetYaxis()->SetLabelSize(0.1);
    h_nom_ratio->GetYaxis()->SetTitleSize(0.1);
    h_nom_ratio->GetYaxis()->SetTitleOffset(0.3);
    h_nom_ratio->GetXaxis()->SetLabelSize(0.1);
    h_nom_ratio->GetXaxis()->SetTitleSize(0.1);
    h_nom_ratio->GetXaxis()->SetTitleOffset(0.7);
    h_var_ratio[i_par][1]->Draw("same");
    h_var_ratio[i_par][0]->Draw("same");

    std::string printname = "PlotReweight_" + all_pars.at(i_par) + ".pdf";
    c1->Print(printname.c_str());



    // std::cout << h_nom_finalP->Integral() << std::endl;
    // std::cout << h_var_finalP[i_par][0]->Integral() << std::endl;
    // std::cout << h_var_finalP[i_par][1]->Integral() << std::endl << std::endl;
    c1->cd();
    pad1->cd();
    h_nom_finalP->SetLabelOffset(5);
    h_nom_finalP->SetTitleOffset(5);
    h_nom_finalP->Draw();
    h_var_finalP[i_par][1]->Draw("same");
    h_var_finalP[i_par][0]->Draw("same");
    leg->Draw();

    pad2->cd();
    // h_nom_ratio->GetXaxis()->SetTitle("Final (non-zero) particle momentum (MeV/c)");
    h_nom_ratio_finalP->GetXaxis()->SetTitle("Final particle momentum (MeV/c)");
    h_nom_ratio_finalP->GetYaxis()->SetRangeUser(0.5,1.5);
    h_nom_ratio_finalP->GetYaxis()->SetNdivisions(010);
    h_nom_ratio_finalP->GetYaxis()->SetLabelSize(0.1);
    h_nom_ratio_finalP->GetYaxis()->SetTitleSize(0.1);
    h_nom_ratio_finalP->GetYaxis()->SetTitleOffset(0.3);
    h_nom_ratio_finalP->GetXaxis()->SetLabelSize(0.1);
    h_nom_ratio_finalP->GetXaxis()->SetTitleSize(0.1);
    h_nom_ratio_finalP->GetXaxis()->SetTitleOffset(0.7);
    h_nom_ratio_finalP->Draw();
    h_var_ratio_finalP[i_par][1]->Draw("same");
    h_var_ratio_finalP[i_par][0]->Draw("same");

    printname = "PlotReweight_finalP_" + all_pars.at(i_par) + ".pdf";
    c1->Print(printname.c_str());



    // std::cout << h_nom_xsec->Integral() << std::endl;
    // std::cout << h_var_xsec[i_par][0]->Integral() << std::endl;
    // std::cout << h_var_xsec[i_par][1]->Integral() << std::endl << std::endl;
    c1->cd();
    pad1->cd();
    h_nom_xsec->SetLabelOffset(5);
    h_nom_xsec->SetTitleOffset(5);
    h_nom_xsec->Draw();
    h_var_xsec[i_par][1]->Draw("same");
    h_var_xsec[i_par][0]->Draw("same");
    leg->Draw();

    pad2->cd();
    // h_nom_ratio->GetXaxis()->SetTitle("Final (non-zero) particle momentum (MeV/c)");
    h_nom_ratio_xsec->GetXaxis()->SetTitle("Final particle momentum (MeV/c)");
    h_nom_ratio_xsec->GetYaxis()->SetRangeUser(0.5,1.5);
    h_nom_ratio_xsec->GetYaxis()->SetNdivisions(010);
    h_nom_ratio_xsec->GetYaxis()->SetLabelSize(0.1);
    h_nom_ratio_xsec->GetYaxis()->SetTitleSize(0.1);
    h_nom_ratio_xsec->GetYaxis()->SetTitleOffset(0.3);
    h_nom_ratio_xsec->GetXaxis()->SetLabelSize(0.1);
    h_nom_ratio_xsec->GetXaxis()->SetTitleSize(0.1);
    h_nom_ratio_xsec->GetXaxis()->SetTitleOffset(0.7);
    h_nom_ratio_xsec->Draw();
    h_var_ratio_xsec[i_par][1]->Draw("same");
    h_var_ratio_xsec[i_par][0]->Draw("same");

    printname = "PlotReweight_xsec_" + all_pars.at(i_par) + ".pdf";
    c1->Print(printname.c_str());

  } // i_par

}
