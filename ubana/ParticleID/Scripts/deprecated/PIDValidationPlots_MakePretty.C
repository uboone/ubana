// Macro to make summary/overlaid plots from the output of the
// ParticleIDValidationPlots module. For now it only makes comparisons
// between true muons and true protons (because that's what we care about most)
// but true pion and kaon information is also available so those plots could be made.
//
// Run using:
// ~ root -b
// root[0] .L PIDValidationPlots_MakePretty.C
// root[1] PIDValidationPlots_MakePretty("particleIDMeta.root")
//
// The output is 1) a set of .png plots produced in the current directory, and
// 2) a root file PIDValidationPlots_out.root containing the same plots, also
// produced in the current directory. The percentage of particles that are
// correctly identified (by which we mean: for which the lowest neg2LL is the
// one corresponding to the correct particle type) is also printed to the screen.
//
// Kirsty Duffy, Fermilab, 11th April 2018

void PIDValidationPlots_MakePretty(std::string inputfile){

  std::vector<std::string> compare_1dplots = { "neglogl_mu",
                                            "neglogl_p",
                                            "neglogl_pi",
                                            "neglogl_K",
                                            "neglogl_MIP",
                                            "neglogl_minmuMIP",
                                            "neglogl_muoverp",
                                            "neglogl_muminusp",
                                            "neglogl_MIPminusp",
                                            "neglogl_minmuMIPminusp",
                                            // "pull_mu",
                                            // "pull_p",
                                            // "pull_pi",
                                            // "pull_K",
                                            // "pull_MIP",
                                            "PIDA" };

  std::vector<std::string> compare_2dplots = { "neglogl_muvsp",
                                              "neglogl_muvspi",
                                              "dEdxtr_len",
                                              "neglogl_MIPvsp",
                                              "neglogl_minmuMIPvsp",
                                              "neglogl_mu_vslength",
                                              "neglogl_MIP_vslength",
                                              "neglogl_minmuMIP_vslength",
                                              "neglogl_p_vslength",
                                              "neglogl_mu_vsangle",
                                              "neglogl_MIP_vsangle",
                                              "neglogl_minmuMIP_vsangle",
                                              "neglogl_p_vsangle",
                                              "neglogl_mu_vsnhits",
                                              "neglogl_MIP_vsnhits",
                                              "neglogl_minmuMIP_vsnhits",
                                              "neglogl_p_vsnhits", };

  std::vector<std::string> categories = {"TrueBragg", "All"};

  // --------------------------------------------------------------------- //
  // Now the actual plotting stuff starts!
  // --------------------------------------------------------------------- //

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);

  TFile *fin = new TFile(inputfile.c_str(),"open");
  TFile *fout = new TFile("PIDValidationPlots_out.root","recreate");

  // --------------------------------------------------------------------- //
  // Comparison plots: 1D
  // --------------------------------------------------------------------- //

  // 1D comparison plots: look through all names in the compare_1dplots vector and overlay true muons with true protons
  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah
  for (int i_name=0; i_name<compare_1dplots.size(); i_name++){
    for (int i_cat=0; i_cat < categories.size(); i_cat++){

      std::string savename = categories.at(i_cat)+"_"+compare_1dplots.at(i_name);
      std::cout << "Making 1D overlay plots for " << savename << std::endl;

      std::string name_mu = "pidvalid/"+categories.at(i_cat)+"_truemu_"+compare_1dplots.at(i_name);
      std::string name_p = "pidvalid/"+categories.at(i_cat)+"_truep_"+compare_1dplots.at(i_name);
      std::string name_pi = "pidvalid/"+categories.at(i_cat)+"_truepi_"+compare_1dplots.at(i_name);
      std::string name_K = "pidvalid/"+categories.at(i_cat)+"_trueK_"+compare_1dplots.at(i_name);

      // Get plots
      fin->cd();
      TH1F *hmu = (TH1F*)fin->Get(name_mu.c_str());
      TH1F *hp = (TH1F*)fin->Get(name_p.c_str());
      TH1F *hpi = (TH1F*)fin->Get(name_pi.c_str());
      TH1F *hK = (TH1F*)fin->Get(name_K.c_str());

      // Style
      hmu->GetYaxis()->SetRangeUser(0,std::max(hmu->GetMaximum(),hp->GetMaximum())*1.15);
      hmu->SetLineColor(kBlue);
      hmu->SetFillColor(kBlue);
      hmu->SetFillStyle(3144);
      hp->SetLineColor(kRed);
      hp->SetFillColor(kRed);
      hp->SetFillStyle(3144);
      if (hpi) hpi->SetLineColor(kGreen+2);
      if (hpi) hpi->SetFillColor(kGreen+2);
      if (hpi) hpi->SetFillStyle(3144);
      if (hK) hK->SetLineColor(kViolet);
      if (hK) hK->SetFillColor(kViolet);
      if (hK) hK->SetFillStyle(3144);

      // Legend
      TLegend *l = new TLegend(0.6,0.6,0.87,0.87);
      l->SetLineColor(kWhite);
      l->SetFillStyle(0);
      l->SetTextFont(132);
      l->AddEntry(hmu,"True Muons","lf");
      l->AddEntry(hp,"True Protons","lf");
      if (hpi) l->AddEntry(hpi,"True Pions","lf");
      if (hK) l->AddEntry(hK,"True Kaons","lf");

      // Draw
      TCanvas *c1 = new TCanvas();
      hmu->Draw();
      hp->Draw("same");
      if (hpi) hpi->Draw("same");
      if (hK)  hK->Draw("same");
      l->Draw();

      // Save
      fout->cd();
      c1->Write(savename.c_str());
      c1->Print(std::string(savename+".png").c_str());

    } // loop over i_cat in categories
  } // loop over i_name in compare_1dplots



  // --------------------------------------------------------------------- //
  // Comparison plots: 2D
  // --------------------------------------------------------------------- //

  // 2D comparison plots: look through all names in the compare_2dplots vector and overlay true muons with true protons
  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah
  for (int i_name=0; i_name<compare_2dplots.size(); i_name++){
    for (int i_cat=0; i_cat < categories.size(); i_cat++){

      std::string savename = categories.at(i_cat)+"_"+compare_2dplots.at(i_name);
      std::cout << "Making 2D overlay plots for " << savename << std::endl;

      std::string name_mu = "pidvalid/"+categories.at(i_cat)+"_truemu_"+compare_2dplots.at(i_name);
      std::string name_p = "pidvalid/"+categories.at(i_cat)+"_truep_"+compare_2dplots.at(i_name);
      std::string name_pi = "pidvalid/"+categories.at(i_cat)+"_truepi_"+compare_2dplots.at(i_name);
      std::string name_K = "pidvalid/"+categories.at(i_cat)+"_trueK_"+compare_2dplots.at(i_name);

      // Get plots
      fin->cd();
      TH2F *hmu = (TH2F*)fin->Get(name_mu.c_str());
      TH2F *hp = (TH2F*)fin->Get(name_p.c_str());
      TH1F *hpi = (TH1F*)fin->Get(name_pi.c_str());
      TH1F *hK = (TH1F*)fin->Get(name_K.c_str());

      // Style
      //hmu->GetYaxis()->SetRangeUser(0,std::max(hmu->GetMaximum(),hp->GetMaximum())*1.15);
      hmu->SetMarkerColor(kBlue);
      hmu->SetMarkerStyle(7);
      hp->SetMarkerColor(kRed);
      hp->SetMarkerStyle(7);
      if (hpi) hpi->SetMarkerColor(kGreen+2);
      if (hpi) hpi->SetMarkerStyle(7);
      if (hK) hK->SetMarkerColor(kViolet);
      if (hK) hK->SetMarkerStyle(7);

      // Legend
      TLegend *l = new TLegend(0.6,0.6,0.87,0.87);
      l->SetLineColor(kWhite);
      l->SetFillStyle(0);
      l->SetTextFont(132);
      l->AddEntry(hmu,"True Muons","p");
      l->AddEntry(hp,"True Protons","p");
      if (hpi) l->AddEntry(hpi,"True Pions","lf");
      if (hK) l->AddEntry(hK,"True Kaons","lf");

      // Draw
      TCanvas *c1 = new TCanvas();
      hmu->Draw();
      hp->Draw("same");
      if (hpi) hpi->Draw("same");
      if (hK)  hK->Draw("same");
      l->Draw();

      // Save
      fout->cd();
      c1->Write(savename.c_str());
      c1->Print(std::string(savename+".png").c_str());

    } // loop over i_cat in categories
  } // loop over i_name in compare_2dplots



  // --------------------------------------------------------------------- //
  // Ratio plots to see how often we get track direction right
  // --------------------------------------------------------------------- //

  // 2D comparison plots: look through all names in the compare_2dplots vector and overlay true muons with true protons
  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah

  // Get plots
  fin->cd();
  TH2F *hTrueBragg_correct = (TH2F*)fin->Get("pidvalid/TrueBragg_correctdirection");
  TH2F *hTrueBragg_incorrect = (TH2F*)fin->Get("pidvalid/TrueBragg_incorrectdirection");
  TH2F *hAll_correct = (TH2F*)fin->Get("pidvalid/All_correctdirection");
  TH2F *hAll_incorrect = (TH2F*)fin->Get("pidvalid/All_incorrectdirection");

  // Make ratio plots

  TH1F *hTrueBragg_correct_ratio = new TH1F("hTrueBragg_correct_ratio","Tracks with true p=0 at end, correct reconstructed direction;PID direction: correct/incorrect;True particle",4,0,4);
  TH1F *hTrueBragg_incorrect_ratio = new TH1F("hTrueBragg_incorrect_ratio","Tracks with true p=0 at end, wrong reconstructed direction;PID direction: correct/incorrect;True particle",4,0,4);
  TH1F *hAll_correct_ratio = new TH1F("hAll_correct_ratio","All tracks, correct reconstructed direction;PID direction: correct/incorrect;True particle",4,0,4);
  TH1F *hAll_incorrect_ratio = new TH1F("hAll_incorrect_ratio","All tracks, wrong reconstructed direction;PID direction: correct/incorrect;True particle",4,0,4);

  const char* particles[4] = {"#mu", "p", "#pi", "K"};
  for (size_t i=0; i<4; i++){
    hTrueBragg_correct_ratio->GetXaxis()->SetBinLabel(i+1,particles[i]);
    hTrueBragg_incorrect_ratio->GetXaxis()->SetBinLabel(i+1,particles[i]);
    hAll_correct_ratio->GetXaxis()->SetBinLabel(i+1,particles[i]);
    hAll_incorrect_ratio->GetXaxis()->SetBinLabel(i+1,particles[i]);

    if (hTrueBragg_correct && hTrueBragg_incorrect){
      hTrueBragg_correct_ratio->SetBinContent(i+1,hTrueBragg_correct->GetBinContent(i+1,2)/hTrueBragg_correct->GetBinContent(i+1,1));
      hTrueBragg_incorrect_ratio->SetBinContent(i+1,hTrueBragg_incorrect->GetBinContent(i+1,2)/hTrueBragg_incorrect->GetBinContent(i+1,1));
    }
    if (hAll_correct && hAll_incorrect){
      hAll_correct_ratio->SetBinContent(i+1,hAll_correct->GetBinContent(i+1,2)/hAll_correct->GetBinContent(i+1,1));
      hAll_incorrect_ratio->SetBinContent(i+1,hAll_incorrect->GetBinContent(i+1,2)/hAll_incorrect->GetBinContent(i+1,1));
    }
  }

  // Draw
  TCanvas *c1 = new TCanvas();
  hTrueBragg_correct_ratio->Draw();
  // Save
  fout->cd();
  c1->Write("hTrueBragg_correct_ratio");
  c1->Print("hTrueBragg_correct_ratio.png");

  // Draw
  hTrueBragg_incorrect_ratio->Draw();
  // Save
  fout->cd();
  c1->Write("hTrueBragg_incorrect_ratio");
  c1->Print("hTrueBragg_incorrect_ratio.png");

  hAll_correct_ratio->Draw();
  // Save
  fout->cd();
  c1->Write("hAll_correct_ratio");
  c1->Print("hAll_correct_ratio.png");

  // Draw
  hAll_incorrect_ratio->Draw();
  // Save
  fout->cd();
  c1->Write("hAll_incorrect_ratio");
  c1->Print("hAll_incorrect_ratio.png");


  // --------------------------------------------------------------------- //
  // Print out how often we got the PID right
  // --------------------------------------------------------------------- //

  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah

  // TrueBragg
  std::cout << " --------- Tracks for which MCParticles have true p=0 at the end (true Bragg peaks) --------- " << std::endl;

  fin->cd();
  TH1F *h = (TH1F*)fin->Get("pidvalid/TrueBragg_truemu_smallest_neglogl");
  std::cout << "--- True muons: " << (h->GetBinContent(1)+h->GetBinContent(5))/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_truep_smallest_neglogl");
  std::cout << "--- True protons: " << h->GetBinContent(2)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_truepi_smallest_neglogl");
  std::cout << "--- True pions: " << (h->GetBinContent(3)+h->GetBinContent(1)+h->GetBinContent(5))/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_trueK_smallest_neglogl");
  std::cout << "--- True kaons: " << h->GetBinContent(4)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;


  // All particles
  std::cout << " --------- All tracks --------- " << std::endl;

  fin->cd();
  h = (TH1F*)fin->Get("pidvalid/All_truemu_smallest_neglogl");
  std::cout << "--- True muons: " << (h->GetBinContent(1)+h->GetBinContent(5))/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_truep_smallest_neglogl");
  std::cout << "--- True protons: " << h->GetBinContent(2)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_truepi_smallest_neglogl");
  std::cout << "--- True pions: " << h->GetBinContent(3)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_trueK_smallest_neglogl");
  std::cout << "--- True kaons: " << h->GetBinContent(4)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as MIP (no Bragg peak): " << h->GetBinContent(5)/h->Integral()*100.0 << "%" << std::endl  << std::endl;


}
