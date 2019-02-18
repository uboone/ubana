//
// KEY: TH1F      AllTracks_neglogl_mu;1
//  KEY: TH1F      AllTracks_neglogl_p;1
//   KEY: TH1F      AllTracks_neglogl_pi;1
//    KEY: TH1F      AllTracks_neglogl_K;1
//     KEY: TH1F      AllTracks_neglogl_MIP;1
//      KEY: TH1F      AllTracks_neglogl_minmuMIP;1
//       KEY: TH1F      AllTracks_neglogl_muoverp;1
//        KEY: TH1F      AllTracks_neglogl_muminusp;1
//         KEY: TH1F      AllTracks_neglogl_MIPminusp;1
//          KEY: TH1F      AllTracks_neglogl_minmuMIPminusp;1
//           KEY: TH2F      AllTracks_neglogl_mu_vslength;1
//            KEY: TH2F      AllTracks_neglogl_p_vslength;1
//             KEY: TH2F      AllTracks_neglogl_pi_vslength;1
//              KEY: TH2F      AllTracks_neglogl_K_vslength;1
//               KEY: TH2F      AllTracks_neglogl_MIP_vslength;1
//                KEY: TH2F      AllTracks_neglogl_minmuMIP_vslength;1
//                 KEY: TH2F      AllTracks_neglogl_mu_vsangle;1
//                  KEY: TH2F      AllTracks_neglogl_p_vsangle;1
//                   KEY: TH2F      AllTracks_neglogl_pi_vsangle;1
//                    KEY: TH2F      AllTracks_neglogl_K_vsangle;1
//                     KEY: TH2F      AllTracks_neglogl_MIP_vsangle;1
//                      KEY: TH2F      AllTracks_neglogl_minmuMIP_vsangle;1
//                       KEY: TH2F      AllTracks_neglogl_mu_vsnhits;1
//                        KEY: TH2F      AllTracks_neglogl_p_vsnhits;1
//                         KEY: TH2F      AllTracks_neglogl_pi_vsnhits;1
//                          KEY: TH2F      AllTracks_neglogl_K_vsnhits;1
//                           KEY: TH2F      AllTracks_neglogl_MIP_vsnhits;1
//                            KEY: TH2F      AllTracks_neglogl_minmuMIP_vsnhits;1
//                             KEY: TH2F      AllTracks_neglogl_muvsp;1
//                              KEY: TH2F      AllTracks_neglogl_muvspi;1
//                               KEY: TH2F      AllTracks_neglogl_MIPvsp;1
//                                KEY: TH2F      AllTracks_neglogl_minmuMIPvsp;1
//                                 KEY: TH2F      AllTracks_dEdxtr_len;1
//
//
//


std::vector<EColor> cols{
  kAzure,
    kRed,
    kGreen,
    kViolet
};

void styleHistogram(TH1D* h, int style){

  h->SetLineWidth(0);
  h->SetFillColor(cols.at(style));

}

void PIDValidationPlots_DataMCComparison_MakePretty(){

  double honPOT = 0.714;
  double hoffPOT= 1;

  /**
   * Input files:
   * inputMC: Simulated data PID root file
   * inputOnBD: "beam-on" physics data
   * inputOffBD: "beam-off" physics data
   */

  TFile *inputMC = 
    new TFile("/uboone/data/users/alister1/particleID/180417-ParticleId/pid_BNBCOS.root", "READ");
  TFile *inputOnBD = 
    new TFile("/uboone/data/users/alister1/particleID/180417-ParticleId/pid_OnBeam.root", "READ");
  TFile *inputOffBD = 
    new TFile("/uboone/data/users/alister1/particleID/180417-ParticleId/pid_OffBeam.root", "READ");

  /** deal with data first */
  std::vector<TString> BDPlotNames{
    "neglogl_mu",
      "neglogl_p",
      "neglogl_pi",
      "neglogl_K",
      "neglogl_MIP",
      "neglogl_minmuMIP",
      "neglogl_muoverp",
      "neglogl_muminusp"
      //"neglogl_MIPminusp",
      //"neglogl_minmuMIPminusp"
  };



  TString prefix("pidvalid/");

  TH1D* hon;
  TH1D* hoff;
  TH1D* honminusoff;
  for (int i = 0; i < BDPlotNames.size(); i++){

    std::cout << "[PLOTMAKER] Getting data plots" << std::endl;

    hon  = (TH1D*)inputOnBD->Get(prefix+"AllTracks_"+BDPlotNames.at(i));
    hoff =  (TH1D*)inputOffBD->Get(prefix+"AllTracks_"+BDPlotNames.at(i));

    hon->RebinX(4);
    hoff->RebinX(4);

    /** rescale histograms and produce onb-offb*/

    std::cout << "hon: " << hon->Integral() << std::endl;
    std::cout << "hoff: " << hoff->Integral() << std::endl;

    std::cout << "hon(2) " << hon->GetBinContent(1) << std::endl;
    std::cout << "hoff(2) " << hoff->GetBinContent(1) << std::endl;

    hon->Sumw2();
    hoff->Sumw2();
    hoff->Scale(honPOT/hoffPOT);

    std::cout << "hoff, post scaling: " << hoff->Integral() << std::endl;
    std::cout << "hoff, post scaling(2): " << hoff->GetBinContent(1) << std::endl;

    honminusoff = (TH1D*)hon->Clone("honminusoff");
    honminusoff->Add(hoff, -1);
    honminusoff->SetMarkerStyle(20);
    honminusoff->SetMarkerSize(0.7);

    std::cout << "honminusoff: " << honminusoff->Integral() << std::endl;
    std::cout << "honminusoff(2) " << honminusoff->GetBinContent(1) << std::endl;

    std::cout << "[PLOTMAKER] Getting simulation plots" << std::endl;

    /** get relevant MC histogramsi and style them */
    TH1D* truemu = (TH1D*)inputMC->Get(prefix+"All_truemu_"+BDPlotNames.at(i));
    TH1D* truep  = (TH1D*)inputMC->Get(prefix+"All_truep_"+BDPlotNames.at(i));
    TH1D* truepi = (TH1D*)inputMC->Get(prefix+"All_truepi_"+BDPlotNames.at(i));
    TH1D* truek  = (TH1D*)inputMC->Get(prefix+"All_trueK_"+BDPlotNames.at(i));
    TH1D* hall = (TH1D*)inputMC->Get(prefix+"AllTracks_"+BDPlotNames.at(i));

    std::cout << "mu: " << truemu->Integral() << std::endl;
    std::cout << "p: " << truep->Integral() << std::endl;
    std::cout << "pi: " << truepi->Integral() << std::endl;
    std::cout << "k: " << truek->Integral() << std::endl;
    std::cout << "all: " << hall->Integral() << std::endl;

    truemu->RebinX(4);
    truep->RebinX(4);
    truepi->RebinX(4);
    truek->RebinX(4);
    hall->RebinX(4);

    styleHistogram(truemu, 0);
    styleHistogram(truep, 1);
    styleHistogram(truepi, 2);
    styleHistogram(truek, 3);

    TH1D* hTotal = (TH1D*)truemu->Clone("hTotal");
    hTotal->Add(truep);
    hTotal->Add(truepi);
    hTotal->Add(truek);

    truemu->Scale(1./hall->Integral());
    truep->Scale(1./hall->Integral());
    truepi->Scale(1./hall->Integral());
    truek->Scale(1./hall->Integral());
    //hall->Scale(1./hall->Integral());

    THStack* hs = new THStack("hs", "");
    hs->Add(truep);
    hs->Add(truemu);
    hs->Add(truepi);
    hs->Add(truek);

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

    hTotal->SetTitle(";"+BDPlotNames.at(i)+";");
    TH1D* hallc = (TH1D*)hall->Clone("hallc");

    hallc->SetFillStyle(3345);
    hallc->SetFillColor(kBlack);
    hall->GetYaxis()->SetRangeUser(0, hTotal->GetMaximum()*1.25);
    hall->DrawNormalized("");
    hs->Draw("same");
    hall->DrawNormalized("same");
    hallc->DrawNormalized("sameE2");
    honminusoff->DrawNormalized("samepE1");

    TLegend *leg = new TLegend(0.15, 0.7, 0.5, .89);
    leg->AddEntry(honminusoff, "OnBeam-OffBeam");
    leg->AddEntry(truek, "Kaon");
    leg->AddEntry(truepi, "Pion");
    leg->AddEntry(truemu, "Muon");
    leg->AddEntry(truep, "Proton");
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->Draw("same");

    c1->SaveAs(BDPlotNames.at(i)+".png");

    c1->Delete();
    hs->Delete();
    truemu->Delete();
    truep->Delete();
    truepi->Delete();
    truek->Delete();
    leg->Delete();
    hon->Reset();
    hoff->Reset();

  }
}
