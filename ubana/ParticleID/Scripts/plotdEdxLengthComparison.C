/**
 * Dear future Adam and Kirsty:
 * I'm sorry this is a mess. I'm sorry that it uses a copy/pasted version of the truncated mean code.
 * I wanted a quick and dirty test of how well the truncated mean/length PID worked.
 *
 * summary looks to be that the separation is *really* excellent but the data/mc agreement is pretty
 * bad, so I don't trust it.
 *
 * Anyway to run, use:
 * root -l bnbcos.root onbeam.root offbeam.root plotdEdxLengthComparison.C
 *
 * You need to update the scale factors in the code but I've put them right at the top so they're
 * easy to find
 */ 

#include "truncatedMean.C"

std::pair<double,double> GetChi2(TH1D *o, TH1D* e){

  double chi2 = 0;
  double dof = 0;

  for (int i = 0; i < o->GetNbinsX(); i++){

    if (e->GetBinContent(i) != 0){

      chi2 += std::pow(o->GetBinContent(i) - e->GetBinContent(i),2)/e->GetBinContent(i);

    }
    if (e->GetBinContent(i) != 0 || o->GetBinContent(i) !=0){

      dof++;

    }

  }

  std::pair<double,double> returner;
  returner.first = chi2;
  returner.second = dof;
  return returner;

}

std::vector<double> templateFit(TH2D* mu, TH2D* p, TH2D* pi, TH2D* k, TH2D* other, TH2D* onbeam, TH2D* offbeam, double offBeamScaling, double potScaling){

  double fracFloat = 0.0;

  // scale offbeam to onbeam
  offbeam->Sumw2();
  offbeam->Scale(offBeamScaling);

  // create onbeam minus offbeam
  TH2D* data = (TH2D*)onbeam->Clone("data");
  data->Add(offbeam, -1);

  // check kaons exist in the sample
  bool isK = true;
  if (k->Integral() == 0) isK = false;

  // get total integral;
  TH2D* total = (TH2D*)mu->Clone("total");
  total->Add(p);
  total->Add(pi);
  total->Add(other);
  total->Add(k);

  // get current fractions
  double muFrac = mu->Integral()/total->Integral();
  double pFrac  = p->Integral()/total->Integral();
  double piFrac  = pi->Integral()/total->Integral();
  double kFrac = 0;
  if (isK) kFrac = k->Integral()/total->Integral();
  double otherFrac  = other->Integral()/total->Integral();

  // scale mchists to POT normally
  mu->Sumw2();
  p->Sumw2();
  pi->Sumw2();
  other->Sumw2();
  if (isK) k->Sumw2();

  mu->Scale(potScaling);
  p->Scale(potScaling);
  pi->Scale(potScaling);
  other->Scale(potScaling);
  if (isK) k->Scale(potScaling);

  // put mchist objects in TObjArray
  TObjArray *mc = new TObjArray(4);
  mc->Add(mu);
  mc->Add(p);
  mc->Add(pi);
  mc->Add(other);
  if (isK) mc->Add(k);

  // call fitting function
  TFractionFitter* fit = new TFractionFitter(data, mc);
  fit->Constrain(0, muFrac * fracFloat, muFrac * (1+fracFloat));
  fit->Constrain(1, pFrac * fracFloat, pFrac * (1+fracFloat));
  fit->Constrain(2, piFrac * fracFloat, piFrac * (1+fracFloat));
  fit->Constrain(3, otherFrac * fracFloat, otherFrac * (1+fracFloat));
  if (isK) fit->Constrain(4, kFrac * fracFloat, kFrac * (1+fracFloat));

  int stat = fit->Fit();

  double resmu = 0;
  double resmuerr = 0;
  fit->GetResult(0, resmu, resmuerr);

  double resp = 0;
  double resperr = 0;
  fit->GetResult(1, resp, resperr);

  double respi = 0;
  double respierr = 0;
  fit->GetResult(2, respi, respierr);

  double resother = 0;
  double resothererr = 0;
  fit->GetResult(3, resother, resothererr);

  double resk = 0;
  double reskerr = 0;
  if (isK) fit->GetResult(4, resk, reskerr);

  std::vector<double> returner;
  returner.push_back(resmu);
  returner.push_back(resp);
  returner.push_back(respi);
  returner.push_back(resk);
  returner.push_back(resother);

  return returner;

}

/**
 * code to convert dQdx to dE/dx. Not actually used in this script
 * except to check that I did the dE/dx to dQ/dx conversion correctly
 */ 
double convertdQdxTodEdx(double dQdx){
  double rho      = 1.383;                    // LAr density in g/cm^3
  double Wion     = 1000./4.237e7;        // 23.6 eV = 1e, Wion in MeV/e
  double E_field  = 0.273;                           // Electric Field in the drift region in KV/cm
  double Beta     = 0.212 / (rho * E_field);
  double Alpha    = 0.930;
  double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;

  return dEdx;
}

/**
 * convert dE/dx (MeV/cm) to dQ/dx (e-/cm). This is done because we don't save dQ/dx because it
 * is a silly variable that nobody should ever use, but this PID makes use of the <dQ/dx>_{tr}
 */ 
double convertdEdxTodQdx(double dEdx){
  double rho      = 1.383;                    // LAr density in g/cm^3
  double Wion     = 1000./4.237e7;        // 23.6 eV = 1e, Wion in MeV/e
  double E_field  = 0.273;                           // Electric Field in the drift region in KV/cm
  double Beta     = 0.212 / (rho * E_field);
  double Alpha    = 0.930;

  double dQdx = (std::log((dEdx*Beta)+Alpha))/(Beta*Wion);

  return dQdx;
}

/**
 * helper function which converts a vector of dE/dx values to a vector of
 * dQ/dx values
 */
std::vector<double> getdQdxVector(std::vector<double>* dEdxVals){

  std::vector<double> dQdxVals;

  for (size_t i = 0; i < dEdxVals->size(); i++){

    dQdxVals.push_back(convertdEdxTodQdx(dEdxVals->at(i)));

  }

  return dQdxVals;

}

/**
 * helper function to convert ADC/cm to e-/cm. This isn't needed because
 * the above function already converts to e-/cm.
 */ 
double convertADCToE(double dQdx){

  return dQdx; //* 196.98;

}

/**
 * main
 */ 
void plotdEdxLengthComparison(){

  // length cuts are used to look at a slice of the
  // 2D distribution in dQ/dx space.
  std::vector<double> track_length_high_cut = {700, 5, 10, 15, 20, 25, 30};
  std::vector<double> track_length_low_cut  = {0  , 0,  5, 10, 15, 20, 25};
  double off_beam_scaling = 0.77;
  double mc_scaling = 2.7;
  double proton_dqdx_scaling = 1.0;
  double proton_normalisation_scaling = 1.0;
  int highval = 300000;
  bool isTemplateFit = true;

  TTree *tree_mc = (TTree*)_file0->Get("pidvalid/pidTree");

  for (int len = 0; len < track_length_low_cut.size(); len++){

    /** BNB+COS MC plots and variables */
    double track_length_mc;
    std::vector<double> *track_dEdx_mc = 0;
    std::vector<double> *track_dEdx_mc_perhit_y = 0;
    int true_PDG;

    tree_mc->SetBranchAddress("track_length", &track_length_mc);
    tree_mc->SetBranchAddress("track_dEdx"  , &track_dEdx_mc);
    tree_mc->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_mc_perhit_y);
    tree_mc->SetBranchAddress("true_PDG"    , &true_PDG);

    TTree *tree_data_onbeam = (TTree*)_file1->Get("pidvalid/pidTree");

    TH2D* h_dEdx_length_proton = new TH2D("h_dEdx_length_proton", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 700);
    TH2D* h_dEdx_length_muon   = new TH2D("h_dEdx_length_muon", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 700);
    TH2D* h_dEdx_length_pion   = new TH2D("h_dEdx_length_pion", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 700);
    TH2D* h_dEdx_length_kaon   = new TH2D("h_dEdx_length_kaon", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 700);
    TH2D* h_dEdx_length_other  = new TH2D("h_dEdx_length_other", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 700);
    TH2D* h_dEdx_length_total  = new TH2D("h_dEdx_length_total", ";<dQ/dx>_{tr} (e-/cm);track length(cm)", 50, -1, highval, 50, 0, 700);

    TH1D* h_truemu = new TH1D("h_truemu", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_truep = new TH1D("h_truep", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_truepi = new TH1D("h_truepi", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_truek = new TH1D("h_truek", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_trueother = new TH1D("h_trueother", ";<dQdx>_{tr};", 50, 0, highval);

    TH1D* h_muID_truemu = new TH1D("h_muID_truemu", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_muID_truep = new TH1D("h_muID_truep", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_muID_truepi = new TH1D("h_muID_truepi", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_muID_truek = new TH1D("h_muID_truek", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_muID_trueother = new TH1D("h_muID_trueother", ";<dQdx>_{tr};", 50, 0, highval);

    TH1D* h_pID_truemu = new TH1D("h_pID_truemu", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_truep = new TH1D("h_pID_truep", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_truepi = new TH1D("h_pID_truepi", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_truek = new TH1D("h_pID_truek", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_trueother = new TH1D("h_pID_trueother", ";<dQdx>_{tr};", 50, 0, highval);

    /** On-beam plots and variables */
    double track_length_data_onbeam;
    std::vector<double> *track_dEdx_data_onbeam = 0;
    std::vector<double> *track_dEdx_data_onbeam_perhit_y = 0;

    tree_data_onbeam->SetBranchAddress("track_length", &track_length_data_onbeam);
    tree_data_onbeam->SetBranchAddress("track_dEdx"  , &track_dEdx_data_onbeam);
    tree_data_onbeam->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_data_onbeam_perhit_y);

    TTree *tree_data_offbeam = (TTree*)_file2->Get("pidvalid/pidTree");

    TH2D* h_dEdx_length_data_onbeam = new TH2D("h_dEdx_length_data_onbeam", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 700);
    TH1D* h_muID_data_onbeam = new TH1D("h_muID_data_onbeam", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_data_onbeam = new TH1D("h_pID_data_onbeam", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_data_onbeam = new TH1D("h_data_onbeam", ";<dQdx>_{tr};", 50, 0, highval);

    /** off-beam plots and variables */
    double track_length_data_offbeam;
    std::vector<double> *track_dEdx_data_offbeam = 0;
    std::vector<double> *track_dEdx_data_offbeam_perhit_y = 0;

    tree_data_offbeam->SetBranchAddress("track_length", &track_length_data_offbeam);
    tree_data_offbeam->SetBranchAddress("track_dEdx"  , &track_dEdx_data_offbeam);
    tree_data_offbeam->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_data_offbeam_perhit_y);

    TH2D* h_dEdx_length_data_offbeam = new TH2D("h_dEdx_length_data_offbeam", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 700);
    TH1D* h_muID_data_offbeam = new TH1D("h_muID_data_offbeam", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_pID_data_offbeam = new TH1D("h_pID_data_offbeam", ";<dQdx>_{tr};", 50, 0, highval);
    TH1D* h_data_offbeam = new TH1D("h_data_offbeam", ";<dQdx>_{tr};", 50, 0, highval);

    /** loop bnb+cos mc */
    for (int i = 0; i < tree_mc->GetEntries(); i++){


      std::fstream infile("dQdxSeparator2.txt");
      tree_mc->GetEntry(i);

      /** convert dEdx values to dQdx values */
      std::vector<double> track_dQdx_mc_perhit_y = getdQdxVector(track_dEdx_mc_perhit_y);

      if (track_dQdx_mc_perhit_y.size() == 0) continue;

      double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_mc_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
      double dQdxInE = convertADCToE(truncatedMeandQdx);

      // get cut value
      std::string s;
      int rnd_trklen = std::round(track_length_mc);
      int cutvalue = 0;
      int j = 1;
      while (std::getline(infile,s)){
        if (j == rnd_trklen){
          cutvalue = std::atof(s.c_str());
          break;
        }
        j++;
      }

      infile.close();

      h_dEdx_length_total->Fill(dQdxInE*proton_dqdx_scaling, track_length_mc);

      if (std::abs(true_PDG) == 2212){
        h_dEdx_length_proton->Fill(dQdxInE*proton_dqdx_scaling, track_length_mc);
        if( track_length_mc > track_length_low_cut.at(len) && track_length_mc < track_length_high_cut.at(len)){
          h_truep->Fill(dQdxInE*proton_dqdx_scaling);
          if (dQdxInE*proton_dqdx_scaling < cutvalue)
            h_muID_truep->Fill(dQdxInE*proton_dqdx_scaling);
          else h_pID_truep->Fill(dQdxInE*proton_dqdx_scaling);
        }

      }
      else if (std::abs(true_PDG) == 13){
        h_dEdx_length_muon->Fill(dQdxInE, track_length_mc);

        if( track_length_mc > track_length_low_cut.at(len) && track_length_mc < track_length_high_cut.at(len)){
          h_truemu->Fill(dQdxInE);
          if (dQdxInE < cutvalue)
            h_muID_truemu->Fill(dQdxInE);
          else h_pID_truemu->Fill(dQdxInE);
        }

      }
      else if (std::abs(true_PDG) == 211){
        h_dEdx_length_pion->Fill(dQdxInE, track_length_mc);

        if( track_length_mc > track_length_low_cut.at(len) && track_length_mc < track_length_high_cut.at(len)){
          h_truepi->Fill(dQdxInE);
          if (dQdxInE < cutvalue)
            h_muID_truepi->Fill(dQdxInE);
          else h_pID_truepi->Fill(dQdxInE);
        }

      }

      else if (std::abs(true_PDG) == 321){
        h_dEdx_length_kaon->Fill(dQdxInE, track_length_mc);

        if( track_length_mc > track_length_low_cut.at(len) && track_length_mc < track_length_high_cut.at(len)){
          h_truek->Fill(dQdxInE);
          if (dQdxInE < cutvalue)
            h_muID_truek->Fill(dQdxInE);
          else h_pID_truek->Fill(dQdxInE);

        }

      }
      else{
        h_dEdx_length_other->Fill(dQdxInE, track_length_mc);

        if( track_length_mc > track_length_low_cut.at(len) && track_length_mc < track_length_high_cut.at(len)){
          h_trueother->Fill(dQdxInE);
          if (dQdxInE < cutvalue)
            h_muID_trueother->Fill(dQdxInE);
          else h_pID_trueother->Fill(dQdxInE);

        }

      }

    }

    /** loop onbeam */
    for (int i = 0; i < tree_data_onbeam->GetEntries() ; i++){

      std::fstream infile("dQdxSeparator2.txt");
      tree_data_onbeam->GetEntry(i);

      if (track_dEdx_data_onbeam_perhit_y->size() == 0)
        continue;

      /** convert dEdx values to dQdx values */
      std::vector<double> track_dQdx_data_onbeam_perhit_y = getdQdxVector(track_dEdx_data_onbeam_perhit_y);

      double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_data_onbeam_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
      double dQdxInE = convertADCToE(truncatedMeandQdx);

      // get cut value
      std::string s;
      int rnd_trklen = std::round(track_length_data_onbeam);
      int cutvalue = 0;
      int j = 1;
      while (std::getline(infile,s)){
        if (j == rnd_trklen){
          cutvalue = std::atof(s.c_str());
          break;
        }
        j++;
      }

      infile.close();

      h_dEdx_length_data_onbeam->Fill(dQdxInE, track_length_data_onbeam);
      if ( track_length_data_onbeam > track_length_low_cut.at(len) && track_length_data_onbeam < track_length_high_cut.at(len)){
        h_data_onbeam->Fill(dQdxInE);
        if (dQdxInE < cutvalue)
          h_muID_data_onbeam->Fill(dQdxInE);
        else h_pID_data_onbeam->Fill(dQdxInE);
      }
    }

    /** loop offbeam */
    for (int i = 0; i < tree_data_offbeam->GetEntries() ; i++){

      std::fstream infile("dQdxSeparator2.txt");
      tree_data_offbeam->GetEntry(i);

      if (track_dEdx_data_offbeam_perhit_y->size() == 0)
        continue;

      /** convert dEdx values to dQdx values */
      std::vector<double> track_dQdx_data_offbeam_perhit_y = getdQdxVector(track_dEdx_data_offbeam_perhit_y);

      double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_data_offbeam_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
      double dQdxInE = convertADCToE(truncatedMeandQdx);

      // get cut value
      std::string s;
      int rnd_trklen = std::round(track_length_data_offbeam);
      int cutvalue = 0;
      int j = 1;
      while (std::getline(infile,s)){
        if (j == rnd_trklen){
          cutvalue = std::atof(s.c_str());
          break;
        }
        j++;
      }

      infile.close();

      h_dEdx_length_data_offbeam->Fill(dQdxInE, track_length_data_offbeam);

      if ( track_length_data_onbeam > track_length_low_cut.at(len) && track_length_data_onbeam < track_length_high_cut.at(len)){
        h_data_offbeam->Fill(dQdxInE);
        if (dQdxInE < cutvalue)
          h_muID_data_offbeam->Fill(dQdxInE);
        else h_pID_data_offbeam->Fill(dQdxInE);
      }
    }

    h_dEdx_length_proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
    h_dEdx_length_muon->SetMarkerColor(TColor::GetColor(8, 64, 129));
    h_dEdx_length_pion->SetMarkerColor(TColor::GetColor(166, 217, 106));
    h_dEdx_length_kaon->SetMarkerColor(TColor::GetColor(133, 1, 98));
    h_dEdx_length_other->SetMarkerColor(TColor::GetColor(197, 197, 197));

    TCanvas *c1 = new TCanvas();

    h_dEdx_length_proton->Draw();
    h_dEdx_length_muon->Draw("same");
    h_dEdx_length_pion->Draw("same");
    h_dEdx_length_kaon->Draw("same");
    h_dEdx_length_other->Draw("same");

    c1->SaveAs("truncateddedxlength_2d.png");

    TCanvas *c2 = new TCanvas();
    h_muID_truep->SetFillColor(TColor::GetColor(215, 48, 39));
    h_muID_truemu->SetFillColor(TColor::GetColor(8,64,129));
    h_muID_truepi->SetFillColor(TColor::GetColor(166,217,106));
    h_muID_truek->SetFillColor(TColor::GetColor(133, 1, 98));
    h_muID_trueother->SetFillColor(TColor::GetColor(197, 197, 197));
    h_muID_data_offbeam->SetFillColor(kBlack);
    h_muID_data_offbeam->SetFillStyle(3345);

    std::vector<double> scalings;

    double premu;
    double prep;
    double prepi;
    double prek;
    double preoth;
    double postmu;
    double postp;
    double postpi;
    double postk;
    double postoth;


    if (isTemplateFit){
      std::cout << "pre muon: "   << h_dEdx_length_muon->Integral()/h_dEdx_length_total->Integral()   << std::endl;
      std::cout << "pre proton: " << h_dEdx_length_proton->Integral()/h_dEdx_length_total->Integral() << std::endl;
      std::cout << "pre pion: "   << h_dEdx_length_pion->Integral()/h_dEdx_length_total->Integral()   << std::endl;
      std::cout << "pre kaon: "   << h_dEdx_length_kaon->Integral()/h_dEdx_length_total->Integral()   << std::endl;
      std::cout << "pre other: "  << h_dEdx_length_other->Integral()/h_dEdx_length_total->Integral()  << std::endl;

      premu  = h_dEdx_length_muon->Integral()/h_dEdx_length_total->Integral();
      prep   = h_dEdx_length_proton->Integral()/h_dEdx_length_total->Integral(); 
      prepi  = h_dEdx_length_pion->Integral()/h_dEdx_length_total->Integral();   
      prek   = h_dEdx_length_kaon->Integral()/h_dEdx_length_total->Integral();   
      preoth = h_dEdx_length_other->Integral()/h_dEdx_length_total->Integral();  

      scalings = templateFit(h_dEdx_length_muon, h_dEdx_length_proton, h_dEdx_length_pion, h_dEdx_length_kaon, h_dEdx_length_other, h_dEdx_length_data_onbeam, h_dEdx_length_data_offbeam, off_beam_scaling, mc_scaling);

      std::cout << "post muon: "   << scalings.at(0) << std::endl;
      std::cout << "post proton: " << scalings.at(1) << std::endl;
      std::cout << "post pion: "   << scalings.at(2) << std::endl;
      std::cout << "post kaon: "   << scalings.at(3) << std::endl;
      std::cout << "post other: "  << scalings.at(4) << std::endl;

      postmu  = scalings.at(0);
      postp   = scalings.at(1);
      postpi  = scalings.at(2);
      postk   = scalings.at(3);
      postoth = scalings.at(4);

    }

    h_dEdx_length_data_offbeam->Sumw2();
    h_muID_data_offbeam->Sumw2();
    h_muID_truemu->Sumw2();
    h_muID_truep->Sumw2();
    h_muID_truepi->Sumw2();
    h_muID_truek->Sumw2();
    h_muID_trueother->Sumw2();

    h_dEdx_length_data_offbeam->Scale(off_beam_scaling);
    h_muID_data_offbeam->Scale(off_beam_scaling);
    double dataInt = h_dEdx_length_data_onbeam->Integral()-h_dEdx_length_data_offbeam->Integral();

    if (isTemplateFit){
      h_muID_truemu->Scale((postmu/premu) * dataInt/h_dEdx_length_total->Integral());
      h_muID_truep->Scale((postp/prep) * dataInt/h_dEdx_length_total->Integral());
      h_muID_truepi->Scale((postpi/prepi) * dataInt/h_dEdx_length_total->Integral());
      if (h_muID_truek->Integral() != 0) h_muID_truek->Scale((postk/prek) * dataInt/h_dEdx_length_total->Integral());
      h_muID_trueother->Scale((postoth/preoth) * dataInt/h_dEdx_length_total->Integral());
    }
    else {
      h_muID_truep->Scale(mc_scaling*proton_normalisation_scaling);
      h_muID_truemu->Scale(mc_scaling);
      h_muID_truepi->Scale(mc_scaling);
      h_muID_truek->Scale(mc_scaling);
      h_muID_trueother->Scale(mc_scaling);
    }

    THStack *hs = new THStack("hs", ";<dQ/dx>_{tr};");
    hs->Add(h_muID_data_offbeam);
    hs->Add(h_muID_truep);
    hs->Add(h_muID_truemu);
    hs->Add(h_muID_truepi);
    hs->Add(h_muID_truek);
    hs->Add(h_muID_trueother);

    TH1D* hsTot = (TH1D*)h_muID_data_offbeam->Clone("hsTot");
    hsTot->Add(h_muID_truep);
    hsTot->Add(h_muID_truemu);
    hsTot->Add(h_muID_truepi);
    hsTot->Add(h_muID_truek);
    hsTot->Add(h_muID_trueother);

    h_muID_data_onbeam->GetYaxis()->SetRangeUser(0, std::max(hs->GetMaximum(),h_muID_data_onbeam->GetMaximum())*1.1);
    h_muID_data_onbeam->Draw("pE0");
    hs->Draw("histsame");
    h_muID_data_onbeam->SetMarkerStyle(20);
    h_muID_data_onbeam->SetMarkerSize(0.6);
    h_muID_data_onbeam->Draw("psameE0");

    std::cout << h_muID_data_onbeam->Integral() << " " << hs->GetHistogram()->Integral() << std::endl;

    std::pair<double,double> chi2_muid = GetChi2(h_muID_data_onbeam, hsTot);

    TPaveText *pt3 = new TPaveText(0.17, 0.87, 0.42, 0.92, "NDC");
    TString chi2string_muid = Form("Chi2/NDF: %.2f/%g", chi2_muid.first, chi2_muid.second);
    pt3->SetFillColor(kWhite);
    pt3->AddText(((std::string)chi2string_muid).c_str());
    pt3->SetTextSize(0.05);

    pt3->Draw("same");

    TString saveString2 = Form("truncdedxlength_muID_%ito%icm.png", (int)track_length_low_cut.at(len), (int)track_length_high_cut.at(len));
    c2->SaveAs(saveString2);

    TCanvas *c3 = new TCanvas();
    h_pID_truep->SetFillColor(TColor::GetColor(215, 48, 39));
    h_pID_truemu->SetFillColor(TColor::GetColor(8,64,129));
    h_pID_truepi->SetFillColor(TColor::GetColor(166,217,106));
    h_pID_truek->SetFillColor(TColor::GetColor(133, 1, 98));
    h_pID_trueother->SetFillColor(TColor::GetColor(197, 197, 197));
    h_pID_data_offbeam->SetFillColor(kBlack);
    h_pID_data_offbeam->SetFillStyle(3345);

    h_pID_data_offbeam->Sumw2();
    h_pID_truemu->Sumw2();
    h_pID_truep->Sumw2();
    h_pID_truepi->Sumw2();
    h_pID_truek->Sumw2();
    h_pID_trueother->Sumw2();

    h_pID_data_offbeam->Scale(off_beam_scaling);

    if (isTemplateFit){
      h_pID_truemu->Scale((postmu/premu) * dataInt/h_dEdx_length_total->Integral());
      h_pID_truep->Scale((postp/prep) * dataInt/h_dEdx_length_total->Integral());
      h_pID_truepi->Scale((postpi/prepi) * dataInt/h_dEdx_length_total->Integral());
      if (h_pID_truek->Integral() !=0) h_pID_truek->Scale((postk/prek) * dataInt/h_dEdx_length_total->Integral());
      h_pID_trueother->Scale((postoth/preoth) * dataInt/h_dEdx_length_total->Integral());
    }
    else{
      h_pID_truep->Scale(mc_scaling*proton_normalisation_scaling);
      h_pID_truemu->Scale(mc_scaling);
      h_pID_truepi->Scale(mc_scaling);
      h_pID_truek->Scale(mc_scaling);
      h_pID_trueother->Scale(mc_scaling);
    }

    THStack *hs2 = new THStack("hs", ";<dQ/dx>_{tr};");
    hs2->Add(h_pID_data_offbeam);
    hs2->Add(h_pID_truep);
    hs2->Add(h_pID_truemu);
    hs2->Add(h_pID_truepi);
    hs2->Add(h_pID_truek);
    hs2->Add(h_pID_trueother);

    TH1D* hs2Tot = (TH1D*)h_pID_data_offbeam->Clone("hs2Tot");
    hs2Tot->Add(h_pID_truep);
    hs2Tot->Add(h_pID_truemu);
    hs2Tot->Add(h_pID_truepi);
    hs2Tot->Add(h_pID_truek);
    hs2Tot->Add(h_pID_trueother);

    h_pID_data_onbeam->GetYaxis()->SetRangeUser(0, std::max(hs2->GetMaximum(),h_pID_data_onbeam->GetMaximum())*1.1);
    h_pID_data_onbeam->Draw("pe0");
    hs2->Draw("histsame");
    h_pID_data_onbeam->SetMarkerStyle(20);
    h_pID_data_onbeam->SetMarkerSize(0.6);
    h_pID_data_onbeam->Draw("psameE0");

    std::cout << "number of protons: " << h_muID_truep->Integral() + h_pID_truep->Integral() << std::endl;
    std::cout << "proton Efficiency: " << h_pID_truep->Integral() / (h_muID_truep->Integral() + h_pID_truep->Integral()) << std::endl;
    std::cout << "proton purity:     " << h_pID_truep->Integral() / (h_pID_truemu->Integral() + h_pID_truepi->Integral() + h_pID_truek->Integral() + h_pID_trueother->Integral() + h_pID_truep->Integral()) << std::endl;

    std::cout << "number of muons: " << h_muID_truemu->Integral() + h_pID_truemu->Integral() << std::endl;
    std::cout << "muon Efficiency: " << h_muID_truemu->Integral() / (h_muID_truemu->Integral() + h_pID_truemu->Integral()) << std::endl;
    std::cout << "muon purity:     " << h_muID_truemu->Integral() / (h_muID_truemu->Integral() + h_muID_truepi->Integral() + h_muID_truek->Integral() + h_muID_trueother->Integral() + h_muID_truep->Integral()) << std::endl;

    std::pair<double,double> chi2_pid = GetChi2(h_pID_data_onbeam, hs2Tot);

    TPaveText *pt2 = new TPaveText(0.17, 0.87, 0.42, 0.92, "NDC");
    TString chi2string_pid = Form("Chi2/NDF: %.2f/%g", chi2_pid.first, chi2_pid.second);
    pt2->SetFillColor(kWhite);
    pt2->AddText(((std::string)chi2string_pid).c_str());
    pt2->SetTextSize(0.05);

    pt2->Draw("same");

    TString saveString3 = Form("truncdedxlength_pID_%ito%icm.png", (int)track_length_low_cut.at(len), (int)track_length_high_cut.at(len));
    c3->SaveAs(saveString3);

    TCanvas *c8 = new TCanvas();
    h_truep->SetFillColor(TColor::GetColor(215, 48, 39));
    h_truemu->SetFillColor(TColor::GetColor(8,64,129));
    h_truepi->SetFillColor(TColor::GetColor(166,217,106));
    h_truek->SetFillColor(TColor::GetColor(133, 1, 98));
    h_trueother->SetFillColor(TColor::GetColor(197, 197, 197));

    h_data_offbeam->Sumw2();
    h_truemu->Sumw2();
    h_truep->Sumw2();
    h_truepi->Sumw2();
    h_truek->Sumw2();
    h_trueother->Sumw2();

    h_data_offbeam->Scale(off_beam_scaling);
    if (isTemplateFit){
      h_truemu->Scale((postmu/premu) * dataInt/h_dEdx_length_total->Integral());
      h_truep->Scale((postp/prep) * dataInt/h_dEdx_length_total->Integral());
      h_truepi->Scale((postpi/prepi) * dataInt/h_dEdx_length_total->Integral());
      if (h_truek->Integral() !=0) h_truek->Scale((postk/prek) * dataInt/h_dEdx_length_total->Integral());
      h_trueother->Scale((postoth/preoth) * dataInt/h_dEdx_length_total->Integral());
    }
    else{
      h_truep->Scale(mc_scaling*proton_normalisation_scaling);
      h_truemu->Scale(mc_scaling);
      h_truepi->Scale(mc_scaling);
      h_truek->Scale(mc_scaling);
      h_trueother->Scale(mc_scaling);
    }

    THStack *hs3 = new THStack("hs3", ";<dQ/dx>_{tr};");
    hs3->Add(h_data_offbeam);
    hs3->Add(h_truep);
    hs3->Add(h_truemu);
    hs3->Add(h_truepi);
    hs3->Add(h_truek);
    hs3->Add(h_trueother);

    TH1D* hs3Tot = (TH1D*)h_muID_data_offbeam->Clone("hs3Tot");
    hs3Tot->Add(h_truep);
    hs3Tot->Add(h_truemu);
    hs3Tot->Add(h_truepi);
    hs3Tot->Add(h_truek);
    hs3Tot->Add(h_trueother);

    h_data_onbeam->GetYaxis()->SetRangeUser(0, std::max(hs3->GetMaximum(),h_data_onbeam->GetMaximum())*1.1);
    h_data_onbeam->Draw("pe0");
    hs3->Draw("histsame");
    h_data_onbeam->SetMarkerStyle(20);
    h_data_onbeam->SetMarkerSize(0.6);
    h_data_onbeam->Draw("psamee0");

    std::pair<double,double> chi2 = GetChi2(h_data_onbeam, hs3Tot);

    TPaveText *pt = new TPaveText(0.17, 0.87, 0.42, 0.92, "NDC");
    TString chi2string = Form("Chi2/NDF: %.2f/%g", chi2.first, chi2.second);
    pt->SetFillColor(kWhite);
    pt->AddText(((std::string)chi2string).c_str());
    pt->SetTextSize(0.05);

    pt->Draw("same");

    TString saveString8 = Form("truncdedxlength_%ito%icm.png", (int)track_length_low_cut.at(len), (int)track_length_high_cut.at(len));
    c8->SaveAs(saveString8);

    TCanvas *c4 = new TCanvas();

    TH2D* hAdded = new TH2D("hAdded", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 700);
    hAdded->Add(h_dEdx_length_muon);
    hAdded->Add(h_dEdx_length_pion);
    hAdded->Add(h_dEdx_length_kaon);
    hAdded->Add(h_dEdx_length_other);
    h_dEdx_length_proton->Sumw2();
    h_dEdx_length_proton->Scale(mc_scaling*proton_normalisation_scaling);
    hAdded->Add(h_dEdx_length_proton);

    hAdded->Sumw2();
    h_dEdx_length_data_onbeam->Sumw2();

    hAdded->Scale(1./hAdded->Integral());
    h_dEdx_length_data_onbeam->Add(h_dEdx_length_data_offbeam, -1.);
    h_dEdx_length_data_onbeam->Scale(1./h_dEdx_length_data_onbeam->Integral());

    h_dEdx_length_data_onbeam->Divide(hAdded);
    h_dEdx_length_data_onbeam->GetZaxis()->SetRangeUser(0,2);
    h_dEdx_length_data_onbeam->Draw("colz");
  }
}
