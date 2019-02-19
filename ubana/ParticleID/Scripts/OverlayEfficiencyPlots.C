std::vector<int> marker_styles = {20,24,34};
std::vector<double> marker_sizes = {0.4,0.4,0.5};

void OverlayEfficiencyPlots(std::vector<std::string> inputfiles, std::vector<std::string> filedescriptions={}){

  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0.);

  // Check size of inputfiles vector
  std::cout << "Overlaying efficiency/purity curves from " << inputfiles.size() << " files." << std::endl;
  if (inputfiles.size() == 0){
    std::cout << "You need to give me some files! Exiting." << std::endl;\
    return;
  }

  TFile *file0 = new TFile(inputfiles.at(0).c_str(),"open");

  // Get list of plots
  // We want one entry in this list for each plane/variable
  // We *do not* want separate entries for efficiency, purity, and efficiency*purity, and we do not want separate entries for muons, pions, and protons
  // So: let's just ask for unique names that start with heff_ and end with _mu
  std::vector<std::string> plotnames;
  TIter next(file0->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    // Get TH1s only
    if (!TString(key->GetClassName()).Contains("TH1")) continue;
    TString name = TString(key->GetName());
    if (!name.Contains("heff_")) continue;
    if (!name.EndsWith("_mu")) continue;
    name.Remove(0,5);
    name.Remove(TString::kTrailing,'u');
    name.Remove(TString::kTrailing,'m');
    name.Remove(TString::kTrailing,'_');
    plotnames.push_back(std::string(name.Data()));
  }

  // for (int i=0; i<plotnames.size(); i++){
  //   std::cout << plotnames.at(i) << std::endl;
  // }

  // Now we have a list of all the plots we want to overlay, let's get going!
  std::vector<std::string> particles = {"_mu","_pi","_p","_k"};
  // Loop over plots
  for (size_t i_plot=0; i_plot<plotnames.size(); i_plot++){

    TCanvas *c1 = new TCanvas();
    c1->Divide(2,2,0.0005,0.0005);

    TLegend *l;

    // Loop over input files
    for (size_t i_file=0; i_file < inputfiles.size(); i_file++){

      TFile *ftmp = TFile::Open(inputfiles.at(i_file).c_str(),"open");
      //TFile *ftmp = new TFile(inputfiles.at(i_file).c_str(),"open");

      // Loop over particles (mu, pi, p)
      for (size_t i_p=0; i_p<particles.size(); i_p++){
        // std::cout << particles.at(i_p) << std::endl;
        // std::cout << inputfiles.at(i_file) << std::endl;
        // std::cout << TString::Format("heff_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data() << std::endl;

        TH1D *heff = (TH1D*)ftmp->Get(TString::Format("heff_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());
        TH1D *hpur = (TH1D*)ftmp->Get(TString::Format("hpur_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());
        TH1D *heffpur = (TH1D*)ftmp->Get(TString::Format("heffpur_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());

        if (!heff) std::cout << "Did not find heff in file " << inputfiles.at(i_file) << std::endl;
        if (!hpur) std::cout << "Did not find hpur in file " << inputfiles.at(i_file) << std::endl;
        if (!heffpur) std::cout << "Did not find heffpur in file " << inputfiles.at(i_file) << std::endl;

        heff->SetDirectory(0);
        hpur->SetDirectory(0);
        heffpur->SetDirectory(0);

        if (i_p==1 && i_file==0){
          double legendy=0.;
          double legendpts=0.;
          for (int i_bin=1; i_bin < heff->GetXaxis()->GetNbins()+1; i_bin++){
            if (heff->GetXaxis()->GetBinCenter(i_bin)>0.7*heff->GetXaxis()->GetBinUpEdge(heff->GetXaxis()->GetLast())){
              legendy+=heff->GetBinContent(i_bin);
              legendpts++;
            }
          }

          legendy = legendy/legendpts;
          if (legendy<0.6) l = new TLegend(0.5,0.6,0.82,0.88);
          else l = new TLegend(0.5,0.3,0.82,0.6);
          l->SetTextFont(132);
          l->SetLineColor(kWhite);
          l->SetFillColor(kWhite);

          TH1D *heff_clone = (TH1D*)heff->Clone("heff_clone");
          TH1D *hpur_clone = (TH1D*)hpur->Clone("hpur_clone");
          TH1D *heffpur_clone = (TH1D*)heffpur->Clone("heffpur_clone");

          heff_clone->SetDirectory(0);
          hpur_clone->SetDirectory(0);
          heffpur_clone->SetDirectory(0);

          l->AddEntry(heff_clone,"Efficiency","l");
          l->AddEntry(hpur_clone,"Purity","l");
          l->AddEntry(heffpur_clone,"Efficiency #times Purity","l");
        }

        int markerstyle = 20;
        if (marker_styles.size() < i_file) markerstyle = marker_styles.at(marker_styles.size()-1)+(i_file - marker_styles.size());
        else markerstyle = marker_styles.at(i_file);

        double markersize = 0.3;
        if (marker_sizes.size() >= i_file) markersize = marker_sizes.at(i_file);

        heff->SetFillStyle(0);
        heff->SetLineColor(kRed);
        heff->SetLineWidth(1);
        heff->SetMarkerColor(kRed);
        heff->SetMarkerStyle(markerstyle);
        heff->SetMarkerSize(markersize);

        hpur->SetFillStyle(0);
        hpur->SetLineColor(kBlue);
        hpur->SetLineWidth(1);
        hpur->SetMarkerColor(kBlue);
        hpur->SetMarkerStyle(markerstyle);
        hpur->SetMarkerSize(markersize);

        heffpur->SetFillStyle(0);
        heffpur->SetLineColor(kBlack);
        heffpur->SetLineWidth(1);
        heffpur->SetMarkerColor(kBlack);
        heffpur->SetMarkerStyle(markerstyle);
        heffpur->SetMarkerSize(markersize);

        heff->GetYaxis()->SetRangeUser(0,1.1);

        c1->cd(i_p+1);
        if (i_file==0){
          heff->Draw("lp");
          hpur->Draw("same lp");
          heffpur->Draw("same lp");
        }
        else{
          heff->Draw("same p");
          hpur->Draw("same p");
          heffpur->Draw("same p");
        }

        if (i_p==1){
          char *legentry = (char*)TString::Format("File %d",(int)i_file).Data();
          if (filedescriptions.size()!=0){
            legentry = (char*)filedescriptions.at(i_file).c_str();
          }

          if (i_file==0){
            l->AddEntry(heffpur,legentry,"lp");
          }
          else{
            l->AddEntry(heffpur,legentry,"p");
          }
        }
      } // end loop over particles (i_p)

      ftmp->Close();

    } // end loop over input files (i_file)

    c1->cd(2);
    l->Draw();
    c1->Print(TString::Format("%s.png",plotnames.at(i_plot).c_str()).Data());
    delete c1;
  } // end loop over plots (i_plot)

}
