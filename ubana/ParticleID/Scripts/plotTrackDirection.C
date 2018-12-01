/**
 * This was intended to show how well the likelihood fits can determine direction but 
 * it's super badly written. Should clean up later.
 */

void plotTrackDirection(){

  int chosenPdg = 211;
  std::string ptype("pi"); // also need to change this when you change above

  TTree *tree = (TTree*)_file0->Get("pidvalid/pidTree");

  TH1D* htot = new TH1D("htot", "", 1, -10, 1100);
  tree->Draw("track_start_z >> htot");

  std::cout << "total: " << htot->Integral() << std::endl;

  TH1D* hpp = new TH1D("hpp", "", 1, -10, 1100);
  TH1D* hpn = new TH1D("hpn", "", 1, -10, 1100);
  TH1D* hnp = new TH1D("hnp", "", 1, -10, 1100);
  TH1D* hnn = new TH1D("hnn", "", 1, -10, 1100);

  tree->Draw("track_start_z >> hpp", ("track_end_z - track_start_z > 0 && true_end_z - true_start_z > 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
  tree->Draw("track_start_z >> hnp", ("track_end_z - track_start_z > 0 && true_end_z - true_start_z < 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
  tree->Draw("track_start_z >> hpn", ("track_end_z - track_start_z < 0 && true_end_z - true_start_z > 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
  tree->Draw("track_start_z >> hnn", ("track_end_z - track_start_z < 0 && true_end_z - true_start_z < 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());

  TH1D* hpp_lc = new TH1D("hpp_lc", "", 1, -10, 1100);
  TH1D* hpp_li = new TH1D("hpp_li", "", 1, -10, 1100);
  TH1D* hnn_lc = new TH1D("hnn_lc", "", 1, -10, 1100);
  TH1D* hnn_li = new TH1D("hnn_li", "", 1, -10, 1100);
   
  tree->Draw("track_start_z >> hpp_lc",  ("track_end_z - track_start_z > 0 && true_end_z - true_start_z > 0 && track_likelihood_fwd_"+ptype+"[2] > track_likelihood_bwd_"+ptype+"[2] && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
   tree->Draw("track_start_z >> hpp_li", ("track_end_z - track_start_z > 0 && true_end_z - true_start_z > 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_"+ptype+"[2] < track_likelihood_bwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
  tree->Draw("track_start_z >> hnn_lc",  ("track_end_z - track_start_z < 0 && true_end_z - true_start_z < 0 && track_likelihood_fwd_"+ptype+"[2] < track_likelihood_bwd_"+ptype+"[2] && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
   tree->Draw("track_start_z >> hnn_li", ("track_end_z - track_start_z < 0 && true_end_z - true_start_z < 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_"+ptype+"[2] > track_likelihood_bwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]" ).c_str());

  TH1D* hnp_lc = new TH1D("hnp_lc", "", 1, -10, 1100);
  TH1D* hnp_li = new TH1D("hnp_li", "", 1, -10, 1100);
  TH1D* hpn_lc = new TH1D("hpn_lc", "", 1, -10, 1100);
  TH1D* hpn_li = new TH1D("hpn_li", "", 1, -10, 1100);
   
  tree->Draw("track_start_z >> hnp_lc",  ("track_end_z - track_start_z < 0 && true_end_z - true_start_z > 0 && track_likelihood_fwd_"+ptype+"[2] > track_likelihood_bwd_"+ptype+"[2] && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
   tree->Draw("track_start_z >> hnp_li", ("track_end_z - track_start_z < 0 && true_end_z - true_start_z > 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_"+ptype+"[2] < track_likelihood_bwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
  tree->Draw("track_start_z >> hpn_lc",  ("track_end_z - track_start_z > 0 && true_end_z - true_start_z < 0 && track_likelihood_fwd_"+ptype+"[2] < track_likelihood_bwd_"+ptype+"[2] && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]").c_str());
   tree->Draw("track_start_z >> hpn_li", ("track_end_z - track_start_z > 0 && true_end_z - true_start_z < 0 && true_PDG == "+std::to_string(chosenPdg)+" && track_likelihood_fwd_"+ptype+"[2] > track_likelihood_bwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_fwd_"+ptype+"[2] && track_likelihood_fwd_mip[2] < track_likelihood_bwd_"+ptype+"[2]" ).c_str());

  std::cout << "correct direction: " << hpp->Integral() + hnn->Integral() << std::endl;
  std::cout << ">> correctly identified: " << hpp_lc->Integral() + hnn_lc->Integral()<< std::endl;
  std::cout << ">> incorrectly identified: " << hpp_li->Integral() + hnn_li->Integral()<< std::endl;

  std::cout << "incorrect direction: " << hpn->Integral() + hnp->Integral() << std::endl;
  std::cout << ">> correctly identified: " << hpn_lc->Integral() + hnp_lc->Integral()<< std::endl;
  std::cout << ">> incorrectly identified: " << hpn_li->Integral() + hnp_li->Integral()<< std::endl;

  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  c1->SetRightMargin(0.12);

  TH2D* h2 = new TH2D("h2", "", 2, 0, 2, 2, 0, 2);
  h2->SetBinContent(1,1,std::round(100*(hpp_lc->Integral() + hnn_lc->Integral())/(hpp->Integral()+hnn->Integral()+hpn->Integral()+hnp->Integral()))/100.);
  h2->SetBinContent(1,2,std::round(100*(hpp_li->Integral() + hnn_li->Integral())/(hpp->Integral()+hnn->Integral()+hpn->Integral()+hnp->Integral()))/100.);
  h2->SetBinContent(2,1,std::round(100*(hpn_lc->Integral() + hnp_lc->Integral())/(hpn->Integral()+hnp->Integral()+hpp->Integral()+hnn->Integral()))/100.);
  h2->SetBinContent(2,2,std::round(100*(hpn_li->Integral() + hnp_li->Integral())/(hpn->Integral()+hnp->Integral()+hpp->Integral()+hnn->Integral()))/100.);
  h2->GetZaxis()->SetRangeUser(0,1);
  gStyle->SetPalette(kLightTemperature);

  const char* x[2] = {"I_{c}", "I_{i}"};
  const char* y[2] = {"L_{c}", "L_{i}"};

  h2->GetXaxis()->SetBinLabel(1, x[0]);
  h2->GetXaxis()->SetBinLabel(2, x[1]);
  h2->GetYaxis()->SetBinLabel(1, y[0]);
  h2->GetYaxis()->SetBinLabel(2, y[1]);

  h2->GetXaxis()->SetLabelSize(0.08);
  h2->GetYaxis()->SetLabelSize(0.08);

  h2->Draw("col");
  h2->DrawClone("textsame");

  c1->SaveAs(("correctDirection_"+ptype+".png").c_str());
}
