/**
 * intended to separate out protons from muons without using likelihoods 
 * so that we don't bias the widths.
 *
 * run with root -l pid_onbeam.root plotProtons.C 
 *
 */ 

// cpp
#include <algorithm>
#include <vector>

// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

// local
#include "../Algorithms/LandauGaussian.h"

bool isInFV(TVector3 vec){

  if (vec.X() > 10 && vec.X() < 246
      && vec.Y() > -105 && vec.Y() < 105
      && vec.Z() > 10 && vec.Z() < 1030){
    return true;
  }
  else return false;

}

Double_t fitFunction(Double_t *x, Double_t *p){
 
  return landauGaussian(x, p) + landauGaussian(x,&p[4]);
  
};

void plotProtons(){

  struct dataFrame {

    double strklen;
    std::vector<double> strkdedx;
    std::vector<double> strkresrg;
    double strk_pfwd;
    double strk_pbwd;
    TVector3 strk_start;
    TVector3 strk_end;

  };

  int ntracks_demanded = 2;
  double low_length = 65;
  double high_length = 70;

  int run;
  int sub_run;
  int event;
  double track_length;
  double track_start_x;
  double track_start_y;
  double track_start_z;
  double track_end_x;
  double track_end_y;
  double track_end_z;
  std::vector<double>* track_dEdx_perhit_y = 0;
  std::vector<double>* track_resrange_perhit_y = 0;
  std::vector<double>* track_likelihood_fwd_p = 0;
  std::vector<double>* track_likelihood_bwd_p = 0;

  TTree* tree = (TTree*)_file0->Get("pidvalid/pidTree");

  tree->SetBranchAddress("track_length", &track_length);
  tree->SetBranchAddress("run", &run); 
  tree->SetBranchAddress("sub_run", &sub_run);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_perhit_y);
  tree->SetBranchAddress("track_resrange_perhit_y", &track_resrange_perhit_y);
  tree->SetBranchAddress("track_likelihood_fwd_p", &track_likelihood_fwd_p);
  tree->SetBranchAddress("track_likelihood_bwd_p", &track_likelihood_bwd_p);
  tree->SetBranchAddress("track_start_x", &track_start_x);
  tree->SetBranchAddress("track_start_y", &track_start_y);
  tree->SetBranchAddress("track_start_z", &track_start_z);
  tree->SetBranchAddress("track_end_x", &track_end_x);
  tree->SetBranchAddress("track_end_y", &track_end_y);
  tree->SetBranchAddress("track_end_z", &track_end_z);

  // do not change
  int prev_run = 0;
  int prev_sub_run = 0;
  int prev_event = 0;
  int ntracks_in_event = 1;

  std::vector<dataFrame> dataFrames;

  TH2D* h_short = new TH2D("h_short", ";Residual Range (cm); dE/dx (MeV/cm)", 100, 0, 200, 50, 0, 10);
  TH2D* h_short_likelihood_identified = new TH2D("h_short_likelihood_identified", ";Residual Range (cm); dE/dx (MeV/cm)", 100, 0, 200, 50, 0, 10);
  TH1D* h_short_likelihood_identified_dedx = new TH1D("h_short_likelihood_identified_dedx", ";Residual Range (cm); dE/dx (MeV/cm)", 50, 0, 10);
  TH1D* h_short_dedx = new TH1D("h_short_dedx", ";dE/dx (MeV/cm);", 50, 0, 10);
  TH2D* h_short_dedx_angle = new TH2D("h_short_dedx_angle", "", 50, 0, 10, 50, 0, 3.15);
  TH1D* h_short_dedx_anglecut = new TH1D("h_short_dedx_anglecut", ";dE/dx (MeV/cm);", 50, 0, 10);
  TH2D* h_long = new TH2D("h_long", ";Residual Range (cm); dE/dx (MeV/cm)", 100, 0, 200, 50, 0, 10);

  for (int i = 0; i < tree->GetEntries(); i++){

    // get entry in TTree
    tree->GetEntry(i);

    // define data frame struct
    dataFrame df;

    // if this event is the same event as the last event then keep
    // counting tracks...
    if (run == prev_run && sub_run == prev_sub_run && event == prev_event){

      ntracks_in_event++;

      df.strklen = track_length;
      df.strkdedx = *track_dEdx_perhit_y;
      df.strkresrg = *track_resrange_perhit_y;
      df.strk_pfwd = track_likelihood_fwd_p->at(2);
      df.strk_pbwd = track_likelihood_bwd_p->at(2);

      TVector3 vect(track_start_x, track_start_y, track_start_z);
      df.strk_start = vect;
      TVector3 vect2(track_end_x, track_end_y, track_end_z);
      df.strk_end = vect2;

      dataFrames.push_back(df);

      prev_run = run;
      prev_sub_run = sub_run;
      prev_event = event;
    }
    // if it's not then first take the dataFrame and fill histogram
    // with track dE/dx values for the shorter track of the pair
    else {

      if (i !=0 && ntracks_in_event == ntracks_demanded){

        // sort dataFrames based on track length
        // post sorting, the first entry should be the shortest
        // track, the last entry should be the longest track
        std::sort(dataFrames.begin(), dataFrames.end(),
            [](const auto &i1, const auto &i2) { return i1.strklen < i2.strklen; } );

       if (dataFrames.at(ntracks_demanded - 1).strklen > 200 && dataFrames.at(0).strklen > 5 && isInFV(dataFrames.at(0).strk_start) == true && isInFV(dataFrames.at(0).strk_end) == true){

          // loop shortest track
          for (size_t j = 0; j < dataFrames.at(0).strkdedx.size(); j++){

            // dE/dx v resrg
            h_short->Fill(dataFrames.at(0).strkresrg.at(j), dataFrames.at(0).strkdedx.at(j));

            // take region between 40 and 45 residual range (to try and fit)
            if (dataFrames.at(0).strkresrg.at(j) > low_length && dataFrames.at(0).strkresrg.at(j) < high_length){

              // 1d representation of above
              h_short_dedx->Fill(dataFrames.at(0).strkdedx.at(j));

              // calulate angle between two tracks
              double angle = (dataFrames.at(0).strk_end - dataFrames.at(0).strk_start).Angle(dataFrames.at(1).strk_end - dataFrames.at(1).strk_start);

              // 2d plot dE/dx v angle
              h_short_dedx_angle->Fill(dataFrames.at(0).strkdedx.at(j), angle);

              // make a generous cut on angle to try and access only protons
              if (angle < 2 && angle > 1)
                h_short_dedx_anglecut->Fill(dataFrames.at(0).strkdedx.at(j));

            }

            // protons identified as such by the likelihood algo
            if(std::max(dataFrames.at(0).strk_pfwd, dataFrames.at(0).strk_pbwd) > 0.18){
              h_short_likelihood_identified->Fill(dataFrames.at(0).strkresrg.at(j), dataFrames.at(0).strkdedx.at(j));
            
                if (dataFrames.at(0).strkresrg.at(j) > low_length
                    && dataFrames.at(0).strkresrg.at(j) < high_length){

                  h_short_likelihood_identified_dedx->Fill(dataFrames.at(0).strkdedx.at(j));

                }


            }

          }

          // loop longest track
          for (size_t j = 0; j < dataFrames.at(1).strkdedx.size(); j++){

            h_long->Fill(dataFrames.at(1).strkresrg.at(j), dataFrames.at(1).strkdedx.at(j));

          }

        }

      }

      // now reset the counters, clear the dataFrame and 
      // begin counting again
      ntracks_in_event = 1;
      prev_run = run;
      prev_sub_run = sub_run;
      prev_event = event;

      dataFrames.clear();

      df.strklen = track_length;
      df.strkdedx = *track_dEdx_perhit_y;
      df.strkresrg = *track_resrange_perhit_y;
      df.strk_pfwd = track_likelihood_fwd_p->at(2);
      df.strk_pbwd = track_likelihood_bwd_p->at(2);
      TVector3 vect(track_start_x, track_start_y, track_start_z);
      df.strk_start = vect;
      TVector3 vect2(track_end_x, track_end_y, track_end_z);
      df.strk_end = vect2;

      dataFrames.push_back(df);

    }

  }

  TCanvas *c_short = new TCanvas();
  h_short->Draw("colz");
  h_short->SetStats(kFALSE);
  c_short->SaveAs("twoTrack_shortTrack_dedxresrg.png");

  TCanvas *c_short_likelihood_identified = new TCanvas();
  h_short_likelihood_identified->Draw("colz");
  c_short_likelihood_identified->SaveAs("twoTrack_shortTrack_likelihoodidentified_dedxresrg.png");


  TCanvas *c_short_likelihood_identified_dedx = new TCanvas();
  h_short_likelihood_identified_dedx->Draw("colz");

  TF1* lg = new TF1("lg", landauGaussian, 0, 10, 4);
  lg->SetParameters(0.13, 1.7, 100, 0.25);
  lg->FixParameter(0,0.13);
  lg->SetParLimits(3,0, 5);
  lg->SetParLimits(1,2.0, 4.0);
  h_short_likelihood_identified_dedx->Fit(lg, "", "", 0.0, 10);
  c_short_likelihood_identified_dedx->SaveAs("twoTrack_shortTrack_likelihoodidentified_dedx.png");
  
  TCanvas *c_short_dedx = new TCanvas();
  h_short_dedx->Draw();

  // now take this histogram and fit the sum of two landau-gaussians
  TF1* fitFcn = new TF1("fitFcn", fitFunction, 0, 10, 8);
  fitFcn->SetParameters(0.09, 1.7, 200, 0.2, 0.09, 3.5, 100, 0.175);
  fitFcn->FixParameter(0, 0.09);
  fitFcn->FixParameter(3, 0.2);
  fitFcn->FixParameter(4, 0.13);
  fitFcn->SetParLimits(5,1.7, 10.0);
  h_short_dedx->Fit(fitFcn, "M", "", 1.2, 20);
  //h_short_dedx->Fit(langaus_p, "+");

  c_short_dedx->SaveAs("twoTrack_shortTrack_dedx.png");

  TCanvas *c_short_dedx_angle = new TCanvas();
  h_short_dedx_angle->Draw("colz");
  c_short_dedx_angle->SaveAs("twoTrack_shortTrack_dedx_angle.png");

  TCanvas *c_short_dedx_anglecut = new TCanvas();
  h_short_dedx_anglecut->Draw();
  c_short_dedx_anglecut->SaveAs("twoTrack_shortTrack_dedx_anglecut.png");

  TCanvas *c_long = new TCanvas();
  h_long->Draw("colz");
  h_long->SetStats(kFALSE);
  c_long->SaveAs("twoTrack_longTrack_dedxresrg.png");

}

int main(){

  plotProtons();

}
