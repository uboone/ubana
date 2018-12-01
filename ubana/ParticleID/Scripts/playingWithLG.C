/**
 * Simple script to check how convolution of a landau with a 
 * gaussian of different widths (as done in LandauGaussian.h)
 * affects the MPV of the convoluted function.
 *
 * This was written to try and determine whether there needed to be an 
 * additional shift for the LG function with respect to theory.
 */

#include "../Algorithms/LandauGaussian.h"

Double_t landauGaussian();

void playingWithLG(){

  double landau_width = 0.19;

  TF1* landau_gaussian = new TF1("landau_gaussian", landauGaussian, 0, 100, 4);
  landau_gaussian->SetNpx(1000);

  TH1D* lan_gaus_mpv = new TH1D("lan_gaus_mpv", "", 1000, 0, 1.0);

  TFile* output = new TFile("output.root", "recreate");
  output->cd();

  int bin_no = 1;
  for (double gaussian_width = 0.001; gaussian_width < 1.0; gaussian_width += 0.001){

    landau_gaussian->SetParameters(landau_width, 10, 1.0, gaussian_width);

    landau_gaussian->SetName(Form("langaus_gauswidth%f", gaussian_width));
    landau_gaussian->Write();

    double max_x = landau_gaussian->GetMaximumX();
    
    lan_gaus_mpv->SetBinContent(bin_no,max_x);

    bin_no++;
  }

  lan_gaus_mpv->Write();
}
