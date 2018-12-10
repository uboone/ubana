#include "LandauGaussianPlot.h"

double muonlandauwidth_mc = 0.06;
double muongaussianwidth_mc = 0.13;
double muonlandauwidth_data = 0.14;
double muongaussianwidth_data = 0.15;

// actually protons
//double muonlandauwidth_mc = 0.07;
//double muongaussianwidth_mc = 0.4;
//double muonlandauwidth_data = 0.16;
//double muongaussianwidth_data = 0.46;

void calculateScalingFactor(){

  // Setup Landau-gaussian
  double likelihoodbefore = 1;
  double likelihoodafter = 1;
 
  TF1 *langaus_mc = new TF1("langaus", landauGaussian, 0, 100, 4);
  langaus_mc->SetNpx(10000);
  langaus_mc->SetParameters(muonlandauwidth_mc, 1.7, 1.0, muongaussianwidth_mc);
  
  // now take the data widths
  TF1 *langaus_data = new TF1("langaus", landauGaussian, 0, 100, 4);
  langaus_data->SetNpx(10000);
  langaus_data->SetParameters(muonlandauwidth_data, 5.0, 1.0, muongaussianwidth_data);

  std::cout << "langaus_mc      : " << langaus_mc->Integral(0, langaus_mc->GetMaximumX()) << " : " << langaus_mc->Integral(langaus_mc->GetMaximumX(), 100) << std::endl;
  std::cout << "langaus_mc mean : " << langaus_mc->Mean(0,1000) << ", Eval@Mean: " << langaus_mc->Eval(langaus_mc->Mean(0,1000)) << std::endl;
  std::cout << "langaus_data    : " << langaus_data->Integral(0, langaus_data->GetMaximumX()) << " : " << langaus_data->Integral(langaus_data->GetMaximumX(), 100) << std::endl;
  std::cout << "langaus_data mean : " << langaus_data->Mean(0,1000) << ", Eval@Mean: " << langaus_data->Eval(langaus_data->Mean(0,1000)) << std::endl;

  std::cout << langaus_mc->Eval(5.2) << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  langaus_mc->SetLineColor(kRed);
  langaus_mc->GetXaxis()->SetRangeUser(0,6);
  langaus_mc->SetTitle(";dE/dx;");
  langaus_mc->Draw();

//  langaus_data->SetLineColor(kBlue);
//  langaus_data->Draw("same");
/*
  double randomPoint = 0;
  for (int i = 0; i < 100; i++){

    // find maximum
    double maxx = langaus_mc->GetMaximumX();

    // choose some points
    gRandom->SetSeed(0);
    double randomPoint = langaus_mc->GetRandom();

    // now calculate the fraction of the pdf encapsulated between the chosen
    // point and the MPV
    double i1;
    if (maxx > randomPoint)
      i1 = langaus_mc->Integral(randomPoint, maxx);
    else i1 = langaus_mc->Integral(maxx, randomPoint);
    likelihoodbefore *= langaus_mc->Eval(randomPoint);
  
    // this is the fraction of the integral below the max peak for the data
    if (i1 > 0.312697) continue;

    std::cout << "-------------" << std::endl;
    std::cout << "randomPoint                                    : " << randomPoint << std::endl;
    std::cout << "langaus_mc evaluated at randomPoint             : " << langaus_mc->Eval(randomPoint) << std::endl;
    std::cout << "langaus_mc integral between maxx and randomPoint: " << i1 << std::endl;
    std::cout << "langaus_mc cumulative likelihood                : " << likelihoodbefore << std::endl;

    double tester = 9999;
    double savedi = 0;
    double integral;
    double testPoint;
    double chosenPoint;
    for (double i = 0.005; i < 5.0; i = i+0.005){

      if (maxx > randomPoint){
        testPoint = maxx-i;
        integral = langaus_data->Integral(testPoint, maxx);
      }
      else {
        testPoint = maxx+i;
        integral = langaus_data->Integral(maxx, testPoint);
      }
      if (std::abs(integral-i1) < tester){
        tester = std::abs(integral-i1);
        chosenPoint = testPoint;
      }


    }

    likelihoodafter *= langaus_data->Eval(chosenPoint);

    std::cout << "langaus_data integral between maxx and chosen point: " << integral << std::endl;
    std::cout << "langaus_data chosen point                          : " << chosenPoint << std::endl;
    std::cout << "langaus_data evaluated at chosen point             : " << langaus_data->Eval(chosenPoint) << std::endl;
    std::cout << "langaus_data cumulative likelihood                 : " << likelihoodafter << std::endl;

  }

  std::cout << "[ScaleFactor] likelihoodbefore/likelihoodafter: " << likelihoodbefore/likelihoodafter << std::endl;
*/
}
