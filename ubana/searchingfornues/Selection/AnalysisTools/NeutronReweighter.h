#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <map>
#include <random>
#include <numeric>
#include <functional>

class NeutronReweighter {
  
  public:
  NeutronReweighter();
  ~NeutronReweighter();
//  void UseStandardUniverses(bool true); // I think not needed
  void CalculateAllWeights(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> endprocess);
  void CalculateAllWeights();
  double CalculateSingleWeight(int univ);
//  double GetUniverseWeight(int univ);

  double GetNominalXsec(double KE);
  double GetUniverseXsec(double KE, int univ);
  double GetCVcorrXsec(double KE);
  
  double CalculateSegmentWeight(double seg_KE, double seg_length, std::string seg_end_process, double curr_xsec, double var_xsec);

  double GetEventWeight(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process, int univ);
  double GetEventWeight(int univ);

  void SetXsecFileName(std::string name);
  void Configure();
  void ConfigureEvent(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process);
//  void 
  void SetRandomSeed(int seed);


  private:
  std::string _xsecFileName;
  TH1D* _nominalG4xsec;
  TH1D* _curr_univ_xsec;
//  int _curr_univ;
//  TH1D* _miniCaptainDataCV;
  TH1D* _scaledG4xsec; // Error bars will be 1-sigma bounds!
  std::map<int,TH1D*> *_all_univ_xsec;
  std::vector<double>* _weights;
  int _numberUniverses;
  std::vector<bool>* _univ_calc;

  bool _configured = false;
  int _randomSeed = 1337;
  std::mt19937 _rng{std::mt19937::result_type(1337)};

  std::vector<double>* _event_KE;
  std::vector<double>* _event_lengths;
  std::vector<std::string>* _event_end_process;

  double CalculateEventWeight(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process, int univ);

  void CalculateUniverseXsecHist(int univ);

  // Global constant for inverse nucleon number density (cm^3 / nucleon)  
  const double _inverse_number_density = 39.95 / (1.394 * 6.023e23);

  bool _useStandardUniverses;

  const double _segmentStepSize=0.1;

  std::vector<double>* _randomNumbers;
};
