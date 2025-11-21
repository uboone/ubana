#include "NeutronReweighter.h"
#include "TH1D.h"
#include "TFile.h"

// Constructor
NeutronReweighter::NeutronReweighter() {
  //default configuration
  _useStandardUniverses=true;
  _numberUniverses=1000;
  
    _weights = new std::vector<double>(_numberUniverses,-999);
  _all_univ_xsec = new std::map<int, TH1D*>;
  _univ_calc = new std::vector<bool>(_numberUniverses,false);
  
    _randomNumbers  = nullptr;
      _nominalG4xsec  = nullptr;
      _scaledG4xsec   = nullptr;
      _curr_univ_xsec = nullptr;
    
//  this->SetXsecFileName("/exp/uboone/data/users/mhernan/neutron_inelastic_cross_section_mod.root");
this->SetXsecFileName("/pnfs/uboone/persistent/users/mhernan/neutron_inelastic_cross_section_mod.root"); 
 return;
}


void NeutronReweighter::SetRandomSeed(int seed){
  _randomSeed = seed;
    if (_randomNumbers) {
        delete _randomNumbers;
      }
  _randomNumbers=new std::vector<double>;
  std::normal_distribution<> gauss(0.0,1.0);
  for (int i(0); i<_numberUniverses; i++){
    std::mt19937 rng(_randomSeed+i);
    _randomNumbers->push_back(gauss(rng));
  }
}

// Set name of file with xsec histograms in
void NeutronReweighter::SetXsecFileName(std::string name){
  _xsecFileName = name;
  return;
}

// Configure the tool - needs to be run before calculating weights
void NeutronReweighter::Configure(){

  _weights = new std::vector<double>(_numberUniverses,-999);
  _randomNumbers=new std::vector<double>;
//  for (int i(0);i<_numberUniverses;i++){
//    _weights->at(i)=1.0;
//  }


  // save a random number per universe.  Use the fixed randomseed
  _randomSeed = 1337;
  std::normal_distribution<> gauss(0.0,1.0);
  for (int i(0); i<_numberUniverses; i++){
    std::mt19937 rng(_randomSeed+i);
    _randomNumbers->push_back(gauss(rng));
  }

//  TFile* file = TFile::Open(_xsecFileName.c_str(), "READ");

//TFile* file = TFile::Open("/exp/uboone/data/users/mhernan/neutron_inelastic_cross_section_mod.root","READ");
TFile* file = TFile::Open("/pnfs/uboone/persistent/users/mhernan/neutron_inelastic_cross_section_mod.root","READ");
  if (!file || file->IsZombie()) {
    //throw std::runtime_error("Failed to open xsec file: " + filename);
   std::cout << "ERROR" << std::endl;
  }
 
	TH1D* h_nom = dynamic_cast<TH1D*>(file->Get("xsec_histogram"));
        TH1D* h_scl = dynamic_cast<TH1D*>(file->Get("scaled_xsec_histogram_withError"));
 
  // Nominal GEANT4 cross section histogram
  _nominalG4xsec = (TH1D*)h_nom->Clone("nominalG4xsec_clone");
  if (!_nominalG4xsec) {
    std::cout << "ERROR" << std::endl;
  }
  
  // MiniCaptain data cross section histogram
  //_miniCaptainDataCV = dynamic_cast<TH1D*>(file->Get("minicaptain_histogram"));
  //if (!_nominalG4xsec) {
  //  //throw std::runtime_error("Failed to find xsec_histogram in file: " + filename);
  //  std::cout << "ERROR" << std::endl;
  //}

  //Scaled CV histogram includes error bar
    _scaledG4xsec  = (TH1D*)h_scl->Clone("scaledG4xsec_clone");
  if (!_scaledG4xsec) {
    std::cout << "ERROR" << std::endl;
  }
  
  _nominalG4xsec->SetDirectory(nullptr);
  _scaledG4xsec->SetDirectory(nullptr);

  _all_univ_xsec->clear();

  for (int i(0); i< _numberUniverses; i++){
	std::string name = "xsec_univ" + std::to_string(i);
    TH1D* tmp = (TH1D*)_scaledG4xsec->Clone(name.c_str());
    tmp->SetDirectory(nullptr);
    (*_all_univ_xsec)[i] = tmp;
  }

  file->Close();
  delete file;
  _configured=true;
  return;
}

// Call on a new event, to make re-use of event weights easy
void NeutronReweighter::ConfigureEvent(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process){
  if (!_configured){this->Configure();}
  for (int i(0); i<_numberUniverses; i++){
    _weights->at(i)=-999; // default weight, which should be reeeeally obvious if you grab it...
  }
  // Maybe I need to loop and fill?
  _event_KE = new std::vector<double>(KE);
  _event_lengths = new std::vector<double>(lengths);
  _event_end_process = new std::vector<std::string>(end_process);
  return;
}


double NeutronReweighter::GetNominalXsec(double KE){
  // Get the nominal xsec used in G4 simulation
  int bin = _nominalG4xsec->FindFixBin(KE);
  // Ensure KE is within valid histogram range
  if (bin <= 0 || bin > _nominalG4xsec->GetNbinsX()) {
    std::cout << "[DEBUG] KE=" << KE << " is outside valid xsec histogram range." << std::endl;
    return 0.0;
  }
  double nom_xsec = _nominalG4xsec->GetBinContent(bin);
  return nom_xsec;
}


double NeutronReweighter::GetCVcorrXsec(double KE){
    // Get the CV scakled xsec we use for the CV
  int bin = _scaledG4xsec->FindFixBin(KE);
  // Ensure KE is within valid histogram range
  if (bin <= 0 || bin > _scaledG4xsec->GetNbinsX()) {
    std::cout << "[DEBUG] KE=" << KE << " is outside valid xsec histogram range." << std::endl;
    return 0.0;
  }
  double cv_corr_xsec = _scaledG4xsec->GetBinContent(bin);
// DEBUG  std::cout << "cv_corr_xsec = " << cv_corr_xsec << std::endl;

  return cv_corr_xsec;
}


double NeutronReweighter::GetUniverseXsec(double KE, int univ){
  // returns cross section in the requested universe.  Universe -1 is the CV correction.
  if (univ==-1){
    double xsec = GetCVcorrXsec(KE);
    return xsec; 
  }
  if (!_univ_calc->at(univ)){
    this->CalculateUniverseXsecHist(univ);
  }
    
  if (!_curr_univ_xsec){
    std::cout << "xsec hist doesn't exist - calculating it" << std::endl;
    this->CalculateUniverseXsecHist(univ);
    //_curr_univ=univ;
  }
  int bin = _curr_univ_xsec->FindFixBin(KE);
  // Ensure KE is within valid histogram range
  if (bin <= 0 || bin > _nominalG4xsec->GetNbinsX()) {
    std::cout << "[DEBUG] KE=" << KE << " is outside valid xsec histogram range." << std::endl;
    return 0.0;
  }
  //this->SetUniverse(univ);
  double var_xsec = _curr_univ_xsec->GetBinContent(bin);
// DEBUG  std::cout << "var_xsec = " << var_xsec << std::endl;
  return var_xsec;
}

void NeutronReweighter::CalculateUniverseXsecHist(int univ){
  // Need to get nominal xsec, MiniCaptain data, etc, and calculate a new histogram for this universe
  if (_univ_calc->at(univ)){
    return;
  }
  //_curr_univ=univ;
  _curr_univ_xsec = _all_univ_xsec->at(univ);
// DEBUG  std::cout << univ << _randomNumbers->at(univ) << std::endl;
  
  for (int i(0); i<_curr_univ_xsec->GetNbinsX(); i++){
    // sigma for this bin - depends on CV scaled xsec and minicaptain data
    double sigma=_scaledG4xsec->GetBinError(i);
    double cv_xsec = _scaledG4xsec->GetBinContent(i);
    double var_xsec = cv_xsec + _randomNumbers->at(univ) * sigma;
    _curr_univ_xsec->SetBinContent(i,var_xsec);
  }
  _univ_calc->at(univ) = true;
  return;
}

// Single weight calculator.  This is the workhorse of the whole class.  Everything else is wrapping...
double NeutronReweighter::CalculateSegmentWeight(double seg_KE, double seg_length, std::string seg_end_process, double nom_xsec, double var_xsec){
  
  // change to xsec
  const double delta_sigma = var_xsec - nom_xsec;
  
  // --- Integrate along the path without a vanishing micro-step ---
  constexpr double tiny = 1e-12;   // cm
  constexpr double epsP = 1e-300;  // guard against 0 denominator

  const int N_full = static_cast<int>(std::floor(seg_length / _segmentStepSize + tiny));
  double used  = static_cast<double>(N_full) * _segmentStepSize;
  double rem   = seg_length - used;
  if (rem < tiny) rem = 0.0;

  // Do survival over (N_full - 1) full steps; apply interaction (or survival) on the final step
  const int    n_surv_steps = std::max(0, N_full - 1);
  const double last_step    = (N_full > 0 ? _segmentStepSize : 0.0) + rem;

  double segment_weight = 1.0;

  // Survival over the first chunk (vectorized)
  if (n_surv_steps > 0) {
    const double f_step = std::exp(-delta_sigma * _segmentStepSize / _inverse_number_density);
    if (std::isfinite(f_step)) segment_weight *= std::pow(f_step, n_surv_steps);
  }

  // Final step
  if (last_step > tiny) {
    const bool interacted = (seg_end_process == "neutronInelastic");

    if (interacted) {
      // Numerically-stable interaction probabilities
      const double x_nom =  nom_xsec   * last_step / _inverse_number_density;
      const double x_var =  var_xsec   * last_step / _inverse_number_density;
      const double P_nom = -std::expm1(-x_nom);   // 1 - exp(-x_nom)
      const double P_var = -std::expm1(-x_var);

      if (!(P_nom == 0.0 && P_var == 0.0)) {
        const double denom = (P_nom == 0.0 ? epsP : P_nom);
        const double w_int = P_var / denom;
        if (std::isfinite(w_int)) segment_weight *= w_int;
        // else: neutralize this step instead of poisoning the weight
      }
      // If both P's are exactly zero, treat as neutral (no ratio info).
    } else {
      const double f_last = std::exp(-delta_sigma * last_step / _inverse_number_density);
      if (std::isfinite(f_last)) segment_weight *= f_last;
    }
  }
// DEBUG  std::cout << "KE = " << seg_KE << ", seg weight = " << segment_weight << std::endl;
  return segment_weight;
}


double NeutronReweighter::CalculateEventWeight(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process, int univ){
  double event_weight=1.0;
  // Check the vectors are all the same length
  //_curr_univ=univ;
  int num_neutrons = KE.size();
  if (lengths.size()!= abs(num_neutrons)){
    std::cout << "ERROR, vectors of KE, lengths, and end processes need to match!  lengths is not the same length as KE!" << std::endl;
    std::cout << "Event weight will be set to 1" << std::endl;
    return 1;
  }
  if (end_process.size()!= abs(num_neutrons)){
    std::cout << "ERROR, vectors of KE, lengths, and end processes need to match!  end_process is not the same length as KE!" << std::endl;
    std::cout << "Event weight will be set to 1" << std::endl;
    return 1;
  }
  if (end_process.size()!=lengths.size()){
    std::cout << "ERROR, vectors of KE, lengths, and end processes need to match!  end_process is not the same length as lengths!" << std::endl;
    std::cout << "Event weight will be set to 1" << std::endl;
    return 1;
  }

  // now loop through segments and multiply event weight by segment weights
  for (int i(0); i<num_neutrons; i++){
    double nom_xsec = GetNominalXsec(KE.at(i));
    double var_xsec = GetUniverseXsec(KE.at(i), univ);
    double seg_weight = CalculateSegmentWeight(KE.at(i), lengths.at(i), end_process.at(i), nom_xsec, var_xsec);
    event_weight*=seg_weight;
  }
  _weights->at(univ) = event_weight;
  return event_weight;
}

// version for anybody to use
double NeutronReweighter::GetEventWeight(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process, int univ){
  if (_weights->size()< abs(univ)){
    std::cout << "ERROR: universe requested is larger than the number calculated" << std::endl;
  }
  if (_weights->at(univ)>-1){
    return _weights->at(univ);
  }
  else{
    double weight = CalculateEventWeight(KE, lengths, end_process, univ);
    _weights->at(univ)=weight;
    return weight;
  }
}

// Version for configured event
double NeutronReweighter::GetEventWeight(int univ){
  if (_weights->at(univ)>-1){
    return _weights->at(univ);
  }
  else{
    double weight = GetEventWeight(*_event_KE, *_event_lengths, *_event_end_process, univ);
    _weights->at(univ)=weight;
    return weight;
  }
}


void NeutronReweighter::CalculateAllWeights(std::vector<double> KE, std::vector<double> lengths, std::vector<std::string> end_process){
  for (int univ(0); univ<_numberUniverses; univ++){
    double weight_univ = GetEventWeight(KE, lengths, end_process, univ);
    _weights->at(univ) = weight_univ;
  }
}

void NeutronReweighter::CalculateAllWeights(){
  this->CalculateAllWeights(*_event_KE, *_event_lengths, *_event_end_process);
}

// Get weight for a given universe
//double NeutronReweighter::GetUniverseWeight(int univ){
//  if (_weights->at(univ)>-1){
//    return _weights->at(univ);
//  }
//  else{
//    weight = this->CalculateEventWeight(univ);
//  }
//}

NeutronReweighter::~NeutronReweighter(){
    if (_weights)       delete _weights;
      if (_randomNumbers) delete _randomNumbers;

      if (_nominalG4xsec) delete _nominalG4xsec;
      if (_scaledG4xsec)  delete _scaledG4xsec;

      if (_all_univ_xsec) {
        for (auto &kv : *_all_univ_xsec) {
          if (kv.second) {
            kv.second->SetDirectory(nullptr);
            delete kv.second;
          }
        }
        delete _all_univ_xsec;
      }

      if (_univ_calc) delete _univ_calc;
    }
