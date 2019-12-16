////////////////////////////////////////////////////////////////////////
// Class:       Arborist
// Plugin Type: analyzer (art v3_01_02)
// File:        Arborist_module.cc
//
// Generated at Mon Mar 18 11:32:48 2019 by Gray Yarbrough using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
// ART includes
#include "art/Framework/Services/Optional/TFileService.h" // used for ROOT file
#include "canvas/Persistency/Common/FindManyP.h" // used for assns

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include "TTree.h"


class Arborist;


class Arborist : public art::EDAnalyzer {
public:
  explicit Arborist(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Arborist(Arborist const&) = delete;
  Arborist(Arborist&&) = delete;
  Arborist& operator=(Arborist const&) = delete;
  Arborist& operator=(Arborist&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void FillEventWeights(art::Event const & e);
  void resetVariables();


private:
  art::ServiceHandle< art::TFileService > tfs;

  art::InputTag fEventWeightInputTag;
  std::string fSplineBugFixLabel;
  std::string fTunedCentralValueLabel;
  std::string fLEESignalLabel;

  TTree* eventweight_tree;
  int run;
  int subrun;
  int event;
  double fSplineBugFixWeight;
  double fTunedCentralValueWeight;
  double fCombinedCentralValueWeight;
  double fLEESignalWeight;
  std::map<std::string, std::vector<double>> fEventWeightMap;

};


Arborist::Arborist(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fEventWeightInputTag(p.get<art::InputTag>("EventWeightInputTag")),
  fSplineBugFixLabel(p.get<std::string>("SplineBugFixLabel")),
  fTunedCentralValueLabel(p.get<std::string>("TunedCentralValueLabel")),
  fLEESignalLabel(p.get<std::string>("LEESignalLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
void Arborist::beginJob()
{
 eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
 eventweight_tree->Branch("run", &run);
 eventweight_tree->Branch("subrun", &subrun);
 eventweight_tree->Branch("event", &event);
 eventweight_tree->Branch("spline_weight", &fSplineBugFixWeight);
 eventweight_tree->Branch("ub_tune_weight", &fTunedCentralValueWeight);
 eventweight_tree->Branch("comb_cv_weight", &fCombinedCentralValueWeight);
 eventweight_tree->Branch("lee_weight", &fLEESignalWeight);
 eventweight_tree->Branch("mcweight", "std::map<std::string, std::vector<double>>", &fEventWeightMap);
}

void Arborist::resetVariables()
{
  run = -1;
  subrun = -1;
  event = -1;
  fSplineBugFixWeight = -1.;
  fTunedCentralValueWeight = -1.;
  fCombinedCentralValueWeight = -1.;
  fLEESignalWeight = -1.;
  fEventWeightMap.clear();
}

void Arborist::FillEventWeights(art::Event const & e){
  art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw =
    e.getValidHandle<std::vector<evwgh::MCEventWeight>>(fEventWeightInputTag);

  std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
  if(ev_evw->size() > 1) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: eventweight has more than one entry\n";
  }

  for ( auto evweight_it = weight_map.begin(); evweight_it != weight_map.end(); evweight_it++ ) {
    std::cout << evweight_it->first << std::endl;
  }

  fEventWeightMap = weight_map;
  
  if ( weight_map.find(fSplineBugFixLabel) != weight_map.end() ) {
    fSplineBugFixWeight = weight_map.find(fSplineBugFixLabel)->second[0];
    fEventWeightMap.erase(fSplineBugFixLabel);
  }
  if ( weight_map.find(fTunedCentralValueLabel) != weight_map.end() ) {
    fTunedCentralValueWeight = weight_map.find(fTunedCentralValueLabel)->second[0];
    fEventWeightMap.erase(fTunedCentralValueLabel);
  }
  if ( weight_map.find(fLEESignalLabel) != weight_map.end() ) {
    fLEESignalWeight = weight_map.find(fLEESignalLabel)->second[0];
    fEventWeightMap.erase(fLEESignalLabel);
  }

  if ( weight_map.find(fSplineBugFixLabel) != weight_map.end() && weight_map.find(fTunedCentralValueLabel) != weight_map.end() )
    fCombinedCentralValueWeight = fSplineBugFixWeight * fTunedCentralValueWeight;
  
}
void Arborist::analyze(art::Event const& e)
{
  resetVariables();

  run = e.run();
  subrun = e.subRun();
  event = e.event();

  FillEventWeights(e);
  
  eventweight_tree->Fill();
}

DEFINE_ART_MODULE(Arborist)
