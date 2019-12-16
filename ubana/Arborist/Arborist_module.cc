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

  TTree* eventweight_tree;
  int run;
  int subrun;
  int event;
  std::map<std::string, std::vector<double>> fmcweight;

};


Arborist::Arborist(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fEventWeightInputTag(p.get<art::InputTag>("EventWeightInputTag"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
void Arborist::beginJob()
{
 eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
 eventweight_tree->Branch("mcweight", "std::map<std::string, std::vector<double>>",&fmcweight);
 eventweight_tree->Branch("run", &run);
 eventweight_tree->Branch("subrun", &subrun);
 eventweight_tree->Branch("event", &event);
}

void Arborist::resetVariables()
{
  run = -1;
  subrun = -1;
  event = -1;
  fmcweight.clear();
}

void Arborist::FillEventWeights(art::Event const & e){
  art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw =
    e.getValidHandle<std::vector<evwgh::MCEventWeight>>(fEventWeightInputTag);

  std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
  if(ev_evw->size() > 1) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: eventweight has more than one entry\n";
  }
  fmcweight=weight_map;
 
}
void Arborist::analyze(art::Event const& e)
{
  resetVariables();
  run = e.run();
  subrun = e.subRun();
  event = e.event();

  FillEventWeights(e);
  eventweight_tree->Fill();

  // Implementation of required member function here.
}

DEFINE_ART_MODULE(Arborist)
