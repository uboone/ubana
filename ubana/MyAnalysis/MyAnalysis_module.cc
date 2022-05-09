////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        MyAnalysis_module.cc
//
// Generated at Sun May  8 19:56:01 2022 by Erin Yandel using cetskelgen
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
#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"


class MyAnalysis;


class MyAnalysis : public art::EDAnalyzer {
public:
  explicit MyAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyAnalysis(MyAnalysis const&) = delete;
  MyAnalysis(MyAnalysis&&) = delete;
  MyAnalysis& operator=(MyAnalysis const&) = delete;
  MyAnalysis& operator=(MyAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  TH1F* fHist; //!< Output histogram
  // Declare member data here.

};


MyAnalysis::MyAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  float maxEnergy = p.get<float>("MaxNuEnergy", 3.0);
  art::ServiceHandle<art::TFileService> tfs;
  fHist = tfs->make<TH1F>("enu", ";E_{#nu} (GeV);Events", 100, 0, maxEnergy);

}

void MyAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<simb::MCTruth> > mctruths;
  e.getByLabel("generator", mctruths);
  for (auto const& truth : *mctruths) {
    const simb::MCNeutrino& mcnu = truth.GetNeutrino();
    const simb::MCParticle& nu = mcnu.Nu();
    float enu = nu.E();
    fHist->Fill(enu);
  }

}

DEFINE_ART_MODULE(MyAnalysis)
