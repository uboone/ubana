////////////////////////////////////////////////////////////////////////
// Class:       BlipVertexAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        BlipVertexAnalysis_module.cc
//
// Generated at Tue Feb  6 00:18:35 2024 by Keng Lin using cetskelgen
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

// ROOT includes
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

// extra includes, trying to fix
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"

// blip reco includes
//#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// for MCParticle.E() and such
#include "nusimdata/SimulationBase/MCTrajectory.h"
// for associations
#include "canvas/Persistency/Common/FindManyP.h" // find associations as pointers
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "lardataobj/RecoBase/Track.h"
// access hits and c++ stuff
#include "lardataobj/RecoBase/Hit.h"
#include <vector> // to store multiple things per event
#include <string> // for product names
#include <set>
#include <string>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"

// this line good for "new" larsoft:
#include "art/Persistency/Common/PtrMaker.h"
// for outdated versions (i.e MCC8) use this line:
//#include "lardata/Utilities/PtrMaker.h"

// additional framework includes
//#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Optional/TFileService.h"

// pot summary
// from /uboone/app/users/afurmans/muonScatteringAnalysis/
//       larsoft/srcs/uboonecode/uboone/AnalysisTree/
//       AnalysisTree_module.cc
#include "larcoreobj/SummaryData/POTSummary.h"


class BlipVertexAnalysis;


class BlipVertexAnalysis : public art::EDAnalyzer {
public:
  explicit BlipVertexAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BlipVertexAnalysis(BlipVertexAnalysis const&) = delete;
  BlipVertexAnalysis(BlipVertexAnalysis&&) = delete;
  BlipVertexAnalysis& operator=(BlipVertexAnalysis const&) = delete;
  BlipVertexAnalysis& operator=(BlipVertexAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void endSubRun(const art::SubRun& sr);


private:

  // Declare member data here.
  // blip
  std::string m_pandoraLabel;
  std::string m_blipLabel;

  TTree * f_output_tree; // output tree

};


BlipVertexAnalysis::BlipVertexAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
//  fhicl::ParameterSet pset_blipalg = p.get<fhicl::ParameterSet>("BlipAlg");
//  fBlipAlg = new blip::BlipRecoAlg( pset_blipalg );
}

void BlipVertexAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  f_output_tree->Fill();
}

void BlipVertexAnalysis::endSubRun(const art::SubRun& sr)
{

}

void BlipVertexAnalysis::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;


  f_output_tree = tfs->make<TTree>("t","t");
//  f_output_tree->Branch("run",&run);
//  f_output_tree->Branch("evt",&evt);

}

void BlipVertexAnalysis::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(BlipVertexAnalysis)
