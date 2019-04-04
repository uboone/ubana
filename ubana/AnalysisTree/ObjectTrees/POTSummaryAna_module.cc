////////////////////////////////////////////////////////////////////////
// Class:       POTSummaryAna
// Plugin Type: analyzer (art v2_11_03)
// File:        POTSummaryAna_module.cc
//
// Generated at Mon Oct 22 18:29:07 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v3_03_01.
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

//include for the TFileService/ROOT
#include "art_root_io/TFileService.h"
#include "TTree.h"

//include the truth objects
#include "larcoreobj/SummaryData/POTSummary.h"

namespace ana {
  class POTSummaryAna;
}


class ana::POTSummaryAna : public art::EDAnalyzer {
public:
  explicit POTSummaryAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  POTSummaryAna(POTSummaryAna const &) = delete;
  POTSummaryAna(POTSummaryAna &&) = delete;
  POTSummaryAna & operator = (POTSummaryAna const &) = delete;
  POTSummaryAna & operator = (POTSummaryAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endSubRun(art::SubRun const&) override;

private:

  art::InputTag fInputTag;

  TTree* fPOTTree;

  unsigned int fRun;
  unsigned int fSubrun;

  double fTotpot;
  double fTotgoodpot;
  int fTotspills;
  int fGoodspills;
  
  void SetupPOTTree();
  void FillSubRunInfo(art::SubRun const&);
  void FillPOTSummaryInfo(sumdata::POTSummary const&);

};

void ana::POTSummaryAna::SetupPOTTree()
{
  fPOTTree->Branch("run",&fRun,"run/i");
  fPOTTree->Branch("subrun",&fSubrun,"subrun/i");

  fPOTTree->Branch("totpot",&fTotpot,"totpot/D");
  fPOTTree->Branch("totgoodpot",&fTotgoodpot,"totgoodpot/D");
  fPOTTree->Branch("totspills",&fTotspills,"totspills/I");
  fPOTTree->Branch("goodspills",&fGoodspills,"goodspills/I");
}

void ana::POTSummaryAna::FillSubRunInfo(art::SubRun const& s)
{
  fRun = s.run();
  fSubrun = s.subRun();
}

void ana::POTSummaryAna::FillPOTSummaryInfo(sumdata::POTSummary const& potsum)
{
  fTotpot = potsum.totpot;
  fTotgoodpot = potsum.totgoodpot;
  fTotspills = potsum.totspills;
  fGoodspills = potsum.goodspills;
}


ana::POTSummaryAna::POTSummaryAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
{
  fInputTag = p.get<art::InputTag>("InputTag");
}

void ana::POTSummaryAna::endSubRun(art::SubRun const& s)
{
  FillSubRunInfo(s);

  auto const& potsum_handle = s.getValidHandle< sumdata::POTSummary >(fInputTag);
  FillPOTSummaryInfo(*potsum_handle);

  fPOTTree->Fill();
}

void ana::POTSummaryAna::analyze(art::Event const & e)
{
}//end analyze event

void ana::POTSummaryAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fPOTTree = tfs->make<TTree>("pot_tree","POT Tree");
  SetupPOTTree();
}

DEFINE_ART_MODULE(ana::POTSummaryAna)
