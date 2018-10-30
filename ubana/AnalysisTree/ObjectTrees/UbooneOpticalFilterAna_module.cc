////////////////////////////////////////////////////////////////////////
// Class:       UbooneOpticalFilterAna
// Plugin Type: analyzer (art v2_11_03)
// File:        UbooneOpticalFilterAna_module.cc
//
// Generated at Mon Oct 22 13:28:33 2018 by Wesley Ketchum using cetskelgen
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

//include the UbooneOpticalFilter object
#include "ubobj/Optical/UbooneOpticalFilter.h"

namespace ana {
  class UbooneOpticalFilterAna;
}


class ana::UbooneOpticalFilterAna : public art::EDAnalyzer {
public:
  explicit UbooneOpticalFilterAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UbooneOpticalFilterAna(UbooneOpticalFilterAna const &) = delete;
  UbooneOpticalFilterAna(UbooneOpticalFilterAna &&) = delete;
  UbooneOpticalFilterAna & operator = (UbooneOpticalFilterAna const &) = delete;
  UbooneOpticalFilterAna & operator = (UbooneOpticalFilterAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fInputTag;

  TTree* fAnaTree;

  unsigned int fRun;
  unsigned int fEvent;
  unsigned int fTime_s;
  unsigned int fTime_ns;

  float fPE_beam;
  float fPE_veto;
  float fMaxFrac;
  float fPE_beam_total;
  float fPE_veto_total;
};


ana::UbooneOpticalFilterAna::UbooneOpticalFilterAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fInputTag = p.get<art::InputTag>("InputTag");
}

void ana::UbooneOpticalFilterAna::analyze(art::Event const & e)
{
  //fill in event info
  fRun = e.run();
  fEvent = e.event();
  fTime_s = e.time().timeHigh();
  fTime_ns = e.time().timeLow();

  //fill optical filter info
  auto const& opfilter_handle = e.getValidHandle<uboone::UbooneOpticalFilter>(fInputTag);
  fPE_beam = opfilter_handle->PE_Beam();
  fPE_veto = opfilter_handle->PE_Veto();
  fMaxFrac = opfilter_handle->PMT_MaxFraction();
  fPE_beam_total = opfilter_handle->PE_Beam_Total();
  fPE_veto_total = opfilter_handle->PE_Veto_Total();

  fAnaTree->Fill();

}

void ana::UbooneOpticalFilterAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fAnaTree = tfs->make<TTree>("opfilter_anatree","UbooneOpticalFilter Ana Tree");

  fAnaTree->Branch("run",&fRun,"run/i");
  fAnaTree->Branch("event",&fEvent,"event/i");
  fAnaTree->Branch("time_s",&fTime_s,"time_s/i");
  fAnaTree->Branch("time_ns",&fTime_ns,"time_ns/i");

  fAnaTree->Branch("pe_beam",&fPE_beam,"pe_beam/F");
  fAnaTree->Branch("pe_veto",&fPE_veto,"pe_veto/F");
  fAnaTree->Branch("maxfrac",&fMaxFrac,"maxfrac/F");
  fAnaTree->Branch("pe_beam_total",&fPE_beam_total,"pe_beam_total/F");
  fAnaTree->Branch("pe_veto_total",&fPE_veto_total,"pe_veto_total/F");
}

DEFINE_ART_MODULE(ana::UbooneOpticalFilterAna)
