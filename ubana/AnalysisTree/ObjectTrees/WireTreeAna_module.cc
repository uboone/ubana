////////////////////////////////////////////////////////////////////////
// Class:       WireTreeAna
// Plugin Type: analyzer (art v2_11_03)
// File:        WireTreeAna_module.cc
//
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

//need the geometry
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()

//include the truth objects
#include "lardataobj/RecoBase/Wire.h"

namespace ana {
  class WireTreeAna;
}


class ana::WireTreeAna : public art::EDAnalyzer {
public:
  explicit WireTreeAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireTreeAna(WireTreeAna const &) = delete;
  WireTreeAna(WireTreeAna &&) = delete;
  WireTreeAna & operator = (WireTreeAna const &) = delete;
  WireTreeAna & operator = (WireTreeAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fInputTag;

  TTree* fWireTree;
  TTree* fROITree;

  unsigned int fRun;
  unsigned int fSubrun;
  unsigned int fEvent;
  
  //wire info
  unsigned int            fWireIndex;
  unsigned int            fChannel;
  int                     fView;
  unsigned int            fNSignal;
  unsigned int            fNROI;
  unsigned int            fWire;
  unsigned int            fPlane;
  unsigned int            fTPC;
  unsigned int            fCryostat;

  //roi info
  unsigned int            fROIIndex;
  unsigned int            fStart_tick;
  unsigned int            fEnd_tick;
  unsigned int            fSize;
  std::vector<float>      fWaveform;
  float                   fQTotal;
  float                   fQPeak;
  float                   fTimeAvg;
  float                   fTimePeak;
  float                   fTimeRMS;


  void SetupWireTree();
  void SetupROITree();

  void FillEventInfo(art::Event const&);
  void FillWireInfo(recob::Wire const&, geo::GeometryCore const&);
  void FillROIInfo(recob::Wire::RegionsOfInterest_t::datarange_t const&);

};

void ana::WireTreeAna::SetupWireTree()
{
  fWireTree->Branch("run",&fRun,"run/i");
  fWireTree->Branch("subrun",&fSubrun,"subrun/i");
  fWireTree->Branch("event",&fEvent,"event/i");

  fWireTree->Branch("wire_index",&fWireIndex,"wire_index/i");
  fWireTree->Branch("channel",&fChannel,"channel/i");
  fWireTree->Branch("view",&fView,"view/I");
  fWireTree->Branch("wire",&fWire,"wire/i");
  fWireTree->Branch("plane",&fPlane,"plane/i");
  fWireTree->Branch("tpc",&fTPC,"tpc/i");
  fWireTree->Branch("cryostat",&fCryostat,"cryostat/i");

  fWireTree->Branch("n_signal",&fNSignal,"n_signal/i");
  fWireTree->Branch("n_roi",&fNROI,"n_roi/i");
}

void ana::WireTreeAna::SetupROITree()
{
  fROITree->Branch("run",&fRun,"run/i");
  fROITree->Branch("subrun",&fSubrun,"subrun/i");
  fROITree->Branch("event",&fEvent,"event/i");

  fROITree->Branch("wire_index",&fWireIndex,"wire_index/i");
  fROITree->Branch("roi_index",&fROIIndex,"roi_index/i");

  fROITree->Branch("roi_start_tick",&fStart_tick,"roi_start_tick/i");
  fROITree->Branch("roi_end_tick",&fEnd_tick,"roi_end_tick/i");
  fROITree->Branch("roi_size",&fSize,"roi_size/i");
  fROITree->Branch("waveform",&fWaveform);

  fROITree->Branch("q_total",&fQTotal,"q_total/F");
  fROITree->Branch("q_peak",&fQPeak,"q_peak/F");
  fROITree->Branch("time_avg",&fTimeAvg,"time_avg/F");
  fROITree->Branch("time_peak",&fTimePeak,"time_peak/F");
  fROITree->Branch("time_rms",&fTimeRMS,"time_rms/F");

}

void ana::WireTreeAna::FillEventInfo(art::Event const& e)
{
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.event();
}

void ana::WireTreeAna::FillWireInfo(recob::Wire const& wire,
				    geo::GeometryCore const& geom)
{

  fChannel = wire.Channel();
  fView = wire.View();

  geo::WireID wid = geom.ChannelToWire(fChannel).at(0);

  fWire = wid.Wire;
  fPlane = wid.Plane;
  fTPC = wid.TPC;
  fCryostat = wid.Cryostat;

  fNSignal = wire.NSignal();
  fNROI = wire.SignalROI().n_ranges();

  for(size_t i_roi=0; i_roi<fNROI; ++i_roi){
    fROIIndex = i_roi;
    FillROIInfo( wire.SignalROI().range(i_roi) );
    fROITree->Fill();
  }

}

void ana::WireTreeAna::FillROIInfo(recob::Wire::RegionsOfInterest_t::datarange_t const& range)
{
  fStart_tick = range.begin_index();
  fEnd_tick = range.end_index();
  fSize = range.size();
  fWaveform = range.data();

  fQTotal=0;
  fQPeak=-999;
  fTimeAvg=0;
  fTimePeak=0;
  fTimeRMS=0;

  for(size_t i_t=0; i_t<fWaveform.size(); ++i_t){
    fQTotal+=fWaveform[i_t];
    fTimeAvg+= (i_t+fStart_tick)*fWaveform[i_t];
    if(fWaveform[i_t]>fQPeak){
      fQPeak=fWaveform[i_t];
      fTimePeak=i_t+fStart_tick;
    }
  }
  fTimeAvg = fTimeAvg/fQTotal;
  for(size_t i_t=0; i_t<fWaveform.size(); ++i_t){
    fTimeRMS += (i_t+fStart_tick-fTimePeak)*(i_t+fStart_tick-fTimePeak)*fWaveform[i_t];
  }
  fTimeRMS= std::sqrt(fTimeRMS/fQTotal);
}

ana::WireTreeAna::WireTreeAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
{
  fInputTag = p.get<art::InputTag>("InputTag");
}

void ana::WireTreeAna::analyze(art::Event const & e)
{
  //first fill the event-level information
  FillEventInfo(e);

  //grab the Hit collection
  auto const& wire_handle = e.getValidHandle< std::vector<recob::Wire> >(fInputTag);
  auto const& wire_vec(*wire_handle);

  //loop over hits and fill the tree
  for(size_t i_w=0; i_w<wire_vec.size(); ++i_w){
    fWireIndex = i_w;
    FillWireInfo(wire_vec[i_w],*(lar::providerFrom<geo::Geometry>()) );
    fWireTree->Fill();
  }//end loop over hits

}//end analyze event

void ana::WireTreeAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fWireTree = tfs->make<TTree>("wire_tree","Wire Tree");
  fROITree = tfs->make<TTree>("roi_tree","ROI Tree");
  SetupWireTree();
  SetupROITree();
}

DEFINE_ART_MODULE(ana::WireTreeAna)
