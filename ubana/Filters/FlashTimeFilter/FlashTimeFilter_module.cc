////////////////////////////////////////////////////////////////////////
// Class:       FlashTimeFilter
// Plugin Type: filter (art v3_01_02)
// File:        FlashTimeFilter_module.cc
//
// Generated at Tue Aug  6 08:16:01 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//#include "art_root_io/TFileService.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TH1.h"
//#include "TTree.h"


class FlashTimeFilter;


class FlashTimeFilter : public art::EDFilter {
public:
  explicit FlashTimeFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.

  FlashTimeFilter(FlashTimeFilter const&) = delete;
  FlashTimeFilter(FlashTimeFilter&&) = delete;
  FlashTimeFilter& operator=(FlashTimeFilter const&) = delete;
  FlashTimeFilter& operator=(FlashTimeFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  art::ServiceHandle<art::TFileService> tfs;
  TH1F *wfstarta;
  TH1F *wfstartb;
  //  detinfo::DetectorClocks const&  detectorClocks;

  // Declare member data here.
  std::string fInputModuleName;
  int fWFcount;
  int fWFTimeStart;
  int fWFTimeEnd;

  bool endSubRun(art::SubRun &subrun) override;
};


FlashTimeFilter::FlashTimeFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  fInputModuleName = p.get< std::string >("InputModule","pmtreadout:OpdetCosmicHighGain" );
  fWFcount = p.get<int>("WFcount",1);
  fWFTimeStart = p.get<int>("WFTimeStart",4300.0);
  fWFTimeEnd = p.get<int>("WFTimeEnd",6400.0);


  // Create histograms                                                                                                      
  wfstarta = tfs->make<TH1F>("wfstarta","wfstarta",500,-3500.,6500.);  // trig time in us                                   
  wfstarta->GetXaxis()->SetTitle("Selected Waveform Start (us)"); 
 
  wfstartb = tfs->make<TH1F>("wfstartb","wfstartb",500,-3500.,6500.);
  wfstartb->GetXaxis()->SetTitle("Rejected Waveform Start (us)"); 

}

bool FlashTimeFilter::filter(art::Event& e)
{

  auto const& detectorClocks(*lar::providerFrom< detinfo::DetectorClocksService >());
    
  art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
  e.getByLabel(fInputModuleName, wvfHandle);

  if(!wvfHandle.isValid()){
    std::cout <<Form("Did not find any waveform") << std::endl;
  }

  int icount =0;
  for(auto const& wvf : (*wvfHandle)){
    auto TimeStart = wvf.TimeStamp();
    double relTime = TimeStart - detectorClocks.TriggerTime();
    // This is always 85086??
    //    int frame = detectorClocks.OpticalClock().Frame(TimeStart);
    //    std::cout << "WF begin: " << relTime << "  " << frame  << std::endl;
    if (relTime>=fWFTimeStart &&  relTime<=fWFTimeEnd) { icount++;
      std::cout << "WF begin: " << relTime << std::endl;
    } 

  }

  std::cout << icount << " " << fWFcount << std::endl;
  if (icount>=fWFcount) {
    std::cout << "kept" << std::endl;
    for(auto const& wvf : (*wvfHandle)){
      auto TimeStart = wvf.TimeStamp();
      double relTime = TimeStart - detectorClocks.TriggerTime();
      wfstarta->Fill(relTime);
    }
    return true;
  }

  // Implementation of required member function here.
  else {
    for(auto const& wvf : (*wvfHandle)){
      auto TimeStart = wvf.TimeStamp();
      double relTime = TimeStart - detectorClocks.TriggerTime();
      wfstartb->Fill(relTime);
    }
  }
  return false;
  
}// end of filter function




bool FlashTimeFilter::endSubRun(art::SubRun &subrun)
{

  return true;
}

void FlashTimeFilter::beginJob()
{
  // Implementation of optional member function here.
}

void FlashTimeFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FlashTimeFilter)
