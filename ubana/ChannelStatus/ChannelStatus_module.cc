////////////////////////////////////////////////////////////////////////
// Class:       ChannelStatus
// Plugin Type: analyzer (art v3_01_00)
// File:        ChannelStatus_module.cc
//
// Generated at Wed Feb 13 09:42:17 2019 by Wanwei Wu using cetskelgen
// from cetlib version v3_05_00.
// Read the channel status from output files of WireCell Toolkit and 
// save the information to a file.
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

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "ubcore/Geometry/UBOpReadoutMap.h"
#include "ubcore/Geometry/UBOpChannelTypes.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/MuCS/MuCSData.h"
#include "ubobj/MuCS/MuCSRecoData.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "TTree.h"
#include "TFile.h"
#include "TTimeStamp.h"

using namespace std;

class ChannelStatus;


class ChannelStatus : public art::EDAnalyzer {
public:
  explicit ChannelStatus(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ChannelStatus(ChannelStatus const&) = delete;
  ChannelStatus(ChannelStatus&&) = delete;
  ChannelStatus& operator=(ChannelStatus const&) = delete;
  ChannelStatus& operator=(ChannelStatus&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void reset();

private:
  // fhicle parameters
  std::string fTPC_badChannelListProducer;
  std::string fTPC_badChannelListLabel;
  bool fSaveTPC_badChannelList;
  //TFile *outputFile;

  // Declare member data here.
  TTree *fBadChannelTree;
  std::string outputFileName;
  
  // Event information
  int fEvent;
  int fRun;
  int fSubRun;
  double fEventTime;

  // TPC bad channel list
  vector<int> fBadChannel;
  vector<int> fBadBegin;
  vector<int> fBadEnd;
  

};


ChannelStatus::ChannelStatus(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
  fTPC_badChannelListProducer (p.get<std::string>("TPC_badChannelListProducer")),
  fTPC_badChannelListLabel (p.get<std::string>("TPC_badChannelListLabel")),
  fSaveTPC_badChannelList (p.get<bool>("SaveTPC_badChannelList"))
  //outputFileName (p.get< std::string >("outputFileName"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  //outputFile = new TFile(outputFileName.c_str(), "recreate");
  
  art::ServiceHandle<art::TFileService> tfs;
  fBadChannelTree = tfs->make<TTree>("badChannelTree", "");
  
  fBadChannelTree->Branch("eventNo", &fEvent);
  fBadChannelTree->Branch("runNo", &fRun);
  fBadChannelTree->Branch("subRunNo", &fSubRun);
  fBadChannelTree->Branch("eventTime", &fEventTime);


  if (fSaveTPC_badChannelList) {
    fBadChannelTree->Branch("badChannel", &fBadChannel);
    fBadChannelTree->Branch("badBegin", &fBadBegin);
    fBadChannelTree->Branch("badEnd", &fBadEnd);
  }
}


void ChannelStatus::reset(){
  if(fSaveTPC_badChannelList){
    fBadChannel.clear();
    fBadBegin.clear();
    fBadEnd.clear();
  }
}


void ChannelStatus::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  reset();
  fEvent = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();
  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fEventTime = tts.AsDouble();

  std::cout << "run number: " << fRun << std::endl;

  // process TPC bad channel list
  if (fSaveTPC_badChannelList) {
    art::Handle< std::vector<int> > badChannelHandle;
    if (!e.getByLabel(fTPC_badChannelListProducer,fTPC_badChannelListLabel, badChannelHandle)) {
      cout << "WARNING: no badchannel list with producer of " << fTPC_badChannelListProducer << ", or label of " << fTPC_badChannelListLabel << endl;
      return;
    }

    vector<int> const &badChannelVector(*badChannelHandle); 
    cout << "size of badChannelVector: " << badChannelVector.size() << endl;
    cout << "number of bad channels: " << (int)badChannelVector.size()/3 << endl;
    int numBadChannels = (int)badChannelVector.size()/3;
    for (int i=0; i<numBadChannels;i++) {
      const int offset = 3*i;
      fBadChannel.push_back(badChannelVector[offset+0]);
      fBadBegin.push_back(badChannelVector[offset+1]);
      fBadEnd.push_back(badChannelVector[offset+2]);
    }

  } // if fSaveTPC_badChannelList

  fBadChannelTree->Fill();

}

DEFINE_ART_MODULE(ChannelStatus)
