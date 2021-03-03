////////////////////////////////////////////////////////////////////////
// Class:       HitTreeAna
// Plugin Type: analyzer (art v2_11_03)
// File:        HitTreeAna_module.cc
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

//include the truth objects
#include "lardataobj/RecoBase/Hit.h"

namespace ana {
  class HitTreeAna;
}


class ana::HitTreeAna : public art::EDAnalyzer {
public:
  explicit HitTreeAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitTreeAna(HitTreeAna const &) = delete;
  HitTreeAna(HitTreeAna &&) = delete;
  HitTreeAna & operator = (HitTreeAna const &) = delete;
  HitTreeAna & operator = (HitTreeAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fInputTag;

  TTree* fHitTree;

  unsigned int fRun;
  unsigned int fSubrun;
  unsigned int fEvent;
  
  //hit info
  unsigned int            fHitIndex;
  unsigned int            fChannel;
  unsigned int            fStart_tick;
  unsigned int            fEnd_tick;
  float                   fPeak_time;
  float                   fSigma_peak_time;
  float                   fRms;
  float                   fPeak_amplitude;
  float                   fSigma_peak_amplitude;
  float                   fSummedADC;
  float                   fHit_integral;
  float                   fHit_sigma_integral;
  int                     fMultiplicity;
  int                     fLocal_index;
  float                   fGoodness_of_fit;
  int                     fDof;
  int                     fView;
  int                     fSignal_type;
  unsigned int            fWire;
  unsigned int            fPlane;
  unsigned int            fTPC;
  unsigned int            fCryostat;

  void SetupHitTree();

  void FillEventInfo(art::Event const&);
  void FillHitInfo(recob::Hit const&);

};

void ana::HitTreeAna::SetupHitTree()
{
  fHitTree->Branch("run",&fRun,"run/i");
  fHitTree->Branch("subrun",&fSubrun,"subrun/i");
  fHitTree->Branch("event",&fEvent,"event/i");

  fHitTree->Branch("hit_index",&fHitIndex,"hit_index/i");
  fHitTree->Branch("channel",&fChannel,"channel/i");
  fHitTree->Branch("view",&fView,"view/I");
  fHitTree->Branch("signal_type",&fSignal_type,"signal_type/I");
  fHitTree->Branch("wire",&fWire,"wire/i");
  fHitTree->Branch("plane",&fPlane,"plane/i");
  fHitTree->Branch("tpc",&fTPC,"tpc/i");
  fHitTree->Branch("cryostat",&fCryostat,"cryostat/i");

  fHitTree->Branch("start_tick",&fStart_tick,"start_tick/i");
  fHitTree->Branch("end_tick",&fEnd_tick,"end_tick/i");
  fHitTree->Branch("peak_time",&fPeak_time,"peak_time/F");
  fHitTree->Branch("sigma_peak_time",&fSigma_peak_time,"sigma_peak_time/F");
  fHitTree->Branch("rms",&fRms,"rms/F");
  fHitTree->Branch("peak_amplitude",&fPeak_amplitude,"peak_ampltiude/F");
  fHitTree->Branch("sigma_peak_amplitude",&fSigma_peak_amplitude,"sigma_peak_ampltiude/F");
  fHitTree->Branch("summedADC",&fSummedADC,"summedADC/F");
  fHitTree->Branch("integral",&fHit_integral,"integral/F");
  fHitTree->Branch("sigma_integral",&fHit_sigma_integral,"sigma_integral/F");
  fHitTree->Branch("multiplicity",&fMultiplicity,"multiplicity/I");
  fHitTree->Branch("local_index",&fLocal_index,"local_index/I");
  fHitTree->Branch("goodness_of_fit",&fGoodness_of_fit,"goodness_of_fit/F");
  fHitTree->Branch("dof",&fDof,"dof/I");
}

void ana::HitTreeAna::FillEventInfo(art::Event const& e)
{
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.event();
}

void ana::HitTreeAna::FillHitInfo(recob::Hit const& hit)
{

  fChannel = hit.Channel();
  fView = hit.View();
  fSignal_type = hit.SignalType();
  fWire = hit.WireID().Wire;
  fPlane = hit.WireID().Plane;
  fTPC = hit.WireID().TPC;
  fCryostat = hit.WireID().Cryostat;

  fStart_tick = hit.StartTick();
  fEnd_tick = hit.EndTick();
  fPeak_time = hit.PeakTime();
  fSigma_peak_time = hit.SigmaPeakTime();
  fRms = hit.RMS();
  fPeak_amplitude = hit.PeakAmplitude();
  fSigma_peak_amplitude = hit.SigmaPeakAmplitude();
  fSummedADC = hit.SummedADC();
  fHit_integral = hit.Integral();
  fHit_sigma_integral = hit.SigmaIntegral();
  fMultiplicity = hit.Multiplicity();
  fLocal_index = hit.LocalIndex();
  fGoodness_of_fit = hit.GoodnessOfFit();
  fDof = hit.DegreesOfFreedom();
}

ana::HitTreeAna::HitTreeAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
{
  fInputTag = p.get<art::InputTag>("InputTag");
}

void ana::HitTreeAna::analyze(art::Event const & e)
{
  //first fill the event-level information
  FillEventInfo(e);

  //grab the Hit collection
  auto const& hit_handle = e.getValidHandle< std::vector<recob::Hit> >(fInputTag);
  auto const& hit_vec(*hit_handle);

  //loop over hits and fill the tree
  for(size_t i_h=0; i_h<hit_vec.size(); ++i_h){
    fHitIndex = i_h;
    FillHitInfo(hit_vec[i_h]);
    fHitTree->Fill();
  }//end loop over hits

}//end analyze event

void ana::HitTreeAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fHitTree = tfs->make<TTree>("hit_tree","Hit Tree");
  SetupHitTree();
}

DEFINE_ART_MODULE(ana::HitTreeAna)
