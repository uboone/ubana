////////////////////////////////////////////////////////////////////////
// Class:       MCS
// Plugin Type: analyzer (art v3_01_02)
// File:        MCS_module.cc
//
// Generated at Wed May  1 16:37:29 2019 by Yifan Chen using cetskelgen
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

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"

#include "FiducialVolume.h"
#include "BackTrackerTruthMatch.h"

class MCS;


class MCS : public art::EDAnalyzer {
public:
  explicit MCS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCS(MCS const&) = delete;
  MCS(MCS&&) = delete;
  MCS& operator=(MCS const&) = delete;
  MCS& operator=(MCS&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  ::ubana::FiducialVolume _fiducial_volume;
 
  //::art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<geo::Geometry> geo;

  art::ServiceHandle<art::TFileService> tfs;
  TTree * my_event_;
  
  void Initialize_event();
 
  std::vector<double> mom_MCS_mu;//MCS momentum of muon track in the every event
  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
  std::vector<double> trk_length;//Range momentum of muon track in the every event

  std::string                         m_hitProducerLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_MCSmuProducerLabel;
  std::string                         Hits_TrackAssLabel;

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


MCS::MCS(fhicl::ParameterSet const& pset)
  : 
  EDAnalyzer{pset},
  m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_MCSmuProducerLabel(pset.get<std::string>("MCSmuProducerLabel")),
  _min_track_len{pset.get<double>("MinTrackLength", 0.1)},
  _trk_mom_calculator{_min_track_len}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _fiducial_volume.Configure(pset.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();
  Hits_TrackAssLabel = pset.get<std::string>("HitsPerTrackAssLabel");
}

void MCS::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  art::Handle<std::vector<recob::Hit> > rawHandle_Hit;
  evt.getByLabel(m_hitProducerLabel, rawHandle_Hit);

  art::Handle<std::vector<recob::Track> > rawHandle_TPCtrack;
  evt.getByLabel(m_trackProducerLabel, rawHandle_TPCtrack);
  std::vector< art::Ptr<recob::Track> > AllTrackCollection;
  art::fill_ptr_vector(AllTrackCollection, rawHandle_TPCtrack);

  art::Handle<std::vector<recob::MCSFitResult> > rawHandle_mcsfitresult_mu;
  evt.getByLabel(m_MCSmuProducerLabel, rawHandle_mcsfitresult_mu);
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
  art::fill_ptr_vector(mcsfitresult_mu_v, rawHandle_mcsfitresult_mu);

  art::FindManyP<recob::Hit> hits_per_track(rawHandle_TPCtrack,evt,Hits_TrackAssLabel);

  //std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(this_track_ptr.key());
  //BackTrackerTruthMatch backtrackertruthmatch;
  //backtrackertruthmatch.MatchToMCParticle(hit_handle,evt,trk_hits_ptrs);

  for(int trk_id = 0; trk_id < (int) AllTrackCollection.size(); trk_id++){
    // Back track matching
    std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
    BackTrackerTruthMatch backtrackertruthmatch;
    backtrackertruthmatch.MatchToMCParticle(rawHandle_Hit,evt,trk_hits_ptrs);
    if(!backtrackertruthmatch.ReturnMCParticle()){
      std::cout<<"MC particle does not exist!"<<std::endl;
    }
    else{
      std::cout<<"Track purity: "<< backtrackertruthmatch.ReturnPurity()<<std::endl;
      std::cout<<"Track ID: "<< backtrackertruthmatch.ReturnMCParticleID()<<std::endl;
      std::cout<<"Track pdg: "<< backtrackertruthmatch.ReturnMCParticle()->PdgCode()<<std::endl;
      //std::cout<<"Track pdg: "<< backtrackertruthmatch.ReturnMCParticle()->type()<<std::endl;
    }
    // Check if the track is contained or not
    bool contained = _fiducial_volume.InFV(AllTrackCollection[trk_id]->Vertex<TVector3>(), AllTrackCollection[trk_id]->End<TVector3>());

    double bestMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum();
    mom_MCS_mu.push_back(bestMCS);
    
    // Track Length
    double trk_len = AllTrackCollection[trk_id]->Length();
    trk_length.push_back(trk_len);

    // Range based momentum (Give muon pdg = 13), track length in [cm] ?
    double RangeMom = _trk_mom_calculator.GetTrackMomentum(trk_len, 13);
    mom_Range_mu.push_back(RangeMom);
    if(contained){
      std::cout<<"MCS: "<< bestMCS<<", range: "<<RangeMom<<std::endl;
    }
    else{
      std::cout<<"Not contained!----MCS: "<< bestMCS<<std::endl;
    }
  }
  my_event_->Fill();
}

void MCS::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;
  // Make a tree to store momentum information
  my_event_ = tfs->make<TTree>("tree","tree");

  my_event_->Branch("mom_MCS_mu", &mom_MCS_mu);
  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("trk_length", &trk_length);
}

void MCS::beginJob()
{
  // Implementation of optional member function here.
  Initialize_event();
}

void MCS::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MCS)
