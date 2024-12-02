////////////////////////////////////////////////////////////////////////
// Class:       NumuCCana
// Plugin Type: analyzer (art v2_05_01)
// File:        NumuCCana_module.cc
//
////////////////////////////////////////////////////////////////////////

/**
 * \brief CC-inclusive analyser
 *
 * \author Thomas Metter
 *
 * \email thomas.mettler@lhep.unibe.ch
 *
 * \notes This module is used as an analyser module for a numu CC inclusive analysis
 *
 */

// Base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h" 
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/EventWeight/Base/WeightManager.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "BackTrackerTruthMatch.h"
#include "PID.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "helpers/PandoraInterfaceHelper.h"
#include "helpers/EnergyHelper.h"
#include "helpers/TrackHelper.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubobj/CRT/CRTTzero.hh"
//#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "TTree.h"
#include <memory>

class NumuCCana;


class NumuCCana : public art::EDAnalyzer {
public:
  explicit NumuCCana(fhicl::ParameterSet const & p);

  NumuCCana(NumuCCana const &) = delete;
  NumuCCana(NumuCCana &&) = delete;
  NumuCCana & operator = (NumuCCana const &) = delete;
  NumuCCana & operator = (NumuCCana &&) = delete;

   void analyze(art::Event const & e) override;
  
  void beginJob() override;
  void endJob() override;
  
  void endSubRun(const art::SubRun &subrun) override;
  
  bool GetVertex(art::Event const &evt);
  bool GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                       const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                       const art::FindManyP<anab::ParticleID> &trackPIDAssn);

  void clearEvent(); // reset all vectors, UPDATE is new vectors are added!
  //void reconfigure(fhicl::ParameterSet const &p);
  
  void initialize_tree();
  void initialize_pot();
  
  void fillMCInfo(art::Event const &evt);
  void fillMCTruth(art::Event const &evt);
  //bool BacktrackDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp);

  
private:
  PandoraInterfaceHelper pandoraInterfaceHelper; // = new PandoraInterfaceHelper();
  EnergyHelper energyHelper;
  TrackHelper trackHelper;
  // LAr Pandora Helper fields
  lar_pandora::PFParticlesToClusters particlesToClusters;
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;
  lar_pandora::ClustersToHits clustersToHits;
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::TrackVector pftracks;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToTracks particlesToTracks;
  lar_pandora::ShowerVector pfshowers;
  lar_pandora::PFParticlesToShowers particlesToShowers;
  
  // Used for reco truth matching
  lar_pandora::PFParticlesToMCParticles matchedParticles;
  std::set<art::Ptr<simb::MCParticle>> matchedMCParticles;
  std::map<art::Ptr<recob::PFParticle>, float> matchedHitFractions;
  
  spacecharge::SpaceCharge const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  
  // MC truth variables:
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::GTruth> > GTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection;
  std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection_G4;
  double EventWeight; // Spine reweight using Steven's tool
  
  // MC neutrino daughter
  //std::vector<int> fTrueNu_DaughterPDG;
  //std::vector<float> fTrueNu_DaughterE;
  //std::vector<bool> fTrueNu_DaughterMatched;

  // fcl producer datalabels
  std::string m_pfp_producer;
  std::string m_hitfinder_producer;
  std::string m_geant_producer;
  std::string m_generatorLabel;
  std::string m_eventweightLabel;
  std::string m_eventweightLabel_G4;
  std::string m_hit_mcp_producer;
  std::string m_spacepointLabel;

  std::string data_label_DAQHeader_;
  std::string data_label_flash_beam_;
  std::string data_label_crthit_;
  std::string data_label_crtmchit_;
  std::string data_label_crtT0asso_;
  std::string data_label_crttricorr_;
  std::string Hits_TrackAssLabel;
  std::string PID_TrackAssLabel;
  
  // fcl parameters
  bool m_isData;
  bool m_hasMCNeutrino;
  bool store_Spacepoints;
  int     verbose_ = 0;
  int     fHardDelay_;
  int     fCRTT0off_;
  double  beam_start_ = 0;
  double  beam_end_ = 0;
  bool    is_data_ = false;
  int     fill_per_track_ = 0;
  double  crt_trig_corr_mean = 0;
  double  crt_trig_corr_med = 0;
  int     use_trig_corr = 0;
  int     no_DAQHeader = 0;
  int     store_weights_all = 0;
  
  art::ServiceHandle<art::TFileService> tfs;
  
  TTree* my_event_;
  
  int event_counter = 0;
  uint32_t fEvtNum;                //Number of current event                       
  uint32_t frunNum;                //Run Number taken from event  
  uint32_t fsubRunNum;             //Subrun Number taken from event 
  
  int has_neutrino_ = -1;
  double NuScore_ = -1;
  double FlashScore_ = -1;
  double FlashScoreTime_ = -1;
  int NuPDG_ = 0;
  int NumPfp_ = -1;
  double Nu_Vx_=-999,Nu_Vy_=-999,Nu_Vz_=-999;
  double Nu_Vx_sce_=-999,Nu_Vy_sce_=-999,Nu_Vz_sce_=-999;
  
  int NuTracks_ = 0;
  int NuShowers_ = 0;
  
  int Nu_NhitsU_ = 0, Nu_NhitsV_ = 0, Nu_NhitsY_ = 0;
  float Nu_CaloU_ = 0, Nu_CaloV_ = 0, Nu_CaloY_ = 0;
  
  // vector for all tracks
  
  std::vector<int> track_key_;
  std::vector<int> TrackPID_;

  std::vector<double> a_crthit_ts0;
  std::vector<double> a_crthit_ts1;
  std::vector<int> a_crthit_s;
  std::vector<double> a_crthit_x_;
  std::vector<double> a_crthit_y_;
  std::vector<double> a_crthit_z_;
  std::vector<int> a_adc_length;
  std::vector<int> a_crt_plane_;
  std::vector<double> a_crt_adc;
  std::vector<int> a_t0_counter;

  std::vector<double> crtt0_time_;
  std::vector<int> crtt0_trig_;
  std::vector<double> crtt0_DCA_;
  std::vector<int> crtt0_plane_;

  std::vector<double> VtxDistance_;
  std::vector<double> VtxDistance_sce_;
  std::vector<double> Vx_;
  std::vector<double> Vy_;
  std::vector<double> Vz_;
  std::vector<double> Vx_sce_;
  std::vector<double> Vy_sce_;
  std::vector<double> Vz_sce_;
  
  std::vector<int> TrackPfp_;
  std::vector<int> ShowerPfp_;
  std::vector<int> isShowerTrack_;
  
  std::vector<double> TrackScore_;
  std::vector<double> TrackScoreGlobal_;
  std::vector<double> ShowerScore_;
  std::vector<double> TrackLength_;
  std::vector<double> TrackPID_chiproton_;
  std::vector<double> TrackPID_chimuon_;
  std::vector<double> TrackPID_chipion_;
  std::vector<double> TrackPID_chikaon_;

  std::vector<float> TrackMomRange_p_;
  std::vector<float> TrackMomRange_mu_;
  std::vector<float> TrackMomMCS_mom_;
  std::vector<float> TrackMomMCS_err_;
  std::vector<float> TrackMomMCS_ll_;

  std::vector<double> TrackStart_x_;
  std::vector<double> TrackStart_y_;
  std::vector<double> TrackStart_z_;
  std::vector<double> TrackStart_x_sce_;
  std::vector<double> TrackStart_y_sce_;
  std::vector<double> TrackStart_z_sce_;
  std::vector<double> TrackEnd_x_;
  std::vector<double> TrackEnd_y_;
  std::vector<double> TrackEnd_z_;
  std::vector<double> TrackEnd_x_sce_;
  std::vector<double> TrackEnd_y_sce_;
  std::vector<double> TrackEnd_z_sce_;
  std::vector<double> TrackDir_x_;
  std::vector<double> TrackDir_y_;
  std::vector<double> TrackDir_z_;
  
  std::vector<double> MCTrackPurity_;
  std::vector<double> MCTrackPDG_;
  std::vector<double> MCTrackEnergy_;
  std::vector<double> MCTrackMomentum_;
  std::vector<double> MCTrackTheta_;
  std::vector<double> MCTrackPhi_;
  std::vector<double> MCTrackLength_;
  
  std::vector<double> MCTrackStart_x_;
  std::vector<double> MCTrackStart_y_;
  std::vector<double> MCTrackStart_z_;
  std::vector<double> MCTrackEnd_x_;
  std::vector<double> MCTrackEnd_y_;
  std::vector<double> MCTrackEnd_z_;
  
  std::vector<double> AllTrack_point_x_;
  std::vector<double> AllTrack_point_y_;
  std::vector<double> AllTrack_point_z_;
  
  std::vector<double> SpacePoint_x_;
  std::vector<double> SpacePoint_y_;
  std::vector<double> SpacePoint_z_;
  
  std::vector<unsigned int> Wire_id_;
  std::vector<unsigned int> Wire_plane_;

  std::vector<double> TrackTheta_;
  std::vector<double> TrackPhi_;
  
  std::vector<int> TrackNHitsU_;
  std::vector<int> TrackNHitsV_;
  std::vector<int> TrackNHitsY_;
  std::vector<float> TrackCaloU_;
  std::vector<float> TrackCaloV_;
  std::vector<float> TrackCaloY_;
  
  std::vector<double> ShowerLength_;
  std::vector<double> ShowerOpenAngle_;
  std::vector<double> ShowerEnergy_;
  std::vector<double> ShowerMIPEnergy_;
  std::vector<double> ShowerDir_x_;
  std::vector<double> ShowerDir_y_;
  std::vector<double> ShowerDir_z_;
  std::vector<float> Shower_dEdxU_;
  std::vector<float> Shower_dEdxV_;
  std::vector<float> Shower_dEdxY_;
  std::vector<int> Shower_dEdxHitsU_;
  std::vector<int> Shower_dEdxHitsV_;
  std::vector<int> Shower_dEdxHitsY_;
  std::vector<float> Shower_dEdxPitchU_;
  std::vector<float> Shower_dEdxPitchV_;
  std::vector<float> Shower_dEdxPitchY_;
  
  std::vector<double> AllShower_point_x_;
  std::vector<double> AllShower_point_y_;
  std::vector<double> AllShower_point_z_;
  
  std::vector<int> ShowerNHitsU_;
  std::vector<int> ShowerNHitsV_;
  std::vector<int> ShowerNHitsY_;
  std::vector<float> ShowerCaloU_;
  std::vector<float> ShowerCaloV_;
  std::vector<float> ShowerCaloY_;
  // end track variables
  int muon_candidate_key = 0;
  int muon_candidate_pfp = 0;
  
  double TriTim_sec_ = 0;          //event trigger time sec
  double TriTim_nsec_ = 0;          //event trigger time ns
  
  // CRT in beam variables
  
  int nr_crthit_ = -1; // # crt hits assigned to a tpc track
  std::vector<double> crthit_ts0_;
  std::vector<double> crthit_ts1_;
  std::vector<int> crthit_s_;
  std::vector<double> crthit_x_;
  std::vector<double> crthit_y_;
  std::vector<double> crthit_z_;
  std::vector<int> adc_length_;
  std::vector<int> crt_plane_;
  std::vector<double> crt_adc_;
  std::vector<int> crtbeam_hit_nr_;
  
  // flash information
  double TimFla_ = -99;
  double flash_PE_ = -99;
  double flash_y_ = -999;
  double flash_z_ = -999;
  
  // MC neutrino info
  uint NuMCnu; // number of MC neutrinos in event, only one gets saved!
  int MCNu_Interaction;
  int MCNu_CCNC;
  int MCNu_PDG;
  float MCNu_Energy;
  float MCNu_leptonPx, MCNu_leptonPy, MCNu_leptonPz;
  float MCNu_LeptonEnergy;
  float MCNu_Px, MCNu_Py, MCNu_Pz;
  float MCNu_leptonTheta;
  float MCNu_time; // time of the true neutrino interaction
  float MCNu_Vx, MCNu_Vy, MCNu_Vz;
  float MCNu_VxSce, MCNu_VySce, MCNu_VzSce;
  //float MCNu_Vx_sce, MCNu_Vy_sce, MCNu_Vz_sce;
  float MCNu_vertexDistance;
  //float MCNu_vertexDistance_sce;
  
  // Matched MCParticle info
  bool MCNU_matched;
  bool MCCosmic_matched;
  //MC truth information
  std::vector<int> MCle_pfp_;
  std::vector<int> MCle_key_;
  std::vector<int> MCle_PDG_;
  std::vector<float> MCle_Energy_;
  std::vector<float> MCle_Px_;
  std::vector<float> MCle_Py_;
  std::vector<float> MCle_Pz_;
  std::vector<float> MCle_Vx_;
  std::vector<float> MCle_Vy_;
  std::vector<float> MCle_Vz_;
  std::vector<float> MCle_Endx_;
  std::vector<float> MCle_Endy_;
  std::vector<float> MCle_Endz_;
  std::vector<float> MCle_length_;
  std::vector<float> MCle_VxSce_;
  std::vector<float> MCle_VySce_;
  std::vector<float> MCle_VzSce_;
  std::vector<float> MCle_Theta_;
  std::vector<float> MCle_Phi_;
  // Genie information
  double Genie_Q2;
  double Genie_q2;
  double Genie_W;
  double Genie_T;
  double Genie_X;
  double Genie_Y;
  int Genie_nNeutron_preFSI;// before FSI 
  int Genie_nProton_preFSI;// before FSI 
  int Genie_nPi0_preFSI;// before FSI 
  int Genie_nPiPlus_preFSI;// before FSI 
  int Genie_nPiMinus_preFSI;// before FSI 
  
  std::vector<double> para[100];
  //char para_name[100][200];

  TTree* _sr_tree;
  int _sr_run = -9999;
  int _sr_subrun = -9999;
  double _sr_begintime = -9999;
  double _sr_endtime = -9999;
  double _sr_pot = -9999;
  //int event_counter = 0;
  
  // variables for cut
  
  //cout counters
  
  int counter_moun = 0;
  int counter_num_nu = 0;
  int counter_trackPIDAssn = 0;
  int counter_neutrino_metadata_vec = 0;
  int counter_getPID = 0;
};/*
void NumuCCana::reconfigure(fhicl::ParameterSet const &p)
{
    m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraConsolidated");
    m_hitfinder_producer = p.get<std::string>("hitfinder_producer", "gaushit");
    m_geant_producer = p.get<std::string>("geant_producer", "largeant");
    m_generatorLabel = p.get<std::string>("GeneratorLabel","generator"),
    m_hit_mcp_producer = p.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
    m_spacepointLabel = p.get<std::string>("SpacePointModule", "pandora");
    energyHelper.reconfigure(p);
}*/

NumuCCana::NumuCCana(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{
    //producer datalabels
    m_isData =          pset.get<bool>("is_data", false);
    m_hasMCNeutrino =   pset.get<bool>("has_MC_neutrino", true);
    store_Spacepoints = pset.get<bool>("store_Spacepoints", false);

    //this->reconfigure(pset);
      
    m_pfp_producer =        pset.get<std::string>("pfp_producer", "pandoraConsolidated");
    m_hitfinder_producer =  pset.get<std::string>("hitfinder_producer", "gaushit");
    m_geant_producer =      pset.get<std::string>("geant_producer", "largeant");
    m_generatorLabel =      pset.get<std::string>("GeneratorLabel","generator"),
    m_eventweightLabel =    pset.get<std::string>("eventweightLabel","eventweight"),
    m_eventweightLabel_G4 =    pset.get<std::string>("eventweightLabel_G4","eventweightreint"),
    m_hit_mcp_producer =    pset.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
    m_spacepointLabel =     pset.get<std::string>("SpacePointModule", "pandora");
    energyHelper.reconfigure(pset);

    data_label_DAQHeader_ =   pset.get<std::string>("data_label_DAQHeader");
    data_label_flash_beam_ =  pset.get<std::string>("data_label_flash_beam");
    data_label_crthit_ =      pset.get<std::string>("data_label_crthit");
    data_label_crtmchit_ =    pset.get<std::string>("data_label_crtmchit");
    data_label_crtT0asso_ =   pset.get<std::string>("data_label_crtT0asso");
    data_label_crttricorr_ =  pset.get<std::string>("data_label_crttricorr");
      
    fHardDelay_ =             pset.get<int>("fHardDelay",40000);
    fCRTT0off_ =              pset.get<int>("fCRTT0off",69000);
    beam_start_ =             pset.get<double>("beam_start",3.2);
    beam_end_ =               pset.get<double>("beam_end",5);
    is_data_ =                pset.get<bool>("is_data",false);
    fill_per_track_ =         pset.get<int>("fill_per_track",0);
    use_trig_corr =           pset.get<int>("use_trig_corr",0);
    no_DAQHeader =            pset.get<int>("no_DAQHeader",0);
    store_weights_all =           pset.get<int>("store_weights_all",0);

    Hits_TrackAssLabel = pset.get<std::string>("HitsPerTrackAssLabel","pandora");
    PID_TrackAssLabel = pset.get<std::string>("PIDTrackAssLabel","pandoracalipidSCE");

    verbose_ = pset.get<int>("verbose");
      
    std::cout.width(45); std::cout << std::left << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] Checking set-up" << std::left << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] pfp_producer: " << std::left << m_pfp_producer << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] hitfinder_producer: " << std::left << m_hitfinder_producer << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] geant_producer: " << std::left << m_geant_producer << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] generatorLabel: " << std::left << m_generatorLabel << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] eventweightLabel: " << std::left << m_eventweightLabel << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] eventweightLabel_G4: " << std::left << m_eventweightLabel_G4 << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] hit_mcp_producer: " << std::left << m_hit_mcp_producer << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] spacepointLabel: " << std::left << m_spacepointLabel << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_DAQHeader: " << std::left << data_label_DAQHeader_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_flash_beam: " << std::left << data_label_flash_beam_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_crthit: " << std::left << data_label_crthit_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_crtmchit: " << std::left << data_label_crtmchit_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_crtT0asso: " << std::left << data_label_crtT0asso_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] data_label_crttricorr: " << std::left << data_label_crttricorr_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] HitsPerTrackAssLabel: " << std::left << Hits_TrackAssLabel << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor]" << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] fHardDelay: " << std::left << fHardDelay_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] fCRTT0off: " << std::left << fCRTT0off_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] beam_start: " << std::left << beam_start_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] beam_end: " << std::left << beam_end_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] is_data: " << std::left << is_data_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] fill_per_track: " << std::left << fill_per_track_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] use_trig_corr: " << std::left << use_trig_corr << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] no_DAQHeader: " << std::left << no_DAQHeader << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] store_weights_all: " << std::left << store_weights_all << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] verbose: " << std::left << verbose_ << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] is_data: " << std::left << m_isData << std::endl;
    std::cout.width(45); std::cout << "[NuCC constructor] has_MC_neutrino: " << std::left << m_hasMCNeutrino << std::endl;
}

void NumuCCana::analyze(art::Event const & evt)
{
  clearEvent();
  std::cout << "[NumuCCana::analyse] Prozessing event nr: " << event_counter << std::endl;
  if(verbose_!=0) std::cout << "[NumuCCana::analyse] Run " << evt.run() << ", subrun " << evt.subRun() << std::endl;
  frunNum    = evt.run();
  fsubRunNum = evt.subRun();
  fEvtNum = evt.event();
  event_counter++;
  // get event time for CRT time calculation
  if(no_DAQHeader !=1){
    art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
    evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
    raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
    art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
    TriTim_sec_ = evtTimeGPS.timeHigh();
    TriTim_nsec_ = evtTimeGPS.timeLow();
  }
  // use only if CRT time has to be corrected... (MC has not to be corrected)
  if(use_trig_corr ==1){
    art::Handle< std::vector<anab::T0> > rawHandle_CRTTriCorr;
    evt.getByLabel(data_label_crttricorr_, rawHandle_CRTTriCorr); //mergerextra
    std::vector<anab::T0> const& CRTTriCorrCollection(*rawHandle_CRTTriCorr);
    for(std::vector<int>::size_type i = 0; i != CRTTriCorrCollection.size(); i++) {
      crt_trig_corr_med = CRTTriCorrCollection.at(i).fTime;
      crt_trig_corr_mean = CRTTriCorrCollection.at(i).fTriggerConfidence;
    }
  }
  //get spacepoints
  if(store_Spacepoints){
    // store space points for event display
    if(verbose_!=0) std::cout << "[NumuCCana::analyse] Storing space point information" << std::endl;
    lar_pandora::SpacePointVector         spacePointVector;
    lar_pandora::SpacePointsToHits        spacePointsToHits;
    larpandora.CollectSpacePoints(evt, m_spacepointLabel, spacePointVector, spacePointsToHits);
    for (auto const spacepoint : spacePointVector){
      auto space_cor = spacepoint->XYZ();
      SpacePoint_x_.push_back(space_cor[0]);
      SpacePoint_y_.push_back(space_cor[1]);
      SpacePoint_z_.push_back(space_cor[2]);
    }
    // store wire hits for eventual cut, only first 10 wires of collection plane
    if(verbose_!=0) std::cout << "[NumuCCana::analyse] Storing wire hit information" << std::endl;
    lar_pandora::HitVector         hitvector;
    larpandora.CollectHits(evt, "gaushit", hitvector);
    for (auto const hit : hitvector){
      auto this_wire = hit->WireID().Wire;
      auto this_plane = hit->WireID().Plane;
      if(this_wire<11 && this_plane==2){
        Wire_id_.push_back(this_wire );
        Wire_plane_.push_back(this_plane );
      }
    }
  }
  // get CRT hits in the beam window
  if(verbose_!=0) std::cout << "[NumuCCana::analyse] Getting CRT Hit information" << std::endl;
  art::Handle< std::vector<crt::CRTHit> > rawHandle_CRTHit;
  evt.getByLabel(data_label_crthit_, rawHandle_CRTHit); //mergerextra
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_CRTHit);
  for(std::vector<int>::size_type i = 0; i != CRTHitCollection.size(); i++) {
    double crthittime = 0;
    //std::cout << "CRT hit: ts0= " << CRTHitCollection.at(i).ts0_ns << " adc= " << CRTHitCollection.at(i).peshit << " with length= " << CRTHitCollection.at(i).pesmap.begin()->second.size() << std::endl;
    //std::cout << "CRT offset: " << fCRTT0off_ << " Tritimensec: " << TriTim_nsec_ << std::endl;
    if(use_trig_corr == 1) crthittime = (CRTHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0;
    else if(CRTHitCollection.at(i).pesmap.begin()->second.size() != 32) crthittime = (CRTHitCollection.at(i).ts0_ns - TriTim_nsec_)/1000.0;
    else crthittime = (CRTHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0;
    //std::cout << "CRT hittime: " << crthittime << " condition test : " << abs( crthittime - (beam_end_ + beam_start_)/2 ) << " < " << (beam_end_-beam_start_)/2 << std::endl;
    if( abs( crthittime - (beam_end_ + beam_start_)/2 ) < (beam_end_-beam_start_)/2 ){
        if(verbose_!=0){
          std::cout << "################################################" << std::endl;
          std::cout << "CRT Hit in beam window:" << std::endl;
          std::cout << "## Time GPS [us]:\t" << crthittime << std::endl;
          std::cout << "################################################" << std::endl;
        }
        // fill tree variables
        //if(use_trig_corr == 1) crthit_ts0_.push_back((double)(CRTHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0);
        //else crthit_ts0_.push_back((double)(CRTHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0);
        crthit_ts0_.push_back((double)(crthittime));
        crthit_ts1_.push_back(((double)CRTHitCollection.at(i).ts1_ns + fHardDelay_)/1000.0);
        crthit_s_.push_back(CRTHitCollection.at(i).ts0_s);
        crt_adc_.push_back(CRTHitCollection.at(i).peshit);
        crthit_x_.push_back(CRTHitCollection.at(i).x_pos);
        crthit_y_.push_back(CRTHitCollection.at(i).y_pos);
        crthit_z_.push_back(CRTHitCollection.at(i).z_pos);
        crt_plane_.push_back(CRTHitCollection.at(i).plane);
        adc_length_.push_back(CRTHitCollection.at(i).pesmap.begin()->second.size());
        nr_crthit_++;
        crtbeam_hit_nr_.push_back(nr_crthit_);
        //std::cout << "CRT got in!!!: x=" << CRTHitCollection.at(i).z_pos << " peshit: " << CRTHitCollection.at(i).peshit << std::endl;
    }
  }
  // Add the following part for V22 or bigger, since the mC CRT hits are not in the crthitcorr collection
  /*if(verbose_!=0) std::cout << "[NumuCCana::analyse] Getting CRT MC Hit information" << std::endl;
  art::Handle< std::vector<crt::CRTHit> > rawHandle_CRTMCHit;
  evt.getByLabel(data_label_crtmchit_, rawHandle_CRTMCHit); //mergerextra
  std::vector<crt::CRTHit> const& CRTMCHitCollection(*rawHandle_CRTMCHit);
  for(std::vector<int>::size_type i = 0; i != CRTMCHitCollection.size(); i++) {
    std::cout << "MC CRT hit: ts0= " << CRTMCHitCollection.at(i).ts0_ns << " adc= " << CRTMCHitCollection.at(i).peshit << " with length= " << CRTMCHitCollection.at(i).pesmap.begin()->second.size() << std::endl;
    std::cout << "MC CRT offset: " << fCRTT0off_ << " Tritimensec: " << TriTim_nsec_ << std::endl;
    
    
    double crthittime = 0;
    if(use_trig_corr == 1) crthittime = (CRTMCHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0;
    else crthittime = (CRTMCHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0;
    std::cout << "MC CRT hittime: " << crthittime << " condition test : " << abs( crthittime - (beam_end_ + beam_start_)/2 ) << " < " << (beam_end_-beam_start_)/2 << std::endl;
    if( abs( crthittime - (beam_end_ + beam_start_)/2 ) < (beam_end_-beam_start_)/2 ){
        if(verbose_!=0){
          std::cout << "################################################" << std::endl;
          std::cout << "CRT MC Hit in beam window:" << std::endl;
          std::cout << "## Time GPS [us]:\t" << crthittime << std::endl;
          std::cout << "################################################" << std::endl;
        }
        // fill tree variables
        if(use_trig_corr == 1) crthit_ts0_.push_back((double)(CRTMCHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0);
        else crthit_ts0_.push_back((double)(CRTMCHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0);
        crthit_ts1_.push_back(((double)CRTMCHitCollection.at(i).ts1_ns + fHardDelay_)/1000.0);
        crt_adc_.push_back(CRTMCHitCollection.at(i).peshit);
        crthit_x_.push_back(CRTMCHitCollection.at(i).x_pos);
        crthit_y_.push_back(CRTMCHitCollection.at(i).y_pos);
        crthit_z_.push_back(CRTMCHitCollection.at(i).z_pos);
        crt_plane_.push_back(CRTMCHitCollection.at(i).plane);
        adc_length_.push_back(CRTMCHitCollection.at(i).pesmap.begin()->second.size());
        nr_crthit_++;
        crtbeam_hit_nr_.push_back(nr_crthit_);
        std::cout << "MC CRT got in!!!: x=" << CRTMCHitCollection.at(i).z_pos << " peshit: " << CRTMCHitCollection.at(i).peshit << std::endl;
    }
  }*/
  
  //get the beam flash
  if(verbose_!=0) std::cout << "[NumuCCana::analyse] Getting OpFlash information" << std::endl;
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlashBeam;
  evt.getByLabel(data_label_flash_beam_, rawHandle_OpFlashBeam);
  std::vector<recob::OpFlash> const& OpFlashCollectionBeam(*rawHandle_OpFlashBeam);
  if(verbose_!=0) std::cout << "There are: " << OpFlashCollectionBeam.size() << " in the beamflash collection" << std::endl;
  for(std::vector<int>::size_type i = 0; i != OpFlashCollectionBeam.size(); i++) {
    if( abs( OpFlashCollectionBeam.at(i).Time() - (beam_end_ + beam_start_)/2 ) < (beam_end_-beam_start_)/2 ){
      if(verbose_!=0){
      std::cout << "################################################" << std::endl;
      std::cout << "Flash in beam window:" << std::endl;
      std::cout << "## Time [us]:\t\t" << OpFlashCollectionBeam.at(i).Time() << std::endl;
      std::cout << "## Total PE []:\t" << OpFlashCollectionBeam.at(i).TotalPE() << std::endl;
      std::cout << "## Y [cm]:\t\t" << OpFlashCollectionBeam.at(i).YCenter() << std::endl;
      std::cout << "## Z [cm]:\t\t" << OpFlashCollectionBeam.at(i).ZCenter() << std::endl;
      std::cout << "################################################" << std::endl;
      }
      //fill tree variables
      TimFla_ = OpFlashCollectionBeam.at(i).Time();
      flash_PE_ = OpFlashCollectionBeam.at(i).TotalPE();
      flash_y_ = OpFlashCollectionBeam.at(i).YCenter();
      flash_z_ = OpFlashCollectionBeam.at(i).ZCenter();
    }
  }
  if(verbose_!=0) std::cout << "[NumuCCana::analyse] Get metadata and build particle map" << std::endl;
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  if (pfparticles.size() == 0){ //this should never happen
    std::cout << "[NumuCCana::analyse] Event failed: No reconstructed PFParticles in event is " << pfparticles.size() << std::endl;
    if (!is_data_) // for all events, for efficiency calculation
      {
        fillMCInfo(evt);
        fillMCTruth(evt);
      }
  }
  else{
    //Get the neutrino candidate info
    if(verbose_!=0) std::cout << "[NumuCCana::analyse] Get neutrino candidate" << std::endl;
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1){ // some events (e.g out of TPC) have no neutrino reconstructed
      if(verbose_!=0) std::cout << "[NumuCCana::analyse] Event failed: Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
      has_neutrino_= 0;
      counter_num_nu++;
      if (!is_data_) // new for >= v08_00_00_16
      {
        fillMCInfo(evt);
        fillMCTruth(evt);
      }
    }
    else{ //check if neutrino candidate is muon neutrino
      has_neutrino_ = 1;
      if (!is_data_)
      {
        fillMCInfo(evt);
        fillMCTruth(evt);
      }
      if(verbose_!=0) std::cout << "[NumuCCana::analyse] Get vertex" << std::endl;
      if(!GetVertex(evt)){//Try to find its daughter particles for further studies
        std::cout << "[NumuCCana::analyse] Event failed: GetVertex not passed" << std::endl;
      }
    }
  }
  if(verbose_!=0) std::cout << "[NumuCCana::analyse] Event finished." << std::endl;
  if(fill_per_track_!=1) my_event_->Fill();
}

void NumuCCana::fillMCTruth(art::Event const &evt){
  if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] Start filling MC truth variables" << std::endl;
  if(!is_data_){
    const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;
    bool MC_beamNeutrino = false;
    if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] MCTruthCollection size: " << MCTruthCollection.size() << std::endl;
    for(unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
      if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] i_mc: " << i_mc << std::endl;
      //if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin) MCNU_matched = true;
      //else MCCosmic_matched = true;
      MCNu_Interaction = MCTruthCollection[i_mc]->GetNeutrino().Mode();
      MCNu_PDG = MCTruthCollection[i_mc]->GetNeutrino().Nu().PdgCode();
      if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] MCNu_PDG: " << MCNu_PDG << std::endl;
      MCNu_CCNC = MCTruthCollection[i_mc]->GetNeutrino().CCNC();

      MCNu_Energy = MCTruthCollection[i_mc]->GetNeutrino().Nu().E();
      MCNu_Vx = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vx();
      MCNu_Vy = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vy();
      MCNu_Vz = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vz();
      MCNu_leptonTheta = MCTruthCollection[i_mc]->GetNeutrino().Theta();
      MCNu_time = MCTruthCollection[i_mc]->GetNeutrino().Nu().T();
      if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin) MC_beamNeutrino = true;
      if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] MC_beamNeutrino: " << MC_beamNeutrino << std::endl;
    }
    for(unsigned int i_mcp = 0; i_mcp < MCParticleCollection.size(); i_mcp++){
      if(MCParticleCollection[i_mcp]->Process() == "primary"){
        // PDG and momemtum of neutrino daughters
        MCle_PDG_.push_back(MCParticleCollection[i_mcp]->PdgCode());
        MCle_Energy_.push_back(MCParticleCollection[i_mcp]->P());
        MCle_Px_.push_back(MCParticleCollection[i_mcp]->Px());
        MCle_Py_.push_back(MCParticleCollection[i_mcp]->Py());
        MCle_Pz_.push_back(MCParticleCollection[i_mcp]->Pz());
        MCle_Vx_.push_back(MCParticleCollection[i_mcp]->Vx());
        MCle_Vy_.push_back(MCParticleCollection[i_mcp]->Vy());
        MCle_Vz_.push_back(MCParticleCollection[i_mcp]->Vz());
        MCle_Endx_.push_back(MCParticleCollection[i_mcp]->EndX());
        MCle_Endy_.push_back(MCParticleCollection[i_mcp]->EndY());
        MCle_Endz_.push_back(MCParticleCollection[i_mcp]->EndZ());
        if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] MCle_PDG: " << MCle_PDG_.back() << std::endl;
        //MCle_length_.push_back(MCParticleCollection[i_mcp]->Length());
        
        TVector3 TrueTrackPos(MCParticleCollection[i_mcp]->Px(), MCParticleCollection[i_mcp]->Py(), MCParticleCollection[i_mcp]->Pz());// The initial momentum represent the angle of true track
        MCle_Theta_.push_back(TrueTrackPos.Theta());
        MCle_Phi_.push_back(TrueTrackPos.Phi());
        TVector3 true_start(MCParticleCollection[i_mcp]->Position().X(), MCParticleCollection[i_mcp]->Position().Y(), MCParticleCollection[i_mcp]->Position().Z());
        TVector3 true_end(MCParticleCollection[i_mcp]->EndPosition().X(), MCParticleCollection[i_mcp]->EndPosition().Y(), MCParticleCollection[i_mcp]->EndPosition().Z());
        MCle_length_.push_back((true_start - true_end).Mag());
        
        //TVector3 ParticleStart, ParticleEnd;
        //ParticleStart.SetXYZ(MCParticleCollection[i_mcp]->Vx(), MCParticleCollection[i_mcp]->Vy(),MCParticleCollection[i_mcp]->Vz());
        //ParticleEnd.SetXYZ(MCParticleCollection[i_mcp]->EndX(), MCParticleCollection[i_mcp]->EndY(),MCParticleCollection[i_mcp]->EndZ());
        //MCle_Theta_.push_back(cos((ParticleEnd - ParticleStart).Theta()));
        //if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] MClepton Theta: " << cos((ParticleEnd - ParticleStart).Theta()) << std::endl;
        //MCle_Phi_.push_back((ParticleEnd - ParticleStart).Phi());

      }// end primary
    } // end loop over mcparticlecollection
    // Get Genie info on how many particles produced
    for(unsigned int i_gn = 0; i_gn < GTruthCollection.size(); i_gn++){
      if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] Genie_Q2: " << GTruthCollection[i_gn]->fgQ2 << std::endl;
      Genie_Q2 = GTruthCollection[i_gn]->fgQ2;
      Genie_q2 = GTruthCollection[i_gn]->fgq2;
      Genie_W = GTruthCollection[i_gn]->fgW;
      Genie_T = GTruthCollection[i_gn]->fgT;
      Genie_X = GTruthCollection[i_gn]->fgX;
      Genie_Y = GTruthCollection[i_gn]->fgY;
      Genie_nNeutron_preFSI = GTruthCollection[i_gn]->fNumNeutron;
      Genie_nProton_preFSI = GTruthCollection[i_gn]->fNumProton;
      Genie_nPi0_preFSI = GTruthCollection[i_gn]->fNumPi0;
      Genie_nPiPlus_preFSI = GTruthCollection[i_gn]->fNumPiPlus;
      Genie_nPiMinus_preFSI = GTruthCollection[i_gn]->fNumPiMinus;
    }
  }
  if(verbose_!=0) std::cout << "[NumuCCana::fillMCTruth] End filling MC truth variables" << std::endl;
}
void NumuCCana::fillMCInfo(art::Event const &evt){
  if(verbose_!=0) std::cout << "[NumuCCana::fillMCInfo] Start filling MC truth vectors" << std::endl;
  if(!is_data_){
    std::cout << "Event counter: " << event_counter << std::endl;
    if(event_counter==1){
      int counter_para = 0;
      art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight;
      if(evt.getByLabel(m_eventweightLabel, Handle_Weight)){
        art::fill_ptr_vector(WeightCollection, Handle_Weight);
        std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(0)->fWeight;
        std::cout << "Vector length: " << evtwgt_map.size() << std::endl;
        for (auto const& evt_weight : evtwgt_map){
          std::cout << "branching parameter: " << counter_para << " name: " << evt_weight.first << " length: " << evt_weight.second.size() << std::endl;
          my_event_->Branch((evt_weight.first).c_str(),             &para[counter_para]);
          counter_para++;
        }
      }
      art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight_G4;
      if(evt.getByLabel(m_eventweightLabel_G4, Handle_Weight_G4)){
        art::fill_ptr_vector(WeightCollection_G4, Handle_Weight_G4);
        std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection_G4.at(0)->fWeight;
        for (auto const& evt_weight : evtwgt_map){
          std::cout << "branching parameter: " << counter_para << " name: " << evt_weight.first << " length: " << evt_weight.second.size() << std::endl;
          //my_event_->Branch((evt_weight.first).c_str(),             &para[counter_para]);
          counter_para++;
        }
      }
    }
    // MC Truth
    art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

    // Genie Truth
    art::Handle< std::vector<simb::GTruth> > Handle_GTruth;
    evt.getByLabel(m_generatorLabel, Handle_GTruth);
    art::fill_ptr_vector(GTruthCollection, Handle_GTruth);
 
    // MC Particle
    art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
    evt.getByLabel(m_geant_producer, Handle_MCParticle);
    art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

    // Reweight
    art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight;
    int counter_para = 0;
    if(evt.getByLabel(m_eventweightLabel, Handle_Weight)){
      art::fill_ptr_vector(WeightCollection, Handle_Weight);
      std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(0)->fWeight;
      
      for (auto const& evt_weight : evtwgt_map){
        if( store_weights_all == 1){
          para[counter_para] = evt_weight.second;
        }
        else if( evt_weight.first == "TunedCentralValue_Genie"){
          para[counter_para] = evt_weight.second;
        }
        else if( evt_weight.first == "splines_general_Spline"){
          para[counter_para] = evt_weight.second;
        }
        else if( store_weights_all == 2 && evt_weight.first == "All_Genie"){
          para[counter_para] = evt_weight.second;
        }
        counter_para++;
        //std::cout << "Type: " << evt_weight.first << " weight: " << std::endl; //<< evt_weight.second << std::endl;
        //std::cout << "Type: " << evt_weight.first << " weight: " << evt_weight.second.size() << std::endl;
        //for(auto const weights : evt_weight.second){
        //  std::cout << weights << " - " ;
        //}
        //std::cout << std::endl;
      }
      const std::vector<double> &weights = evtwgt_map.at("splines_general_Spline");
      EventWeight = weights.front();
    }
    else{
      EventWeight = 1;
    }
    art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight_G4;
    if(evt.getByLabel(m_eventweightLabel_G4, Handle_Weight_G4)){
      art::fill_ptr_vector(WeightCollection_G4, Handle_Weight_G4);
      std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection_G4.at(0)->fWeight;
      
      for (auto const& evt_weight : evtwgt_map){
        if( store_weights_all == 1){
          para[counter_para] = evt_weight.second;
        }
        counter_para++;
      }
    }
  }
  if(verbose_!=0) std::cout << "[NumuCCana::fillMCInfo] End filling MC truth vectors" << std::endl;
}

bool NumuCCana::GetVertex(art::Event const &evt)
{
  if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Get vertex" << std::endl;
  //NumPfp_ = pfparticles.size();
  // Get vertex information
  lar_pandora::VertexVector vertexVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);
  lar_pandora::ClusterVector clusterVector_dummy;
  larpandora.CollectClusters(evt, m_pfp_producer, clusterVector_dummy, clustersToHits);
  lar_pandora::PFParticleVector particleVector_dummy;
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToClusters);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToSpacePoints);
  
  // get information of downstream tracks
  larpandora.CollectTracks(evt, m_pfp_producer, pftracks, particlesToTracks);
  larpandora.CollectShowers(evt, m_pfp_producer, pfshowers, particlesToShowers);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, evt, "pandoracalipidSCE");
  if (!trackPIDAssn.isValid()){
    std::cout << "[NumuCCana::getVertex] Event failed: PID is invalid" << std::endl;
    counter_trackPIDAssn++;
  }
  //check the neutrino information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  NuPDG_ = pfnu->PdgCode(); // has to be 14
  
  art::Handle< std::vector<recob::PFParticle> > theParticles;
  evt.getByLabel(m_pfp_producer, theParticles);
  // get Associated flash and then the flash score
  art::FindManyP<anab::T0> nuFlashScoreAsso(theParticles, evt, "flashmatch");
  const std::vector<art::Ptr<anab::T0> > T0_v = nuFlashScoreAsso.at( pfnu.key() );
  if(T0_v.size()==1){
    FlashScoreTime_ = T0_v.at(0)->Time();
    FlashScore_ = T0_v.at(0)->TriggerConfidence();
    if(verbose_!=0) std::cout << "[NumuCCana::getVertex] got flash score" << std::endl;
  }
  else std::cout << "[NumuCCana::getVertex] Flash score invalid" << std::endl;
  
  // get the metadata
  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  //check if there is a neutrino vertex reconstructed and only one
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1){
    std::cout << "[NumuCCana::getVertex] Event failed: Neutrino association failed" << std::endl;
    counter_neutrino_metadata_vec++;
    return false;
  }
  else{
    //check the topological score ov the event
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    NuScore_ = neutrino_properties.at("NuScore"); // nuscore $$
    if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Got NuScore: " << NuScore_ << std::endl;
    if(verbose_ !=0) for (auto& t : neutrino_properties) std::cout << "[NumuCCana::getVertex] neutrino properties: " << t.first << " " << t.second << std::endl;
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    Nu_Vx_ = neutrino_vtx.X();
    Nu_Vy_ = neutrino_vtx.Y();
    Nu_Vz_ = neutrino_vtx.Z();
    auto NuVtx_offset = SCE->GetCalPosOffsets(geo::Point_t(neutrino_vtx.X(), neutrino_vtx.Y(), neutrino_vtx.Z()), 0);
    Nu_Vx_sce_ = (neutrino_vtx.X() - NuVtx_offset.X() );
    Nu_Vy_sce_ = (neutrino_vtx.Y() + NuVtx_offset.Y() );
    Nu_Vz_sce_ = (neutrino_vtx.Z() + NuVtx_offset.Z() );
  }
  if (!is_data_){
    MCNu_vertexDistance = pandoraInterfaceHelper.Distance3D(MCNu_VxSce, MCNu_VySce, MCNu_VzSce, MCNu_Vx, MCNu_Vx, MCNu_Vx);
  }
  //now get the muon candidate track for further investigtion (aka the longest track)
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  double max_track_length = 0; //max muon length
  int muon_pfp_key = -1; // which pfp has longest track
  int pfp_counter = 0; // pfp counter
  // Implement here smart muon track id selection
  NuShowers_ = 0;
  NuTracks_ = 0;
  art::FindMany<anab::T0> trk_T0_assn_v(trackHandle, evt, data_label_crtT0asso_);
  art::FindMany<crt::CRTHit> trk_crt_assn_v(trackHandle, evt, data_label_crtT0asso_);
  //Hits per track
  art::FindManyP<recob::Hit> hits_per_track(trackHandle,evt,Hits_TrackAssLabel);
  for (auto const pfp : pfdaughters){  //loop over all daughter pfparticles
    if (!pfp->IsPrimary()){ 
      pfp_counter++;
      if(verbose_!=0) std::cout << "[NumuCCana::getVertex] looking at pfparticle" << std::endl;
      const lar_pandora::ClusterVector cluster_vec = particlesToClusters.at(pfp);
      std::vector<uint> nHits;
      std::vector<float> pfenergy;
      energyHelper.energy_from_hits(cluster_vec, nHits, pfenergy);
      Nu_NhitsU_ += nHits[0];
      Nu_NhitsV_ += nHits[1];
      Nu_NhitsY_ += nHits[2];
      Nu_CaloU_ += pfenergy[0];
      Nu_CaloV_ += pfenergy[1];
      Nu_CaloY_ += pfenergy[2];
      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
      TrackScoreGlobal_.push_back(pfp_properties.at("TrackScore"));
      int isShowerTrack = 0;
      
      if (particlesToTracks.find(pfp) != particlesToTracks.end() ){ // get the track like pfp
        isShowerTrack+=100;
        NuTracks_++;
        TrackPfp_.push_back(pfp->Self());
        if(verbose_!=0) std::cout << "[NumuCCana::getVertex] pfparticle is track like" << std::endl;
        TrackNHitsU_.push_back(nHits[0]);
        TrackNHitsV_.push_back(nHits[1]);
        TrackNHitsY_.push_back(nHits[2]);
        TrackCaloU_.push_back(pfenergy[0]);
        TrackCaloV_.push_back(pfenergy[1]);
        TrackCaloY_.push_back(pfenergy[2]);
        
        const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
        track_key_.push_back(this_track.key());
        
        TrackScore_.push_back(pfp_properties.at("TrackScore"));
        
        const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
        Vx_.push_back(pfp_vtx.X());//  vtx_distance < 5cm $$
        Vy_.push_back(pfp_vtx.Y());
        Vz_.push_back(pfp_vtx.Z());
        VtxDistance_.push_back(pandoraInterfaceHelper.Distance3D(pfp_vtx.X(), pfp_vtx.Y(), pfp_vtx.Z(), Nu_Vx_, Nu_Vy_, Nu_Vz_));
        
        auto Vtx_offset = SCE->GetCalPosOffsets(geo::Point_t(pfp_vtx.X(), pfp_vtx.Y(), pfp_vtx.Z()), 0);
        Vx_sce_.push_back(pfp_vtx.X() - Vtx_offset.X() );
        Vy_sce_.push_back(pfp_vtx.Y() + Vtx_offset.Y() );
        Vz_sce_.push_back(pfp_vtx.Z() + Vtx_offset.Z() );
        VtxDistance_sce_.push_back(pandoraInterfaceHelper.Distance3D(Vx_sce_.back(), Vy_sce_.back(), Vz_sce_.back(), Nu_Vx_sce_, Nu_Vy_sce_, Nu_Vz_sce_));
        
         const recob::Track::Vector_t &track_dir = this_track->StartDirection();
        //check track length
        TrackPID_.push_back(this_track->ParticleId());
        TrackLength_.push_back(this_track->Length());
        TrackDir_x_.push_back(track_dir.X());
        TrackDir_y_.push_back(track_dir.Y());
        TrackDir_z_.push_back(track_dir.Z());
        TrackStart_x_.push_back(this_track->Start().X());
        TrackStart_y_.push_back(this_track->Start().Y());
        TrackStart_z_.push_back(this_track->Start().Z());
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(this_track->Start().X(), this_track->Start().Y(), this_track->Start().Z()), 0);
        TrackStart_x_sce_.push_back(this_track->Start().X() - Trk_start_offset.X() );
        TrackStart_y_sce_.push_back(this_track->Start().Y() + Trk_start_offset.Y() );
        TrackStart_z_sce_.push_back(this_track->Start().Z() + Trk_start_offset.Z() );
        TrackEnd_x_.push_back(this_track->End().X());
        TrackEnd_y_.push_back(this_track->End().Y());
        TrackEnd_z_.push_back(this_track->End().Z());
        auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(this_track->End().X(), this_track->End().Y(), this_track->End().Z()), 0);
        TrackEnd_x_sce_.push_back(this_track->End().X() - Trk_end_offset.X() );
        TrackEnd_y_sce_.push_back(this_track->End().Y() + Trk_end_offset.Y() );
        TrackEnd_z_sce_.push_back(this_track->End().Z() + Trk_end_offset.Z() );
        TrackTheta_.push_back(this_track->Theta());
        TrackPhi_.push_back(this_track->Phi());
        
        TVector3 Trk_start_SCEcorr;
        TVector3 Trk_end_SCEcorr;
        Trk_start_SCEcorr.SetX(TrackStart_x_sce_.back());
        Trk_start_SCEcorr.SetY(TrackStart_y_sce_.back());
        Trk_start_SCEcorr.SetZ(TrackStart_z_sce_.back());
        Trk_end_SCEcorr.SetX(TrackEnd_x_sce_.back());
        Trk_end_SCEcorr.SetY(TrackEnd_y_sce_.back());
        Trk_end_SCEcorr.SetZ(TrackEnd_z_sce_.back());
        
        //get PID
        PID pid;
        pid.Chi2(trackPIDAssn,this_track, Trk_start_SCEcorr, Trk_end_SCEcorr);
        TrackPID_chiproton_.push_back(pid.PID_Chi2P_3pl);  // chi squqre cuts
        TrackPID_chimuon_.push_back(pid.PID_Chi2Mu_3pl);
        TrackPID_chipion_.push_back(pid.PID_Chi2Pi_3pl);
        TrackPID_chikaon_.push_back(pid.PID_Chi2K_3pl);
        std::cout << "[NumuCCana::getVertex] Got track PID (muon 3pl)  : " << pid.PID_Chi2Mu_3pl << std::endl;
        std::cout << "[NumuCCana::getVertex] Got track PID (proton 3pl): " << pid.PID_Chi2P_3pl << std::endl;
        // get track MCS fit results
        const recob::MCSFitResult &mcsMu = MCSMu_handle->at(this_track.key());
        TrackMomMCS_mom_.push_back(mcsMu.fwdMomentum());
        TrackMomMCS_err_.push_back(mcsMu.fwdMomUncertainty());
        TrackMomMCS_ll_.push_back(mcsMu.fwdLogLikelihood());
        // get track momentum by range results
        float track_mom_range_p, track_mom_range_mu;
        trackHelper.getRangeMomentum(this_track->Length(), track_mom_range_p, track_mom_range_mu);
        TrackMomRange_p_.push_back(track_mom_range_p);
        TrackMomRange_mu_.push_back(track_mom_range_mu);
        
        // get track PID chi2 results (old)
        //std::map<std::string, float> pid_map;
        //if(trackHelper.getPID(pid_map, this_track, trackPIDAssn)){
          //TrackPID_chiproton_.push_back(pid_map.at("chi2_proton"));  // chi squqre cuts
          //TrackPID_chimuon_.push_back(pid_map.at("chi2_muon"));
          //TrackPID_chipion_.push_back(pid_map.at("chi2_pion"));
          //TrackPID_chikaon_.push_back(pid_map.at("chi2_kaon"));
        //}
        
        // get the maximal length of all tracks for later
        if( this_track->Length() > max_track_length){ //take the longest track as muon candidate
          max_track_length = this_track->Length();
          muon_pfp_key = pfp_counter;
        }
        // for Overlay: get the back tracked MC particles
        if(!is_data_){
          std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(this_track.key());
          BackTrackerTruthMatch backtrackertruthmatch;
          art::Handle<std::vector<recob::Hit> > Handle_Hit;
          evt.getByLabel("gaushit", Handle_Hit);
          backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
          auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
          if(MCparticle){
            // MC_granddau_pdg stores all granddaughters of all the primary pfparticle(s)
           // MC_granddau_pdg.push_back(MCparticle->PdgCode());
            std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully with PDG: " << MCparticle->PdgCode() << std::endl;
            TVector3 TrueTrackPos(MCparticle->Px(), MCparticle->Py(), MCparticle->Pz());
            TVector3 true_start(MCparticle->Position().X(), MCparticle->Position().Y(), MCparticle->Position().Z());
            TVector3 true_end(MCparticle->EndPosition().X(), MCparticle->EndPosition().Y(), MCparticle->EndPosition().Z());
            
            MCTrackPDG_.push_back(MCparticle->PdgCode());
            MCTrackEnergy_.push_back(-1);
            MCTrackMomentum_.push_back(MCparticle->P());
            MCTrackTheta_.push_back(TrueTrackPos.Theta());
            MCTrackPhi_.push_back(TrueTrackPos.Phi());
            MCTrackLength_.push_back((true_start - true_end).Mag());

            MCTrackStart_x_.push_back(MCparticle->Position().X());
            MCTrackStart_y_.push_back(MCparticle->Position().Y());
            MCTrackStart_z_.push_back(MCparticle->Position().Z());
            MCTrackEnd_x_.push_back(MCparticle->EndPosition().X());
            MCTrackEnd_y_.push_back(MCparticle->EndPosition().Y());
            MCTrackEnd_z_.push_back(MCparticle->EndPosition().Y());
          }
          else{
            MCTrackPDG_.push_back(-1);
            MCTrackEnergy_.push_back(-1);
            MCTrackMomentum_.push_back(-1);
            MCTrackTheta_.push_back(-1);
            MCTrackPhi_.push_back(-1);
            MCTrackLength_.push_back(-1);

            MCTrackStart_x_.push_back(-1);
            MCTrackStart_y_.push_back(-1);
            MCTrackStart_z_.push_back(-1);
            MCTrackEnd_x_.push_back(-1);
            MCTrackEnd_y_.push_back(-1);
            MCTrackEnd_z_.push_back(-1);
          }
          MCTrackPurity_.push_back( backtrackertruthmatch.ReturnPurity() );
          std::cout << "[NumuCCana::getVertex] Matched MC daughter with purity: " << backtrackertruthmatch.ReturnPurity() << std::endl;
        }
        // get associated CRT T0 if any
        const std::vector<const anab::T0*>& T0crt = trk_T0_assn_v.at(this_track.key() );
        if(T0crt.size()!=0){
          if(verbose_!=0){
            std::cout << "################################################" << std::endl;
            std::cout << "Found T0 object from crthit - track association:" << std::endl;
            std::cout << "## Time [us]:\t\t" << T0crt.at(0)->Time() << std::endl;
            std::cout << "## Is from CRT []:\t" << T0crt.at(0)->TriggerType() << std::endl;
            std::cout << "## DCA [cm]:\t\t" << T0crt.at(0)->TriggerConfidence() << std::endl;
            std::cout << "## plane:\t\t" << T0crt.at(0)->TriggerBits() << std::endl;
            std::cout << "################################################" << std::endl;
          }
          //fill tree variables
          crtt0_time_.push_back(T0crt.at(0)->Time());
          crtt0_trig_.push_back(T0crt.at(0)->TriggerType());
          crtt0_DCA_.push_back(T0crt.at(0)->TriggerConfidence());
          crtt0_plane_.push_back(T0crt.at(0)->TriggerBits());
        }
        else{
          crtt0_time_.push_back(-1);
          crtt0_trig_.push_back(-1);
          crtt0_DCA_.push_back(-1);
          crtt0_plane_.push_back(-1);
        }
        // get the corresponding CRT hit itself
        const std::vector<const crt::CRTHit*>& CRTHit_v = trk_crt_assn_v.at(this_track.key());
        if(CRTHit_v.size()!=0){
          a_crthit_ts0.push_back(CRTHit_v.at(0)->ts0_ns);
          a_crthit_ts1.push_back(CRTHit_v.at(0)->ts1_ns);
          a_crthit_s.push_back(CRTHit_v.at(0)->ts0_s);
          a_crthit_x_.push_back(CRTHit_v.at(0)->x_pos);
          a_crthit_y_.push_back(CRTHit_v.at(0)->y_pos);
          a_crthit_z_.push_back(CRTHit_v.at(0)->z_pos);
          a_adc_length.push_back(CRTHit_v.at(0)->pesmap.begin()->second.size());
          a_crt_adc.push_back(CRTHit_v.at(0)->peshit);
          a_crt_plane_.push_back(CRTHit_v.at(0)->plane);
          a_t0_counter.push_back(CRTHit_v.size());
        }
        else{
          a_crthit_ts0.push_back(-1);
          a_crthit_ts1.push_back(-1);
          a_crthit_s.push_back(-1);
          a_adc_length.push_back(-1);
          a_crt_adc.push_back(-1);
          a_t0_counter.push_back(0);
        }
        
        if (particlesToSpacePoints.find(pfp) == particlesToSpacePoints.end())
        {
          // If a daughter has no associated spacepoints, count the hits to contribute to the total, but dont save the daughter
          std::cout << "[NuCCana::getVertex] Daughter had no associated spacepoints." << std::endl;
          //return false;
        }
        for (unsigned int p = 0; p < this_track->NumberTrajectoryPoints(); ++p)
        {
            const auto position(this_track->LocationAtPoint(p));
            AllTrack_point_x_.push_back(position.x());
            AllTrack_point_y_.push_back(position.y());
            AllTrack_point_z_.push_back(position.z());
        }
      }
      if(verbose_!=0) std::cout << "[NumuCCana::getVertex] End of loop over track like particles" << std::endl;
      
      if (particlesToShowers.find(pfp) != particlesToShowers.end()){
        isShowerTrack+=10;
        NuShowers_++;
        ShowerScore_.push_back(pfp_properties.at("TrackScore"));
        ShowerPfp_.push_back(pfp->Self());
        ShowerNHitsU_.push_back(nHits[0]);
        ShowerNHitsV_.push_back(nHits[1]);
        ShowerNHitsY_.push_back(nHits[2]);
        ShowerCaloU_.push_back(pfenergy[0]);
        ShowerCaloV_.push_back(pfenergy[1]);
        ShowerCaloY_.push_back(pfenergy[2]);
        
        if(verbose_!=0) std::cout << "[NumuCCana::getVertex] pfparticle is shower like, Shower nr: " << NuShowers_ << std::endl;
        const art::Ptr<recob::Shower> this_shower = particlesToShowers.at(pfp).front();
        
        if (this_shower->has_length() && this_shower->has_open_angle() )
        {
          if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Fill Shower variables" << std::endl;
          const TVector3 &shower_dir = this_shower->Direction();
          ShowerLength_.push_back(this_shower->Length());
          ShowerOpenAngle_.push_back(this_shower->OpenAngle());
          //auto energy_tmp = this_shower->Energy();//.at(this_shower->best_plane());
          //std::cout << "Energy Size:  Best plane: " << this_shower->best_plane() << std::endl;
          //if(abs(this_shower->best_plane()) < 4) ShowerEnergy_.push_back(this_shower->Energy().at(0));
          //if(energy_tmp.size()!=0) ShowerEnergy_.push_back(energy_tmp.at(this_shower->best_plane()));
          //auto energyMIP_tmp = this_shower->MIPEnergy();
          //if(energyMIP_tmp.size()!=0) ShowerMIPEnergy_.push_back(energyMIP_tmp.at(this_shower->best_plane()));
          ShowerDir_x_.push_back(shower_dir.X());
          ShowerDir_y_.push_back(shower_dir.Y());
          ShowerDir_z_.push_back(shower_dir.Z());
          
          if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Get Shower variables Energy" << std::endl;
          std::vector<float> pitches(3, std::numeric_limits<float>::lowest());
          std::vector<float> dqdx(3, std::numeric_limits<float>::lowest());
          std::vector<std::vector<float>> dqdx_hits(3, std::vector<float>());
          if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Call energy helper function" << std::endl;
          //energyHelper.dQdx(shower_dir, cluster_vec, clustersToHits, dqdx, dqdx_hits, pitches);
          //std::vector<float> dedx = energyHelper.dEdx_from_dQdx(dqdx);

          if(verbose_!=0) std::cout << "[NumuCCana::getVertex] Fill Shower variables Energy" << std::endl;
          /*Shower_dEdxU_.push_back(dedx[0]);
          Shower_dEdxV_.push_back(dedx[1]);
          Shower_dEdxY_.push_back(dedx[2]);
          Shower_dEdxHitsU_.push_back(dqdx_hits[0].size());
          Shower_dEdxHitsV_.push_back(dqdx_hits[1].size());
          Shower_dEdxHitsY_.push_back(dqdx_hits[2].size());
          Shower_dEdxPitchU_.push_back(pitches[0]);
          Shower_dEdxPitchV_.push_back(pitches[1]);
          Shower_dEdxPitchY_.push_back(pitches[2]);*/
          
        }
        else{
          std::cout << "[NumuCCana::getVertex] Shower has no length or opening angle!" << std::endl;
        }
      }
      isShowerTrack_.push_back(isShowerTrack);
    }
    
  }
  NumPfp_ = pfp_counter;
  if(verbose_!=0) std::cout << "[NumuCCana::getVertex] looking into the tracks" << std::endl;
  /*
  if(fill_per_track_!=1){
    if(verbose_!=0) std::cout << "[NumuCCana::getVertex] take only one track" << std::endl;
    if(NuTracks_>1){
      if(verbose_!=0) std::cout << "[NumuCCana::getVertex] loop over tracks" << std::endl;
      lar_pandora::PFParticleVector pass_proton_cut;
      for (auto const pfp : pfdaughters){  //loop over all daughter pfparticles
        if (!pfp->IsPrimary()){ 
          if (particlesToTracks.find(pfp) != particlesToTracks.end() ){ // get the track like pfp
            //NuTracks_++;
            double this_chiproton = -1;
            double this_chimuon_ = -1;
            const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
            std::map<std::string, float> pid_map;
            if(trackHelper.getPID(pid_map, this_track, trackPIDAssn)){
              this_chiproton = pid_map.at("chi2_proton"); 
              this_chimuon_ = pid_map.at("chi2_muon");
              
            }
            if(this_chimuon_!=-1 && (this_chimuon_/this_chiproton)<0.168 ){
              pass_proton_cut.push_back(pfp);
            }
          }
        }
      }
      if(pass_proton_cut.size()>1){
        lar_pandora::PFParticleVector pass_pion_cut;
        for (auto const pfp : pass_proton_cut){  //loop over all daughter pfparticles
          if (!pfp->IsPrimary()){ 
            if (particlesToTracks.find(pfp) != particlesToTracks.end() ){ // get the track like pfp
              //NuTracks_++;
              double this_chipion = -1;
              double this_chimuon_ = -1;
              const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
              std::map<std::string, float> pid_map;
              if(trackHelper.getPID(pid_map, this_track, trackPIDAssn)){
                this_chipion = pid_map.at("chi2_pion"); 
                this_chimuon_ = pid_map.at("chi2_muon");

              }
              if(this_chimuon_!=-1 && (this_chimuon_/this_chipion)<1.06 ){
                pass_pion_cut.push_back(pfp);
              }
            }
          }
        }
        if(pass_pion_cut.size()>1){// take longest here
          double max_length = 0;
          int longest_index = 0;
          int longest_counter = 0;
          for (auto const pfp : pass_pion_cut){  //loop over all daughter pfparticles
            const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
            if( this_track->Length() < max_length){
              max_length = this_track->Length();
              longest_index = longest_counter;
            }
            longest_counter++;
          }
          if (!GetMuon(pass_pion_cut.at(longest_index), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
            std::cout << "[NumuCCana::getVertex] Event failed: Daughter investigation unsuccessful" << std::endl;
            return false;
          }
          else if (1){//MatchDaughter(evt, pass_pion_cut.at(longest_index))){
            std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully" << std::endl;
            //fill tree
           // my_event_->Fill();
            return true;
          }
          else return false;
        }
        else if(pass_pion_cut.size()==1){ // take this track
          if (!GetMuon(pass_pion_cut.at(0), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
            std::cout << "[NumuCCana::getVertex] Event failed: Daughter investigation unsuccessful" << std::endl;
            return false;
          }
          else if (1){//MatchDaughter(evt, pass_pion_cut.at(0))){
            std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully" << std::endl;
            //fill tree
            //my_event_->Fill();
            return true;
          }
          else return false;
        }
        else{ //take longest from proton cut
          double max_length = 0;
          int longest_index = 0;
          int longest_counter = 0;
          for (auto const pfp : pass_proton_cut){  //loop over all daughter pfparticles
            const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
            if( this_track->Length() < max_length){
              max_length = this_track->Length();
              longest_index = longest_counter;
            }
            longest_counter++;
          }
          if (!GetMuon(pass_proton_cut.at(longest_index), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
            std::cout << "[NumuCCana::getVertex] Event failed: Daughter investigation unsuccessful" << std::endl;
            return false;
          }
          else if (1){//MatchDaughter(evt, pass_proton_cut.at(longest_index))){
            std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully" << std::endl;
            //fill tree
            //my_event_->Fill();
            return true;
          }
          else return false;
        }
      }
      else if(pass_proton_cut.size()==1){
        if (!GetMuon(pfdaughters.at(muon_pfp_key), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
          std::cout << "[NumuCCana::getVertex] Event failed: Daughter investigation unsuccessful" << std::endl;
          return false;
        }
        else if (1){//MatchDaughter(evt, pfdaughters.at(muon_pfp_key))){
            std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully" << std::endl;
            //fill tree
            //my_event_->Fill();
            return true;
          }
          else return false;
      }
    }
    else if(NuTracks_ == 1){
      if (!GetMuon(pfdaughters.at(muon_pfp_key), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
        std::cout << "[NumuCCana::getVertex] Event failed: Daughter investigation unsuccessful" << std::endl;
        return false;
      }
      else if (1){//MatchDaughter(evt, pfdaughters.at(muon_pfp_key))){
        std::cout << "[NumuCCana::getVertex] Matched MC daughter successfully" << std::endl;
        //fill tree
        //my_event_->Fill();
        return true;
      }
      else return false;
    }
    else{ // No track like object
      std::cout << "[NumuCCana::getVertex] Event failed: No Daughter is muon/tracklike" << std::endl;
      counter_moun++;
      return false;
    }
  }*/
  if(verbose_!=0) std::cout << "blabla" << muon_pfp_key << std::endl;
  return true;
}

bool NumuCCana::GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                         const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                         const art::FindManyP<anab::ParticleID> &trackPIDAssn){
  //check the track score value
  if(verbose_!=0) std::cout << "[NumuCCana::GetMuon] Start with track" << std::endl;

  if (particlesToTracks.find(pfp) != particlesToTracks.end()){ // always true
    // get the recob track
    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
    muon_candidate_key = this_track.key();
    muon_candidate_pfp = pfp->Self();
  }
  // all requirements passed, to take this event with this muon candidate
  if(verbose_!=0) std::cout << "[NumuCCana::GetMuon] End with track" << std::endl;
  return true;
}

void NumuCCana::clearEvent(){ // reset all pfp related vectors
    pfparticles.clear();
    pfneutrinos.clear();
    pfdaughters.clear();
    pftracks.clear();
    particleMap.clear();
    particlesToMetadata.clear();
    particlesToVertices.clear();
    particlesToTracks.clear();
    particlesToClusters.clear();
    particlesToSpacePoints.clear();
  //mc
  matchedParticles.clear();
  matchedMCParticles.clear();
  
  MCTruthCollection.clear();
  GTruthCollection.clear();
  MCParticleCollection.clear();
  WeightCollection.clear();
  WeightCollection_G4.clear();
  
  EventWeight = 1;
  
  has_neutrino_ = -1;
  track_key_.clear();
  NuScore_ = -1;
  FlashScore_ = -1;
  FlashScoreTime_ = -1;
  NuPDG_ = 0;
  NumPfp_ = -1;
  Nu_Vx_=-999,Nu_Vy_=-999,Nu_Vz_=-999;
  Nu_Vx_sce_=-999,Nu_Vy_sce_=-999,Nu_Vz_sce_=-999;
  NuTracks_ = -1;
  NuShowers_ = -1;
  muon_candidate_key = -1;
  muon_candidate_pfp = -1;
  
  Nu_NhitsU_ = 0, Nu_NhitsV_ = 0, Nu_NhitsY_ = 0;
  Nu_CaloU_ = 0, Nu_CaloV_ = 0, Nu_CaloY_ = 0;
  
  Vx_.clear();
  Vy_.clear();
  Vz_.clear();
  Vx_sce_.clear();
  Vy_sce_.clear();
  Vz_sce_.clear();
  VtxDistance_.clear();
  VtxDistance_sce_.clear();
  TrackScore_.clear();
  ShowerScore_.clear();
  TrackScoreGlobal_.clear();
  TrackLength_.clear();
  TrackPID_.clear();
  
  TrackPfp_.clear();
  ShowerPfp_.clear();
  isShowerTrack_.clear();
  
  TrackStart_x_.clear();
  TrackStart_y_.clear();
  TrackStart_z_.clear();
  TrackStart_x_sce_.clear();
  TrackStart_y_sce_.clear();
  TrackStart_z_sce_.clear();
  TrackEnd_x_.clear();
  TrackEnd_y_.clear();
  TrackEnd_z_.clear();
  TrackEnd_x_sce_.clear();
  TrackEnd_y_sce_.clear();
  TrackEnd_z_sce_.clear();
  TrackDir_x_.clear();
  TrackDir_y_.clear();
  TrackDir_z_.clear();
  TrackTheta_.clear();
  TrackPhi_.clear();
  
  MCTrackPurity_.clear();
  MCTrackPDG_.clear();
  MCTrackEnergy_.clear();
  MCTrackMomentum_.clear();
  MCTrackTheta_.clear();
  MCTrackPhi_.clear();
  MCTrackLength_.clear();

  MCTrackStart_x_.clear();
  MCTrackStart_y_.clear();
  MCTrackStart_z_.clear();
  MCTrackEnd_x_.clear();
  MCTrackEnd_y_.clear();
  MCTrackEnd_z_.clear();
  
  SpacePoint_x_.clear();
  SpacePoint_y_.clear();
  SpacePoint_z_.clear();
  
  Wire_id_.clear();
  Wire_plane_.clear();
  
  TrackMomRange_p_.clear();
  TrackMomRange_mu_.clear();
  TrackMomMCS_mom_.clear();
  TrackMomMCS_err_.clear();
  TrackMomMCS_ll_.clear();
  
  TrackPID_chiproton_.clear();
  TrackPID_chimuon_.clear();
  TrackPID_chipion_.clear();
  TrackPID_chikaon_.clear();
  
  TrackNHitsU_.clear();
  TrackNHitsV_.clear();
  TrackNHitsY_.clear();
  TrackCaloU_.clear();
  TrackCaloV_.clear();
  TrackCaloY_.clear();
  
  AllTrack_point_x_.clear();
  AllTrack_point_y_.clear();
  AllTrack_point_z_.clear();


  ShowerLength_.clear();
  ShowerOpenAngle_.clear();
  ShowerEnergy_.clear();
  ShowerMIPEnergy_.clear();
  ShowerDir_x_.clear();
  ShowerDir_y_.clear();
  ShowerDir_z_.clear();

  Shower_dEdxU_.clear();
  Shower_dEdxV_.clear();
  Shower_dEdxY_.clear();
  Shower_dEdxHitsU_.clear();
  Shower_dEdxHitsV_.clear();
  Shower_dEdxHitsY_.clear();
  Shower_dEdxPitchU_.clear();
  Shower_dEdxPitchV_.clear();
  Shower_dEdxPitchY_.clear();

  ShowerNHitsU_.clear();
  ShowerNHitsV_.clear();
  ShowerNHitsY_.clear();
  ShowerCaloU_.clear();
  ShowerCaloV_.clear();
  ShowerCaloY_.clear();
  
  crtt0_time_.clear();
  crtt0_trig_.clear();
  crtt0_DCA_.clear();
  crtt0_plane_.clear();
  
  a_crthit_ts0.clear();
  a_crthit_ts1.clear();
  a_crthit_s.clear();
  a_crthit_x_.clear();
  a_crthit_y_.clear();
  a_crthit_z_.clear();
  a_crt_plane_.clear();
  a_adc_length.clear();
  a_crt_adc.clear();
  a_t0_counter.clear();
  
  nr_crthit_ = -1; // # crt hits assigned to a tpc track
  crthit_ts0_.clear();
  crthit_ts1_.clear();
  crthit_s_.clear();
  crthit_x_.clear();
  crthit_y_.clear();
  crthit_z_.clear();
  crt_plane_.clear();
  adc_length_.clear();
  crt_adc_.clear();
  
  TimFla_ = -99;
  flash_PE_ = -99;
  flash_y_ = -999;
  flash_z_ = -999;
  
  crt_trig_corr_mean = 0;
  crt_trig_corr_med = 0;
  
  if(!is_data_){
  NuMCnu = -1; 
  MCNu_Interaction = -1;
  MCNu_CCNC = -1;
  MCNu_PDG = -1;
  MCNu_Energy = -1;
  MCNu_leptonPx = -9, MCNu_leptonPy = -9, MCNu_leptonPz = -9;
  MCNu_LeptonEnergy = -9;
  MCNu_Px = -9, MCNu_Py = -9, MCNu_Pz = -9;
  MCNu_leptonTheta = -9;
  MCNu_time = -99; 
  MCNu_Vx = -999, MCNu_Vy = -999, MCNu_Vz = -999;
  MCNu_VxSce = -999, MCNu_VySce = -999, MCNu_VzSce = -999;
  MCNu_vertexDistance = -9;
    
  Genie_Q2 = -1;
  Genie_q2 = -1;
  Genie_W = -1;
  Genie_T = -1;
  Genie_X = -1;
  Genie_Y = -1;
  Genie_nNeutron_preFSI = -1;
  Genie_nProton_preFSI = -1;
  Genie_nPi0_preFSI = -1;
  Genie_nPiPlus_preFSI = -1; 
  Genie_nPiMinus_preFSI = -1;

  MCNU_matched = false;
  MCCosmic_matched = false;
    
  MCle_pfp_.clear();
  MCle_key_.clear();
  MCle_PDG_.clear();
  MCle_Energy_.clear();
  MCle_Px_.clear();
  MCle_Py_.clear();
  MCle_Pz_.clear();
  MCle_Vx_.clear();
  MCle_Vy_.clear(); 
  MCle_Vz_.clear();
  MCle_Endx_.clear();
  MCle_Endy_.clear();
  MCle_Endz_.clear();
  MCle_length_.clear();
  MCle_VxSce_.clear();
  MCle_VySce_.clear();
  MCle_VzSce_.clear();
  MCle_Theta_.clear();
  MCle_Phi_.clear();
  }
}

void NumuCCana::initialize_tree()
{
  // Implementation of required member function here.
  std::cout << "Initialize variables and histograms for event tree" << std::endl;
  //tree stuff for hits: //////////////////////////////////////////////////////////////////////////////////
  my_event_ = tfs->make<TTree>("event","numuCC event tree");
  my_event_->Branch("event_counter", &event_counter, "event_counter/I");
  my_event_->Branch("frunNum", &frunNum, "Run Number/i");
  my_event_->Branch("fsubRunNum", &fsubRunNum, "SubRun Number/i");
  my_event_->Branch("fEvtNum", &fEvtNum, "Event Number/i");
  my_event_->Branch("EventWeight", &EventWeight, "EventWeight/D");
  
  my_event_->Branch("TriTim_sec",        &TriTim_sec_,            "TriTim_sec/D");
  my_event_->Branch("TriTim_nsec",       &TriTim_nsec_,           "TriTim_nsec/D");
  
  my_event_->Branch("NuScore",           &NuScore_,               "NuScore/D");
  my_event_->Branch("FlashScore",        &FlashScore_,            "FlashScore/D");
  my_event_->Branch("FlashScoreTime",    &FlashScoreTime_,        "FlashScoreTime/D");
  //my_event_->Branch("NuPDG",             &NuPDG_,                 "NuPDG/I");
  my_event_->Branch("NumPfp",            &NumPfp_,                "NumPfp/I");
  my_event_->Branch("Nu_Vx",             &Nu_Vx_,                 "Nu_Vx/D");
  my_event_->Branch("Nu_Vy",             &Nu_Vy_,                 "Nu_Vy/D");
  my_event_->Branch("Nu_Vz",             &Nu_Vz_,                 "Nu_Vz/D");
  my_event_->Branch("Nu_Vx_sce",             &Nu_Vx_sce_,                 "Nu_Vx_sce/D");
  my_event_->Branch("Nu_Vy_sce",             &Nu_Vy_sce_,                 "Nu_Vy_sce/D");
  my_event_->Branch("Nu_Vz_sce",             &Nu_Vz_sce_,                 "Nu_Vz_sce/D");
  
  my_event_->Branch("NuTracks",          &NuTracks_,              "NuTracks/I");
  //my_event_->Branch("NuShowers",         &NuShowers_,             "NuShowers/I");
  /* avoid too much information
  my_event_->Branch("Nu_NhitsU",         &Nu_NhitsU_,             "Nu_NhitsU/I");
  my_event_->Branch("Nu_NhitsV",         &Nu_NhitsV_,             "Nu_NhitsV/I");
  my_event_->Branch("Nu_NhitsY",         &Nu_NhitsY_,             "Nu_NhitsY/I");
  my_event_->Branch("Nu_CaloU",         &Nu_CaloU_,             "Nu_CaloU/F");
  my_event_->Branch("Nu_CaloV",         &Nu_CaloV_,             "Nu_CaloV/F");
  my_event_->Branch("Nu_CaloY",         &Nu_CaloY_,             "Nu_CaloY/F");
  */
  my_event_->Branch("muon_candidate_key",&muon_candidate_key,     "muon_candidate_key/I");
  my_event_->Branch("muon_candidate_pfp",&muon_candidate_pfp,     "muon_candidate_pfp/I");
  
  my_event_->Branch("track_key",        &track_key_);
  my_event_->Branch("TrackPID",        &TrackPID_);
  
  my_event_->Branch("Vx",                &Vx_);
  my_event_->Branch("Vy",                &Vy_);
  my_event_->Branch("Vz",                &Vz_);
  my_event_->Branch("Vx_sce",                &Vx_sce_);
  my_event_->Branch("Vy_sce",                &Vy_sce_);
  my_event_->Branch("Vz_sce",                &Vz_sce_);
  
  my_event_->Branch("TrackScore",        &TrackScore_);
  //my_event_->Branch("ShowerScore",        &ShowerScore_); // avoid too much information
  my_event_->Branch("TrackScoreGlobal",        &TrackScoreGlobal_);
  my_event_->Branch("VtxDistance",       &VtxDistance_);
  my_event_->Branch("VtxDistance_sce",       &VtxDistance_sce_);
  my_event_->Branch("TrackLength",       &TrackLength_);
  my_event_->Branch("TrackPfp",       &TrackPfp_);
  //my_event_->Branch("isShowerTrack",       &isShowerTrack_); //avoid too much information
  
  my_event_->Branch("TrackMomRange_p",   &TrackMomRange_p_);
  my_event_->Branch("TrackMomRange_mu",  &TrackMomRange_mu_);
  my_event_->Branch("TrackMomMCS_mom",   &TrackMomMCS_mom_);
  my_event_->Branch("TrackMomMCS_err",   &TrackMomMCS_err_);
  my_event_->Branch("TrackMomMCS_ll",    &TrackMomMCS_ll_);
  
  my_event_->Branch("TrackStart_x",      &TrackStart_x_);
  my_event_->Branch("TrackStart_y",      &TrackStart_y_);
  my_event_->Branch("TrackStart_z",      &TrackStart_z_);
  my_event_->Branch("TrackStart_x_sce",      &TrackStart_x_sce_);
  my_event_->Branch("TrackStart_y_sce",      &TrackStart_y_sce_);
  my_event_->Branch("TrackStart_z_sce",      &TrackStart_z_sce_);
  my_event_->Branch("TrackEnd_x",        &TrackEnd_x_);
  my_event_->Branch("TrackEnd_y",        &TrackEnd_y_);
  my_event_->Branch("TrackEnd_z",        &TrackEnd_z_);
  my_event_->Branch("TrackEnd_x_sce",        &TrackEnd_x_sce_);
  my_event_->Branch("TrackEnd_y_sce",        &TrackEnd_y_sce_);
  my_event_->Branch("TrackEnd_z_sce",        &TrackEnd_z_sce_);
  my_event_->Branch("TrackDir_x",        &TrackDir_x_);
  my_event_->Branch("TrackDir_y",        &TrackDir_y_);
  my_event_->Branch("TrackDir_z",        &TrackDir_z_);
  my_event_->Branch("TrackTheta",        &TrackTheta_);
  my_event_->Branch("TrackPhi",          &TrackPhi_);
  
  /* Avoid getting too much information
  my_event_->Branch("SpacePoint_x",          &SpacePoint_x_);
  my_event_->Branch("SpacePoint_y",          &SpacePoint_y_);
  my_event_->Branch("SpacePoint_z",          &SpacePoint_z_);
  
  my_event_->Branch("Wire_id",          &Wire_id_);
  my_event_->Branch("Wire_plane",          &Wire_plane_);
  
  my_event_->Branch("AllTrack_point_x",          &AllTrack_point_x_);
  my_event_->Branch("AllTrack_point_y",          &AllTrack_point_y_);
  my_event_->Branch("AllTrack_point_z",          &AllTrack_point_z_);
  
  my_event_->Branch("TrackNHitsU",       &TrackNHitsU_);
  my_event_->Branch("TrackNHitsV",       &TrackNHitsV_);
  my_event_->Branch("TrackNHitsY",       &TrackNHitsY_);
  my_event_->Branch("TrackCaloU",        &TrackCaloU_);
  my_event_->Branch("TrackCaloV",        &TrackCaloV_);
  my_event_->Branch("TrackCaloY",        &TrackCaloY_);
  
  my_event_->Branch("ShowerPfp",      &ShowerPfp_);
  my_event_->Branch("ShowerLength",      &ShowerLength_);
  my_event_->Branch("ShowerOpenAngle",   &ShowerOpenAngle_);
  my_event_->Branch("ShowerEnergy",   &ShowerEnergy_);
  my_event_->Branch("ShowerMIPEnergy",   &ShowerMIPEnergy_);
  my_event_->Branch("ShowerDir_x",       &ShowerDir_x_);
  my_event_->Branch("ShowerDir_y",       &ShowerDir_y_);
  my_event_->Branch("ShowerDir_z",       &ShowerDir_z_);
  
  my_event_->Branch("Shower_dEdxU",      &Shower_dEdxU_);
  my_event_->Branch("Shower_dEdxV",      &Shower_dEdxV_);
  my_event_->Branch("Shower_dEdxY",      &Shower_dEdxY_);
  
  my_event_->Branch("Shower_dEdxHitsU",  &Shower_dEdxHitsU_);
  my_event_->Branch("Shower_dEdxHitsV",  &Shower_dEdxHitsV_);
  my_event_->Branch("Shower_dEdxHitsY",  &Shower_dEdxHitsY_);
  
  my_event_->Branch("Shower_dEdxPitchU", &Shower_dEdxPitchU_);
  my_event_->Branch("Shower_dEdxPitchV", &Shower_dEdxPitchV_);
  my_event_->Branch("Shower_dEdxPitchY", &Shower_dEdxPitchY_);
  
  my_event_->Branch("ShowerNHitsU",      &ShowerNHitsU_);
  my_event_->Branch("ShowerNHitsV",      &ShowerNHitsV_);
  my_event_->Branch("ShowerNHitsY",      &ShowerNHitsY_);
  
  my_event_->Branch("ShowerCaloU",       &ShowerCaloU_);
  my_event_->Branch("ShowerCaloV",       &ShowerCaloV_);
  my_event_->Branch("ShowerCaloY",       &ShowerCaloY_);
  */
  my_event_->Branch("crtt0_time",        &crtt0_time_);
  my_event_->Branch("crtt0_trig",        &crtt0_trig_);
  my_event_->Branch("crtt0_DCA",         &crtt0_DCA_);
  my_event_->Branch("crtt0_plane",       &crtt0_plane_);
  
  my_event_->Branch("TrackPID_chiproton",&TrackPID_chiproton_);
  my_event_->Branch("TrackPID_chipion",  &TrackPID_chipion_);
  my_event_->Branch("TrackPID_chikaon",  &TrackPID_chikaon_);
  my_event_->Branch("TrackPID_chimuon",  &TrackPID_chimuon_);

  my_event_->Branch("nr_crthit",         &nr_crthit_,              "nr_crthit_/I");
  my_event_->Branch("crthit_ts0",        &crthit_ts0_);
  my_event_->Branch("crthit_ts1",        &crthit_ts1_);
  my_event_->Branch("crthit_s",        &crthit_s_);
  my_event_->Branch("crthit_x",        &crthit_x_);
  my_event_->Branch("crthit_y",        &crthit_y_);
  my_event_->Branch("crthit_z",        &crthit_z_);
  my_event_->Branch("crt_plane",        &crt_plane_);
  my_event_->Branch("adc_length",        &adc_length_);
  my_event_->Branch("crt_adc",           &crt_adc_);
  my_event_->Branch("crtbeam_hit_nr",   &crtbeam_hit_nr_);

  my_event_->Branch("a_crthit_ts0",      &a_crthit_ts0);
  my_event_->Branch("a_crthit_ts1",      &a_crthit_ts1);
  my_event_->Branch("a_crthit_s",      &a_crthit_s);
  my_event_->Branch("a_crthit_x",        &a_crthit_x_);
  my_event_->Branch("a_crthit_y",        &a_crthit_y_);
  my_event_->Branch("a_crthit_z",        &a_crthit_z_);
  my_event_->Branch("a_crt_plane",        &a_crt_plane_);
  my_event_->Branch("a_adc_length",      &a_adc_length);
  my_event_->Branch("a_crt_adc",         &a_crt_adc);
  my_event_->Branch("a_t0_counter",      &a_t0_counter);
  
  /* avoid getting to much infos
  my_event_->Branch("TimFla",            &TimFla_,                "TimFla/D");
  my_event_->Branch("flash_PE",          &flash_PE_,              "flash_PE/D");
  my_event_->Branch("flash_y",           &flash_y_,               "flash_y/D");
  my_event_->Branch("flash_z",           &flash_z_,               "flash_z/D");
  */
  my_event_->Branch("crt_trig_corr_mean",&crt_trig_corr_mean,     "crt_trig_corr_mean/D");
  my_event_->Branch("crt_trig_corr_med", &crt_trig_corr_med,      "crt_trig_corr_med/D");
  
  if(!is_data_){
    my_event_->Branch("MCNu_Interaction",&MCNu_Interaction,"MCNu_Interaction/I");
    my_event_->Branch("MCNu_CCNC",           &MCNu_CCNC,           "MCNu_CCNC/I");
    my_event_->Branch("MCNu_PDG",            &MCNu_PDG,            "MCNu_PDG/I");
    my_event_->Branch("MCNu_Energy",         &MCNu_Energy,         "MCNu_Energy/F");
    my_event_->Branch("MCNu_leptonPx",       &MCNu_leptonPx,       "MCNu_leptonPx/F");
    my_event_->Branch("MCNu_leptonPy",       &MCNu_leptonPy,       "MCNu_leptonPy/F");
    my_event_->Branch("MCNu_leptonPz",       &MCNu_leptonPz,       "MCNu_leptonPz/F");
    my_event_->Branch("MCNu_LeptonEnergy",   &MCNu_LeptonEnergy,   "MCNu_LeptonEnergy/F");
    my_event_->Branch("MCNu_Px",             &MCNu_Px,             "MCNu_Px/F");
    my_event_->Branch("MCNu_Py",             &MCNu_Py,             "MCNu_Py/F");
    my_event_->Branch("MCNu_Pz",             &MCNu_Pz,             "MCNu_Pz/F");
    my_event_->Branch("MCNu_leptonTheta",    &MCNu_leptonTheta,    "MCNu_leptonTheta/F");
    my_event_->Branch("MCNu_time",           &MCNu_time,           "MCNu_time/F");
    my_event_->Branch("MCNu_Vx",             &MCNu_Vx,             "MCNu_Vx/F");
    my_event_->Branch("MCNu_Vy",             &MCNu_Vy,             "MCNu_Vy/F");
    my_event_->Branch("MCNu_Vz",             &MCNu_Vz,             "MCNu_Vz/F");
    my_event_->Branch("MCNu_VxSce",          &MCNu_VxSce,          "MCNu_VxSce/F");
    my_event_->Branch("MCNu_VySce",          &MCNu_VySce,          "MCNu_VySce/F");
    my_event_->Branch("MCNu_VzSce",          &MCNu_VzSce,          "MCNu_VzSce/F");
    //my_event_->Branch("MCNu_Vx_sce",          &MCNu_Vx_sce,          "MCNu_Vx_sce/F");
    //my_event_->Branch("MCNu_Vy_sce",          &MCNu_Vy_sce,          "MCNu_Vy_sce/F");
    //my_event_->Branch("MCNu_Vz_sce",          &MCNu_Vz_sce,          "MCNu_Vz_sce/F");
    my_event_->Branch("MCNu_vertexDistance",    &MCNu_vertexDistance,    "MCNu_vertexDistance/F");
    
    my_event_->Branch("MCle_pfp",               &MCle_pfp_);
    my_event_->Branch("MCle_key",               &MCle_key_);
    my_event_->Branch("MCle_PDG",               &MCle_PDG_);
    //my_event_->Branch("MCle_purity",            &MCle_purity_);
    my_event_->Branch("MCle_Energy",            &MCle_Energy_);
    my_event_->Branch("MCle_Px",                &MCle_Px_);
    my_event_->Branch("MCle_Py",                &MCle_Py_);
    my_event_->Branch("MCle_Pz",                &MCle_Pz_);
    my_event_->Branch("MCle_Vx",                &MCle_Vx_);
    my_event_->Branch("MCle_Vy",                &MCle_Vy_);
    my_event_->Branch("MCle_Vz",                &MCle_Vz_);
    my_event_->Branch("MCle_Endx",                &MCle_Endx_);
    my_event_->Branch("MCle_Endy",                &MCle_Endy_);
    my_event_->Branch("MCle_Endz",                &MCle_Endz_);
    my_event_->Branch("MCle_length",            &MCle_length_);
    my_event_->Branch("MCle_VxSce",             &MCle_VxSce_);
    my_event_->Branch("MCle_VySce",             &MCle_VySce_);
    my_event_->Branch("MCle_VzSce",             &MCle_VzSce_);
    my_event_->Branch("MCle_Theta",             &MCle_Theta_);
    my_event_->Branch("MCle_Phi",             &MCle_Phi_);
    
    my_event_->Branch("MCTrackPurity",             &MCTrackPurity_);
    my_event_->Branch("MCTrackPDG",             &MCTrackPDG_);
    my_event_->Branch("MCTrackEnergy",             &MCTrackEnergy_);
    my_event_->Branch("MCTrackMomentum",             &MCTrackMomentum_);
    my_event_->Branch("MCTrackTheta",             &MCTrackTheta_);
    my_event_->Branch("MCTrackPhi",             &MCTrackPhi_);
    my_event_->Branch("MCTrackLength",             &MCTrackLength_);
    
    my_event_->Branch("MCTrackStart_x",             &MCTrackStart_x_);
    my_event_->Branch("MCTrackStart_y",             &MCTrackStart_y_);
    my_event_->Branch("MCTrackStart_z",             &MCTrackStart_z_);
    my_event_->Branch("MCTrackEnd_x",             &MCTrackEnd_x_);
    my_event_->Branch("MCTrackEnd_y",             &MCTrackEnd_y_);
    my_event_->Branch("MCTrackEnd_z",             &MCTrackEnd_z_);
    
    my_event_->Branch("Genie_Q2",&Genie_Q2,"Genie_Q2/D");
    my_event_->Branch("Genie_q2",&Genie_q2,"Genie_q2/D");
    my_event_->Branch("Genie_W",&Genie_W,"Genie_W/D");
    my_event_->Branch("Genie_T",&Genie_T,"Genie_T/D");
    my_event_->Branch("Genie_X",&Genie_X,"Genie_X/D");
    my_event_->Branch("Genie_Y",&Genie_Y,"Genie_Y/D");
    my_event_->Branch("Genie_nNeutron_preFSI",&Genie_nNeutron_preFSI,"Genie_nNeutron_preFSI/I");
    my_event_->Branch("Genie_nProton_preFSI",&Genie_nProton_preFSI,"Genie_nProton_preFSI/I");
    my_event_->Branch("Genie_nPi0_preFSI",&Genie_nPi0_preFSI,"Genie_nPi0_preFSI/I");
    my_event_->Branch("Genie_nPiPlus_preFSI",&Genie_nPiPlus_preFSI,"Genie_nPiPlus_preFSI/I");
    my_event_->Branch("Genie_nPiMinus_preFSI",&Genie_nPiMinus_preFSI,"Genie_nPiMinus_preFSI/I");

  }
  
  
  
}
void NumuCCana::initialize_pot()
{
  // Implementation of required member function here.
  std::cout << "Initialize variables and histograms for pot tree" << std::endl;
  //tree stuff for hits: //////////////////////////////////////////////////////////////////////////////////
  _sr_tree = tfs->make<TTree>("pottree","pottree");
  _sr_tree->Branch("run",                &_sr_run,                "run/I");
  _sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/D");
 // _sr_tree->Branch("event_count",        &event_counter,          "event count/I");
  /*std::cout << "-------Using the following fcl parameters:-------" << std::endl;
  std::cout << "Pandora label:\t\t" << m_pfp_producer << std::endl;
  std::cout << "hitfinder label:\t\t" << m_hitfinder_producer << std::endl;
  std::cout << "geant producer label:\t\t" << m_geant_producer << std::endl;
  std::cout << "mcp producer label:\t\t" << m_hit_mcp_producer << std::endl;
  std::cout << "DAQ label:\t\t" << data_label_DAQHeader_ << std::endl;
  std::cout << "Beam flash label:\t" << data_label_flash_beam_ << std::endl;
  std::cout << "mcp producer label:\t\t" << m_hit_mcp_producer << std::endl;
  std::cout << "CRT hit label:\t\t" << data_label_crthit_ << std::endl;
  std::cout << "CRT T0 asso label:\t" << data_label_crtT0asso_ << std::endl;
  
  std::cout << "fHardDelay:\t\t" << fHardDelay_ << std::endl;
  std::cout << "fCRTT0off:\t\t" << fCRTT0off_ << std::endl;
  std::cout << "Beam start:\t\t" << beam_start_ << std::endl;
  std::cout << "Beam end:\t\t" << beam_end_ << std::endl;

  std::cout << "is_data:\t\t\t" << is_data_ << std::endl;
  std::cout << "verbose:\t\t\t" << verbose_ << std::endl;
  
  std::cout << "m_isData:\t\t\t" << m_isData << std::endl;
  std::cout << "m_hasMCNeutrino:\t\t\t" << m_hasMCNeutrino << std::endl;
  
  
  std::cout << "------------end fcl parameters-------------------" << std::endl;
  */
  
}
void NumuCCana::endSubRun(art::SubRun const &sr){
  
  if( verbose_ !=0) std::cout << "Write POT infos in to tree..." << std::endl;
  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;
  if ( sr . getByLabel ( "generator" , potsum_h ) ) {
    _sr_pot = potsum_h->totpot;
    //std::cout << "Got POT info" << std::endl;
  }
  _sr_tree->Fill();
  
}

void NumuCCana::beginJob()
{
  // Implementation of optional member function here.
  initialize_tree();
  initialize_pot();
}
void NumuCCana::endJob(){
  // Implementation of optional member function here.
  std::cout.width(35); std::cout << "###### NumuCCincl summary  ######################" << std::left << std::endl;
  std::cout.width(35); std::cout << "Failed due to has one neutrino: " << std::left << counter_num_nu << std::endl;
  std::cout.width(35); std::cout << "Failed due to no muon candidate: " << std::left << counter_moun << std::endl;
  std::cout.width(35); std::cout << "Failed due to PID asso: " << std::left << counter_trackPIDAssn << std::endl;
  std::cout.width(35); std::cout << "Failed due to neutrino meta data: " << std::left << counter_neutrino_metadata_vec << std::endl;
  std::cout.width(35); std::cout << "Failed due to getPID: " << std::left << counter_getPID << std::endl;
  std::cout.width(35); std::cout << "###### End summary  #########################" << std::endl;
}

DEFINE_ART_MODULE(NumuCCana)
