//////////////////////////////////////////////////////////////////////
// Class:       CC2p
// Plugin Type: analyzer (art v3_01_02)
// File:        CC2p_module.cc
//
// Generated at Tue Feb  4 13:05:22 2020 by Samantha Sword-fehlberg using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

//Art Includes
/////////////////
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft Includes
///////////////////
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

//Uboone Includes
///////////////////
#include "ubobj/CRT/CRTHit.hh"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "ubana/CC2p/Algorithms/BackTrackerTruthMatch.h"
#include "ubana/CC2p/Algorithms/TrackFeatures.h"

#include "ubana/CC2p/Algorithms/PID.h"

#include "ubana/LEEPhotonAnalysis/ParticleAssociations.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

//Nusim Includes
////////////////
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//ROOT Includes
///////////////////
#include "TTree.h" 
#include "Objects/CartesianVector.h"

class CC2p;

class CC2p : public art::EDAnalyzer {
public:
  explicit CC2p(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  CC2p(CC2p const&) = delete;
  CC2p(CC2p&&) = delete;
  CC2p& operator=(CC2p const&) = delete;
  CC2p& operator=(CC2p&&) = delete;

  //Define all your functions
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endSubRun(art::SubRun const& sr) override;
  void FillGENIETruth(art::Event const& e);
  bool IsContained(float st_x, float st_y, float st_z, float end_x, float end_y, float end_z);
  bool IsNearVertex(float st_x, float st_y, float st_z, float end_x, float end_y, float end_z);
  float KEfromLength(float trk_length);
  void ClearLocalData();
  //void endJob() override;

private:

  //Define any specific services you might need later
  art::ServiceHandle<art::TFileService> tfs; //tfile service   

  // Declare member data here. i.e. declare the default values for all of the things you will use throughout the codes
  std::string m_pfp_producer; //pfp producer
  std::string m_hit_producer; //hit producer
  std::string m_crt_producer; //crt producer
  std::string m_geant_producer; //geant producer
  std::string m_hit_mcp_producer; //hit MCParticle producer
  std::string m_genie_producer; //GENIE producer
  std::string _hitassoclabel; //hit association label
  std::string _tpcobject_producer; //TPC Object Producer
  std::string m_potsum_producer; //POT Sum Producer
  std::string m_potsum_instance; //POT Sum Instance
  std::string fMCParticleModuleLabel; //MCParticle Module Label
  std::string fHitModuleLabel; //Hit module label
  std::string HitsPerTrackAssModuleLabel; //hits per track association label
  std::string fMCParticleToMCTruthAssModuleLabel; //MCParticle to MCTruth Association Module Label
  std::string m_flash_producer; //flassh producer
  double m_beamwindow_low; //left edge of the beam window
  double m_beamwindow_high; //right edge of the beam window
  double m_trklen_low; //
  double m_trklen_high; //
  double m_p_muon; //
  double m_p_proton; //
  double m_p_pionpm; //
  double m_p_pion0; //
  double m_p_electron; //
  double m_p_neutron; //
  double m_vtx_dis; //
  bool _debug=true; // ALWAYS DEBUG!
  bool m_doReweighting;  // Should I do the MC reweighting?
  bool m_isOverlay; //Is this overlay I am looking at?
  bool m_isData; //is this data I am looking at?
  bool m_CRT; // Should I use CRT Information?

  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles,pfparticles1,pfparticles2;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfnotneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::VertexVector pfvertices;
  lar_pandora::ShowerVector pfshowers;
  lar_pandora::TrackVector pftracks,trackVector2;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::TracksToHits tracksToHits;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToClusters particlesToClusters;
  lar_pandora::PFParticlesToShowers particlesToShowers;
  lar_pandora::PFParticlesToTracks particlesToTracks;
  lar_pandora::PFParticlesToHits pfParticlesToHits;
  lar_pandora::HitsToPFParticles hitsToPfParticles;
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;
  lar_pandora::ClustersToHits clustersToHits;
  lar_pandora::SpacePointVector spacePoints;
  lar_pandora::HitsToSpacePoints hitsToSpacePoints;
  lar_pandora::SpacePointsToHits spacePointsToHits;
  lar_pandora::PFParticlesToVertices  pfParticleToVertexMap;
  
  // Tree Information & Basic Variables
  ////////////////////////////////////////
  TTree* _tree1; //main tree
  TTree* _sr_tree; //pot tree
  int _run, _subrun, _event; //duh
  int _n_pfp_per_event=-9999; //number of pfp per event
  int _n_trk_per_event=-9999; //number of track per evet
  int _n_nu_per_event=-9999; //number of neutrinos per event
  int _n_nu_pfp_per_event=-9999; //number or nu pfp per event
  int _nu_PDG_per_event=-9999; //pdg code of the nu's
  int _n_shower_per_event=-9999; //number of showers per event
  int _n_neutrinos_per_event=-9999; //number of neutrinos per event
  float _sr_pot; //total POT /1.0E16

  //GENIE Information
  //////////////////////////
  int _mc_ccnc=-9999, _mc_mode=-9999, _mc_interactiontype=-9999,  _mc_hitnuc=-9999, _mc_nupdg=-9999; //genie variables!
  int _mc_n_muon=-9999; //number of muons in MC
  int _mc_n_proton=-9999; //number of protons in MC
  int _mc_n_electron=-9999; //number of electrons in MC
  int _mc_n_pionpm=-9999; //number of pion+- in MC
  int _mc_n_pion0=-9999; //number of pion0 in MC
  int _mc_n_neutron=-9999; //number of neutrons in MC
  int _mc_n_threshold_muon=-9999; //number of muons above threshold in MC
  int _mc_n_threshold_proton=-9999; //number of protons above threshold in MC   
  int _mc_n_threshold_electron=-9999; //number of electrons above threshold in MC   
  int _mc_n_threshold_pionpm=-9999; //number of pi+-s above threshold in MC   
  int _mc_n_threshold_pion0=-9999; //number of pi0s above threshold in MC   
  int _mc_n_threshold_neutron=-9999; //number of neutronss above threshold in MC   
  std::vector<float> _mc_g4_mom_all,_mc_g4_mom_muon,_mc_g4_mom_proton,_mc_g4_mom_electron,_mc_g4_mom_pionpm,_mc_g4_mom_pion0,_mc_g4_mom_neutron; //mc_geant4 momentum of various particles
  std::vector<float> _mc_g4_E, _mc_g4_p,_mc_g4_mass, _mc_g4_phi, _mc_g4_theta,_mc_g4_start_x,_mc_g4_start_y,_mc_g4_start_z,_mc_g4_end_x,_mc_g4_end_y,_mc_g4_end_z,_mc_g4_start_x_sce,_mc_g4_start_y_sce,_mc_g4_start_z_sce,_mc_g4_end_x_sce,_mc_g4_end_y_sce,_mc_g4_end_z_sce; //mc_geant4 energy, mass, phi, theta, start, end, start sce, end sce
  std::vector<float> _mc_length, _mc_start_x,_mc_start_y, _mc_start_z, _mc_end_x, _mc_end_y, _mc_end_z, _mc_theta, _mc_phi, _mc_ke, _mc_mom; //mc length, start, end, theta, phi, ke, mom
  std::vector<float>_mc_start_x_sce,_mc_start_y_sce, _mc_start_z_sce, _mc_end_x_sce, _mc_end_y_sce, _mc_end_z_sce; //mc start, end
  std::vector<int > _mc_pdg,  _mc_primary, _mc_origin; //mc pdg, if primary particle, origin of the particle
  std::vector<int> _mc_g4_pdg; //the mc g4 pdg
  float _mc_q2=-9999., _mc_nu_vtxx, _mc_nu_vtxy, _mc_nu_vtxz, _mc_nu_vtxx_sce, _mc_nu_vtxy_sce, _mc_nu_vtxz_sce, _mc_enu; //mc q2, vertex, sce vertex, energy of the neutrino
  float _mc_X=-9999., _mc_Y=-9999., _mc_Pt=-9999.;

  //ALL DA GENIE WEIGHTS FOR SYSTEMATICS
  float _mc_wgt=1;
  float _mc_wgt_cv=1;
  float _mc_wgt_0_EtaNCEL =1;
  float _mc_wgt_0_FrAbs_N =1;
  float _mc_wgt_0_FrAbs_pi =1;
  float _mc_wgt_0_FrCEx_N =1;
  float _mc_wgt_0_FrCEx_pi =1;
  float _mc_wgt_0_FrInel_N =1;
  float _mc_wgt_0_FrInel_pi =1;
  float _mc_wgt_0_FrPiProd_N =1;
  float _mc_wgt_0_FrPiProd_pi =1;
  float _mc_wgt_0_MFP_N =1;
  float _mc_wgt_0_MFP_pi =1;
  float _mc_wgt_0_MaCCQE =1;
  float _mc_wgt_0_MaCCRES =1;
  float _mc_wgt_0_MaNCEL =1;
  float _mc_wgt_0_MaNCRES =1;
  float _mc_wgt_0_MvCCRES =1;
  float _mc_wgt_0_MvNCRES =1;
  float _mc_wgt_0_NonRESBGvnCC1pi =1;
  float _mc_wgt_0_NonRESBGvnCC2pi =1;
  float _mc_wgt_0_NonRESBGvnNC1pi =1;
  float _mc_wgt_0_NonRESBGvnNC2pi =1;
  float _mc_wgt_0_NonRESBGvpCC1pi =1;
  float _mc_wgt_0_NonRESBGvpCC2pi =1;
  float _mc_wgt_0_NonRESBGvpNC1pi =1;
  float _mc_wgt_0_NonRESBGvpNC2pi =1;
  float _mc_wgt_0_NormCCMEC =1;
  float _mc_wgt_0_NormNCMEC =1;
  float _mc_wgt_1_EtaNCEL =1;
  float _mc_wgt_1_FrAbs_N =1;
  float _mc_wgt_1_FrAbs_pi =1;
  float _mc_wgt_1_FrCEx_N =1;
  float _mc_wgt_1_FrCEx_pi =1;
  float _mc_wgt_1_FrInel_N =1;
  float _mc_wgt_1_FrInel_pi =1;
  float _mc_wgt_1_FrPiProd_N =1;
  float _mc_wgt_1_FrPiProd_pi =1;
  float _mc_wgt_1_MFP_N =1;
  float _mc_wgt_1_MFP_pi =1;
  float _mc_wgt_1_MaCCQE =1;
  float _mc_wgt_1_MaCCRES =1;
  float _mc_wgt_1_MaNCEL =1;
  float _mc_wgt_1_MaNCRES =1;
  float _mc_wgt_1_MvCCRES =1;
  float _mc_wgt_1_MvNCRES =1;
  float _mc_wgt_1_NonRESBGvnCC1pi =1;
  float _mc_wgt_1_NonRESBGvnCC2pi =1;
  float _mc_wgt_1_NonRESBGvnNC1pi =1;
  float _mc_wgt_1_NonRESBGvnNC2pi =1;
  float _mc_wgt_1_NonRESBGvpCC1pi =1;
  float _mc_wgt_1_NonRESBGvpCC2pi =1;
  float _mc_wgt_1_NonRESBGvpNC1pi =1;
  float _mc_wgt_1_NonRESBGvpNC2pi =1;
  float _mc_wgt_1_NormCCMEC =1;
  float _mc_wgt_1_NormNCMEC =1;

  //GENIE Weight Block Borrowed From UBXSec
  /////////////////////////////////////////
  int evtwgt_genie_pm1_nfunc; ///< Number of functions used for GENIE reweighting (pm1sigma)                                                                        
  std::vector<std::string> evtwgt_genie_pm1_funcname; ///< Names of the functions used for GENIE reweighting (pm1sigma)                                             
  std::vector<int> evtwgt_genie_pm1_nweight; ///< Number of weights per function name used for GENIE reweighting (pm1sigma)                                         
  std::vector<std::vector<double>> evtwgt_genie_pm1_weight; ///< Weights per function name used for GENIE reweighting (pm1sigma)                                    
  int evtwgt_genie_multisim_nfunc; ///< Number of functions used for GENIE reweighting (multisim)                                                                   
  std::vector<std::string> evtwgt_genie_multisim_funcname; ///< Names of the functions used for GENIE reweighting (multisim)                                        
  std::vector<int> evtwgt_genie_multisim_nweight; ///< Number of weights per function name used for GENIE reweighting (multisim)                                    
  std::vector<std::vector<double>> evtwgt_genie_multisim_weight; ///< Weights per function name used for GENIE reweighting (multisim)                               
  int evtwgt_flux_multisim_nfunc; ///< Number of functions used for FLUX reweighting (multisim)                                                                     
  std::vector<std::string> evtwgt_flux_multisim_funcname; ///< Names of the functions used for FLUX reweighting (multisim)                                          
  std::vector<int> evtwgt_flux_multisim_nweight; ///< Number of weights per function name used for FLUX reweighting (multisim)                                      
  std::vector<std::vector<double>> evtwgt_flux_multisim_weight; ///< Weights per function name used for FLUX reweighting (multisim)                                 

  //Reconstructed Innformation
  ////////////////////////////////////
  float _reco_nu_vtxx=-9999., _reco_nu_vtxy=-9999.,  _reco_nu_vtxz=-9999.; //reco vertex
  std::vector<float> _reco_q2, _reco_length, _reco_start_x, _reco_start_y, _reco_start_z, _reco_end_x, _reco_end_y, _reco_end_z, _reco_theta, _reco_phi, _reco_ke, _reco_mom; //reco q2, length, start, end, theta, phi, ke, mom
  std::vector<float> _reco_mom_muon, _reco_mom_proton, _reco_mom_pion; //reco mom of muon, proton, pion
  std::vector<bool> _is_primary; //is it primary?
  std::vector<bool> _is_contained; //is it contained?
  std::vector<bool> _is_from_nu_slice; //did it come from the neutrino slice?
  std::vector<bool> _has_shower; //does it have a shower?
  std::vector<int> _n_pfp; //number of pfp particles
  std::vector<int> _n_trk; //number of tacks
  std::vector<int> _n_shower; //number of showers
  std::vector<int> _id_pfp; //id of the pfp
  std::vector<int> _n_daughters; //number of daughter particles
  std::vector<int> _parentPDG; //pdg of the parent
  std::vector<float> _track_score; //track score
  std::vector<float> _dislen_ratio; //conpare 3D calculated length to the reco length
  std::vector<float> _KE_len; //kinetic energy from length
  std::vector<float> _top_score; //neutrino score
  std::vector<float> _flash_score; //calculating the flash score
  int _vtx_n_pfp; //number of tracks near a vertex
 
  //Calorimetry and PID                                                                                                                                             
  //////////////////////////////////////////////////////// 
  std::vector<int> _nhits_0, _nhits_1, _nhits_2;
  std::vector<float> _chi2p_3D, _chi2mu_3D, _chi2pi_3D, _chi2K_3D; //the 3D chi2 from thomas 
  std::vector<float> _chi2_p_0, _chi2_p_1, _chi2_p_2;
  std::vector<float> _chi2_mu_0, _chi2_mu_1, _chi2_mu_2;
  std::vector<float> _LL3;
  std::vector<float> _LL_p_0, _LL_p_1, _LL_p_2;
  std::vector<float> _LL_mip_0, _LL_mip_1, _LL_mip_2;
  std::vector<float> _LL_mu_0, _LL_mu_1, _LL_mu_2;
  std::vector<float> _LL_back_p_0, _LL_back_p_1, _LL_back_p_2;
  std::vector<float> _LL_back_mu_0, _LL_back_mu_1, _LL_back_mu_2;
  std::vector<float> _TM_dedx_0, _TM_dedx_1, _TM_dedx_2;
  std::vector<float> _PIDA_0, _PIDA_1, _PIDA_2;
  std::vector<float> _start_dedx_0, _start_dedx_1, _start_dedx_2;
  std::vector<float> _end_dedx_0, _end_dedx_1, _end_dedx_2;
  std::vector<float> _ratio_dedx_0, _ratio_dedx_1, _ratio_dedx_2;
  std::vector<float> _avg_dedx_0, _avg_dedx_1, _avg_dedx_2;
  std::vector<float> _total_dedx_0, _total_dedx_1, _total_dedx_2;
  std::vector<float> _KE_calo_0, _KE_calo_1, _KE_calo_2;
  ::trkf::TrackMomentumCalculator _trk_mom_calculator;

};


CC2p::CC2p(fhicl::ParameterSet const& p): EDAnalyzer{p}  
{
  //Default fhicl parameters values. These can be modified by the fhicl
  /////////////////////////////////////////////////////////////////////
  m_pfp_producer      = p.get<std::string>("pfp_producer", "pandora");
  m_hit_producer      = p.get<std::string>("hit_producer", "gaushit");
  m_geant_producer    = p.get<std::string>("geant_producer", "largeant");
  m_hit_mcp_producer  = p.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
  m_flash_producer    = p.get<std::string>("flash_producer", "simpleFlashBeam");
  m_crt_producer      = p.get<std::string>("crt_producer", "crtveto");
  HitsPerTrackAssModuleLabel = p.get<std::string>("hit_trackass","pandora");
  m_beamwindow_low    = p.get<double>("Beam_low", 3.57); //beam window low
  m_beamwindow_high   = p.get<double>("Beam_high", 5.25); //beam window high
  m_trklen_low        = p.get<double>("Trklen_low", 0.5); //track length lower limit
  m_trklen_high       = p.get<double>("Trklen_high", 200); //track length upper limit
  m_vtx_dis           = p.get<double>("Vertex_distance", 5); //distance between  vertex and start/end of track
  m_p_muon            = p.get<double>("P_muon", 0.1); //momentum of muon
  m_p_proton          = p.get<double>("P_proton", 0.1); //momentum of proton
  m_p_pionpm          = p.get<double>("P_pionpm", 0.1); //momentum of  pion
  m_p_pion0           = p.get<double>("P_pion0", 0.1); //momentum of pio0
  m_p_electron        = p.get<double>("P_electron", 0.1); //momentum of electron
  m_p_neutron         = p.get<double>("P_neutron", 0.1); //momentum  of neutron
  m_doReweighting     = p.get<bool>("Reweighting", false);
  m_isOverlay         = p.get<bool>("Overlay", false);
  m_CRT               = p.get<bool>("CRTVeto", false);
  _debug              = p.get<bool>("DebugMode","false");
  m_genie_producer    = p.get<std::string>("genie_producer","generator");
  m_potsum_producer   = p.get<std::string>("POTSummaryProducer","generator");

  fMCParticleToMCTruthAssModuleLabel = p.get< std::string >("MCParticleToMCTruthAssModuleLabel","largeant");

  //Now to define my trees and their various branches
  //////////////////////////////////////////////////////
  _tree1 = tfs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("n_pfp_per_event", &_n_pfp_per_event, "n_pfp_per_event/I");
  _tree1->Branch("n_trk_per_event", &_n_trk_per_event, "n_trk_per_event/I");
  _tree1->Branch("n_shower_per_event", &_n_shower_per_event, "n_shower_per_event/I");
  _tree1->Branch("n_neutrinos_per_event",&_n_neutrinos_per_event,"n_neutrinos_per_event/I");
  _tree1->Branch("n_nu_per_event", &_n_nu_per_event, "n_nu_per_event/I");
  _tree1->Branch("n_nu_pfp_per_event", &_n_nu_pfp_per_event, "n_nu_pfp_per_event/I");
  _tree1->Branch("nu_PDG_per_event", &_nu_PDG_per_event, "nu_PDG_per_event/I");
  _tree1->Branch("mc_ccnc",&_mc_ccnc,"mc_ccnc/I");
  _tree1->Branch("mc_mode",&_mc_mode,"mc_mode/I");
  _tree1->Branch("mc_interactiontype",&_mc_interactiontype,"mc_interactiontype/I");
  _tree1->Branch("mc_hitnuc",&_mc_hitnuc,"mc_hitnuc/I");
  _tree1->Branch("mc_q2",&_mc_q2,"mc_q2/F");
  _tree1->Branch("mc_X",&_mc_X,"mc_X/F");
  _tree1->Branch("mc_Y",&_mc_Y,"mc_Y/F");
  _tree1->Branch("mc_Pt",&_mc_Pt,"mc_Pt/F");
  _tree1->Branch("mc_nu_vtxx",&_mc_nu_vtxx,"mc_nu_vtxx/F");
  _tree1->Branch("mc_nu_vtxy",&_mc_nu_vtxy,"mc_nu_vtxy/F");
  _tree1->Branch("mc_nu_vtxz",&_mc_nu_vtxz,"mc_nu_vtxz/F");
  _tree1->Branch("mc_nu_vtxx_sce",&_mc_nu_vtxx_sce,"mc_nu_vtxx_sce/F");
  _tree1->Branch("mc_nu_vtxy_sce",&_mc_nu_vtxy_sce,"mc_nu_vtxy_sce/F");
  _tree1->Branch("mc_nu_vtxz_sce",&_mc_nu_vtxz_sce,"mc_nu_vtxz_sce/F");
  _tree1->Branch("mc_enu",&_mc_enu,"mc_enu/F");
  _tree1->Branch("mc_wgt",&_mc_wgt,"mc_wgt/F");
  _tree1->Branch("mc_wgt_cv",&_mc_wgt_cv,"mc_wgt_cv/F");
  _tree1->Branch("mc_wgt_0_EtaNCEL",&_mc_wgt_0_EtaNCEL,"mc_wgt_0_EtaNCEL/F");
  _tree1->Branch("mc_wgt_0_FrAbs_N",&_mc_wgt_0_FrAbs_N,"mc_wgt_0_FrAbs_N/F");
  _tree1->Branch("mc_wgt_0_FrAbs_pi",&_mc_wgt_0_FrAbs_pi,"mc_wgt_0_FrAbs_pi/F");
  _tree1->Branch("mc_wgt_0_FrCEx_N",&_mc_wgt_0_FrCEx_N,"mc_wgt_0_FrCEx_N/F");
  _tree1->Branch("mc_wgt_0_FrCEx_pi",&_mc_wgt_0_FrCEx_pi,"mc_wgt_0_FrCEx_pi/F");
  _tree1->Branch("mc_wgt_0_FrInel_N",&_mc_wgt_0_FrInel_N,"mc_wgt_0_FrInel_N/F");
  _tree1->Branch("mc_wgt_0_FrInel_pi",&_mc_wgt_0_FrInel_pi,"mc_wgt_0_FrInel_pi/F");
  _tree1->Branch("mc_wgt_0_FrPiProd_N",&_mc_wgt_0_FrPiProd_N,"mc_wgt_0_FrPiProd_N/F");
  _tree1->Branch("mc_wgt_0_FrPiProd_pi",&_mc_wgt_0_FrPiProd_pi,"mc_wgt_0_FrPiProd_pi/F");
  _tree1->Branch("mc_wgt_0_MFP_N",&_mc_wgt_0_MFP_N,"mc_wgt_0_MFP_N/F");
  _tree1->Branch("mc_wgt_0_MFP_pi",&_mc_wgt_0_MFP_pi,"mc_wgt_0_MFP_pi/F");
  _tree1->Branch("mc_wgt_0_MaCCQE",&_mc_wgt_0_MaCCQE,"mc_wgt_0_MaCCQE/F");
  _tree1->Branch("mc_wgt_0_MaCCRES",&_mc_wgt_0_MaCCRES,"mc_wgt_0_MaCCRES/F");
  _tree1->Branch("mc_wgt_0_MaNCEL",&_mc_wgt_0_MaNCEL,"mc_wgt_0_MaNCEL/F");
  _tree1->Branch("mc_wgt_0_MaNCRES",&_mc_wgt_0_MaNCRES,"mc_wgt_0_MaNCRES/F");
  _tree1->Branch("mc_wgt_0_MvCCRES",&_mc_wgt_0_MvCCRES,"mc_wgt_0_MvCCRES/F");
  _tree1->Branch("mc_wgt_0_MvNCRES",&_mc_wgt_0_MvNCRES,"mc_wgt_0_MvNCRES/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvnCC1pi",&_mc_wgt_0_NonRESBGvnCC1pi,"mc_wgt_0_NonRESBGvnCC1pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvnCC2pi",&_mc_wgt_0_NonRESBGvnCC2pi,"mc_wgt_0_NonRESBGvnCC2pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvnNC1pi",&_mc_wgt_0_NonRESBGvnNC1pi,"mc_wgt_0_NonRESBGvnNC1pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvnNC2pi",&_mc_wgt_0_NonRESBGvnNC2pi,"mc_wgt_0_NonRESBGvnNC2pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvpCC1pi",&_mc_wgt_0_NonRESBGvpCC1pi,"mc_wgt_0_NonRESBGvpCC1pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvpCC2pi",&_mc_wgt_0_NonRESBGvpCC2pi,"mc_wgt_0_NonRESBGvpCC2pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvpNC1pi",&_mc_wgt_0_NonRESBGvpNC1pi,"mc_wgt_0_NonRESBGvpNC1pi/F");
  _tree1->Branch("mc_wgt_0_NonRESBGvpNC2pi",&_mc_wgt_0_NonRESBGvpNC2pi,"mc_wgt_0_NonRESBGvpNC2pi/F");
  _tree1->Branch("mc_wgt_0_NormCCMEC",&_mc_wgt_0_NormCCMEC,"mc_wgt_0_NormCCMEC/F");
  _tree1->Branch("mc_wgt_0_NormNCMEC",&_mc_wgt_0_NormNCMEC,"mc_wgt_0_NormNCMEC/F");
  _tree1->Branch("mc_wgt_1_EtaNCEL",&_mc_wgt_1_EtaNCEL,"mc_wgt_1_EtaNCEL/F");
  _tree1->Branch("mc_wgt_1_FrAbs_N",&_mc_wgt_1_FrAbs_N,"mc_wgt_1_FrAbs_N/F");
  _tree1->Branch("mc_wgt_1_FrAbs_pi",&_mc_wgt_1_FrAbs_pi,"mc_wgt_1_FrAbs_pi/F");
  _tree1->Branch("mc_wgt_1_FrCEx_N",&_mc_wgt_1_FrCEx_N,"mc_wgt_1_FrCEx_N/F");
  _tree1->Branch("mc_wgt_1_FrCEx_pi",&_mc_wgt_1_FrCEx_pi,"mc_wgt_1_FrCEx_pi/F");
  _tree1->Branch("mc_wgt_1_FrInel_N",&_mc_wgt_1_FrInel_N,"mc_wgt_1_FrInel_N/F");
  _tree1->Branch("mc_wgt_1_FrInel_pi",&_mc_wgt_1_FrInel_pi,"mc_wgt_1_FrInel_pi/F");
  _tree1->Branch("mc_wgt_1_FrPiProd_N",&_mc_wgt_1_FrPiProd_N,"mc_wgt_1_FrPiProd_N/F");
  _tree1->Branch("mc_wgt_1_FrPiProd_pi",&_mc_wgt_1_FrPiProd_pi,"mc_wgt_1_FrPiProd_pi/F");
  _tree1->Branch("mc_wgt_1_MFP_N",&_mc_wgt_1_MFP_N,"mc_wgt_1_MFP_N/F");
  _tree1->Branch("mc_wgt_1_MFP_pi",&_mc_wgt_1_MFP_pi,"mc_wgt_1_MFP_pi/F");
  _tree1->Branch("mc_wgt_1_MaCCQE",&_mc_wgt_1_MaCCQE,"mc_wgt_1_MaCCQE/F");
  _tree1->Branch("mc_wgt_1_MaCCRES",&_mc_wgt_1_MaCCRES,"mc_wgt_1_MaCCRES/F");
  _tree1->Branch("mc_wgt_1_MaNCEL",&_mc_wgt_1_MaNCEL,"mc_wgt_1_MaNCEL/F");
  _tree1->Branch("mc_wgt_1_MaNCRES",&_mc_wgt_1_MaNCRES,"mc_wgt_1_MaNCRES/F");
  _tree1->Branch("mc_wgt_1_MvCCRES",&_mc_wgt_1_MvCCRES,"mc_wgt_1_MvCCRES/F");
  _tree1->Branch("mc_wgt_1_MvNCRES",&_mc_wgt_1_MvNCRES,"mc_wgt_1_MvNCRES/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvnCC1pi",&_mc_wgt_1_NonRESBGvnCC1pi,"mc_wgt_1_NonRESBGvnCC1pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvnCC2pi",&_mc_wgt_1_NonRESBGvnCC2pi,"mc_wgt_1_NonRESBGvnCC2pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvnNC1pi",&_mc_wgt_1_NonRESBGvnNC1pi,"mc_wgt_1_NonRESBGvnNC1pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvnNC2pi",&_mc_wgt_1_NonRESBGvnNC2pi,"mc_wgt_1_NonRESBGvnNC2pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvpCC1pi",&_mc_wgt_1_NonRESBGvpCC1pi,"mc_wgt_1_NonRESBGvpCC1pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvpCC2pi",&_mc_wgt_1_NonRESBGvpCC2pi,"mc_wgt_1_NonRESBGvpCC2pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvpNC1pi",&_mc_wgt_1_NonRESBGvpNC1pi,"mc_wgt_1_NonRESBGvpNC1pi/F");
  _tree1->Branch("mc_wgt_1_NonRESBGvpNC2pi",&_mc_wgt_1_NonRESBGvpNC2pi,"mc_wgt_1_NonRESBGvpNC2pi/F");
  _tree1->Branch("mc_wgt_1_NormCCMEC",&_mc_wgt_1_NormCCMEC,"mc_wgt_1_NormCCMEC/F");
  _tree1->Branch("mc_wgt_1_NormNCMEC",&_mc_wgt_1_NormNCMEC,"mc_wgt_1_NormNCMEC/F");
  _tree1->Branch("evtwgt_genie_pm1_nfunc",&evtwgt_genie_pm1_nfunc,"evtwgt_genie_pm1_nfunc/I");
  _tree1->Branch("evtwgt_genie_pm1_funcname",&evtwgt_genie_pm1_funcname);
  _tree1->Branch("evtwgt_genie_pm1_nweight",&evtwgt_genie_pm1_nweight);
  _tree1->Branch("evtwgt_genie_pm1_weight",&evtwgt_genie_pm1_weight);
  _tree1->Branch("evtwgt_genie_multisim_nfunc",&evtwgt_genie_multisim_nfunc,"evtwgt_genie_multisim_nfunc/I");
  _tree1->Branch("evtwgt_genie_multisim_funcname",&evtwgt_genie_multisim_funcname);
  _tree1->Branch("evtwgt_genie_multisim_nweight",&evtwgt_genie_multisim_nweight);
  _tree1->Branch("evtwgt_genie_multisim_weight",&evtwgt_genie_multisim_weight);
  _tree1->Branch("evtwgt_flux_multisim_nfunc",&evtwgt_flux_multisim_nfunc,"evtwgt_flux_multisim_nfunc/I");
  _tree1->Branch("evtwgt_flux_multisim_funcname",&evtwgt_flux_multisim_funcname);
  _tree1->Branch("evtwgt_flux_multisim_nweight",&evtwgt_flux_multisim_nweight);
  _tree1->Branch("evtwgt_flux_multisim_weight",&evtwgt_flux_multisim_weight);
  _tree1->Branch("mc_nupdg",&_mc_nupdg,"mc_nupdg/I");
  _tree1->Branch("mc_n_muon",&_mc_n_muon,"mc_n_muon/I");
  _tree1->Branch("mc_n_proton",&_mc_n_proton,"mc_n_proton/I");
  _tree1->Branch("mc_n_pionpm",&_mc_n_pionpm,"mc_n_pionpm/I");
  _tree1->Branch("mc_n_pion0",&_mc_n_pion0,"mc_n_pion0/I");
  _tree1->Branch("mc_n_electron",&_mc_n_electron,"mc_n_electron/I");
  _tree1->Branch("mc_n_neutron",&_mc_n_neutron,"mc_n_neutron/I");
  _tree1->Branch("mc_n_threshold_muon",&_mc_n_threshold_muon,"mc_n_threshold_muon/I");
  _tree1->Branch("mc_n_threshold_proton",&_mc_n_threshold_proton,"mc_n_threshold_proton/I");
  _tree1->Branch("mc_n_threshold_pionpm",&_mc_n_threshold_pionpm,"mc_n_threshold_pionpm/I");
  _tree1->Branch("mc_n_threshold_pion0",&_mc_n_threshold_pion0,"mc_n_threshold_pion0/I");
  _tree1->Branch("mc_n_threshold_electron",&_mc_n_threshold_electron,"mc_n_threshold_electron/I");
  _tree1->Branch("mc_n_threshold_neutron",&_mc_n_threshold_neutron,"mc_n_threshold_neutron/I");
  _tree1->Branch("mc_g4_mom_all",&_mc_g4_mom_all);
  _tree1->Branch("mc_g4_mom_muon",&_mc_g4_mom_muon);
  _tree1->Branch("mc_g4_mom_proton",&_mc_g4_mom_proton);
  _tree1->Branch("mc_g4_mom_electron",&_mc_g4_mom_electron);
  _tree1->Branch("mc_g4_mom_pionpm",&_mc_g4_mom_pionpm);
  _tree1->Branch("mc_g4_mom_pion0",&_mc_g4_mom_pion0);
  _tree1->Branch("mc_g4_mom_neutron",&_mc_g4_mom_neutron);
  _tree1->Branch("mc_g4_E",&_mc_g4_E);
  _tree1->Branch("mc_g4_p",&_mc_g4_p);
  _tree1->Branch("mc_g4_mass",&_mc_g4_mass);
  _tree1->Branch("mc_g4_phi",&_mc_g4_phi);
  _tree1->Branch("mc_g4_theta",&_mc_g4_theta);
  _tree1->Branch("mc_g4_pdg",&_mc_g4_pdg);
  _tree1->Branch("mc_g4_start_x",&_mc_g4_start_x);
  _tree1->Branch("mc_g4_start_y",&_mc_g4_start_y);
  _tree1->Branch("mc_g4_start_z",&_mc_g4_start_z);
  _tree1->Branch("mc_g4_end_x",&_mc_g4_end_x);
  _tree1->Branch("mc_g4_end_y",&_mc_g4_end_y);
  _tree1->Branch("mc_g4_end_z",&_mc_g4_end_z);
  _tree1->Branch("mc_g4_start_x_sce",&_mc_g4_start_x_sce);
  _tree1->Branch("mc_g4_start_y_sce",&_mc_g4_start_y_sce);
  _tree1->Branch("mc_g4_start_z_sce",&_mc_g4_start_z_sce);
  _tree1->Branch("mc_g4_end_x_sce",&_mc_g4_end_x_sce);
  _tree1->Branch("mc_g4_end_y_sce",&_mc_g4_end_y_sce);
  _tree1->Branch("mc_g4_end_z_sce",&_mc_g4_end_z_sce);
  _tree1->Branch("is_from_nu_slice",&_is_from_nu_slice);
  _tree1->Branch("is_primary",&_is_primary);
  _tree1->Branch("is_contained",&_is_contained);
  _tree1->Branch("mc_pdg",&_mc_pdg);
  _tree1->Branch("mc_primary",&_mc_primary);
  _tree1->Branch("mc_origin",&_mc_origin);
  _tree1->Branch("mc_length",&_mc_length);
  _tree1->Branch("mc_start_x",&_mc_start_x);
  _tree1->Branch("mc_start_y",&_mc_start_y);
  _tree1->Branch("mc_start_z",&_mc_start_z);
  _tree1->Branch("mc_end_x",&_mc_end_x);
  _tree1->Branch("mc_end_y",&_mc_end_y);
  _tree1->Branch("mc_end_z",&_mc_end_z);
  _tree1->Branch("mc_start_x_sce",&_mc_start_x_sce);
  _tree1->Branch("mc_start_y_sce",&_mc_start_y_sce);
  _tree1->Branch("mc_start_z_sce",&_mc_start_z_sce);
  _tree1->Branch("mc_end_x_sce",&_mc_end_x_sce);
  _tree1->Branch("mc_end_y_sce",&_mc_end_y_sce);
  _tree1->Branch("mc_end_z_sce",&_mc_end_z_sce);
  _tree1->Branch("mc_theta",&_mc_theta);
  _tree1->Branch("mc_phi",&_mc_phi);
  _tree1->Branch("mc_ke",&_mc_ke);
  _tree1->Branch("mc_mom",&_mc_mom);
  _tree1->Branch("n_pfp",&_n_pfp);
  _tree1->Branch("n_trk",&_n_trk);
  _tree1->Branch("id_pfp",&_id_pfp);
  _tree1->Branch("n_shower",&_n_shower);
  _tree1->Branch("parentPDG",&_parentPDG);
  _tree1->Branch("trk_score",&_track_score);
  _tree1->Branch("KE_len",&_KE_len);
  _tree1->Branch("dislen_ratio",&_dislen_ratio);
  _tree1->Branch("reco_q2",&_reco_q2);
  _tree1->Branch("top_score",&_top_score);
  _tree1->Branch("n_daughters",&_n_daughters);
  _tree1->Branch("vtx_n_pfp",&_vtx_n_pfp,"vtx_n_pfp/I");
  _tree1->Branch("has_shower",&_has_shower);
  _tree1->Branch("reco_nu_vtxx",&_reco_nu_vtxx,"reco_nu_vtxx/F");
  _tree1->Branch("reco_nu_vtxy",&_reco_nu_vtxy,"reco_nu_vtxy/F");
  _tree1->Branch("reco_nu_vtxz",&_reco_nu_vtxz,"reco_nu_vtxz/F");
  _tree1->Branch("reco_length",&_reco_length);
  _tree1->Branch("reco_start_x",&_reco_start_x);
  _tree1->Branch("reco_start_y",&_reco_start_y);
  _tree1->Branch("reco_start_z",&_reco_start_z);
  _tree1->Branch("reco_end_x",&_reco_end_x);
  _tree1->Branch("reco_end_y",&_reco_end_y);
  _tree1->Branch("reco_end_z",&_reco_end_z);
  _tree1->Branch("reco_theta",&_reco_theta);
  _tree1->Branch("reco_phi",&_reco_phi);
  _tree1->Branch("reco_ke",&_reco_ke);
  _tree1->Branch("reco_mom",&_reco_mom);
  _tree1->Branch("reco_mom_muon",&_reco_mom_muon);
  _tree1->Branch("reco_mom_proton",&_reco_mom_proton);
  _tree1->Branch("reco_mom_pion",&_reco_mom_pion);
  _tree1->Branch("nhits_0",&_nhits_0);
  _tree1->Branch("nhits_1",&_nhits_1);
  _tree1->Branch("nhits_2",&_nhits_2);
  _tree1->Branch("chi2p_3D",&_chi2p_3D);
  _tree1->Branch("chi2mu_3D",&_chi2mu_3D);
  _tree1->Branch("chi2pi_3D",&_chi2pi_3D);
  _tree1->Branch("chi2K_3D",&_chi2K_3D);
  _tree1->Branch("chi2_p_0",&_chi2_p_0);
  _tree1->Branch("chi2_p_1",&_chi2_p_1);
  _tree1->Branch("chi2_p_2",&_chi2_p_2);
  _tree1->Branch("chi2_mu_0",&_chi2_mu_0);
  _tree1->Branch("chi2_mu_1",&_chi2_mu_1);
  _tree1->Branch("chi2_mu_2",&_chi2_mu_2);
  _tree1->Branch("LL3",&_LL3);
  _tree1->Branch("LL_p_0",&_LL_p_0);
  _tree1->Branch("LL_p_1",&_LL_p_1);
  _tree1->Branch("LL_p_2",&_LL_p_2);
  _tree1->Branch("LL_mip_0",&_LL_mip_0);
  _tree1->Branch("LL_mip_1",&_LL_mip_1);
  _tree1->Branch("LL_mip_2",&_LL_mip_2);
  _tree1->Branch("LL_mu_0",&_LL_mu_0);
  _tree1->Branch("LL_mu_1",&_LL_mu_1);
  _tree1->Branch("LL_mu_2",&_LL_mu_2);
  _tree1->Branch("LL_back_p_0",&_LL_back_p_0);
  _tree1->Branch("LL_back_p_1",&_LL_back_p_1);
  _tree1->Branch("LL_back_p_2",&_LL_back_p_2);
  _tree1->Branch("LL_back_mu_0",&_LL_back_mu_0);
  _tree1->Branch("LL_back_mu_1",&_LL_back_mu_1);
  _tree1->Branch("LL_back_mu_2",&_LL_back_mu_2);
  _tree1->Branch("TM_dedx_0",&_TM_dedx_0);
  _tree1->Branch("TM_dedx_1",&_TM_dedx_1);
  _tree1->Branch("TM_dedx_2",&_TM_dedx_2);
  _tree1->Branch("PIDA_0",&_PIDA_0);
  _tree1->Branch("PIDA_1",&_PIDA_1);
  _tree1->Branch("PIDA_2",&_PIDA_2);
  _tree1->Branch("start_dedx_0",&_start_dedx_0);
  _tree1->Branch("start_dedx_1",&_start_dedx_1);
  _tree1->Branch("start_dedx_2",&_start_dedx_2);
  _tree1->Branch("end_dedx_0",&_end_dedx_0);
  _tree1->Branch("end_dedx_1",&_end_dedx_1);
  _tree1->Branch("end_dedx_2",&_end_dedx_2);
  _tree1->Branch("ratio_dedx_0",&_ratio_dedx_0);
  _tree1->Branch("ratio_dedx_1",&_ratio_dedx_1);
  _tree1->Branch("ratio_dedx_2",&_ratio_dedx_2);
  _tree1->Branch("avg_dedx_0",&_avg_dedx_0);
  _tree1->Branch("avg_dedx_1",&_avg_dedx_1);
  _tree1->Branch("avg_dedx_2",&_avg_dedx_2);
  _tree1->Branch("total_dedx_0",&_total_dedx_0);
  _tree1->Branch("total_dedx_1",&_total_dedx_1);
  _tree1->Branch("total_dedx_2",&_total_dedx_2);
  _tree1->Branch("KE_calo_0",&_KE_calo_0);
  _tree1->Branch("KE_calo_1",&_KE_calo_1);
  _tree1->Branch("KE_calo_2",&_KE_calo_2);

  _sr_tree = tfs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                &_run,                "run/I");
  _sr_tree->Branch("subrun",             &_subrun,             "subrun/I");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/F");

}

void CC2p::analyze(art::Event const& e)
{

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  m_isData = e.isRealData();
  
  if(_debug) std::cout << "=======Start of 2 Proton Module for  Run/Subrun/Event: "<<_run<<" / "<<_subrun<<" / "<<_event<<"======="<<std::endl;
  
  //Truth Per Event
  //////////////////
  FillGENIETruth(e);

  //Now to start everything else
  //////////////////////////////
  larpandora.CollectPFParticleMetadata(e, m_pfp_producer, pfparticles, particlesToMetadata); //Get the particle meta data using the event and the pfp_producer
  larpandora.BuildPFParticleMap(pfparticles, particleMap); //build the pfparticle map
  larpandora.CollectTracks(e,m_pfp_producer, trackVector2, tracksToHits); //collect the tracks
  larpandora.CollectTracks(e, m_pfp_producer, pftracks, particlesToTracks);
  larpandora.CollectShowers(e, m_pfp_producer, pfshowers, particlesToShowers); //collect the showers
  larpandora.CollectVertices(e, m_pfp_producer, pfvertices, particlesToVertices); //collect all the vertices
  larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos); //Get the neutrino PFPs
  larpandora.CollectPFParticles(e, m_pfp_producer, pfparticles1, particlesToSpacePoints); //Collect your PFPs
  larpandora.CollectSpacePoints(e, m_pfp_producer, spacePoints, spacePointsToHits); //Collect the space points

  if(_debug) std::cout << "[McPfpMatch] RecoNeutrinos: " << pfneutrinos.size() << std::endl;
      
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  art::FindManyP<anab::Calorimetry> trackCaloAssn(trackHandle, e, "pandoracaliSCE");
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, "pandoracalipidSCE");
  art::Handle<std::vector<recob::PFParticle>> pfparticles_handle;
  e.getByLabel(m_pfp_producer, pfparticles_handle);
  art::FindManyP<anab::T0> nuFlashScoreAsso(pfparticles_handle, e, "flashmatch");
  art::FindManyP<recob::Track> PFPTrackAsso(pfparticles_handle, e, "pandora"); //grabbing the track assocaitons
  art::FindManyP<recob::Shower> PFPShowerAsso(pfparticles_handle, e, "pandora"); //grabbing the shower assocations

  if (!trackPIDAssn.isValid()){
    if(_debug) 
      throw cet::exception("[Numu0pi2p]")<< "trackPIDAAssn in not valid" << std::endl;
  } 
  
  if (pfparticles.size() == 0){
    if(_debug)
      throw cet::exception("[Numu0pi2p]")<< "No Reconstructed PFParticles in the event." << std::endl;
  }
  else {//if there are reconstructed pfparticles continue on your merry way

    if(_debug) std::cout << "[Numu0pi2p] There are " << pfparticles.size()<< " PFPs" << std::endl; //how many pfps in the event
    if(_debug) std::cout << "[Numu0pi2p] There are " << pftracks.size()<< " tracks" << std::endl; //better be exactly 3 for all events
    if(_debug) std::cout << "[Numu0pi2p] There are " << pfshowers.size()<< " showers" << std::endl; //better be exactly 0 for all events
    if(_debug) std::cout << "[Numu0pi2p] There are " << pfneutrinos.size()<< " neutrino candidates" << std::endl; //better be neutrinos!

    _n_pfp_per_event=pfparticles.size();
    _n_trk_per_event=pftracks.size();
    _n_shower_per_event=pfshowers.size();
    _n_neutrinos_per_event=pfneutrinos.size();
    
    if(pfneutrinos.size()==1){ //when there is exactly 1 neutrino candidate in the event
     
        if(_debug) std::cout<< "[Numu0pi2p] Starting to Fill In the Vertex" <<std::endl;

	art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
	lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
	const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();

	_reco_nu_vtxx = neutrino_vtx.X();
	_reco_nu_vtxy = neutrino_vtx.Y();
	_reco_nu_vtxz = neutrino_vtx.Z();

	if(_debug) std::cout << "[Numu0pi2p] Location of the Vertex:" << _reco_nu_vtxx << " " << _reco_nu_vtxy << " " << _reco_nu_vtxz << std::endl;
	if(_debug) std::cout << "[Numu0pi2p] The neutrino candidate's PDG code is   " << pfnu->PdgCode() << std::endl;
	
	_nu_PDG_per_event = pfnu->PdgCode();

	const std::vector<size_t> &daughterIDs = pfnu->Daughters();
	_n_nu_pfp_per_event=daughterIDs.size();

	} //end of if loop concerning events where the pfneutrinos is exactly 1

    _vtx_n_pfp=0;
    
    // Now define all variables first, so I wouldn't screw up later....
    for (unsigned int n = 0; n < pfparticles.size(); ++n)// loop over all PFParticles
      {
	float p_mc_length=-9999;
	float p_mc_start_x=-9999;
	float p_mc_start_y=-9999;
	float p_mc_start_z=-9999;
	float p_mc_end_x=-9999;
	float p_mc_end_y=-9999;
	float p_mc_end_z=-9999;
	float p_mc_start_x_sce=-9999;
	float p_mc_start_y_sce=-9999;
	float p_mc_start_z_sce=-9999;
	float p_mc_end_x_sce=-9999;
	float p_mc_end_y_sce=-9999;
	float p_mc_end_z_sce=-9999;
	float p_mc_theta=-9999;
	float p_mc_phi=-9999;
	float p_mc_ke=-9999;
	float p_mc_mom=-9999;
	int p_mc_pdg=-9999;
	int p_mc_primary=-9999;
	int p_mc_origin=-9999;
	bool p_is_primary=false;
	bool p_is_contained=false;
	bool p_is_from_nu_slice=false;
	int p_n_pfp=-9999;
	int p_n_trk=-9999;
	int p_n_shower=-9999;
	int p_parentPDG=-9999;
	float p_track_score=-9999;
	float p_dislen_ratio=-9999;
	float p_KE_len=-9999;

	float p_reco_q2=-9999;
	float p_reco_length=-9999;
	float p_reco_start_x=-9999;
	float p_reco_start_y=-9999;
	float p_reco_start_z=-9999;
	float p_reco_end_x=-9999;
	float p_reco_end_y=-9999;
	float p_reco_end_z=-9999;
	float p_reco_theta=-9999;
	float p_reco_phi=-9999;
	float p_reco_ke=-9999;
	float p_reco_mom=-9999;
	float p_reco_mom_muon=-9999;
	float p_reco_mom_proton=-9999;
	float p_reco_mom_pion=-9999;

	int p_id_pfp=-9999;
	bool p_has_shower=false;
	int p_n_daughters=-9999;
	int p_nhits_0=-9999, p_nhits_1=-9999, p_nhits_2=-9999;
	float p_chi2p_3D=-9999,p_chi2mu_3D=-9999,p_chi2pi_3D=-9999,p_chi2K_3D=-9999;
	float p_chi2_p_0=-9999, p_chi2_p_1=-9999, p_chi2_p_2=-9999;
	float p_chi2_mu_0=-9999, p_chi2_mu_1=-9999, p_chi2_mu_2=-9999;
	float p_LL3=-9999;
	float p_LL_p_0=-9999, p_LL_p_1=-9999, p_LL_p_2=-9999;
	float p_LL_mip_0=-9999, p_LL_mip_1=-9999, p_LL_mip_2=-9999;
	float p_LL_mu_0=-9999, p_LL_mu_1=-9999, p_LL_mu_2=-9999;
	float p_LL_back_p_0=-9999, p_LL_back_p_1=-9999, p_LL_back_p_2=-9999;
	float p_LL_back_mu_0=-9999, p_LL_back_mu_1=-9999, p_LL_back_mu_2=-9999;
	float p_TM_dedx_0=-9999, p_TM_dedx_1=-9999, p_TM_dedx_2=-9999;
	float p_PIDA_0=-9999, p_PIDA_1=-9999, p_PIDA_2=-9999;
	float p_start_dedx_0=-9999, p_start_dedx_1=-9999, p_start_dedx_2=-9999;
	float p_end_dedx_0=-9999, p_end_dedx_1=-9999, p_end_dedx_2=-9999;
	float p_ratio_dedx_0=-9999, p_ratio_dedx_1=-9999, p_ratio_dedx_2=-9999;
	float p_avg_dedx_0=-9999, p_avg_dedx_1=-9999, p_avg_dedx_2=-9999;
	float p_total_dedx_0=-9999, p_total_dedx_1=-9999, p_total_dedx_2=-9999;
	float p_KE_calo_0=-9999, p_KE_calo_1=-9999, p_KE_calo_2=-9999;
	p_id_pfp=n;

	const art::Ptr<recob::PFParticle> particle = pfparticles.at(n); //define an individual particle to be the nth particle of the pfparticles
	p_is_primary=particle->IsPrimary(); //is the particle the primary?
	
	if(!particle->IsPrimary()){ //if the particle is not the primary (i.e. not the neutrino)

	  p_parentPDG=larpandora.GetParentNeutrino(particleMap, particle); //grab the parent PDG
	  const auto parentIterator = particleMap.find(particle->Parent());  //find the parent of the particle from the Map
	  const int parentPDG = std::abs(parentIterator->second->PdgCode()); //what is the particles parent pdg

	  if (abs(parentPDG) == 14) p_is_from_nu_slice=true; //how we define the particle being from the neutrino slice

	  if (parentIterator == particleMap.end()) continue;
	  if (!parentIterator->second->IsPrimary()) continue;

	  const std::vector<size_t> &daughterIDs = parentIterator->second->Daughters(); //grab the PFP ids of the daughters
	  p_n_pfp=daughterIDs.size(); //number of daughters of our parent particle
	  p_n_trk=0; //initalize
	  p_n_shower=0; //initialize
	  
	  for(int j = 0; j< p_n_pfp; j++)
	    {
	      auto Iterator = particleMap.find(parentIterator->second->Daughters().at(j));
	      auto this_pfp = Iterator->second;
	      auto assoTrack = PFPTrackAsso.at(this_pfp.key());
	      auto assoShower = PFPShowerAsso.at(this_pfp.key());
	      
	      if(assoTrack.size()==1){
		p_n_trk++; //number of daughter tracks
		if(_debug) std::cout<< "This_PFP: "<<this_pfp<<" "<<this_pfp.key()<<std::endl; //it is returning the IDs of the good PFP's
		}
	      if(assoShower.size()==1){
		p_n_shower++; //number of daughter showers
		}
	    }

	  art::Ptr<recob::PFParticle> pfnu1 =parentIterator->second; //the parents of our particles
	  if(_debug) std::cout << "[Numu0pi2p]  The neutrino candidate's PDG code is   " << pfnu1->PdgCode() << std::endl;

	}// end the if loop for when the particle is NOT the primary

        
	lar_pandora::PFParticlesToTracks::const_iterator trkIter = particlesToTracks.find(particle);
       
	if (particlesToTracks.end() != trkIter){
	    const lar_pandora::TrackVector &pftracks = trkIter->second;
	    if (!pftracks.empty()) {
	      //want to ensure that each particle has only one associated track
	      if (pftracks.size() !=1){    
		  std::cout << " Warning: Found particle with more than one associated track "<<pftracks.size() << std::endl;
		    continue;
		  }

		const art::Ptr<recob::Track> track = *(pftracks.begin());
		const auto &trackVtxPosition = track->Vertex();
		const auto &trackEndPosition = track->End();
		const auto &track_direction  = track->StartDirection();

		p_reco_start_x=trackVtxPosition.x(); //start of the particle in x
		p_reco_start_y=trackVtxPosition.y(); //start of the particle in y
		p_reco_start_z=trackVtxPosition.z(); //start of the particle in z
		  
		p_reco_end_x=trackEndPosition.x(); //end of the track in x
		p_reco_end_y=trackEndPosition.y(); //end of the track in y
		p_reco_end_z=trackEndPosition.z(); //end of the track in z
		  
		p_reco_length=track->Length(); //length of the track
		p_reco_mom_muon=_trk_mom_calculator.GetTrackMomentum(p_reco_length, 13); //gets the momentum of the track using muon hypothesis
		p_reco_mom_proton=_trk_mom_calculator.GetTrackMomentum(p_reco_length, 2212); //gets the momentum of the track using the proton hypothesis
		p_reco_mom_pion=_trk_mom_calculator.GetTrackMomentum(p_reco_length, 211); //gets the momentum of the track using the pion hypothesis

		p_reco_theta=track_direction.Theta(); //theta of the track
		p_reco_phi=track_direction.Phi(); //phi of the track
		p_dislen_ratio=TMath::Sqrt( TMath::Power(p_reco_end_x-p_reco_start_x,2) + TMath::Power(p_reco_end_y-p_reco_start_y,2)  + TMath::Power(p_reco_end_z-p_reco_start_z,2))/ p_reco_length; //3d length of the track
		  
		float kelen=KEfromLength(track->Length()); //get the kinetic energy of the track from the length
		p_KE_len=kelen;
		p_is_contained=IsContained(p_reco_start_x,p_reco_start_y,p_reco_start_z, p_reco_end_x, p_reco_end_y, p_reco_end_z); //is the track contained?

		if(_debug) std::cout<<"Starting to do the vertex check"<<std::endl;

		if(!particle->IsPrimary()&&p_n_pfp>0&&(IsNearVertex(_reco_nu_vtxx,_reco_nu_vtxy,_reco_nu_vtxz,p_reco_start_x,p_reco_start_y,p_reco_start_z)||IsNearVertex(_reco_nu_vtxx,_reco_nu_vtxy,_reco_nu_vtxz,p_reco_end_x,p_reco_end_y,p_reco_end_z))){
		    _vtx_n_pfp++;
		} //the vertex check loop


		//NOW TO DO SOME PID
		/////////////////////
		if(_debug) std::cout<<"This is me trying to implement Thomas' 3D PID stuff"<<std::endl;
		PID pid;

		TVector3 track_start(trackVtxPosition.x(),trackVtxPosition.y(),trackVtxPosition.z());
		TVector3 track_end(trackEndPosition.x(),trackEndPosition.y(),trackEndPosition.z());

		pid.Chi2(trackPIDAssn,track, track_start, track_end);
		p_chi2p_3D = pid.PID_Chi2P_3pl;  // chi squqre cuts
		p_chi2mu_3D = pid.PID_Chi2Mu_3pl;
		p_chi2pi_3D = pid.PID_Chi2Pi_3pl;
		p_chi2K_3D = pid.PID_Chi2K_3pl;
		std::cout << "[1mu2p] Got track PID (muon 3pl)  : " << pid.PID_Chi2Mu_3pl << std::endl;
		std::cout << "[1mu2p] Got track PID (proton 3pl): " << pid.PID_Chi2P_3pl << std::endl;

		//NOW TO DO SOME PID: More traditional algorithms
		/////////////////////////////////////////////////
		if(_debug) std::cout<<"Now we are starting to do some PID"<<std::endl;
		std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track.key()); //grabbing all the track PIDs
		
		if (trackPID.size() != 0){
		    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
		    for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
		      {
			anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
			int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
			  
			/*  std::cout << "\n ParticleIDAlg " << AlgScore.fAlgName
			    << "\n -- Variable type: " << AlgScore.fVariableType
			    << "\n -- Track direction: " << AlgScore.fTrackDir
			    << "\n -- Assuming PDG: " << AlgScore.fAssumedPdg
			    << "\n -- Number of degrees of freedom: " << AlgScore.fNdf
			    << "\n -- Value: " << AlgScore.fValue
			    << "\n -- Using planeID: " << planenum << std::endl;
			*/



			if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
			    if (AlgScore.fAlgName == "Chi2" && TMath::Abs(AlgScore.fAssumedPdg) == 2212){
				if(planenum==0)    p_chi2_p_0=AlgScore.fValue;
				if(planenum==1)    p_chi2_p_1=AlgScore.fValue;
				if(planenum==2)    p_chi2_p_2=AlgScore.fValue;
			      }
			    if (AlgScore.fAlgName == "Chi2" && TMath::Abs(AlgScore.fAssumedPdg) == 13){
				if(planenum==0)    p_chi2_mu_0=AlgScore.fValue;
				if(planenum==1)    p_chi2_mu_1=AlgScore.fValue;
				if(planenum==2)    p_chi2_mu_2=AlgScore.fValue;
			      }
			  }
			if (AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood ){
			    if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward&&TMath::Abs(AlgScore.fAssumedPdg) == 0){
				if(planenum==0)    p_LL_mip_0=AlgScore.fValue;
				if(planenum==1)    p_LL_mip_1=AlgScore.fValue;
				if(planenum==2)    p_LL_mip_2=AlgScore.fValue;
			      }
			          
			    if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward && TMath::Abs(AlgScore.fAssumedPdg) == 2212){
				if(planenum==0)    p_LL_p_0=AlgScore.fValue;
				if(planenum==1)    p_LL_p_1=AlgScore.fValue;
				if(planenum==2)    p_LL_p_2=AlgScore.fValue;
			      }
			    if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward && TMath::Abs(AlgScore.fAssumedPdg) == 13){
				if(planenum==0)    p_LL_mu_0=AlgScore.fValue;
				if(planenum==1)    p_LL_mu_1=AlgScore.fValue;
				if(planenum==2)    p_LL_mu_2=AlgScore.fValue;
			      }
			    if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward && TMath::Abs(AlgScore.fAssumedPdg) == 2212){
				if(planenum==0)    p_LL_back_p_0=AlgScore.fValue;
				if(planenum==1)    p_LL_back_p_1=AlgScore.fValue;
				if(planenum==2)    p_LL_back_p_2=AlgScore.fValue;
			      }
			    if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward && TMath::Abs(AlgScore.fAssumedPdg) == 13){
				if(planenum==0)    p_LL_back_mu_0=AlgScore.fValue;
				if(planenum==1)    p_LL_back_mu_1=AlgScore.fValue;
				if(planenum==2)    p_LL_back_mu_2=AlgScore.fValue;
			    }
			} //BraggpeakLLh
			if (AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
			    if(planenum==0)    p_PIDA_0=AlgScore.fValue;
			    if(planenum==1)    p_PIDA_1=AlgScore.fValue;
			    if(planenum==2)    p_PIDA_2=AlgScore.fValue;
			  }
			if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward && TMath::Abs(AlgScore.fAssumedPdg) == 2212){
			    if(planenum==2)   p_LL3=AlgScore.fValue;
			}
		      }
		}//end of pid info

		//Calorimetry, dEdx, etc.
		////////////////////////////////////////////
		if(trackCaloAssn.isValid()){

		    std::vector<art::Ptr<anab::Calorimetry> > calos = trackCaloAssn.at(track.key());
		    p_nhits_0=calos.at(0)->dEdx().size();
		    p_nhits_1=calos.at(1)->dEdx().size();
		    p_nhits_2=calos.at(2)->dEdx().size();
		    p_KE_calo_0=calos.at(0)->KineticEnergy();
		    p_KE_calo_1=calos.at(1)->KineticEnergy();
		    p_KE_calo_2=calos.at(2)->KineticEnergy();
		    float f_total_dedx_0=0;
		    float f_total_dedx_1=0;
		    float f_total_dedx_2=0;
		    float f_start_dedx_0=0;
		    float f_start_dedx_1=0;
		    float f_start_dedx_2=0;
		    float f_end_dedx_0=0;
		    float f_end_dedx_1=0;
		    float f_end_dedx_2=0;
		    int ii0=-1;
		    int jj0=calos.at(0)->dEdx().size();
		    int ii1=-1;
		    int jj1=calos.at(1)->dEdx().size();
		    int ii2=-1;
		    int jj2=calos.at(2)->dEdx().size();
		    int nhits0=calos.at(0)->dEdx().size();
		    int nhits1=calos.at(1)->dEdx().size();
		    int nhits2=calos.at(2)->dEdx().size();
		    int nhit0_min=6;
		    int nhit1_min=6;
		    int nhit2_min=6;

		    if(nhits0<6){
			nhit0_min=TMath::Floor((float)nhits0/2);
		      }
		    if(nhits1<6){
			nhit1_min=TMath::Floor((float)nhits1/2);
		      }
		    if(nhits2<6){
			nhit2_min=TMath::Floor((float)nhits2/2);
		      }
		    if(calos.at(0)->dEdx().size()>0){
			for( auto& de : calos.at(0)->dEdx()){
			    f_total_dedx_0 += de;
			    ii0++;
			    jj0--;
			    if(ii0<nhit0_min){f_start_dedx_0+=de;}
			    if(jj0<nhit0_min){f_end_dedx_0+=de;}
			  }

			p_avg_dedx_0=(float)f_total_dedx_0 / (float)nhits0;
			p_start_dedx_0=f_start_dedx_0;
			p_end_dedx_0=f_end_dedx_0;
			p_total_dedx_0=f_total_dedx_0;

			if(f_start_dedx_0>0){
			    p_ratio_dedx_0=f_end_dedx_0/f_start_dedx_0;
			  }
		      }//plane 0
		    if(calos.at(1)->dEdx().size()>0){
			for( auto& de : calos.at(1)->dEdx()){
			    f_total_dedx_1 += de;
			    ii1++;
			    jj1--;
			    if(ii1<nhit1_min){f_start_dedx_1+=de;}
			    if(jj1<nhit1_min){f_end_dedx_1+=de;}
			  }
			
			p_avg_dedx_1= (float)f_total_dedx_1 / (float)nhits1;
			p_start_dedx_1=f_start_dedx_1;
			p_end_dedx_1=f_end_dedx_1;
			p_total_dedx_1=f_total_dedx_1;
			
			if(f_start_dedx_1>0){
			    p_ratio_dedx_1=f_end_dedx_1/f_start_dedx_1;
			  }
		      }//plane 1
		    
		    if(calos.at(2)->dEdx().size()>0){
			for( auto& de : calos.at(2)->dEdx()){
			    f_total_dedx_2 += de;
			    ii2++;
			    jj2--;
			    if(ii2<nhit2_min){f_start_dedx_2+=de;}
			    if(jj2<nhit2_min){f_end_dedx_2+=de;}
			  }
			
			p_avg_dedx_2= (float)f_total_dedx_2 / (float)nhits2;
			p_start_dedx_2=f_start_dedx_2;
			p_end_dedx_2=f_end_dedx_2;
			p_total_dedx_2=f_total_dedx_2;

			if(f_start_dedx_2>0){
			    p_ratio_dedx_2=f_end_dedx_2/f_start_dedx_2;
			  }
		      }//plane 2
		  }//end of calo-related stuffs
		
		
		//Fill in the truth information for the track if you are looking at Overlay
		////////////////////////////////////////////////////////////////////////////                 
		if(m_isOverlay){
		    
		    lar_pandora::TracksToHits::const_iterator trkIter2 = tracksToHits.find(track); //map of track and hit objects
		    std::vector<art::Ptr<recob::Hit>> hitVec; //vector of recob hits
		    hitVec.clear(); //clear the vector
		    hitVec= trkIter2->second; //hits for a particular track
		    
		    if (_debug)   std::cout << "[Numu0pi2p] DEBUGGING 0 hitVec size "<<hitVec.size()<< std::endl;
		    art::Handle<std::vector<simb::MCParticle>> mcparticle_handle; //handle of the vector of simb::MCparticles
		    std::vector<art::Ptr<simb::MCParticle> > largeant_vec; //actual vector pointer for MCParticle
		    
		    if (e.getByLabel(m_geant_producer,mcparticle_handle)) { art::fill_ptr_vector(largeant_vec,mcparticle_handle); }
		    art::Handle<std::vector<recob::Hit > > hit_handle;
		    std::vector<art::Ptr<recob::Hit> > hit_vec;
		    hit_vec.clear();
		    
		    if(e.getByLabel(m_hit_producer,hit_handle)) { art::fill_ptr_vector(hit_vec,hit_handle); }
		    if (_debug)   std::cout << "[Numu0pi2p] DEBUGGING 0 hit_vec size "<<hit_vec.size()<< std::endl;

		    //art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,m_geant_producer); //I wonder if this is correct
		    //These lines come directly from the backtracker documentation in the overlay stuff
		    //art:;Ptr<recob::Track> CurrentTrack = trk_vec.at(i_t);
		    art::FindManyP<recob::Hit> hits_per_track(trackHandle,e,HitsPerTrackAssModuleLabel);
		    std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(track.key());
		    
		    //this matches exactly from Afro's CCQE code
		    art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,fMCParticleToMCTruthAssModuleLabel);
		    BackTrackerTruthMatch backtrackertruthmatch;
		    backtrackertruthmatch.MatchToMCParticle(hit_handle,e,trk_hits_ptrs);
		    art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();
		    
		    //What lu was using before
		    //BackTrackerTruthMatch backtrackertruthmatch;
		    //backtrackertruthmatch.MatchToMCParticle(hit_handle,e,hitVec);
		    // art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();

		    if (_debug) std::cout<<"Purity for a given track: "<<backtrackertruthmatch.ReturnPurity()<<std::endl;
		    if (_debug) std::cout<<"Completeness for a given track: "<<backtrackertruthmatch.ReturnCompleteness()<<std::endl;

		    if(!(maxp_me.isNull()))
		      {
			if (_debug) std::cout << "Track matched by the BackTracker" << std::endl;
			const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruth.at(maxp_me.key());
			p_mc_pdg=maxp_me->PdgCode();
			p_mc_primary=float(maxp_me->Process() == "primary");
			p_mc_origin=mctruth->Origin();
			p_mc_length=maxp_me->Trajectory().TotalLength();
			p_mc_start_x= maxp_me->Vx();
			p_mc_start_y=maxp_me->Vy();
			p_mc_start_z=maxp_me->Vz();
			p_mc_end_x=maxp_me->EndX();
			p_mc_end_y=maxp_me->EndY();
			p_mc_end_z=maxp_me->EndZ();
			p_mc_theta=maxp_me->Momentum().Theta();
			p_mc_phi=maxp_me->Momentum().Phi();
			p_mc_ke=maxp_me->E() - maxp_me->Mass();
			p_mc_mom=maxp_me->P();

			if(_debug) std::cout<<"MCTruth Origin: "<<mctruth->Origin()<<std::endl;
			if(_debug) std::cout<<"Particle's MC PDG Code: "<<p_mc_pdg<<std::endl;

			//Make sure to do space charge correction
			///////////////////////////////////////                                          
			auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); //get the space charge correction                
			auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); //Get time offset for X space charge correction
			auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
			auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
			auto gen = mct_h->at(0);
			double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
			double xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(p_mc_start_x,p_mc_start_y,p_mc_start_z));
			auto sce_offset_end = SCE->GetPosOffsets(geo::Point_t(p_mc_end_x,p_mc_end_y,p_mc_end_z));

			p_mc_start_x_sce = (p_mc_start_x - sce_offset.X() + xtimeoffset) + 0.6;
			p_mc_start_y_sce = p_mc_start_y + sce_offset.Y();
			p_mc_start_z_sce = p_mc_start_z + sce_offset.Z();
			p_mc_end_x_sce = (p_mc_end_x - sce_offset_end.X() + xtimeoffset) + 0.6;
			p_mc_end_y_sce = p_mc_end_y + sce_offset_end.Y();
			p_mc_end_z_sce = p_mc_end_z + sce_offset_end.Z();

		      } else {
			if (_debug) std::cout << "Track not matched by the BackTracker- Cosmic!" << std::endl;
		      }
		    hitVec.clear();
		    largeant_vec.clear();
		} //End of the Overlay loop
		
		if(_debug)
		  {
		    std::cout<<"[Numu0pi2p] track score "<<p_track_score<<std::endl;
		    std::cout<<"[Numu0pi2p] reco vtx x, y, z: ("<<p_reco_start_x<<", "<<p_reco_start_y<<", "<<p_reco_start_z<<" )"<<std::endl;
		    std::cout<<"[Numu0pi2p] reco end x, y, z: ("<<p_reco_end_x<<", "<<p_reco_end_y<<", "<<p_reco_end_z<<" )"<<std::endl;
		    std::cout<<"[Numu0pi2p] reco length: "<< p_reco_length <<std::endl;
		    std::cout<<"[Numu0pi2p] reco theta: "<< p_reco_theta <<std::endl;
		    std::cout<<"[Numu0pi2p] reco phi: "<< p_reco_phi <<std::endl;
		    std::cout<<"[Numu0pi2p] chi2 proton 0, 1, 2: "<<p_chi2_p_0<<", "<<p_chi2_p_1<<", "<<p_chi2_p_2<<std::endl;
		    std::cout<<"[Numu0pi2p] chi2 muon 0, 1, 2: "<<p_chi2_mu_0<<", "<<p_chi2_mu_1<<", "<<p_chi2_mu_2<<std::endl;
		    std::cout<<"[Numu0pi2p] chi2 3D: proton, muon, pion, kaon "<<p_chi2p_3D<<", "<<p_chi2mu_3D<<", "<<p_chi2pi_3D<<", "<<p_chi2K_3D<<std::endl;
		    std::cout<<"[Numu0pi2p] 3-plane LL: "<<p_LL3<<std::endl;
		    std::cout<<"[Numu0pi2p] KE 0 "<<p_KE_calo_0<<std::endl;
		    std::cout<<"[Numu0pi2p] KE 1 "<<p_KE_calo_1<<std::endl;
		    std::cout<<"[Numu0pi2p] KE 2 "<<p_KE_calo_2<<std::endl;
		    std::cout<<"[Numu0pi2p] NHits 0 "<<p_nhits_0<<std::endl;
		    std::cout<<"[Numu0pi2p] NHits 1 "<<p_nhits_1<<std::endl;
		    std::cout<<"[Numu0pi2p] NHits 2 "<<p_nhits_2<<std::endl;
		    std::cout<<"[Numu0pi2p] true PDG: "<< p_mc_pdg<<std::endl;
		    std::cout<<"[Numu0pi2p] true origin: "<< p_mc_origin<<std::endl;
		    std::cout<<"[Numu0pi2p] true K.E.: "<< p_mc_ke<<std::endl;
		    std::cout<<"[Numu0pi2p] true start x, y, z: ("<<p_mc_start_x<<", "<<p_mc_start_y<<", "<<p_mc_start_z<<" )"<<std::endl;
		    std::cout<<"[Numu0pi2p] true end x, y, z: ("<< p_mc_end_x<<", "<< p_mc_end_y<<", "<< p_mc_end_z<<" )"<<std::endl;
		    std::cout<<"[Numu0pi2p] true length: "<< p_mc_length<<std::endl;
		    std::cout<<"[Numu0pi2p] true theta: "<< p_mc_theta<<std::endl;
		    std::cout<<"[Numu0pi2p] true phi: "<< p_mc_phi<<std::endl;
		    std::cout<<"[Numu0pi2p] --------------end of PFP  "<<n<<" ----------------"<<std::endl;
		  }//if debug 

	    }//if trackvector is not empty                                                                                                                            
	}//if there is a track 
	
	_mc_length.push_back(p_mc_length);
	_mc_start_x.push_back(p_mc_start_x);
	_mc_start_y.push_back(p_mc_start_y);
	_mc_start_z.push_back(p_mc_start_z);
	_mc_end_x.push_back(p_mc_end_x);
	_mc_end_y.push_back(p_mc_end_y);
	_mc_end_z.push_back(p_mc_end_z);
	_mc_start_x_sce.push_back(p_mc_start_x_sce);
	_mc_start_y_sce.push_back(p_mc_start_y_sce);
	_mc_start_z_sce.push_back(p_mc_start_z_sce);
	_mc_end_x_sce.push_back(p_mc_end_x_sce);
	_mc_end_y_sce.push_back(p_mc_end_y_sce);
	_mc_end_z_sce.push_back(p_mc_end_z_sce);
	_mc_theta.push_back(p_mc_theta);
	_mc_phi.push_back(p_mc_phi);
	_mc_ke.push_back(p_mc_ke);
	_mc_mom.push_back(p_mc_mom);
	_mc_pdg.push_back(p_mc_pdg);
	_mc_primary.push_back(p_mc_primary);
	_mc_origin.push_back(p_mc_origin);
	_is_contained.push_back(p_is_contained);
	_is_primary.push_back(p_is_primary);
	_is_from_nu_slice.push_back(p_is_from_nu_slice);
	_n_pfp.push_back(p_n_pfp);
	_id_pfp.push_back(p_id_pfp);
	_n_trk.push_back(p_n_trk);
	_n_shower.push_back(p_n_shower);
	_parentPDG.push_back(p_parentPDG);
	_track_score.push_back(p_track_score);
	_dislen_ratio.push_back(p_dislen_ratio);
	_KE_len.push_back(p_KE_len);
	_reco_q2.push_back(p_reco_q2);
	_reco_length.push_back(p_reco_length);
	_reco_start_x.push_back(p_reco_start_x);
	_reco_start_y.push_back(p_reco_start_y);
	_reco_start_z.push_back(p_reco_start_z);
	_reco_end_x.push_back(p_reco_end_x);
	_reco_end_y.push_back(p_reco_end_y);
	_reco_end_z.push_back(p_reco_end_z);
	_reco_theta.push_back(p_reco_theta);
	_reco_phi.push_back(p_reco_phi);
	_reco_ke.push_back(p_reco_ke);
	_reco_mom.push_back(p_reco_mom);
	_reco_mom_muon.push_back(p_reco_mom_muon);
	_reco_mom_proton.push_back(p_reco_mom_proton);
	_reco_mom_pion.push_back(p_reco_mom_pion);
	_nhits_0.push_back(p_nhits_0);
	_nhits_1.push_back(p_nhits_1);
	_nhits_2.push_back( p_nhits_2);
	_chi2_p_0.push_back(p_chi2_p_0);
	_chi2_p_1.push_back( p_chi2_p_1);
	_chi2_p_2.push_back( p_chi2_p_2);
	_chi2_mu_0.push_back(p_chi2_mu_0);
	_chi2_mu_1.push_back( p_chi2_mu_1);
	_chi2_mu_2.push_back( p_chi2_mu_2);
	_chi2p_3D.push_back(p_chi2p_3D);
	_chi2mu_3D.push_back(p_chi2mu_3D);
	_chi2pi_3D.push_back(p_chi2pi_3D);
	_chi2K_3D.push_back(p_chi2K_3D);
	_has_shower.push_back(p_has_shower);
	_n_daughters.push_back(p_n_daughters);
	_LL3.push_back(p_LL3);
	_LL_p_0.push_back(p_LL_p_0);
	_LL_p_1.push_back(p_LL_p_1);
	_LL_p_2.push_back(p_LL_p_2);
	_LL_mip_0.push_back(p_LL_mip_0);
	_LL_mip_1.push_back(p_LL_mip_1);
	_LL_mip_2.push_back(p_LL_mip_2);
	_LL_mu_0.push_back(p_LL_mu_0);
	_LL_mu_1.push_back(p_LL_mu_1);
	_LL_mu_2.push_back(p_LL_mu_2);
	_LL_back_p_0.push_back(p_LL_back_p_0);
	_LL_back_p_1.push_back( p_LL_back_p_1);
	_LL_back_p_2.push_back( p_LL_back_p_2);
	_LL_back_mu_0.push_back(p_LL_back_mu_0);
	_LL_back_mu_1.push_back( p_LL_back_mu_1);
	_LL_back_mu_2.push_back( p_LL_back_mu_2);
	_TM_dedx_0.push_back(p_TM_dedx_0);
	_TM_dedx_1.push_back( p_TM_dedx_1);
	_TM_dedx_2.push_back( p_TM_dedx_2);
	_PIDA_0.push_back(p_PIDA_0);
	_PIDA_1.push_back( p_PIDA_1);
	_PIDA_2.push_back( p_PIDA_2);
	_start_dedx_0.push_back(p_start_dedx_0);
	_start_dedx_1.push_back( p_start_dedx_1);
	_start_dedx_2.push_back( p_start_dedx_2);
	_end_dedx_0.push_back(p_end_dedx_0);
	_end_dedx_1.push_back( p_end_dedx_1);
	_end_dedx_2.push_back( p_end_dedx_2);
	_ratio_dedx_0.push_back(p_ratio_dedx_0);
	_ratio_dedx_1.push_back(p_ratio_dedx_1);
	_ratio_dedx_2.push_back( p_ratio_dedx_2);
	_avg_dedx_0.push_back(p_avg_dedx_0);
	_avg_dedx_1.push_back(p_avg_dedx_1);
	_avg_dedx_2.push_back( p_avg_dedx_2);
	_total_dedx_0.push_back(p_total_dedx_0);
	_total_dedx_1.push_back( p_total_dedx_1);
	_total_dedx_2.push_back( p_total_dedx_2);
	_KE_calo_0.push_back(p_KE_calo_0);
	_KE_calo_1.push_back( p_KE_calo_1);
	_KE_calo_2.push_back(p_KE_calo_2);
	                                                                              
      }//end of looping over all PFParticles                                                                 
  }// end of everything that has more than 0 PFParticles                                                   
_tree1->Fill();
ClearLocalData();

}//end of the CC2p analyze

//Here is where the rest of the functions are defined
/////////////////////////////////////////////////////
void CC2p::endSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.                                                                                                  
  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  _run       = sr.run();
  _subrun    = sr.subRun();

  art::Handle<sumdata::POTSummary> potsum_h;

  // MC                                                                                                                                                
  if (m_isOverlay) {
    if (_debug) std::cout << "[NCE::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(m_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[NCE::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot/1.0e16;
      //std::cout <<"MC pot:  "<<_sr_pot << std::endl;                                                                                             
    } 
    else
      _sr_pot = 0.;
    if (_debug) std::cout <<"MC pot:  "<<_sr_pot << std::endl;
  }
  _sr_tree->Fill();
}

void CC2p::ClearLocalData()
{
  pfparticles.clear();
  pfneutrinos.clear();
  pfnotneutrinos.clear();
  pfdaughters.clear();
  pfshowers.clear();
  pftracks.clear();
  particleMap.clear();
  particlesToMetadata.clear();
  particlesToVertices.clear();
  particlesToClusters.clear();
  particlesToSpacePoints.clear();
  particlesToShowers.clear();
  particlesToTracks.clear();
  clustersToHits.clear();
  hitsToSpacePoints.clear();
  spacePointsToHits.clear();
  tracksToHits.clear();
  _mc_ccnc=-9999; 
  _mc_mode=-9999; 
  _mc_nupdg=-9999; 
  _mc_wgt=1; 
  _mc_wgt_cv=1;
  _mc_interactiontype=-9999; 
  _mc_hitnuc=-9999; 
  _mc_nupdg=-9999;
  _mc_q2=-9999.;
  _mc_X=-9999.;
  _mc_Y=-9999.;
  _mc_Pt=-9999.;
  _nhits_0.clear(); 
  _nhits_1.clear(); 
  _nhits_2.clear();
  _mc_length.clear();
  _mc_start_x.clear();
  _mc_start_y.clear();
  _mc_start_z.clear();
  _mc_end_x.clear();
  _mc_end_y.clear(); 
  _mc_end_z.clear();
  _mc_start_x_sce.clear();
  _mc_start_y_sce.clear();
  _mc_start_z_sce.clear();
  _mc_end_x_sce.clear();
  _mc_end_y_sce.clear(); 
  _mc_end_z_sce.clear();
  _mc_theta.clear();
  _mc_phi.clear();
  _mc_ke.clear();
  _mc_mom.clear();
  _mc_pdg.clear();
  _mc_primary.clear(); 
  _mc_origin.clear(); 
  _is_contained.clear();
  _n_pfp.clear();
  _id_pfp.clear();
  _n_trk.clear();
  _n_shower.clear();
  _mc_g4_start_x.clear();
  _mc_g4_start_y.clear();
  _mc_g4_start_z.clear();
  _mc_g4_end_x.clear();
  _mc_g4_end_y.clear(); 
  _mc_g4_end_z.clear();
  _mc_g4_start_x_sce.clear();
  _mc_g4_start_y_sce.clear();
  _mc_g4_start_z_sce.clear();
  _mc_g4_end_x_sce.clear();
  _mc_g4_end_y_sce.clear(); 
  _mc_g4_end_z_sce.clear();
  _mc_g4_mom_all.clear();
  _mc_g4_mom_muon.clear();
  _mc_g4_mom_proton.clear();
  _mc_g4_mom_electron.clear();
  _mc_g4_mom_pionpm.clear();
  _mc_g4_mom_pion0.clear();
  _mc_g4_mom_neutron.clear();
  _mc_g4_E.clear();
  _mc_g4_p.clear();
  _mc_g4_mass.clear();
  _mc_g4_phi.clear();
  _mc_g4_theta.clear();
  _mc_g4_pdg.clear();
  _parentPDG.clear();
  _is_from_nu_slice.clear();
  _is_primary.clear();
  _dislen_ratio.clear();
  _KE_len.clear();
  _has_shower.clear();
  _n_daughters.clear();
  _vtx_n_pfp=-9999;
  _reco_q2.clear();
  _reco_nu_vtxx=-9999.;
  _reco_nu_vtxy=-9999.;
  _reco_nu_vtxz=-9999.;
  _reco_length.clear(); 
  _reco_start_x.clear(); 
  _reco_start_y.clear(); 
  _reco_start_z.clear(); 
  _reco_end_x.clear(); 
  _reco_end_y.clear(); 
  _reco_end_z.clear(); 
  _reco_theta.clear(); 
  _reco_phi.clear();
  _reco_ke.clear(); 
  _reco_mom.clear();
  _reco_mom_muon.clear();
  _reco_mom_proton.clear();
  _reco_mom_pion.clear();
  _chi2_p_0.clear(); 
  _chi2_p_1.clear(); 
  _chi2_p_2.clear();
  _chi2_mu_0.clear(); 
  _chi2_mu_1.clear(); 
  _chi2_mu_2.clear();
  _chi2p_3D.clear();
  _chi2mu_3D.clear();
  _chi2pi_3D.clear();
  _chi2K_3D.clear();
  _LL3.clear();
  _LL_p_0.clear();
  _LL_p_1.clear();
  _LL_p_2.clear();
  _LL_mip_0.clear();
  _LL_mip_1.clear();
  _LL_mip_2.clear();
  _LL_mu_0.clear();
  _LL_mu_1.clear();
  _LL_mu_2.clear();
  _LL_back_p_0.clear();
  _LL_back_p_1.clear();
  _LL_back_p_2.clear();
  _LL_back_mu_0.clear();
  _LL_back_mu_1.clear();
  _LL_back_mu_2.clear();
  _TM_dedx_0.clear();
  _TM_dedx_1.clear();
  _TM_dedx_2.clear();
  _PIDA_0.clear();
  _PIDA_1.clear();
  _PIDA_2.clear();
  _start_dedx_0.clear();
  _start_dedx_1.clear();
  _start_dedx_2.clear();
  _end_dedx_0.clear();
  _end_dedx_1.clear();
  _end_dedx_2.clear();
  _ratio_dedx_0.clear();
  _ratio_dedx_1.clear();
  _ratio_dedx_2.clear();
   _avg_dedx_0.clear();
  _avg_dedx_1.clear();
  _avg_dedx_2.clear();
   _KE_calo_0.clear();
  _KE_calo_1.clear();
  _KE_calo_2.clear();
  _total_dedx_0.clear();
  _total_dedx_1.clear();
  _total_dedx_2.clear();
  _start_dedx_0.clear();
  _start_dedx_1.clear();
  _start_dedx_2.clear();
  evtwgt_genie_pm1_funcname.clear();
  evtwgt_genie_pm1_weight.clear();
  evtwgt_genie_pm1_nweight.clear();
  evtwgt_genie_multisim_funcname.clear();
  evtwgt_genie_multisim_weight.clear();
  evtwgt_genie_multisim_nweight.clear();
  evtwgt_flux_multisim_funcname.clear();
  evtwgt_flux_multisim_weight.clear();
  evtwgt_flux_multisim_nweight.clear();
}

//Fill in all of the GENIE Truth
/////////////////////////////////
void CC2p::FillGENIETruth(art::Event const& e)
{
  if(m_isOverlay)
    {
      if(m_doReweighting)
	{
	  //Flux
	  art::InputTag flux_eventweight_tag("eventweightFluxJan29");
	  art::Handle<std::vector<evwgh::MCEventWeight>> fluxeventweight_h;
	  e.getByLabel(flux_eventweight_tag, fluxeventweight_h);
	  if(!fluxeventweight_h.isValid())
	    {
	      std::cout << "[NumuCC0pi2p] MCEventWeight for FLUX reweight multisim, product not found..." << std::endl;
	    } 
	  else 
	    {
	      std::vector<art::Ptr<evwgh::MCEventWeight>> fluxeventweight_v;
	      art::fill_ptr_vector(fluxeventweight_v, fluxeventweight_h);
	      if (fluxeventweight_v.size() > 0)
		{
		  art::Ptr<evwgh::MCEventWeight> evt_wgt = fluxeventweight_v.at(0); // Just for the first nu interaction
		  std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
		  int countFunc = 0;
		  // loop over the map and save the name of the function and the vector of weights for each function
		  for(auto it : evtwgt_map) {
		    std::string func_name = it.first;
		    std::vector<double> weight_v = it.second;
		    evtwgt_flux_multisim_funcname.push_back(func_name);
		    evtwgt_flux_multisim_weight.push_back(weight_v);
		    evtwgt_flux_multisim_nweight.push_back(weight_v.size());
		    countFunc++;
		    if(_debug)
		      {
			for ( size_t u = 0; u < weight_v.size(); ++u ) 
			  {
			    std::cout<<"name  "<<func_name<<" universe "<<u<<" :  "<<weight_v.at(u)<<std::endl;
			  }
		      }
		  }
		  evtwgt_flux_multisim_nfunc = countFunc;
		}
	    }
	  
	  //GENIE multisim
	  art::InputTag genie_eventweight_tag("eventweightGenieALL");
	  art::Handle<std::vector<evwgh::MCEventWeight>> genieeventweight_h;
	  e.getByLabel(genie_eventweight_tag, genieeventweight_h);
	  if(!genieeventweight_h.isValid())
	    {
	      std::cout << "[NumuCC0pi2p] MCEventWeight for GENIE reweight multisim, product not found..." << std::endl;
	    } 
	  else 
	    {
	      std::vector<art::Ptr<evwgh::MCEventWeight>> genieeventweight_v;
	      art::fill_ptr_vector(genieeventweight_v, genieeventweight_h);
	      if (genieeventweight_v.size() > 0)
		{
		  art::Ptr<evwgh::MCEventWeight> evt_wgt = genieeventweight_v.at(0); // Just for the first nu interaction
		  std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
		  int countFunc = 0;
		  // loop over the map and save the name of the function and the vector of weights for each function
		  for(auto it : evtwgt_map) 
		    {
		      std::string func_name = it.first;
		      std::vector<double> weight_v = it.second;
		      evtwgt_genie_multisim_funcname.push_back(func_name);
		      evtwgt_genie_multisim_weight.push_back(weight_v);
		      evtwgt_genie_multisim_nweight.push_back(weight_v.size());
		      countFunc++;
		      if(_debug)
			{
			  for ( size_t u = 0; u < weight_v.size(); ++u ) 
			    {
			      std::cout<<"name  "<<func_name<<" universe "<<u<<" :  "<<weight_v.at(u)<<std::endl;
			    }
			}
		    }
		  evtwgt_genie_multisim_nfunc = countFunc;
		}
	    }
	  
	  
	  //GENIE PM and cv tune and v4fix wgt
	  art::InputTag eventweight_tag("eventweightGeniePM");
	  art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
	  if (e.getByLabel(eventweight_tag, eventweights_handle))
	    {
	      std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
	      art::fill_ptr_vector(eventweights, eventweights_handle);
	      std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
	      int countFunc = 0;
	      for(auto it : evtwgt_map) 
		{
		  std::string func_name = it.first;
		  std::vector<double> weight_v = it.second; 
		  evtwgt_genie_pm1_funcname.push_back(func_name);
		  evtwgt_genie_pm1_weight.push_back(weight_v);
		  evtwgt_genie_pm1_nweight.push_back(weight_v.size());
		  countFunc++;

		  if(_debug)
		    {
		      for ( size_t u = 0; u < weight_v.size(); ++u ) 
			{
			  std::cout<<"name  "<<func_name<<" universe "<<u<<" :  "<<weight_v.at(u)<<std::endl;
			}
		    }
		}
	      evtwgt_genie_pm1_nfunc = countFunc;
	      
	      const std::vector<double> &weights_cv = evtwgt_map.at("TunedCentralValue_Genie");
	      const std::vector<double> &weights = evtwgt_map.at("splines_general_Spline");
	      const std::vector<double> &weights_EtaNCEL = evtwgt_map.at("EtaNCEL_Genie");
	      const std::vector<double> &weights_FrAbs_N  = evtwgt_map.at("FrAbs_N_Genie");
	      const std::vector<double> &weights_FrAbs_pi  = evtwgt_map.at("FrAbs_pi_Genie");
	      const std::vector<double> &weights_FrCEx_N  = evtwgt_map.at("FrCEx_N_Genie");
	      const std::vector<double> &weights_FrCEx_pi  = evtwgt_map.at("FrCEx_pi_Genie");
	      const std::vector<double> &weights_FrInel_N  = evtwgt_map.at("FrInel_N_Genie");
	      const std::vector<double> &weights_FrInel_pi  = evtwgt_map.at("FrInel_pi_Genie");
	      const std::vector<double> &weights_FrPiProd_N  = evtwgt_map.at("FrPiProd_N_Genie");
	      const std::vector<double> &weights_FrPiProd_pi  = evtwgt_map.at("FrPiProd_pi_Genie");
	      const std::vector<double> &weights_MFP_N  = evtwgt_map.at("MFP_N_Genie");
	      const std::vector<double> &weights_MFP_pi  = evtwgt_map.at("MFP_pi_Genie");
	      const std::vector<double> &weights_MaCCQE  = evtwgt_map.at("MaCCQE_Genie");
	      const std::vector<double> &weights_MaCCRES  = evtwgt_map.at("MaCCRES_Genie");
	      const std::vector<double> &weights_MaNCEL  = evtwgt_map.at("MaNCEL_Genie");
	      const std::vector<double> &weights_MaNCRES  = evtwgt_map.at("MaNCRES_Genie");
	      const std::vector<double> &weights_MvCCRES  = evtwgt_map.at("MvCCRES_Genie");
	      const std::vector<double> &weights_MvNCRES  = evtwgt_map.at("MvNCRES_Genie");
	      const std::vector<double> &weights_NonRESBGvnCC1pi  = evtwgt_map.at("NonRESBGvnCC1pi_Genie");
	      const std::vector<double> &weights_NonRESBGvnCC2pi  = evtwgt_map.at("NonRESBGvnCC2pi_Genie");
	      const std::vector<double> &weights_NonRESBGvnNC1pi  = evtwgt_map.at("NonRESBGvnNC1pi_Genie");
	      const std::vector<double> &weights_NonRESBGvnNC2pi  = evtwgt_map.at("NonRESBGvnNC2pi_Genie");
	      const std::vector<double> &weights_NonRESBGvpCC1pi  = evtwgt_map.at("NonRESBGvpCC1pi_Genie");
	      const std::vector<double> &weights_NonRESBGvpCC2pi  = evtwgt_map.at("NonRESBGvpCC2pi_Genie");
	      const std::vector<double> &weights_NonRESBGvpNC1pi  = evtwgt_map.at("NonRESBGvpNC1pi_Genie");
	      const std::vector<double> &weights_NonRESBGvpNC2pi  = evtwgt_map.at("NonRESBGvpNC2pi_Genie");
	      const std::vector<double> &weights_NormCCMEC  = evtwgt_map.at("NormCCMEC_Genie");
	      const std::vector<double> &weights_NormNCMEC  = evtwgt_map.at("NormNCMEC_Genie");
	      
	      _mc_wgt= weights.front();
	      _mc_wgt_cv= weights_cv.front();

	      std::cout<<"Value of mc_wgt_cv: "<<_mc_wgt_cv<<std::endl;

	      _mc_wgt_0_EtaNCEL = weights_EtaNCEL.at( 0 );
	      _mc_wgt_1_EtaNCEL = weights_EtaNCEL.at( 1 );
	      _mc_wgt_0_FrAbs_N = weights_FrAbs_N.at( 0 );
	      _mc_wgt_1_FrAbs_N = weights_FrAbs_N.at( 1 );
	      _mc_wgt_0_FrAbs_pi = weights_FrAbs_pi.at( 0 );
	      _mc_wgt_1_FrAbs_pi = weights_FrAbs_pi.at( 1 );
	      _mc_wgt_0_FrCEx_N = weights_FrCEx_N.at( 0 );
	      _mc_wgt_1_FrCEx_N = weights_FrCEx_N.at( 1 );
	      _mc_wgt_0_FrCEx_pi = weights_FrCEx_pi.at( 0 );
	      _mc_wgt_1_FrCEx_pi = weights_FrCEx_pi.at( 1 );
	      _mc_wgt_0_FrInel_N = weights_FrInel_N.at( 0 );
	      _mc_wgt_1_FrInel_N = weights_FrInel_N.at( 1 );
	      _mc_wgt_0_FrInel_pi = weights_FrInel_pi.at( 0 );
	      _mc_wgt_1_FrInel_pi = weights_FrInel_pi.at( 1 );
	      _mc_wgt_0_FrPiProd_N = weights_FrPiProd_N.at( 0 );
	      _mc_wgt_1_FrPiProd_N = weights_FrPiProd_N.at( 1 );
	      _mc_wgt_0_FrPiProd_pi = weights_FrPiProd_pi.at( 0 );
	      _mc_wgt_1_FrPiProd_pi = weights_FrPiProd_pi.at( 1 );
	      _mc_wgt_0_MFP_N = weights_MFP_N.at( 0 );
	      _mc_wgt_1_MFP_N = weights_MFP_N.at( 1 );
	      _mc_wgt_0_MFP_pi = weights_MFP_pi.at( 0 );
	      _mc_wgt_1_MFP_pi = weights_MFP_pi.at( 1 );
	      _mc_wgt_0_MaCCQE = weights_MaCCQE.at( 0 );
	      _mc_wgt_1_MaCCQE = weights_MaCCQE.at( 1 );
	      _mc_wgt_0_MaCCRES = weights_MaCCRES.at( 0 );
	      _mc_wgt_1_MaCCRES = weights_MaCCRES.at( 1 );
	      _mc_wgt_0_MaNCEL = weights_MaNCEL.at( 0 );
	      _mc_wgt_1_MaNCEL = weights_MaNCEL.at( 1 );
	      _mc_wgt_0_MaNCRES = weights_MaNCRES.at( 0 );
	      _mc_wgt_1_MaNCRES = weights_MaNCRES.at( 1 );
	      _mc_wgt_0_MvCCRES = weights_MvCCRES.at( 0 );
	      _mc_wgt_1_MvCCRES = weights_MvCCRES.at( 1 );
	      _mc_wgt_0_MvNCRES = weights_MvNCRES.at( 0 );
	      _mc_wgt_1_MvNCRES = weights_MvNCRES.at( 1 );
	      _mc_wgt_0_NonRESBGvnCC1pi = weights_NonRESBGvnCC1pi.at( 0 );
	      _mc_wgt_1_NonRESBGvnCC1pi = weights_NonRESBGvnCC1pi.at( 1 );
	      _mc_wgt_0_NonRESBGvnCC2pi = weights_NonRESBGvnCC2pi.at( 0 );
	      _mc_wgt_1_NonRESBGvnCC2pi = weights_NonRESBGvnCC2pi.at( 1 );
	      _mc_wgt_0_NonRESBGvnNC1pi = weights_NonRESBGvnNC1pi.at( 0 );
	      _mc_wgt_1_NonRESBGvnNC1pi = weights_NonRESBGvnNC1pi.at( 1 );
	      _mc_wgt_0_NonRESBGvnNC2pi = weights_NonRESBGvnNC2pi.at( 0 );
	      _mc_wgt_1_NonRESBGvnNC2pi = weights_NonRESBGvnNC2pi.at( 1 );
	      _mc_wgt_0_NonRESBGvpCC1pi = weights_NonRESBGvpCC1pi.at( 0 );
	      _mc_wgt_1_NonRESBGvpCC1pi = weights_NonRESBGvpCC1pi.at( 1 );
	      _mc_wgt_0_NonRESBGvpCC2pi = weights_NonRESBGvpCC2pi.at( 0 );
	      _mc_wgt_1_NonRESBGvpCC2pi = weights_NonRESBGvpCC2pi.at( 1 );
	      _mc_wgt_0_NonRESBGvpNC1pi = weights_NonRESBGvpNC1pi.at( 0 );
	      _mc_wgt_1_NonRESBGvpNC1pi = weights_NonRESBGvpNC1pi.at( 1 );
	      _mc_wgt_0_NonRESBGvpNC2pi = weights_NonRESBGvpNC2pi.at( 0 );
	      _mc_wgt_1_NonRESBGvpNC2pi = weights_NonRESBGvpNC2pi.at( 1 );
	      _mc_wgt_0_NormCCMEC = weights_NormCCMEC.at( 0 );
	      _mc_wgt_1_NormCCMEC = weights_NormCCMEC.at( 1 );
	      _mc_wgt_0_NormNCMEC = weights_NormNCMEC.at( 0 );
	      _mc_wgt_1_NormNCMEC = weights_NormNCMEC.at( 1 );
	      
	      
	      if(_debug) 
		{
		  std::cout << "[Numu0pi2p] v17 Event Weight:  " << _mc_wgt << std::endl;
		  std::cout << "[Numu0pi2p] cv:  " << weights_cv.front() << std::endl;
		  std::cout << "[Numu0pi2p] Eta:  " << weights_EtaNCEL.front() << std::endl;
		  for ( size_t u = 0; u < weights_EtaNCEL.size(); ++u ) 
		    {
		      double w = weights_EtaNCEL.at( u );
		      std::cout << "    universe #" << u << " has weight = " << w << '\n';
		    }	  
		}
	    }
	  else
	    {
	      _mc_wgt=1;
	      if(_debug) std::cout << "[Numu0pi2p] v17 Event Weight: Failed obtaining eventweight" << std::endl;
	    }
	}
      
      art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
      std::vector<art::Ptr<simb::MCTruth> > mclist;
      std::vector<art::Ptr<simb::MCParticle> > mcparticle;
      if (e.getByLabel(m_genie_producer,mctruthListHandle))
	{
	  art::fill_ptr_vector(mclist, mctruthListHandle);
	}
      if (e.getByLabel(m_geant_producer,mcparticleHandle))
	{
	  art::fill_ptr_vector(mcparticle, mcparticleHandle);
	}
      
      // Get space charge correction
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      
      // Get time offset for x space charge correction
      auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
      auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
      auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
      auto gen = mct_h->at(0);
      double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
      double xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
      if(!mclist.empty())
	{
	  if(mclist[0]->NeutrinoSet())
	    {
	      art::Ptr<simb::MCTruth> mctruth = mclist[0];
	      _mc_ccnc   = mctruth->GetNeutrino().CCNC();
	      _mc_mode   = mctruth->GetNeutrino().Mode();
	      _mc_hitnuc = mctruth->GetNeutrino().HitNuc();
	      _mc_q2     = mctruth->GetNeutrino().QSqr();
	      _mc_X = mctruth->GetNeutrino().X();
	      _mc_Y = mctruth->GetNeutrino().Y();
	      _mc_Pt = mctruth->GetNeutrino().Pt();
	      _mc_nu_vtxx    = mctruth->GetNeutrino().Nu().Vx();
	      _mc_nu_vtxy    = mctruth->GetNeutrino().Nu().Vy();
	      _mc_nu_vtxz    = mctruth->GetNeutrino().Nu().Vz();
	      _mc_enu    = mctruth->GetNeutrino().Nu().E();
	      _mc_interactiontype= mctruth->GetNeutrino().InteractionType();
	      _mc_nupdg  = mctruth->GetNeutrino().Nu().PdgCode();
	      
	      
	      auto sce_offset = SCE->GetPosOffsets(geo::Point_t(_mc_nu_vtxx,_mc_nu_vtxy,_mc_nu_vtxz));

	      _mc_nu_vtxx_sce = (_mc_nu_vtxx - sce_offset.X() + xtimeoffset) + 0.6;
	      _mc_nu_vtxy_sce = _mc_nu_vtxy + sce_offset.Y();
	      _mc_nu_vtxz_sce = _mc_nu_vtxz + sce_offset.Z();
	      
	      if(_debug)
		{
		  std::cout<<"[Numu0pi2p] CCNC: "<<_mc_ccnc<<std::endl;
		  std::cout<<"[Numu0pi2p] Mode: "<<_mc_mode<<std::endl;
		  std::cout<<"[Numu0pi2p] Hitnuc: "<<_mc_hitnuc<<std::endl;
		  std::cout<<"[Numu0pi2p] Q2: "<<_mc_q2<<std::endl;
		  std::cout<<"[Numu0pi2p] Bjorken X"<<_mc_X<<std::endl;
		  std::cout<<"[Numu0pi2p] Bjorken X"<<_mc_Y<<std::endl;
		  std::cout<<"[Numu0pi2p] Transverse Momentum"<<_mc_Pt<<std::endl;
		  std::cout<<"[Numu0pi2p] nu (x,y,z):  ("<<_mc_nu_vtxx<<", "<<_mc_nu_vtxy<<", "<<_mc_nu_vtxz<<")"<<std::endl;
		  std::cout<<"[Numu0pi2p] Ev:  "<<_mc_enu<<std::endl;
		}//if debug
	    }//if neutrino
	}//if mclist not empty
      /// count particles
      int p_mc_n_muon=0;
      int p_mc_n_proton=0;
      int p_mc_n_electron=0;
      int p_mc_n_pionpm=0;
      int p_mc_n_pion0=0;
      int p_mc_n_neutron=0;
      int p_mc_n_threshold_muon=0;
      int p_mc_n_threshold_proton=0;
      int p_mc_n_threshold_electron=0;
      int p_mc_n_threshold_pionpm=0;
      int p_mc_n_threshold_pion0=0;
      int p_mc_n_threshold_neutron=0;
      for(int i_mcp = 0; i_mcp < (int) mcparticle.size(); i_mcp++)
	{
	  if(mcparticle[i_mcp]->Process() == "primary")
	    {
	      _mc_g4_pdg.push_back(mcparticle[i_mcp]->PdgCode());
	      _mc_g4_E.push_back(mcparticle[i_mcp]->E());
	      _mc_g4_p.push_back(mcparticle[i_mcp]->P());
	      _mc_g4_mass.push_back(mcparticle[i_mcp]->Mass());
	      _mc_g4_phi.push_back(mcparticle[i_mcp]->Momentum().Phi());
	      _mc_g4_theta.push_back(mcparticle[i_mcp]->Momentum().Theta());
	      _mc_g4_mom_all.push_back(mcparticle[i_mcp]->P());
	      _mc_g4_start_x.push_back(mcparticle[i_mcp]->Vx());
	      _mc_g4_start_y.push_back(mcparticle[i_mcp]->Vy());
	      _mc_g4_start_z.push_back(mcparticle[i_mcp]->Vz());
	      _mc_g4_end_x.push_back(mcparticle[i_mcp]->EndX());
	      _mc_g4_end_y.push_back(mcparticle[i_mcp]->EndY());
	      _mc_g4_end_z.push_back(mcparticle[i_mcp]->EndZ());
	      auto start_sce_offset = SCE->GetPosOffsets(geo::Point_t(mcparticle[i_mcp]->Vx(),mcparticle[i_mcp]->Vy(),mcparticle[i_mcp]->Vz()));
              auto end_sce_offset = SCE->GetPosOffsets(geo::Point_t(mcparticle[i_mcp]->EndX(),mcparticle[i_mcp]->EndY(),mcparticle[i_mcp]->EndZ()));
	      _mc_g4_start_x_sce.push_back(mcparticle[i_mcp]->Vx() - start_sce_offset.X() + xtimeoffset + 0.6);
	      _mc_g4_start_y_sce.push_back(mcparticle[i_mcp]->Vy() + start_sce_offset.Y());
	      _mc_g4_start_z_sce.push_back(mcparticle[i_mcp]->Vz() + start_sce_offset.Z());
	      _mc_g4_end_x_sce.push_back(mcparticle[i_mcp]->EndX()  - end_sce_offset.X() + xtimeoffset + 0.6);
	      _mc_g4_end_y_sce.push_back(mcparticle[i_mcp]->EndY() + end_sce_offset.Y());
	      _mc_g4_end_z_sce.push_back(mcparticle[i_mcp]->EndZ() + end_sce_offset.Z());

             
	      if(abs(mcparticle[i_mcp]->PdgCode()) == 13) 
		{
		  p_mc_n_muon++;
		  _mc_g4_mom_muon.push_back(mcparticle[i_mcp]->P());
		  if(mcparticle[i_mcp]->P() >= m_p_muon)
		    {
		      p_mc_n_threshold_muon++;
		    }
		}
	      if(abs(mcparticle[i_mcp]->PdgCode()) == 11) 
		{
		  p_mc_n_electron++;
		  _mc_g4_mom_electron.push_back(mcparticle[i_mcp]->P());
		   if(mcparticle[i_mcp]->P() >= m_p_electron)
		    {
		      p_mc_n_threshold_electron++;
		    }
		}
	      if(abs(mcparticle[i_mcp]->PdgCode()) == 2112)
		{
		  p_mc_n_neutron++;
		  _mc_g4_mom_neutron.push_back(mcparticle[i_mcp]->P());
		  if(mcparticle[i_mcp]->P() >= m_p_neutron)
		    {
		      p_mc_n_threshold_neutron++;
		    }
		}

	      if(abs(mcparticle[i_mcp]->PdgCode()) == 2212)
		{
		  p_mc_n_proton++;
		   _mc_g4_mom_proton.push_back(mcparticle[i_mcp]->P());
		    if(mcparticle[i_mcp]->P() >= m_p_proton)
		    {
		      p_mc_n_threshold_proton++;
		    }
		}
	      if(abs(mcparticle[i_mcp]->PdgCode()) == 111) 
		{
		  p_mc_n_pion0++;
		   _mc_g4_mom_pion0.push_back(mcparticle[i_mcp]->P());
		    if(mcparticle[i_mcp]->P() >= m_p_pion0)
		    {
		      p_mc_n_threshold_pion0++;
		    }
		}
	      if(abs(mcparticle[i_mcp]->PdgCode()) == 211)
		{
		  p_mc_n_pionpm++;
		   _mc_g4_mom_pionpm.push_back(mcparticle[i_mcp]->P());
		   if(mcparticle[i_mcp]->P() >= m_p_pionpm)
		     {
		       p_mc_n_threshold_pionpm++;
		     }
		}
	    }// is primary mcparticle
	}//loop over mcparticle
      ///end of counting particles
      _mc_n_muon= p_mc_n_muon;
      _mc_n_proton= p_mc_n_proton;
      _mc_n_electron= p_mc_n_electron;
      _mc_n_pionpm= p_mc_n_pionpm;
      _mc_n_pion0= p_mc_n_pion0;
      _mc_n_neutron= p_mc_n_neutron;
      _mc_n_threshold_muon= p_mc_n_threshold_muon;
      _mc_n_threshold_proton= p_mc_n_threshold_proton;
      _mc_n_threshold_electron= p_mc_n_threshold_electron;
      _mc_n_threshold_pionpm= p_mc_n_threshold_pionpm;
      _mc_n_threshold_pion0= p_mc_n_threshold_pion0;
      _mc_n_threshold_neutron= p_mc_n_threshold_neutron;
       if(_debug)
	 {
	   std::cout<<"[Numu0pi2p] # of muon: "<<_mc_n_muon<<std::endl;
	   std::cout<<"[Numu0pi2p] # of proton: "<<_mc_n_proton<<std::endl;
	   std::cout<<"[Numu0pi2p] # of electron: "<<_mc_n_electron<<std::endl;
	   std::cout<<"[Numu0pi2p] # of pion+-: "<<_mc_n_pionpm<<std::endl;
	   std::cout<<"[Numu0pi2p] # of pion0: "<<_mc_n_pion0<<std::endl;
	   std::cout<<"[Numu0pi2p] # of neutron: "<<_mc_n_neutron<<std::endl;
	   
	   std::cout<<"[Numu0pi2p] # of threshold muon: "<<_mc_n_threshold_muon<<std::endl;
	   std::cout<<"[Numu0pi2p] # of threshold proton: "<<_mc_n_threshold_proton<<std::endl;
	   std::cout<<"[Numu0pi2p] # of threshold electron: "<<_mc_n_threshold_electron<<std::endl;
	   std::cout<<"[Numu0pi2p] # of threshold pion+-: "<<_mc_n_threshold_pionpm<<std::endl;
	   std::cout<<"[Numu0pi2p] # of threshold pion0: "<<_mc_n_threshold_pion0<<std::endl;
	   std::cout<<"[Numu0pi2p] # of threshold neutron: "<<_mc_n_threshold_neutron<<std::endl;


	 }//if debug
    }
} //End of Fill GENIE Truth

bool CC2p::IsNearVertex(float st_x, float st_y, float st_z, float end_x, float end_y, float end_z)
{
  if((st_x - end_x)*(st_x - end_x)+(st_y - end_y)*(st_y - end_y)+(st_z - end_z)*(st_z - end_z)< 25) return 1;
  else return 0;
}

float CC2p::KEfromLength(float trk_length){
   const double A=31.3;
   const double bp1=0.578;
   float KE= A* TMath::Power(trk_length,bp1);
   return KE;
 }

bool CC2p::IsContained(float st_x, float st_y, float st_z, float end_x, float end_y, float end_z)
{
  const float xmin=0;
  const float xmax=256.35;
  const float ymin=-116.35;
  const float ymax=116.35;
  const float zmin=0;
  const float zmax=1036.8;
  const float border_xorz=10;
  const float border_y=20;
  
  if(st_x < xmin+ border_xorz) return false;
  if(end_x< xmin + border_xorz) return false;
  if(st_y< ymin + border_y) return false;
  if(end_y<ymin + border_y) return false;
  if(st_z< zmin + border_xorz) return false;
  if(end_z< zmin + border_xorz) return false;
  
  if(st_x > xmax - border_xorz) return false;
  if(end_x > xmax - border_xorz) return false;
  if(st_y > ymax - border_y) return false;
  if(end_y > ymax - border_y) return false;
  if(st_z > zmax - border_xorz) return false;
  if(end_z > zmax - border_xorz) return false;
  return true;
}


void CC2p::beginJob()
{
  // Implementation of optional member function here.
}

//void CC2p::endJob()
//{
  // Implementation of optional member function here.
//}

DEFINE_ART_MODULE(CC2p)
