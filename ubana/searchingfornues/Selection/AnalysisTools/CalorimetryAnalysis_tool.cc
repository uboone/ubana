#ifndef ANALYSIS_CALORIMETRYANALYSIS_CXX
#define ANALYSIS_CALORIMETRYANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Typedefs.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"

#include "lardataobj/AnalysisBase/T0.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CalorimetryAnalysis
// File:        CalorimetryAnalysis.cc
//
//              A basic analysis example
//
// Configuration parameters:
//
// TBD
//
// Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
//
////////////////////////////////////////////////////////////////////////

class CalorimetryAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  CalorimetryAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~CalorimetryAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief Fill Default info for event associated to neutrino
     */
  void fillDefault();

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;


private:

  // function that given a track and its calo info fills NTuple variables
  void FillCalorimetry(art::Event const &e,
           const searchingfornues::ProxyPfpElem_t pfp,
           const searchingfornues::ProxyCaloColl_t calo_proxy,
           const searchingfornues::ProxyPIDColl_t pid_proxy,
           const searchingfornues::ProxyClusColl_t clus_proxy,
           const TVector3 nu_vtx,
           const bool fData,
           const bool fShrFit,
           const float fEnergyThresholdForMCHits,
           const std::vector<searchingfornues::BtPart> btparts_v,
           const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart);


  const trkf::TrackMomentumCalculator _trkmom;
  const trkf::TrajectoryMCSFitter _mcsfitter;

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fPFPproducer;
  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fCLSproducer;
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  float fEnergyThresholdForMCHits;
  art::InputTag fMCRproducer;
  art::InputTag fMCTproducer;
  art::InputTag fT0producer;
  art::InputTag fMCPproducer;

  bool fShrFit; // use shower track-fitter info?

  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]

  bool fGetCaloID; // get the index of the calorimetry object manually. Needs to be true unless object produced by Pandora hierarchy
  // (This is at least what I foudn empirically)

  TTree* _calo_tree;

  int _run, _sub, _evt;

  float _true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z;
  float _true_nu_vtx_sce_x, _true_nu_vtx_sce_y, _true_nu_vtx_sce_z;
  float _reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z;
  float _reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z;

  // backtracking information
  int _backtracked_pdg;            // PDG code of backtracked particle
  float _backtracked_e;            // energy of backtracked particle
  float _backtracked_purity;       // purity of backtracking
  float _backtracked_completeness; // completeness of backtracking
  float _backtracked_overlay_purity; // purity of overlay

  float _backtracked_px;
  float _backtracked_py;
  float _backtracked_pz;

  float _backtracked_start_x;
  float _backtracked_start_y;
  float _backtracked_start_z;
  float _backtracked_start_t;
  float _backtracked_start_U;
  float _backtracked_start_V;
  float _backtracked_start_Y;
  float _backtracked_sce_start_x;
  float _backtracked_sce_start_y;
  float _backtracked_sce_start_z;
  float _backtracked_sce_start_U;
  float _backtracked_sce_start_V;
  float _backtracked_sce_start_Y;

  std::string _backtracked_end_process;
  bool _backtracked_end_in_tpc;

  uint _generation;    // generation, 1 is primary
  uint _shr_daughters; // number of shower daughters
  uint _trk_daughters; // number of track daughters

  // track information
  int _nplanehits_U;
  int _nplanehits_V;
  int _nplanehits_Y;
  float _trk_score;

  float _trk_theta;
  float _trk_phi;

  float _trk_dir_x;
  float _trk_dir_y;
  float _trk_dir_z;

  float _trk_len;
  float _trk_distance;

  float _trk_start_x;
  float _trk_start_y;
  float _trk_start_z;

  float _trk_sce_start_x;
  float _trk_sce_start_y;
  float _trk_sce_start_z;

  float _trk_end_x;
  float _trk_end_y;
  float _trk_end_z;

  float _trk_sce_end_x;
  float _trk_sce_end_y;
  float _trk_sce_end_z;

  float _trk_bragg_p_y;
  float _trk_bragg_mu_y;
  float _trk_bragg_mip_y;
  float _trk_pid_chipr_y;
  float _trk_pid_chika_y;
  float _trk_pid_chipi_y;
  float _trk_pid_chimu_y;
  float _trk_pida_y;

  float _trk_bragg_p_u;
  float _trk_bragg_mu_u;
  float _trk_bragg_mip_u;
  float _trk_pid_chipr_u;
  float _trk_pid_chika_u;
  float _trk_pid_chipi_u;
  float _trk_pid_chimu_u;
  float _trk_pida_u;

  float _trk_bragg_p_v;
  float _trk_bragg_mu_v;
  float _trk_bragg_mip_v;
  float _trk_pid_chipr_v;
  float _trk_pid_chika_v;
  float _trk_pid_chipi_v;
  float _trk_pid_chimu_v;
  float _trk_pida_v;

  float _trk_bragg_p_three_planes;

  float _trk_mcs_muon_mom;
  float _trk_range_muon_mom;
  float _trk_energy_proton;
  float _trk_energy_muon;

  int _longest; // longest track in slice?

  // dedx vector
  std::vector<float> _dqdx_u;
  std::vector<float> _dqdx_v;
  std::vector<float> _dqdx_y;

  std::vector<float> _dedx_u;
  std::vector<float> _dedx_v;
  std::vector<float> _dedx_y;

  std::vector<float> _rr_u;
  std::vector<float> _rr_v;
  std::vector<float> _rr_y;

  std::vector<float> _pitch_u;
  std::vector<float> _pitch_v;
  std::vector<float> _pitch_y;

  std::vector<float> _x_u;
  std::vector<float> _x_v;
  std::vector<float> _x_y;

  std::vector<float> _y_u;
  std::vector<float> _y_v;
  std::vector<float> _y_y;

  std::vector<float> _z_u;
  std::vector<float> _z_v;
  std::vector<float> _z_y;

  std::vector<float> _dir_x_u;
  std::vector<float> _dir_x_v;
  std::vector<float> _dir_x_y;

  std::vector<float> _dir_y_u;
  std::vector<float> _dir_y_v;
  std::vector<float> _dir_y_y;

  std::vector<float> _dir_z_u;
  std::vector<float> _dir_z_v;
  std::vector<float> _dir_z_y;

  std::vector<bool> _is_hit_montecarlo_u;
  std::vector<bool> _is_hit_montecarlo_v;
  std::vector<bool> _is_hit_montecarlo_y;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CalorimetryAnalysis::CalorimetryAnalysis(const fhicl::ParameterSet &p) :
_mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(p.get<fhicl::ParameterSet>("mcsfitmu")))
{
  fPFPproducer  = p.get< art::InputTag > ("PFPproducer","pandora");
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );
  fT0producer  = p.get< art::InputTag > ("T0producer", "" );

  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
  fHproducer = p.get<art::InputTag>("Hproducer");
  fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fMCTproducer = p.get<art::InputTag>("MCTproducer");

  fShrFit    = p.get<bool>("ShrFit"   , false);
  fGetCaloID = p.get<bool>("GetCaloID", false);

  fADCtoE = p.get<std::vector<float>>("ADCtoE");

  art::ServiceHandle<art::TFileService> tfs;

  _calo_tree = tfs->make<TTree>("CalorimetryAnalyzer", "Calo Tree");
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CalorimetryAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CalorimetryAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[CalorimetryAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: "<< _evt << std::endl;

  // if flag to store cosmic muon info is on -> save ACPT tagged tracks if any
  if (fT0producer == "")
    return;

  // grab PFParticles in event
  searchingfornues::ProxyPfpColl_t const &pfp_proxy_v = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                               proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                               proxy::withAssociated<recob::Cluster>(fPFPproducer),
                               proxy::withAssociated<recob::Slice>(fPFPproducer),
                               proxy::withAssociated<recob::Track>(fPFPproducer),
                               proxy::withAssociated<recob::Vertex>(fPFPproducer),
                               proxy::withAssociated<recob::PCAxis>(fPFPproducer),
                               proxy::withAssociated<recob::Shower>(fPFPproducer),
                               proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

  // load clusters
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));

  // load backtrack information
  std::vector<searchingfornues::BtPart> btparts_v;
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

  if (fData)
  {
    const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
    const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
    /*
    auto const &hit_mcpart_assn = *(e.getValidHandle<art::Assns<recob::Hit,simb::MCParticle,anab::BackTrackerHitMatchingData>>(fBacktrackTag));

    std::unordered_set<size_t> my_hit_key_mc;

    int ctr = 0.;
    std::cout << "there are " << hit_mcpart_assn.size() << " hit_mcpart_assn" << std::endl;
    for(auto &pair : hit_mcpart_assn){

      std::cout << "hit " << pair.first.key() << " is associated to track-id " << pair.second->TrackId() << " " << std::endl;
      my_hit_key_mc.insert(pair.first.key());
      ctr += 1;

    }
    std::cout << " ctr is " << ctr << std::endl;

    for (size_t h=0; h < inputHits->size(); h++) {
      if(my_hit_key_mc.count(h) > 0)
  std::cout << "hit backtracked!" << std::endl;
    }
    */
  }

  _true_nu_vtx_t = std::numeric_limits<float>::lowest();
  _true_nu_vtx_x = std::numeric_limits<float>::lowest();
  _true_nu_vtx_y = std::numeric_limits<float>::lowest();
  _true_nu_vtx_z = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();



  // load calorimetry
  searchingfornues::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const& pid_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::ParticleID>(fPIDproducer));

  // grab tracks in the event
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track> >(fPFPproducer);
  // grab anab::T0 for ACPT in-time associated to Tracks
  art::FindManyP<anab::T0> trk_t0_assn_v(trk_h, e, fT0producer);

  for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy_v) {

    fillDefault();

    // grab associated tracks
    auto trk_v = pfp_pxy.get<recob::Track>();
    if (trk_v.size() != 1)
      continue;

    auto trk = trk_v.at(0);

    // grab track key to find anab::T0 association, if present
    auto const T0_v = trk_t0_assn_v.at( trk.key() );
    if (T0_v.size() == 1)
    {
      TVector3 nuvtx;
      nuvtx.SetXYZ(0, 0, 0);
      FillCalorimetry(e,
          pfp_pxy,
          calo_proxy,
          pid_proxy,
          clus_proxy,
          nuvtx,
          fData,
          fShrFit,
          fEnergyThresholdForMCHits,
          btparts_v,
    assocMCPart);

    }//  if T0 assocaition is found
  }// for all PFParticles
}// analyzeEvent

void CalorimetryAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  // if fT0producer was specified -> means we don't want to run on neutrino candidate tracks -> skip
  if (fT0producer != "")
    return;

  // load clusters
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));

  // load backtrack information
  std::vector<searchingfornues::BtPart> btparts_v;
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
  if (!fData)
  {
    const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
    const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);

    // save truth
    // load MCTruth
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

    auto mct = mct_h->at(0);
    auto neutrino = mct.GetNeutrino();
    auto nu = neutrino.Nu();

    _true_nu_vtx_t = nu.T();
    _true_nu_vtx_x = nu.Vx();
    _true_nu_vtx_y = nu.Vy();
    _true_nu_vtx_z = nu.Vz();

    float _true_nu_vtx_sce[3];
    searchingfornues::True2RecoMappingXYZ(_true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z, _true_nu_vtx_sce);

    _true_nu_vtx_sce_x = _true_nu_vtx_sce[0];
    _true_nu_vtx_sce_y = _true_nu_vtx_sce[1];
    _true_nu_vtx_sce_z = _true_nu_vtx_sce[2];
  }
  else
  {
    _true_nu_vtx_t = std::numeric_limits<float>::lowest();
    _true_nu_vtx_x = std::numeric_limits<float>::lowest();
    _true_nu_vtx_y = std::numeric_limits<float>::lowest();
    _true_nu_vtx_z = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
  }

  // load calorimetry
  searchingfornues::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const& pid_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::ParticleID>(fPIDproducer));


  // Build larpandora info:
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleMap particleMap;
  larpandora.CollectPFParticles(e, "pandora", pfparticles);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  //find vertex
  TVector3 nuvtx;
  for (auto pfp : slice_pfp_v)
  {
    if (pfp->IsPrimary())
    {
      double xyz[3] = {};
      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      else
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);

        _reco_nu_vtx_x = nuvtx.X();
        _reco_nu_vtx_y = nuvtx.Y();
        _reco_nu_vtx_z = nuvtx.Z();

        float _reco_nu_vtx_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(_reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z, _reco_nu_vtx_sce);
        _reco_nu_vtx_sce_x = _reco_nu_vtx_sce[0];
        _reco_nu_vtx_sce_y = _reco_nu_vtx_sce[1];
        _reco_nu_vtx_sce_z = _reco_nu_vtx_sce[2];
      }
      break;
    }
  }


  // find longest track
  float  lenmax = 0;
  size_t idxmax = 0;
  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto pfp = slice_pfp_v.at(i_pfp);

    if ( pfp->IsPrimary() )
      continue;

    auto trk_v = pfp.get<recob::Track>();
    if (trk_v.size() != 1)
      continue;

    auto trk = trk_v.at(0);

    float len = trk->Length();
    if (len > lenmax)
    {
      lenmax = len;
      idxmax = i_pfp;
    }
  }// for all pfps


  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    fillDefault();
    auto pfp = slice_pfp_v[i_pfp];

    // Hieracrchy information:
    _generation = larpandora.GetGeneration(particleMap, particleMap.at(pfp->Self()));
    uint this_num_trk_d = 0;
    uint this_num_shr_d = 0;
    for (size_t daughter : pfp->Daughters())
    {
      if (larpandora.IsTrack(particleMap.at(daughter)))
      {
        this_num_trk_d++; // Track daughter
      }
      else
      {
        this_num_shr_d++; // Shower daughter
      }
    }
    _shr_daughters = this_num_shr_d;
    _trk_daughters = this_num_trk_d;

    _longest = 0;
    if ((i_pfp == idxmax) && (lenmax != 0))
      _longest = 1;

    FillCalorimetry(e,
        pfp,
        calo_proxy,
        pid_proxy,
        clus_proxy,
        nuvtx,
        fData,
        fShrFit,
        fEnergyThresholdForMCHits,
        btparts_v,
        assocMCPart);
  } // for all PFParticles
}

void CalorimetryAnalysis::fillDefault()
{
  // backtracking information
  _backtracked_pdg = std::numeric_limits<int>::lowest();            // PDG code of backtracked particle
  _backtracked_e = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_purity = std::numeric_limits<float>::lowest();       // purity of backtracking
  _backtracked_completeness = std::numeric_limits<float>::lowest(); // completeness of backtracking
  _backtracked_overlay_purity = std::numeric_limits<float>::lowest(); // purity of overlay
  _backtracked_end_process = "";
  _backtracked_end_in_tpc = false;

  _backtracked_px = std::numeric_limits<float>::lowest();
  _backtracked_py = std::numeric_limits<float>::lowest();
  _backtracked_pz = std::numeric_limits<float>::lowest();

  _backtracked_start_x = std::numeric_limits<float>::lowest();
  _backtracked_start_y = std::numeric_limits<float>::lowest();
  _backtracked_start_z = std::numeric_limits<float>::lowest();
  _backtracked_start_t = std::numeric_limits<float>::lowest();
  _backtracked_start_U = std::numeric_limits<float>::lowest();
  _backtracked_start_V = std::numeric_limits<float>::lowest();
  _backtracked_start_Y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_x = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_z = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_U = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_V = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_Y = std::numeric_limits<float>::lowest();

  _generation = std::numeric_limits<uint>::lowest();
  _shr_daughters = std::numeric_limits<uint>::lowest();
  _trk_daughters = std::numeric_limits<uint>::lowest();

  // track information
  _nplanehits_U = std::numeric_limits<int>::lowest();
  _nplanehits_V = std::numeric_limits<int>::lowest();
  _nplanehits_Y = std::numeric_limits<int>::lowest();
  _trk_score = std::numeric_limits<float>::lowest();

  _trk_theta = std::numeric_limits<float>::lowest();
  _trk_phi = std::numeric_limits<float>::lowest();
  _trk_len = std::numeric_limits<float>::lowest();
  _trk_distance = std::numeric_limits<float>::lowest();

  _trk_dir_x = std::numeric_limits<float>::lowest();
  _trk_dir_y = std::numeric_limits<float>::lowest();
  _trk_dir_z = std::numeric_limits<float>::lowest();

  _trk_start_x = std::numeric_limits<float>::lowest();
  _trk_start_y = std::numeric_limits<float>::lowest();
  _trk_start_z = std::numeric_limits<float>::lowest();

  _trk_sce_start_x = std::numeric_limits<float>::lowest();
  _trk_sce_start_y = std::numeric_limits<float>::lowest();
  _trk_sce_start_z = std::numeric_limits<float>::lowest();

  _trk_end_x = std::numeric_limits<float>::lowest();
  _trk_end_y = std::numeric_limits<float>::lowest();
  _trk_end_z = std::numeric_limits<float>::lowest();

  _trk_sce_end_x = std::numeric_limits<float>::lowest();
  _trk_sce_end_y = std::numeric_limits<float>::lowest();
  _trk_sce_end_z = std::numeric_limits<float>::lowest();

  _trk_bragg_p_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_y = std::numeric_limits<float>::lowest();
  _trk_pid_chika_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_y = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_y = std::numeric_limits<float>::lowest();
  _trk_pida_y = std::numeric_limits<float>::lowest();

  _trk_bragg_p_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_u = std::numeric_limits<float>::lowest();
  _trk_pid_chika_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_u = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_u = std::numeric_limits<float>::lowest();
  _trk_pida_u = std::numeric_limits<float>::lowest();

  _trk_bragg_p_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_v = std::numeric_limits<float>::lowest();
  _trk_pid_chika_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_v = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_v = std::numeric_limits<float>::lowest();
  _trk_pida_v = std::numeric_limits<float>::lowest();

  _trk_bragg_p_three_planes = std::numeric_limits<float>::lowest();

  _longest = std::numeric_limits<int>::lowest();

  _trk_mcs_muon_mom = std::numeric_limits<float>::lowest();
  _trk_range_muon_mom = std::numeric_limits<float>::lowest();
  _trk_energy_proton = std::numeric_limits<float>::lowest();
  _trk_energy_muon = std::numeric_limits<float>::lowest();

  // dedx vector
  _dqdx_u.clear();
  _dqdx_v.clear();
  _dqdx_y.clear();

  _dedx_u.clear();
  _dedx_v.clear();
  _dedx_y.clear();

  _rr_u.clear();
  _rr_v.clear();
  _rr_y.clear();

  _pitch_u.clear();
  _pitch_v.clear();
  _pitch_y.clear();

  _x_u.clear();
  _x_v.clear();
  _x_y.clear();

  _y_u.clear();
  _y_v.clear();
  _y_y.clear();

  _z_u.clear();
  _z_v.clear();
  _z_y.clear();

  _dir_x_u.clear();
  _dir_x_v.clear();
  _dir_x_y.clear();

  _dir_y_u.clear();
  _dir_y_v.clear();
  _dir_y_y.clear();

  _dir_z_u.clear();
  _dir_z_v.clear();
  _dir_z_y.clear();

  _is_hit_montecarlo_u.clear();
  _is_hit_montecarlo_v.clear();
  _is_hit_montecarlo_y.clear();
}

void CalorimetryAnalysis::setBranches(TTree *_tree)
{
  _calo_tree->Branch("run", &_run, "run/i");
  _calo_tree->Branch("sub", &_sub, "sub/i");
  _calo_tree->Branch("evt", &_evt, "evt/i");

  _calo_tree->Branch("true_nu_vtx_t", &_true_nu_vtx_t, "true_nu_vtx_t/F");
  _calo_tree->Branch("true_nu_vtx_x", &_true_nu_vtx_x, "true_nu_vtx_x/F");
  _calo_tree->Branch("true_nu_vtx_y", &_true_nu_vtx_y, "true_nu_vtx_y/F");
  _calo_tree->Branch("true_nu_vtx_z", &_true_nu_vtx_z, "true_nu_vtx_z/F");
  _calo_tree->Branch("true_nu_vtx_sce_x", &_true_nu_vtx_sce_x, "true_nu_vtx_sce_x/F");
  _calo_tree->Branch("true_nu_vtx_sce_y", &_true_nu_vtx_sce_y, "true_nu_vtx_sce_y/F");
  _calo_tree->Branch("true_nu_vtx_sce_z", &_true_nu_vtx_sce_z, "true_nu_vtx_sce_z/F");
  _calo_tree->Branch("reco_nu_vtx_x", &_reco_nu_vtx_x, "reco_nu_vtx_x/F");
  _calo_tree->Branch("reco_nu_vtx_y", &_reco_nu_vtx_y, "reco_nu_vtx_y/F");
  _calo_tree->Branch("reco_nu_vtx_z", &_reco_nu_vtx_z, "reco_nu_vtx_z/F");
  _calo_tree->Branch("reco_nu_vtx_sce_x", &_reco_nu_vtx_sce_x, "reco_nu_vtx_sce_x/F");
  _calo_tree->Branch("reco_nu_vtx_sce_y", &_reco_nu_vtx_sce_y, "reco_nu_vtx_sce_y/F");
  _calo_tree->Branch("reco_nu_vtx_sce_z", &_reco_nu_vtx_sce_z, "reco_nu_vtx_sce_z/F");
  // backtracking information
  _calo_tree->Branch("backtracked_pdg", &_backtracked_pdg, "backtracked_pdg/I");            // PDG code of backtracked particle
  _calo_tree->Branch("backtracked_e", &_backtracked_e, "backtracked_e/f");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_purity", &_backtracked_purity, "backtracked_purity/f");       // purity of backtracking
  _calo_tree->Branch("backtracked_completeness", &_backtracked_completeness, "backtracked_completeness/f"); // completeness of backtracking
  _calo_tree->Branch("backtracked_overlay_purity", &_backtracked_overlay_purity, "backtracked_overlay_purity/f"); // purity of overlay

  _calo_tree->Branch("backtracked_start_x", &_backtracked_start_x, "backtracked_start_x/f");
  _calo_tree->Branch("backtracked_start_y", &_backtracked_start_y, "backtracked_start_y/f");
  _calo_tree->Branch("backtracked_start_z", &_backtracked_start_z, "backtracked_start_z/f");
  _calo_tree->Branch("backtracked_start_t", &_backtracked_start_t, "backtracked_start_t/f");
  // _calo_tree->Branch("backtracked_start_U", &_backtracked_start_U, "backtracked_start_U/f");
  // _calo_tree->Branch("backtracked_start_V", &_backtracked_start_V, "backtracked_start_V/f");
  // _calo_tree->Branch("backtracked_start_Y", &_backtracked_start_Y, "backtracked_start_Y/f");
  _calo_tree->Branch("backtracked_sce_start_x", &_backtracked_sce_start_x, "backtracked_sce_start_x/f");
  _calo_tree->Branch("backtracked_sce_start_y", &_backtracked_sce_start_y, "backtracked_sce_start_y/f");
  _calo_tree->Branch("backtracked_sce_start_z", &_backtracked_sce_start_z, "backtracked_sce_start_z/f");
  // _calo_tree->Branch("backtracked_sce_start_U", &_backtracked_sce_start_U, "backtracked_sce_start_U/f");
  // _calo_tree->Branch("backtracked_sce_start_V", &_backtracked_sce_start_V, "backtracked_sce_start_V/f");
  // _calo_tree->Branch("backtracked_sce_start_Y", &_backtracked_sce_start_Y, "backtracked_sce_start_Y/f");

  _calo_tree->Branch("backtracked_end_process", &_backtracked_end_process);
  _calo_tree->Branch("backtracked_end_in_tpc", &_backtracked_end_in_tpc, "backtracked_end_in_tpc/O");

  _calo_tree->Branch("generation", &_generation, "generation/i");
  _calo_tree->Branch("trk_daughters", &_trk_daughters, "trk_daughters/i");
  _calo_tree->Branch("shr_daughters", &_shr_daughters, "shr_daughters/i");

  // track information
  _calo_tree->Branch("nplanehits_U", &_nplanehits_U, "nplanehits_U/I");
  _calo_tree->Branch("nplanehits_V", &_nplanehits_V, "nplanehits_V/I");
  _calo_tree->Branch("nplanehits_Y", &_nplanehits_Y, "nplanehits_Y/I");
  _calo_tree->Branch("trk_score", &_trk_score, "trk_score/f");

  _calo_tree->Branch("trk_theta", &_trk_theta, "trk_theta/f");
  _calo_tree->Branch("trk_phi", &_trk_phi, "trk_phi/f");
  _calo_tree->Branch("trk_len", &_trk_len, "trk_len/f");
  _calo_tree->Branch("trk_distance", &_trk_distance, "trk_distance/f");

  _calo_tree->Branch("trk_dir_x", &_trk_dir_x, "trk_dir_x/f");
  _calo_tree->Branch("trk_dir_y", &_trk_dir_y, "trk_dir_y/f");
  _calo_tree->Branch("trk_dir_z", &_trk_dir_z, "trk_dir_z/f");

  _calo_tree->Branch("longest", &_longest, "longest/I");

  _calo_tree->Branch("trk_start_x", &_trk_start_x, "trk_start_x/f");
  _calo_tree->Branch("trk_start_y", &_trk_start_y, "trk_start_y/f");
  _calo_tree->Branch("trk_start_z", &_trk_start_z, "trk_start_z/f");

  _calo_tree->Branch("trk_sce_start_x", &_trk_sce_start_x, "trk_sce_start_x/f");
  _calo_tree->Branch("trk_sce_start_y", &_trk_sce_start_y, "trk_sce_start_y/f");
  _calo_tree->Branch("trk_sce_start_z", &_trk_sce_start_z, "trk_sce_start_z/f");

  _calo_tree->Branch("trk_end_x", &_trk_end_x, "trk_end_x/f");
  _calo_tree->Branch("trk_end_y", &_trk_end_y, "trk_end_y/f");
  _calo_tree->Branch("trk_end_z", &_trk_end_z, "trk_end_z/f");

  _calo_tree->Branch("trk_sce_end_x", &_trk_sce_end_x, "trk_sce_end_x/f");
  _calo_tree->Branch("trk_sce_end_y", &_trk_sce_end_y, "trk_sce_end_y/f");
  _calo_tree->Branch("trk_sce_end_z", &_trk_sce_end_z, "trk_sce_end_z/f");

  _calo_tree->Branch("trk_bragg_p_u", &_trk_bragg_p_u, "trk_bragg_p_u/f");
  _calo_tree->Branch("trk_bragg_mu_u", &_trk_bragg_mu_u, "trk_bragg_mu_u/f");
  _calo_tree->Branch("trk_bragg_mip_u", &_trk_bragg_mip_u, "trk_bragg_mip_u/f");
  _calo_tree->Branch("trk_pid_chipr_u", &_trk_pid_chipr_u, "trk_pid_chipr_u/f");
  _calo_tree->Branch("trk_pid_chika_u", &_trk_pid_chika_u, "trk_pid_chika_u/f");
  _calo_tree->Branch("trk_pid_chipi_u", &_trk_pid_chipi_u, "trk_pid_chipi_u/f");
  _calo_tree->Branch("trk_pid_chimu_u", &_trk_pid_chimu_u, "trk_pid_chimu_u/f");
  _calo_tree->Branch("trk_pida_u", &_trk_pida_u, "trk_pida_u/f");

  _calo_tree->Branch("trk_bragg_p_v", &_trk_bragg_p_v, "trk_bragg_p_v/f");
  _calo_tree->Branch("trk_bragg_mu_v", &_trk_bragg_mu_v, "trk_bragg_mu_v/f");
  _calo_tree->Branch("trk_bragg_mip_v", &_trk_bragg_mip_v, "trk_bragg_mip_v/f");
  _calo_tree->Branch("trk_pid_chipr_v", &_trk_pid_chipr_v, "trk_pid_chipr_v/f");
  _calo_tree->Branch("trk_pid_chika_v", &_trk_pid_chika_v, "trk_pid_chika_v/f");
  _calo_tree->Branch("trk_pid_chipi_v", &_trk_pid_chipi_v, "trk_pid_chipi_v/f");
  _calo_tree->Branch("trk_pid_chimu_v", &_trk_pid_chimu_v, "trk_pid_chimu_v/f");
  _calo_tree->Branch("trk_pida_v", &_trk_pida_v, "trk_pida_v/f");

  _calo_tree->Branch("trk_bragg_p_y", &_trk_bragg_p_y, "trk_bragg_p_y/f");
  _calo_tree->Branch("trk_bragg_mu_y", &_trk_bragg_mu_y, "trk_bragg_mu_y/f");
  _calo_tree->Branch("trk_bragg_mip_y", &_trk_bragg_mip_y, "trk_bragg_mip_y/f");
  _calo_tree->Branch("trk_pid_chipr_y", &_trk_pid_chipr_y, "trk_pid_chipr_y/f");
  _calo_tree->Branch("trk_pid_chika_y", &_trk_pid_chika_y, "trk_pid_chika_y/f");
  _calo_tree->Branch("trk_pid_chipi_y", &_trk_pid_chipi_y, "trk_pid_chipi_y/f");
  _calo_tree->Branch("trk_pid_chimu_y", &_trk_pid_chimu_y, "trk_pid_chimu_y/f");
  _calo_tree->Branch("trk_pida_y", &_trk_pida_y, "trk_pida_y/f");

  _calo_tree->Branch("trk_bragg_p_three_planes", &_trk_bragg_p_three_planes, "trk_bragg_p_three_planes/f");

  _calo_tree->Branch("trk_mcs_muon_mom", &_trk_mcs_muon_mom, "trk_mcs_muon_mom/f");
  _calo_tree->Branch("trk_range_muon_mom", &_trk_range_muon_mom, "trk_range_muon_mom/f");
  _calo_tree->Branch("trk_energy_proton", &_trk_energy_proton, "trk_energy_proton/f");
  _calo_tree->Branch("trk_energy_muon", &_trk_energy_muon, "trk_energy_muon/f");

  // dedx vector
  _calo_tree->Branch("dqdx_u", "std::vector<float>", &_dqdx_u);
  _calo_tree->Branch("dqdx_v", "std::vector<float>", &_dqdx_v);
  _calo_tree->Branch("dqdx_y", "std::vector<float>", &_dqdx_y);

  _calo_tree->Branch("dedx_u", "std::vector<float>", &_dedx_u);
  _calo_tree->Branch("dedx_v", "std::vector<float>", &_dedx_v);
  _calo_tree->Branch("dedx_y", "std::vector<float>", &_dedx_y);

  _calo_tree->Branch("rr_u", "std::vector<float>", &_rr_u);
  _calo_tree->Branch("rr_v", "std::vector<float>", &_rr_v);
  _calo_tree->Branch("rr_y", "std::vector<float>", &_rr_y);

  _calo_tree->Branch("pitch_u", "std::vector<float>", &_pitch_u);
  _calo_tree->Branch("pitch_v", "std::vector<float>", &_pitch_v);
  _calo_tree->Branch("pitch_y", "std::vector<float>", &_pitch_y);

  _calo_tree->Branch("x_u", "std::vector<float>", &_x_u);
  _calo_tree->Branch("x_v", "std::vector<float>", &_x_v);
  _calo_tree->Branch("x_y", "std::vector<float>", &_x_y);

  _calo_tree->Branch("y_u", "std::vector<float>", &_y_u);
  _calo_tree->Branch("y_v", "std::vector<float>", &_y_v);
  _calo_tree->Branch("y_y", "std::vector<float>", &_y_y);

  _calo_tree->Branch("z_u", "std::vector<float>", &_z_u);
  _calo_tree->Branch("z_v", "std::vector<float>", &_z_v);
  _calo_tree->Branch("z_y", "std::vector<float>", &_z_y);

  _calo_tree->Branch("dir_x_u", "std::vector<float>", &_dir_x_u);
  _calo_tree->Branch("dir_x_v", "std::vector<float>", &_dir_x_v);
  _calo_tree->Branch("dir_x_y", "std::vector<float>", &_dir_x_y);

  _calo_tree->Branch("dir_y_u", "std::vector<float>", &_dir_y_u);
  _calo_tree->Branch("dir_y_v", "std::vector<float>", &_dir_y_v);
  _calo_tree->Branch("dir_y_y", "std::vector<float>", &_dir_y_y);

  _calo_tree->Branch("dir_z_u", "std::vector<float>", &_dir_z_u);
  _calo_tree->Branch("dir_z_v", "std::vector<float>", &_dir_z_v);
  _calo_tree->Branch("dir_z_y", "std::vector<float>", &_dir_z_y);

  _calo_tree->Branch("is_hit_montecarlo_u", "std::vector<bool>", &_is_hit_montecarlo_u);
  _calo_tree->Branch("is_hit_montecarlo_v", "std::vector<bool>", &_is_hit_montecarlo_v);
  _calo_tree->Branch("is_hit_montecarlo_y", "std::vector<bool>", &_is_hit_montecarlo_y);
}

void CalorimetryAnalysis::resetTTree(TTree *_tree)
{
}


void CalorimetryAnalysis::FillCalorimetry(art::Event const &e,
             const searchingfornues::ProxyPfpElem_t pfp,
            const searchingfornues::ProxyCaloColl_t calo_proxy,
            const searchingfornues::ProxyPIDColl_t pid_proxy,
            const searchingfornues::ProxyClusColl_t clus_proxy,
            const TVector3 nu_vtx,
            const bool fData,
            const bool fShrFit,
            const float fEnergyThresholdForMCHits,
            const std::vector<searchingfornues::BtPart> btparts_v,
            const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart)
{

  if ( (pfp->IsPrimary()) && (fT0producer == ""))
    return;

  auto trk_v = pfp.get<recob::Track>();
  if (trk_v.size() != 1)
    return;

  auto trk = trk_v.at(0);

  //add is primary daughter ?

  // store pfp information
  _trk_score = searchingfornues::GetTrackShowerScore(pfp);

  // get hits associated to this PFParticle through the clusters
  std::vector<art::Ptr<recob::Hit>> hit_v;
  auto clus_pxy_v = pfp.get<recob::Cluster>();
  for (auto ass_clus : clus_pxy_v)
  {
    // get cluster proxy
    const auto &clus = clus_proxy[ass_clus.key()];
    auto clus_hit_v = clus.get<recob::Hit>();
    auto nhits = clus_hit_v.size();

    _nplanehits_U = 0;
    _nplanehits_V = 0;
    _nplanehits_Y = 0;

    if (clus->Plane().Plane == 0)
    {
      _nplanehits_U = nhits;
    }
        else if (clus->Plane().Plane == 1)
    {
      _nplanehits_V = nhits;
    }
        else if (clus->Plane().Plane == 2)
    {
      _nplanehits_Y = nhits;
    }
        for (const auto &hit : clus_hit_v)
    {
      hit_v.push_back(hit);
    }
  } // for all clusters associated to PFP

  // store Backtracking
  if (!fData)
  {
    if (clus_pxy_v.size() != 0)
    {
      float purity = 0., completeness = 0., overlay_purity = 0.;
      _backtracked_pdg = 0;
      int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
      if (ibt >= 0)
      {
        auto &mcp = btparts_v[ibt];
        _backtracked_e = mcp.e;
        _backtracked_pdg = mcp.pdg;
        _backtracked_purity = purity;
        _backtracked_completeness = completeness;
        _backtracked_overlay_purity = overlay_purity;

        _backtracked_px = mcp.px;
        _backtracked_py = mcp.py;
        _backtracked_pz = mcp.pz;
        _backtracked_start_x = mcp.start_x;
        _backtracked_start_y = mcp.start_y;
        _backtracked_start_z = mcp.start_z;
        _backtracked_start_t = mcp.start_t;

        _backtracked_start_U = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0);
        _backtracked_start_V = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1);
        _backtracked_start_Y = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2);

        float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};

        if (mcp.pdg == 11 || mcp.pdg == 22)
        {
          reco_st[0] += searchingfornues::x_offset(mcp.start_t);
        }
        else
        {
          searchingfornues::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
        }
        _backtracked_sce_start_x = reco_st[0];
        _backtracked_sce_start_y = reco_st[1];
        _backtracked_sce_start_z = reco_st[2];

        _backtracked_sce_start_U = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0);
        _backtracked_sce_start_V = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1);
        _backtracked_sce_start_Y = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2);

        auto const &mcparticles_v = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
        for (size_t p = 0; p < mcparticles_v->size(); p++)
        {
          auto mcp = mcparticles_v->at(p);
          if ((mcp.StatusCode() == 1) && ((mcp.Momentum(0).E() - _backtracked_e) < 0.0001))
          {
            _backtracked_end_process = mcp.EndProcess();
            // stops in the detector?
            art::ServiceHandle<geo::Geometry> geo;
            geo::TPCGeo const &thisTPC = geo->TPC();
            geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
            if ((mcp.EndPosition().X() > theTpcGeo.MinX()) && (mcp.EndPosition().X() < theTpcGeo.MaxX()) &&
                (mcp.EndPosition().Y() > theTpcGeo.MinY()) && (mcp.EndPosition().Y() < theTpcGeo.MaxY()) &&
                (mcp.EndPosition().Z() > theTpcGeo.MinZ()) && (mcp.EndPosition().Z() < theTpcGeo.MaxZ()))
              _backtracked_end_in_tpc = true;
          }
        }
      }
    }
  }

  /*
  // get trk proxy in order to fetch PID
  auto trkpxy2 = pid_proxy[trk.key()];
  auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

  //collection plane
  _trk_bragg_p_y = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
        searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));
  _trk_bragg_mu_y = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
         searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));
  _trk_bragg_mip_y = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
  _trk_pid_chipr_y = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 2);
  _trk_pid_chimu_y = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 2);
  _trk_pid_chipi_y = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 2);
  _trk_pid_chika_y = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 2);
  _trk_pida_y = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

  //u plane
  _trk_bragg_p_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 0),
          searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 0));
  _trk_bragg_mu_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 0),
           searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 0));
  _trk_bragg_mip_u = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 0);
  _trk_pid_chipr_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
  _trk_pid_chimu_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
  _trk_pid_chipi_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
  _trk_pid_chika_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);
  _trk_pida_u = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 0);

  //v plane
  _trk_bragg_p_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 1),
          searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 1));
  _trk_bragg_mu_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 1),
           searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 1));
  _trk_bragg_mip_v = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 1);
  _trk_pid_chipr_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 1);
  _trk_pid_chimu_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 1);
  _trk_pid_chipi_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 1);
  _trk_pid_chika_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 1);
  _trk_pida_v = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 1);

  //three plane pid from Pip
  _trk_bragg_p_three_planes = searchingfornues::PID(pidpxy_v[0], "ThreePlaneProtonPID", anab::kLikelihood, anab::kForward, 2212, -1);
  */


  // Kinetic energy using tabulated stopping power (GeV)
  _trk_mcs_muon_mom = _mcsfitter.fitMcs(trk->Trajectory(), 13).bestMomentum();
  _trk_range_muon_mom = _trkmom.GetTrackMomentum(searchingfornues::GetSCECorrTrackLength(trk), 13);
  _trk_energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(searchingfornues::GetSCECorrTrackLength(trk), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
  _trk_energy_muon = std::sqrt(std::pow(_trk_mcs_muon_mom, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

  _trk_theta = trk->Theta();
  _trk_phi = trk->Phi();
  _trk_len = searchingfornues::GetSCECorrTrackLength(trk);

  TVector3 trk_vtx_v;
  trk_vtx_v.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
  trk_vtx_v -= nu_vtx;
  _trk_distance = trk_vtx_v.Mag();

  _trk_dir_x = trk->StartDirection().X();
  _trk_dir_y = trk->StartDirection().Y();
  _trk_dir_z = trk->StartDirection().Z();

  _trk_start_x = trk->Start().X();
  _trk_start_y = trk->Start().Y();
  _trk_start_z = trk->Start().Z();

  _trk_end_x = trk->End().X();
  _trk_end_y = trk->End().Y();
  _trk_end_z = trk->End().Z();

  float _trk_start_sce[3];
  searchingfornues::ApplySCECorrectionXYZ(_trk_start_x, _trk_start_y, _trk_start_z, _trk_start_sce);
  _trk_sce_start_x = _trk_start_sce[0];
  _trk_sce_start_y = _trk_start_sce[1];
  _trk_sce_start_z = _trk_start_sce[2];

  float _trk_end_sce[3];
  searchingfornues::ApplySCECorrectionXYZ(_trk_end_x, _trk_end_y, _trk_end_z, _trk_end_sce);
  _trk_sce_end_x = _trk_end_sce[0];
  _trk_sce_end_y = _trk_end_sce[1];
  _trk_sce_end_z = _trk_end_sce[2];

  // fill Calorimetry

  int key = trk.key();

  if (fGetCaloID)
  {
    int caloctr = 0;
    for (const searchingfornues::ProxyCaloElem_t tkcalo : calo_proxy)
    {
      // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
      if (tkcalo->ID() == int(pfp.index()))
      {
        key = caloctr;
        break;
      }
      caloctr += 1;
    }// for all calo-proxy objects
  }// if we want to get the calo object index manually

  auto calo_v = calo_proxy[key].get<anab::Calorimetry>();

  for (auto const& calo : calo_v)
  {
    auto const& plane = calo->PlaneID().Plane;
    if (plane>2)
    {
      continue;
    }
    auto const& xyz_v = calo->XYZ();

    // fill vector of boolean to determine if hit has to be corrected or not
    std::vector<bool> is_hit_montecarlo;

    if (!fData)
    {
      art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
      std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
      assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));

      const std::vector< size_t > &tp_indices = calo->TpIndices();
      for (size_t i = 0; i < tp_indices.size(); i++)
      {
        size_t tp_index = tp_indices[i];
        is_hit_montecarlo.push_back(searchingfornues::isHitBtMonteCarlo(tp_index, assocMCPart, fEnergyThresholdForMCHits));
      }
    }
    else
    {
      for (size_t i = 0; i < xyz_v.size(); i++)
      {
        is_hit_montecarlo.push_back(false);
      }
    }

    if (plane == 0)
    {
      _dqdx_u = calo->dQdx();
      _rr_u = calo->ResidualRange();
      _pitch_u = calo->TrkPitchVec();
      _is_hit_montecarlo_u = is_hit_montecarlo;
      for (auto xyz : xyz_v)
      {
        _x_u.push_back(xyz.X());
        _y_u.push_back(xyz.Y());
        _z_u.push_back(xyz.Z());

        float _dir_u[3];
        searchingfornues::TrkDirectionAtXYZ(trk.value(), xyz.X(), xyz.Y(), xyz.Z(), _dir_u);
        _dir_x_u.push_back(_dir_u[0]);
        _dir_y_u.push_back(_dir_u[1]);
        _dir_z_u.push_back(_dir_u[2]);
      }
      if (fShrFit)
      {
        _dedx_u = searchingfornues::GetdEdxfromdQdx(_dqdx_u, _x_u, _y_u, _z_u, 2.1, fADCtoE[plane]);
      }
      else
      {
        _dedx_u = calo->dEdx();
      }
    }
    else if (plane == 1)
    {
      _dqdx_v = calo->dQdx();
      _rr_v = calo->ResidualRange();
      _pitch_v = calo->TrkPitchVec();
      _is_hit_montecarlo_v = is_hit_montecarlo;
      for (auto xyz : xyz_v)
      {
        _x_v.push_back(xyz.X());
        _y_v.push_back(xyz.Y());
        _z_v.push_back(xyz.Z());

        float _dir_v[3];
        searchingfornues::TrkDirectionAtXYZ(trk.value(), xyz.X(), xyz.Y(), xyz.Z(), _dir_v);
        _dir_x_v.push_back(_dir_v[0]);
        _dir_y_v.push_back(_dir_v[1]);
        _dir_z_v.push_back(_dir_v[2]);
      }
      if (fShrFit)
      {
        _dedx_v = searchingfornues::GetdEdxfromdQdx(_dqdx_v, _x_v, _y_v, _z_v, 2.1, fADCtoE[plane]);
      }
      else
      {
        _dedx_v = calo->dEdx();
      }
    }
    else if (plane == 2) //collection
    {
      _dqdx_y = calo->dQdx();
      _rr_y = calo->ResidualRange();
      _pitch_y = calo->TrkPitchVec();
      _is_hit_montecarlo_y = is_hit_montecarlo;
      for (auto xyz : xyz_v)
      {
        _x_y.push_back(xyz.X());
        _y_y.push_back(xyz.Y());
        _z_y.push_back(xyz.Z());

        float _dir_y[3];
        searchingfornues::TrkDirectionAtXYZ(trk.value(), xyz.X(), xyz.Y(), xyz.Z(), _dir_y);
        _dir_x_y.push_back(_dir_y[0]);
        _dir_y_y.push_back(_dir_y[1]);
        _dir_z_y.push_back(_dir_y[2]);
      }
      if (fShrFit)
      {
        _dedx_y = searchingfornues::GetdEdxfromdQdx(_dqdx_y, _x_y, _y_y, _z_y, 2.1, fADCtoE[plane]);
      }
      else
      {
        _dedx_y = calo->dEdx();
      }
    }
  }

  _calo_tree->Fill();
}


DEFINE_ART_CLASS_TOOL(CalorimetryAnalysis)
} // namespace analysis

#endif
