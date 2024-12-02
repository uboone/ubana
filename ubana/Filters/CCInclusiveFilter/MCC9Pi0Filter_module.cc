////////////////////////////////////////////////////////////////////////
// Class:       MCC9Pi0Filter
// Plugin Type: filter (art v3_01_02)
// File:        MCC9Pi0Filter_module.cc
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

#include "art_root_io/TFileService.h"
#include "TTree.h"

using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle>>(std::declval<art::Event>(), std::declval<art::InputTag>(),
										     proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Track>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::PCAxis>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Shower>(std::declval<art::InputTag>())));
using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;


class MCC9Pi0Filter;


class MCC9Pi0Filter : public art::EDFilter {
public:
  explicit MCC9Pi0Filter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCC9Pi0Filter(MCC9Pi0Filter const&) = delete;
  MCC9Pi0Filter(MCC9Pi0Filter&&) = delete;
  MCC9Pi0Filter& operator=(MCC9Pi0Filter const&) = delete;
  MCC9Pi0Filter& operator=(MCC9Pi0Filter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_pxy : parent pfparticle proxy for which to add daughters
   * @input pfp_pxy_col : evnt PFP proxy collection
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                    const ProxyPfpColl_t &pfp_pxy_col,
                    std::vector<ProxyPfpElem_t> &slice_v);

  std::pair<double,double> VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir);

  float GetTrackShowerScore(const ProxyPfpElem_t &pfp_pxy);

  // select pi0 event
  bool selectEvent(const std::vector<ProxyPfpElem_t>& pfp_pxy_v);

  void SaveTruth(art::Event const &e);

  bool endSubRun(art::SubRun &subrun) override;

  void ResetTTree();

  // input module labels
  art::InputTag fPANDORAproducer;
  art::InputTag fSHRproducer;
  art::InputTag fMCTproducer;
  art::InputTag fMCPproducer;

  // variables with which to apply cuts
  float fDotMin, fTrkShrScore, fE1Min, fE2Min, fGammaDotMax, fRMin;
  bool fData;

  // TTree variables
  TTree* _tree;
  // event information
  int _run, _sub, _evt;
  // reco variables
  float _e1, _e2; // energy of leading / subleading photon
  float _r1, _r2; // conversion distances
  float _dedx1, _dedx2; // dedx
  float _gammadot; // opening angle between photons
  float _d1, _d2;   // dot product between shower dir and vtx -> shr vector.
  float _shrscore1, _shrscore2; // shower scores
  float _mass; // pi0 reco mass
  int _nslice;
  // truth variables
  int _nu_pdg;                  /**< neutrino PDG code */
  int _ccnc;                    /**< CC or NC tag from GENIE */
  int _interaction;             /**< Interaction code from GENIE */
  float _vtx_x, _vtx_y, _vtx_z; /**< neutrino interaction vertex coordinates [cm] */
  float _vtx_t;                 /**< neutrino generation time */
  float _nu_e;
  float _lep_e;
  // final state particle information
  int _nmuon;                            /**< is there a final-state muon from the neutrino? [1=yes 0=no] */
  float _muon_e;                         /**< energy, purity, completeness. */
  int _nelec;                            /**< is there a final-state electron from the neutrino? [1=yes 0=no] */
  float _elec_e;                         /**< energy, purity, completeness. */
  int _npi0;                             /**< how many pi0s are there? */
  //int _pi0;                              /**< is there a final-state pi0 from the neutrino? [1=yes 0=no] */
  float _pi0_e;                          /**< energy, purity, completeness. */
  int _nneutron;                         /**< how many neutrons are there? */
  int _nproton;                          /**< how many protons are there? */
  //int _proton;                           /**< is there a final-state proton from the neutrino? [1=yes 0=no] */
  float _proton_e;                       /**< energy, purity, completeness. */
  int _npion;                            /**< how many pions are there? */
  //int _pion;                             /**< is there a final-state charged pion from the neutrino? [1=yes 0=no] */
  float _pion_e;                         /**< energy, purity, completeness. */

  TTree *_subrun_tree;
  int _run_sr; // The run number
  int _sub_sr; // The subRun number
  float _pot;  // The total amount of POT for the current sub run



};


MCC9Pi0Filter::MCC9Pi0Filter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  fPANDORAproducer = p.get<art::InputTag>("PANDORAproducer");
  fSHRproducer     = p.get<art::InputTag>("SHRproducer");
  fMCTproducer     = p.get<art::InputTag>("MCTproducer");
  fMCPproducer     = p.get<art::InputTag>("MCPproducer");
  fData            = p.get<bool>("IsData");
  fRMin            = p.get<float>("RMin");
  fDotMin          = p.get<float>("DotMin");
  fTrkShrScore     = p.get<float>("TrkShrScore");
  fE1Min           = p.get<float>("E1Min");
  fE2Min           = p.get<float>("E2Min");
  fGammaDotMax     = p.get<float>("GammaDotMax");

  // setup TTree
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("MCC9Pi0Filter", "Pi0 Filter TTree");
  // event information
  _tree->Branch("evt",&_evt,"evt/I");
  _tree->Branch("sub",&_sub,"sub/I");
  _tree->Branch("run",&_run,"run/I");
  // reco information
  _tree->Branch("nslice",&_nslice,"nslice/I");
  _tree->Branch("e1",&_e1,"e1/F");
  _tree->Branch("e2",&_e2,"e2/F");
  _tree->Branch("r1",&_r1,"r1/F");
  _tree->Branch("r2",&_r2,"r2/F");
  _tree->Branch("dedx1",&_dedx1,"dedx1/F");
  _tree->Branch("dedx2",&_dedx2,"dedx2/F");
  _tree->Branch("shrscore1",&_shrscore1,"shrscore1/F");
  _tree->Branch("shrscore2",&_shrscore2,"shrscore2/F");
  _tree->Branch("d1",&_d1,"d1/F");
  _tree->Branch("d2",&_d2,"d2/F");
  _tree->Branch("gammadot",&_gammadot,"gammadot/F");
  _tree->Branch("mass",&_mass,"mass/F");
  // truth information
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("ccnc", &_ccnc, "ccnc/I");
  _tree->Branch("interaction", &_interaction, "interaction/I");
  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("lep_e", &_lep_e, "lep_e/F");
  // neutrino vertex
  _tree->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _tree->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _tree->Branch("vtx_z", &_vtx_z, "vtx_z/F");
  // muon
  _tree->Branch("nmuon", &_nmuon, "nmuon/I");
  _tree->Branch("muon_e", &_muon_e, "muon_e/F");
  // electron
  _tree->Branch("nelec", &_nelec, "nelec/I");
  _tree->Branch("elec_e", &_elec_e, "elec_e/F");
  // pi0
  _tree->Branch("npi0", &_npi0, "npi0/I");
  _tree->Branch("pi0_e", &_pi0_e, "pi0_e/F");
  // neutrons
  _tree->Branch("nneutron", &_nneutron, "nneutron/I");
  // first [highest momentum] proton
  _tree->Branch("nproton", &_nproton, "nproton/I");
  _tree->Branch("proton_e", &_proton_e, "proton_e/F");
  // charged pions
  _tree->Branch("npion", &_npion, "npion/I");
  _tree->Branch("pion_e", &_pion_e, "pion_e/F");

  _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
  _subrun_tree->Branch("run", &_run_sr, "run/I");
  _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

  if (!fData)
    _subrun_tree->Branch("pot", &_pot, "pot/F");

}

bool MCC9Pi0Filter::filter(art::Event& e)
{

  ResetTTree();
  
  // grab PFParticles in event
  ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPANDORAproducer,
											 proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPANDORAproducer),
											 proxy::withAssociated<recob::Cluster>(fPANDORAproducer),
											 proxy::withAssociated<recob::Slice>(fPANDORAproducer),
											 proxy::withAssociated<recob::Track>(fPANDORAproducer),
											 proxy::withAssociated<recob::Vertex>(fPANDORAproducer),
											 proxy::withAssociated<recob::PCAxis>(fPANDORAproducer),
											 proxy::withAssociated<recob::Shower>(fSHRproducer));

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();
  
  BuildPFPMap(pfp_proxy);

  if (!fData) 
    SaveTruth(e); 

  bool selectedpi0 = false;
  
  // select neutrino slice if it exists
  for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
    
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false)
      continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      _nslice += 1;

      // collect PFParticle hierarchy originating from this neutrino candidate
      std::vector<ProxyPfpElem_t> slice_pfp_v;
      AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
      
      // apply pi0 selection
      selectedpi0 = selectEvent(slice_pfp_v);

    }// if the neutrino pfparticle

  // Implementation of required member function here.
  }// for all PFParticles in the proxy

  _tree->Fill();

  if (selectedpi0) return true;

  return false;
  
}// end of filter function

void MCC9Pi0Filter::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
			     const ProxyPfpColl_t &pfp_pxy_col,
			     std::vector<ProxyPfpElem_t> &slice_v)
{
  
  auto daughters = pfp_pxy->Daughters();
  
  slice_v.push_back(pfp_pxy);
  
  for (auto const &daughterid : daughters)
    {
      
      if (_pfpmap.find(daughterid) == _pfpmap.end())
	{
	  std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
	  continue;
	}
      
      // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
      auto pfp_pxy2 = pfp_pxy_col.begin();
      for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
	++pfp_pxy2;
      // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));
      
      AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
      
    } // for all daughters
  
  return;
} // AddDaughters

void MCC9Pi0Filter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col) {
  
  _pfpmap.clear();
  
  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
    {
      _pfpmap[pfp_pxy->Self()] = p;
      p++;
    }
  
  return;
} // BuildPFPMap


bool MCC9Pi0Filter::selectEvent(const std::vector<ProxyPfpElem_t>& pfp_pxy_v) {
  
  TVector3 nuvtx;
  Double_t xyz[3] = {};
  
  TVector3 gammadir1, gammadir2;
  
  // vector of pfp indices for showers
  std::vector<float> shr_energy_v;
  // vector of shower energies
  std::vector<short> pfp_idx_v;
  
  short pfp_ctr = 0;
  
  // loop over all PFPs and get metadata, find the neutrino and save its vertex
  for (const auto& pfp_pxy : pfp_pxy_v) {
    
    pfp_ctr += 1;
    
    auto PDG = fabs(pfp_pxy->PdgCode());
    
    if ( (PDG == 12) || (PDG == 14) ) {
      
      // grab vertex
      auto vtx = pfp_pxy.get<recob::Vertex>();
      if (vtx.size() != 1) {
	return false;
      }
      
      // save vertex to array
      vtx.at(0)->XYZ(xyz);
      nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);
      
    }// if neutrino PFP
    
    // grab shower/track score
    auto trkshrscore = GetTrackShowerScore(pfp_pxy);
    
    auto nshr = pfp_pxy.get<recob::Shower>().size();
    //auto ntrk = pfp_pxy.get<recob::Track>().size();
    
    // 1 -> track-like
    if (trkshrscore > fTrkShrScore)  continue;
    
    if (nshr != 1) continue;
    
    auto const& shr = pfp_pxy.get<recob::Shower>().at(0);
    
    auto energy = shr->Energy()[2];
    
    auto vtxcompat = VtxCompatibility(nuvtx, shr->ShowerStart(), shr->Direction());
    
    // if blank result, continue
    if ( (vtxcompat.first == -1) && (vtxcompat.second == -1) ) continue;
    // if too close to vertex or too mis-algined, continue
    if ( (vtxcompat.second < fRMin) || (vtxcompat.first < fDotMin) ) continue;
    
    shr_energy_v.push_back( energy );
    pfp_idx_v.push_back( pfp_ctr - 1 );
    
  }// for all PFP
  
  // if we did not find any shower, return
  if (shr_energy_v.size() == 0) return false;
  
  // if a single shower, backtrack and quit
  if (shr_energy_v.size() == 1) return false;
  
  // if two or more, sort by energy
  double e1 = 0; // energy of highest energy reco shower
  short  i1 = 0; // pfp index for highest energy reco shower
  double e2 = 0; // energy of 2nd highest energy reco shower
  short  i2 = 0; // pfp index for 2nd highest energy reco shower
  
  // find highest energy shower
  for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
    if (shr_energy_v[n0] > e1) {
      e1 = shr_energy_v[n0];
      i1 = pfp_idx_v[n0];
    } 
  }
  // find second highest energy shower
  for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
    if (pfp_idx_v[n0] == i1) continue;
    if (shr_energy_v[n0] > e2) {
      e2 = shr_energy_v[n0];
      i2 = pfp_idx_v[n0];
    }
  }
  
  auto shr1 = pfp_pxy_v.at(i1).get<recob::Shower>().at(0);
  auto shr2 = pfp_pxy_v.at(i2).get<recob::Shower>().at(0);
  
  _e1 = shr1->Energy()[2];
  _e2 = shr2->Energy()[2];

  _shrscore1 = GetTrackShowerScore(pfp_pxy_v.at(i1));
  _shrscore2 = GetTrackShowerScore(pfp_pxy_v.at(i2));
  
  auto vtxcompat1 = VtxCompatibility(nuvtx, shr1->ShowerStart(), shr1->Direction());
  
  _r1    = vtxcompat1.second;
  _d1    = vtxcompat1.first;
  _dedx1 = shr1->dEdx()[2];
  
  auto vtxcompat2 = VtxCompatibility(nuvtx, shr2->ShowerStart(), shr2->Direction());
  
  _r2    = vtxcompat2.second;
  _d2    = vtxcompat2.first;
  _dedx2 = shr2->dEdx()[2];

  _gammadot = shr1->Direction().Dot(shr2->Direction());

  _mass = sqrt( 2 * _e1 * _e2 * (1 - _gammadot ) );

  if (_r1 < fRMin) return false;
  if (_r2 < fRMin) return false;
  if (_e1 < fE1Min) return false;
  if (_e2 < fE2Min) return false;
  if (_gammadot > fGammaDotMax) return false;
  
  return true;
}

std::pair<double,double> MCC9Pi0Filter::VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir) {
  
  // grab shower start point and direction
  //auto shrvtx = shr->ShowerStart();
  //auto shrdir = shr->Direction();
  //gammadir = shrdir;
  
  // assess compatibility
  auto nuvtx2shrvtx = (shrvtx - nuvtx).Unit();
  auto shrdirnormed = shrdir.Unit();
  
  double dot  = nuvtx2shrvtx.Dot(shrdirnormed);
  double dist = (nuvtx-shrvtx).Mag();
  
  return std::make_pair(dot,dist);
  
}// end of vertex compatibility

float MCC9Pi0Filter::GetTrackShowerScore(const ProxyPfpElem_t &pfp_pxy)
{
  
  const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();
  
  if (pfParticleMetadataList.size() == 0)
    return 1;
  
  for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
    {

    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
    if (!pfParticlePropertiesMap.empty())
    {
      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
      {
        if (it->first == "TrackScore")
          return it->second;
      } // for map elements
    }   // if pfp metadata map not empty
  }     // for list

  return 1;
}

void MCC9Pi0Filter::SaveTruth(art::Event const &e)
{

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu = neutrino.Nu();
  auto lep = neutrino.Lepton();

  _ccnc = neutrino.CCNC();
  _interaction = neutrino.Mode();
  _nu_pdg = nu.PdgCode();
  _nu_e = -1;
  if ( (nu.Trajectory().size() == 0) || (lep.Trajectory().size() == 0) )
    return;

  _nu_e = nu.Trajectory().E(0);
  _lep_e = lep.Trajectory().E(0);
  _vtx_x = nu.EndX();
  _vtx_y = nu.EndY();
  _vtx_z = nu.EndZ();
  _vtx_t = nu.T();

  _nelec = 0;
  _nmuon = 0;
  _npi0 = 0;
  _nproton = 0;
  _npion = 0;
  _nneutron = 0;


  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++)
  {

    auto const &part = mct.GetParticle(i);
    if (part.StatusCode() != 1)
    {
      continue;
    }

    // if muon
    if ((std::abs(part.PdgCode()) == 11) and (part.StatusCode() == 1))
    {
      _nmuon += 1;
      
      if (part.Momentum(0).E() > _muon_e)
        _muon_e = part.Momentum(0).E();
    } // if muon

    // if electron
    else if ((std::abs(part.PdgCode()) == 13) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() > 0.020)
      {
        _nelec += 1;
      }
      if (part.Momentum(0).E() > _elec_e)
        _elec_e = part.Momentum(0).E();
    } // if electron

    // if pi0
    else if ((part.PdgCode() == 111) and (part.StatusCode() == 1))
    {
      _npi0 += 1;
      if (part.Momentum(0).E() > _pi0_e)
        _pi0_e = part.Momentum(0).E();
    } // if pi0

    // if proton
    else if ((part.PdgCode() == 2212) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() - 0.938 > 0.040)
      {
        _nproton += 1;
      }
      if (part.Momentum(0).E() > _proton_e)
        _proton_e = part.Momentum(0).E();

    } // if proton

    // if neutron
    else if ((part.PdgCode() == 2112) and (part.StatusCode() == 1))
    {
      _nneutron += 1;
    }

    // if pion
    else if ((std::abs(part.PdgCode()) == 211) and (part.StatusCode() == 1))
    {
      _npion += 1;
      if (part.Momentum(0).E() > _pion_e)
        _pion_e = part.Momentum(0).E();
    } // if pion

  } // for all MCParticles

  return;
}

bool MCC9Pi0Filter::endSubRun(art::SubRun &subrun)
{
  if (!fData)
  {
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    _pot = subrun.getByLabel(fMCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    std::cout << "[LArPandoraExternalEventBuilding::endSubRun] Storing POT info!" << std::endl;
  }

  _run_sr = subrun.run();
  _sub_sr = subrun.subRun();
  _subrun_tree->Fill();
  return true;
}

void MCC9Pi0Filter::ResetTTree() {

  _nslice = 0;
  _e1 = _e2 = 0;
  _r1 = _r2 = 0;
  _dedx1 = _dedx2 = 0;
  _shrscore1 = _shrscore2 = 0;
  _d1 = _d2 = 0;
  _gammadot = 0;
  _mass = 0;
  _nu_pdg = -1;
  _ccnc = -1;
  _interaction = -1;
  _nu_e = -1;
  _lep_e = -1;
  _vtx_x = _vtx_y = _vtx_z = -1e3;
  

}

void MCC9Pi0Filter::beginJob()
{
  // Implementation of optional member function here.
}

void MCC9Pi0Filter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MCC9Pi0Filter)
