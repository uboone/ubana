////////////////////////////////////////////////////////////////////////
// Class:       Pi0Filter
// Plugin Type: filter (art v3_01_02)
// File:        Pi0Filter_module.cc
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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle>>(
										     std::declval<art::Event>(), std::declval<art::InputTag>(),
										     proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Track>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::PCAxis>(std::declval<art::InputTag>()),
										     proxy::withAssociated<recob::Shower>(std::declval<art::InputTag>())));
using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;


class Pi0Filter;


class Pi0Filter : public art::EDFilter {
public:
  explicit Pi0Filter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Filter(Pi0Filter const&) = delete;
  Pi0Filter(Pi0Filter&&) = delete;
  Pi0Filter& operator=(Pi0Filter const&) = delete;
  Pi0Filter& operator=(Pi0Filter&&) = delete;

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


  // input module labels
  art::InputTag fPANDORAproducer;
  art::InputTag fSHRproducer;

  // variables with which to apply cuts
  float _dmin, _dotmin, _trkshrscore, _e1min, _e2min, _gammadotmax;
 

};


Pi0Filter::Pi0Filter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  fPANDORAproducer = p.get<art::InputTag>("PANDORAproducer");
  fSHRproducer     = p.get<art::InputTag>("SHRproducer");

}

bool Pi0Filter::filter(art::Event& e)
{

  // grab PFParticles in event
  ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPANDORAproducer,
											 proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPANDORAproducer),
											 proxy::withAssociated<recob::Cluster>(fPANDORAproducer),
											 proxy::withAssociated<recob::Slice>(fPANDORAproducer),
											 proxy::withAssociated<recob::Track>(fPANDORAproducer),
											 proxy::withAssociated<recob::Vertex>(fPANDORAproducer),
											 proxy::withAssociated<recob::PCAxis>(fPANDORAproducer),
											 proxy::withAssociated<recob::Shower>(fSHRproducer));
  
  BuildPFPMap(pfp_proxy);
  
  
  // select neutrino slice if it exists
  for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
    
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false)
      continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      // collect PFParticle hierarchy originating from this neutrino candidate
      std::vector<ProxyPfpElem_t> slice_pfp_v;
      AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
      
      // apply pi0 selection
      bool selectedpi0 = selectEvent(slice_pfp_v);

      if (selectedpi0) return true;

    }// if the neutrino pfparticle

  // Implementation of required member function here.
  }// for all PFParticles in the proxy

  return false;
  
}// end of filter function

void Pi0Filter::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
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

void Pi0Filter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col) {
  
  _pfpmap.clear();
  
  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
    {
      _pfpmap[pfp_pxy->Self()] = p;
      p++;
    }
  
  return;
} // BuildPFPMap


bool Pi0Filter::selectEvent(const std::vector<ProxyPfpElem_t>& pfp_pxy_v) {
  
  TVector3 nuvtx;
  Double_t xyz[3] = {};
  
  TVector3 gammadir1, gammadir2;
  
  // vector of pfp indices for showers
  std::vector<float> shr_energy_v;
  // vector of shower energies
  std::vector<short> pfp_idx_v;
  
  short pfp_ctr = 0;
  
  // tagged shower reco variables
  float _energy1_Y, _energy2_Y; // _dot1, _dot2, _radlen1, _radlen2, _energy1_Y, _energy2_Y, _dedx1_Y, _dedx2_Y;
  // pi0 reco variables
  float _gammadot;
  
  // loop over all PFPs and get metadata, find the neutrino and save its vertex
  for (const auto& pfp_pxy : pfp_pxy_v) {
    
    pfp_ctr += 1;
    
    auto PDG = fabs(pfp_pxy->PdgCode());
    
    if ( (PDG == 12) || (PDG == 14) ) {
      
      // grab vertex
      auto vtx = pfp_pxy.get<recob::Vertex>();
      if (vtx.size() != 1) {
	std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
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
    if (trkshrscore > _trkshrscore)  continue;
    
    if (nshr != 1) continue;
    
    auto const& shr = pfp_pxy.get<recob::Shower>().at(0);
    
    auto energy = shr->Energy()[2];
    
    auto vtxcompat = VtxCompatibility(nuvtx, shr->ShowerStart(), shr->Direction());
    
    // if blank result, continue
    if ( (vtxcompat.first == -1) && (vtxcompat.second == -1) ) continue;
    // if too close to vertex or too mis-algined, continue
    if ( (vtxcompat.second < _dmin) || (vtxcompat.first < _dotmin) ) continue;
    
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
  
  //_shrscore1 = GetTrackShowerScore(pfp_pxy_v.at(i1));
  //_shrscore2 = GetTrackShowerScore(pfp_pxy_v.at(i2));
  
  //auto vtxcompat1 = VtxCompatibility(nuvtx, shr1->ShowerStart(), shr1->Direction());
  
  //_radlen1    = vtxcompat1.second;
  //_dot1       = vtxcompat1.first;
  _energy1_Y  = shr1->Energy()[2];
  //_dedx1_Y    = shr1->dEdx()[2];
  
  //auto vtxcompat2 = VtxCompatibility(nuvtx, shr2->ShowerStart(), shr2->Direction());
  
  //_radlen2   = vtxcompat2.second;
  //_dot2      = vtxcompat2.first;
  _energy2_Y = shr2->Energy()[2];
  //_dedx2_Y   = shr2->dEdx()[2];
  
  _gammadot = shr1->Direction().Dot(shr2->Direction());

  if (_energy1_Y < _e1min) return false;
  if (_energy2_Y < _e2min) return false;
  if (_gammadot > _gammadotmax) return false;
  
  return true;
}

std::pair<double,double> Pi0Filter::VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir) {
  
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

float Pi0Filter::GetTrackShowerScore(const ProxyPfpElem_t &pfp_pxy)
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

void Pi0Filter::beginJob()
{
  // Implementation of optional member function here.
}

void Pi0Filter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Pi0Filter)
