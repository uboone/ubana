////////////////////////////////////////////////////////////////////////
// Class:       MCC9SliceIDFilter
// Plugin Type: filter (art v3_01_02)
// File:        MCC9SliceIDFilter_module.cc
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


class MCC9SliceIDFilter;


class MCC9SliceIDFilter : public art::EDFilter {
public:
  explicit MCC9SliceIDFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCC9SliceIDFilter(MCC9SliceIDFilter const&) = delete;
  MCC9SliceIDFilter(MCC9SliceIDFilter&&) = delete;
  MCC9SliceIDFilter& operator=(MCC9SliceIDFilter const&) = delete;
  MCC9SliceIDFilter& operator=(MCC9SliceIDFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // Declare member data here.

  bool endSubRun(art::SubRun &subrun) override;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

  // input module labels
  art::InputTag fPANDORAproducer;
};


MCC9SliceIDFilter::MCC9SliceIDFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  fPANDORAproducer = p.get<art::InputTag>("PANDORAproducer");

}

bool MCC9SliceIDFilter::filter(art::Event& e)
{

  // grab PFParticles in event
  ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPANDORAproducer,
											 proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPANDORAproducer),
											 proxy::withAssociated<recob::Cluster>(fPANDORAproducer),
											 proxy::withAssociated<recob::Slice>(fPANDORAproducer),
											 proxy::withAssociated<recob::Track>(fPANDORAproducer),
											 proxy::withAssociated<recob::Vertex>(fPANDORAproducer),
											 proxy::withAssociated<recob::PCAxis>(fPANDORAproducer),
											 proxy::withAssociated<recob::Shower>(fPANDORAproducer));

  
  BuildPFPMap(pfp_proxy);

  // select neutrino slice if it exists
  for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
    
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false)
      continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      return true;
      
    }// if the neutrino pfparticle
    
    // Implementation of required member function here.
  }// for all PFParticles in the proxy
  
  return false;
  
}// end of filter function


void MCC9SliceIDFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col) {
  
  _pfpmap.clear();
  
  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
    {
      _pfpmap[pfp_pxy->Self()] = p;
      p++;
    }
  
  return;
} // BuildPFPMap


bool MCC9SliceIDFilter::endSubRun(art::SubRun &subrun)
{

  return true;
}

void MCC9SliceIDFilter::beginJob()
{
  // Implementation of optional member function here.
}

void MCC9SliceIDFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MCC9SliceIDFilter)
