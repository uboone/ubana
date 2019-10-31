////////////////////////////////////////////////////////////////////////
// Class:       NuCCfilter
// Plugin Type: filter (art v3_01_02)
// File:        NuCCfilter_module.cc
//
// Generated at Wed Jul 31 15:59:55 2019 by Wouter Van de pontseele using cetskelgen
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include <memory>

class NuCCfilter;

class NuCCfilter : public art::EDFilter
{
public:
  explicit NuCCfilter(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuCCfilter(NuCCfilter const &) = delete;
  NuCCfilter(NuCCfilter &&) = delete;
  NuCCfilter &operator=(NuCCfilter const &) = delete;
  NuCCfilter &operator=(NuCCfilter &&) = delete;

  // Required functions.
  bool filter(art::Event &e) override;

private:
  std::string m_pfp_producer;
  std::string m_muon_producer;
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleMap particleMap;
};

NuCCfilter::NuCCfilter(fhicl::ParameterSet const &p)
    : EDFilter{p}
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandora");
  m_muon_producer = p.get<std::string>("muon_producer", "NuCCproducer");
}

bool NuCCfilter::filter(art::Event &evt)
{
  // Reset fields for next event
  pfparticles.clear();
  pfneutrinos.clear();
  particleMap.clear();
  larpandora.CollectPFParticles(evt, m_pfp_producer, pfparticles);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  bool selected = false;

  if (pfparticles.size() == 0)
  {
    std::cout << "[NuCCfilter::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  }
  else
  {
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1)
    {
      std::cout << "[NuCCfilter::FillReconstructed] Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
    }
    else // We have a reconstructed neutrino
    {
      art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

      art::Handle<std::vector<recob::PFParticle>> pfparticles_handle;
      evt.getByLabel(m_pfp_producer, pfparticles_handle);
      art::FindManyP<anab::T0> pfp_muon_assn(pfparticles_handle, evt, m_muon_producer);

      for (size_t daughter_id : pfnu->Daughters())
      {
        const std::vector<art::Ptr<anab::T0>> T0_muon = pfp_muon_assn.at(particleMap.at(daughter_id).key());
        if (T0_muon.size() != 0)
        {
          std::cout << "[NuCCfilter::filter] Muon neutrino daughter found! Event passed filter." << std::endl;
          selected = true;
        }
      }
    }
  }
  return selected;
}
DEFINE_ART_MODULE(NuCCfilter)