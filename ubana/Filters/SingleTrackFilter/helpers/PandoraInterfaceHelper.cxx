#ifndef PANDORAINTERFACEHELPER_CXX
#define PANDORAINTERFACEHELPER_CXX

#include "PandoraInterfaceHelper.h"

PandoraInterfaceHelper::PandoraInterfaceHelper()
{
  m_configured = false;
}

void PandoraInterfaceHelper::Configure(art::Event const &e,
                                       std::string m_pfp_producer,
                                       std::string m_spacepoint_producer,
                                       std::string m_hitfinder_producer,
                                       std::string m_geant_producer,
                                       std::string m_hit_mcp_producer)
{

  lar_pandora::LArPandoraHelper::DaughterMode daughterMode = lar_pandora::LArPandoraHelper::kAddDaughters;

  // Collect hits

  lar_pandora::HitVector hitVector;
  //lar_pandora::LArPandoraHelper::CollectHits(e, m_hitfinder_producer, hitVector);
  art::Handle<std::vector<recob::Hit>> hit_h;

  e.getByLabel(m_hitfinder_producer, hit_h);
  if (!hit_h.isValid())
  {
    std::cout << "[McPfpMatch] Hit Handle is not valid." << std::endl;
    // throw std::exception();
  }

  art::fill_ptr_vector(hitVector, hit_h);

  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcps_from_hit(hit_h, e, m_hit_mcp_producer);

  // Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector recoParticleVector;
  lar_pandora::PFParticleVector recoNeutrinoVector;
  lar_pandora::PFParticlesToHits pfp_to_hits_map;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, m_pfp_producer, recoParticleVector);

  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);

  try
  {
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e,
                                                          m_pfp_producer,
                                                          m_spacepoint_producer,
                                                          pfp_to_hits_map,
                                                          recoHitsToParticles,
                                                          daughterMode,
                                                          true); // Use clusters to go from pfp to hits
  }
  catch (...)
  {
    std::cout << "[McPfpMatch] BuildPFParticleHitMaps error" << std::endl;
  }

  m_verbose = false;
  if (m_verbose)
  {
    std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
  }

  // Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits trueParticlesToHits;
  lar_pandora::HitsToMCParticles hit_to_mcps_map;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, trueParticleVector);
  lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);

  // Construct a Particle Map (trackID to MCParticle)
  lar_pandora::MCParticleMap particleMap;

  for (lar_pandora::MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end(); iter1 != iterEnd1; ++iter1)
  {
    const lar_pandora::MCParticleVector &particleVector = iter1->second;
    for (lar_pandora::MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
    {
      const art::Ptr<simb::MCParticle> particle = *iter2;
      particleMap[particle->TrackId()] = particle;
    }
  }

  // Loop over the hits, get the ass MCP, and then tru to link

  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  std::vector<anab::BackTrackerHitMatchingData const *> match_v;
  for (auto hit : hitVector)
  {

    mcp_v.clear();
    match_v.clear();
    mcps_from_hit.get(hit.key(), mcp_v, match_v);

    double max_energy = -1;
    int best_match_id = -1;
    for (size_t m = 0; m < match_v.size(); m++)
    {
      double this_energy = match_v[m]->energy;
      if (this_energy > max_energy)
      {
        best_match_id = m;
        max_energy = this_energy;
      }
    }

    if (best_match_id > -1)
    {
      try
      {
        const art::Ptr<simb::MCParticle> thisParticle = mcp_v.at(best_match_id);
        const art::Ptr<simb::MCParticle> primaryParticle(lar_pandora::LArPandoraHelper::GetFinalStateMCParticle(particleMap, thisParticle));
        const art::Ptr<simb::MCParticle> selectedParticle((lar_pandora::LArPandoraHelper::kAddDaughters == daughterMode) ? primaryParticle : thisParticle);

        if ((lar_pandora::LArPandoraHelper::kIgnoreDaughters == daughterMode) && (selectedParticle != primaryParticle))
          continue;

        if (!(lar_pandora::LArPandoraHelper::IsVisible(selectedParticle)))
          continue;

        hit_to_mcps_map[hit] = selectedParticle;
      }
      catch (...)
      {
        std::cout << "[PandoraInterfaceHelper] "
                  << "Error in the loop of the hits" << std::endl;
        continue;
      }
    }
  }

  // Now set the things we need for the future
  m_hit_to_mcps_map = hit_to_mcps_map;
  m_pfp_to_hits_map = pfp_to_hits_map;

  // std::cout << "hit_to_mcps_map size " << hit_to_mcps_map.size() << std::endl;

  m_configured = true;
}

art::Ptr<simb::MCTruth> PandoraInterfaceHelper::TrackIDToMCTruth(art::Event const &e, std::string m_geant_producer, int geant_track_id)
{

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);

  for (auto iter : particlesToTruth)
  {
    if (iter.first->TrackId() == geant_track_id)
    {
      return iter.second;
    }
  }
  art::Ptr<simb::MCTruth> null_ptr;
  return null_ptr;
}

void PandoraInterfaceHelper::GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles)
{
  bool m_debug = false;

  if (!m_configured)
  {
    std::cout << "Call to " << __PRETTY_FUNCTION__ << " whitout having done configuration. Abort." << std::endl;
    // throw std::exception();
  }

  // Loop over the reco particles
  for (auto iter1 : m_pfp_to_hits_map)
  {

    // The PFParticle
    const art::Ptr<recob::PFParticle> recoParticle = iter1.first;

    if (m_debug)
      std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] Looking at PFP with ID " << recoParticle->Self() << std::endl;

    // The PFParticle's hits
    const lar_pandora::HitVector &hitVector = iter1.second;

    if (m_debug)
      std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t This PFP has " << hitVector.size() << " hits." << std::endl;

    lar_pandora::MCParticlesToHits truthContributionMap;
    // Loop over all the hits associated to this reco particle
    for (auto hit : hitVector)
    {

      // Find the MCParticle that share this same hit (if any)
      auto iter3 = m_hit_to_mcps_map.find(hit);
      if (m_hit_to_mcps_map.end() == iter3)
        continue;

      // If exists, get the MCParticle
      const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

      if (m_debug)
        std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t Found a hit shared with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

      // This map will contain all the true particles that match some or all of the hits of the reco particle
      truthContributionMap[trueParticle].push_back(hit);
    }

    // Now we want to find the true particle that has more hits in common with this reco particle than the others
    lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

    for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
         iter4 != iterEnd4; ++iter4)
    {
      if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
      {
        mIter = iter4;
      }
    }

    if (truthContributionMap.end() != mIter)
    {
      const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

      if (m_debug)
        std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t >>> Match found with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

      // Emplace into the output map
      matchedParticles[recoParticle] = trueParticle;
    }

  } // m_pfp_to_hits_map loop ends
}

void PandoraInterfaceHelper::CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                                          const art::Ptr<recob::PFParticle> &particle,
                                                          lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
    downstreamPFParticles.push_back(particle);

  for (const auto &daughterId : particle->Daughters())
  {
    const auto iter(pfParticleMap.find(daughterId));
    if (iter == pfParticleMap.end())
      throw cet::exception("PandoraInterfaceHelper::CollectDownstreamPFParticles") << "Scrambled PFParticle IDs" << std::endl;

    this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
  }
}

float PandoraInterfaceHelper::Distance3D(float x1, float y1, float z1, float x2, float y2, float z2)
{
  return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

void PandoraInterfaceHelper::SCE(const float &x,
                                 const float &y,
                                 const float &z,
                                 const float &nu_time,
                                 float &x_out, float &y_out, float &z_out)
{
  auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
  auto sce_start = SCE.GetPosOffsets(geo::Point_t(x, y, z));

  y_out = y + sce_start.Y();
  z_out = z + sce_start.Z();

  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  float g4Ticks = detClocks->TPCG4Time2Tick(nu_time) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  float xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
  
  x_out = (x + xtimeoffset - sce_start.X())+ 0.6;
}

#endif