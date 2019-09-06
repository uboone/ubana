#include "NuCCproducer.h"

void NuCCproducer::produce(art::Event &evt)
{
  clearEvent();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  std::cout << "[NuCCproducer::produce]: Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent << std::endl;

  std::unique_ptr<std::vector<anab::T0>> T0_v(new std::vector<anab::T0>);
  std::unique_ptr<art::Assns<recob::PFParticle, anab::T0>> pfp_t0_assn_v(new art::Assns<recob::PFParticle, anab::T0>);

  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  if (pfparticles.size() == 0)
    std::cout << "[NuCCproducer::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  else
  {
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1)
      std::cout << "[NuCCproducer::FillReconstructed] Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
    else // We have a reconstructed neutrino
    {
      FillReconstructed(evt);
      // After all the fields are filled, do the selection and create association.
      fIsNuMuCC = IsNuMuCC(evt, T0_v, pfp_t0_assn_v);
    }
  }
  evt.put(std::move(T0_v));
  evt.put(std::move(pfp_t0_assn_v));
}

void NuCCproducer::FillReconstructed(art::Event const &evt)
{
  // Load associations and collections
  lar_pandora::VertexVector vertexVector_dummy;
  lar_pandora::PFParticleVector particleVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);
  larpandora.CollectTracks(evt, m_pfp_producer, pftracks, particlesToTracks);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track>>(m_pfp_producer);
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, evt, "pandoracalipidSCE");
  art::Handle<std::vector<recob::PFParticle>> pfparticles_handle;
  evt.getByLabel(m_pfp_producer, pfparticles_handle);
  art::FindManyP<anab::T0> nuFlashScoreAsso(pfparticles_handle, evt, "flashmatch");

  if (!trackPIDAssn.isValid())
  {
    std::cout << "[NuCCproducer::FillReconstructed] trackPIDAssn.isValid() == false" << std::endl;
  }

  // Start filling information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  fNu_PDG = pfnu->PdgCode();

  const std::vector<art::Ptr<anab::T0>> T0_flashchi_v = nuFlashScoreAsso.at(pfnu.key());
  if (T0_flashchi_v.size() == 1)
  {
    fNu_FlashChi2 = T0_flashchi_v.at(0)->TriggerConfidence();
  }

  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1)
  {
    std::cout << "[NuCCproducer::FillReconstructed] Neutrino association problem." << std::endl;
  }
  else
  {
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    fNu_Score = neutrino_properties.at("NuScore");
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    fNu_Vx = neutrino_vtx.X();
    fNu_Vy = neutrino_vtx.Y();
    fNu_Vz = neutrino_vtx.Z();
    std::vector<float> fid_vtx_v = {m_vtx_fid_x_start, m_vtx_fid_y_start, m_vtx_fid_z_start,
                                    m_vtx_fid_x_end, m_vtx_fid_y_end, m_vtx_fid_z_end};
    fNu_Contained = IsContained(fNu_Vx, fNu_Vy, fNu_Vz, fid_vtx_v);
  }
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);

  for (auto const pfp : pfdaughters)
  {
    if (!pfp->IsPrimary())
    {
      FillDaughters(pfp, trackPIDAssn);
    }
  }

  // Store the obvious cosmic with the lowest score:
  fBestObviousCosmic_FlashChi2 = std::numeric_limits<float>::max();
  for (auto const pfp : pfparticles)
  {
    // Only look at obvious cosmics:
    lar_pandora::MetadataVector pfp_metadata_vec = particlesToMetadata.at(pfp);
    const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = pfp_metadata_vec.front()->GetPropertiesMap();

    if (pfp_properties.count("IsClearCosmic"))
    {
      if (pfp_properties.at("IsClearCosmic") && pfp->IsPrimary())
      {
        const std::vector<art::Ptr<anab::T0>> T0_v = nuFlashScoreAsso.at(pfp.key());
        if (T0_v.size() == 1)
        {
          if (fBestObviousCosmic_FlashChi2 > T0_v.at(0)->TriggerConfidence())
          {
            fBestObviousCosmic_FlashChi2 = T0_v.at(0)->TriggerConfidence();
          }
        }
      }
    }
  }
}

bool NuCCproducer::FillDaughters(const art::Ptr<recob::PFParticle> &pfp,
                                 const art::FindManyP<anab::ParticleID> &trackPIDAssn)
{
  clearDaughter();
  if (particlesToVertices.find(pfp) == particlesToVertices.end())
  {
    // If a daughter has no associated vertex, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NuCCproducer::FillDaughters] Daughter had no associated vertex." << std::endl;
    return false;
  }

  if (particlesToMetadata.at(pfp).size() != 1 || particlesToVertices.at(pfp).size() != 1)
  {
    std::cout << "[NuCCproducer::FillDaughters] Daughter association problem." << std::endl;
    return false;
  }
  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  fVx = pfp_vtx.X();
  fVy = pfp_vtx.Y();
  fVz = pfp_vtx.Z();
  std::vector<float> pfp_start_fid_v(6, m_pfp_start_border);
  fStartContained = IsContained(fVx, fVy, fVz, pfp_start_fid_v);
  if (!fStartContained)
  {
    fDaughtersStartContained = false;
  }
  fGeneration = larpandora.GetGeneration(particleMap, pfp);

  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");
  fVtxDistance = pandoraInterfaceHelper.Distance3D(fVx, fVy, fVz, fNu_Vx, fNu_Vy, fNu_Vz);

  // Track-like fields
  if (particlesToTracks.find(pfp) != particlesToTracks.end())
  {
    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
    fTrackLength = this_track->Length();

    // PID information:
    std::map<std::string, float> pid_map;
    if (trackHelper.getPID(pid_map, this_track, trackPIDAssn))
    {
      fTrackPID_chiproton = pid_map.at("chi2_proton");
      fTrackPID_chimuon = pid_map.at("chi2_muon");
    }
    else
    {
      std::cout << "[NuCCproducer::FillDaughters] Track has no PID attached to it" << std::endl;
    }

    if (IsMuonCandidate())
    {
      // add pfp pointer to vector
      m_muon_candidates.push_back(pfp);
    }
  }
  return true;
}

bool NuCCproducer::IsContained(float x, float y, float z, const std::vector<float> &borders) const
{
  float fidvolXstart = borders[0];
  float fidvolYstart = borders[1];
  float fidvolZstart = borders[2];
  float fidvolXend = borders[3];
  float fidvolYend = borders[4];
  float fidvolZend = borders[5];

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x = x > (bnd[0] + fidvolXstart) && x < (bnd[1] - fidvolXend);
  bool is_y = y > (bnd[2] + fidvolYstart) && y < (bnd[3] - fidvolYend);
  bool is_z = z > (bnd[4] + fidvolZstart) && z < (bnd[5] - fidvolZend);

  return is_x && is_y && is_z;
}

bool NuCCproducer::IsMuonCandidate()
{
  fIsMuonCandidate = fGeneration == 2 &&
                     m_muon_cut_trackscore < fTrackScore &&
                     m_muon_cut_vtxdistance > fVtxDistance &&
                     m_muon_cut_protonchi2 < fTrackPID_chiproton &&
                     m_muon_cut_muonchi2 > fTrackPID_chimuon &&
                     m_muon_cut_length < fTrackLength &&
                     m_muon_cut_chiratio < (fTrackPID_chiproton / fTrackPID_chimuon);

  return fIsMuonCandidate;
}

bool NuCCproducer::IsNuMuCC(art::Event &e,
                            std::unique_ptr<std::vector<anab::T0>> &T0_v,
                            std::unique_ptr<art::Assns<recob::PFParticle, anab::T0>> &pfp_t0_assn_v)
{
  if (fNu_PDG == 14 &&
      fNu_Score > m_event_cut_nuscore_hard &&
      (fNu_Score > m_event_cut_nuscore_soft || fNu_FlashChi2 < m_event_cut_flashchi2) &&
      fNu_Contained &&
      fDaughtersStartContained &&
      ((fNu_FlashChi2 / fBestObviousCosmic_FlashChi2) < m_event_cut_flashchi2_ratio))
  {
    // Check if there is a muon candidate:
    if (m_muon_candidates.size() == 0)
    {
      std::cout << "[NuCCproducer::IsNuMuCC] Passed cosmic rejection but no muon candidate found." << std::endl;
      return false;
    }
    else
    {
      float max_length = 0;
      art::Ptr<recob::PFParticle> muon_candidate;

      for (const art::Ptr<recob::PFParticle> candidate : m_muon_candidates)
      {
        const art::Ptr<recob::Track> this_track = particlesToTracks.at(candidate).front();
        if (this_track->Length() > max_length)
        {
          max_length = this_track->Length();
          muon_candidate = candidate;
        }
      }
      if (max_length > m_event_cut_length)
      {
        std::cout << "[NuCCproducer::IsNuMuCC] Passed cosmic rejection and muon candidates found, the longest was picked! SELECTED" << std::endl;
        anab::T0 t0(1, 1, 1);
        T0_v->emplace_back(t0);
        util::CreateAssn(*this, e, *T0_v, muon_candidate, *pfp_t0_assn_v);

        return true;
      }
      else // Muon candidate is not long enough
      {
        return false;
      }
    }
  }
  else // Fails background rejection cuts
  {
    return false;
  }
}
