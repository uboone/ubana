////////////////////////////////////////////////////////////////////////
// Class:       NuCCproducer
// Plugin Type: producer (art v3_01_02)
// File:        NuCCproducer_module.cc
//
// Generated at Wed Jul 31 09:22:00 2019 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#ifndef NUCCPRODUCER_H
#define NUCCPRODUCER_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "helpers/TrackHelper.h"
#include "helpers/PandoraInterfaceHelper.h"

class NuCCproducer;

class NuCCproducer : public art::EDProducer
{
public:
  explicit NuCCproducer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuCCproducer(NuCCproducer const &) = delete;
  NuCCproducer(NuCCproducer &&) = delete;
  NuCCproducer &operator=(NuCCproducer const &) = delete;
  NuCCproducer &operator=(NuCCproducer &&) = delete;

  // Required functions.
  void produce(art::Event &e) override;

  void reconfigure(fhicl::ParameterSet const &p);
  void clearEvent();
  void clearDaughter();

  /**
     *  @brief  Collect and fill the reconstructed information.
     *
     *  @param  e Art event
     */
  void FillReconstructed(art::Event const &e);

  /**
     *  @brief  Fill the tree for every daughter of the neutrino candidate.
     *
     *  @param  pfparticle ptr The Pfparticle corresponding to the daughter.
     *  @return 1 if succesful, 0 if failure.
     */
  bool FillDaughters(const art::Ptr<recob::PFParticle> &pfp,
                     const art::FindManyP<anab::ParticleID> &trackPIDAssn);
  /**
     *  @brief  Tag a daughter as a muon candidate
     *
     *  @return 1 if succesfully, 0 if not.
     */
  bool IsMuonCandidate();

  /**
     *  @brief  Tag event as NuMuCC, if so, create association
     *
     *  @param Event e, to create the association.
     *  @return 1 if succesfully, 0 if not.
     */
  bool IsNuMuCC(art::Event &e, 
                std::unique_ptr<std::vector<anab::T0>> &T0_v,
                std::unique_ptr<art::Assns<recob::PFParticle, anab::T0>> &pfp_t0_assn_v);

  /**
     *  @brief  Returns if point is inside a fiducial volume
     *
     *  @param fiducial volume tolerance: -x,+x,-y,+y,-z,+z
     *  @return 1 if succesfully, 0 if not.
     */
  bool IsContained(float x, float y, float z, const std::vector<float> &borders) const;

private:
  // Fields needed for the analyser
  std::string m_pfp_producer;
  std::string m_hitfinder_producer;

  float m_vtx_fid_x_start;
  float m_vtx_fid_y_start;
  float m_vtx_fid_z_start;
  float m_vtx_fid_x_end;
  float m_vtx_fid_y_end;
  float m_vtx_fid_z_end;
  float m_pfp_start_border;

  float m_muon_cut_trackscore;
  float m_muon_cut_vtxdistance;
  float m_muon_cut_protonchi2;
  float m_muon_cut_muonchi2;
  float m_muon_cut_chiratio;
  float m_muon_cut_length;

  float m_event_cut_flashchi2;
  float m_event_cut_nuscore_soft;
  float m_event_cut_nuscore_hard;
  float m_event_cut_flashchi2_ratio;
  float m_event_cut_length;

  TrackHelper trackHelper;
  PandoraInterfaceHelper pandoraInterfaceHelper;

  // Store the pfps that qualify as muon candidates
  std::vector<art::Ptr<recob::PFParticle>> m_muon_candidates;

  // LAr Pandora Helper fields
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::TrackVector pftracks;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToTracks particlesToTracks;

  uint fRun, fSubrun, fEvent;
  // Reco candidate info
  int fNu_PDG;
  float fNu_Score;
  float fNu_Vx, fNu_Vy, fNu_Vz;
  bool fNu_Contained;
  bool fDaughtersStartContained;
  float fNu_FlashChi2;
  float fBestObviousCosmic_FlashChi2;
  bool fIsNuMuCC;

  // Reco daughter info
  uint fGeneration;
  float fTrackScore;
  float fVx, fVy, fVz;
  bool fStartContained;
  float fVtxDistance;
  // Track info
  float fTrackLength;
  float fTrackPID_chiproton;
  float fTrackPID_chimuon;
  bool fIsMuonCandidate;
};

void NuCCproducer::reconfigure(fhicl::ParameterSet const &p)
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandora");
  m_hitfinder_producer = p.get<std::string>("hitfinder_producer", "gaushit");

  m_vtx_fid_x_start = p.get<float>("vtx_fid_x_start", 10);
  m_vtx_fid_y_start = p.get<float>("vtx_fid_y_start", 10);
  m_vtx_fid_z_start = p.get<float>("vtx_fid_z_start", 10);
  m_vtx_fid_x_end = p.get<float>("vtx_fid_x_end", 10);
  m_vtx_fid_y_end = p.get<float>("vtx_fid_y_end", 10);
  m_vtx_fid_z_end = p.get<float>("vtx_fid_z_end", 50);
  m_pfp_start_border = p.get<float>("pfp_start_border", 10);

  m_muon_cut_trackscore = p.get<float>("muon_cut_trackscore", 0.8);
  m_muon_cut_vtxdistance = p.get<float>("muon_cut_vtxdistance", 4.0);
  m_muon_cut_protonchi2 = p.get<float>("muon_cut_protonchi2", 60);
  m_muon_cut_muonchi2 = p.get<float>("muon_cut_muonchi2", 30);
  m_muon_cut_chiratio = p.get<float>("muon_cut_chiratio", 7);
  m_muon_cut_length = p.get<float>("muon_cut_length", 5);

  m_event_cut_flashchi2 = p.get<float>("event_cut_flashchi2", 10);
  m_event_cut_nuscore_soft = p.get<float>("event_cut_nuscore_soft", 0.25);
  m_event_cut_nuscore_hard = p.get<float>("event_cut_nuscore_hard", 0.06);
  m_event_cut_flashchi2_ratio = p.get<float>("event_cut_flashchi2_ratio", 5);
  m_event_cut_length = p.get<float>("event_cut_length", 20);
}

NuCCproducer::NuCCproducer(fhicl::ParameterSet const &p)
    : EDProducer{p}
{
  this->reconfigure(p);

  produces<std::vector<anab::T0>>();
  produces<art::Assns<recob::PFParticle, anab::T0>>();
}

void NuCCproducer::clearEvent()
{
  fNu_PDG = 0; // if 0, no neutrinocandidate was found, only look at truth information.
  fIsNuMuCC = false;
  fDaughtersStartContained = true;
  fNu_FlashChi2 = 0;
  fBestObviousCosmic_FlashChi2 = 0;

  pfparticles.clear();
  pfneutrinos.clear();
  pfdaughters.clear();
  pftracks.clear();
  particleMap.clear();
  particlesToMetadata.clear();
  particlesToVertices.clear();
  particlesToTracks.clear();

  m_muon_candidates.clear();
}

void NuCCproducer::clearDaughter()
{
  // Track info
  fTrackLength = 0;
  fTrackPID_chiproton = 0;
  fTrackPID_chimuon = 0;
  fIsMuonCandidate = false;
}

DEFINE_ART_MODULE(NuCCproducer)
#endif // NUCCPRODUCER_H
