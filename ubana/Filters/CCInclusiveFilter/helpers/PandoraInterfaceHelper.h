#ifndef PANDORAINTERFACEHELPER_H
#define PANDORAINTERFACEHELPER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include <math.h>

namespace lar_pandora
{
typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle>> PFParticlesToMCParticles;
}

typedef std::map<art::Ptr<recob::PFParticle>, unsigned int> RecoParticleToNMatchedHits;
typedef std::map<art::Ptr<simb::MCParticle>, RecoParticleToNMatchedHits> ParticleMatchingMap;
typedef std::set<art::Ptr<recob::PFParticle>> PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle>> MCParticleSet;

class PandoraInterfaceHelper
{
public:
  PandoraInterfaceHelper();
  ~PandoraInterfaceHelper() = default;

  /**
    *  @brief Returns matching between true and reconstructed particles
    *
    *  @param matchedParticles the output matches between reconstructed and true particles
    */
  void GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles,
                            std::map<art::Ptr<recob::PFParticle>, float> &matchedHitFractions);

  /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param m_pfp_producer the PFParticle producer label
     *  @param m_spacepoint_producer the SpacePoint producer label
     *  @param m_hitfinder_producer the Hit producer label
     *  @param m_geant_producer The Geant4 producer label
     */
  void Configure(art::Event const &e,
                 std::string m_pfp_producer,
                 std::string m_spacepoint_producer,
                 std::string m_hitfinder_producer,
                 std::string m_geant_producer,
                 std::string m_hit_mcp_producer);

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id);

  /**
 *  @brief  Collect all downstream particles of those in the input vector
 *
 *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
 *  @param  parentPFParticles the input vector of PFParticles
 *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
 */
  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                    const art::Ptr<recob::PFParticle> &particle,
                                    lar_pandora::PFParticleVector &downstreamPFParticles) const;

  /**
 *  @brief  Distance between 2 point in 3D
 *  @return  distance 
 */
  float Distance3D(float x1, float y1, float z1, float x2, float y2, float z2);

/**
 *  @brief  Applies SCE on a point (MCC9)
 */
  void SCE(const float &x, 
            const float &y, 
            const float &z, 
            const float &nu_time, 
            float &x_out, float &y_out, float &z_out);

protected:
  lar_pandora::HitsToMCParticles m_hit_to_mcps_map; ///< A map from recon hits to MCParticles
  lar_pandora::PFParticlesToHits m_pfp_to_hits_map; ///< A map from PFParticles to recon hits

private:
  bool m_configured = false;
  bool m_debug = false;
  bool m_verbose = false;
};

#endif