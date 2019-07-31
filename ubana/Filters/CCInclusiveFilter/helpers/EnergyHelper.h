#ifndef ENERGYHELPER_H
#define ENERGYHELPER_H

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TPrincipal.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

class EnergyHelper
{
public:
  explicit EnergyHelper();
  ~EnergyHelper() = default;

  /**
    * @brief Configure all of the parameters of this class
    *
    * @param p fcl parameter set
    */
  void reconfigure(fhicl::ParameterSet const &pset);


  /**
   * @brief      Measure calorimetric energy for a reconstructed object
   *
   * @param[in]  clusters          Pointer to the vector of reconstructed clusters
   * @param[out] nHits             Address of the vector of the number of hits per plane
   * @param[out] pfenergy          Address of the vector of reconstructed energy per plane
   */
  void energy_from_hits(const lar_pandora::ClusterVector &clusters,
                        std::vector<uint> &nHits,
                        std::vector<float> &pfenergy);


  /**
   * @brief      Measure the dQdx of a shower
   *
   * @param[in]  pfp_dir             shower pfp direction
   * @param[in]  clusters            clusters
   * @param[in]  hits_per_cluster    hits-cluster map
   * @param[out] dqdx                Vector of the dQ/dx median values per plane
   * @param[out] dqdx_hits           Vector of the dQ/dx hits values per plane
   * @param[out] pitches             Vector of the pitches per plane
   */
  void dQdx(const TVector3 &pfp_dir,
            const lar_pandora::ClusterVector &clusters,
            const lar_pandora::ClustersToHits &hits_per_cluster,
            std::vector<float> &dqdx,
            std::vector<std::vector<float>> &dqdx_hits,
            std::vector<float> &pitches);

  /**
   * @brief      Convert dQ/dx vector into dE/dx vector (in MeV)
   *
   * @return dedx          Address of the dE/dx vector
   * @param[in]  dqdx          dQ/dx vector
   */
  std::vector<float> dEdx_from_dQdx(std::vector<float> dqdx);

  /**
   * @brief      Returns the 3D effective pitch give a direction and a plane
   *
   * @param[in]  direction  The direction
   * @param[in]  pl         Plane of interest (0, 1, 2)
   *
   * @return     The pitch.
   */
  float getPitch(const TVector3 &direction, const int &pl) const;

  /**
   * @brief      Builds a rectangle.
   *
   * @param[in]  length  The length
   * @param[in]  width   The width
   * @param      start   The start
   * @param      axis    The axis
   * @return      points  The points
   */
  std::vector<std::vector<float>> buildRectangle(float length, float width, std::vector<float> &start,
                                                 std::vector<float> &axis);

  /**
   * @brief Determine if a point is within a rectangle (in 2D)
   *
   * @param P   vector containing the point's coordinates
   * @param V   vector of vectors containing rectangle's vertices
   * @return    True if the points is inside, False otherwise
   * https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
   */
  bool isInside(std::vector<float> P, std::vector<std::vector<float>> V);

private:
  std::vector<float> _data_gain = {238.4, 238.4, 238.4}; // DocDB 20227
  std::vector<float> _mc_gain = {248.2, 248.2, 248.2};   // Plane 0, plane 1, plane 2
  std::vector<float> _gain;
  const detinfo::DetectorProperties *_detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  float _drift = _detprop->DriftVelocity() * 1e-3;
  float _readout_window = 4.8;
  float _from_tick_to_ns = _readout_window / _detprop->ReadOutWindowSize() * 1e6;
  float _wire_spacing = 0.3;
  float _work_function = 23 / 1e6;
  float _betap;
  float _alpha;
  float _recombination_factor;
  float _dQdx_rectangle_length;
  float _dQdx_rectangle_width;
  bool m_isData;
};

#endif
