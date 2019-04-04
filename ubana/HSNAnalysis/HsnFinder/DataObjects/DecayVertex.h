/******************************************************************************
 * @file DecayVertex.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.cxx
 * ****************************************************************************/

#ifndef DECAYVERTEX_H
#define DECAYVERTEX_H

#include "TMatrixD.h"
#include "TVector3.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// #include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace AuxVertex
{

  // Decay vertex class and functions
  class DecayVertex
  {
  public:
    // Constructor and destructor
    DecayVertex();
    virtual ~DecayVertex();

    DecayVertex(
            const art::Ptr<recob::Vertex> &nuVertex,
            const art::Ptr<recob::Vertex> &t1Vertex,
            const art::Ptr<recob::Vertex> &t2Vertex,
            const art::Ptr<recob::Track> &t1Track,
            const art::Ptr<recob::Track> &t2Track,
            const std::vector<art::Ptr<recob::Hit>> &t1Hits,
            const std::vector<art::Ptr<recob::Hit>> &t2Hits,
            const art::Ptr<recob::MCSFitResult> &t1Mcs,
            const art::Ptr<recob::MCSFitResult> &t2Mcs);

    // New Getters
    art::Ptr<recob::Vertex> GetNuVertex() const;
    art::Ptr<recob::Vertex> GetProngVertex(int prong) const;
    art::Ptr<recob::Track> GetProngTrack(int prong) const;
    std::vector<art::Ptr<recob::Hit>> GetProngHits(int prong) const;
    std::vector<art::Ptr<recob::Hit>> GetTotHits() const;

    // Setters
    void SetDetectorCoordinates(
      const std::vector<double>& minTpcBound,
      const std::vector<double>& maxTpcBound,
      geo::GeometryCore const* geometry,
      detinfo::DetectorProperties const* detectorProperties);
    void SetChannelLoc(int channel0, int channel1, int channel2);
    void SetTickLoc(float tick0, float tick1, float tick2);
    void SetProngChannelLoc(int par, int channel0, int channel1, int channel2);
    void SetProngTickLoc(int par, float tick0, float tick1, float tick2);
    void SetProngXYZ(int par, float x, float y, float z);
    void SetIsInsideTPC(bool val);
    void SetIsDetLocAssigned(bool val);
    void SetTotHits(std::vector<art::Ptr<recob::Hit>> totHitsInMaxRadius);
    void SetHypothesisLabels();
    void SetMomentumQuantities_ByRange();
    void SetMomentumQuantities_ByMCS();


    // Printers
    void PrintInformation() const;

    // Data products pointers
    art::Ptr<recob::Vertex> fNuVertex;
    std::vector<art::Ptr<recob::Vertex>> fProngVertex;
    std::vector<art::Ptr<recob::Track>> fProngTrack;
    std::vector<art::Ptr<recob::MCSFitResult>> fProngMcs;
    std::vector<std::vector<art::Ptr<recob::Hit>>> fProngHits;
    std::vector<art::Ptr<recob::Hit>> fTotHitsInMaxRadius;

    // Coordinates of the pandora neutrino recob::Vertex object
    float fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
    std::vector<int> fChannelLoc; // Nearest channel in each plane.
    std::vector<float> fTickLoc; // Nearest time tick in each plane.
    /**/
    // Coordinates of the two prongs recob::Vertex objects
    std::vector<float> fProngX, fProngY, fProngZ; // Spatial coordinates of the vertex of the track inside the detector.
    std::vector<std::vector<int>> fProngChannelLoc; // Nearest channel in each plane for the vertex parent.
    std::vector<std::vector<float>> fProngTickLoc; // Nearest time tick in each plane for the vertex parent.
    /**/
    // Coordinates of the two prongs start and end points for the recob::Track objects
    std::vector<float> fProngStartX, fProngStartY, fProngStartZ; // Spatial coordinates for the start of the track.
    std::vector<float> fProngEndX, fProngEndY, fProngEndZ; // Spatial coordinates for the end of the track.

    // Track direction (no calorimetry data)
    std::vector<float> fProngDirX, fProngDirY, fProngDirZ; // Direction of momentum for each track.
    std::vector<float> fProngTheta, fProngPhi; // Direction angles of each prong.

    // Hypotheses for pdg code and masses
    std::vector<int> fProngPdgCode_h1, fProngPdgCode_h2;
    std::vector<float> fProngMass_h1, fProngMass_h2;

    // Momentum information (measured by range and using two different hypotheses for pdg)
    std::vector<float> fProngMom_ByRange_h1_X, fProngMom_ByRange_h1_Y, fProngMom_ByRange_h1_Z; // Components of momentum.
    std::vector<float> fProngMom_ByRange_h2_X, fProngMom_ByRange_h2_Y, fProngMom_ByRange_h2_Z; // Components of momentum.
    std::vector<float> fProngMomMag_ByRange_h1, fProngEnergy_ByRange_h1; // Momentum and energy of each prong.
    std::vector<float> fProngMomMag_ByRange_h2, fProngEnergy_ByRange_h2; // Momentum and energy of each prong.
    /**/
    float fTotMom_ByRange_h1_X, fTotMom_ByRange_h1_Y, fTotMom_ByRange_h1_Z; // Momentum component for neutrino.
    float fTotMom_ByRange_h2_X, fTotMom_ByRange_h2_Y, fTotMom_ByRange_h2_Z; // Momentum component for neutrino.
    float fTotDir_ByRange_h1_X, fTotDir_ByRange_h1_Y, fTotDir_ByRange_h1_Z; // Direction components of neutrino
    float fTotDir_ByRange_h2_X, fTotDir_ByRange_h2_Y, fTotDir_ByRange_h2_Z; // Direction components of neutrino
    float fTotTheta_ByRange_h1, fTotPhi_ByRange_h1; // Direction angles of neutrino
    float fTotTheta_ByRange_h2, fTotPhi_ByRange_h2; // Direction angles of neutrino
    float fTotMomMag_ByRange_h1, fTotEnergy_ByRange_h1, fInvMass_ByRange_h1; // Total momentum, total energy and invariant mass.
    float fTotMomMag_ByRange_h2, fTotEnergy_ByRange_h2, fInvMass_ByRange_h2; // Total momentum, total energy and invariant mass.

    // Momentum information (measured by MCS and using two different hypotheses for Pdg)
    std::vector<int> fProngPdgCodeHypothesis_ByMcs;
    std::vector<bool> fProngIsBestFwd_ByMcs;
    std::vector<float> fProngMomMag_ByMcs_fwd_h1, fProngMomMag_ByMcs_best_h1; // Momentum of prongs
    std::vector<float> fProngMomMag_ByMcs_fwd_h2, fProngMomMag_ByMcs_best_h2; // Momentum of prongs
    std::vector<float> fProngMom_ByMcs_fwd_h1_X, fProngMom_ByMcs_fwd_h1_Y, fProngMom_ByMcs_fwd_h1_Z; // Components of momentum.
    std::vector<float> fProngMom_ByMcs_best_h1_X, fProngMom_ByMcs_best_h1_Y, fProngMom_ByMcs_best_h1_Z; // Components of momentum.
    std::vector<float> fProngMom_ByMcs_fwd_h2_X, fProngMom_ByMcs_fwd_h2_Y, fProngMom_ByMcs_fwd_h2_Z; // Components of momentum.
    std::vector<float> fProngMom_ByMcs_best_h2_X, fProngMom_ByMcs_best_h2_Y, fProngMom_ByMcs_best_h2_Z; // Components of momentum.
    /**/
    std::vector<float> fProngEnergy_ByMcs_fwd_h1, fProngEnergy_ByMcs_best_h1; // Energy of each prong.
    std::vector<float> fProngEnergy_ByMcs_fwd_h2, fProngEnergy_ByMcs_best_h2; // Energy of each prong.
    /**/
    float fTotMom_ByMcs_fwd_h1_X, fTotMom_ByMcs_fwd_h1_Y, fTotMom_ByMcs_fwd_h1_Z; // Momentum component for neutrino.
    float fTotMom_ByMcs_best_h1_X, fTotMom_ByMcs_best_h1_Y, fTotMom_ByMcs_best_h1_Z; // Momentum component for neutrino.
    float fTotDir_ByMcs_fwd_h1_X, fTotDir_ByMcs_fwd_h1_Y, fTotDir_ByMcs_fwd_h1_Z; // Direction components of neutrino
    float fTotDir_ByMcs_best_h1_X, fTotDir_ByMcs_best_h1_Y, fTotDir_ByMcs_best_h1_Z; // Direction components of neutrino
    float fTotTheta_ByMcs_fwd_h1, fTotPhi_ByMcs_fwd_h1; // Direction angles of neutrino
    float fTotTheta_ByMcs_best_h1, fTotPhi_ByMcs_best_h1; // Direction angles of neutrino
    float fTotMomMag_ByMcs_fwd_h1, fTotEnergy_ByMcs_fwd_h1, fInvMass_ByMcs_fwd_h1; // Total momentum, total energy and invariant mass.
    float fTotMomMag_ByMcs_best_h1, fTotEnergy_ByMcs_best_h1, fInvMass_ByMcs_best_h1; // Total momentum, total energy and invariant mass.

    float fTotMom_ByMcs_fwd_h2_X, fTotMom_ByMcs_fwd_h2_Y, fTotMom_ByMcs_fwd_h2_Z; // Momentum component for neutrino.
    float fTotMom_ByMcs_best_h2_X, fTotMom_ByMcs_best_h2_Y, fTotMom_ByMcs_best_h2_Z; // Momentum component for neutrino.
    float fTotDir_ByMcs_fwd_h2_X, fTotDir_ByMcs_fwd_h2_Y, fTotDir_ByMcs_fwd_h2_Z; // Direction components of neutrino
    float fTotDir_ByMcs_best_h2_X, fTotDir_ByMcs_best_h2_Y, fTotDir_ByMcs_best_h2_Z; // Direction components of neutrino
    float fTotTheta_ByMcs_fwd_h2, fTotPhi_ByMcs_fwd_h2; // Direction angles of neutrino
    float fTotTheta_ByMcs_best_h2, fTotPhi_ByMcs_best_h2; // Direction angles of neutrino
    float fTotMomMag_ByMcs_fwd_h2, fTotEnergy_ByMcs_fwd_h2, fInvMass_ByMcs_fwd_h2; // Total momentum, total energy and invariant mass.
    float fTotMomMag_ByMcs_best_h2, fTotEnergy_ByMcs_best_h2, fInvMass_ByMcs_best_h2; // Total momentum, total energy and invariant mass.

    // Other variables
    std::vector<float> fProngLength; // Length of each prong
    float fOpeningAngle; // Opening angle between the two prongs
    std::vector<float> fProngStartToNeutrinoDistance; // Distance from start point to neutrino vertex for each prong
    std::vector<int> fProngNumHits; // Number of hits associated with each prong

    // Status
    bool fIsInsideTPC; // Whetehr the vertex is inside the TPC.
    bool fIsDetLocAssigned; // Whether channel/tick coordinates have been determined.
  };
  
} //END namespace AuxVertex

#endif