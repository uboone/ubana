/******************************************************************************
 * @file CandidateTreeFiller.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  CandidateTreeFiller.cxx
 * ****************************************************************************/

#ifndef CandidateTreeFiller_H
#define CandidateTreeFiller_H

// C++ standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubana/HELEPHANTE/HelephantFinder/DataObjects/DecayVertex.h"
#include "EventTreeFiller.h"


namespace AuxEvent
{

  // CandidateTreeFiller class and functions
  class CandidateTreeFiller
  {
  public:
    // Constructor and destructor
    CandidateTreeFiller();
    virtual ~CandidateTreeFiller();

    void Initialize(AuxEvent::EventTreeFiller & etf, int i_hsnID, const AuxVertex::DecayVertex & decayVertex, std::vector<double> centerCoordinates);

    // General
    int run;
    int subrun;
    int event;
    int hsnID;
    int nHsnCandidatesInSameEvent;
    // Cheat reco-truth
    float recoTruthDistance;
    bool isClosestToTruth;
    std::vector<float> truthCoordinates;
    // Coordinates
    float geo_nuPosX, geo_nuPosY, geo_nuPosZ;
    std::vector<float> geo_prongPosX, geo_prongPosY, geo_prongPosZ;
    std::vector<float> geo_prongStartPosX, geo_prongStartPosY, geo_prongStartPosZ;
    std::vector<float> geo_prongEndPosX, geo_prongEndPosY, geo_prongEndPosZ;
    std::vector<float> geo_prongLength;
    float geo_openingAngle;
    // Direction
    std::vector<float> geo_prongDirX, geo_prongDirY, geo_prongDirZ;
    std::vector<float> geo_prongTheta, geo_prongPhi;
    // Hypothesis information
    std::vector<int> hypo_prongPdgCode_h1, hypo_prongPdgCode_h2;
    std::vector<float> hypo_prongMass_h1, hypo_prongMass_h2;
    // Prong momentum (by range, assuming both muons)
    std::vector<float> range_prongMomMag_h1, range_prongEnergy_h1;
    std::vector<float> range_prongMom_h1_X, range_prongMom_h1_Y, range_prongMom_h1_Z;
    std::vector<float> range_prongMomMag_h2, range_prongEnergy_h2;
    std::vector<float> range_prongMom_h2_X, range_prongMom_h2_Y, range_prongMom_h2_Z;
    // Tot momentum (by range, assuming both muons)
    float range_totMomMag_h1, range_totEnergy_h1, range_invariantMass_h1;
    float range_totMom_h1_X, range_totMom_h1_Y, range_totMom_h1_Z;
    float range_totMomMag_h2, range_totEnergy_h2, range_invariantMass_h2;
    float range_totMom_h2_X, range_totMom_h2_Y, range_totMom_h2_Z;
    // Tot momentum direction (by range, assuming both muons)
    float range_totTheta_h1, range_totPhi_h1;
    float range_totDir_h1_X, range_totDir_h1_Y, range_totDir_h1_Z;
    float range_totTheta_h2, range_totPhi_h2;
    float range_totDir_h2_X, range_totDir_h2_Y, range_totDir_h2_Z;
    // MCS variables
    std::vector<int> mcs_prongPdgCodeHypothesis;
    std::vector<bool>  mcs_prongIsBestFwd;
    // Prong Momentum (By MCS, forward)
    std::vector<float> mcs_prongMomMag_fwd_h1, mcs_prongEnergy_fwd_h1;
    std::vector<float> mcs_prongMom_fwd_h1_X, mcs_prongMom_fwd_h1_Y, mcs_prongMom_fwd_h1_Z;
    std::vector<float> mcs_prongMomMag_fwd_h2, mcs_prongEnergy_fwd_h2;
    std::vector<float> mcs_prongMom_fwd_h2_X, mcs_prongMom_fwd_h2_Y, mcs_prongMom_fwd_h2_Z;
    // Tot momentum (by MCS, assuming both muons, forward)
    float mcs_totMomMag_fwd_h1, mcs_totEnergy_fwd_h1, mcs_invariantMass_fwd_h1;
    float mcs_totMom_fwd_h1_X, mcs_totMom_fwd_h1_Y, _totMom_fwd_h1_Z;
    float mcs_totMomMag_fwd_h2, mcs_totEnergy_fwd_h2, mcs_invariantMass_fwd_h2;
    float mcs_totMom_fwd_h2_X, mcs_totMom_fwd_h2_Y, mcs_totMom_fwd_h2_Z;
    // Tot momentum direction (by MCS, assuming both muons, forward)
    float mcs_totTheta_fwd_h1, mcs_totPhi_fwd_h1;
    float mcs_totDir_fwd_h1_X, mcs_totDir_fwd_h1_Y, mcs_totDir_fwd_h1_Z;
    float mcs_totTheta_fwd_h2, mcs_totPhi_fwd_h2;
    float mcs_totDir_fwd_h2_X, mcs_totDir_fwd_h2_Y, mcs_totDir_fwd_h2_Z;
    // Prong Momentum (By MCS, best)
    std::vector<float> mcs_prongMomMag_best_h1, mcs_prongEnergy_best_h1;
    std::vector<float> mcs_prongMom_best_h1_X, mcs_prongMom_best_h1_Y, mcs_prongMom_best_h1_Z;
    std::vector<float> mcs_prongMomMag_best_h2, mcs_prongEnergy_best_h2;
    std::vector<float> mcs_prongMom_best_h2_X, mcs_prongMom_best_h2_Y, mcs_prongMom_best_h2_Z;
    // Tot momentum (by MCS, assuming both muons, best)
    float mcs_totMomMag_best_h1, mcs_totEnergy_best_h1, mcs_invariantMass_best_h1;
    float mcs_totMom_best_h1_X, mcs_totMom_best_h1_Y, mcs_totMom_best_h1_Z;
    float mcs_totMomMag_best_h2, mcs_totEnergy_best_h2, mcs_invariantMass_best_h2;
    float mcs_totMom_best_h2_X, mcs_totMom_best_h2_Y, mcs_totMom_best_h2_Z;
    // Tot momentum direction (by MCS, assuming both muons, best)
    float mcs_totTheta_best_h1, mcs_totPhi_best_h1;
    float mcs_totDir_best_h1_X, mcs_totDir_best_h1_Y, mcs_totDir_best_h1_Z;
    float mcs_totTheta_best_h2, mcs_totPhi_best_h2;
    float mcs_totDir_best_h2_X, mcs_totDir_best_h2_Y, mcs_totDir_best_h2_Z;
    // Extra
    std::vector<float> prongStartToNeutrinoDistance;
    std::vector<int> prongNumHits;
    float maxEndPointX, maxEndPointY, maxEndPointZ;
    float deltaPhi, deltaTheta;
    float maxStartToNeutrinoDistance;
    float lengthDiff, lengthRatio;
    // Pandora statusnostic
    int status_nuWithMissingAssociatedVertex;
    int status_nuWithMissingAssociatedTrack;
    int status_nuProngWithMissingAssociatedHits;
    // Pandora calo
    std::vector<std::vector<float>> calo_totChargeInRadius;
    std::vector<std::vector<float>> calo_prong1ChargeInRadius;
    std::vector<std::vector<float>> calo_prong2ChargeInRadius;
    std::vector<std::vector<float>> calo_caloRatio;
    // Flash info
    float flash_flashDistance, flash_flashPE;
    // Trigger info
    bool trigger_passedBNB, trigger_passedBNB5PE, trigger_passedEXT, trigger_passedEXT5PE, trigger_passedHSN, trigger_passedHSNEXT;
    // G4 info
    int g4_pi_pdg, g4_mu_pdg;
    float g4_pi_mom, g4_mu_mom;
    float g4_pi_endMom, g4_mu_endMom;
    float g4_reco_pi_momRange, g4_reco_mu_momRange;
    float g4_reco_pi_momMCS, g4_reco_mu_momMCS;
    float g4_pi_length, g4_mu_length;
    int g4_pi_nScatterings, g4_mu_nScatterings;
    std::vector<std::vector<float>> g4_recoTruth_distancesScore;
    std::string g4_pi_declaredEndProcess, g4_mu_declaredEndProcess;
    std::string g4_pi_declaredScatteredEndProcess, g4_mu_declaredScatteredEndProcess;
    std::string g4_pi_endProcess, g4_mu_endProcess;
    // Calo info
    // std::vector<std::vector<bool>> calo_prongIsTrackFlipped = {{false,false,false},{false,false,false}};
    std::vector<std::vector<bool>> cali_prongIsTrackFlipped = {{false,false,false},{false,false,false}};
    std::vector<std::vector<int>> cali_prongCaloPlane = {{-1,-1,-1},{-1,-1,-1}};

  };


} //END namespace AuxEvent

#endif
