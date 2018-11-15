/******************************************************************************
 * @file CandidateTreeFiller.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  CandidateTreeFiller.h
 * ****************************************************************************/

// Decay vertex header
#include "CandidateTreeFiller.h"

namespace AuxEvent
{
  CandidateTreeFiller::CandidateTreeFiller()
  {}
  CandidateTreeFiller::~CandidateTreeFiller()
  {}

  void CandidateTreeFiller::Initialize( AuxEvent::EventTreeFiller & etf, int i_hsnID, const AuxVertex::DecayVertex & dv, std::vector<double> centerCoordinates)
  {
    // General
    run = etf.run;
    subrun = etf.subrun;
    event = etf.event;
    nHsnCandidatesInSameEvent = etf.nHsnCandidates;
    hsnID = i_hsnID;
    // Cheat reco-truth
    if (etf.recoTruthDistances.size() == 0)
    {
      recoTruthDistance = -999;
      isClosestToTruth = false;
      truthCoordinates = {-999,-999,-999};
    }
    else
    {
      truthCoordinates = {etf.truth_vx,etf.truth_vy,etf.truth_vz};
      recoTruthDistance = etf.recoTruthDistances[hsnID];
      isClosestToTruth = etf.isClosestToTruth[hsnID];
    }
    // Coordinates
    geo_nuPosX = dv.fX;
    geo_nuPosY = dv.fY;
    geo_nuPosZ = dv.fZ;
    geo_prongPosX = {dv.fProngX[0],dv.fProngX[1]};
    geo_prongPosY = {dv.fProngY[0],dv.fProngY[1]};
    geo_prongPosZ = {dv.fProngZ[0],dv.fProngZ[1]};
    geo_prongStartPosX = {dv.fProngStartX[0],dv.fProngStartX[1]};
    geo_prongStartPosY = {dv.fProngStartY[0],dv.fProngStartY[1]};
    geo_prongStartPosZ = {dv.fProngStartZ[0],dv.fProngStartZ[1]};
    geo_prongEndPosX = {dv.fProngEndX[0],dv.fProngEndX[1]};
    geo_prongEndPosY = {dv.fProngEndY[0],dv.fProngEndY[1]};
    geo_prongEndPosZ = {dv.fProngEndZ[0],dv.fProngEndZ[1]};
    geo_prongLength = {dv.fProngLength[0],dv.fProngLength[1]};
    geo_openingAngle = dv.fOpeningAngle;
    // Direction
    geo_prongDirX = {dv.fProngDirX[0],dv.fProngDirX[1]};
    geo_prongDirY = {dv.fProngDirY[0],dv.fProngDirY[1]};
    geo_prongDirZ = {dv.fProngDirZ[0],dv.fProngDirZ[1]};
    geo_prongTheta = {dv.fProngTheta[0],dv.fProngTheta[1]};
    geo_prongPhi = {dv.fProngPhi[0],dv.fProngPhi[1]};
    // Prong momentum (by range, assuming both muons)
    hypo_prongPdgCode_h1 = {dv.fProngPdgCode_h1[0],dv.fProngPdgCode_h1[1]};
    hypo_prongMass_h1 = {dv.fProngMass_h1[0],dv.fProngMass_h1[1]};
    range_prongMomMag_h1 = {dv.fProngMomMag_ByRange_h1[0],dv.fProngMomMag_ByRange_h1[1]};
    range_prongEnergy_h1 = {dv.fProngEnergy_ByRange_h1[0],dv.fProngEnergy_ByRange_h1[1]};
    range_prongMom_h1_X = {dv.fProngMom_ByRange_h1_X[0],dv.fProngMom_ByRange_h1_X[1]};
    range_prongMom_h1_Y = {dv.fProngMom_ByRange_h1_Y[0],dv.fProngMom_ByRange_h1_Y[1]};
    range_prongMom_h1_Z = {dv.fProngMom_ByRange_h1_Z[0],dv.fProngMom_ByRange_h1_Z[1]};
    //
    hypo_prongPdgCode_h2 = {dv.fProngPdgCode_h2[0],dv.fProngPdgCode_h2[1]};
    hypo_prongMass_h2 = {dv.fProngMass_h2[0],dv.fProngMass_h2[1]};
    range_prongMomMag_h2 = {dv.fProngMomMag_ByRange_h2[0],dv.fProngMomMag_ByRange_h2[1]};
    range_prongEnergy_h2 = {dv.fProngEnergy_ByRange_h2[0],dv.fProngEnergy_ByRange_h2[1]};
    range_prongMom_h2_X = {dv.fProngMom_ByRange_h2_X[0],dv.fProngMom_ByRange_h2_X[1]};
    range_prongMom_h2_Y = {dv.fProngMom_ByRange_h2_Y[0],dv.fProngMom_ByRange_h2_Y[1]};
    range_prongMom_h2_Z = {dv.fProngMom_ByRange_h2_Z[0],dv.fProngMom_ByRange_h2_Z[1]};
    // Tot momentum (by direction, assuming both muons)
    range_totMomMag_h1 = dv.fTotMomMag_ByRange_h1;
    range_totEnergy_h1 = dv.fTotEnergy_ByRange_h1;
    range_invariantMass_h1 = dv.fInvMass_ByRange_h1;
    range_totMom_h1_X = dv.fTotMom_ByRange_h1_X;
    range_totMom_h1_Y = dv.fTotMom_ByRange_h1_Y;
    range_totMom_h1_Z = dv.fTotMom_ByRange_h1_Z;
    //
    range_totMomMag_h2 = dv.fTotMomMag_ByRange_h2;
    range_totEnergy_h2 = dv.fTotEnergy_ByRange_h2;
    range_invariantMass_h2 = dv.fInvMass_ByRange_h2;
    range_totMom_h2_X = dv.fTotMom_ByRange_h2_X;
    range_totMom_h2_Y = dv.fTotMom_ByRange_h2_Y;
    range_totMom_h2_Z = dv.fTotMom_ByRange_h2_Z;
    // Tot momentum direction (by direction, assuming both muons)
    range_totTheta_h1 = dv.fTotTheta_ByRange_h1;
    range_totPhi_h1 = dv.fTotPhi_ByRange_h1;
    range_totDir_h1_X = dv.fTotDir_ByRange_h1_X;
    range_totDir_h1_Y = dv.fTotDir_ByRange_h1_Y;
    range_totDir_h1_Z = dv.fTotDir_ByRange_h1_Z;
    //
    range_totTheta_h2 = dv.fTotTheta_ByRange_h2;
    range_totPhi_h2 = dv.fTotPhi_ByRange_h2;
    range_totDir_h2_X = dv.fTotDir_ByRange_h2_X;
    range_totDir_h2_Y = dv.fTotDir_ByRange_h2_Y;
    range_totDir_h2_Z = dv.fTotDir_ByRange_h2_Z;
    // Momentum (By Mcs)
    mcs_prongPdgCodeHypothesis = dv.fProngPdgCodeHypothesis_ByMcs;
    mcs_prongIsBestFwd = dv.fProngIsBestFwd_ByMcs;
    // Prong Momentum (By Mcs, best)
    mcs_prongMomMag_best_h1 = dv.fProngMomMag_ByMcs_best_h1;
    mcs_prongEnergy_best_h1 = dv.fProngEnergy_ByMcs_best_h1;
    mcs_prongMom_best_h1_X = dv.fProngMom_ByMcs_best_h1_X;
    mcs_prongMom_best_h1_Y = dv.fProngMom_ByMcs_best_h1_Y;
    mcs_prongMom_best_h1_Z = dv.fProngMom_ByMcs_best_h1_Z;
    mcs_prongMomMag_best_h2 = dv.fProngMomMag_ByMcs_best_h2;
    mcs_prongEnergy_best_h2 = dv.fProngEnergy_ByMcs_best_h2;
    mcs_prongMom_best_h2_X = dv.fProngMom_ByMcs_best_h2_X;
    mcs_prongMom_best_h2_Y = dv.fProngMom_ByMcs_best_h2_Y;
    mcs_prongMom_best_h2_Z = dv.fProngMom_ByMcs_best_h2_Z;
    // Tot momentum (by range, assuming both muons, best)
    mcs_totMomMag_best_h1 = dv.fTotMomMag_ByMcs_best_h1;
    mcs_totEnergy_best_h1 = dv.fTotEnergy_ByMcs_best_h1;
    mcs_invariantMass_best_h1 = dv.fInvMass_ByMcs_best_h1;
    mcs_totMom_best_h1_X = dv.fTotMom_ByMcs_best_h1_X;
    mcs_totMom_best_h1_Y = dv.fTotMom_ByMcs_best_h1_Y;
    mcs_totMom_best_h1_Z = dv.fTotMom_ByMcs_best_h1_Z;
    mcs_totMomMag_best_h2 = dv.fTotMomMag_ByMcs_best_h2;
    mcs_totEnergy_best_h2 = dv.fTotEnergy_ByMcs_best_h2;
    mcs_invariantMass_best_h2 = dv.fInvMass_ByMcs_best_h2;
    mcs_totMom_best_h2_X = dv.fTotMom_ByMcs_best_h2_X;
    mcs_totMom_best_h2_Y = dv.fTotMom_ByMcs_best_h2_Y;
    mcs_totMom_best_h2_Z = dv.fTotMom_ByMcs_best_h2_Z;
    // Tot momentum direction (by range, assuming both muons, best)
    mcs_totTheta_best_h1 = dv.fTotTheta_ByMcs_best_h1;
    mcs_totPhi_best_h1 = dv.fTotPhi_ByMcs_best_h1;
    mcs_totDir_best_h1_X = dv.fTotDir_ByMcs_best_h1_X;
    mcs_totDir_best_h1_Y = dv.fTotDir_ByMcs_best_h1_Y;
    mcs_totDir_best_h1_Z = dv.fTotDir_ByMcs_best_h1_Z;
    mcs_totTheta_best_h2 = dv.fTotTheta_ByMcs_best_h2;
    mcs_totPhi_best_h2 = dv.fTotPhi_ByMcs_best_h2;
    mcs_totDir_best_h2_X = dv.fTotDir_ByMcs_best_h2_X;
    mcs_totDir_best_h2_Y = dv.fTotDir_ByMcs_best_h2_Y;
    mcs_totDir_best_h2_Z = dv.fTotDir_ByMcs_best_h2_Z;
    // Extra
    prongNumHits = {dv.fProngNumHits[0],dv.fProngNumHits[1]};
    prongStartToNeutrinoDistance = {dv.fProngStartToNeutrinoDistance[0],dv.fProngStartToNeutrinoDistance[1]};

    // Calculate extra quantities
    // Calculate end points
    float disp1 = abs(geo_prongEndPosX[0] - centerCoordinates[0]);
    float disp2 = abs(geo_prongEndPosX[1] - centerCoordinates[0]);
    if (disp1 > disp2) maxEndPointX = geo_prongEndPosX[0];
    else maxEndPointX = geo_prongEndPosX[1];

    disp1 = abs(geo_prongEndPosY[0] - centerCoordinates[1]);
    disp2 = abs(geo_prongEndPosY[1] - centerCoordinates[1]);
    if (disp1 > disp2) maxEndPointY = geo_prongEndPosY[0];
    else maxEndPointY = geo_prongEndPosY[1];

    disp1 = abs(geo_prongEndPosZ[0] - centerCoordinates[2]);
    disp2 = abs(geo_prongEndPosZ[1] - centerCoordinates[2]);
    if (disp1 > disp2) maxEndPointZ = geo_prongEndPosZ[0];
    else maxEndPointZ = geo_prongEndPosZ[1];

    // Calculate delta phi and theta
    deltaPhi = abs(geo_prongPhi[0] - geo_prongPhi[1]);
    deltaTheta = abs(geo_prongTheta[0] - geo_prongTheta[1]);

    // Calculate length ratio and diff, and maxDist
    lengthDiff = abs(geo_prongLength[0] - geo_prongLength[1]);
    lengthRatio = std::min(geo_prongLength[0],geo_prongLength[1])/float(std::max(geo_prongLength[0],geo_prongLength[1]));
    maxStartToNeutrinoDistance = std::max(prongStartToNeutrinoDistance[0],prongStartToNeutrinoDistance[1]);


  } // END function Initialize
} // END namespace CandidateTreeFiller 
