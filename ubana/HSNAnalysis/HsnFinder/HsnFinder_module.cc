#ifndef HSNFINDER_MODULE
#define HSNFINDER_MODULE

#include "HsnFinder.h"

HsnFinder::HsnFinder(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fFindPandoraVertexAlg(pset),
    fCalorimetryRadiusAlg(pset),
    fExtractTruthInformationAlg(pset),
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fCenterCoordinates(pset.get<std::vector<double>>("CenterCoordinates")),
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fMcsLabel(pset.get<std::string>("McsLabel")),
    fRadiusProfileLimits(pset.get<std::vector<double>>("RadiusProfileLimits")),
    fRadiusProfileBins(pset.get<int>("RadiusProfileBins")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fVerbose(pset.get<bool>("VerboseMode")),
    fSaveDrawTree(pset.get<bool>("SaveDrawTree") ),
    fSaveTruthDrawTree(pset.get<bool>("SaveTruthDrawTree") ),
    fUseTruthDistanceMetric(pset.get<bool>("UseTruthDistanceMetric")),
    fMcTrackLabel(pset.get<std::string>("McTrackLabel")),
    fIsHSN(pset.get<bool>("IsHSN"))
{
  // Get geometry and detector services
  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // Determine profile ticks
  double profileStep = (fRadiusProfileLimits[1] - fRadiusProfileLimits[0]) / float(fRadiusProfileBins);
  double currTick = fRadiusProfileLimits[0];
  for (int i=0; i<fRadiusProfileBins; i++)
  {
    currTick += profileStep;
    profileTicks.push_back(currTick);
  }
} // END constructor HsnFinder

HsnFinder::~HsnFinder()
{} // END destructor HsnFinder

void HsnFinder::beginJob()
{
  // Declare file service handle
  art::ServiceHandle< art::TFileService > tfs;

  // Meta tree containing fcl file parameters
  metaTree = tfs->make<TTree>("MetaData","");
  metaTree->Branch("instanceName",&fInstanceName);
  metaTree->Branch("iteration",&fIteration,"iteration/I");
  metaTree->Branch("minTpcBound",&fMinTpcBound);
  metaTree->Branch("maxTpcBound",&fMaxTpcBound);
  metaTree->Branch("pfpLabel",&fPfpLabel);
  metaTree->Branch("hitLabel",&fHitLabel);
  metaTree->Branch("mcsLabel",&fMcsLabel);
  metaTree->Branch("radiusProfileLimits",&fRadiusProfileLimits);
  metaTree->Branch("radiusProfileBins",&fRadiusProfileBins);
  metaTree->Branch("profileTicks",&profileTicks);
  metaTree->Branch("channelNorm",&fChannelNorm,"channelNorm/D");
  metaTree->Branch("tickNorm",&fTickNorm,"tickNorm/D");
  metaTree->Branch("saveDrawTree",&fSaveDrawTree,"saveDrawTree/O");
  metaTree->Fill();

  // Tree containing data about current event
  eventTree = tfs->make<TTree>("EventData","");
  eventTree->Branch("run",&etf.run);
  eventTree->Branch("subrun",&etf.subrun);
  eventTree->Branch("event",&etf.event);
  eventTree->Branch("nNeutrinos",&etf.nNeutrinos);
  eventTree->Branch("neutrinoPdgCode",&etf.neutrinoPdgCode);
  eventTree->Branch("neutrinoNumDaughters",&etf.neutrinoNumDaughters);
  eventTree->Branch("neutrinoNumTracks",&etf.neutrinoNumTracks);
  eventTree->Branch("neutrinoNumShowers",&etf.neutrinoNumShowers);
  eventTree->Branch("nTwoProngedNeutrinos",&etf.nTwoProngedNeutrinos);
  eventTree->Branch("nContainedTwoProngedNeutrinos",&etf.nContainedTwoProngedNeutrinos);
  eventTree->Branch("nHsnCandidates",&etf.nHsnCandidates);
  eventTree->Branch("nHsnCandidates",&etf.nHsnCandidates);
  eventTree->Branch("truth_vx",&etf.truth_vx);
  eventTree->Branch("truth_vy",&etf.truth_vy);
  eventTree->Branch("truth_vz",&etf.truth_vz);
  eventTree->Branch("recoTruthDistances",&etf.recoTruthDistances);
  eventTree->Branch("isClosestToTruth",&etf.isClosestToTruth);

  // Tree containing data about current HSN candidate
  candidateTree = tfs->make<TTree>("CandidateData","");
  // HSN ID
  candidateTree->Branch("run",&ctf.run);
  candidateTree->Branch("subrun",&ctf.subrun);
  candidateTree->Branch("event",&ctf.event);
  candidateTree->Branch("hsnID",&ctf.hsnID);
  candidateTree->Branch("nHsnCandidatesInSameEvent",&ctf.nHsnCandidatesInSameEvent);
  // Cheat reco-truth
  if ( fUseTruthDistanceMetric )
  {
    candidateTree->Branch("recoTruthDistance",&ctf.recoTruthDistance);
    candidateTree->Branch("isClosestToTruth",&ctf.isClosestToTruth);
    candidateTree->Branch("truthCoordinates",&ctf.truthCoordinates);
  }
  // Coordinates
  candidateTree->Branch("geo_nuPositionX",&ctf.geo_nuPosX);
  candidateTree->Branch("geo_nuPositionY",&ctf.geo_nuPosY);
  candidateTree->Branch("geo_nuPositionZ",&ctf.geo_nuPosZ);
  candidateTree->Branch("geo_prongPositionX",&ctf.geo_prongPosX);
  candidateTree->Branch("geo_prongPositionY",&ctf.geo_prongPosY);
  candidateTree->Branch("geo_prongPositionZ",&ctf.geo_prongPosZ);
  candidateTree->Branch("geo_prongStartPositionX",&ctf.geo_prongStartPosX);
  candidateTree->Branch("geo_prongStartPositionY",&ctf.geo_prongStartPosY);
  candidateTree->Branch("geo_prongStartPositionZ",&ctf.geo_prongStartPosZ);
  candidateTree->Branch("geo_prongEndPositionX",&ctf.geo_prongEndPosX);
  candidateTree->Branch("geo_prongEndPositionY",&ctf.geo_prongEndPosY);
  candidateTree->Branch("geo_prongEndPositionZ",&ctf.geo_prongEndPosZ);
  candidateTree->Branch("geo_prongLength",&ctf.geo_prongLength);
  candidateTree->Branch("geo_openingAngle",&ctf.geo_openingAngle);
  // Direction
  candidateTree->Branch("geo_prongDirectionX",&ctf.geo_prongDirX);
  candidateTree->Branch("geo_prongDirectionY",&ctf.geo_prongDirY);
  candidateTree->Branch("geo_prongDirectionZ",&ctf.geo_prongDirZ);
  candidateTree->Branch("geo_prongTheta",&ctf.geo_prongTheta);
  candidateTree->Branch("geo_prongPhi",&ctf.geo_prongPhi);
  // Hypothesis info
  candidateTree->Branch("hypo_prongPdgCode_h1",&ctf.hypo_prongPdgCode_h1);
  candidateTree->Branch("hypo_prongPdgCode_h2",&ctf.hypo_prongPdgCode_h2);
  candidateTree->Branch("hypo_prongMass_h1",&ctf.hypo_prongMass_h1);
  candidateTree->Branch("hypo_prongMass_h2",&ctf.hypo_prongMass_h2);
  // Prong momentum (by range, assuming h1)
  candidateTree->Branch("range_prongEnergy_h1",&ctf.range_prongEnergy_h1);
  candidateTree->Branch("range_prongMomMag_h1",&ctf.range_prongMomMag_h1);
  candidateTree->Branch("range_prongMom_h1_X",&ctf.range_prongMom_h1_X);
  candidateTree->Branch("range_prongMom_h1_Y",&ctf.range_prongMom_h1_Y);
  candidateTree->Branch("range_prongMom_h1_Z",&ctf.range_prongMom_h1_Z);
  // Tot momentum (by range, assuming h1)
  candidateTree->Branch("range_invariantMass_h1",&ctf.range_invariantMass_h1);
  candidateTree->Branch("range_totEnergy_h1",&ctf.range_totEnergy_h1);
  candidateTree->Branch("range_totMomMag_h1",&ctf.range_totMomMag_h1);
  candidateTree->Branch("range_totMom_h1_X",&ctf.range_totMom_h1_X);
  candidateTree->Branch("range_totMom_h1_Y",&ctf.range_totMom_h1_Y);
  candidateTree->Branch("range_totMom_h1_Z",&ctf.range_totMom_h1_Z);
  // Tot momentum direction (by range, assuming h1)
  candidateTree->Branch("range_totDirection_h1_X",&ctf.range_totDir_h1_X);
  candidateTree->Branch("range_totDirection_h1_Y",&ctf.range_totDir_h1_Y);
  candidateTree->Branch("range_totDirection_h1_Z",&ctf.range_totDir_h1_Z);
  candidateTree->Branch("range_totTheta_h1",&ctf.range_totTheta_h1);
  candidateTree->Branch("range_totPhi_h1",&ctf.range_totPhi_h1);
  // Prong momentum (by range, assuming h2)
  candidateTree->Branch("range_prongEnergy_h2",&ctf.range_prongEnergy_h2);
  candidateTree->Branch("range_prongMomMag_h2",&ctf.range_prongMomMag_h2);
  candidateTree->Branch("range_prongMom_h2_X",&ctf.range_prongMom_h2_X);
  candidateTree->Branch("range_prongMom_h2_Y",&ctf.range_prongMom_h2_Y);
  candidateTree->Branch("range_prongMom_h2_Z",&ctf.range_prongMom_h2_Z);
  // Tot momentum (by range, assuming h2)
  candidateTree->Branch("range_invariantMass_h2",&ctf.range_invariantMass_h2);
  candidateTree->Branch("range_totEnergy_h2",&ctf.range_totEnergy_h2);
  candidateTree->Branch("range_totMomMag_h2",&ctf.range_totMomMag_h2);
  candidateTree->Branch("range_totMom_h2_X",&ctf.range_totMom_h2_X);
  candidateTree->Branch("range_totMom_h2_Y",&ctf.range_totMom_h2_Y);
  candidateTree->Branch("range_totMom_h2_Z",&ctf.range_totMom_h2_Z);
  // Tot momentum direction (by range, assuming h2)
  candidateTree->Branch("range_totDirection_h2_X",&ctf.range_totDir_h2_X);
  candidateTree->Branch("range_totDirection_h2_Y",&ctf.range_totDir_h2_Y);
  candidateTree->Branch("range_totDirection_h2_Z",&ctf.range_totDir_h2_Z);
  candidateTree->Branch("range_totTheta_h2",&ctf.range_totTheta_h2);
  candidateTree->Branch("range_totPhi_h2",&ctf.range_totPhi_h2);
  // Momentum (By Mcs)
  candidateTree->Branch("mcs_prongPdgCodeHypothesis",&ctf.mcs_prongPdgCodeHypothesis);
  candidateTree->Branch("mcs_prongIsBestFwd",&ctf.mcs_prongIsBestFwd);
  // Prong Momentum (By Mcs, best)
  candidateTree->Branch("mcs_prongMomMag_best_h1",&ctf.mcs_prongMomMag_best_h1);
  candidateTree->Branch("mcs_prongEnergy_best_h1",&ctf.mcs_prongEnergy_best_h1);
  candidateTree->Branch("mcs_prongMom_best_h1_X",&ctf.mcs_prongMom_best_h1_X);
  candidateTree->Branch("mcs_prongMom_best_h1_Y",&ctf.mcs_prongMom_best_h1_Y);
  candidateTree->Branch("mcs_prongMom_best_h1_Z",&ctf.mcs_prongMom_best_h1_Z);
  candidateTree->Branch("mcs_prongMomMag_best_h2",&ctf.mcs_prongMomMag_best_h2);
  candidateTree->Branch("mcs_prongEnergy_best_h2",&ctf.mcs_prongEnergy_best_h2);
  candidateTree->Branch("mcs_prongMom_best_h2_X",&ctf.mcs_prongMom_best_h2_X);
  candidateTree->Branch("mcs_prongMom_best_h2_Y",&ctf.mcs_prongMom_best_h2_Y);
  candidateTree->Branch("mcs_prongMom_best_h2_Z",&ctf.mcs_prongMom_best_h2_Z);
  // Tot momentum (by range, assuming both muons, best)
  candidateTree->Branch("mcs_totMomMag_best_h1",&ctf.mcs_totMomMag_best_h1);
  candidateTree->Branch("mcs_totEnergy_best_h1",&ctf.mcs_totEnergy_best_h1);
  candidateTree->Branch("mcs_invariantMass_best_h1",&ctf.mcs_invariantMass_best_h1);
  candidateTree->Branch("mcs_totMom_best_h1_X",&ctf.mcs_totMom_best_h1_X);
  candidateTree->Branch("mcs_totMom_best_h1_Y",&ctf.mcs_totMom_best_h1_Y);
  candidateTree->Branch("mcs_totMom_best_h1_Z",&ctf.mcs_totMom_best_h1_Z);
  candidateTree->Branch("mcs_totMomMag_best_h2",&ctf.mcs_totMomMag_best_h2);
  candidateTree->Branch("mcs_totEnergy_best_h2",&ctf.mcs_totEnergy_best_h2);
  candidateTree->Branch("mcs_invariantMass_best_h2",&ctf.mcs_invariantMass_best_h2);
  candidateTree->Branch("mcs_totMom_best_h2_X",&ctf.mcs_totMom_best_h2_X);
  candidateTree->Branch("mcs_totMom_best_h2_Y",&ctf.mcs_totMom_best_h2_Y);
  candidateTree->Branch("mcs_totMom_best_h2_Z",&ctf.mcs_totMom_best_h2_Z);
  // Tot momentum direction (by range, assuming both muons, best)
  candidateTree->Branch("mcs_totTheta_best_h1",&ctf.mcs_totTheta_best_h1);
  candidateTree->Branch("mcs_totPhi_best_h1",&ctf.mcs_totPhi_best_h1);
  candidateTree->Branch("mcs_totDir_best_h1_X",&ctf.mcs_totDir_best_h1_X);
  candidateTree->Branch("mcs_totDir_best_h1_Y",&ctf.mcs_totDir_best_h1_Y);
  candidateTree->Branch("mcs_totDir_best_h1_Z",&ctf.mcs_totDir_best_h1_Z);
  candidateTree->Branch("mcs_totTheta_best_h2",&ctf.mcs_totTheta_best_h2);
  candidateTree->Branch("mcs_totPhi_best_h2",&ctf.mcs_totPhi_best_h2);
  candidateTree->Branch("mcs_totDir_best_h2_X",&ctf.mcs_totDir_best_h2_X);
  candidateTree->Branch("mcs_totDir_best_h2_Y",&ctf.mcs_totDir_best_h2_Y);
  candidateTree->Branch("mcs_totDir_best_h2_Z",&ctf.mcs_totDir_best_h2_Z);
  // Others
  candidateTree->Branch("prongStartToNeutrinoDistance",&ctf.prongStartToNeutrinoDistance);
  candidateTree->Branch("prongNumHits",&ctf.prongNumHits);
  candidateTree->Branch("maxEndPointX",&ctf.maxEndPointX);
  candidateTree->Branch("maxEndPointY",&ctf.maxEndPointY);
  candidateTree->Branch("maxEndPointZ",&ctf.maxEndPointZ);
  candidateTree->Branch("deltaPhi",&ctf.deltaPhi);
  candidateTree->Branch("deltaTheta",&ctf.deltaTheta);
  candidateTree->Branch("lengthDiff",&ctf.lengthDiff);
  candidateTree->Branch("lengthRatio",&ctf.lengthRatio);
  candidateTree->Branch("maxStartToNeutrinoDistance",&ctf.maxStartToNeutrinoDistance);
  // // Calorimetry
  // candidateTree->Branch("calo_totChargeInRadius",&ctf.calo_totChargeInRadius);
  // candidateTree->Branch("calo_prong1ChargeInRadius",&ctf.calo_prong1ChargeInRadius);
  // candidateTree->Branch("calo_prong2ChargeInRadius",&ctf.calo_prong2ChargeInRadius);
  // candidateTree->Branch("calo_caloRatio",&ctf.calo_caloRatio);
  // // Status
  // candidateTree->Branch("status_nuWithMissingAssociatedVertex",&ctf.status_nuWithMissingAssociatedVertex);
  // candidateTree->Branch("status_nuWithMissingAssociatedTrack",&ctf.status_nuWithMissingAssociatedTrack);
  // candidateTree->Branch("status_nuProngWithMissingAssociatedHits",&ctf.status_nuProngWithMissingAssociatedHits);


  if (fSaveDrawTree)
  {
    drawTree = tfs->make<TTree>("DrawData","");
    drawTree->Branch("run",&dtf.run);
    drawTree->Branch("subrun",&dtf.subrun);
    drawTree->Branch("event",&dtf.event);
    drawTree->Branch("hsnID",&dtf.hsnID);
    drawTree->Branch("nHsnCandidatesInSameEvent",&dtf.nHsnCandidatesInSameEvent);
    drawTree->Branch("dv_p0_wireCoordinates",&dtf.dv_p0_wireCoordinates);
    drawTree->Branch("dv_p0_tickCoordinates",&dtf.dv_p0_tickCoordinates);
    drawTree->Branch("dv_p1_wireCoordinates",&dtf.dv_p1_wireCoordinates);
    drawTree->Branch("dv_p1_tickCoordinates",&dtf.dv_p1_tickCoordinates);
    drawTree->Branch("dv_p2_wireCoordinates",&dtf.dv_p2_wireCoordinates);
    drawTree->Branch("dv_p2_tickCoordinates",&dtf.dv_p2_tickCoordinates);
    drawTree->Branch("prong1_p0_wireCoordinates",&dtf.prong1_p0_wireCoordinates);
    drawTree->Branch("prong1_p0_tickCoordinates",&dtf.prong1_p0_tickCoordinates);
    drawTree->Branch("prong1_p1_wireCoordinates",&dtf.prong1_p1_wireCoordinates);
    drawTree->Branch("prong1_p1_tickCoordinates",&dtf.prong1_p1_tickCoordinates);
    drawTree->Branch("prong1_p2_wireCoordinates",&dtf.prong1_p2_wireCoordinates);
    drawTree->Branch("prong1_p2_tickCoordinates",&dtf.prong1_p2_tickCoordinates);
    drawTree->Branch("prong2_p0_wireCoordinates",&dtf.prong2_p0_wireCoordinates);
    drawTree->Branch("prong2_p0_tickCoordinates",&dtf.prong2_p0_tickCoordinates);
    drawTree->Branch("prong2_p1_wireCoordinates",&dtf.prong2_p1_wireCoordinates);
    drawTree->Branch("prong2_p1_tickCoordinates",&dtf.prong2_p1_tickCoordinates);
    drawTree->Branch("prong2_p2_wireCoordinates",&dtf.prong2_p2_wireCoordinates);
    drawTree->Branch("prong2_p2_tickCoordinates",&dtf.prong2_p2_tickCoordinates);
    drawTree->Branch("prong1_hits_p0_wireCoordinates",&dtf.prong1_hits_p0_wireCoordinates);
    drawTree->Branch("prong1_hits_p0_tickCoordinates",&dtf.prong1_hits_p0_tickCoordinates);
    drawTree->Branch("prong1_hits_p1_wireCoordinates",&dtf.prong1_hits_p1_wireCoordinates);
    drawTree->Branch("prong1_hits_p1_tickCoordinates",&dtf.prong1_hits_p1_tickCoordinates);
    drawTree->Branch("prong1_hits_p2_wireCoordinates",&dtf.prong1_hits_p2_wireCoordinates);
    drawTree->Branch("prong1_hits_p2_tickCoordinates",&dtf.prong1_hits_p2_tickCoordinates);
    drawTree->Branch("prong2_hits_p0_wireCoordinates",&dtf.prong2_hits_p0_wireCoordinates);
    drawTree->Branch("prong2_hits_p0_tickCoordinates",&dtf.prong2_hits_p0_tickCoordinates);
    drawTree->Branch("prong2_hits_p1_wireCoordinates",&dtf.prong2_hits_p1_wireCoordinates);
    drawTree->Branch("prong2_hits_p1_tickCoordinates",&dtf.prong2_hits_p1_tickCoordinates);
    drawTree->Branch("prong2_hits_p2_wireCoordinates",&dtf.prong2_hits_p2_wireCoordinates);
    drawTree->Branch("prong2_hits_p2_tickCoordinates",&dtf.prong2_hits_p2_tickCoordinates);
    drawTree->Branch("tot_hits_p0_wireCoordinates",&dtf.tot_hits_p0_wireCoordinates);
    drawTree->Branch("tot_hits_p0_tickCoordinates",&dtf.tot_hits_p0_tickCoordinates);
    drawTree->Branch("tot_hits_p1_wireCoordinates",&dtf.tot_hits_p1_wireCoordinates);
    drawTree->Branch("tot_hits_p1_tickCoordinates",&dtf.tot_hits_p1_tickCoordinates);
    drawTree->Branch("tot_hits_p2_wireCoordinates",&dtf.tot_hits_p2_wireCoordinates);
    drawTree->Branch("tot_hits_p2_tickCoordinates",&dtf.tot_hits_p2_tickCoordinates);
    if (fSaveTruthDrawTree)
    {
      drawTree->Branch("truth_dv_p0_wireCoordinates",&dtf.truth_dv_p0_wireCoordinates);
      drawTree->Branch("truth_dv_p0_tickCoordinates",&dtf.truth_dv_p0_tickCoordinates);
      drawTree->Branch("truth_dv_p1_wireCoordinates",&dtf.truth_dv_p1_wireCoordinates);
      drawTree->Branch("truth_dv_p1_tickCoordinates",&dtf.truth_dv_p1_tickCoordinates);
      drawTree->Branch("truth_dv_p2_wireCoordinates",&dtf.truth_dv_p2_wireCoordinates);
      drawTree->Branch("truth_dv_p2_tickCoordinates",&dtf.truth_dv_p2_tickCoordinates);
      drawTree->Branch("truth_nPrimaryTracks",&dtf.truth_nPrimaryTracks);
      drawTree->Branch("truth_nSecondaryTracks",&dtf.truth_nSecondaryTracks);
      drawTree->Branch("truth_primaryTracks_start_p0_wireCoordinates",&dtf.truth_primaryTracks_start_p0_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_start_p0_tickCoordinates",&dtf.truth_primaryTracks_start_p0_tickCoordinates);
      drawTree->Branch("truth_primaryTracks_start_p1_wireCoordinates",&dtf.truth_primaryTracks_start_p1_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_start_p1_tickCoordinates",&dtf.truth_primaryTracks_start_p1_tickCoordinates);
      drawTree->Branch("truth_primaryTracks_start_p2_wireCoordinates",&dtf.truth_primaryTracks_start_p2_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_start_p2_tickCoordinates",&dtf.truth_primaryTracks_start_p2_tickCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p0_wireCoordinates",&dtf.truth_primaryTracks_end_p0_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p0_tickCoordinates",&dtf.truth_primaryTracks_end_p0_tickCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p1_wireCoordinates",&dtf.truth_primaryTracks_end_p1_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p1_tickCoordinates",&dtf.truth_primaryTracks_end_p1_tickCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p2_wireCoordinates",&dtf.truth_primaryTracks_end_p2_wireCoordinates);
      drawTree->Branch("truth_primaryTracks_end_p2_tickCoordinates",&dtf.truth_primaryTracks_end_p2_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p0_wireCoordinates",&dtf.truth_secondaryTracks_start_p0_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p0_tickCoordinates",&dtf.truth_secondaryTracks_start_p0_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p1_wireCoordinates",&dtf.truth_secondaryTracks_start_p1_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p1_tickCoordinates",&dtf.truth_secondaryTracks_start_p1_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p2_wireCoordinates",&dtf.truth_secondaryTracks_start_p2_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_start_p2_tickCoordinates",&dtf.truth_secondaryTracks_start_p2_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p0_wireCoordinates",&dtf.truth_secondaryTracks_end_p0_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p0_tickCoordinates",&dtf.truth_secondaryTracks_end_p0_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p1_wireCoordinates",&dtf.truth_secondaryTracks_end_p1_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p1_tickCoordinates",&dtf.truth_secondaryTracks_end_p1_tickCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p2_wireCoordinates",&dtf.truth_secondaryTracks_end_p2_wireCoordinates);
      drawTree->Branch("truth_secondaryTracks_end_p2_tickCoordinates",&dtf.truth_secondaryTracks_end_p2_tickCoordinates);
      drawTree->Branch("truth_nPrimaryShowers",&dtf.truth_nPrimaryShowers);
      drawTree->Branch("truth_nSecondaryShowers",&dtf.truth_nSecondaryShowers);
      drawTree->Branch("truth_primaryShowers_start_p0_wireCoordinates",&dtf.truth_primaryShowers_start_p0_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_start_p0_tickCoordinates",&dtf.truth_primaryShowers_start_p0_tickCoordinates);
      drawTree->Branch("truth_primaryShowers_start_p1_wireCoordinates",&dtf.truth_primaryShowers_start_p1_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_start_p1_tickCoordinates",&dtf.truth_primaryShowers_start_p1_tickCoordinates);
      drawTree->Branch("truth_primaryShowers_start_p2_wireCoordinates",&dtf.truth_primaryShowers_start_p2_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_start_p2_tickCoordinates",&dtf.truth_primaryShowers_start_p2_tickCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p0_wireCoordinates",&dtf.truth_primaryShowers_end_p0_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p0_tickCoordinates",&dtf.truth_primaryShowers_end_p0_tickCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p1_wireCoordinates",&dtf.truth_primaryShowers_end_p1_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p1_tickCoordinates",&dtf.truth_primaryShowers_end_p1_tickCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p2_wireCoordinates",&dtf.truth_primaryShowers_end_p2_wireCoordinates);
      drawTree->Branch("truth_primaryShowers_end_p2_tickCoordinates",&dtf.truth_primaryShowers_end_p2_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p0_wireCoordinates",&dtf.truth_secondaryShowers_start_p0_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p0_tickCoordinates",&dtf.truth_secondaryShowers_start_p0_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p1_wireCoordinates",&dtf.truth_secondaryShowers_start_p1_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p1_tickCoordinates",&dtf.truth_secondaryShowers_start_p1_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p2_wireCoordinates",&dtf.truth_secondaryShowers_start_p2_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_start_p2_tickCoordinates",&dtf.truth_secondaryShowers_start_p2_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p0_wireCoordinates",&dtf.truth_secondaryShowers_end_p0_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p0_tickCoordinates",&dtf.truth_secondaryShowers_end_p0_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p1_wireCoordinates",&dtf.truth_secondaryShowers_end_p1_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p1_tickCoordinates",&dtf.truth_secondaryShowers_end_p1_tickCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p2_wireCoordinates",&dtf.truth_secondaryShowers_end_p2_wireCoordinates);
      drawTree->Branch("truth_secondaryShowers_end_p2_tickCoordinates",&dtf.truth_secondaryShowers_end_p2_tickCoordinates);
    }
  }

} // END function beginJob

void HsnFinder::endJob()
{} // END function endJob

void HsnFinder::ClearData()
{} // END function ClearData


// Core analysis. This is where all the functions are executed. Gets repeated event by event.
void HsnFinder::analyze(art::Event const & evt)
{
  if (fVerbose) {printf("\n\n\n---------------------------------------------------\n");}
  if (fVerbose) {printf("||HSN FINDER MODULE: EVENT %i [RUN %i, SUBRUN %i]||\n", evt.id().event(), evt.id().subRun(), evt.id().run());}

  // Determine event ID and initialize event tree filler.
  // The event tree filler is a special class in which we fill all the information we want to know about the current event.
  // At the end of the event loop the information is taken from the event tree filler and filled into the anatree.
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();
  etf.Initialize(run,subrun,event);

  // Search among pfparticles and get vector of potential neutrino pfps with only two tracks. Return vectors of pfps for neutrinos, tracks and showers in event and decay vertices, which contain information about neutrino vertices with exctly two tracks.
  std::vector<AuxVertex::DecayVertex> ana_decayVertices;
  fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evt, etf, ana_decayVertices);
  etf.nHsnCandidates = ana_decayVertices.size();

  // Now, IF there are any candidates, go on. Otherwise you can stop here
  if (ana_decayVertices.size() == 0)
  {
    printf("No clean vertex candidates found. Moving to next event...\n");
    eventTree->Fill();
    return;
  }
  else
  {
    // IF there are candidate, continue with analysis
    // Perform calorimetry analysis (TO BE REVIEWED)
    // fCalorimetryRadiusAlg.PerformCalorimetry(evt, etf, ana_decayVertices);

    // If want to use reco-truth distance as a metric for finding best HSN candidate, do it here.
    if ( fUseTruthDistanceMetric ) { fExtractTruthInformationAlg.FillEventTreeWithTruth(evt,etf,ana_decayVertices);}

    // Now loop for each candidate and fill the tree
    for (std::vector<int>::size_type i=0; i!=ana_decayVertices.size(); i++)
    {
      // The candidate tree filler is a special class in which we fill all the information we want to know about the current HSN candidate.
      // It is filled multiple times in each event.
      ctf.Initialize(etf,i,ana_decayVertices[i],fCenterCoordinates);
      // Fill tree
      candidateTree->Fill();

      // If requested, do the same for the draw tree
      if (fSaveDrawTree)
      {
        dtf.Initialize(etf,i,ana_decayVertices[i]);
        if (fSaveTruthDrawTree)
        {
          fExtractTruthInformationAlg.FillDrawTreeWithTruth(evt,dtf);
        }
        drawTree->Fill();
      }
    } // END FOR loop for each candidate
    eventTree->Fill();
  } // END IF there are any candidates
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(HsnFinder)

#endif // END def HsnFinder_module
