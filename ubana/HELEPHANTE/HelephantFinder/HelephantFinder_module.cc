#ifndef HELEPHANTFINDER_MODULE
#define HELEPHANTFINDER_MODULE

#include "HelephantFinder.h"

HelephantFinder::HelephantFinder(fhicl::ParameterSet const & pset) :
    // Algorithms
    EDAnalyzer(pset),
    fFindPandoraVertexAlg(pset),
    fCalorimetryRadiusAlg(pset),
    fExtractTruthInformationAlg(pset),
    fFlashMatchingAlg(pset),
    fTriggerInformationAlg(pset),
    fG4InformationAlg(pset),
    fdEdXInformationAlg(pset),
    // Meta information
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fRadiusProfileLimits(pset.get<std::vector<double>>("RadiusProfileLimits")),
    fRadiusProfileBins(pset.get<int>("RadiusProfileBins")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fVerbose(pset.get<bool>("VerboseMode")),
    // Geometry settings
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fCenterCoordinates(pset.get<std::vector<double>>("CenterCoordinates")),
    // Data product labels settings
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fMcsLabel(pset.get<std::string>("McsLabel")),
    fMcTrackLabel(pset.get<std::string>("McTrackLabel")),
    fFlashLabel(pset.get<std::string>("FlashLabel")),
    fCalorimetryLabel(pset.get<std::string>("CalorimetryLabel")),
    fCalibratedLabel(pset.get<std::string>("CalibratedLabel")),
    // Additional tree information
    fSaveDrawTree(pset.get<bool>("SaveDrawTree")),
    fSaveTruthDrawTree(pset.get<bool>("SaveTruthDrawTree")),
    fSaveTriggerInformation(pset.get<bool>("SaveTriggerInformation")),
    fSaveG4Information(pset.get<bool>("SaveG4Information")),
    fSavedEdXInformation(pset.get<bool>("SavedEdXInformation")),
    // Data type dependent settings
    fUseTruthDistanceMetric(pset.get<bool>("UseTruthDistanceMetric")),
    fIsHSN(pset.get<bool>("IsHSN")),
    // Additional settings (finer tuning)
    fSaveAlsoUncalibrateddEdXInformation(pset.get<bool>("SaveAlsoUncalibrateddEdXInformation"))
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
} // END constructor HelephantFinder

HelephantFinder::~HelephantFinder()
{} // END destructor HelephantFinder

void HelephantFinder::beginJob()
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

  // Tree containing data at the event level
  // (each row is for a single event)
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

  // Tree containing data at the candidate level
  // (each row is for a possible HSN candidate, there can be multiple rows corresponding to the same event!)
  candidateTree = tfs->make<TTree>("CandidateData","");
  // HSN ID (used to put in correlation with event tree)
  candidateTree->Branch("run",&ctf.run);
  candidateTree->Branch("subrun",&ctf.subrun);
  candidateTree->Branch("event",&ctf.event);
  candidateTree->Branch("hsnID",&ctf.hsnID);
  candidateTree->Branch("nHsnCandidatesInSameEvent",&ctf.nHsnCandidatesInSameEvent);
  // Reco-truth matching information (only for MC)
  if ( fUseTruthDistanceMetric )
  {
    candidateTree->Branch("recoTruthDistance",&ctf.recoTruthDistance);
    candidateTree->Branch("isClosestToTruth",&ctf.isClosestToTruth);
    candidateTree->Branch("truthCoordinates",&ctf.truthCoordinates);
  }
  // Candidate coordinates
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
  // Candidate direction
  candidateTree->Branch("geo_prongDirectionX",&ctf.geo_prongDirX);
  candidateTree->Branch("geo_prongDirectionY",&ctf.geo_prongDirY);
  candidateTree->Branch("geo_prongDirectionZ",&ctf.geo_prongDirZ);
  candidateTree->Branch("geo_prongTheta",&ctf.geo_prongTheta);
  candidateTree->Branch("geo_prongPhi",&ctf.geo_prongPhi);
  // Candidate hypothesis (which pdg and mass we are using for each h)
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
  // Flash
  candidateTree->Branch("flash_flashDistance",&ctf.flash_flashDistance);
  candidateTree->Branch("flash_flashPE",&ctf.flash_flashPE);
  // Trigger
  if (fSaveTriggerInformation)
  {
    candidateTree->Branch("trigger_passedBNB",&ctf.trigger_passedBNB);
    candidateTree->Branch("trigger_passedBNB5PE",&ctf.trigger_passedBNB5PE);
    candidateTree->Branch("trigger_passedEXT",&ctf.trigger_passedEXT);
    candidateTree->Branch("trigger_passedEXT5PE",&ctf.trigger_passedEXT5PE);
    candidateTree->Branch("trigger_passedHSN",&ctf.trigger_passedHSN);
    candidateTree->Branch("trigger_passedHSNEXT",&ctf.trigger_passedHSNEXT);
  }
  // G4 Information
  if (fSaveG4Information)
  {
    candidateTree->Branch("g4_pi_pdg",&ctf.g4_pi_pdg);
    candidateTree->Branch("g4_mu_pdg",&ctf.g4_mu_pdg);
    candidateTree->Branch("g4_pi_mom",&ctf.g4_pi_mom);
    candidateTree->Branch("g4_mu_mom",&ctf.g4_mu_mom);
    candidateTree->Branch("g4_pi_endMom",&ctf.g4_pi_endMom);
    candidateTree->Branch("g4_mu_endMom",&ctf.g4_mu_endMom);
    candidateTree->Branch("g4_reco_pi_momRange",&ctf.g4_reco_pi_momRange);
    candidateTree->Branch("g4_reco_mu_momRange",&ctf.g4_reco_mu_momRange);
    candidateTree->Branch("g4_reco_pi_momMCS",&ctf.g4_reco_pi_momMCS);
    candidateTree->Branch("g4_reco_mu_momMCS",&ctf.g4_reco_mu_momMCS);
    candidateTree->Branch("g4_pi_length",&ctf.g4_pi_length);
    candidateTree->Branch("g4_mu_length",&ctf.g4_mu_length);
    candidateTree->Branch("g4_pi_nScatterings",&ctf.g4_pi_nScatterings);
    candidateTree->Branch("g4_mu_nScatterings",&ctf.g4_mu_nScatterings);
    candidateTree->Branch("g4_pi_declaredEndProcess",&ctf.g4_pi_declaredEndProcess);
    candidateTree->Branch("g4_mu_declaredEndProcess",&ctf.g4_mu_declaredEndProcess);
    candidateTree->Branch("g4_pi_declaredScatteredEndProcess",&ctf.g4_pi_declaredScatteredEndProcess);
    candidateTree->Branch("g4_mu_declaredScatteredEndProcess",&ctf.g4_mu_declaredScatteredEndProcess);
    candidateTree->Branch("g4_pi_endProcess",&ctf.g4_pi_endProcess);
    candidateTree->Branch("g4_mu_endProcess",&ctf.g4_mu_endProcess);
    candidateTree->Branch("g4_recoTruth_distancesScore",&ctf.g4_recoTruth_distancesScore);
  }

  if (fSavedEdXInformation)
  {
    candidateTree->Branch("cali_prongCaloPlane",&ctf.cali_prongCaloPlane);
  }

  // Separate tree for custom event display
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
    // Additional information to draw also MCTruth in the custom event display
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
    } // END if save truth in draw tree
  } //  END if save draw tree
} // END function beginJob

void HelephantFinder::endJob()
{} // END function endJob

void HelephantFinder::ClearData()
{} // END function ClearData


// Core analysis. This is where all the functions are executed. Gets repeated event by event.
void HelephantFinder::analyze(art::Event const & evt)
{
  // Verbose information
  if (fVerbose) {printf("\n\n\n---------------------------------------------------\n");}
  if (fVerbose) {printf("||HSN FINDER MODULE: EVENT %i [RUN %i, SUBRUN %i]||\n", evt.id().event(), evt.id().subRun(), evt.id().run());}

  // Determine event ID and initialize event tree filler (etf).
  // The event tree filler (etf) is a special class in which we store temporarily all the information we want to know about the current event, at the event level (not for each HSN candidate). At the end of the event loop the information from the event tree filler is taken and saved into the anatree.
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();
  etf.Initialize(run,subrun,event);

  // Bulk of the module, doing the heavy lifting. Search among all the reconstructed pfparticles and find a HSN candidate: A neutrino pfparticle that has two and only two tracks. Then store all this information in a special class, AuxVertex::DecayVertex, which contains pointers to reco tracks, hits and reco vertex representing the current HSN candidate.
  std::vector<AuxVertex::DecayVertex> ana_decayVertices;
  fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evt, etf, ana_decayVertices);
  etf.nHsnCandidates = ana_decayVertices.size();

  // Once all the possible candidates in this event have been found (there can be more than one), continue with the analysis, obtaining additional information. Of course, if there's not candidate in this event you can just move on to the next. Now, IF there are candidates, keep performing analysis.
  if (ana_decayVertices.size() != 0)
  {
    // Before we start dealing with each candidate one by one, we need to consider information at the event level (which might contain multiple candidates).
    // Run the reco-truth distance metric here, at the event level. Why? The next algorithm picks up the best candidate based on the shortest reco-truth distance metric. Since it needs to compare candidates with each other this needs to be done at the event level (the for loop below is for each individual candidates, which doesn't know anything about the other candidates).
    if ( fUseTruthDistanceMetric ) { fExtractTruthInformationAlg.FillEventTreeWithTruth(evt,etf,ana_decayVertices);}

    // Finally, let's loop for each candidate and run analysis algorithms on them.
    for (std::vector<int>::size_type i=0; i!=ana_decayVertices.size(); i++)
    {
      AuxVertex::DecayVertex currentVertex = ana_decayVertices[i];
      // First of all, we initialize the candidate tree filler which is a special class in which we fill all the information we want to know about the current HSN candidate.
      // It is filled multiple times in each event.
      ctf.Initialize(etf,i,currentVertex,fCenterCoordinates);

      // Run flash matching for each candidate
      fFlashMatchingAlg.AddFlashMatching(evt,ctf,currentVertex);

      // If requested, save also trigger information regarding the event the candidates is in.
      // (technically this should be at event level, but it ends up being more useful here).
      if (fSaveTriggerInformation) fTriggerInformationAlg.AddTriggerInformation(evt,ctf);

      // If requested, save G4-level information regarding the truth of the event the current candidate is in.
      // (technically this should be at event level, but it ends up being more useful here).
      if (fSaveG4Information) fG4InformationAlg.AddG4Information(evt,ctf,currentVertex);

      // If requested, save information regarding dEdX for the current vertex
      if (fSavedEdXInformation) fdEdXInformationAlg.AdddEdXInformation(evt,ctf,currentVertex);

      // If requested, save information for the custom event display
      if (fSaveDrawTree)
      {
        // Initialize ALSO FILLS the tree (that's why there's no function)
        dtf.Initialize(etf,i,currentVertex);
        // If requested save also the truth information in the event display.
        // This is optional, because DATA shouldn't run this function (unlike MC, which can).
        if (fSaveTruthDrawTree) fExtractTruthInformationAlg.FillDrawTreeWithTruth(evt,dtf);
        drawTree->Fill();
      }
      
      // Fill candidate-level tree for each candidate
      candidateTree->Fill();
    } // END FOR loop for each candidate
  } // END IF there are any candidates

  // IF there aren't any candidates, print some information and move on.
  else printf("No clean vertex candidates found. Moving to next event...\n");

  // Fill event-level tree for each event
  eventTree->Fill();
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(HelephantFinder)

#endif // END def HelephantFinder_module
