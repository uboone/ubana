#ifndef EXPLORETRUTH_MODULE
#define EXPLORETRUTH_MODULE

#include "SingleParticleMomentumStudy.h"

SingleParticleMomentumStudy::SingleParticleMomentumStudy(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fDetailed_dEdX(pset.get<bool>("detailed_dEdX"))
{} // END constructor SingleParticleMomentumStudy

SingleParticleMomentumStudy::~SingleParticleMomentumStudy()
{} // END destructor SingleParticleMomentumStudy

void SingleParticleMomentumStudy::beginJob()
{
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
  tDataTree->Branch("truth_pdgCode",&truth_pdgCode);
  tDataTree->Branch("truth_nTracks",&truth_nTracks);
  tDataTree->Branch("truth_length",&truth_length);
  tDataTree->Branch("declaredEndProcess",&declaredEndProcess);
  tDataTree->Branch("daughterProcess",&daughterProcess);
  tDataTree->Branch("daughterPdg",&daughterPdg);
  tDataTree->Branch("daughterEnergy",&daughterEnergy);
  tDataTree->Branch("reco_nTracks",&reco_nTracks);
  tDataTree->Branch("reco_length",&reco_length);
  // Coordinates
  tDataTree->Branch("truth_start",&truth_start);
  tDataTree->Branch("truth_end",&truth_end);
  tDataTree->Branch("reco_start",&reco_start);
  tDataTree->Branch("reco_end",&reco_end); 
  // Truth momentum
  tDataTree->Branch("truth_momentum",&truth_momentum);
  // Range momentum
  tDataTree->Branch("rangeTruth_muon_momentum",&rangeTruth_muon_momentum);
  tDataTree->Branch("rangeTruth_pion_momentum",&rangeTruth_pion_momentum);
  tDataTree->Branch("range_muon_momentum",&range_muon_momentum);
  tDataTree->Branch("range_pion_momentum",&range_pion_momentum);
  tDataTree->Branch("range_pion_momentum2",&range_pion_momentum2);
  tDataTree->Branch("range_proton_momentum2",&range_pion_momentum2);
  // MCS momentum
  tDataTree->Branch("mcs_momentum_best",&mcs_momentum_best);
  tDataTree->Branch("mcs_momentum_forward",&mcs_momentum_forward);
  // MCS test
  tDataTree->Branch("mcs_momentum_mu1",&mcs_momentum_mu1); 
  tDataTree->Branch("mcs_momentum_mu2",&mcs_momentum_mu2); 
  tDataTree->Branch("mcs_momentum_mu3",&mcs_momentum_mu3); 
  tDataTree->Branch("mcs_momentum_mu4",&mcs_momentum_mu4); 
  tDataTree->Branch("mcs_momentum_pi1",&mcs_momentum_pi1); 
  tDataTree->Branch("mcs_momentum_pi2",&mcs_momentum_pi2); 
  tDataTree->Branch("mcs_momentum_pi3",&mcs_momentum_pi3); 
  tDataTree->Branch("mcs_momentum_pi4",&mcs_momentum_pi4); 
  // Calorimetry momentum
  tDataTree->Branch("cali_caloPlane",&cali_caloPlane); 
  tDataTree->Branch("cali_isTrackFlipped",&cali_isTrackFlipped); 
  tDataTree->Branch("cali_kinEnergy",&cali_kinEnergy); 
  tDataTree->Branch("cali_range",&cali_range);
  if (fDetailed_dEdX)
  {
    tDataTree->Branch("cali_dEdX",&cali_dEdX); 
    tDataTree->Branch("cali_dQdX",&cali_dQdX); 
    tDataTree->Branch("cali_resRange",&cali_resRange);
  }
} // END function beginJob

void SingleParticleMomentumStudy::endJob()
{
} // END function endJob



// FUNCTION
// Fill end process for McParticle
void SingleParticleMomentumStudy::FillEndProcess(
        art::Ptr<simb::MCParticle> const & mcp,
        std::map< int,art::Ptr<simb::MCParticle> > mcpMap)
{
  daughterProcess.clear();
  daughterPdg.clear();
  daughterEnergy.clear();
  declaredEndProcess = mcp->EndProcess();

  for(int j=0; j!=mcp->NumberDaughters(); j++)
  {
    int daughterId = mcp->Daughter(j);
    art::Ptr<simb::MCParticle> daughter;
    if (mcpMap.count(daughterId)) daughter = mcpMap[daughterId];
    else continue;

    std::string dProcess = daughter->Process();
    // There's too many electrons/photons coming from the ionization we don't care about
    if (dProcess!=std::string("hIoni") && dProcess!=std::string("muIoni"))
    {
      daughterProcess.push_back(dProcess);
      daughterPdg.push_back(daughter->PdgCode());
      daughterEnergy.push_back(daughter->E());
    }
  } // END loop through all daughters of current primary
} // END function FillEndProcess



// FUNCTION
// Build map of trackId - McParticle pointers
std::map< int,art::Ptr<simb::MCParticle> > SingleParticleMomentumStudy::BuildMcpMap(art::ValidHandle< std::vector<simb::MCParticle> > const & mcpHandle)
{
  std::map< int,art::Ptr<simb::MCParticle> > mcpMap;
  for(int i=0; i!=int((*mcpHandle).size()); i++)
  {
    art::Ptr<simb::MCParticle> map_mcp(mcpHandle,i);
    int trackId = map_mcp->TrackId();
    mcpMap[trackId] = map_mcp;
  }
  return mcpMap;
} // END function BuildMcpMap



// FUNCTION
// Build map of trackId - McTrack pointers
std::map< int,art::Ptr<sim::MCTrack> > SingleParticleMomentumStudy::BuildMctMap(art::ValidHandle< std::vector<sim::MCTrack> > const & mctHandle)
{
  std::map< int,art::Ptr<sim::MCTrack> > mctMap;
  for(int i=0; i!=int((*mctHandle).size()); i++)
  {
    art::Ptr<sim::MCTrack> map_mct(mctHandle,i);
    int trackId = map_mct->TrackID();
    mctMap[trackId] = map_mct;
  }
  return mctMap;
} // END function BuildMcpMap



// FUNCTION
// Associate the correct reco track to the primary
art::Ptr<recob::Track> SingleParticleMomentumStudy::FindCorrectRecoTrack(
        art::ValidHandle< std::vector<recob::Track> > const & recoHandle,
        art::Ptr<sim::MCTrack> const & primary_mct)
{
  art::Ptr<recob::Track> correct_recot(recoHandle,0);
  float minScore = 1e30;
  for(int i=0; i!=int((*recoHandle).size()); i++)
  {
    art::Ptr<recob::Track> recot(recoHandle,i);
    float recoTruthScore = RecoTruthScore(primary_mct, recot);
    if (recoTruthScore < minScore)
    {
      minScore = recoTruthScore;
      correct_recot = recot;
    }
  }
  return correct_recot;
} // END function FindCorrectRecoTrack



// FUNCTION
// Determine the "score" between truth and reco tracks based 
float SingleParticleMomentumStudy::RecoTruthScore(
        art::Ptr<sim::MCTrack> const & primary_mct,
        art::Ptr<recob::Track> const & reco_candidate)
{ 
  float x1, x2, y1, y2, z1, z2;
  float startScore, endScore, totScore_fwd, totScore_bwd, bestTotScore;
  // Forward score
  x1 = primary_mct->Start().X();
  y1 = primary_mct->Start().Y();
  z1 = primary_mct->Start().Z();
  x2 = reco_candidate->Start().X();
  y2 = reco_candidate->Start().Y();
  z2 = reco_candidate->Start().Z();
  startScore = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
  x1 = primary_mct->End().X();
  y1 = primary_mct->End().Y();
  z1 = primary_mct->End().Z();
  x2 = reco_candidate->End()[0];
  y2 = reco_candidate->End()[1];
  z2 = reco_candidate->End()[2];
  endScore = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
  totScore_fwd = sqrt(pow(startScore,2.) + pow(endScore,2.));

  // Backward score
  x1 = primary_mct->Start().X();
  y1 = primary_mct->Start().Y();
  z1 = primary_mct->Start().Z();
  x2 = reco_candidate->End()[0];
  y2 = reco_candidate->End()[1];
  z2 = reco_candidate->End()[2];
  startScore = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
  x1 = primary_mct->End().X();
  y1 = primary_mct->End().Y();
  z1 = primary_mct->End().Z();
  x2 = reco_candidate->Start().X();
  y2 = reco_candidate->Start().Y();
  z2 = reco_candidate->Start().Z();
  endScore = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
  totScore_bwd = sqrt(pow(startScore,2.) + pow(endScore,2.));

  bestTotScore = std::min(totScore_fwd, totScore_bwd);
  return bestTotScore;
} //  END function RecoTruthSCore



// FUNCTION
// Determine MCTrack length
float SingleParticleMomentumStudy::MCTrackLength(art::Ptr<sim::MCTrack> const & mct)
{
  float x1, x2, y1, y2, z1, z2;
  x1 = mct->Start().X();
  y1 = mct->Start().Y();
  z1 = mct->Start().Z();
  x2 = mct->End().X();
  y2 = mct->End().Y();
  z2 = mct->End().Z();
  float length = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));  
  return length;
} // END function MCTrackLength



// FUNCTION
// Determine Track length
float SingleParticleMomentumStudy::TrackLength(art::Ptr<recob::Track> const & recot)
{
  float x1, x2, y1, y2, z1, z2;
  x1 = recot->Start().X();
  y1 = recot->Start().Y();
  z1 = recot->Start().Z();
  x2 = recot->End()[0];
  y2 = recot->End()[1];
  z2 = recot->End()[2];
  float length = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
  // float length = recot->Length();
  return length;
} // END function TrackLength



// FUNCTION
// Determine whether track has been flipped by making sure that dqdx is larger at the end than at the start.
bool SingleParticleMomentumStudy::IsTrackFlipped(art::Ptr<anab::Calorimetry> const & calo)
{
  bool isTrackFlipped;
  if ((calo->dQdx()).size()>5)
  {
    float dqdx_start = (calo->dQdx())[0] +  (calo->dQdx())[1] +  (calo->dQdx())[2];
    float dqdx_end = (calo->dQdx())[-1] +  (calo->dQdx())[-2] +  (calo->dQdx())[-3];
    isTrackFlipped = dqdx_start < dqdx_end;
  }
  else isTrackFlipped = false;
  return isTrackFlipped;
} //  END function IsTrackFlipped



// MAIN FUNCTION ANALYZE
void SingleParticleMomentumStudy::analyze(art::Event const & evt)
{
  // Clear vectors
  truth_start.clear();
  truth_end.clear();
  reco_start.clear();
  reco_end.clear();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  // Get McTruth handle
  art::InputTag mcpTag {std::string("largeant")};
  art::ValidHandle< std::vector<simb::MCParticle> > const & mcpHandle = evt.getValidHandle< std::vector<simb::MCParticle> >(mcpTag);
  // Get McTrack handle
  art::InputTag mctTag {std::string("mcreco")};
  art::ValidHandle< std::vector<sim::MCTrack> > const & mctHandle = evt.getValidHandle< std::vector<sim::MCTrack> >(mctTag);
  // Get Reco handle
  art::InputTag recoTag {std::string("pandoraNu")};
  art::ValidHandle< std::vector<recob::Track> > const & recoHandle = evt.getValidHandle< std::vector<recob::Track> >(recoTag);
  // Get MCSFIT handle
  art::InputTag mcsTag {std::string("pandoraNuMCSMu")};
  art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag);
  // MCS test
  // art::InputTag mcsTag_mu1 {std::string("pandoraNuMCSMuShort1")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_mu1 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_mu1);
  // art::InputTag mcsTag_mu2 {std::string("pandoraNuMCSMuShort2")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_mu2 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_mu2);
  // art::InputTag mcsTag_mu3 {std::string("pandoraNuMCSMuShort3")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_mu3 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_mu3);
  // art::InputTag mcsTag_mu4 {std::string("pandoraNuMCSMuShort4")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_mu4 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_mu4);

  // art::InputTag mcsTag_pi1 {std::string("pandoraNuMCSPiShort1")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_pi1 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_pi1);
  // art::InputTag mcsTag_pi2 {std::string("pandoraNuMCSPiShort2")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_pi2 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_pi2);
  // art::InputTag mcsTag_pi3 {std::string("pandoraNuMCSPiShort3")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_pi3 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_pi3);
  // art::InputTag mcsTag_pi4 {std::string("pandoraNuMCSPiShort4")};
  // art::ValidHandle< std::vector<recob::MCSFitResult> > const & mcsHandle_pi4 = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag_pi4);

  if((*recoHandle).size()==0)
  {
    // printf("No tracks in current event. Skipping.");
    return;
  }

  // Build map from TrackId to McParticles and McTracks
  std::map< int,art::Ptr<simb::MCParticle> > mcpMap = BuildMcpMap(mcpHandle);
  std::map< int,art::Ptr<sim::MCTrack> > mctMap = BuildMctMap(mctHandle);
  // Find the only primary among McParticles
  std::vector<art::Ptr<simb::MCParticle>> primaries_mcp;
  for(int i=0; i!=int((*mcpHandle).size()); i++)
  {
    art::Ptr<simb::MCParticle> mcp(mcpHandle,i);
    if(mcp->Process() == std::string("primary")) primaries_mcp.push_back(mcp);
  }

  // Determine right MCParticle and MCTrack for only primary
  art::Ptr<simb::MCParticle> primary_mcp = primaries_mcp[0];
  art::Ptr<sim::MCTrack> primary_mct = mctMap[primary_mcp->TrackId()];
  // Find the right Track associated with the MCParticle and its MCSFitResult
  art::Ptr<recob::Track> primary_reco = FindCorrectRecoTrack(recoHandle, primary_mct);
  art::Ptr<recob::MCSFitResult> primary_mcs(mcsHandle, primary_reco.key());
  // Now fill all information relevant to the tree
  // General
  FillEndProcess(primary_mcp,mcpMap);
  truth_pdgCode = primary_mcp->PdgCode();
  truth_nTracks = int((*mctHandle).size());
  reco_nTracks = int((*recoHandle).size());
  truth_length = MCTrackLength(primary_mct);
  reco_length = TrackLength(primary_reco);
  // Coordinates
  float x1, x2, y1, y2, z1, z2;
  x1 = primary_mct->Start().X();
  y1 = primary_mct->Start().Y();
  z1 = primary_mct->Start().Z();
  x2 = primary_mct->End().X();
  y2 = primary_mct->End().Y();
  z2 = primary_mct->End().Z();
  truth_start = {x1,y1,z1};
  truth_end = {x2,y2,z2};
  x1 = primary_reco->Start().X();
  y1 = primary_reco->Start().Y();
  z1 = primary_reco->Start().Z();
  x2 = primary_reco->End()[0];
  y2 = primary_reco->End()[1];
  z2 = primary_reco->End()[2];
  reco_start = {x1,y1,z1};
  reco_end = {x2,y2,z2};
  // Truth momentum
  float truth_momentum_x = primary_mct->Start().Px();
  float truth_momentum_y = primary_mct->Start().Py();
  float truth_momentum_z = primary_mct->Start().Pz();
  truth_momentum = sqrt(pow(truth_momentum_x,2.) + pow(truth_momentum_y,2.) + pow(truth_momentum_z,2.));
  // printf("MCTrack momentum: [ %.2f ], [ %.2f, %.2f, %.2f ]\n", primary_mcp->P(), primary_mcp->Px(), primary_mcp->Py() , primary_mcp->Pz() );
  // printf("MCParticle momentum: [ %.2f ], [ %.2f, %.2f, %.2f ]\n", truth_momentum, truth_momentum_x, truth_momentum_y, truth_momentum_z);
  // Range momentum
  float muon_mass = 0.105658;
  float pion_mass = 0.139570;
  trkf::TrackMomentumCalculator tmc;
  rangeTruth_muon_momentum = tmc.GetTrackMomentum(truth_length,13);
  rangeTruth_pion_momentum = tmc.GetTrackMomentum(truth_length*(muon_mass/pion_mass),13)*(pion_mass/muon_mass);
  range_muon_momentum = tmc.GetTrackMomentum(reco_length,13);
  range_pion_momentum = tmc.GetTrackMomentum(reco_length*(muon_mass/pion_mass),13)*(pion_mass/muon_mass);
  // MCS momentum
  mcs_momentum_best = primary_mcs->bestMomentum();
  mcs_momentum_forward = primary_mcs->fwdMomentum();
  // printf("Track momentum: [ %.2f ]\n", rangeTruth_muon_momentum);
  // MCS test
  // art::Ptr<recob::MCSFitResult> primary_mcs_mu1(mcsHandle_mu1, primary_reco.key());
  // mcs_momentum_mu1 = primary_mcs_mu1->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_mu2(mcsHandle_mu2, primary_reco.key());
  // mcs_momentum_mu2 = primary_mcs_mu2->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_mu3(mcsHandle_mu3, primary_reco.key());
  // mcs_momentum_mu3 = primary_mcs_mu3->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_mu4(mcsHandle_mu4, primary_reco.key());
  // mcs_momentum_mu4 = primary_mcs_mu4->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_pi1(mcsHandle_pi1, primary_reco.key());
  // mcs_momentum_pi1 = primary_mcs_pi1->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_pi2(mcsHandle_pi2, primary_reco.key());
  // mcs_momentum_pi2 = primary_mcs_pi2->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_pi3(mcsHandle_pi3, primary_reco.key());
  // mcs_momentum_pi3 = primary_mcs_pi3->bestMomentum();
  // art::Ptr<recob::MCSFitResult> primary_mcs_pi4(mcsHandle_pi4, primary_reco.key());
  // mcs_momentum_pi4 = primary_mcs_pi4->bestMomentum();
  //  Calo momentum
  art::InputTag caliTag {std::string("pandoraNucali")};
  std::vector<art::Ptr<recob::Track>> tracks;
  tracks.push_back(primary_reco);
  art::FindManyP<anab::Calorimetry> caliTrackAss(tracks,evt,caliTag);
  std::vector<art::Ptr<anab::Calorimetry>> calis = caliTrackAss.at(0);
  for(int j=0; j!=3; j++)
  {    
    art::Ptr<anab::Calorimetry> cali = calis[j];
    int planeID = cali->PlaneID().Plane;
    if ((planeID < 0) || (planeID > 2)) continue;
    cali_caloPlane[planeID] = planeID;
    cali_isTrackFlipped[planeID] = IsTrackFlipped(cali);
    cali_kinEnergy[planeID] = cali->KineticEnergy();
    cali_range[planeID] = cali->Range();
    if (fDetailed_dEdX)
    {
      // Get the dEdX
      cali_dEdX[planeID].resize(cali->dEdx().size(),-999.);
      cali_dQdX[planeID].resize(cali->dQdx().size(),-999.);
      cali_resRange[planeID].resize(cali->ResidualRange().size(),-999.);
      for(int k=0; k<int(cali->dEdx().size()); k++)
      {
        cali_dEdX[planeID][k] = (cali->dEdx())[k];
        cali_dQdX[planeID][k] = (cali->dQdx())[k];
        cali_resRange[planeID][k] = (cali->ResidualRange())[k];
      } //  END loop for each hit
    } // END if detailed dEdX required
  } // END loop for each plane
  tDataTree->Fill();
} // END function analyze



// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(SingleParticleMomentumStudy)

#endif // END def SingleParticleMomentumStudy_module
