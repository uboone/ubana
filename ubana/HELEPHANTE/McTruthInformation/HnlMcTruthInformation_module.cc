#ifndef HSNMCTRUTHINFORMATION_MODULE
#define HSNMCTRUTHINFORMATION_MODULE
#include "HnlMcTruthInformation.h"

HnlMcTruthInformation::HnlMcTruthInformation(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor HnlMcTruthInformation

HnlMcTruthInformation::~HnlMcTruthInformation()
{} // END destructor HnlMcTruthInformation

void HnlMcTruthInformation::beginJob()
{
  art::ServiceHandle< art::TFileService > tfs;
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
  // NU + GENERIC
  tDataTree->Branch("X_nu",&X_nu);
  tDataTree->Branch("Y_nu",&Y_nu);
  tDataTree->Branch("Z_nu",&Z_nu);
  tDataTree->Branch("T_nu",&T_nu);
  tDataTree->Branch("pX_nu",&pX_nu);
  tDataTree->Branch("pY_nu",&pY_nu);
  tDataTree->Branch("pZ_nu",&pZ_nu);
  tDataTree->Branch("p_nu",&p_nu);
  tDataTree->Branch("e_nu",&e_nu);
  tDataTree->Branch("theta_nu",&theta_nu);
  tDataTree->Branch("phi_nu",&phi_nu);
  tDataTree->Branch("openingAngle",&openingAngle);
  tDataTree->Branch("invariantMass",&invariantMass);
  tDataTree->Branch("isContained",&isContained);
  // MUON MCPARTICLE
  tDataTree->Branch("process_mu",&process_mu);
  tDataTree->Branch("endProcess_mu",&endProcess_mu);
  tDataTree->Branch("pdgCode_mu",&pdgCode_mu);
  tDataTree->Branch("start_X_mu",&start_X_mu);
  tDataTree->Branch("start_Y_mu",&start_Y_mu);
  tDataTree->Branch("start_Z_mu",&start_Z_mu);
  tDataTree->Branch("start_T_mu",&start_T_mu);
  tDataTree->Branch("end_X_mu",&end_X_mu);
  tDataTree->Branch("end_Y_mu",&end_Y_mu);
  tDataTree->Branch("end_Z_mu",&end_Z_mu);
  tDataTree->Branch("end_T_mu",&end_T_mu);
  tDataTree->Branch("start_pX_mu",&start_pX_mu);
  tDataTree->Branch("start_pY_mu",&start_pY_mu);
  tDataTree->Branch("start_pZ_mu",&start_pZ_mu);
  tDataTree->Branch("start_p_mu",&start_p_mu);
  tDataTree->Branch("start_e_mu",&start_e_mu);
  tDataTree->Branch("end_pX_mu",&end_pX_mu);
  tDataTree->Branch("end_pY_mu",&end_pY_mu);
  tDataTree->Branch("end_pZ_mu",&end_pZ_mu);
  tDataTree->Branch("end_p_mu",&end_p_mu);
  tDataTree->Branch("end_e_mu",&end_e_mu);
  tDataTree->Branch("length_mu",&length_mu);
  tDataTree->Branch("theta_mu",&theta_mu);
  tDataTree->Branch("phi_mu",&phi_mu);
  // MUON MCTRACK
  tDataTree->Branch("track_start_X_mu",&track_start_X_mu);
  tDataTree->Branch("track_start_Y_mu",&track_start_Y_mu);
  tDataTree->Branch("track_start_Z_mu",&track_start_Z_mu);
  tDataTree->Branch("track_start_T_mu",&track_start_T_mu);
  tDataTree->Branch("track_end_X_mu",&track_end_X_mu);
  tDataTree->Branch("track_end_Y_mu",&track_end_Y_mu);
  tDataTree->Branch("track_end_Z_mu",&track_end_Z_mu);
  tDataTree->Branch("track_end_T_mu",&track_end_T_mu);
  tDataTree->Branch("track_start_pX_mu",&track_start_pX_mu);
  tDataTree->Branch("track_start_pY_mu",&track_start_pY_mu);
  tDataTree->Branch("track_start_pZ_mu",&track_start_pZ_mu);
  tDataTree->Branch("track_start_p_mu",&track_start_p_mu);
  tDataTree->Branch("track_start_e_mu",&track_start_e_mu);
  tDataTree->Branch("track_end_pX_mu",&track_end_pX_mu);
  tDataTree->Branch("track_end_pY_mu",&track_end_pY_mu);
  tDataTree->Branch("track_end_pZ_mu",&track_end_pZ_mu);
  tDataTree->Branch("track_end_p_mu",&track_end_p_mu);
  tDataTree->Branch("track_end_e_mu",&track_end_e_mu);
  tDataTree->Branch("track_length_mu",&track_length_mu);
  tDataTree->Branch("track_theta_mu",&track_theta_mu);
  tDataTree->Branch("track_phi_mu",&track_phi_mu);
  // MUON OTHER
  tDataTree->Branch("hasMCParticle_mu",&hasMCParticle_mu);
  tDataTree->Branch("hasMCTrack_mu",&hasMCTrack_mu);
  // PION MCPARTICLE
  tDataTree->Branch("process_pi",&process_pi);
  tDataTree->Branch("endProcess_pi",&endProcess_pi);
  tDataTree->Branch("pdgCode_pi",&pdgCode_pi);
  tDataTree->Branch("start_X_pi",&start_X_pi);
  tDataTree->Branch("start_Y_pi",&start_Y_pi);
  tDataTree->Branch("start_Z_pi",&start_Z_pi);
  tDataTree->Branch("start_T_pi",&start_T_pi);
  tDataTree->Branch("end_X_pi",&end_X_pi);
  tDataTree->Branch("end_Y_pi",&end_Y_pi);
  tDataTree->Branch("end_Z_pi",&end_Z_pi);
  tDataTree->Branch("end_T_pi",&end_T_pi);
  tDataTree->Branch("start_pX_pi",&start_pX_pi);
  tDataTree->Branch("start_pY_pi",&start_pY_pi);
  tDataTree->Branch("start_pZ_pi",&start_pZ_pi);
  tDataTree->Branch("start_p_pi",&start_p_pi);
  tDataTree->Branch("start_e_pi",&start_e_pi);
  tDataTree->Branch("end_pX_pi",&end_pX_pi);
  tDataTree->Branch("end_pY_pi",&end_pY_pi);
  tDataTree->Branch("end_pZ_pi",&end_pZ_pi);
  tDataTree->Branch("end_p_pi",&end_p_pi);
  tDataTree->Branch("end_e_pi",&end_e_pi);
  tDataTree->Branch("length_pi",&length_pi);
  tDataTree->Branch("theta_pi",&theta_pi);
  tDataTree->Branch("phi_pi",&phi_pi);
  // PION MCTRACK
  tDataTree->Branch("track_start_X_pi",&track_start_X_pi);
  tDataTree->Branch("track_start_Y_pi",&track_start_Y_pi);
  tDataTree->Branch("track_start_Z_pi",&track_start_Z_pi);
  tDataTree->Branch("track_start_T_pi",&track_start_T_pi);
  tDataTree->Branch("track_end_X_pi",&track_end_X_pi);
  tDataTree->Branch("track_end_Y_pi",&track_end_Y_pi);
  tDataTree->Branch("track_end_Z_pi",&track_end_Z_pi);
  tDataTree->Branch("track_end_T_pi",&track_end_T_pi);
  tDataTree->Branch("track_start_pX_pi",&track_start_pX_pi);
  tDataTree->Branch("track_start_pY_pi",&track_start_pY_pi);
  tDataTree->Branch("track_start_pZ_pi",&track_start_pZ_pi);
  tDataTree->Branch("track_start_p_pi",&track_start_p_pi);
  tDataTree->Branch("track_start_e_pi",&track_start_e_pi);
  tDataTree->Branch("track_end_pX_pi",&track_end_pX_pi);
  tDataTree->Branch("track_end_pY_pi",&track_end_pY_pi);
  tDataTree->Branch("track_end_pZ_pi",&track_end_pZ_pi);
  tDataTree->Branch("track_end_p_pi",&track_end_p_pi);
  tDataTree->Branch("track_end_e_pi",&track_end_e_pi);
  tDataTree->Branch("track_length_pi",&track_length_pi);
  tDataTree->Branch("track_theta_pi",&track_theta_pi);
  tDataTree->Branch("track_phi_pi",&track_phi_pi);
  // PION OTHER
  tDataTree->Branch("hasMCParticle_pi",&hasMCParticle_pi);
  tDataTree->Branch("hasMCTrack_pi",&hasMCTrack_pi);
} // END function beginJob

void HnlMcTruthInformation::endJob()
{
} // END function endJob

void HnlMcTruthInformation::ClearData()
{
  // MUON MCPARTICLE
  pdgCode_mu = -999;
  start_X_mu = -999;
  start_Y_mu = -999;
  start_Z_mu = -999;
  start_T_mu = -999;
  end_X_mu = -999;
  end_Y_mu = -999;
  end_Z_mu = -999;
  end_T_mu= -999;
  start_pX_mu = -999;
  start_pY_mu = -999;
  start_pZ_mu = -999;
  start_p_mu = -999;
  start_e_mu = -999;
  end_pX_mu = -999;
  end_pY_mu = -999;
  end_pZ_mu = -999;
  end_p_mu = -999;
  end_e_mu = -999;
  length_mu = -999;
  theta_mu = -999;
  phi_mu = -999;
  // MUON MCTRACK
  track_start_X_mu = -999;
  track_start_Y_mu = -999;
  track_start_Z_mu = -999;
  track_start_T_mu = -999;
  track_end_X_mu = -999;
  track_end_Y_mu = -999;
  track_end_Z_mu = -999;
  track_end_T_mu= -999;
  track_start_pX_mu = -999;
  track_start_pY_mu = -999;
  track_start_pZ_mu = -999;
  track_start_p_mu = -999;
  track_start_e_mu = -999;
  track_end_pX_mu = -999;
  track_end_pY_mu = -999;
  track_end_pZ_mu = -999;
  track_end_p_mu = -999;
  track_end_e_mu = -999;
  track_length_mu = -999;
  track_theta_mu = -999;
  track_phi_mu = -999;
  // PION MCPARTICLE
  pdgCode_pi = -999;
  start_X_pi = -999;
  start_Y_pi = -999;
  start_Z_pi = -999;
  start_T_pi = -999;
  end_X_pi = -999;
  end_Y_pi = -999;
  end_Z_pi = -999;
  end_T_pi = -999;
  start_pX_pi = -999;
  start_pY_pi = -999;
  start_pZ_pi = -999;
  start_p_pi = -999;
  start_e_pi = -999;
  end_pX_pi = -999;
  end_pY_pi = -999;
  end_pZ_pi = -999;
  end_p_pi = -999;
  end_e_pi = -999;
  length_pi = -999;
  theta_pi = -999;
  phi_pi = -999;
  // PION MCTRACK
  track_start_X_pi = -999;
  track_start_Y_pi = -999;
  track_start_Z_pi = -999;
  track_start_T_pi = -999;
  track_end_X_pi = -999;
  track_end_Y_pi = -999;
  track_end_Z_pi = -999;
  track_end_T_pi = -999;
  track_start_pX_pi = -999;
  track_start_pY_pi = -999;
  track_start_pZ_pi = -999;
  track_start_p_pi = -999;
  track_start_e_pi = -999;
  track_end_pX_pi = -999;
  track_end_pY_pi = -999;
  track_end_pZ_pi = -999;
  track_end_p_pi = -999;
  track_end_e_pi = -999;
  track_length_pi = -999;
  track_theta_pi = -999;
  track_phi_pi = -999;
  // NEUTRINO
  X_nu = -999;
  Y_nu = -999;
  Z_nu = -999;
  T_nu = -999;
  pX_nu = -999;
  pY_nu = -999;
  pZ_nu = -999;
  p_nu = -999;
  e_nu = -999;
  theta_nu = -999;
  phi_nu = -999;
  openingAngle = -999;
  invariantMass = -999;
  // OTHERS
  isContained = false;
  hasMCParticle_pi = false;
  hasMCParticle_mu = false;
  hasMCTrack_pi = false;
  hasMCTrack_mu = false;
} // END function ClearData



// FUNCTION
// Build map of trackId - McParticle pointers
std::map< int,art::Ptr<simb::MCParticle> > HnlMcTruthInformation::BuildMcpMap(art::ValidHandle< std::vector<simb::MCParticle> > const & mcpHandle)
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
std::map< int,art::Ptr<sim::MCTrack> > HnlMcTruthInformation::BuildMctMap(art::ValidHandle< std::vector<sim::MCTrack> > const & mctHandle)
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
// Determine MCParticle length
float HnlMcTruthInformation::MCParticleLength(art::Ptr<simb::MCParticle> const & mcp)
{
  float x1, x2, y1, y2, z1, z2;
  x1 = mcp->Vx();
  y1 = mcp->Vy();
  z1 = mcp->Vz();
  x2 = mcp->EndX();
  y2 = mcp->EndY();
  z2 = mcp->EndZ();
  float length = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));  
  return length;
}



// FUNCTION
// Determine Track length
float HnlMcTruthInformation::MCTrackLength(art::Ptr<sim::MCTrack> const & mct)
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
}



void HnlMcTruthInformation::analyze(art::Event const & evt)
{
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  // Prepare handle labels
  // art::InputTag mcTruthTag {std::string("generator")};
  // art::ValidHandle<std::vector<simb::MCTruth>> mcTruthHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTruthTag);
  art::InputTag mcParticleTag {std::string("largeant")};
  art::ValidHandle< std::vector<simb::MCParticle>> mcParticleHandle = evt.getValidHandle< std::vector<simb::MCParticle> >(mcParticleTag);
  art::InputTag mcTrackTag {std::string("mcreco")};
  art::ValidHandle< std::vector<sim::MCTrack>> mcTrackHandle = evt.getValidHandle< std::vector<sim::MCTrack> >(mcTrackTag);
  // Prepare maps of TrackID -> McParticle/McTrack
  std::map< int,art::Ptr<simb::MCParticle> > mcpMap = BuildMcpMap(mcParticleHandle);
  std::map< int,art::Ptr<sim::MCTrack> > mctMap = BuildMctMap(mcTrackHandle);

  // Find muon and pion
  art::Ptr<simb::MCParticle> mcp_mu, mcp_pi;
  art::Ptr<sim::MCTrack> mct_mu, mct_pi;
  int trackID_mu, trackID_pi;
  int nPrimary_mu, nPrimary_pi;
  nPrimary_mu = 0;
  nPrimary_pi = 0;
  for(std::vector<int>::size_type i=0; i!=(*mcParticleHandle).size(); i++)
  {
    art::Ptr<simb::MCParticle> mcp(mcParticleHandle,i);
    if(mcp->Process() == std::string("primary"))
    {
      printf("-------- LOOOOK HERE!!!!\n");
      printf("PDG CODE %i\n", mcp->PdgCode());
      if((mcp->PdgCode() == 13) || (mcp->PdgCode() == -13))
      {
        trackID_mu = mcp->TrackId();
        nPrimary_mu += 1;
        mcp_mu = mcpMap[trackID_mu];
        hasMCParticle_mu = true;
      } // END if particle is muon
      if((mcp->PdgCode() == 211) || (mcp->PdgCode() == -211))
      {
        trackID_pi = mcp->TrackId();
        nPrimary_pi += 1;
        mcp_pi = mcpMap[trackID_pi];
        hasMCParticle_pi = true;
      } // END if particle is pion
    } // END if particle is primary
  } // END loop through all MCParticles

  // Check that we have only one type of primary for each (muon or pion). If we have more, something went wrong
  if ((nPrimary_mu > 1) || (nPrimary_pi > 1))
  {
    printf("ERROR!!!\nThere is more than one pion or more than one muon among the primary MCParticles.\nThis is incompatible with a HNL event. This module is written only to analyze mu-pi decays of HNL.\n");
      std::exit(1);
  }

  // Now that we're sure we got the right products, let's find also the MCTracks
  if (mctMap.count(trackID_mu))
  {
    mct_mu = mctMap[trackID_mu];
    hasMCTrack_mu = true;
  }
  if (mctMap.count(trackID_pi))
  {
    mct_pi = mctMap[trackID_pi];
    hasMCTrack_pi = true;
  }

  // We have all MCParticle and MCTracks. Let's derive the quantities we need
  // DO IT FOR MUON
  if(hasMCParticle_mu)
  {
    pdgCode_mu = mcp_mu->PdgCode();
    // START
    process_mu = mcp_mu->Process();
    start_X_mu = mcp_mu->Vx();
    start_Y_mu = mcp_mu->Vy();
    start_Z_mu = mcp_mu->Vz();
    start_T_mu = mcp_mu->T();
    start_pX_mu = mcp_mu->Px();
    start_pY_mu = mcp_mu->Py();
    start_pZ_mu = mcp_mu->Pz();
    start_p_mu = mcp_mu->P();
    start_e_mu = mcp_mu->E();
    // END
    endProcess_mu = mcp_mu->EndProcess();
    end_X_mu = mcp_mu->EndX();
    end_Y_mu = mcp_mu->EndY();
    end_Z_mu = mcp_mu->EndZ();
    end_T_mu = mcp_mu->EndT();
    end_pX_mu = mcp_mu->EndPx();
    end_pY_mu = mcp_mu->EndPy();
    end_pZ_mu = mcp_mu->EndPz();
    end_p_mu = sqrt(pow(end_pX_mu,2.)+pow(end_pY_mu,2.)+pow(end_pZ_mu,2.));
    end_e_mu = mcp_mu->EndE();
    length_mu = MCParticleLength(mcp_mu);
    theta_mu = acos(start_pZ_mu/start_p_mu);
    phi_mu = atan2(start_pY_mu,start_pX_mu);
    if(hasMCTrack_mu)
    {
      // START
      track_start_X_mu = mct_mu->Start().X();
      track_start_Y_mu = mct_mu->Start().Y();
      track_start_Z_mu = mct_mu->Start().Z();
      track_start_T_mu = mct_mu->Start().T();
      track_start_pX_mu = mct_mu->Start().Px();
      track_start_pY_mu = mct_mu->Start().Py();
      track_start_pZ_mu = mct_mu->Start().Pz();
      track_start_p_mu = sqrt(pow(track_start_pX_mu,2.)+pow(track_start_pY_mu,2.)+pow(track_start_pZ_mu,2.));
      track_start_e_mu = mct_mu->Start().E();
      // END
      track_end_X_mu = mct_mu->End().X();
      track_end_Y_mu = mct_mu->End().Y();
      track_end_Z_mu = mct_mu->End().Z();
      track_end_T_mu = mct_mu->End().T();
      track_end_pX_mu = mct_mu->End().Px();
      track_end_pY_mu = mct_mu->End().Py();
      track_end_pZ_mu = mct_mu->End().Pz();
      track_end_p_mu = sqrt(pow(track_end_pX_mu,2.)+pow(track_end_pY_mu,2.)+pow(track_end_pZ_mu,2.));
      track_end_e_mu = mct_mu->End().E();
      // OTHERS
      track_length_mu = MCTrackLength(mct_mu);
      track_theta_mu = acos(track_start_pZ_mu/track_start_p_mu);
      track_phi_mu = atan2(track_start_pY_mu,track_start_pX_mu);
    } // END if muon has McTrack
  } // END if muon has MCParticle

  // DO IT FOR PION
  if(hasMCParticle_pi)
  {
    pdgCode_pi = mcp_pi->PdgCode();
    // START
    process_pi = mcp_pi->Process();
    start_X_pi = mcp_pi->Vx();
    start_Y_pi = mcp_pi->Vy();
    start_Z_pi = mcp_pi->Vz();
    start_T_pi = mcp_pi->T();
    start_pX_pi = mcp_pi->Px();
    start_pY_pi = mcp_pi->Py();
    start_pZ_pi = mcp_pi->Pz();
    start_p_pi = mcp_pi->P();
    start_e_pi = mcp_pi->E();
    // END
    endProcess_pi = mcp_pi->EndProcess();
    end_X_pi = mcp_pi->EndX();
    end_Y_pi = mcp_pi->EndY();
    end_Z_pi = mcp_pi->EndZ();
    end_T_pi = mcp_pi->EndT();
    end_pX_pi = mcp_pi->EndPx();
    end_pY_pi = mcp_pi->EndPy();
    end_pZ_pi = mcp_pi->EndPz();
    end_p_pi = sqrt(pow(end_pX_pi,2.)+pow(end_pY_pi,2.)+pow(end_pZ_pi,2.));
    end_e_pi = mcp_pi->EndE();
    length_pi = MCParticleLength(mcp_pi);
    theta_pi = acos(start_pZ_pi/start_p_pi);
    phi_pi = atan2(start_pY_pi,start_pX_pi);
    if(hasMCTrack_pi)
    {
      // START
      track_start_X_pi = mct_pi->Start().X();
      track_start_Y_pi = mct_pi->Start().Y();
      track_start_Z_pi = mct_pi->Start().Z();
      track_start_T_pi = mct_pi->Start().T();
      track_start_pX_pi = mct_pi->Start().Px();
      track_start_pY_pi = mct_pi->Start().Py();
      track_start_pZ_pi = mct_pi->Start().Pz();
      track_start_p_pi = sqrt(pow(track_start_pX_pi,2.)+pow(track_start_pY_pi,2.)+pow(track_start_pZ_pi,2.));
      track_start_e_pi = mct_pi->Start().E();
      // END
      track_end_X_pi = mct_pi->End().X();
      track_end_Y_pi = mct_pi->End().Y();
      track_end_Z_pi = mct_pi->End().Z();
      track_end_T_pi = mct_pi->End().T();
      track_end_pX_pi = mct_pi->End().Px();
      track_end_pY_pi = mct_pi->End().Py();
      track_end_pZ_pi = mct_pi->End().Pz();
      track_end_p_pi = sqrt(pow(track_end_pX_pi,2.)+pow(track_end_pY_pi,2.)+pow(track_end_pZ_pi,2.));
      track_end_e_pi = mct_pi->End().E();
      // OTHERS
      track_length_pi = MCTrackLength(mct_pi);
      track_theta_pi = acos(track_start_pZ_pi/track_start_p_pi);
      track_phi_pi = atan2(track_start_pY_pi,track_start_pX_pi);
    } // END if pion has McTrack
  } // END if pion has MCParticle

  // DO IT FOR NEUTRINO
  if (hasMCParticle_pi && hasMCParticle_mu)
  {
    X_nu = start_X_mu;
    Y_nu = start_Y_mu;
    Z_nu = start_Z_mu;
    T_nu = start_T_mu;
    pX_nu = start_pX_mu + start_pX_pi;
    pY_nu = start_pY_mu + start_pY_pi;
    pZ_nu = start_pZ_mu + start_pZ_pi;
    // Calculate opening angle
    float dotProduct = start_pX_mu*start_pX_pi + start_pY_mu*start_pY_pi + start_pZ_mu*start_pZ_pi;
    openingAngle = acos(dotProduct / float(start_p_mu*start_p_pi));
    // Calculate invariant mass
    float eTerm = pow((start_e_mu + start_e_pi),2.);
    float pTerm = pow(start_p_mu,2.) + pow(start_p_pi,2.) + 2.*dotProduct;
    invariantMass = sqrt(eTerm - pTerm);
    p_nu = sqrt(pow(pX_nu,2.) + pow(pY_nu,2.) + pow(pZ_nu,2.));
    e_nu = sqrt(pow(invariantMass,2.) + pow(p_nu,2.));
    theta_nu = acos(pZ_nu/p_nu);
    phi_nu = atan2(pY_nu,pX_nu);
  } // END if there's MCParticles for both 

  // Fill tree and finish event loop
  tDataTree->Fill();
}


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(HnlMcTruthInformation)

#endif // END def HnlMcTruthInformation_module
