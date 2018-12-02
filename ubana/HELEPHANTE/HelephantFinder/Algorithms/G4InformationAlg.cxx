/******************************************************************************
 * @file G4InformationAlg.cxx
 * @brief Retrieve G4 information regarding the current event (e.g. truth track start-end points, if they decay/interact/absorb) and compare G4 Truth momentum with reconstructed momentum using a very basic reco-truth matching
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see G4InformationAlg.h
 * ****************************************************************************/
#include "G4InformationAlg.h"

namespace G4Information
{
  // Constructor/destructor
  G4InformationAlg::G4InformationAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
  }
  G4InformationAlg::~G4InformationAlg()
  {}
  void G4InformationAlg::reconfigure(fhicl::ParameterSet const & pset)
  {}

  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void G4InformationAlg::AddG4Information(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf,
            AuxVertex::DecayVertex const & dv)
  {
    printf("DV1 stuff: [ %.2f, %.2f ]\n", dv.fProngMomMag_ByRange_h1[0], dv.fProngMomMag_ByRange_h1[1]);
    // Get all the handles to McTruth objects.
    std::string partLabel("largeant");
    art::InputTag partTag {partLabel};
    const auto& partHandle = evt.getValidHandle< std::vector<simb::MCParticle> >(partTag);

    // TRACKID - MCPARTICLE MAP
    // There's no way to go from trackID number to the actual McParticle it refers to,
    // so you'll need to build a map of trackID-pointers first.
    std::map< int,art::Ptr<simb::MCParticle> > mcpMap;
    std::vector< art::Ptr<simb::MCParticle> > primaries;
    for(int i=0; i!=int((*partHandle).size()); i++)
    {
      art::Ptr<simb::MCParticle> map_mcp(partHandle,i);
      int trackId = map_mcp->TrackId();
      mcpMap[trackId] = map_mcp;
      if(map_mcp->Process().c_str() == std::string("primary"))
      {
        primaries.push_back(map_mcp);
      }
    }

    // FIND INTERACTION / DECAY / ABSORPTION
    // Loop through each primary
    for(int i=0; i!=int(primaries.size()); i++)
    {
      art::Ptr<simb::MCParticle> primary = primaries[i];
      std::string originalDeclaredEndProcess = primary->EndProcess();

      // ELASTIC SCATTERING CASE
      // Code I'm not particularly proud of. I can write clearer than that.
      // Check if any daughters has same pdg as primary and replace primary with daughter
      // Keep doing until there are no more elastic scatterings
      bool stillElastic = true;
      int numberElastic = 0;
      art::Ptr<simb::MCParticle> tempPrimary;
      while (stillElastic)
      {
        stillElastic = false;
        int numberDaughters = primary->NumberDaughters();
        for(int j=0; j!=numberDaughters; j++)
        {
          int daughterId = primary->Daughter(j);
          art::Ptr<simb::MCParticle> daughter;
          if (mcpMap.count(daughterId)) daughter = mcpMap[daughterId];
          else continue;

          if (daughter->PdgCode() == primary->PdgCode())
          {
            printf("Primary: %i, Daughter: %i\n", primary->TrackId(), daughter->TrackId());
            printf("\n");
            numberElastic += 1;
            tempPrimary = daughter;
            stillElastic = true;
          } // END comparing iterative pdgCode
        } // END iterative for loop through daughters of daughters
      if (stillElastic==true) primary = tempPrimary;
      } // END iterative while loop for elastic scattering


      std::string declaredEndProcess = primary->EndProcess();
      std::string otherEndProcess = std::string("Other");

      // DECAY / INTERACTION STUDY
      // Loop through the daughters and look for case by case 
      for(int j=0; j!=primary->NumberDaughters(); j++)
      {
        int daughterId = primary->Daughter(j);
        art::Ptr<simb::MCParticle> daughter;
        if (mcpMap.count(daughterId)) daughter = mcpMap[daughterId];
        else continue;

        if (primary->PdgCode() == 211 || primary->PdgCode() == -211)
        {
          if (daughter->Process() == std::string("Decay")) otherEndProcess = std::string("Decay");
          if (daughter->Process() == std::string("pi+Inelastic")) otherEndProcess = std::string("pi+Inelastic");
          if (daughter->Process() == std::string("pi-Inelastic")) otherEndProcess = std::string("pi-Inelastic");
          else if (daughter->PdgCode() == 111) otherEndProcess = std::string("Pi0");
          // Now you have all the information you need. Fill the candidate tree
          float piMass = 0.139570;
          ctf.g4_pi_pdg = primary->PdgCode();
          ctf.g4_pi_mom = primary->P();
          ctf.g4_pi_endMom = sqrt(pow(primary->EndE(),2.) - pow(piMass,2.));
          ctf.g4_pi_length = McpLength(primary);
          ctf.g4_pi_nScatterings = numberElastic;
          ctf.g4_pi_declaredEndProcess = originalDeclaredEndProcess;
          ctf.g4_pi_declaredScatteredEndProcess = declaredEndProcess;
          ctf.g4_pi_endProcess = otherEndProcess;
        } // END case study for pion  antipion
        if (primary->PdgCode() == 13 || primary->PdgCode() == -13)
        {
          if (daughter->Process() == std::string("muMinusCaptureAtRest")) otherEndProcess = std::string("muMinusCaptureAtRest");
          if (daughter->Process() == std::string("muPlusCaptureAtRest")) otherEndProcess = std::string("muPlusCaptureAtRest");
          if (daughter->Process() == std::string("Decay")) otherEndProcess = std::string("Decay");
          // Now you have all the information you need. Fill the candidate tree
          float muMass = 0.105658;
          ctf.g4_mu_pdg = primary->PdgCode();
          ctf.g4_mu_mom = primary->P();
          ctf.g4_mu_endMom = sqrt(pow(primary->EndE(),2.) - pow(muMass,2.));
          ctf.g4_mu_length = McpLength(primary);
          ctf.g4_mu_nScatterings = numberElastic;
          ctf.g4_mu_declaredEndProcess = originalDeclaredEndProcess;
          ctf.g4_mu_declaredScatteredEndProcess = declaredEndProcess;
          ctf.g4_mu_endProcess = otherEndProcess;
        } // END case study for muon / antimuon
      } // END loop through all daughters of current primary
    } // END loop through all primaries in event

    // Perform a very basic truth-reco matching, assign the g4 track with one of the two decay vertex prongs based on proximity of end points. Then, store the momenta from truth, range and mcs.
    printf("DV2 stuff: [ %.2f, %.2f ]\n", dv.fProngMomMag_ByRange_h1[0], dv.fProngMomMag_ByRange_h1[1]);
    BadTruthRecoMatching(primaries, dv, ctf);
  } // END function AddG4Information


  float G4InformationAlg::McpLength(art::Ptr<simb::MCParticle> mcp)
  {
    float x1 = mcp->Vx();
    float y1 = mcp->Vy();
    float z1 = mcp->Vz();
    float x2 = mcp->EndX();
    float y2 = mcp->EndY();
    float z2 = mcp->EndZ();
    float length = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
    return length;
  }

  void G4InformationAlg::BadTruthRecoMatching(
            std::vector< art::Ptr<simb::MCParticle> > const & mcps, 
            AuxVertex::DecayVertex const & dv,
            AuxEvent::CandidateTreeFiller & ctf)
  {
    art::Ptr<simb::MCParticle> g4muon, g4pion;
    for(int i=0; i!=int(mcps.size()); i++)
    {
      art::Ptr<simb::MCParticle> mcp = mcps[i];
      if ((mcp->PdgCode()==13) || (mcp->PdgCode()==-13)) g4muon = mcp;
      if ((mcp->PdgCode()==211) || (mcp->PdgCode()==-211)) g4pion = mcp;
    }

    // Hypo1 (0:mu, 1:pi)
    float h1Dist1 = End2EndDistance(0,g4muon,dv);
    float h1Dist2 = End2EndDistance(1,g4pion,dv);
    float h1Score = sqrt(pow(h1Dist1,2.)+pow(h1Dist2,2.));
    // Hypo2 (1:mu, 0:pi)
    float h2Dist1 = End2EndDistance(1,g4muon,dv);
    float h2Dist2 = End2EndDistance(0,g4pion,dv);
    float h2Score = sqrt(pow(h2Dist1,2.)+pow(h2Dist2,2.));

    // Diagnostics
    printf("HYPOTHESIS 1:\n");
    printf("|_Distances: [ %.1f, %.1f ]\n", h1Dist1, h1Dist2);
    printf("|_Score: %.1f\n", h1Score);
    printf("HYPOTHESIS 2:\n");
    printf("|_Distances: [ %.1f, %.1f ]\n", h2Dist1, h2Dist2);
    printf("|_Score: %.1f\n", h2Score);
    printf("\n");

    // Store scores and distances
    std::vector<float> v1{h1Dist1,h1Dist2, h1Score};
    std::vector<float> v2{h2Dist1,h2Dist2, h2Score};
    std::vector<std::vector<float>> distScores{v1,v2};
    ctf.g4_recoTruth_distancesScore = distScores;
    if (h1Score<h2Score)
    {
      ctf.g4_reco_mu_momRange = dv.fProngMomMag_ByRange_h1[0];
      ctf.g4_reco_pi_momRange = dv.fProngMomMag_ByRange_h1[1];
      ctf.g4_reco_mu_momMCS = dv.fProngMomMag_ByMcs_best_h1[0];
      ctf.g4_reco_pi_momMCS = dv.fProngMomMag_ByMcs_best_h1[1];
    }
    else
    {
      ctf.g4_reco_mu_momRange = dv.fProngMomMag_ByRange_h1[1];
      ctf.g4_reco_pi_momRange = dv.fProngMomMag_ByRange_h1[0];
      ctf.g4_reco_mu_momMCS = dv.fProngMomMag_ByMcs_best_h1[1];
      ctf.g4_reco_pi_momMCS = dv.fProngMomMag_ByMcs_best_h1[0];
    }
  }



  float G4InformationAlg::End2EndDistance(
            int prongID,
            art::Ptr<simb::MCParticle> const & mcp,
            AuxVertex::DecayVertex const & dv)
  {
    float x1 = dv.fProngEndX[prongID];
    float y1 = dv.fProngEndY[prongID];
    float z1 = dv.fProngEndZ[prongID];
    float x2 = mcp->EndX();
    float y2 = mcp->EndY();
    float z2 = mcp->EndZ();
    float distance = sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)+pow(z1-z2,2.));
    return distance; 
  }
} // END namespace G4Information
