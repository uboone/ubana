#ifndef HSNMCTRUTHINFORMATION_MODULE
#define HSNMCTRUTHINFORMATION_MODULE

#include "HsnMcTruthInformation.h"

HsnMcTruthInformation::HsnMcTruthInformation(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMcTruthLabel(pset.get<std::string>("mcTruthLabel")),
    fMcTrackLabel(pset.get<std::string>("mcTrackLabel"))
{} // END constructor HsnMcTruthInformation

HsnMcTruthInformation::~HsnMcTruthInformation()
{} // END destructor HsnMcTruthInformation

void HsnMcTruthInformation::beginJob()
{
  art::ServiceHandle< art::TFileService > tfs;
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("pdgCode",&pdgCode);
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
  tDataTree->Branch("Vx",&Vx);
  tDataTree->Branch("Vy",&Vy);
  tDataTree->Branch("Vz",&Vz);
  tDataTree->Branch("T",&T);
  tDataTree->Branch("EndX",&EndX);
  tDataTree->Branch("EndY",&EndY);
  tDataTree->Branch("EndZ",&EndZ);
  tDataTree->Branch("EndT",&EndT);
  tDataTree->Branch("Px",&Px);
  tDataTree->Branch("Py",&Py);
  tDataTree->Branch("Pz",&Pz);
  tDataTree->Branch("E",&E);
  tDataTree->Branch("P",&P);
  tDataTree->Branch("Pt",&Pt);
  tDataTree->Branch("Length",&Length);
  tDataTree->Branch("Theta",&Theta);
  tDataTree->Branch("Phi",&Phi);
  tDataTree->Branch("Nu_Px",&Nu_Px);
  tDataTree->Branch("Nu_Py",&Nu_Py);
  tDataTree->Branch("Nu_Pz",&Nu_Pz);
  tDataTree->Branch("Nu_P",&Nu_P);
  tDataTree->Branch("Nu_E",&Nu_E);
  tDataTree->Branch("Nu_Theta",&Nu_Theta);
  tDataTree->Branch("Nu_Phi",&Nu_Phi);
  tDataTree->Branch("OpeningAngle",&OpeningAngle);
  tDataTree->Branch("InvariantMass",&InvariantMass);
  tDataTree->Branch("Contained",&Contained);

} // END function beginJob

void HsnMcTruthInformation::endJob()
{
} // END function endJob

void HsnMcTruthInformation::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  OpeningAngle = -999;
  InvariantMass = -999;
  pdgCode.clear();
  Vx.clear();
  Vy.clear();
  Vz.clear();
  T.clear();
  EndX.clear();
  EndY.clear();
  EndZ.clear();
  EndT.clear();
  Length.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  E.clear();
  P.clear();
  Pt.clear();
  Theta.clear();
  Phi.clear();
} // END function ClearData

float HsnMcTruthInformation::TrackLength(std::vector<float> start, std::vector<float> end)
{
  return sqrt(pow(start[0] - end[0],2.) + pow(start[1] - end[1],2.) + pow(start[2] - end[2],2.));

}

void HsnMcTruthInformation::GetTruthParticles(art::Event const & evt)
{
  // Prepare handle labels
  art::InputTag mcTruthTag {fMcTruthLabel};
  const auto& mcTruthHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTruthTag);

  art::InputTag mcTrackTag {fMcTrackLabel};
  const auto& mcTrackHandle = evt.getValidHandle< std::vector<sim::MCTrack> >(mcTrackTag);

  // Find mcTruth
  for(std::vector<int>::size_type i=0; i!=(*mcTruthHandle).size(); i++)
  {
    art::Ptr<simb::MCTruth> mcTruth(mcTruthHandle,i);
    int nParticles = mcTruth->NParticles();
    printf("|_Number of MCTruth: %i\n", nParticles);
    printf("|_Number of MCTracks: %i\n", (int) (*mcTrackHandle).size());
    printf("|\n");

    for (int j=0; j<nParticles; j++)
    {
      const simb::MCParticle & mcPart = mcTruth->GetParticle(j);
      art::Ptr<sim::MCTrack> mcTrack; 
      printf("|_Found MCPart (%i of %i) | PDG: %i\n", j+1, nParticles, mcPart.PdgCode());

      // Find mcTrack object associated with this mcPart. mcPart doesn't have simulation of interaction in argon,
      // so we can't recover the end points of tracks (and thus the lengths)
      // There is no mcPart/mcTrack association so at the moment we look for primary particles
      // that have the same PDG code. This immediately fails if an interaction contains two particles with
      // same pdg coming out of nucleus (albeit unlikely), but it should be kept to mind.
      printf("| |_Looping through candidates.\n");
      bool matchFound = false;
      for(std::vector<int>::size_type k=0; k!=(*mcTrackHandle).size(); k++)
      {
        art::Ptr<sim::MCTrack> potMcTrack(mcTrackHandle,k);
        bool mctIsPrimary = (potMcTrack->Process()=="primary");
        bool mctHasSamePdgCode = (potMcTrack->PdgCode()==mcPart.PdgCode());
        printf("| | |_Examining MCTrack candidate %i of %i | PDG: %i | Process: %s\n", (int) k+1, (int) (*mcTrackHandle).size(), potMcTrack->PdgCode(), potMcTrack->Process().c_str());
        if (mctIsPrimary && mctHasSamePdgCode)
        {
          printf("| | | |_Particle matches!\n");
          mcTrack = potMcTrack;
          matchFound = true;
        }
      }


      pdgCode.push_back(mcPart.PdgCode());
      Vx.push_back((float) mcPart.Vx());
      Vy.push_back((float) mcPart.Vy());
      Vz.push_back((float) mcPart.Vz());
      T.push_back((float) mcPart.T());
      if (matchFound)
      {
        EndX.push_back((float) mcTrack->End().X());
        EndY.push_back((float) mcTrack->End().Y());
        EndZ.push_back((float) mcTrack->End().Z());
        EndT.push_back((float) mcTrack->End().T());
        std::vector<float> start = {(float) mcPart.Vx(), (float) mcPart.Vy(), (float) mcPart.Vz()};
        std::vector<float> end = {(float) mcTrack->End().X(), (float) mcTrack->End().Y(), (float) mcTrack->End().Z()};
        Length.push_back(TrackLength(start, end));
      }
      else
      {
        EndX.push_back(-999999);
        EndY.push_back(-999999);
        EndZ.push_back(-999999);
        EndT.push_back(-999999);
        Length.push_back(-999999);
        Contained = false;
        printf("|_BAD EVENT! Not all matches found.\n");
      }

      Px.push_back((float) mcPart.Px());
      Py.push_back((float) mcPart.Py());
      Pz.push_back((float) mcPart.Pz());
      E.push_back((float) mcPart.E());
      P.push_back((float) mcPart.P());
      Pt.push_back((float) mcPart.Pt());

      // Calculate other quantities

      Theta.push_back((float) (acos(mcPart.Pz()/mcPart.P())));
      Phi.push_back((float) (atan2(mcPart.Py(),mcPart.Px())));
    }


    // Calculate combined quantities: opening angle and invariant mass
    if (nParticles==2)
    {
      // Calculate opening angle
      float dotProduct = Px[0]*Px[1] + Py[0]*Py[1] + Pz[0]*Pz[1];
      OpeningAngle = acos(dotProduct / float(P[0]*P[1]));
      // Calculate invariant mass
      float eTerm = pow((E[0] + E[1]),2.);
      float pTerm = pow(P[0],2.) + pow(P[1],2.) + 2.*dotProduct;
      InvariantMass = sqrt(eTerm - pTerm);

      // Calculate neutrino quantities
      Nu_Px = Px[0] + Px[1];
      Nu_Py = Py[0] + Py[1];
      Nu_Pz = Pz[0] + Pz[1];
      Nu_P = sqrt(pow(Nu_Px,2.) + pow(Nu_Py,2.) + pow(Nu_Pz,2.));
      Nu_E = sqrt(pow(InvariantMass,2.) + pow(Nu_P,2.));
      Nu_Theta = acos(Nu_Pz/Nu_P);
      Nu_Phi = atan2(Nu_Py,Nu_Px);
    }
    else
    {
      OpeningAngle = -999;
      InvariantMass = -999;
    }

  } // End of pfp loop

  return;
} // END function GetTruthParticles


void HsnMcTruthInformation::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  printf("\n-------------------------------------------------------\n");

  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  // Get vector of primaries and secondaries pfps
  GetTruthParticles(evt);

  // Fill tree and finish event loop
  tDataTree->Fill();
  printf("-------------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(HsnMcTruthInformation)

#endif // END def HsnMcTruthInformation_module
