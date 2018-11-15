#ifndef EXPLORETRUTH_MODULE
#define EXPLORETRUTH_MODULE

#include "ExploreTruth.h"

ExploreTruth::ExploreTruth(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor ExploreTruth

ExploreTruth::~ExploreTruth()
{} // END destructor ExploreTruth

void ExploreTruth::beginJob()
{
} // END function beginJob

void ExploreTruth::endJob()
{
} // END function endJob

void ExploreTruth::analyze(art::Event const & evt)
{
  // Determine event ID
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

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
    PrintInfo(primary,"Primary","");


    for(int j=0; j!=primary->NumberDaughters(); j++)
    {
      int daughterId = primary->Daughter(j);
      art::Ptr<simb::MCParticle> daughter;
      if (mcpMap.count(daughterId)) daughter = mcpMap[daughterId];
      else continue;

      std::string dProcess = daughter->Process();
      if (dProcess!=std::string("hIoni") && dProcess!=std::string("muIoni"))
      {
        PrintInfo(daughter,"Daughter","  ");
      }
    } // END loop through all daughters
  } // END loop through all primaries
  printf("\n\n");
} // END function analyze


void ExploreTruth::PrintInfo(art::Ptr<simb::MCParticle> mcp, std::string type, std::string spaces)
{
  printf("%s|_%s %i | PDG: %i, Process: %s\n", spaces.c_str(), type.c_str(), mcp->TrackId(), mcp->PdgCode(), mcp->Process().c_str());
  printf("%s  |_End process: %s\n", spaces.c_str(), mcp->EndProcess().c_str());
  printf("%s  |_Number of daughters: %i\n", spaces.c_str(), mcp->NumberDaughters());
  // printf("%s  |_Status code: %i\n",spaces.c_str(), mcp->StatusCode());
  printf("%s  |_Start coordinates: [ %.1f, %.1f, %.1f, %.1f ]\n", spaces.c_str(), mcp->Vx(), mcp->Vy(), mcp->Vz(), mcp->T());
  printf("%s  |_End coordinates: [ %.1f, %.1f, %.1f, %.1f ]\n", spaces.c_str(), mcp->EndX(), mcp->EndY(), mcp->EndZ(), mcp->EndT());
  printf("%s  |_Start/End energies: [ %.3f, %.3f ]\n", spaces.c_str(), mcp->E(), mcp->EndE());
}

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ExploreTruth)

#endif // END def ExploreTruth_module
