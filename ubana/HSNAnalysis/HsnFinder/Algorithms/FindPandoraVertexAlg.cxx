#include "FindPandoraVertexAlg.h"

namespace FindPandoraVertex
{
  // Constructor/destructor
  FindPandoraVertexAlg::FindPandoraVertexAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
  }
  FindPandoraVertexAlg::~FindPandoraVertexAlg()
  {}

  void FindPandoraVertexAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fMinTpcBound = pset.get<std::vector<double>>("MinTpcBound");
    fMaxTpcBound = pset.get<std::vector<double>>("MaxTpcBound");
    fPfpLabel = pset.get<std::string>("PfpLabel");
    fMcsLabel = pset.get<std::string>("McsLabel");
    fVerbose = pset.get<bool>("VerboseMode");
  }


  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void FindPandoraVertexAlg::GetPotentialNeutrinoVertices(
            art::Event const & evt,
            AuxEvent::EventTreeFiller & etf,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices)
  {
    if (fVerbose) printf("\n--- GetPotentialNeutrinoVertices message ---\n");

    // Clear stuff that will be modified by function
    etf.nNeutrinos = 0;
    etf.nTwoProngedNeutrinos = 0;
    etf.nContainedTwoProngedNeutrinos = 0;
    etf.neutrinoPdgCode.clear();
    etf.neutrinoNumDaughters.clear();
    etf.neutrinoNumTracks.clear();
    etf.neutrinoNumShowers.clear();
    etf.status_nuWithMissingAssociatedVertex = 0;
    etf.status_nuWithMissingAssociatedTrack = 0;
    etf.status_nuProngWithMissingAssociatedHits = 0;

    //Prepare the pfp handle
    art::InputTag pfpTag {fPfpLabel};
    art::InputTag mcsTag {fMcsLabel};
    const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
    // const auto& trackHandle = evt.getValidHandle< std::vector<recob::Track> >(pfpTag);
    const auto& mcsHandle = evt.getValidHandle< std::vector<recob::MCSFitResult> >(mcsTag);

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt);
    // Loop through each pfp
    for(std::vector<int>::size_type i=0; i!=(*pfpHandle).size(); i++)
    {
      art::Ptr<recob::PFParticle> pfp(pfpHandle,i);
      // Find out if this pfp is a neutrino
      if (pfp->IsPrimary())
      {
        // Fill useful variables for the tree
        etf.nNeutrinos += 1;
        etf.neutrinoPdgCode.push_back(pfp->PdgCode());
        etf.neutrinoNumDaughters.push_back(pfp->NumDaughters());
        int thisNeutrino_numTracks = 0;
        int thisNeutrino_numShowers = 0;
        // ID of current neutrino and vectors where we will fill pointers to daughters (only tracks for now)
        size_t nuID = pfp->Self();
        std::vector<art::Ptr<recob::PFParticle>> thisNeutrino_pfpTrackPointers;

        // Diagnostic message
        if (fVerbose)
        {
          printf("Neutrino %i (ID: %i, PDG: %i)\n", etf.nNeutrinos, (int) nuID, pfp->PdgCode());
          printf("|_Number of daughters: %i (ID:", pfp->NumDaughters());

          // Prepare vector of ID of neutrino daughters
          auto nuDaughtersID = pfp->Daughters();
          // Loop through each daughter and print their ID
          for (std::vector<int>::size_type j=0; j!=nuDaughtersID.size(); j++)
          {
            printf(" %i", (int) nuDaughtersID[j]);
          }
          printf(" )\n");
        }

        // Loop through each pfp and check if their parent was the neutrino we are currently looping through
        for(std::vector<int>::size_type j=0; j!=(*pfpHandle).size(); j++)
        {
          art::Ptr<recob::PFParticle> daughter_pfp(pfpHandle,j);
          // Find correct pfp corresponding to daughter we are analyzing
          if (daughter_pfp->Parent()==nuID)
          {
            // Separate in track and shower pfps and save their pointers to corresponding vectors
            if (daughter_pfp->PdgCode()==13)
            {
              if (fVerbose) printf("| |_Found track with ID: %i\n", (int) daughter_pfp->Self());
              thisNeutrino_numTracks += 1;
              thisNeutrino_pfpTrackPointers.push_back(daughter_pfp);
            }
            if (daughter_pfp->PdgCode()==11)
            {
              if (fVerbose) printf("| |_Found shower with ID: %i\n", (int) daughter_pfp->Self());
              thisNeutrino_numShowers += 1;
            }
          }
        }
        // Save to vector number of tracks and showers found
        etf.neutrinoNumTracks.push_back(thisNeutrino_numTracks);
        etf.neutrinoNumShowers.push_back(thisNeutrino_numShowers);

        // Diagnostic message
        if (fVerbose)
        {
          // Cross check, loop through each pointer, make sure their ID is correct and their parent is as well
          for (auto const& pfpTrack : thisNeutrino_pfpTrackPointers)
          {
            printf("| |_Checking saved daughter with ID %i and parent ID %i\n", (int) pfpTrack->Self(), (int) pfpTrack->Parent());
          }
        }
        // Diagnostic message
        if (fVerbose)
        {
          printf("|_Summary: %i daughters, %i tracks and %i showers.\n", pfp->NumDaughters(),thisNeutrino_numTracks, thisNeutrino_numShowers);
        }

        // If this neutrino contains two and only two tracks we can create a specific decay vertex for it (to use later for calorimetry), but first we have to make sure we have all the associations we need.
        if (thisNeutrino_numTracks==2)
        {
          etf.nTwoProngedNeutrinos += 1;
          if (fVerbose) printf("|_Neutrino is potential candidate n. %i in event.\n", etf.nTwoProngedNeutrinos);

          // Creating a stupid vector with a pointer to the current neutrino, all other methods won't work so I have to go through this stupid way of retrieving associations
          std::vector<art::Ptr<recob::PFParticle>> thisNeutrino_pfpNeutrinoPointer;
          thisNeutrino_pfpNeutrinoPointer.push_back(pfp);

          // Association to vertex and tracks objects of our pfpTracks
          art::FindOneP<recob::Vertex> nu_pva(thisNeutrino_pfpNeutrinoPointer,evt,pfpTag);
          art::FindOneP<recob::Vertex> pva(thisNeutrino_pfpTrackPointers,evt,pfpTag);
          art::FindOneP<recob::Track> pta(thisNeutrino_pfpTrackPointers,evt,pfpTag);

          // Make sure we have the necessary vertices associated with the pfp
          art::Ptr<recob::Vertex> nuVertex, t1Vertex, t2Vertex;
          nu_pva.get(0,nuVertex);
          pva.get(0,t1Vertex);
          pva.get(1,t2Vertex);
          bool rightNumVertices = (nuVertex.isNonnull() && t1Vertex.isNonnull() && t2Vertex.isNonnull());

          // Make sure we have the necessary tracks associated with the pfps
          art::Ptr<recob::Track> t1Track, t2Track;
          pta.get(0,t1Track);
          pta.get(1,t2Track);
          bool rightNumTracks = (t1Track.isNonnull() && t2Track.isNonnull());

          if (!rightNumVertices) etf.status_nuWithMissingAssociatedVertex += 1;
          if (!rightNumTracks) etf.status_nuWithMissingAssociatedTrack += 1;
          if (rightNumVertices && rightNumTracks)
          {
            if (fVerbose) printf("| | |_Neutrino has correct number of vertex and tracks associated to PFP.\n");
            std::vector<art::Ptr<recob::Track>> thisNu_tracks = {t1Track, t2Track};

            // Make sure also we have the necessary hits associated to tracks
            art::FindManyP<recob::Hit> tha(thisNu_tracks,evt,pfpTag);
            std::vector<art::Ptr<recob::Hit>> t1Hits, t2Hits;
            tha.get(0,t1Hits);
            tha.get(1,t2Hits);
            bool rightNumHits = (t1Hits.size()>1 && t2Hits.size()>1);
            std::cout << "| |_Track 1: There are " << t1Hits.size() << " associated hits." << std::endl;
            std::cout << "| |_Track 2: There are " << t2Hits.size() << " associated hits." << std::endl;

            if (!rightNumHits) etf.status_nuProngWithMissingAssociatedHits += 1;
            else
            {
              if (fVerbose) printf("| | |_Neutrino has correct number of hits vectors associated to tracks.\n");

              // For each track, find in the mcsHandle the MCS fit result with the same index (they don't have associations unfortunately but they should be paired by same index, so you can retrieve them this way).
              art::Ptr<recob::MCSFitResult> t1Mcs(mcsHandle,t1Track.key());
              art::Ptr<recob::MCSFitResult> t2Mcs(mcsHandle,t2Track.key());

              // Time to dump all associations in the neutrino vertex
              AuxVertex::DecayVertex nuV(nuVertex,t1Vertex,t2Vertex,t1Track,t2Track,t1Hits,t2Hits,t1Mcs,t2Mcs);
              nuV.SetDetectorCoordinates(fMinTpcBound,fMaxTpcBound,fGeometry,detProp);
              nuV.PrintInformation();
              if (nuV.fIsInsideTPC)
              {
                etf.nContainedTwoProngedNeutrinos += 1;
                ana_decayVertices.push_back(nuV);
              }
            } // END if right number of hits associations
          } // END if right number of vertices / tracks associations
        } // END if neutrino has 2 tracks
      } // END if pfp is a neutrino
    } // END loop for each pfp
  } // END function GetOrderedPFParticles

} // END namespace FindPandoraVertex
