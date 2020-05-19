#include "RecoTruthDistanceAlg.h"

namespace RecoTruthDistance
{
  // Constructor/destructor
  RecoTruthDistanceAlg::RecoTruthDistanceAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
  }
  RecoTruthDistanceAlg::~RecoTruthDistanceAlg()
  {}
  void RecoTruthDistanceAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fVerbose = pset.get<bool>("VerboseMode");
  }


  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void RecoTruthDistanceAlg::DetermineRecoTruthDistance(
            art::Event const & evt,
            AuxEvent::EventTreeFiller & etf,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices)
  {
    // Prepare handle labels
    std::string mcTruthLabel = "generator";
    art::InputTag mcTruthTag {mcTruthLabel};
    const auto& mcTruthHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTruthTag);
    // Convention valid only for HSN! (first mcParticle in first mcTruth is either pi or mu from decay).
    art::Ptr<simb::MCTruth> mcTruth(mcTruthHandle,0);
    const simb::MCParticle & mcPart = mcTruth->GetParticle(0);
    etf.truth_vx = mcPart.Vx();
    etf.truth_vy = mcPart.Vy();
    etf.truth_vz = mcPart.Vz();

    // Loop through each decay vertex and find reco-truth distance
    float minDist = 1e10;
    int minDistInd = -1;
    for (std::vector<int>::size_type i=0; i!=ana_decayVertices.size(); i++)
    {
      AuxVertex::DecayVertex dv = ana_decayVertices[i];
      float distance = sqrt( pow((etf.truth_vx - dv.fX),2.) + pow((etf.truth_vy - dv.fY),2.) + pow((etf.truth_vz - dv.fZ),2.) );
      if ( distance < minDist )
      {
        minDist = distance;
        minDistInd = i;
      }
      etf.recoTruthDistances.push_back(distance);
      etf.isClosestToTruth.push_back(0);
    }
    etf.isClosestToTruth[minDistInd] = 1;
  } // END function DetermineRecoTruthDistance

} // END namespace RecoTruthDistance
