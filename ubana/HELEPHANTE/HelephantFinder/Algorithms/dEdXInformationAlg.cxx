/******************************************************************************
 * @file dEdXInformationAlg.cxx
 * @brief Retrieve dEdX information regarding the current vertex in order to perform calorimetry.
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see dEdXInformationAlg.h
 * ****************************************************************************/
#include "dEdXInformationAlg.h"

namespace dEdXInformation
{
  // Constructor/destructor
  dEdXInformationAlg::dEdXInformationAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
  }
  dEdXInformationAlg::~dEdXInformationAlg(){}
  void dEdXInformationAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCalorimetryLabel = pset.get<std::string>("CalorimetryLabel");
    fCalibratedLabel = pset.get<std::string>("CalibratedLabel");
    fSaveAlsoUncalibrateddEdXInformation = pset.get<bool>("SaveAlsoUncalibrateddEdXInformation");
  }

  // Determine whether track has been flipped by making sure that dqdx is larger at the end than at the start.
  bool dEdXInformationAlg::IsTrackFlipped(art::Ptr<anab::Calorimetry> const & calo)
  {
    float dqdx_start = (calo->dQdx())[0] +  (calo->dQdx())[0] +  (calo->dQdx())[0];
    float dqdx_end = (calo->dQdx())[-1] +  (calo->dQdx())[-2] +  (calo->dQdx())[-3];
    bool isTrackFlipped = dqdx_start > dqdx_end;
    return isTrackFlipped;
  }

  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void dEdXInformationAlg::AdddEdXInformation(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf,
            AuxVertex::DecayVertex const & dv)
  {

    // Retrieve tracks from decay vertex and associated calorimetric objects
    std::vector<art::Ptr<recob::Track>> tracks;
    tracks.push_back(dv.GetProngTrack(0));
    tracks.push_back(dv.GetProngTrack(1));
    art::FindManyP<anab::Calorimetry> caliTrackAss(tracks,evt,fCalibratedLabel);
    for(int trackID=0; trackID!=int(tracks.size()); trackID++)
    {
      // art::Ptr<recob::Track> track = tracks[i];
      std::vector<art::Ptr<anab::Calorimetry>> calis = caliTrackAss.at(trackID);
      for(int j=0; j!=3; j++)
      {
        art::Ptr<anab::Calorimetry> cali = calis[j];
        int planeID = cali->PlaneID().Plane;
        if ((planeID < 0) || (planeID > 2)) continue;
        ctf.cali_prongCaloPlane[trackID][planeID] = planeID;
        ctf.cali_prongIsTrackFlipped[trackID][planeID] = IsTrackFlipped(cali);
      }

    } //  END loop for each dv track
  } // END function AdddEdX



} // END namespace dEdX
