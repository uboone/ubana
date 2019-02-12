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
    // Retrieve calorimetric information for each dv track
    for(int trackID=0; trackID!=int(tracks.size()); trackID++)
    {
      // Get calo/cali object associated to current track
      std::vector<art::Ptr<anab::Calorimetry>> calis = caliTrackAss.at(trackID);
      // Loop through each plane
      for(int j=0; j!=3; j++)
      {
        art::Ptr<anab::Calorimetry> cal = calis[j];
        int planeID = cal->PlaneID().Plane;
        printf("%i.%i Got plane %i \n", trackID,j, planeID);
        if ((planeID < 0) || (planeID > 2)) continue;
        std::vector<double> dqdx_v = cal->dQdx();
        int vSize = dqdx_v.size();
        float start, end;
        if (vSize>6)
        {
          start = Start_dqdx(cal);
          end = End_dqdx(cal);
        }
        else
        {
          start = -9999.;
          end = -9999.;
        }
        // Assign calorimetric quantities
        bool enoughHits = (vSize>6);
        ctf.cali_prongEnoughHits[trackID][planeID] = enoughHits;
        ctf.cali_prongCaloPlane[trackID][planeID] = planeID;
        ctf.cali_prongKinEnergy[trackID][planeID] = cal->KineticEnergy();
        ctf.cali_prongRange[trackID][planeID] = cal->Range();
        ctf.cali_prongTruncMean[trackID][planeID] = utilsHNLAlg.GetTruncatedMean(dqdx_v);
        ctf.cali_prongdqdxMeanStart[trackID][planeID] = start;
        ctf.cali_prongdqdxMeanEnd[trackID][planeID] = end;
        ctf.cali_prongIsTrackFlipped[trackID][planeID] = IsTrackFlipped(start,end);
      } // END loop for each plane
    } //  END loop for each dv track
    
    // Now do the same for uncalibrated object if required
    if (fSaveAlsoUncalibrateddEdXInformation)
    {
      art::FindManyP<anab::Calorimetry> caloTrackAss(tracks,evt,fCalorimetryLabel);
      // Retrieve calorimetric information for each dv track
      for(int trackID=0; trackID!=int(tracks.size()); trackID++)
      {
        // Get calo/calo object associated to current track
        std::vector<art::Ptr<anab::Calorimetry>> calos = caloTrackAss.at(trackID);
        // Loop through each plane
        for(int j=0; j!=3; j++)
        {
          art::Ptr<anab::Calorimetry> cal = calos[j];
          int planeID = cal->PlaneID().Plane;
          if ((planeID < 0) || (planeID > 2)) continue;
          std::vector<double> dqdx_v = cal->dQdx();
          int vSize = dqdx_v.size();
          float start, end;
          if (vSize>6)
          {
            start = Start_dqdx(cal);
            end = End_dqdx(cal);
          }
          else
          {
            start = -9999.;
            end = -9999.;
          }
          // Assign calorimetric quantities
          bool enoughHits = (vSize>6);
          ctf.calo_prongEnoughHits[trackID][planeID] = enoughHits;
          ctf.calo_prongCaloPlane[trackID][planeID] = planeID;
          ctf.calo_prongKinEnergy[trackID][planeID] = cal->KineticEnergy();
          ctf.calo_prongRange[trackID][planeID] = cal->Range();
          ctf.calo_prongTruncMean[trackID][planeID] = utilsHNLAlg.GetTruncatedMean(dqdx_v);
          ctf.calo_prongdqdxMeanStart[trackID][planeID] = start;
          ctf.calo_prongdqdxMeanEnd[trackID][planeID] = end;
          ctf.calo_prongIsTrackFlipped[trackID][planeID] = IsTrackFlipped(start,end);
        } // END loop for each plane
      } //  END loop for each dv track
    } // END if save also uncalibrated information

  } // END function AdddEdX


  //_________________________________________________________________________________
  // AUXILIARY FUNCTIONS
  //_________________________________________________________________________________
  // Determine whether track has been flipped by making sure that dqdx is larger at the end than at the start.
  float dEdXInformationAlg::Start_dqdx(art::Ptr<anab::Calorimetry> const & calo)
  {
    return ((calo->dQdx())[0] +  (calo->dQdx())[1] +  (calo->dQdx())[2])/3.;
  }
  float dEdXInformationAlg::End_dqdx(art::Ptr<anab::Calorimetry> const & calo)
  {
    return ((calo->dQdx())[-1] +  (calo->dQdx())[-2] +  (calo->dQdx())[-3])/3.;
  }
  bool dEdXInformationAlg::IsTrackFlipped(float const & start_dqdx, float const & end_dqdx)
  {
    return start_dqdx > end_dqdx;
  }

} // END namespace dEdX
