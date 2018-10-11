/******************************************************************************
 * @file DrawTreeFiller.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeFiller.cxx
 * ****************************************************************************/

#ifndef DrawTreeFiller_H
#define DrawTreeFiller_H

// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <exception>
// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraph.h"
// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
// art includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
// larsoft object includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// HSN finder includes
#include "ubana/HSNAnalysis/HsnFinder/DataObjects/DecayVertex.h"
#include "EventTreeFiller.h"


namespace AuxEvent
{

  // DrawTreeFiller class and functions
  class DrawTreeFiller
  {
  public:
    // Constructor and destructor
    DrawTreeFiller();
    virtual ~DrawTreeFiller();
    void Initialize(AuxEvent::EventTreeFiller & etf, int i_hsnID, const AuxVertex::DecayVertex & decayVertex);

    // General
    int run;
    int subrun;
    int event;
    int hsnID;
    int nHsnCandidatesInSameEvent;

    // Declare drawTree variables
    // P0
    int dv_p0_wireCoordinates;
    int dv_p0_tickCoordinates;

    int prong1_p0_wireCoordinates;
    int prong1_p0_tickCoordinates;

    int prong2_p0_wireCoordinates;
    int prong2_p0_tickCoordinates;

    std::vector<int> prong1_hits_p0_wireCoordinates;
    std::vector<float> prong1_hits_p0_tickCoordinates;
    std::vector<int> prong2_hits_p0_wireCoordinates;
    std::vector<float> prong2_hits_p0_tickCoordinates;
    std::vector<int> tot_hits_p0_wireCoordinates;
    std::vector<float> tot_hits_p0_tickCoordinates;

    // P1
    int dv_p1_wireCoordinates;
    int dv_p1_tickCoordinates;

    int prong1_p1_wireCoordinates;
    int prong1_p1_tickCoordinates;

    int prong2_p1_wireCoordinates;
    int prong2_p1_tickCoordinates;

    std::vector<int> prong1_hits_p1_wireCoordinates;
    std::vector<float> prong1_hits_p1_tickCoordinates;
    std::vector<int> prong2_hits_p1_wireCoordinates;
    std::vector<float> prong2_hits_p1_tickCoordinates;
    std::vector<int> tot_hits_p1_wireCoordinates;
    std::vector<float> tot_hits_p1_tickCoordinates;

    // P2
    int dv_p2_wireCoordinates;
    int dv_p2_tickCoordinates;

    int prong1_p2_wireCoordinates;
    int prong1_p2_tickCoordinates;

    int prong2_p2_wireCoordinates;
    int prong2_p2_tickCoordinates;

    std::vector<int> prong1_hits_p2_wireCoordinates;
    std::vector<float> prong1_hits_p2_tickCoordinates;
    std::vector<int> prong2_hits_p2_wireCoordinates;
    std::vector<float> prong2_hits_p2_tickCoordinates;
    std::vector<int> tot_hits_p2_wireCoordinates;
    std::vector<float> tot_hits_p2_tickCoordinates;

    // Edges
    float p0_maxTick, p0_minTick;
    float p1_maxTick, p1_minTick;
    float p2_maxTick, p2_minTick;
    int p0_maxWire, p0_minWire;
    int p1_maxWire, p1_minWire;
    int p2_maxWire, p2_minWire;

    // Truth information
    int truth_dv_p0_wireCoordinates;
    float truth_dv_p0_tickCoordinates;
    int truth_dv_p1_wireCoordinates;
    float truth_dv_p1_tickCoordinates;
    int truth_dv_p2_wireCoordinates;
    float truth_dv_p2_tickCoordinates;

    int truth_nPrimaryTracks, truth_nSecondaryTracks;
    std::vector<int> truth_primaryTracks_start_p0_wireCoordinates;
    std::vector<float> truth_primaryTracks_start_p0_tickCoordinates;
    std::vector<int> truth_primaryTracks_start_p1_wireCoordinates;
    std::vector<float> truth_primaryTracks_start_p1_tickCoordinates;
    std::vector<int> truth_primaryTracks_start_p2_wireCoordinates;
    std::vector<float> truth_primaryTracks_start_p2_tickCoordinates;
    std::vector<int> truth_primaryTracks_end_p0_wireCoordinates;
    std::vector<float> truth_primaryTracks_end_p0_tickCoordinates;
    std::vector<int> truth_primaryTracks_end_p1_wireCoordinates;
    std::vector<float> truth_primaryTracks_end_p1_tickCoordinates;
    std::vector<int> truth_primaryTracks_end_p2_wireCoordinates;
    std::vector<float> truth_primaryTracks_end_p2_tickCoordinates;
    std::vector<int> truth_secondaryTracks_start_p0_wireCoordinates;
    std::vector<float> truth_secondaryTracks_start_p0_tickCoordinates;
    std::vector<int> truth_secondaryTracks_start_p1_wireCoordinates;
    std::vector<float> truth_secondaryTracks_start_p1_tickCoordinates;
    std::vector<int> truth_secondaryTracks_start_p2_wireCoordinates;
    std::vector<float> truth_secondaryTracks_start_p2_tickCoordinates;
    std::vector<int> truth_secondaryTracks_end_p0_wireCoordinates;
    std::vector<float> truth_secondaryTracks_end_p0_tickCoordinates;
    std::vector<int> truth_secondaryTracks_end_p1_wireCoordinates;
    std::vector<float> truth_secondaryTracks_end_p1_tickCoordinates;
    std::vector<int> truth_secondaryTracks_end_p2_wireCoordinates;
    std::vector<float> truth_secondaryTracks_end_p2_tickCoordinates;

    int truth_nPrimaryShowers, truth_nSecondaryShowers;
    std::vector<int> truth_primaryShowers_start_p0_wireCoordinates;
    std::vector<float> truth_primaryShowers_start_p0_tickCoordinates;
    std::vector<int> truth_primaryShowers_start_p1_wireCoordinates;
    std::vector<float> truth_primaryShowers_start_p1_tickCoordinates;
    std::vector<int> truth_primaryShowers_start_p2_wireCoordinates;
    std::vector<float> truth_primaryShowers_start_p2_tickCoordinates;
    std::vector<int> truth_primaryShowers_end_p0_wireCoordinates;
    std::vector<float> truth_primaryShowers_end_p0_tickCoordinates;
    std::vector<int> truth_primaryShowers_end_p1_wireCoordinates;
    std::vector<float> truth_primaryShowers_end_p1_tickCoordinates;
    std::vector<int> truth_primaryShowers_end_p2_wireCoordinates;
    std::vector<float> truth_primaryShowers_end_p2_tickCoordinates;
    std::vector<int> truth_secondaryShowers_start_p0_wireCoordinates;
    std::vector<float> truth_secondaryShowers_start_p0_tickCoordinates;
    std::vector<int> truth_secondaryShowers_start_p1_wireCoordinates;
    std::vector<float> truth_secondaryShowers_start_p1_tickCoordinates;
    std::vector<int> truth_secondaryShowers_start_p2_wireCoordinates;
    std::vector<float> truth_secondaryShowers_start_p2_tickCoordinates;
    std::vector<int> truth_secondaryShowers_end_p0_wireCoordinates;
    std::vector<float> truth_secondaryShowers_end_p0_tickCoordinates;
    std::vector<int> truth_secondaryShowers_end_p1_wireCoordinates;
    std::vector<float> truth_secondaryShowers_end_p1_tickCoordinates;
    std::vector<int> truth_secondaryShowers_end_p2_wireCoordinates;
    std::vector<float> truth_secondaryShowers_end_p2_tickCoordinates;

  }; // END class AuxEvent
} //END namespace AuxEvent

#endif