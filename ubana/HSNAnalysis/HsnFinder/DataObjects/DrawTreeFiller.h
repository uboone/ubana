/******************************************************************************
 * @file DrawTreeFiller.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeFiller.cxx
 * ****************************************************************************/

#ifndef DrawTreeFiller_H
#define DrawTreeFiller_H

// C++ standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
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
    // Geometric coordinates
    std::vector<float> dv_xyzCoordinates; 
    std::vector<float> prong1_xyzCoordinates; 
    std::vector<float> prong2_xyzCoordinates; 

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
  }; // END class AuxEvent
} //END namespace AuxEvent

#endif