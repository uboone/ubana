/******************************************************************************
 * @file EventTreeFiller.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  EventTreeFiller.cxx
 * ****************************************************************************/

#ifndef EventTreeFiller_H
#define EventTreeFiller_H

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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
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
#include "ubana/HELEPHANT/HelephantFinder/DataObjects/DecayVertex.h"

namespace AuxEvent
{

  // EventTreeFiller class and functions
  class EventTreeFiller
  {
  public:
    // Constructor and destructor
    EventTreeFiller();
    virtual ~EventTreeFiller();

    void Initialize(int i_run, int i_subrun, int i_event);

    // General
    int run;
    int subrun;
    int event;
    int hsnID;

    // Event variables search
    int nNeutrinos;
    std::vector<int> neutrinoPdgCode;
    std::vector<int> neutrinoNumDaughters;
    std::vector<int> neutrinoNumTracks;
    std::vector<int> neutrinoNumShowers;
    int nTwoProngedNeutrinos;
    int nContainedTwoProngedNeutrinos;
    int nHsnCandidates;
    // Status
    int status_nuWithMissingAssociatedVertex;
    int status_nuWithMissingAssociatedTrack;
    int status_nuProngWithMissingAssociatedHits;

    // Truth distance metric (find best HSN candidate in event by proximity to truth)
    float truth_vx, truth_vy, truth_vz;
    std::vector<float> recoTruthDistances;
    std::vector<bool> isClosestToTruth;
  };


} //END namespace AuxEvent

#endif
