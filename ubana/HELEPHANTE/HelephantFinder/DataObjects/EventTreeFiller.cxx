/******************************************************************************
 * @file EventTreeFiller.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  EventTreeFiller.h
 * ****************************************************************************/

// Decay vertex header
#include "EventTreeFiller.h"

namespace AuxEvent
{
  EventTreeFiller::EventTreeFiller()
  {}
  EventTreeFiller::~EventTreeFiller()
  {}

  void EventTreeFiller::Initialize(int i_run, int i_subrun, int i_event)
  {
    // Event
    run = i_run;
    subrun = i_subrun;
    event = i_event;
    nNeutrinos = -999;
    nTwoProngedNeutrinos = -999;
    nContainedTwoProngedNeutrinos = -999;
    nHsnCandidates = -999;
    neutrinoPdgCode.clear();
    neutrinoNumDaughters.clear();
    neutrinoNumTracks.clear();
    neutrinoNumShowers.clear();
    // Status
    status_nuWithMissingAssociatedVertex = -999;
    status_nuWithMissingAssociatedTrack = -999;
    status_nuProngWithMissingAssociatedHits = -999;
    // Truth distance metric (find best HSN candidate in event by proximity to truth)
    truth_vx = -999;
    truth_vy = -999;
    truth_vz = -999;
    recoTruthDistances.clear();
    isClosestToTruth.clear();
  } // END function Initialize


} // END namespace EventTreeFiller 
