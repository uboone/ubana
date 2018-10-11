/******************************************************************************
 * @file DrawTreeFiller.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeFiller.h
 * ****************************************************************************/

// Decay vertex header
#include "DrawTreeFiller.h"

namespace AuxEvent
{
  DrawTreeFiller::DrawTreeFiller()
  {}
  DrawTreeFiller::~DrawTreeFiller()
  {}

  void DrawTreeFiller::Initialize(AuxEvent::EventTreeFiller & etf, int i_hsnID, const AuxVertex::DecayVertex & decayVertex)
  {
    // General
    run = etf.run;
    subrun = etf.subrun;
    event = etf.event;
    nHsnCandidatesInSameEvent = etf.nHsnCandidates;
    hsnID = i_hsnID;

    // Drawing edges
    p0_maxTick = -1e10;
    p0_minTick = 1e10;
    p1_maxTick = -1e10;
    p1_minTick = 1e10;
    p2_maxTick = -1e10;
    p2_minTick = 1e10;
    p0_maxWire = -1e6;
    p0_minWire = 1e6;
    p1_maxWire = -1e6;
    p1_minWire = 1e6;
    p2_maxWire = -1e6;
    p2_minWire = 1e6;

    // Get vector of hits
    std::vector<art::Ptr<recob::Hit>> prong1_hits = decayVertex.GetProngHits(0);
    std::vector<art::Ptr<recob::Hit>> prong2_hits = decayVertex.GetProngHits(1);
    std::vector<art::Ptr<recob::Hit>> thisTot_hits = decayVertex.GetTotHits();

    // Fill prong1
    prong1_hits_p0_wireCoordinates.clear();
    prong1_hits_p0_tickCoordinates.clear();
    prong1_hits_p1_wireCoordinates.clear();
    prong1_hits_p1_tickCoordinates.clear();
    prong1_hits_p2_wireCoordinates.clear();
    prong1_hits_p2_tickCoordinates.clear();
    for (auto hit : prong1_hits)
    {
      float meanTick = (hit->StartTick() + hit->EndTick())/2.;
      int channel = hit->Channel();
      if (hit->View() == 0) {
        prong1_hits_p0_wireCoordinates.push_back(channel);
        prong1_hits_p0_tickCoordinates.push_back(meanTick);
        if (channel<p0_minWire) p0_minWire = channel;
        if (channel>p0_maxWire) p0_maxWire = channel;
        if (meanTick<p0_minTick) p0_minTick = meanTick;
        if (meanTick>p0_maxTick) p0_maxTick = meanTick;
      }
      if (hit->View() == 1) {
        prong1_hits_p1_wireCoordinates.push_back(channel);
        prong1_hits_p1_tickCoordinates.push_back(meanTick);
        if (channel<p1_minWire) p1_minWire = channel;
        if (channel>p1_maxWire) p1_maxWire = channel;
        if (meanTick<p1_minTick) p1_minTick = meanTick;
        if (meanTick>p1_maxTick) p1_maxTick = meanTick;
      }
      if (hit->View() == 2) {
        prong1_hits_p2_wireCoordinates.push_back(channel);
        prong1_hits_p2_tickCoordinates.push_back(meanTick);
        if (channel<p2_minWire) p2_minWire = channel;
        if (channel>p2_maxWire) p2_maxWire = channel;
        if (meanTick<p2_minTick) p2_minTick = meanTick;
        if (meanTick>p2_maxTick) p2_maxTick = meanTick;
      }
    }

    // Fill prong2
    prong2_hits_p0_wireCoordinates.clear();
    prong2_hits_p0_tickCoordinates.clear();
    prong2_hits_p1_wireCoordinates.clear();
    prong2_hits_p1_tickCoordinates.clear();
    prong2_hits_p2_wireCoordinates.clear();
    prong2_hits_p2_tickCoordinates.clear();
    for (auto hit : prong2_hits)
    {
      float meanTick = (hit->StartTick() + hit->EndTick())/2.;
      int channel = hit->Channel();
      if (hit->View() == 0) {
        prong2_hits_p0_wireCoordinates.push_back(channel);
        prong2_hits_p0_tickCoordinates.push_back(meanTick);
        if (channel<p0_minWire) p0_minWire = channel;
        if (channel>p0_maxWire) p0_maxWire = channel;
        if (meanTick<p0_minTick) p0_minTick = meanTick;
        if (meanTick>p0_maxTick) p0_maxTick = meanTick;
      }
      if (hit->View() == 1) {
        prong2_hits_p1_wireCoordinates.push_back(channel);
        prong2_hits_p1_tickCoordinates.push_back(meanTick);
        if (channel<p1_minWire) p1_minWire = channel;
        if (channel>p1_maxWire) p1_maxWire = channel;
        if (meanTick<p1_minTick) p1_minTick = meanTick;
        if (meanTick>p1_maxTick) p1_maxTick = meanTick;
      }
      if (hit->View() == 2) {
        prong2_hits_p2_wireCoordinates.push_back(channel);
        prong2_hits_p2_tickCoordinates.push_back(meanTick);
        if (channel<p2_minWire) p2_minWire = channel;
        if (channel>p2_maxWire) p2_maxWire = channel;
        if (meanTick<p2_minTick) p2_minTick = meanTick;
        if (meanTick>p2_maxTick) p2_maxTick = meanTick;
      }
    }

    // Fill totHits
    tot_hits_p0_wireCoordinates.clear();
    tot_hits_p0_tickCoordinates.clear();
    tot_hits_p1_wireCoordinates.clear();
    tot_hits_p1_tickCoordinates.clear();
    tot_hits_p2_wireCoordinates.clear();
    tot_hits_p2_tickCoordinates.clear();
    for (auto hit : thisTot_hits)
    {
      float meanTick = (hit->StartTick() + hit->EndTick())/2.;
      int channel = hit->Channel();
      if (hit->View() == 0) {
        tot_hits_p0_wireCoordinates.push_back(channel);
        tot_hits_p0_tickCoordinates.push_back(meanTick);
        if (channel<p0_minWire) p0_minWire = channel;
        if (channel>p0_maxWire) p0_maxWire = channel;
        if (meanTick<p0_minTick) p0_minTick = meanTick;
        if (meanTick>p0_maxTick) p0_maxTick = meanTick;
      }
      if (hit->View() == 1) {
        tot_hits_p1_wireCoordinates.push_back(channel);
        tot_hits_p1_tickCoordinates.push_back(meanTick);
        if (channel<p1_minWire) p1_minWire = channel;
        if (channel>p1_maxWire) p1_maxWire = channel;
        if (meanTick<p1_minTick) p1_minTick = meanTick;
        if (meanTick>p1_maxTick) p1_maxTick = meanTick;
      }
      if (hit->View() == 2) {
        tot_hits_p2_wireCoordinates.push_back(channel);
        tot_hits_p2_tickCoordinates.push_back(meanTick);
        if (channel<p2_minWire) p2_minWire = channel;
        if (channel>p2_maxWire) p2_maxWire = channel;
        if (meanTick<p2_minTick) p2_minTick = meanTick;
        if (meanTick>p2_maxTick) p2_maxTick = meanTick;
      }
    }

    // Get coordinates
    dv_p0_wireCoordinates = decayVertex.fChannelLoc[0];
    dv_p0_tickCoordinates = decayVertex.fTickLoc[0];
    dv_p1_wireCoordinates = decayVertex.fChannelLoc[1];
    dv_p1_tickCoordinates = decayVertex.fTickLoc[1];
    dv_p2_wireCoordinates = decayVertex.fChannelLoc[2];
    dv_p2_tickCoordinates = decayVertex.fTickLoc[2];
    if (dv_p0_wireCoordinates<p0_minWire) p0_minWire = dv_p0_wireCoordinates;
    if (dv_p0_wireCoordinates>p0_maxWire) p0_maxWire = dv_p0_wireCoordinates;
    if (dv_p1_wireCoordinates<p1_minWire) p1_minWire = dv_p1_wireCoordinates;
    if (dv_p1_wireCoordinates>p1_maxWire) p1_maxWire = dv_p1_wireCoordinates;
    if (dv_p2_wireCoordinates<p2_minWire) p2_minWire = dv_p2_wireCoordinates;
    if (dv_p2_wireCoordinates>p2_maxWire) p2_maxWire = dv_p2_wireCoordinates;
    if (dv_p0_tickCoordinates<p0_minTick) p0_minTick = dv_p0_tickCoordinates;
    if (dv_p0_tickCoordinates>p0_maxTick) p0_maxTick = dv_p0_tickCoordinates;
    if (dv_p1_tickCoordinates<p1_minTick) p1_minTick = dv_p1_tickCoordinates;
    if (dv_p1_tickCoordinates>p1_maxTick) p1_maxTick = dv_p1_tickCoordinates;
    if (dv_p2_tickCoordinates<p2_minTick) p2_minTick = dv_p2_tickCoordinates;
    if (dv_p2_tickCoordinates>p2_maxTick) p2_maxTick = dv_p2_tickCoordinates;

    prong1_p0_wireCoordinates = decayVertex.fProngChannelLoc[0][0];
    prong1_p0_tickCoordinates = decayVertex.fProngTickLoc[0][0];
    prong1_p1_wireCoordinates = decayVertex.fProngChannelLoc[0][1];
    prong1_p1_tickCoordinates = decayVertex.fProngTickLoc[0][1];
    prong1_p2_wireCoordinates = decayVertex.fProngChannelLoc[0][2];
    prong1_p2_tickCoordinates = decayVertex.fProngTickLoc[0][2];
    if (prong1_p0_wireCoordinates<p0_minWire) p0_minWire = prong1_p0_wireCoordinates;
    if (prong1_p0_wireCoordinates>p0_maxWire) p0_maxWire = prong1_p0_wireCoordinates;
    if (prong1_p1_wireCoordinates<p1_minWire) p1_minWire = prong1_p1_wireCoordinates;
    if (prong1_p1_wireCoordinates>p1_maxWire) p1_maxWire = prong1_p1_wireCoordinates;
    if (prong1_p2_wireCoordinates<p2_minWire) p2_minWire = prong1_p2_wireCoordinates;
    if (prong1_p2_wireCoordinates>p2_maxWire) p2_maxWire = prong1_p2_wireCoordinates;
    if (prong1_p0_tickCoordinates<p0_minTick) p0_minTick = prong1_p0_tickCoordinates;
    if (prong1_p0_tickCoordinates>p0_maxTick) p0_maxTick = prong1_p0_tickCoordinates;
    if (prong1_p1_tickCoordinates<p1_minTick) p1_minTick = prong1_p1_tickCoordinates;
    if (prong1_p1_tickCoordinates>p1_maxTick) p1_maxTick = prong1_p1_tickCoordinates;
    if (prong1_p2_tickCoordinates<p2_minTick) p2_minTick = prong1_p2_tickCoordinates;
    if (prong1_p2_tickCoordinates>p2_maxTick) p2_maxTick = prong1_p2_tickCoordinates;

    prong2_p0_wireCoordinates = decayVertex.fProngChannelLoc[1][0];
    prong2_p0_tickCoordinates = decayVertex.fProngTickLoc[1][0];
    prong2_p1_wireCoordinates = decayVertex.fProngChannelLoc[1][1];
    prong2_p1_tickCoordinates = decayVertex.fProngTickLoc[1][1];
    prong2_p2_wireCoordinates = decayVertex.fProngChannelLoc[1][2];
    prong2_p2_tickCoordinates = decayVertex.fProngTickLoc[1][2];
    if (prong2_p0_wireCoordinates<p0_minWire) p0_minWire = prong2_p0_wireCoordinates;
    if (prong2_p0_wireCoordinates>p0_maxWire) p0_maxWire = prong2_p0_wireCoordinates;
    if (prong2_p1_wireCoordinates<p1_minWire) p1_minWire = prong2_p1_wireCoordinates;
    if (prong2_p1_wireCoordinates>p1_maxWire) p1_maxWire = prong2_p1_wireCoordinates;
    if (prong2_p2_wireCoordinates<p2_minWire) p2_minWire = prong2_p2_wireCoordinates;
    if (prong2_p2_wireCoordinates>p2_maxWire) p2_maxWire = prong2_p2_wireCoordinates;
    if (prong2_p0_tickCoordinates<p0_minTick) p0_minTick = prong2_p0_tickCoordinates;
    if (prong2_p0_tickCoordinates>p0_maxTick) p0_maxTick = prong2_p0_tickCoordinates;
    if (prong2_p1_tickCoordinates<p1_minTick) p1_minTick = prong2_p1_tickCoordinates;
    if (prong2_p1_tickCoordinates>p1_maxTick) p1_maxTick = prong2_p1_tickCoordinates;
    if (prong2_p2_tickCoordinates<p2_minTick) p2_minTick = prong2_p2_tickCoordinates;
    if (prong2_p2_tickCoordinates>p2_maxTick) p2_maxTick = prong2_p2_tickCoordinates;

  } // END function initialize
} // END namespace DrawTreeFiller 
