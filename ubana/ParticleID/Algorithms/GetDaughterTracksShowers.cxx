#ifndef DAUGHTERTRACKSSHOWERS_CXX
#define DAUGHTERTRACKSSHOWERS_CXX

#include "GetDaughterTracksShowers.h"

// ---------------------- This function returns the number of daughter tracks from a reco track ------------------------ //

int GetNDaughterTracks(std::vector<recob::Track> trk_handle, int trkID, double cutdist, double cutfrac)
{
  // Get the track under consideration
  auto orig_track = trk_handle.at(trkID);

  // Get its end point
  double orig_endX = orig_track.End().X();
  double orig_endY = orig_track.End().Y();
  double orig_endZ = orig_track.End().Z();

  // Get track length and compare it to the cut distance
  // Actually put a cut at radius cutdist *or* x% track length, whichever is smaller
  // x% comes from cutfrac
  // This is to avoid counting tracks coming from the original vertex
  // (eg. muon produced with a short proton track)
  double orig_trklen = orig_track.Length();
  double cut_r;
  if (cutdist > cutfrac*orig_trklen){ cut_r = cutfrac*orig_trklen; }
  else{ cut_r = cutdist; }
  //std::cout << "Searching within distance " << cut_r << " cm" << std::endl;

  // int to count number of daughter tracks
  int ndaughters = 0;

  // Look through all other tracks in the event, record the number
  // with start *or* end point within cut_r of orig_trkendpoint
  for (auto const& other_trk : trk_handle){
    
    // Skip if it's this track (don't double count)
    if (other_trk.ID() == trkID){ continue; }

    // Check start point of other track
    double other_startX = other_trk.Start().X();
    double other_startY = other_trk.Start().Y();
    double other_startZ = other_trk.Start().Z();

    double startdist = pow((other_startX-orig_endX),2) + pow((other_startY-orig_endY),2) + pow((other_startZ-orig_endZ),2);
    startdist = TMath::Sqrt(startdist);

    if (startdist < cut_r){
      ndaughters++;
      continue; // don't check end -- count each track only once
    }

    // If the start point of the track isn't near the end point of our original one, check
    // the end point (in case it got flipped)
    double other_endX = other_trk.End().X();
    double other_endY = other_trk.End().Y();
    double other_endZ = other_trk.End().Z();

    double enddist = pow((other_endX-orig_endX),2) + pow((other_endY-orig_endY),2) + pow((other_endZ-orig_endZ),2);
    enddist = TMath::Sqrt(enddist);

    if (enddist < cut_r){
      ndaughters++;
      continue; // don't check end -- count each track only once
    }
   
  } // loop over other tracks

  return ndaughters;
}


// ---------------------- This function returns the number of daughter showers from a reco track ------------------------ //

int GetNDaughterShowers(std::vector<recob::Track> trk_handle, int trkID, std::vector<recob::Shower> shwr_handle, double cutdist, double cutfrac)
{
  // Get the track under consideration
  auto orig_track = trk_handle.at(trkID);

  // Get its end point
  double orig_endX = orig_track.End().X();
  double orig_endY = orig_track.End().Y();
  double orig_endZ = orig_track.End().Z();

  // Get track length and compare it to the cut distance
  // Actually put a cut at radius cutdist *or* x% track length, whichever is smaller
  // x% comes from cutfrac
  // This is to avoid counting tracks coming from the original vertex
  // (eg. muon produced with a short proton track)
  double orig_trklen = orig_track.Length();
  double cut_r;
  if (cutdist > cutfrac*orig_trklen){ cut_r = cutfrac*orig_trklen; }
  else{ cut_r = cutdist; }
  //std::cout << "Searching within distance " << cut_r << " cm" << std::endl;

  // int to count number of daughter showers
  int ndaughters = 0;

  // Look through all showers in the event, record the number
  // with start *or* end point within cut_r of orig_trkendpoint
  for (auto const& shower : shwr_handle){
    
    // Check start point of shower
    double other_startX = shower.ShowerStart().X();
    double other_startY = shower.ShowerStart().Y();
    double other_startZ = shower.ShowerStart().Z();
    
    double startdist = pow((other_startX-orig_endX),2) + pow((other_startY-orig_endY),2) + pow((other_startZ-orig_endZ),2);
    startdist = TMath::Sqrt(startdist);

    if (startdist < cut_r){
      ndaughters++;
    }
   
  } // loop over showers

  return ndaughters;
}

#endif
