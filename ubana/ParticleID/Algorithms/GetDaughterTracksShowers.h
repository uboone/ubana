/* --------------------------------------------------------------------------------------------------
 * Algorithm to return the number of "daughters" of a reconstsructed track (defined as tracks/showers
 * that start within some distance of the ending point of this track)
 *
 * Two variables define the distance within which you search for "daughters":
 * -- cutdist = defines a sphere of radius cutdist around the end point of your track
 * -- cutfrac = a fraction (as a decimal, e.g. 0.5 = 50%) of the track length (of the initial
 *              track, not daughters)
 * The algorithm will take whichever is smaller, cutdist or cutfrac*track_length, and search in
 * a sphere with that radius for new tracks or showers.
 *
 * GetNDaughterTracks will return the number of tracks with their start *or* end point within this
 * distance (both start and end are checked in case of flipped tracks, but each track is only counted
 * once, even if it starts and ends in this distance)
 *
 * GetNDaughterShowers will return the number of showers with their start point within this distance
 * (since the end point of a shower is not well defined)
 *
 * Other variables that you need to give the function:
 * -- std::vector<recob::Track> trk_handle: vector of all tracks in the event
 * -- int trkID: track ID for the track you want to find the daughters of
 * -- std::vector<recob::Shower> shwr_handle: vector of all showers in the event
 *
 * Kirsty Duffy (kduffy@fnal.gov), Fermilab, Jan 31 2018
 * -------------------------------------------------------------------------------------------------- */

#ifndef DAUGHTERTRACKSSHOWERS_H
#define DAUGHTERTRACKSSHOWERS_H

// art/canvas/gallery includes
#include "canvas/Persistency/Common/FindManyP.h"

// larsoft datatypes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "TMath.h"

int GetNDaughterTracks(std::vector<recob::Track> trk_handle, int trkID, float cutdist, float cutfrac);
int GetNDaughterShowers(std::vector<recob::Track> trk_handle, int trkID, std::vector<recob::Shower> shwr_handle, float cutdist, float cutfrac);

#endif
