#ifndef DETECTOROBJECTS_H
#define DETECTOROBJECTS_H

#include <map>

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "art/Framework/Principal/Handle.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoVector.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoSphere.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoCone.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"

#include "../SinglePhoton_module.h"

namespace single_photon{
//DetectorObject is a container for labeling tracks and showers with IDs, and for each track/shower, the corresponding details are also stored.
struct DetectorObject {

  size_t const fid;//track/shower IDs
  size_t const foriginal_index;//PFParticle IDs, IDToPFParticleMap[fid].first?;
  int const freco_type;
  bool fis_associated;//true when an object is associated.
  
  DetectorObject(size_t const id, size_t const original_index, int const reco_type) :
    fid(id),
    foriginal_index(original_index),
    freco_type(reco_type),
    fis_associated(false) {}
  
  virtual ~DetectorObject(){}
  
};


struct Track : public DetectorObject {
  
  geoalgo::Trajectory ftrajectory;
  
  Track(size_t const id, size_t const original_index, int const reco_type, recob::Track const & t) :
    DetectorObject(id, original_index, reco_type) {
    for(size_t i = 0; i < t.NumberTrajectoryPoints(); ++i)
      ftrajectory.push_back(t.LocationAtPoint<TVector3>(i));
  }
  
};


struct Shower : public DetectorObject {
  
  geoalgo::Cone fcone;
  
  Shower(size_t const id, size_t const original_index, int const reco_type, recob::Shower const & s) :
    DetectorObject(id, original_index, reco_type) {
    fcone = geoalgo::Cone(s.ShowerStart(),
			  s.Direction(),
			  s.Length(),
			  0);
  }
  
};

//collectinos of the DetectorObject??
class DetectorObjects_all {

	friend class SinglePhoton;
  std::map<size_t, DetectorObject *> fobject_m;//index , DetectorObject
  std::vector<size_t> ftrack_index_v;
  std::vector<size_t> fshower_index_v;
  size_t fobject_id;

//  std::map<size_t, size_t> foriginal_track_index_m;
//  std::map<size_t, size_t> foriginal_shower_index_m;

public:

  int const ftrack_reco_type;
  int const fshower_reco_type;

  DetectorObjects_all();

  ~DetectorObjects_all() {
    for(std::pair<size_t, DetectorObject *> const & p : fobject_m) delete p.second;
  }

//  void AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t, bool const track_original_indices = false);
//  void AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s, bool const track_original_indices = false);
  void AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t, bool const track_original_indices = false);
  void AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_s, bool const track_original_indices = false);

  void SetAssociated(size_t const i);
  
  /*******
   *
   * Return the freco_type of an fobject_m whose index is i.
   *
   *	fobject_m is a map <index, DetectorObject>
   *	freco_type - 2 for track, and 1 for shower;
   *
   *
   * ******/
  int GetRecoType(size_t const i) const;
  std::vector<size_t> const & GetTrackIndices() const {return ftrack_index_v;}
  std::vector<size_t> const & GetShowerIndices() const {return fshower_index_v;}

  DetectorObject const & GetDetectorObject(size_t const i) const;

  Track & GetTrack(size_t const i);
  Track const & GetTrack(size_t const i) const;

  Shower & GetShower(size_t const i);
  Shower const & GetShower(size_t const i) const;
  
//  size_t GetTrackIndexFromOriginalIndex(size_t const i) const;
//  size_t GetShowerIndexFromOriginalIndex(size_t const i) const;

};











//THE HEADER FILE IS THE ABOVE








DetectorObjects_all::DetectorObjects_all() :
  fobject_id(0),
  ftrack_reco_type(2),
  fshower_reco_type(1){}


//void DetectorObjects_all::AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t,
//DetectorObjects_all::AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t,
//				bool const track_original_indices) {
void DetectorObjects_all::AddTracks(
	std::vector<art::Ptr<recob::Track>> const & ev_t, 
	bool const track_original_indices){
  for(size_t i = 0; i < ev_t.size(); ++i) {
    recob::Track const & t = *(ev_t.at(i));
    fobject_m.emplace(fobject_id, new Track(fobject_id, i, ftrack_reco_type, t)); 
    ftrack_index_v.push_back(fobject_id);
//    if(track_original_indices) foriginal_track_index_m.emplace(i, fobject_id);
    ++fobject_id;
  }
}

//void DetectorObjects_all::AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s,
void DetectorObjects_all::AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_s, bool const track_original_indices) {
  for(size_t i = 0; i < ev_s.size(); ++i) {
    recob::Shower const & s = *(ev_s.at(i));
    fobject_m.emplace(fobject_id, new Shower(fobject_id, i, fshower_reco_type, s));
    fshower_index_v.push_back(fobject_id);
//    if(track_original_indices) foriginal_shower_index_m.emplace(i, fobject_id);
    ++fobject_id;
  }
}  


void DetectorObjects_all::SetAssociated(size_t const i) {

  auto om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  om_it->second->fis_associated = true;
  
}


int DetectorObjects_all::GetRecoType(size_t const i) const {

  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }

  return om_it->second->freco_type;

}


DetectorObject const & DetectorObjects_all::GetDetectorObject(size_t const i) const {

  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  return *om_it->second;

}


Track & DetectorObjects_all::GetTrack(size_t const i) {
  
  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  Track * t = dynamic_cast<Track *>(fobject_m.find(i)->second);
  
  if(!t) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not convert: " << i << std::endl;
    exit(1);
  }
  
  return *t;

}


Track const & DetectorObjects_all::GetTrack(size_t const i) const {
  
  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  Track * t = dynamic_cast<Track *>(fobject_m.find(i)->second);
  
  if(!t) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not convert: " << i << std::endl;
    exit(1);
  }
  
  return *t;

}


Shower & DetectorObjects_all::GetShower(size_t const i) {
  
  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  Shower * s = dynamic_cast<Shower *>(fobject_m.find(i)->second);
  
  if(!s) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find convert: " << i << std::endl;
    exit(1);
  }
  
  return *s;

}


Shower const & DetectorObjects_all::GetShower(size_t const i) const {
  
  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  Shower * s = dynamic_cast<Shower *>(fobject_m.find(i)->second);
  
  if(!s) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find convert: " << i << std::endl;
    exit(1);
  }
  
  return *s;

}


//size_t DetectorObjects_all::GetTrackIndexFromOriginalIndex(size_t const i) const {
//
//  auto om_it = foriginal_track_index_m.find(i);
//  
//  if(om_it == foriginal_track_index_m.end()) {
//    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object associated to original index: " << i << std::endl;
//    exit(1);
//  }
//
//  return om_it->second;
//
//}
//
//
//size_t DetectorObjects_all::GetShowerIndexFromOriginalIndex(size_t const i) const {
//
//  auto om_it = foriginal_shower_index_m.find(i);
//  
//  if(om_it == foriginal_shower_index_m.end()) {
//    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object associated to original index: " << i << std::endl;
//    exit(1);
//  }
//
//  return om_it->second;
//
//}
}
#endif
