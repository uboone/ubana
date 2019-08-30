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
struct DetectorObject {

  size_t const fid;
  size_t const foriginal_index;
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


class DetectorObjects {

	friend class SinglePhoton;
  std::map<size_t, DetectorObject *> fobject_m;
  std::vector<size_t> ftrack_index_v;
  std::vector<size_t> fshower_index_v;
  size_t fobject_id;

  std::map<size_t, size_t> foriginal_track_index_m;
  std::map<size_t, size_t> foriginal_shower_index_m;

public:

  int const ftrack_reco_type;
  int const fshower_reco_type;

  DetectorObjects();

  ~DetectorObjects() {
    for(std::pair<size_t, DetectorObject *> const & p : fobject_m) delete p.second;
  }

//  void AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t, bool const track_original_indices = false);
//  void AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s, bool const track_original_indices = false);
  void AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t, bool const track_original_indices = false);
  void AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_s, bool const track_original_indices = false);

  void SetAssociated(size_t const i);
  
  int GetRecoType(size_t const i) const;
  std::vector<size_t> const & GetTrackIndices() const {return ftrack_index_v;}
  std::vector<size_t> const & GetShowerIndices() const {return fshower_index_v;}

  DetectorObject const & GetDetectorObject(size_t const i) const;

  Track & GetTrack(size_t const i);
  Track const & GetTrack(size_t const i) const;

  Shower & GetShower(size_t const i);
  Shower const & GetShower(size_t const i) const;
  
  size_t GetTrackIndexFromOriginalIndex(size_t const i) const;
  size_t GetShowerIndexFromOriginalIndex(size_t const i) const;

};











//HEADER FILE ARE ABOVE








DetectorObjects::DetectorObjects() :
  fobject_id(0),
  ftrack_reco_type(2),
  fshower_reco_type(1){}


//void DetectorObjects::AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t,
//DetectorObjects::AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t,
//				bool const track_original_indices) {
void DetectorObjects::AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t, bool const track_original_indices){
  for(size_t i = 0; i < ev_t.size(); ++i) {
    recob::Track const & t = *(ev_t.at(i));
    fobject_m.emplace(fobject_id, new Track(fobject_id, i, ftrack_reco_type, t)); 
    ftrack_index_v.push_back(fobject_id);
    if(track_original_indices) foriginal_track_index_m.emplace(i, fobject_id);
    ++fobject_id;
  }
}

//void DetectorObjects::AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s,
void DetectorObjects::AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_s, bool const track_original_indices) {
  for(size_t i = 0; i < ev_s.size(); ++i) {
    recob::Shower const & s = *(ev_s.at(i));
    fobject_m.emplace(fobject_id, new Shower(fobject_id, i, fshower_reco_type, s));
    fshower_index_v.push_back(fobject_id);
    if(track_original_indices) foriginal_shower_index_m.emplace(i, fobject_id);
    ++fobject_id;
  }
}  


void DetectorObjects::SetAssociated(size_t const i) {

  auto om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  om_it->second->fis_associated = true;
  
}


int DetectorObjects::GetRecoType(size_t const i) const {

  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }

  return om_it->second->freco_type;

}


DetectorObject const & DetectorObjects::GetDetectorObject(size_t const i) const {

  auto const om_it = fobject_m.find(i);
  
  if(om_it == fobject_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object: " << i << std::endl;
    exit(1);
  }
  
  return *om_it->second;

}


Track & DetectorObjects::GetTrack(size_t const i) {
  
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


Track const & DetectorObjects::GetTrack(size_t const i) const {
  
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


Shower & DetectorObjects::GetShower(size_t const i) {
  
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


Shower const & DetectorObjects::GetShower(size_t const i) const {
  
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


size_t DetectorObjects::GetTrackIndexFromOriginalIndex(size_t const i) const {

  auto om_it = foriginal_track_index_m.find(i);
  
  if(om_it == foriginal_track_index_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object associated to original index: " << i << std::endl;
    exit(1);
  }

  return om_it->second;

}


size_t DetectorObjects::GetShowerIndexFromOriginalIndex(size_t const i) const {

  auto om_it = foriginal_shower_index_m.find(i);
  
  if(om_it == foriginal_shower_index_m.end()) {
    std::cout << __PRETTY_FUNCTION__ << "\nCould not find object associated to original index: " << i << std::endl;
    exit(1);
  }

  return om_it->second;

}
}
#endif
