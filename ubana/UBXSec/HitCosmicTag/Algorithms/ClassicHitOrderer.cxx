#ifndef CLASSICHITORDERER_CXX
#define CLASSICHITORDERER_CXX

#include "ClassicHitOrderer.h"


namespace cosmictag {

  static ClassicHitOrdererFactory __global_ClassicHitOrdererFactory__;

  ClassicHitOrderer::ClassicHitOrderer(const std::string name)
    : BaseHitOrdererAlgo(name)
  {}

  void ClassicHitOrderer::_Configure_(const Config_t &pset)
  {
    _max_allowed_hit_distance = pset.get<double>("MaxAllowedHitDistance");
    _max_allowed_hit_distance_coll_cop = pset.get<double>("MaxAllowedHitDistanceCollectionCoplanar");
    _slope_threshold          = pset.get<double>("SlopeThreshold"); 
  }

  void ClassicHitOrderer::CollectionCoplanar(bool status)
  { _collection_coplanar = status; }
  
  int ClassicHitOrderer::OrderHits(SimpleCluster& cluster) const
  {

    int                    & _start_index  = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v      = cluster._s_hit_v;
    std::vector<double>    & _ds_v         = cluster._ds_v;
    bool                   & _hits_ordered = cluster._hits_ordered;

    _hits_ordered = false;

    double max_allowed_hit_distance = _max_allowed_hit_distance;

    if (_collection_coplanar) {
      CT_DEBUG() << "Collection Coplanar" << std::endl;
      max_allowed_hit_distance = _max_allowed_hit_distance_coll_cop;
    }

    if (_start_index < 0) {
      CT_NORMAL() << "Start hit not set." << std::endl;
      return 0;
    }

    if ((size_t)_start_index >= _s_hit_v.size()) {
      CT_NORMAL() << "Start hit not compatible with current hit vector." 
       << "Start hit index is " << _start_index 
       << ", while hit vector size is " << _s_hit_v.size() << std::endl;
      return 0;
    }


    std::vector<SimpleHit> new_vector;
    new_vector.clear();
    new_vector.reserve(_s_hit_v.size());

    _ds_v.clear();
    _ds_v.reserve(_s_hit_v.size());

    new_vector.push_back(_s_hit_v.at(_start_index));

    //  for (auto h : _s_hit_v) {
    //    std::cout << "BEFORE: " << h.wire << ", " << h.time*4 << std::endl;
    //  }

    _s_hit_v.erase(_s_hit_v.begin() + _start_index);

    double min_dist = 1e9; 
    int min_index = -1;

    while (_s_hit_v.size() != 0) {

      min_dist = 1e9;
      min_index = -1; 

      for (size_t i = 0; i < _s_hit_v.size(); i++){

        TVector3 pt1(new_vector.back().time, new_vector.back().wire, 0);
        TVector3 pt2(_s_hit_v.at(i).time,    _s_hit_v.at(i).wire,    0);
        double dist = (pt1 - pt2).Mag();

        if (dist < min_dist) {
          min_index = i;
          min_dist = dist;
        }
      }

      if (min_index < 0) {
         CT_CRITICAL() << "Logic fail." << std::endl;
         throw std::exception();
      }

      // Emplace the next hit in the new vector...
      if (min_dist < max_allowed_hit_distance) {
        CT_DEBUG()  << "min_dist: " << min_dist <<std::endl;
        new_vector.push_back(_s_hit_v.at(min_index));
        _ds_v.push_back(min_dist);
      } else if (_s_hit_v.at(min_index).wire == new_vector.back().wire && min_dist < 50) {
        CT_DEBUG()  << "min_dist: " << min_dist << " => but on same wire, so continue"<< std::endl;
        new_vector.push_back(_s_hit_v.at(min_index));
        _ds_v.push_back(min_dist);
      } else if (new_vector.size() > 5){
 
        // Calculate previous slope
        auto iter = new_vector.end();
        auto sh_3 = _s_hit_v.at(min_index);
        auto sh_2 = *(--iter);
        auto sh_1 = *(iter-5); // go 5 hits back
        double slope = (sh_2.time - sh_1.time) / (sh_2.wire - sh_1.wire);

        // Calculate the new slope
        double slope_new = (sh_3.time - sh_2.time) / (sh_3.wire - sh_2.wire);

        CT_DEBUG()  << "sh_1.wire: "<<sh_1.wire<<", sh_1.time: "<< sh_1.time << std::endl;
        CT_DEBUG()  << "sh_2.wire: "<<sh_2.wire<<", sh_2.time: "<< sh_2.time << std::endl;
        CT_DEBUG()  << "sh_3.wire: "<<sh_3.wire<<", sh_3.time: "<< sh_3.time << std::endl;
        CT_DEBUG()  << "Current slope : " << slope 
                    << " New slope: " << slope_new 
                    << " Diff: " << slope_new - slope << std::endl;

        // Check the next hit will be in a consecutive wire
        bool progressive_order = false;

        if (sh_1.wire < sh_2.wire) {
          if (sh_3.wire > sh_2.wire) {
            progressive_order = true;
          }
        }
        if (sh_2.wire < sh_1.wire) {
          if (sh_3.wire < sh_2.wire) {
            progressive_order = true;
          }
        }

        CT_DEBUG() << "Progressive order? " << (progressive_order ? "YES" : "NO") << std::endl;

        // If the two slopes are close, than there is 
        // probably a dead region between the point.
        // If so, increase the min distance by half a meter
        // and add the hit.
        if (std::abs(slope_new - slope) < _slope_threshold &&
            min_dist < max_allowed_hit_distance + 50 &&
            progressive_order) {

          new_vector.push_back(_s_hit_v.at(min_index)); 
          _ds_v.push_back(min_dist);

        } 

      }

      // ...and delete it from the old vector
      _s_hit_v.erase(_s_hit_v.begin() + min_index);

    }

    // For the last hit, use the last values of min_dist
    _ds_v.push_back(min_dist);

    // Now that the vector is ordered, reassing to the original one
    _s_hit_v = new_vector;

    //for (auto h : _s_hit_v) {
      //std::cout << "Ordered hits: " << h.wire << ", " << h.time*4 << std::endl;
    //}

    _hits_ordered = true;

    return _s_hit_v.size();
  }
}
#endif
