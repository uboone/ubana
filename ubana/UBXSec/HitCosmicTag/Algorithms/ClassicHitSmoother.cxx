#ifndef CLASSICHitSmoother_CXX
#define CLASSICHitSmoother_CXX

#include "ClassicHitSmoother.h"


namespace cosmictag {

  static ClassicHitSmootherFactory __global_ClassicHitSmootherFactory__;

  ClassicHitSmoother::ClassicHitSmoother(const std::string name)
    : BaseHitSmootherAlgo(name)
  {}

  void ClassicHitSmoother::_Configure_(const Config_t &pset)
  {
     _slider_window = pset.get<size_t>("SliderWindow", 4);
  }
  
  int ClassicHitSmoother::Smooth(SimpleCluster& cluster) const
  {

    std::vector<SimpleHit> & _s_hit_v      = cluster._s_hit_v;
    std::vector<double>    & _ds_v         = cluster._ds_v;
    bool                   & _hits_ordered = cluster._hits_ordered;

    if (!_hits_ordered) {
      CT_CRITICAL() << ": need to order hits first." << std::endl;
      throw HitCosmicTagException();
    }

    if (_s_hit_v.size() < _slider_window) {
      CT_NORMAL() << ": _s_hit_v size is less than " << _slider_window << "." << std::endl;
      return 0;
    }

    std::vector<SimpleHit> new_vector;
    new_vector.clear(); 

    std::vector<double> new_vector_ds;
    new_vector_ds.clear();

    std::vector<double> mean_v;
    mean_v.clear();

    std::vector<double> wire_v;

    for(const auto& window : get_windows(_s_hit_v, _slider_window)) {

      for (auto s_h : window) {
        wire_v.push_back(s_h.wire);
      }

      mean_v.push_back(mean(wire_v));

      wire_v.clear();
    }

    new_vector.push_back(_s_hit_v.at(0));
    new_vector.push_back(_s_hit_v.at(1));
    new_vector_ds.push_back(_ds_v.at(0));
    new_vector_ds.push_back(_ds_v.at(1));

    for (size_t i = 2; i < mean_v.size()-1; i++) {

      //if (_debug) std::cout << "i: " << i 
      //                      << " wire : " << _s_hit_v.at(i).wire 
      //                      << " time : " << _s_hit_v.at(i).time 
      //                      << " mean: " << mean_v.at(i) << std::endl;

      if (std::abs(mean_v.at(i-1) - mean_v.at(i)) < 1    && 
          _s_hit_v.at(i-1).wire !=  _s_hit_v.at(i).wire  &&
          std::abs(mean_v.at(i)   - mean_v.at(i+1)) < 1  &&
          _s_hit_v.at(i).wire !=  _s_hit_v.at(i+1).wire    ) {

        //std::cout << ">>>" << std::endl;
        if (_s_hit_v.at(i).integral > _s_hit_v.at(i+1).integral) {
          new_vector.push_back(_s_hit_v.at(i));
          new_vector_ds.push_back(_ds_v.at(i));
        } else {
          new_vector.push_back(_s_hit_v.at(i+1));
          new_vector_ds.push_back(_ds_v.at(i+1)); 
        }

        i++;

      } else {

        new_vector.push_back(_s_hit_v.at(i));
        new_vector_ds.push_back(_ds_v.at(i));
      }

    }

    new_vector.push_back(_s_hit_v.at(_s_hit_v.size()-1));
    new_vector_ds.push_back(_ds_v.at(_ds_v.size()-1)); 

    std::swap(new_vector, _s_hit_v);
    std::swap(new_vector_ds, _ds_v); 

    if (_s_hit_v.size() != _ds_v.size()) {
      CT_CRITICAL() << "_s_hit_v and _ds_v size mismatch" << std::endl;
      throw std::exception();
    }

    return _s_hit_v.size();
  }
}
#endif
