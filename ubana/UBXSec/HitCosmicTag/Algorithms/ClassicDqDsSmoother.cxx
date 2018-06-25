#ifndef CLASSICSDQDSSMOOTHER_CXX
#define CLASSICSDQDSSMOOTHER_CXX

#include "ClassicDqDsSmoother.h"

namespace cosmictag {

  static ClassicDqDsSmootherFactory __global_ClassicDqDsSmootherFactory__;

  ClassicDqDsSmoother::ClassicDqDsSmoother(const std::string name)
    : BaseDqDsSmootherAlgo(name)
  {}

  void ClassicDqDsSmoother::_Configure_(const Config_t &pset)
  {
    _slider_window = pset.get<size_t>("SliderWindow");
  }
  
  int ClassicDqDsSmoother::SmoothDqDs(SimpleCluster& cluster) const
  {

    //int                    & _start_index      = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    std::vector<double>    & _dqds_v           = cluster._dqds_v;
    std::vector<double>    & _dqds_slider      = cluster._dqds_slider;
    //bool                   & _hits_ordered     = cluster._hits_ordered;
    //std::vector<double>    & _ds_v        = cluster._ds_v;

     if (_dqds_v.size() != _s_hit_v.size()) {
      CT_ERROR() << " dqds size is different than hit vector size" << std::endl;
      throw HitCosmicTagException();
    }

    if (_dqds_v.size() < _slider_window) {
      CT_NORMAL() <<  "Not enough hits" << std::endl;
      return 0;
    }

    //size_t window = _slider_window;

    if (_slider_window % 2 != 0) {
      CT_ERROR() << "_slider_window has to be even." << std::endl; 
      throw HitCosmicTagException();
    }

    _dqds_slider.clear();
    _dqds_slider.reserve(_dqds_v.size());

    std::vector<double> dqds_vector = _dqds_v;

    for(const auto& window : get_windows(dqds_vector, _slider_window)) {

      double median_dqds = get_smooth_trunc_median(window);
      _dqds_slider.push_back(median_dqds);
      //if (_debug) std::cout << "dqds_slider value " << median_dqds << std::endl;

    }

    return _s_hit_v.size();
  }

}
#endif
