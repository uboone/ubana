#ifndef CLASSICSDQDSCALCULATOR_CXX
#define CLASSICSDQDSCALCULATOR_CXX

#include "ClassicDqDsCalculator.h"


namespace cosmictag {

  static ClassicDqDsCalculatorFactory __global_ClassicDqDsCalculatorFactory__;

  ClassicDqDsCalculator::ClassicDqDsCalculator(const std::string name)
    : BaseDqDsCalculatorAlgo(name)
  {}

  void ClassicDqDsCalculator::_Configure_(const Config_t &pset)
  {
    _w2cm = pset.get<double>("WireToCmConstant");
    _t2cm = pset.get<double>("TimeToCmConstant");
    _dqds_calib = pset.get<double>("GainCalib");
  }
  
  int ClassicDqDsCalculator::CalculateDqDs(SimpleCluster& cluster) const
  {

    //int                    & _start_index      = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    std::vector<double>    & _dqds_v           = cluster._dqds_v;
    std::vector<double>    & _ds_v             = cluster._ds_v;
    bool                   & _hits_ordered     = cluster._hits_ordered;
    

     if (!_hits_ordered) {
      CT_ERROR() << "Need to order hits first." << std::endl;
      throw HitCosmicTagException();
    }
 
    if (_ds_v.size() != _s_hit_v.size()) {
      CT_ERROR() << ": ds size is different than hit vector size" << std::endl;
      throw HitCosmicTagException();
    }
 
    _dqds_v.clear();
    _dqds_v.reserve(_s_hit_v.size());

    double ds = 1.;

    for (size_t i = 0; i < _s_hit_v.size()-1; i++) {

      TVector3 this_point(_s_hit_v.at(i).wire * _w2cm, _s_hit_v.at(i).time*4*_t2cm, 0);
      TVector3 next_point(_s_hit_v.at(i+1).wire * _w2cm, _s_hit_v.at(i+1).time*4*_t2cm, 0);
      ds = (this_point - next_point).Mag();

      _dqds_v.emplace_back(_s_hit_v.at(i).integral/ds * _dqds_calib);
      CT_DEBUG() << "_dqds_v.back() " << _dqds_v.back() << std::endl;

    }

    // Finish with the last point
    _dqds_v.emplace_back(_s_hit_v.at(_s_hit_v.size()-1).integral/ds * _dqds_calib);

    return _s_hit_v.size();
  }
}
#endif
