#ifndef CLASSICLocalLinearityCalculator_CXX
#define CLASSICLocalLinearityCalculator_CXX

#include "ClassicLocalLinearityCalculator.h"


namespace cosmictag {

  static ClassicLocalLinearityCalculatorFactory __global_ClassicLocalLinearityCalculatorFactory__;

  ClassicLocalLinearityCalculator::ClassicLocalLinearityCalculator(const std::string name)
    : BaseLocalLinearityCalculatorAlgo(name)
  {}

  void ClassicLocalLinearityCalculator::_Configure_(const Config_t &pset)
  {
    _max_allowed_hit_distance = pset.get<double>("MaxAllowedHitDistance");
    _slider_window = pset.get<size_t>("SliderWindow", 20); 
  }
  
  int ClassicLocalLinearityCalculator::CalculateLocalLinearity(SimpleCluster& cluster) const
  {

    //int                    & _start_index      = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    std::vector<double>    & _linearity_v      = cluster._linearity_v;
    bool                   & _linearity_is_set = cluster._linearity_is_set;


    if (_s_hit_v.size() < 2) {
      CT_DEBUG() << "Cannot calculate linearity if less than 2 hits" << std::endl;
      return 0;
    }

    std::vector<double> R;
    R.reserve(_s_hit_v.size());

    std::vector<double> X;
    std::vector<double> Y;
    X.reserve(_slider_window);
    Y.reserve(_slider_window);

    for(const auto& window : get_windows(_s_hit_v, _slider_window)) {

      for(const auto& s_hit : window) {
        X.push_back(s_hit.wire); 
        Y.push_back(s_hit.time);
      }

      auto c  = cov(X,Y);
      auto sX = stdev(X);
      auto sY = stdev(Y);
      auto r  = std::abs(c/(sX * sY));

      //if(_debug) {
        //std::cout << "c: "  << c << std::endl
        //          << "sX: " << sX <<  std::endl
        //          << "sY: " << sY <<  std::endl
        //          << "r: "  << r <<  std::endl;
        //std::cout << "Local Linearity: " << r << std::endl;
      //}

      if(std::isnan(r)) r = 0.0; 
      R.push_back(r);
      
      X.clear(); Y.clear();
    }   
 
    //first and last points will be nan. Lets set them equal to the points just above and below
    R.at(0)            = R.at(1);
    R.at(R.size() - 1) = R.at(R.size() - 2);
    
    _linearity_v = R;

    _linearity_is_set = true;

    return _s_hit_v.size();
  }
}
#endif
