#ifndef COSMICSIMPLEMIP_CXX
#define COSMICSIMPLEMIP_CXX

#include "CosmicSimpleMIP.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"

namespace cosmictag {

  static CosmicSimpleMIPFactory __global_CosmicSimpleMIPFactory__;

  CosmicSimpleMIP::CosmicSimpleMIP(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
  {}

  void CosmicSimpleMIP::_Configure_(const Config_t &pset)
  {
    _local_linearity_threshold   = pset.get< double > ( "LocalLinerityThreshold", 0.9 );
  }

  void CosmicSimpleMIP::PrintConfig() const 
  {
    std::cout << "--- CosmicSimpleMIP Config" << std::endl;
    std::cout << "--- _local_linearity_threshold = " << _local_linearity_threshold << std::endl;
  }

 
  bool CosmicSimpleMIP::IsCosmicSimpleMIP(const cosmictag::SimpleCluster & cluster) const {


    //const int                    & _start_index      = cluster._start_index;
    //const std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    //const bool                   & _start_hit_is_set = cluster._start_hit_is_set;
    const std::vector<double>    & _dqds_v           = cluster._dqds_v;
    const std::vector<double>    & _dqds_slider      = cluster._dqds_slider;
    const std::vector<double>    & _linearity_v      = cluster._linearity_v;


       // This algo excludes the first and last hit

    // Check that the local linearity in never below threshold
    for (size_t i = 1; i < _linearity_v.size() - 1; i++) {

      if (_linearity_v.at(i) < _local_linearity_threshold) {

        CT_DEBUG() << "Local linearity at hit " << i 
                   << " is " << _linearity_v.at(i) 
                   << " which is below threshold (" 
                   << _local_linearity_threshold << ")." << std::endl;
        return false;

      }
    }

    // Check compatibility
    if (_dqds_v.size() != _dqds_slider.size()) {
      CT_DEBUG() << "_dqds_v vector size is " << _dqds_v.size() 
                 << " which is different to _dqds_slider vector size, which is "
                 << _dqds_slider.size() << "." << std::endl;
      throw std::exception();
    }

    // Check that we have enough hits
    if (_dqds_v.size() < 6 + 1 + 6 + 1) {
      CT_DEBUG() << "Not enough hits." << std::endl;
      return false;
    }

    // Now verify the first and last hits are flat in dqds
    std::vector<double> start_v (_dqds_v.begin() + 1, _dqds_v.begin() + 6);
    std::vector<double> end_v (_dqds_v.end() - 6, _dqds_v.end() - 1);

    //double std_start = stdev(start_v);
    //double std_end = stdev(end_v);
    //std::cout << "std_start: " << std_start << ", std_end: " << std_end << std::endl;

    bool good_start = true;
    bool good_end = true;

    for (auto q : start_v) {
      if (q < 40000 || q > 75000) {
        good_start = false;
      }
    }

    for (auto q : end_v) { 
      if (q < 40000 || q > 75000) {
        good_end = false;
      }
    }

    CT_DEBUG() << "Start is " << (good_start ? "GOOD" : "BAD") << std::endl;
    CT_DEBUG() << "End is " << (good_end ? "GOOD" : "BAD") << std::endl;


    // Now use the smoothed dqds to evaluate the start and 
    // end dqds average value
    double start_mean = 0, end_mean = 0;

    start_mean = std::accumulate(_dqds_slider.begin() + 1, 
                                 _dqds_slider.begin() + 6, 0);
    start_mean /= 5.;

    end_mean = std::accumulate(_dqds_slider.end() - 6,
                               _dqds_slider.end() - 1, 0);
    end_mean /= 5.;

    double perc_diff = (start_mean - end_mean) / start_mean * 100.; 

    CT_DEBUG() << "Start mean: " << start_mean
               << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl; 

    if (good_start && good_end && std::abs(perc_diff) < 10) {
      return true;
    }

    return false;
  }
}


#endif
