#ifndef HITCOSMICTAG_TOOLS_H
#define HITCOSMICTAG_TOOLS_H

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>

namespace cosmictag {

  /// ?
  template<typename T>
  std::vector<std::vector<T>> get_windows(const std::vector<T>& the_thing,
                                          const size_t window_size);
  
}
#include "Tools.tcxx" // for template functions

namespace cosmictag {

  /// ?
  double get_smooth_trunc_median(std::vector<double> v);

  
  /// ?
  double mean(const std::vector<double>& data);
  

  /// ?
  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2);
  

  /// ?
  double stdev(const std::vector<double>& data);
  
}
#endif
