#ifndef HITCOSMICTAG_TOOLS_CXX
#define HITCOSMICTAG_TOOLS_CXX

#include "Tools.h"

#include <vector>
#include <numeric>
#include <string>
#include <cmath>

namespace cosmictag {


    double get_smooth_trunc_median(std::vector<double> v) {

    if (v.size() > 2) {
      // Find and erase max element
      auto it_max = std::max_element(v.begin(), v.end());
      v.erase(it_max);

      // Find and erase min element
      auto it_min = std::min_element(v.begin(), v.end());
      v.erase(it_min);
    }

    double median = -1;

    size_t size = v.size();
    std::sort(v.begin(), v.end());
    if (size % 2 == 0){
      median = (v[size/2 - 1] + v[size/2]) / 2;
    }
    else{
      median = v[size/2];
    }

    return median;
  }



  double mean(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to mean" << std::endl;
        
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }





  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2)
  {
    if(data1.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;
    if(data2.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }




  double stdev(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to stdev" << std::endl;

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return std::sqrt(result/((double)data.size()));
  }




  

}
#endif
