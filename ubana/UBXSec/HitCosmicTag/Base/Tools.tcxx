#ifndef HITCOSMICTAG_TOOLS_IMPL_H
#define HITCOSMICTAG_TOOLS_IMPL_H

//#include "Tools.h"

#include <vector>
#include <numeric>
#include <string>

namespace cosmictag {

  /// Enumerator for different types of algorithm
  template<typename T>
  std::vector<std::vector<T>> get_windows(const std::vector<T>& the_thing,
                                          const size_t window_size)
  {

    // given a vector of values return a vector of the same length
    // with each element being a vector of the values of the local neighbors
    // of the element at position i in the original vector
    // input  : [0,1,2,3,4,5,6,...,...,N-3,N-2,N-1] (input vector of size N)
    // output  (assuming a value of 'w' below == 3):
    // 0th element: [0]
    // 1st element: [0,1,2]
    // 2nd element: [0,1,2,3,4]
    // jth element: [j-w,j-w+1,..,j+w-2,j+w-1]
    
    std::vector<std::vector<T>> data;
    
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();
    
    data.reserve(num);
    
    for(size_t i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      // if we are at the beginning of the vector (and risk accessing -1 elements)
      if(i < w)
        {
          for(size_t j = 0; j < 2 * (i%w) - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are at the end of the vector (and risk going past it)
      else if (i > num - w + 1)
        {
          for(size_t j = num - 2*((num - i)%w)-1 ; j < num; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are in the middle of the waveform
      else
        {
          for(size_t j = i - w; j < i + w - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      data.emplace_back(inner);
    }

    return data;
  
  }
}

#endif
