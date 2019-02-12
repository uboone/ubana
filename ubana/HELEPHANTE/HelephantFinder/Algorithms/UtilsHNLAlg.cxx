/******************************************************************************
 * @file UtilsHNLAlg.cxx
 * @brief Useful functions for other HNL algorithms
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see UtilsHNLAlg.h
 * ****************************************************************************/
#include "UtilsHNLAlg.h"

namespace UtilsHNL
{
  // Constructor/destructor
  UtilsHNLAlg::UtilsHNLAlg(){}
  UtilsHNLAlg::~UtilsHNLAlg(){}


  //_________________________________________________________________________________
  // Get mean of a vector
  float UtilsHNLAlg::GetMean(std::vector<double> const & vec)
  {
    float mean = -9999;
    size_t size = vec.size();
    if (size == 0) return mean;

    float sum = 0;
    for (auto v : vec)
    {
      sum += v;
    }
    mean = sum/(float)size;
    return mean;
  } // END function GetMean


  //_________________________________________________________________________________
  // Get median of a vector
  float UtilsHNLAlg::GetMedian(std::vector<double> const & vec)
  {
    float median = -9999;
    size_t size = vec.size();
    if (size == 0) return median;

    std::vector<double> sorted_vec = vec;
    std::sort(sorted_vec.begin(), sorted_vec.end());
    if (size % 2 == 0) median = (sorted_vec[size/2 - 1] + sorted_vec[size/2]) / 2;
    else median = sorted_vec[size/2];

    return median;
  } // END function GetMedian


  //_________________________________________________________________________________
  // Get standard deviation of a vector
  float UtilsHNLAlg::GetVariance(std::vector<double> const & vec)
  {
    float variance = -1;
    float sum = 0;
    float sum2 = 0;
    size_t size = vec.size();
    if (size == 0) return variance;

    for (auto value : vec)
    {
      sum  += value;
      sum2 += value*value;
    }  

    variance = sum2/(float)size - (sum/(float)size)*(sum/(float)size);

    return variance;
  } // END function GetVariance


  //_________________________________________________________________________________
  // Get standard deviation of a vector
  float UtilsHNLAlg::GetSTD(std::vector<double> const & vec)
  {
    if (vec.size() == 0) return -9999;
    float variance = GetVariance(vec);
    
    if (variance > 0) return std::sqrt(variance);
    else return -9999;
  } // END function GetSTD


  //_________________________________________________________________________________
  // Get truncated mean of a vector
  float UtilsHNLAlg::GetTruncatedMean(std::vector<double> const & vec)
  {
    float result = -9999;
    float n = 1.;
    if (vec.size() == 0) return result;

    float median = GetMedian(vec);
    float std = GetSTD(vec);

    std::vector<double> vec_trimmed;
    vec_trimmed.clear();

    for (auto v : vec)
    {
      if (v > median - n * std && v < median + n * std)
      {
        vec_trimmed.emplace_back(v);
      }
    }
    result = GetMean(vec_trimmed);
    return result;
  } // END function GetTruncatedMean




} // END namespace dEdX
