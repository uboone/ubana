#ifndef STOPMUBRAGG_CXX
#define STOPMUBRAGG_CXX

#include "StopMuBragg.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"

namespace cosmictag {

  static StopMuBraggFactory __global_StopMuBraggFactory__;

  StopMuBragg::StopMuBragg(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
  {}

  void StopMuBragg::_Configure_(const Config_t &pset)
  {
    _hits_to_remove              = pset.get< int > ( "HitsToRemove", 3 );
    _pre_post_window             = pset.get< int > ( "PrePostWindow", 5 ); 
    _max_muon_hits               = pset.get< int > ( "MaxMuonHits", 2000 );
    _min_muon_hits               = pset.get< int > ( "MinMuonHits", 20 );
    _local_linearity_threshold   = pset.get< double > ( "LocalLinerityThreshold", 0.85 );
    _perc_diff_cut               = pset.get< double > ( "PercDiffCut", 20 );
  }

  void StopMuBragg::PrintConfig() const 
  {
    std::cout << "--- StopMuBragg Config" << std::endl;
    std::cout << "--- _hits_to_remove = " << _hits_to_remove << std::endl;
    std::cout << "--- _pre_post_window = " << _pre_post_window << std::endl;
    std::cout << "--- _max_muon_hits = " << _max_muon_hits << std::endl;
    std::cout << "--- _min_muon_hits = " << _min_muon_hits << std::endl;
    std::cout << "--- _local_linearity_threshold = " << _local_linearity_threshold << std::endl;
    std::cout << "--- _perc_diff_cut = " << _perc_diff_cut << std::endl;
  }

 
  bool StopMuBragg::IsStopMuBragg(const cosmictag::SimpleCluster & cluster) const {


    const std::vector<double>    & _dqds_slider      = cluster._dqds_slider;
    const std::vector<double>    & _linearity_v      = cluster._linearity_v;

    if (_dqds_slider.size() < (unsigned int) (_pre_post_window * 2)) {
      CT_DEBUG() << "Can't make decision, number of simple hits is " << _dqds_slider.size() 
                 << ", which is less then " << _pre_post_window * 2 << std::endl;
      return false;
    }

    // Find the hits with the maximum dqds, that one will be the hit
    // where the Bragg peak is
    size_t bragg_index;
    auto it_max = std::max_element(_dqds_slider.begin(), _dqds_slider.end());
    bragg_index = it_max - _dqds_slider.begin();
    double bragg_dqds = *it_max;
    for (bool flag = true; flag && it_max != _dqds_slider.end(); it_max++) {
      if (*it_max < bragg_dqds) {
        bragg_index = --it_max - _dqds_slider.begin();
        flag = false;
      }
    }


    CT_DEBUG() << "Bragg peak hit index is " << bragg_index << std::endl;

    // Check that the number of muon hits are below the maximum allowed
    int n_muon_hits = bragg_index + 1;
    if (n_muon_hits > _max_muon_hits) {
      CT_DEBUG() << "Number of muon hits is " << n_muon_hits
                 << " which is above maximum allowed (" << _max_muon_hits << ")" << std::endl;
      return false;
    }

    // Check that the number of muon hits are above the minimum allowed
    if (n_muon_hits < _min_muon_hits) {
      CT_DEBUG() << "Number of muon hits is " << n_muon_hits
                 << " which is below minimum allowed (" << _min_muon_hits << ")" << std::endl;
      return false;
    }

    // In this case we are looking for events that don't have a Michel, 
    // so we want to ensure that the local linearity is not below threshold
    // in the Bragg region
    /*
    double bragg_local_linearity = _linearity_v.at(bragg_index);
    if (bragg_local_linearity < _local_linearity_threshold) {
      if (_debug) std::cout << "[IsStopMuBragg] Local linearity is " << bragg_local_linearity
                            << " which is less than threshold (" << _local_linearity_threshold << ")" << std::endl;
      return false;
    }
    */

    // We actually want that there is no kink in this cluster,
    // as we just want the muon to stop. But exclude first hits 
    // as things can get funny at the beginning
    for (size_t l = _hits_to_remove; l < _linearity_v.size() - _hits_to_remove; l++) {
      if (_linearity_v.at(l) < _local_linearity_threshold) {
        CT_DEBUG() << "Local linearity at hit " << l << " is " << _linearity_v.at(l)
                   << " which is less than threshold (" << _local_linearity_threshold << ")" << std::endl;
        return false;
      }
    }


    // Take firsts hits, then lasts hits
    bool _use_mean = false; 
    double start_mean, end_mean;

    if (_use_mean) { 

      start_mean = std::accumulate(_dqds_slider.begin(), _dqds_slider.begin() + _pre_post_window, 0);
      start_mean /= _pre_post_window;
      end_mean = std::accumulate(_dqds_slider.end() - _pre_post_window, _dqds_slider.end(), 0);
      end_mean /= _pre_post_window;

    } else {

      std::vector<double>::const_iterator first = _dqds_slider.begin();
      std::vector<double>::const_iterator last  = _dqds_slider.begin() + _pre_post_window;
      std::vector<double> temp(first, last);
      std::sort(temp.begin(), temp.end());

      double _n_hits_remove = 7.;

      start_mean = std::accumulate(temp.begin(), temp.begin() + _n_hits_remove, 0);
      start_mean /= _n_hits_remove;

      first = _dqds_slider.end() - _pre_post_window;
      last  = _dqds_slider.end();
      std::vector<double> temp2(first, last);
      std::sort(temp2.begin(), temp2.end());

      end_mean = std::accumulate(temp2.begin(), temp2.begin() + _n_hits_remove, 0);
      end_mean /= _n_hits_remove;
    }

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    CT_DEBUG() << "Start mean: " << start_mean 
               << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff < -_perc_diff_cut) {
      return true;
    }

    return false;
  }
}


#endif
