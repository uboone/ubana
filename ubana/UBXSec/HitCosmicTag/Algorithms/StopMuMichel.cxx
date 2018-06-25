#ifndef STOPMUMICHEL_CXX
#define STOPMUMICHEL_CXX

#include "StopMuMichel.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"

namespace cosmictag {

  static StopMuMichelFactory __global_StopMuMichelFactory__;

  StopMuMichel::StopMuMichel(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
  {}

  void StopMuMichel::_Configure_(const Config_t &pset)
  {
    _hits_to_remove              = pset.get< int > ( "HitsToRemove", 3 );
    _pre_post_window             = pset.get< int > ( "PrePostWindow", 5 ); 
    _max_muon_hits               = pset.get< int > ( "MaxMuonHits", 2000 );
    _min_muon_hits               = pset.get< int > ( "MinMuonHits", 20 );
    _max_michel_hits             = pset.get< int > ( "MaxMichelHits", 70 );
    _min_michel_hits             = pset.get< int > ( "MinMichelHits", 2 );
    _max_end_hits                = pset.get< int > ( "MaxEndHits", 120 ); 
    _local_linearity_threshold   = pset.get< double > ( "LocalLinerityThreshold", 0.85 );
    _perc_diff_cut               = pset.get< double > ( "PercDiffCut", 20 );
  }

  void StopMuMichel::PrintConfig()
  {
    std::cout << "--- StopMuMichel Config" << std::endl;
    std::cout << "--- _hits_to_remove = " << _hits_to_remove << std::endl;
    std::cout << "--- _pre_post_window = " << _pre_post_window << std::endl;
    std::cout << "--- _max_muon_hits = " << _max_muon_hits << std::endl;
    std::cout << "--- _min_muon_hits = " << _min_muon_hits << std::endl;
    std::cout << "--- _max_michel_hits = " << _max_michel_hits << std::endl;
    std::cout << "--- _min_michel_hits = " << _min_michel_hits << std::endl;
    std::cout << "--- _max_end_hits = " << _max_end_hits << std::endl;
    std::cout << "--- _local_linearity_threshold = " << _local_linearity_threshold << std::endl;
    std::cout << "--- _perc_diff_cut = " << _perc_diff_cut << std::endl;

  }

 
  bool StopMuMichel::IsStopMuMichel(const cosmictag::SimpleCluster & cluster) {

    const std::vector<double>    & _dqds_slider      = cluster._dqds_slider;
    const std::vector<double>    & _linearity_v      = cluster._linearity_v;


    if (_dqds_slider.size() < (unsigned int) (_hits_to_remove * 2 + _pre_post_window * 2 + 6)) {
      CT_DEBUG() << "Can't make decision, number of simple hits is " << _dqds_slider.size() 
                 << ", which is less then " << _hits_to_remove * 2 + _pre_post_window * 2 + 6<< std::endl;
      return false;
    }


    // Vertex must not be in the FV
    //if (_fv.InFV(_vertex))
    //  return false;


    // Find the hits with the maximum dqds, that one will be the hit
    // where the Bragg peak is. If this is a big cluster, is likely that 
    // there'll be a big delta ray that may fake a Bragg peak. In this case, 
    // look only at the end of the cluster to find the Bragg peak.
    int offset = 1; // at least start from 1 as first hit is usually bad
    if (_dqds_slider.size() > (size_t)_max_end_hits) {
      CT_DEBUG() << "Many hits in this cluster ("
                 << _dqds_slider.size() << "), finding Bragg in the last "
                 << _max_end_hits << " hits." << std::endl;
      offset = _dqds_slider.size() - _max_end_hits;
    }
    int bragg_index;
    auto it_max = std::max_element(_dqds_slider.begin()+offset, _dqds_slider.end());
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

    // Check that the local linearity is less than threshold in the Bragg region
    double bragg_local_linearity = _linearity_v.at(bragg_index);

    if (bragg_local_linearity > _local_linearity_threshold) {
      CT_DEBUG() << "Local linearity is " << bragg_local_linearity 
                 << " which is above threshold (" << _local_linearity_threshold << ")" << std::endl;
      return false;
    }

    // Check that the local linearity in the muon region (before Bragg) 
    // is above threshold. Exclude first hits as things can get funny
    // at the beginning. Also exclude last hits as we expect the muon 
    // to curve in the last region
    /*
    for (size_t i = offset; i < (size_t) bragg_index - _pre_post_window; i++) {
      if (_linearity_v.at(i) < _local_linearity_threshold) {
        if (_debug) std::cout << "[IsStopMuMichel] Local linearity at hit " << i << " (before Bragg) is " << _linearity_v.at(i)
                              << " which is below threshold (" << _local_linearity_threshold << ")" << std::endl;
        return false;
      }
    }
    */

    // Check that the photon hits are below the maximum allowed
    int n_michel_hits = _dqds_slider.size() - bragg_index - 1;
    if (n_michel_hits > _max_michel_hits) {
      CT_DEBUG() << "Number of Michel hits is " << n_michel_hits
                 << " which is above maximum allowed (" << _max_michel_hits << ")" << std::endl;
      return false;
    }

    // Check that the photon hits are above the minimum allowed
    if (n_michel_hits < _min_michel_hits) {
      CT_DEBUG() << "Number of Michel hits is " << n_michel_hits
                 << " which is below the minimum allowed (" << _min_michel_hits << ")" << std::endl;
      return false;
    }

    // Get mean of first and last hits
    //dqds_end.clear();
    dqds_end = _dqds_slider;

    dqds_end.erase(dqds_end.begin(), dqds_end.begin() + bragg_index - (_pre_post_window + 5));
    dqds_end.erase(dqds_end.end() - _hits_to_remove, dqds_end.end());

    bragg_index = (_pre_post_window + 5);

    if (dqds_end.size() <= (size_t)bragg_index) {
      CT_DEBUG() << "Not enough hits." << std::endl;
      return false;
    }

    for (size_t i = 0; i < dqds_end.size(); i++) CT_DEBUG() << i << ": dqds_end = " << dqds_end.at(i) << std::endl;

    double start_mean = std::accumulate(dqds_end.begin(), dqds_end.begin() + _pre_post_window, 0);
    start_mean /= _pre_post_window;

    double end_mean = std::accumulate(dqds_end.end() - _pre_post_window, dqds_end.end(), 0);
    end_mean /= _pre_post_window;

    int edge = bragg_index + 5;

    if (dqds_end.size() - edge < (size_t) _pre_post_window) {
      int vector_size = (int) dqds_end.size();
      CT_DEBUG() << "Few Michel hits, calculating average only on "
                 << vector_size - edge << " hits." << std::endl;
      end_mean = std::accumulate(dqds_end.end() - (vector_size - edge), dqds_end.end(), 0);
      end_mean /= (double) vector_size - edge;
    }

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    CT_DEBUG() << "Start mean: " << start_mean 
               << ", end mean " << end_mean 
               << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff > _perc_diff_cut) {
      return true;
    } 

    return false;   


  }
}


#endif
