#ifndef QHI2MATCH_CXX
#define QHI2MATCH_CXX

#include "Chi2Match.h"
#include "ubana/LLSelectionTool/OpT0Finder/Base/OpT0FinderException.h"
#include <TMinuit.h>
#include <cmath>
#include <numeric>
#include <TMath.h>

namespace flashana {

  static Chi2MatchFactory __global_Chi2MatchFactory__;

  void MIN_vtx_qll(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
  
  Chi2Match::Chi2Match(const std::string name)
    : BaseFlashMatch(name), _normalize(false)
  { _current_llhd = _current_chi2 = -1.0; }

  Chi2Match::Chi2Match()
  { throw OpT0FinderException("Use the other ctor"); }
  
  void Chi2Match::_Configure_(const Config_t &pset) {
    _normalize = pset.get<bool>("NormalizeHypothesis");
    _penalty_threshold_v = pset.get<std::vector<double> >("PEPenaltyThreshold");
    _penalty_value_v = pset.get<std::vector<double> >("PEPenaltyValue");

    _recox_penalty_threshold = pset.get<double>("XPenaltyThreshold");
    _recoz_penalty_threshold = pset.get<double>("ZPenaltyThreshold");

    _onepmt_score_threshold = pset.get<double>("OnePMTScoreThreshold");
    _onepmt_xdiff_threshold = pset.get<double>("OnePMTXDiffThreshold");
    _onepmt_pesum_threshold = pset.get<double>("OnePMTPESumThreshold");
    _onepmt_pefrac_threshold = pset.get<double>("OnePMTPEFracThreshold");
  }
  
  FlashMatch_t Chi2Match::Match(const QCluster_t &pt_v, const Flash_t &flash) {
    
    //
    // Prepare TPC
    //
    _raw_trk.resize(pt_v.size());
    double min_x = 1e9;
    double max_x = 0;
    for (size_t i = 0; i < pt_v.size(); ++i) {
      auto const &pt = pt_v[i];
      _raw_trk[i] = pt;
      if (pt.x < min_x) { min_x = pt.x; _raw_xmin_pt = pt; }
      if (pt.x > max_x) { max_x = pt.x; _raw_xmax_pt = pt; }
    }
    //for (auto &pt : _raw_trk) pt.x -= min_x;

    //
    // Prepare Flash
    //
    if (_hypothesis.pe_v.empty()) _hypothesis.pe_v.resize(NOpDets(), 0.);
    if (_hypothesis.pe_v.size() != NOpDets()) {
      throw OpT0FinderException("Hypothesis vector length != PMT count");
    }
    
    for (auto &v : _hypothesis.pe_v) v = 0;
    

    FillEstimate(_raw_trk, _hypothesis);

    if (_normalize) {
      double qsum = std::accumulate(std::begin(_hypothesis.pe_v),
				    std::end(_hypothesis.pe_v),
				    0.0);
      for (auto &v : _hypothesis.pe_v) v /= qsum;
    }

    FlashMatch_t res;

    res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;
    
    double weight = 0;
    
    for (size_t pmt_index = 0; pmt_index < NOpDets(); ++pmt_index) {

      res.tpc_point.y += OpDetY(pmt_index) * _hypothesis.pe_v[pmt_index];
      res.tpc_point.z += OpDetZ(pmt_index) * _hypothesis.pe_v[pmt_index];
      
      weight += _hypothesis.pe_v[pmt_index];
    }
    
    res.tpc_point.y /= weight;
    res.tpc_point.z /= weight;

    res.tpc_point.x = _reco_x_offset;
    res.tpc_point_err.x = _reco_x_offset_err;
    
    res.hypothesis  = _hypothesis.pe_v;

    //
    // Compute score
    //
    res.score = Chi2(_hypothesis, flash);
    
    return res;

    /*
    auto res1 = PESpectrumMatch(pt_v,flash,true);
    auto res2 = PESpectrumMatch(pt_v,flash,false);

    std::cout << "Using   mid-x-init: " << res1.tpc_point.x << " [cm] @ " << res1.score << std::endl;
    std::cout << "Without mid-x-init: " << res2.tpc_point.x << " [cm] @ " << res2.score << std::endl;


    auto res = (res1.score > res2.score ? res1 : res2);

    if(res.score < _onepmt_score_threshold) {

      auto res_onepmt = OnePMTMatch(flash);

      if(res_onepmt.score >= 0.)
	return res_onepmt;
    }

    return res;
    */
  }
  
  
  double Chi2Match::Chi2(const Flash_t &hypothesis,
			 const Flash_t &measurement) {
    
    double nvalid_pmt = 0;
    
    double PEtot_Hyp = 0;
    for (auto const &pe : hypothesis.pe_v)
      PEtot_Hyp += pe;
    double PEtot_Obs = 0;
    for (auto const &pe : measurement.pe_v)
      PEtot_Obs += pe;
    
    _current_chi2 = 0.;
    
    if (measurement.pe_v.size() != hypothesis.pe_v.size())
      throw OpT0FinderException("Cannot compute QLL for unmatched length!");
    
    double O, H, Error;
    for (size_t pmt_index = 0; pmt_index < hypothesis.pe_v.size(); ++pmt_index) {
      
      O = measurement.pe_v[pmt_index]; // observation
      H = hypothesis.pe_v[pmt_index];  // hypothesis

      //std::cout << "@ PMT index " << pmt_index << " O : " << O << " H : " << H;
      
      if( H < 0 ) throw OpT0FinderException("Cannot have hypothesis value < 0!");
      
      if(O < 0) {
	if(H < _penalty_threshold_v[pmt_index]) continue;
	O = _penalty_value_v[pmt_index];
      }
      
      nvalid_pmt += 1;

      Error = O;
      if( Error < 1.0 ) Error = 1.0;

      auto chilocal = std::pow((O - H), 2) / (Error);

      //std::cout << "\t chi local : " << chilocal << std::endl;

      _current_chi2 += chilocal;
      
    }// for all PMTs
    
    _current_chi2 /= nvalid_pmt;

    return _current_chi2;
  }
  
  
}
#endif
