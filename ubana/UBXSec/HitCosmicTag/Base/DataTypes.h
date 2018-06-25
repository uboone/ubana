#ifndef HITCOSMICTAG_DATATYPES_H
#define HITCOSMICTAG_DATATYPES_H

#include <vector>
#include <numeric>
#include <string>

namespace cosmictag {

  /// Enumerator for different types of algorithm
  enum AlgoType{
    kAlgoUnknown = -1,  ///< Algorithm type to f?
    kStartHitFinder = 0,
    kHitOrderer,
    kHitSmoother,
    kDqDsCalculator,
    kDqDsSmoother,
    kLocalLinearity,
    kCustomAlgo,
    kAlgorithmTypeMax   ///< enum flag for algorithm type count & invalid type
  };

  /// Struct to represent a hit
   struct SimpleHit {

    public: 

    double time;     ///< in cm, basically x
    int wire;        ///< in cm, given the wire pitch
    int plane;       ///< plane number
    double integral; ///< integral (~charge)

    double t;        ///< in time ticks
    int w;           ///< in wire number

    SimpleHit () {
      time = t = -1;
      wire = w = -1;
      plane = -1;
    }

    bool operator==(const SimpleHit& x) const {
      return (time == x.time) && 
             (wire == x.wire) &&
             (plane == x.plane);
    }
  };

  /// Struct to represent a cluster of SimpleHits
   struct SimpleCluster {

    public: 

    std::vector<SimpleHit> _s_hit_v;
    int _start_index = -1;
    std::vector<double> _dqds_v;
    std::vector<double> _ds_v;
    std::vector<double> _dqds_slider;
    std::vector<double> _linearity_v;

    bool _linearity_is_set = false;
    bool _hits_ordered = false;
    bool _start_hit_is_set = false;

    SimpleCluster () {
      _s_hit_v.clear();
      _dqds_v.clear();
      _ds_v.clear();
      _dqds_slider.clear();
      _linearity_v.clear();

      _linearity_is_set = false;
      _hits_ordered = false;
      _start_hit_is_set = false;

    }

    SimpleCluster (std::vector<SimpleHit> v) {
      _s_hit_v.clear();
      _s_hit_v = v;
      _dqds_v.clear();
      _ds_v.clear();
      _dqds_slider.clear();
      _linearity_v.clear();

      _linearity_is_set = false;
      _hits_ordered = false;
      _start_hit_is_set = false;

    }
  };

 
  

  
  
  namespace msg {
    /// Verbosity message level
    enum Level_t {
      kDEBUG,
      kINFO,
      kNORMAL,
      kWARNING,
      kERROR,
      kEXCEPTION,
      kMSG_TYPE_MAX
    };
    
    const std::string kStringPrefix[kMSG_TYPE_MAX] =
      {
	" [DEBUG] ", ///< DEBUG message prefix
	" [INFO] ", ///< INFO message prefix
	" [NORMAL] ", ///< NORMAL message prefix
	" [WARNING] ", ///< WARNING message prefix
	" [ERROR] ", ///< ERROR message prefix
	" [EXCEPTION] "  ///< CRITICAL message prefix
      };
    ///< Prefix of message
  }
}
#endif
