/**
 * \file StopMuMichel.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class StopMuMichel
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/
#ifndef STOPMUMICHEL_H
#define STOPMUMICHEL_H

#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>
#include "../Base/BaseAlgorithm.h"
#include "../Base/CustomAlgoFactory.h"

namespace cosmictag{
/**
   \class StopMuMichel
   User defined class StopMuMichel ... these comments are used to generate
   doxygen documentation!
 */

  class StopMuMichel : public cosmictag::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    StopMuMichel(const std::string name="StopMuMichel");
    
    /// Default destructor
    ~StopMuMichel(){}

    bool IsStopMuMichel(const cosmictag::SimpleCluster & cluster);

    void PrintConfig();


  protected:

    void _Configure_(const Config_t &pset);
    
    int _hits_to_remove;
    int _pre_post_window;
    int _max_muon_hits;
    int _min_muon_hits;
    int _max_michel_hits;
    int _min_michel_hits;
    int _max_end_hits;
    double _local_linearity_threshold;
    double _perc_diff_cut;

    std::vector<double> dqds_end;

  };
  
  /**
     \class cosmictag::StopMuMichelFactory
  */
  class StopMuMichelFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    StopMuMichelFactory() { CustomAlgoFactory::get().add_factory("StopMuMichel",this); }
    /// dtor
    ~StopMuMichelFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new StopMuMichel(instance_name); }
  };
} 

#endif
/** @} */ // end of doxygen group 

