/**
 * \file StopMuBragg.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class StopMuBragg
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/
#ifndef STOPMUBRAGG_H
#define STOPMUBRAGG_H

#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>
#include "../Base/BaseAlgorithm.h"
#include "../Base/CustomAlgoFactory.h"

namespace cosmictag{
/**
   \class StopMuBragg
   User defined class StopMuBragg ... these comments are used to generate
   doxygen documentation!
 */

  class StopMuBragg : public cosmictag::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    StopMuBragg(const std::string name="StopMuBragg");
    
    /// Default destructor
    ~StopMuBragg(){}

    bool IsStopMuBragg(const cosmictag::SimpleCluster & cluster) const;

    void PrintConfig() const;


  protected:

    void _Configure_(const Config_t &pset);
    
    int _hits_to_remove;
    int _pre_post_window;
    int _max_muon_hits;
    int _min_muon_hits;
    double _local_linearity_threshold;
    double _perc_diff_cut;

  };
  
  /**
     \class cosmictag::StopMuBraggFactory
  */
  class StopMuBraggFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    StopMuBraggFactory() { CustomAlgoFactory::get().add_factory("StopMuBragg",this); }
    /// dtor
    ~StopMuBraggFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new StopMuBragg(instance_name); }
  };
} 

#endif
/** @} */ // end of doxygen group 

