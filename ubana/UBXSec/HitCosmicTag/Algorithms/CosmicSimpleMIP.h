/**
 * \file CosmicSimpleMIP.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class CosmicSimpleMIP
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/
#ifndef COSMICSIMPLEMIP_H
#define COSMICSIMPLEMIP_H

#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>
#include "../Base/BaseAlgorithm.h"
#include "../Base/CustomAlgoFactory.h"

namespace cosmictag{
/**
   \class CosmicSimpleMIP
   User defined class CosmicSimpleMIP ... these comments are used to generate
   doxygen documentation!
 */

  class CosmicSimpleMIP : public cosmictag::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    CosmicSimpleMIP(const std::string name="CosmicSimpleMIP");
    
    /// Default destructor
    ~CosmicSimpleMIP(){}

    bool IsCosmicSimpleMIP(const cosmictag::SimpleCluster & cluster) const;

    void PrintConfig() const;


  protected:

    void _Configure_(const Config_t &pset);
    
    double _local_linearity_threshold;

  };
  
  /**
     \class cosmictag::CosmicSimpleMIPFactory
  */
  class CosmicSimpleMIPFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    CosmicSimpleMIPFactory() { CustomAlgoFactory::get().add_factory("CosmicSimpleMIP",this); }
    /// dtor
    ~CosmicSimpleMIPFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new CosmicSimpleMIP(instance_name); }
  };
} 

#endif
/** @} */ // end of doxygen group 

