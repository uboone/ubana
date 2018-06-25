/**
 * \file BaseHitOrdererAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseHitOrdererAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEHITORDERERALGO_H
#define HITCOSMICTAG_BASEHITORDERERALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseHitOrdererAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseHitOrdererAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseHitOrdererAlgo(const std::string name="noname") 
       : BaseAlgorithm(kHitOrderer,name)
    {}
 
    /// Default destructor
    virtual ~BaseHitOrdererAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int OrderHits(SimpleCluster&) const = 0;

    virtual void CollectionCoplanar(bool status) = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

