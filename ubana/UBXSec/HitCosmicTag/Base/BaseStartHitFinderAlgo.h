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
#ifndef HITCOSMICTAG_BASESTARTHITFINDERALGO_H
#define HITCOSMICTAG_BASESTARTHITFINDERALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseHitOrdererAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseStartHitFinderAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseStartHitFinderAlgo(const std::string name="noname") 
       : BaseAlgorithm(kStartHitFinder,name)
    {}
 
    /// Default destructor
    virtual ~BaseStartHitFinderAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int FindStartHit(SimpleCluster&, SimpleHit&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

