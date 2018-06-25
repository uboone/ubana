/**
 * \file BaseHitSmootherAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseHitSmootherAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEHitSmootherALGO_H
#define HITCOSMICTAG_BASEHitSmootherALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseHitSmootherAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseHitSmootherAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseHitSmootherAlgo(const std::string name="noname") 
       : BaseAlgorithm(kHitSmoother,name)
    {}
 
    /// Default destructor
    virtual ~BaseHitSmootherAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int Smooth(SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

