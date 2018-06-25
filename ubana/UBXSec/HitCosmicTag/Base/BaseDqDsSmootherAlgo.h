/**
 * \file BaseDqDsSmootherAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseDqDsSmootherAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEDQDSSMOOTHERALGO_H
#define HITCOSMICTAG_BASEDQDSSMOOTHERALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseDqDsSmootherAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseDqDsSmootherAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseDqDsSmootherAlgo(const std::string name="noname") 
       : BaseAlgorithm(kDqDsSmoother,name)
    {}
 
    /// Default destructor
    virtual ~BaseDqDsSmootherAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int SmoothDqDs(SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

