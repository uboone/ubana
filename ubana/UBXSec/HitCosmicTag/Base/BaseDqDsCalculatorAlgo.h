/**
 * \file BaseDqDsCalculatorAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseDqDsCalculatorAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEDQDSCALCULATORALGO_H
#define HITCOSMICTAG_BASEDQDSCALCULATORALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseDqDsCalculatorAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseDqDsCalculatorAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseDqDsCalculatorAlgo(const std::string name="noname") 
       : BaseAlgorithm(kDqDsCalculator,name)
    {}
 
    /// Default destructor
    virtual ~BaseDqDsCalculatorAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int CalculateDqDs(SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

