/**
 * \file LocalLinearityCalculatorAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class LocalLinearityCalculatorAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASELOCALLINEARITYCALCULATORALGO_H
#define HITCOSMICTAG_BASELOCALLINEARITYCALCULATORALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class LocalLinearityCalculatorAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseLocalLinearityCalculatorAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseLocalLinearityCalculatorAlgo(const std::string name="noname") 
       : BaseAlgorithm(kLocalLinearity,name)
    {}
 
    /// Default destructor
    virtual ~BaseLocalLinearityCalculatorAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int CalculateLocalLinearity(SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

