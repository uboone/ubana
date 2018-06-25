/**
 * \file BaseAlgorithm.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseAlgorithm
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEALGORITHM_H
#define HITCOSMICTAG_BASEALGORITHM_H

#include "DataTypes.h"
#include "HitCosmicTagFMWKInterface.h"
#include "LoggerFeature.h"

#include <vector>

namespace cosmictag {

  class CosmicTagManager;
  
  /**
     \class BaseAlgorithm
  */
  class BaseAlgorithm : public LoggerFeature {

    friend class CosmicTagManager;
    
  public:
    
    /// Default constructor
    BaseAlgorithm(const AlgoType type, const std::string name);
    
    /// Default destructor
    ~BaseAlgorithm(){}

    /// Function to accept configuration
    void Configure(const Config_t &pset);

    /// Algorithm type
    AlgoType AlgorithmType() const;

    /// Algorithm name
    const std::string& AlgorithmName() const;
    
  protected:

    virtual void _Configure_(const Config_t &pset) = 0;

  private:
    
    AlgoType _type; ///< Algorithm type
    std::string _name; ///< Algorithm name
   
  };
}
#endif
/** @} */ // end of doxygen group 

