/**
 * \file ClassicHitOrderer.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class ClassicHitOrderer
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/

#ifndef CLASSICSDQDSCALCULATOR_H
#define CLASSICSDQDSCALCULATOR_H

#include <iostream>
#include "../Base/BaseDqDsCalculatorAlgo.h"
#include "../Base/DqDsCalculatorFactory.h"
#include "../Base/HitCosmicTagException.h"

#include <TVector3.h>

namespace cosmictag {
  /**
     \class ClassicHitOrderer
     User custom analysis class
   */
  class ClassicDqDsCalculator : public BaseDqDsCalculatorAlgo {
  
  public:

    /// Default constructor
    ClassicDqDsCalculator(const std::string name="ClassicDqDsCalculator");

    /// Default destructor
    virtual ~ClassicDqDsCalculator(){}

    /// ?
    int CalculateDqDs(SimpleCluster&) const;

  protected:

    void _Configure_(const Config_t &pset);
    double _w2cm;                             ///< ?
    double _t2cm;                             ///< ?
    double _dqds_calib;                       ///< ?
  };
  
  /**
     \class cosmictag::ClassicHitOrdererFactory
  */
  class ClassicDqDsCalculatorFactory : public DqDsCalculatorFactoryBase {
  public:
    /// ctor
    ClassicDqDsCalculatorFactory() { DqDsCalculatorFactory::get().add_factory("ClassicDqDsCalculator",this); }
    /// dtor
    ~ClassicDqDsCalculatorFactory() {}
    /// creation method
    BaseDqDsCalculatorAlgo* create(const std::string instance_name) { return new ClassicDqDsCalculator(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group 
