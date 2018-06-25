/**
 * \file ClassicLocalLinearityCalculator.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class ClassicLocalLinearityCalculator
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/

#ifndef ClassicLocalLinearityCalculator_H
#define ClassicLocalLinearityCalculator_H

#include <iostream>
#include "../Base/BaseLocalLinearityCalculatorAlgo.h"
#include "../Base/LocalLinearityCalculatorFactory.h"
#include "../Base/HitCosmicTagException.h"
#include "../Base/Tools.h"

#include <TVector3.h>

namespace cosmictag {
  /**
     \class ClassicLocalLinearityCalculator
     User custom analysis class
   */
  class ClassicLocalLinearityCalculator : public BaseLocalLinearityCalculatorAlgo {
  
  public:

    /// Default constructor
    ClassicLocalLinearityCalculator(const std::string name="ClassicLocalLinearityCalculator");

    /// Default destructor
    virtual ~ClassicLocalLinearityCalculator(){}

    /// ?
    int CalculateLocalLinearity(SimpleCluster&) const;

  protected:

    void _Configure_(const Config_t &pset);
    double _max_allowed_hit_distance;         ///< ?
    size_t _slider_window; 
  };
  
  /**
     \class cosmictag::ClassicLocalLinearityCalculatorFactory
  */
  class ClassicLocalLinearityCalculatorFactory : public LocalLinearityCalculatorFactoryBase {
  public:
    /// ctor
    ClassicLocalLinearityCalculatorFactory() { LocalLinearityCalculatorFactory::get().add_factory("ClassicLocalLinearityCalculator",this); }
    /// dtor
    ~ClassicLocalLinearityCalculatorFactory() {}
    /// creation method
    BaseLocalLinearityCalculatorAlgo* create(const std::string instance_name) { return new ClassicLocalLinearityCalculator(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group 
