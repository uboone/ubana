/**
 * \file ClassicDqDsSmoother.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class ClassicDqDsSmoother
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/

#ifndef CLASSICSDQDSSMOOTHER_H
#define CLASSICSDQDSSMOOTHER_H

#include <iostream>
#include "../Base/BaseDqDsSmootherAlgo.h"
#include "../Base/DqDsSmootherFactory.h"
#include "../Base/HitCosmicTagException.h"
#include "../Base/Tools.h"

#include <TVector3.h>

namespace cosmictag {
  /**
     \class ClassicDqDsSmoother
     User custom analysis class
   */
  class ClassicDqDsSmoother : public BaseDqDsSmootherAlgo {
  
  public:

    /// Default constructor
    ClassicDqDsSmoother(const std::string name="ClassicDqDsSmoother");

    /// Default destructor
    virtual ~ClassicDqDsSmoother(){}

    /// ?
    int SmoothDqDs(SimpleCluster&) const;

  protected:

    void _Configure_(const Config_t &pset);
    size_t _slider_window;                             ///< ?

  };
  
  /**
     \class cosmictag::ClassicHitOrdererFactory
  */
  class ClassicDqDsSmootherFactory : public DqDsSmootherFactoryBase {
  public:
    /// ctor
    ClassicDqDsSmootherFactory() { DqDsSmootherFactory::get().add_factory("ClassicDqDsSmoother",this); }
    /// dtor
    ~ClassicDqDsSmootherFactory() {}
    /// creation method
    BaseDqDsSmootherAlgo* create(const std::string instance_name) { return new ClassicDqDsSmoother(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group 
