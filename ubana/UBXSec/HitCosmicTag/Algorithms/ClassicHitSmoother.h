/**
 * \file ClassicHitSmoother.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class ClassicHitSmoother
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/

#ifndef ClassicHitSmoother_H
#define ClassicHitSmoother_H

#include <iostream>
#include "../Base/BaseHitSmootherAlgo.h"
#include "../Base/HitSmootherFactory.h"
#include "../Base/HitCosmicTagException.h"
#include "../Base/Tools.h"

#include <TVector3.h>

namespace cosmictag {
  /**
     \class ClassicHitSmoother
     User custom analysis class
   */
  class ClassicHitSmoother : public BaseHitSmootherAlgo {
  
  public:

    /// Default constructor
    ClassicHitSmoother(const std::string name="ClassicHitSmoother");

    /// Default destructor
    virtual ~ClassicHitSmoother(){}

    /// ?
    int Smooth(SimpleCluster&) const;

  protected:

    void _Configure_(const Config_t &pset);
    size_t _slider_window;                             ///< ?
  };
  
  /**
     \class cosmictag::ClassicHitSmootherFactory
  */
  class ClassicHitSmootherFactory : public HitSmootherFactoryBase {
  public:
    /// ctor
    ClassicHitSmootherFactory() { HitSmootherFactory::get().add_factory("ClassicHitSmoother",this); }
    /// dtor
    ~ClassicHitSmootherFactory() {}
    /// creation method
    BaseHitSmootherAlgo* create(const std::string instance_name) { return new ClassicHitSmoother(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group 
