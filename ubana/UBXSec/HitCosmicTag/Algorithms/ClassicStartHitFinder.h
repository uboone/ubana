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

#ifndef CLASSICSTARTHITFINDER_H
#define CLASSICSTARTHITFINDER_H

#include <iostream>
#include "../Base/BaseStartHitFinderAlgo.h"
#include "../Base/StartHitFinderFactory.h"
#include "../Base/HitCosmicTagException.h"

#include <TVector3.h>

namespace cosmictag {
  /**
     \class ClassicHitOrderer
     User custom analysis class
   */
  class ClassicStartHitFinder : public BaseStartHitFinderAlgo {
  
  public:

    /// Default constructor
    ClassicStartHitFinder(const std::string name="ClassicStartHitFinder");

    /// Default destructor
    virtual ~ClassicStartHitFinder(){}

    /// ?
    int FindStartHit(SimpleCluster&, SimpleHit&) const;

    void PrintConfig() const;

  protected:

    void _Configure_(const Config_t &pset);

    double _max_allowed_hit_distance;         ///< ?
  };
  
  /**
     \class cosmictag::ClassicHitOrdererFactory
  */
  class ClassicStartHitFinderFactory : public StartHitFinderFactoryBase {
  public:
    /// ctor
    ClassicStartHitFinderFactory() { StartHitFinderFactory::get().add_factory("ClassicStartHitFinder",this); }
    /// dtor
    ~ClassicStartHitFinderFactory() {}
    /// creation method
    BaseStartHitFinderAlgo* create(const std::string instance_name) { return new ClassicStartHitFinder(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group 
