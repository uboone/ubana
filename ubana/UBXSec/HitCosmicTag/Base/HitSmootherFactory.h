/**
 * \file HitSmootherFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class HitSmootherFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HitSmootherFACTORY_H
#define HitSmootherFACTORY_H

#include <iostream>
#include <map>
#include "BaseHitSmootherAlgo.h"


namespace cosmictag {

  /**
     \class HitSmootherFactoryBase
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class HitSmootherFactoryBase {
  public:
    /// Default ctor
    HitSmootherFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~HitSmootherFactoryBase(){}
    /// Abstract constructor method
    virtual BaseHitSmootherAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class HitSmootherFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class HitSmootherFactory {
  private:
    /// Default ctor, shouldn't be used
    HitSmootherFactory() {}
  public:
    /// Default dtor
    ~HitSmootherFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static HitSmootherFactory& get()
    { if(!_me) _me = new HitSmootherFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::HitSmootherFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseHitSmootherAlgo* create(const std::string name, const std::string instance_name) {
      auto iter = _factory_map.find(name);
      if(iter == _factory_map.end() || !((*iter).second)) {
	      std::cerr << "Found no registered class " << name << std::endl;
	      return nullptr;
      }
      auto ptr = (*iter).second->create(instance_name);
      return ptr;
    }

  private:
    /// Static factory container
    std::map<std::string,cosmictag::HitSmootherFactoryBase*> _factory_map;
    /// Static self
    static HitSmootherFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

