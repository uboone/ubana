/**
 * \file HitOrdererFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class HitOrdererFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITORDERERFACTORY_H
#define HITORDERERFACTORY_H

#include <iostream>
#include <map>
#include "BaseHitOrdererAlgo.h"


namespace cosmictag {

  /**
     \class HitOrdererFactoryBase
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class HitOrdererFactoryBase {
  public:
    /// Default ctor
    HitOrdererFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~HitOrdererFactoryBase(){}
    /// Abstract constructor method
    virtual BaseHitOrdererAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class HitOrdererFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class HitOrdererFactory {
  private:
    /// Default ctor, shouldn't be used
    HitOrdererFactory() {}
  public:
    /// Default dtor
    ~HitOrdererFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static HitOrdererFactory& get()
    { if(!_me) _me = new HitOrdererFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::HitOrdererFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseHitOrdererAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,cosmictag::HitOrdererFactoryBase*> _factory_map;
    /// Static self
    static HitOrdererFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

