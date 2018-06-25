/**
 * \file StartHitFinderFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class StartHitFinderFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef STARTHITFINDERFACTORY_H
#define STARTHITFINDERFACTORY_H

#include <iostream>
#include <map>
#include "BaseStartHitFinderAlgo.h"


namespace cosmictag {

  /**
     \class BaseStartHitFinderAlgo
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class StartHitFinderFactoryBase {
  public:
    /// Default ctor
    StartHitFinderFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~StartHitFinderFactoryBase(){}
    /// Abstract constructor method
    virtual BaseStartHitFinderAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class StartHitFinderFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class StartHitFinderFactory {
  private:
    /// Default ctor, shouldn't be used
    StartHitFinderFactory() {}
  public:
    /// Default dtor
    ~StartHitFinderFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static StartHitFinderFactory& get()
    { if(!_me) _me = new StartHitFinderFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::StartHitFinderFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseStartHitFinderAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,cosmictag::StartHitFinderFactoryBase*> _factory_map;
    /// Static self
    static StartHitFinderFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

