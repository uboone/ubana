/**
 * \file LocalLinearityCalculatorFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class LocalLinearityCalculatorFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef LOCALLINEARITYCALCULATORFACTORY_H
#define LOCALLINEARITYCALCULATORFACTORY_H

#include <iostream>
#include <map>
#include "BaseLocalLinearityCalculatorAlgo.h"


namespace cosmictag {

  /**
     \class LocalLinearityCalculatorFactoryBase
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class LocalLinearityCalculatorFactoryBase {
  public:
    /// Default ctor
    LocalLinearityCalculatorFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~LocalLinearityCalculatorFactoryBase(){}
    /// Abstract constructor method
    virtual BaseLocalLinearityCalculatorAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class LocalLinearityCalculatorFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class LocalLinearityCalculatorFactory {
  private:
    /// Default ctor, shouldn't be used
    LocalLinearityCalculatorFactory() {}
  public:
    /// Default dtor
    ~LocalLinearityCalculatorFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static LocalLinearityCalculatorFactory& get()
    { if(!_me) _me = new LocalLinearityCalculatorFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::LocalLinearityCalculatorFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseLocalLinearityCalculatorAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,cosmictag::LocalLinearityCalculatorFactoryBase*> _factory_map;
    /// Static self
    static LocalLinearityCalculatorFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

