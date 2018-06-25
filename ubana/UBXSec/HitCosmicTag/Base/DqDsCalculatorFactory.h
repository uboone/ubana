/**
 * \file DqDsCalculatorFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class DqDsCalculatorFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef DQDSCALCULATORFACTORY_H
#define DQDSCALCULATORFACTORY_H

#include <iostream>
#include <map>
#include "BaseDqDsCalculatorAlgo.h"


namespace cosmictag {

  /**
     \class DqDsCalculatorFactoryBase
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class DqDsCalculatorFactoryBase {
  public:
    /// Default ctor
    DqDsCalculatorFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~DqDsCalculatorFactoryBase(){}
    /// Abstract constructor method
    virtual BaseDqDsCalculatorAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class DqDsCalculatorFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class DqDsCalculatorFactory {
  private:
    /// Default ctor, shouldn't be used
    DqDsCalculatorFactory() {}
  public:
    /// Default dtor
    ~DqDsCalculatorFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static DqDsCalculatorFactory& get()
    { if(!_me) _me = new DqDsCalculatorFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::DqDsCalculatorFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseDqDsCalculatorAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,cosmictag::DqDsCalculatorFactoryBase*> _factory_map;
    /// Static self
    static DqDsCalculatorFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

