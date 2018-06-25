/**
 * \file DqDsSmootherFactoryy.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class DqDsSmootherFactory
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef DQDSSMOOTHERFACTORYFACTORY_H
#define DQDSSMOOTHERFACTORYFACTORY_H

#include <iostream>
#include <map>
#include "BaseDqDsSmootherAlgo.h"


namespace cosmictag {

  /**
     \class DqDsSmootherFactoryBase
     \brief Abstract base class for factory (to be implemented per cluster)
  */
  class DqDsSmootherFactoryBase {
  public:
    /// Default ctor
    DqDsSmootherFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~DqDsSmootherFactoryBase(){}
    /// Abstract constructor method
    virtual BaseDqDsSmootherAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class DqDsSmootherFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class DqDsSmootherFactory {
  private:
    /// Default ctor, shouldn't be used
    DqDsSmootherFactory() {}
  public:
    /// Default dtor
    ~DqDsSmootherFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static DqDsSmootherFactory& get()
    { if(!_me) _me = new DqDsSmootherFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, cosmictag::DqDsSmootherFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseDqDsSmootherAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,cosmictag::DqDsSmootherFactoryBase*> _factory_map;
    /// Static self
    static DqDsSmootherFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

