#ifndef HITCOSMICTAG_BASEALGORITHM_CXX
#define HITCOSMICTAG_BASEALGORITHM_CXX

#include "BaseAlgorithm.h"
#include "HitCosmicTagException.h"
namespace cosmictag {

  BaseAlgorithm::BaseAlgorithm(const AlgoType type,const std::string name)
    : _type(type)
    , _name(name)
  {}

  AlgoType BaseAlgorithm::AlgorithmType() const
  { return _type; }

  void BaseAlgorithm::Configure(const Config_t &pset)
  {
    //this->set_verbosity((msg::Level_t)(pset.get<unsigned int>("Verbosity",(unsigned int)(msg::kNORMAL))));
    this->_Configure_(pset);
  }


  const std::string& BaseAlgorithm::AlgorithmName() const
  { return _name; }
  
}

#endif
