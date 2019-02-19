#ifndef CHISQ_H
#define CHISQ_H

#include <vector>
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "larana/ParticleIdentification/Chi2PIDAlg.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

namespace particleid{

  class Chisquare{

    public:
      
      std::vector<double> getChisq( art::Ptr< anab::Calorimetry > caloObj, fhicl::ParameterSet const &p);

  };

}

#endif
