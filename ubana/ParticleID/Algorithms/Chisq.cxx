/**
 * Handle to Chisquare code
 */

#include "Chisq.h"

namespace particleid{

  std::vector<float> Chisquare::getChisq( art::Ptr< anab::Calorimetry > caloObj, fhicl::ParameterSet const &p){

    std::vector<float> ChisqValues;

    //pid::Chi2PIDAlg fChiAlg;
    pid::Chi2PIDAlg fChiAlg(p);

    anab::ParticleID particleIdObj;

    fChiAlg.DoParticleID(caloObj, particleIdObj);
 
    ChisqValues.push_back(particleIdObj.Chi2Muon());
    ChisqValues.push_back(particleIdObj.Chi2Proton());
    ChisqValues.push_back(particleIdObj.Chi2Pion());
    ChisqValues.push_back(particleIdObj.Chi2Kaon());

    return ChisqValues;

  }
}
