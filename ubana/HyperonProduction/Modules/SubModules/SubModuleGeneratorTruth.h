#ifndef _SubModuleGeneratorTruth_h_
#define _SubModuleGeneratorTruth_h_

#include <string>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

namespace hyperon {

struct GeneratorTruth {

   double Weight = 1.0;
   int NMCTruths;
   std::string Mode;
   std::vector<SimParticle> Neutrino;
   TVector3 TruePrimaryVertex;
   std::string CCNC;

};

class SubModuleGeneratorTruth {

public:

   SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset);

   GeneratorTruth GetGeneratorTruth();

private:

   art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
   std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;

   GeneratorTruth theTruth;

};

}

#endif
