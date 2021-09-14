#ifndef _SubModuleGeneratorTruth_cxx_
#define _SubModuleGeneratorTruth_cxx_

#include "SubModuleGeneratorTruth.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleGeneratorTruth::SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset){

   if(!e.getByLabel(pset.get<std::string>("GeneratorModuleLabel","generator"),Handle_MCTruth))  
      throw cet::exception("GeneratorTruthSubModule") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

GeneratorTruth SubModuleGeneratorTruth::GetGeneratorTruth(){

   if(!Vect_MCTruth.size()){
      std::cout << "MCTruth vector is empty" << std::endl;
      return theTruth;
   }

   theTruth.NMCTruths = Vect_MCTruth.size();

   art::Ptr<simb::MCTruth> theMCTruth = Vect_MCTruth.at(0);

   simb::MCNeutrino Nu = theMCTruth->GetNeutrino();

   int mode = Nu.Mode();
   int ccnc = Nu.CCNC();

   if(ccnc == 0) theTruth.CCNC = "CC";
   else theTruth.CCNC = "NC";

   // NOTE: Hyperon events produced in GENIE use 0, NuWro uses 1095
   if(mode == 0) theTruth.Mode = "QEL";
   else if(mode == 1) theTruth.Mode = "RES";
   else if(mode == 2) theTruth.Mode = "DIS";
   else if(mode == 3) theTruth.Mode = "COH";
   else if(mode == 5) theTruth.Mode = "ElectronScattering";
   else if(mode == 10) theTruth.Mode = "MEC";
   else if(mode == 11) theTruth.Mode = "Diffractive";
   else if(mode == 1095) theTruth.Mode = "HYP";
   else theTruth.Mode = "Other";	

   for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){

      simb::MCParticle Part = theMCTruth->GetParticle(k_particles);

      if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) 
         theTruth.TruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());

      if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
         SimParticle P = MakeSimParticle(Part);
         P.Origin = 0;
         theTruth.Neutrino.push_back(P);
      }

      // If there is a hyperon i the final state in a QEL event, change mode to HYP
      if(isHyperon(Part.PdgCode()) && Part.StatusCode() == 1 && mode == 0) theTruth.Mode = "HYP";
   }

   return theTruth;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
