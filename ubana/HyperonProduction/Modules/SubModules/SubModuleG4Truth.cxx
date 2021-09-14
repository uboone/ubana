#ifndef _SubModuleG4Truth_cxx_
#define _SubModuleG4Truth_cxx_

#include "SubModuleG4Truth.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,std::string g4label){

   if(!e.getByLabel(g4label,Handle_G4)) 
      throw cet::exception("SubModuleG4Truth") << "No Geant4 Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_G4,Handle_G4);

   // Create map between particle ID and Particles
   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,fhicl::ParameterSet pset) :
 SubModuleG4Truth(e,pset.get<std::string>("G4ModuleLabel","largeant"))
{
   SetNeutronScatterThresholds(pset.get<double>("NeutronScatterProtonThresh",0.15),pset.get<double>("NeutronScatterPionThresh",0.05));
   SetDecayThresholds(pset.get<double>("DecayProtonThresh",0.0),pset.get<double>("DecayPionThresh",0.0));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetParticleLists(){

   Primary_IDs.clear();
   Daughter_IDs.clear();
   SigmaZero_Daughter_IDs.clear();
   Kaon_Daughter_IDs.clear();

   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){

      // Skip any particles not produced at PV
      if(g4p->Mother() != 0)  continue;

      Primary_IDs.push_back(g4p->TrackId());

      std::vector<int> IDs = GetChildIDs(g4p);

      if(isHyperon(g4p->PdgCode())){

         if(g4p->PdgCode() == 3212 && g4p->EndProcess() == "Decay"){
            SigmaZero_Daughter_IDs.insert(SigmaZero_Daughter_IDs.begin(),IDs.begin(),IDs.end());                         

         }

         else if(g4p->EndProcess() == "Decay")
            Daughter_IDs.insert(Daughter_IDs.begin(),IDs.begin(),IDs.end());                

      }		

      else if(isKaon(g4p->PdgCode()) && g4p->EndProcess() == "Decay")
         Kaon_Daughter_IDs.insert(Kaon_Daughter_IDs.begin(),IDs.begin(),IDs.end());                     

   }


   // Get the SigmaZero daughters  
   for(size_t i_d=0;i_d<SigmaZero_Daughter_IDs.size();i_d++){

      if(partByID.find(SigmaZero_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[SigmaZero_Daughter_IDs.at(i_d)];

      if(part->PdgCode() == 3122){
         std::vector<int> IDs = GetChildIDs(part);
         Daughter_IDs.insert(Daughter_IDs.begin(),IDs.begin(),IDs.end());                     
      }     

   } 

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

G4Truth SubModuleG4Truth::GetG4Info(){

   GetParticleLists();

   // Clear everything
   theTruth.IsHyperon = false;
   theTruth.IsLambda = false;
   theTruth.IsLambdaCharged = false;
   theTruth.Lepton.clear();
   theTruth.Hyperon.clear();
   theTruth.PrimaryNucleon.clear();
   theTruth.PrimaryPion.clear();
   theTruth.PrimaryKaon.clear();
   theTruth.Decay.clear();
   theTruth.KaonDecay.clear();
   theTruth.SigmaZeroDecayPhoton.clear();
   theTruth.SigmaZeroDecayLambda.clear();

   GetPrimaryParticles();

   if(SigmaZero_Daughter_IDs.size()) GetSigmaZeroDecay(); 
   if(theTruth.Hyperon.size() || theTruth.SigmaZeroDecayLambda.size()) GetHyperonDecay();
   if(theTruth.PrimaryKaon.size()) GetKaonDecay();


   if(theTruth.Hyperon.size()) theTruth.IsHyperon = true;
   if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3122) theTruth.IsLambda = true;
   if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3212) theTruth.IsSigmaZero = true;
   if(theTruth.Hyperon.size() && theTruth.PrimaryKaon.size()) theTruth.IsAssociatedHyperon = true;

   if(theTruth.Decay.size() == 2){
      if(theTruth.IsLambda && theTruth.Decay.at(0).PDG == 2212 && theTruth.Decay.at(1).PDG == -211) theTruth.IsLambdaCharged = true;
      if(theTruth.IsLambda && theTruth.Decay.at(1).PDG == 2212 && theTruth.Decay.at(0).PDG == -211) theTruth.IsLambdaCharged = true; 
   }

   theTruth.HasNeutronScatter = FindNeutronScatter();   

   return theTruth;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetPrimaryParticles(){

   for(size_t i_p=0;i_p<Primary_IDs.size();i_p++){

      // Geant does not always keep everything it simulates, make sure particle is in list of IDs (crashes otherwise!)
      if(partByID.find(Primary_IDs[i_p]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs.at(i_p)];

      // Anything with very large pdg code is a nucleus, skip these
      if(part->PdgCode() > 10000) continue;

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 1;

      //hyperon produced at primary vertex
      if(isHyperon(part->PdgCode())){
         theTruth.Hyperon.push_back(P);
         theTruth.IsHyperon = true;	
      }

      if(isLepton(part->PdgCode()) || isNeutrino(part->PdgCode())) theTruth.Lepton.push_back(P);
      if(isNucleon(part->PdgCode())) theTruth.PrimaryNucleon.push_back(P);
      if(isPion(part->PdgCode())) theTruth.PrimaryPion.push_back(P);
      if(isKaon(part->PdgCode())) theTruth.PrimaryKaon.push_back(P);

   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetHyperonDecay(){

   if(!Daughter_IDs.size()) return;

   for(size_t i_d=0;i_d<Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 2;
      theTruth.Decay.push_back(P);     

   }

   theTruth.DecayVertex = TVector3(theTruth.Decay.at(0).StartX,theTruth.Decay.at(0).StartY,theTruth.Decay.at(0).StartZ);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetKaonDecay(){

   for(size_t i_d=0;i_d<Kaon_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Kaon_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Kaon_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 4;
      theTruth.Decay.push_back(P);     
   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetSigmaZeroDecay(){

   for(size_t i_d=0;i_d<SigmaZero_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(SigmaZero_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[SigmaZero_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 5;

      if(part->PdgCode() == 3122)
         theTruth.SigmaZeroDecayLambda.push_back(P);     
      else if(part->PdgCode() == 22)
         theTruth.SigmaZeroDecayPhoton.push_back(P);     

   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> SubModuleG4Truth::GetChildIDs(art::Ptr<simb::MCParticle> g4p,bool IsNeutron){

   std::vector<int> DecayProduct_IDs;

   if(g4p->EndProcess() != "Decay" && !IsNeutron) return DecayProduct_IDs;

   for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[g4p->Daughter(i_d)];

      if(part->PdgCode() > 10000) continue;

      double X = part->Position().X();
      double Y = part->Position().Y();
      double Z = part->Position().Z();

      double EndX = g4p->EndPosition().X();
      double EndY = g4p->EndPosition().Y();
      double EndZ = g4p->EndPosition().Z();

      if(X != EndX || Y != EndY || Z != EndZ) continue;

      DecayProduct_IDs.push_back(g4p->Daughter(i_d));
   }

   return DecayProduct_IDs;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int SubModuleG4Truth::GetOrigin(int trackid){

   if(std::find(Primary_IDs.begin(),Primary_IDs.end(),trackid) != Primary_IDs.end()) return 1;
   else if(std::find(Daughter_IDs.begin(),Daughter_IDs.end(),trackid) != Daughter_IDs.end()) return 2;
   else if(std::find(Kaon_Daughter_IDs.begin(),Kaon_Daughter_IDs.end(),trackid) != Kaon_Daughter_IDs.end()) return 4;
   else if(std::find(SigmaZero_Daughter_IDs.begin(),SigmaZero_Daughter_IDs.end(),trackid) != SigmaZero_Daughter_IDs.end()) return 5;
   else return 3;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::FindNeutronScatter(){

   std::vector<int> Neutron_IDs;

   for(size_t i=0;i<Primary_IDs.size();i++){

      if(partByID.find(Primary_IDs[i]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs[i]];

      if(part->PdgCode() != 2112) continue;

      if(part->EndProcess() == "neutronInelastic") Neutron_IDs.push_back(part->TrackId()); 

   }

   for(size_t i=0;i<Neutron_IDs.size();i++){
      art::Ptr<simb::MCParticle> part = partByID[Neutron_IDs[i]];
      std::vector<int> NeutronDaughter_IDs = GetChildIDs(part,true);      

      int nProtons = 0;
      int nPions = 0;

      for(size_t i_d=0;i_d<NeutronDaughter_IDs.size();i_d++){

         // Look for protons and charged pions
         art::Ptr<simb::MCParticle> part2 = partByID[NeutronDaughter_IDs[i_d]];


         double P = sqrt(part2->Px()*part2->Px() + part2->Py()*part2->Py() + part2->Pz()*part2->Pz());

         if(part2->PdgCode() == 2212 && P > NeutronScatterProtonThresh) nProtons++; 
         if(abs(part2->PdgCode()) == 211 && P > NeutronScatterPionThresh) nPions++; 

      }

      // Flag event as containing neutron scatter if there are either 2+ protons, 2+ charged pions or 1+ of each
      if(nProtons >= 2 || nPions >= 2 || (nProtons >= 1 && nPions >= 1)) return true;

   }

   return false;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetNeutronScatterThresholds(double neutronscatterprotonthresh,double neutronscatterpionthresh){

   NeutronScatterProtonThresh = neutronscatterprotonthresh;
   NeutronScatterPionThresh = neutronscatterpionthresh;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetDecayThresholds(double decayprotonthresh,double decaypionthresh){

   DecayProtonThresh = decayprotonthresh;
   DecayPionThresh = decaypionthresh;
}

#endif
