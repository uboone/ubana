////////////////////////////////////////////////////////////////////////
// Class:       HyperonNtuples
// Plugin Type: analyzer (art v3_03_01)
// File:        HyperonNtuples_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

//root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//local includes

//objects and helpers
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

//algorithms
#include "ubana/HyperonProduction/Alg/ConnectednessHelper.h"

//submodules
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"

namespace hyperon {
   class HyperonNtuples;
}


class hyperon::HyperonNtuples : public art::EDAnalyzer {
   public:
      explicit HyperonNtuples(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HyperonNtuples(HyperonNtuples const&) = delete;
      HyperonNtuples(HyperonNtuples&&) = delete;
      HyperonNtuples& operator=(HyperonNtuples const&) = delete;
      HyperonNtuples& operator=(HyperonNtuples&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void FinishEvent();

      //check if event contains a reco'd muon, proton and pion from Lambda decay
      //records their positions in track vector if they exist
      //void StoreTrackTruth();

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      // Output trees
      TTree * OutputTree;
      TTree * MetaTree;

      // Basic event info
      unsigned int fEventID;
      int run,subrun,event;

      double t_Weight=1.0;

      // Generator/Geant4 truth info

      int t_NMCTruths=0;	

      std::string t_Mode; //interaction mode
      std::string t_CCNC; //charged current/neutral current

      bool t_IsLambda;
      bool t_IsHyperon;
      bool t_IsSigmaZero; 		
      bool t_IsLambdaCharged;
      bool t_IsAssociatedHyperon;
      bool t_HasNeutronScatter;
      bool t_IsSignal;	
      bool t_GoodReco;

      std::vector<SimParticle> t_Neutrino;
      std::vector<SimParticle> t_Lepton;
      std::vector<SimParticle> t_Hyperon;
      std::vector<SimParticle> t_PrimaryNucleon;
      std::vector<SimParticle> t_PrimaryPion;
      std::vector<SimParticle> t_PrimaryKaon; 
      std::vector<SimParticle> t_Decay; 
      std::vector<SimParticle> t_SigmaZeroDecayPhoton;
      std::vector<SimParticle> t_SigmaZeroDecayLambda;
      std::vector<SimParticle> t_KaonDecay;

      TVector3 t_TruePrimaryVertex;
      TVector3 t_DecayVertex;

      int t_NPrimaryDaughters;
      int t_NPrimaryTrackDaughters;
      int t_NPrimaryShowerDaughters;

      std::vector<RecoParticle> t_TrackPrimaryDaughters;
      std::vector<RecoParticle> t_ShowerPrimaryDaughters;   

      TVector3 t_RecoPrimaryVertex;

      ////////////////////////////
      //   Connectedness test   //
      ////////////////////////////

      std::vector<std::vector<int>> t_Conn_SeedIndexes_Plane0;
      std::vector<std::vector<int>> t_Conn_OutputIndexes_Plane0;
      std::vector<std::vector<int>> t_Conn_OutputSizes_Plane0;
      std::vector<std::vector<int>> t_Conn_SeedChannels_Plane0;
      std::vector<std::vector<int>> t_Conn_SeedTicks_Plane0;

      std::vector<std::vector<int>> t_Conn_SeedIndexes_Plane1;
      std::vector<std::vector<int>> t_Conn_OutputIndexes_Plane1;
      std::vector<std::vector<int>> t_Conn_OutputSizes_Plane1;
      std::vector<std::vector<int>> t_Conn_SeedChannels_Plane1;
      std::vector<std::vector<int>> t_Conn_SeedTicks_Plane1;

      std::vector<std::vector<int>> t_Conn_SeedIndexes_Plane2;
      std::vector<std::vector<int>> t_Conn_OutputIndexes_Plane2;
      std::vector<std::vector<int>> t_Conn_OutputSizes_Plane2;
      std::vector<std::vector<int>> t_Conn_SeedChannels_Plane2;
      std::vector<std::vector<int>> t_Conn_SeedTicks_Plane2;

      std::vector<std::string> t_SysDials;
      std::vector<std::vector<double>> t_SysWeights;

      /////////////////////////
      // Metadata for sample //
      /////////////////////////

      int fNEvents;
      int fNHyperons;
      int fNSignal;      
      int fNGoodReco;

      double fPOT = 0; //total POT of the sample

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;
      fhicl::ParameterSet f_Reco;
      std::string fWireLabel;

      std::string fWeightLabel;

      bool fIsData;

      std::string fPOTSummaryLabel;

      bool fDebug = false;

      ///////////////////////
      //      Objects      //
      ///////////////////////

      ConnectednessHelper Conn_Helper;

};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::HyperonNtuples::HyperonNtuples(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Reco(p.get<fhicl::ParameterSet>("Reco")),
   fWeightLabel(p.get<std::string>("WeightLabel","None")),
   Conn_Helper(p.get<bool>("DrawConnectedness",false))   // ,
{
   fIsData = p.get<bool>("IsData");

   fWireLabel = p.get<std::string>("WireLabel");
   fPOTSummaryLabel = p.get<std::string>("POTSummaryLabel");

   fDebug = p.get<bool>("Debug",false);
}

void hyperon::HyperonNtuples::analyze(art::Event const& e)
{
   if(fDebug) std::cout << "New Event" << std::endl;

   //begin by resetting everything

   t_Weight = 1.0;
   t_Mode = "NONE";
   t_NMCTruths = 0;
   t_IsHyperon = false;
   t_IsSigmaZero = false;
   t_IsLambda = false;
   t_IsLambdaCharged = false;
   t_IsAssociatedHyperon = false;
   t_HasNeutronScatter = false;
   t_IsSignal = false;	
   t_GoodReco = false;

   t_Neutrino.clear();
   t_Lepton.clear();
   t_Hyperon.clear();
   t_PrimaryNucleon.clear();
   t_PrimaryPion.clear();
   t_PrimaryKaon.clear();
   t_Decay.clear();
   t_KaonDecay.clear();
   t_SigmaZeroDecayPhoton.clear();
   t_SigmaZeroDecayLambda.clear();

   t_TruePrimaryVertex.SetXYZ(-1000,-1000,-1000);
   t_DecayVertex.SetXYZ(-1000,-1000,-1000);
 
   t_NPrimaryDaughters = 0; //number of primary daughters
   t_NPrimaryTrackDaughters=0; //num of track like primary daughters
   t_NPrimaryShowerDaughters=0; //num of shower like primary daughters

   t_TrackPrimaryDaughters.clear();
   t_ShowerPrimaryDaughters.clear();

   t_RecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

   t_Conn_SeedIndexes_Plane0.clear();
   t_Conn_OutputIndexes_Plane0.clear();
   t_Conn_OutputSizes_Plane0.clear();
   t_Conn_SeedChannels_Plane0.clear();
   t_Conn_SeedTicks_Plane0.clear();

   t_Conn_SeedIndexes_Plane1.clear();
   t_Conn_OutputIndexes_Plane1.clear();
   t_Conn_OutputSizes_Plane1.clear();
   t_Conn_SeedChannels_Plane1.clear();
   t_Conn_SeedTicks_Plane1.clear();

   t_Conn_SeedIndexes_Plane2.clear();
   t_Conn_OutputIndexes_Plane2.clear();
   t_Conn_OutputSizes_Plane2.clear();
   t_Conn_SeedChannels_Plane2.clear();
   t_Conn_SeedTicks_Plane2.clear();
 
   // General Event Info

   fEventID = e.id().event();
   run = e.run();
   subrun = e.subRun();
   event = e.event();

   // Event Generator Info

   if(!fIsData){

      if(fDebug) std::cout << "Getting EG Info" << std::endl;

      SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();
      
      t_Weight *= GenT.Weight;
      t_CCNC = GenT.CCNC;
      t_Mode = GenT.Mode;
      t_NMCTruths = GenT.NMCTruths;
      t_Neutrino = GenT.Neutrino;
      t_TruePrimaryVertex = GenT.TruePrimaryVertex;
      delete Generator_SM;
   }

   // G4 Info

   if(!fIsData){

      if(fDebug) std::cout << "Getting G4 Info" << std::endl;

      SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4);
      G4Truth G4T = G4_SM->GetG4Info();

      t_IsHyperon = G4T.IsHyperon;
      t_IsSigmaZero = G4T.IsSigmaZero;
      t_IsLambda = G4T.IsLambda;
      t_IsLambdaCharged = G4T.IsLambdaCharged;
      t_IsAssociatedHyperon = G4T.IsAssociatedHyperon;
      t_HasNeutronScatter = G4T.HasNeutronScatter;       
      t_Weight *= G4T.Weight;
      t_Lepton = G4T.Lepton;
      t_Hyperon = G4T.Hyperon;
      t_PrimaryNucleon = G4T.PrimaryNucleon;
      t_PrimaryPion = G4T.PrimaryPion;
      t_PrimaryKaon = G4T.PrimaryKaon;
      t_Decay = G4T.Decay;
      t_KaonDecay = G4T.KaonDecay;
      t_SigmaZeroDecayPhoton = G4T.SigmaZeroDecayPhoton;
      t_SigmaZeroDecayLambda = G4T.SigmaZeroDecayLambda;
      t_DecayVertex = G4T.DecayVertex;

      t_IsSignal = t_Neutrino.size() == 1 && t_Mode == "HYP" && t_Neutrino.at(0).PDG == -14 && t_IsLambdaCharged;
        
      delete G4_SM;
   }

   // Reconstructed Info

   if(fDebug) std::cout << "Getting Reconstructed Info" << std::endl;

   SubModuleReco* Reco_SM = new SubModuleReco(e,fIsData,f_Reco);
   Reco_SM->PrepareInfo();
   Reco_SM->SetIndices(t_IsSignal);
   RecoData RecoD =  Reco_SM->GetInfo();   

   t_NPrimaryDaughters = RecoD.NPrimaryDaughters;
   t_NPrimaryTrackDaughters = RecoD.NPrimaryTrackDaughters;
   t_NPrimaryShowerDaughters = RecoD.NPrimaryShowerDaughters;
   t_TrackPrimaryDaughters = RecoD.TrackPrimaryDaughters;
   t_ShowerPrimaryDaughters = RecoD.ShowerPrimaryDaughters;
   t_RecoPrimaryVertex = RecoD.RecoPrimaryVertex;

   // Results of connectedness test on different combinations of tracks

   if(fDebug) std::cout << "Performing Connectedness Tests" << std::endl;

   CTOutcome ConnData = Conn_Helper.PrepareAndTestEvent(e,fWireLabel,RecoD.TrackStarts);   

   t_Conn_SeedIndexes_Plane0 = ConnData.SeedIndexes_Plane0;
   t_Conn_OutputIndexes_Plane0 = ConnData.OutputIndexes_Plane0;
   t_Conn_OutputSizes_Plane0 = ConnData.OutputSizes_Plane0;
   t_Conn_SeedChannels_Plane0 = ConnData.SeedChannels_Plane0;
   t_Conn_SeedTicks_Plane0 = ConnData.SeedTicks_Plane0;
   t_Conn_SeedIndexes_Plane1 = ConnData.SeedIndexes_Plane1;
   t_Conn_OutputIndexes_Plane1 = ConnData.OutputIndexes_Plane1;
   t_Conn_OutputSizes_Plane1 = ConnData.OutputSizes_Plane1;
   t_Conn_SeedChannels_Plane1 = ConnData.SeedChannels_Plane1;
   t_Conn_SeedTicks_Plane1 = ConnData.SeedTicks_Plane1;
   t_Conn_SeedIndexes_Plane2 = ConnData.SeedIndexes_Plane2;
   t_Conn_OutputIndexes_Plane2 = ConnData.OutputIndexes_Plane2;
   t_Conn_OutputSizes_Plane2 = ConnData.OutputSizes_Plane2;
   t_Conn_SeedChannels_Plane2 = ConnData.SeedChannels_Plane2;
   t_Conn_SeedTicks_Plane2 = ConnData.SeedTicks_Plane2;

   // Systematics weights if requested

   if(fWeightLabel != "None"){

      // Try to get some systematics info
      art::Handle<std::vector<evwgh::MCEventWeight>> Handle_EventWeight;
      std::vector<art::Ptr<evwgh::MCEventWeight>> Vect_EventWeight;

      if(!e.getByLabel(fWeightLabel,Handle_EventWeight)) 
         throw cet::exception("HyperonNtuples") << "No EventWeight Found!" << std::endl;

      art::fill_ptr_vector(Vect_EventWeight,Handle_EventWeight);

      if(!Vect_EventWeight.size())
         throw cet::exception("HyperonNtuples") << "Weight vector empty!" << std::endl;

      std::map<std::string,std::vector<double>> theWeights = Vect_EventWeight.at(0)->fWeight;
      std::map<std::string,std::vector<double>>::iterator it;

      for(it = theWeights.begin(); it != theWeights.end();it++){
         t_SysDials.push_back(it->first);
         t_SysWeights.push_back(it->second);
      }

   }

   FinishEvent();

   delete Reco_SM;
}

///////////////////////////////////////////////////////////////	
// Finished processing event - update Metadata and fill tree //
///////////////////////////////////////////////////////////////

void hyperon::HyperonNtuples::FinishEvent(){

   if(fDebug) std::cout << "Finishing Event" << std::endl;

   OutputTree->Fill();

   fNEvents++;
   if(t_IsHyperon) fNHyperons++;
   if(t_IsSignal) fNSignal++;
   if(t_GoodReco) fNGoodReco++;

}

///////////////////////////////////////////////////////////////	

void hyperon::HyperonNtuples::beginJob(){

   if(fDebug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //             Output Tree	           //
   //////////////////////////////////////////

   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");

   OutputTree->Branch("IsData",&fIsData);
   OutputTree->Branch("EventID",&fEventID);
   OutputTree->Branch("run",&run);
   OutputTree->Branch("subrun",&subrun);
   OutputTree->Branch("event",&event);

   OutputTree->Branch("Weight",&t_Weight);
   OutputTree->Branch("Mode",&t_Mode);
   OutputTree->Branch("CCNC",&t_CCNC);
   OutputTree->Branch("NMCTruths",&t_NMCTruths);
   OutputTree->Branch("IsHyperon",&t_IsHyperon);
   OutputTree->Branch("IsSigmaZero",&t_IsSigmaZero);
   OutputTree->Branch("IsLambda",&t_IsLambda);
   OutputTree->Branch("IsLambdaCharged",&t_IsLambdaCharged);
   OutputTree->Branch("IsSignal",&t_IsSignal);
   OutputTree->Branch("GoodReco",&t_GoodReco);
   OutputTree->Branch("IsAssociatedHyperon",&t_IsAssociatedHyperon);
   OutputTree->Branch("HasNeutronScatter",&t_HasNeutronScatter);
   OutputTree->Branch("Neutrino","vector<SimParticle>",&t_Neutrino);

   OutputTree->Branch("Lepton","vector<SimParticle>",&t_Lepton);
   OutputTree->Branch("Hyperon","vector<SimParticle>",&t_Hyperon);
   OutputTree->Branch("PrimaryNucleon","vector<SimParticle>",&t_PrimaryNucleon);
   OutputTree->Branch("PrimaryPion","vector<SimParticle>",&t_PrimaryPion);
   OutputTree->Branch("PrimaryKaon","vector<SimParticle>",&t_PrimaryKaon);
   OutputTree->Branch("Decay","vector<SimParticle>",&t_Decay);
   OutputTree->Branch("SigmaZeroDecayPhoton","vector<SimParticle>",&t_SigmaZeroDecayPhoton);
   OutputTree->Branch("SigmaZeroDecayLambda","vector<SimParticle>",&t_SigmaZeroDecayLambda);
   OutputTree->Branch("KaonDecay","vector<SimParticle>",&t_KaonDecay);
   OutputTree->Branch("TruePrimaryVertex","TVector3",&t_TruePrimaryVertex);
   OutputTree->Branch("DecayVertex","TVector3",&t_DecayVertex); //position ot_ hyperon decay vertex

   OutputTree->Branch("RecoPrimaryVertex","TVector3",&t_RecoPrimaryVertex);
   OutputTree->Branch("NPrimaryTrackDaughters",&t_NPrimaryTrackDaughters);
   OutputTree->Branch("NPrimaryShowerDaughters",&t_NPrimaryShowerDaughters);
   OutputTree->Branch("TracklikePrimaryDaughters","vector<RecoParticle>",&t_TrackPrimaryDaughters);
   OutputTree->Branch("ShowerlikePrimaryDaughters","vector<RecoParticle>",&t_ShowerPrimaryDaughters);
   
   OutputTree->Branch("ConnSeedIndexes_Plane0",&t_Conn_SeedIndexes_Plane0);
   OutputTree->Branch("ConnOutputIndexes_Plane0",&t_Conn_OutputIndexes_Plane0);
   OutputTree->Branch("ConnOutputSizes_Plane0",&t_Conn_OutputSizes_Plane0);
   OutputTree->Branch("ConnSeedChannels_Plane0",&t_Conn_SeedChannels_Plane0);
   OutputTree->Branch("ConnSeedTicks_Plane0",&t_Conn_SeedTicks_Plane0);
   OutputTree->Branch("ConnSeedIndexes_Plane1",&t_Conn_SeedIndexes_Plane1);
   OutputTree->Branch("ConnOutputIndexes_Plane1",&t_Conn_OutputIndexes_Plane1);
   OutputTree->Branch("ConnOutputSizes_Plane1",&t_Conn_OutputSizes_Plane1);
   OutputTree->Branch("ConnSeedChannels_Plane1",&t_Conn_SeedChannels_Plane1);
   OutputTree->Branch("ConnSeedTicks_Plane1",&t_Conn_SeedTicks_Plane1);
   OutputTree->Branch("ConnSeedIndexes_Plane2",&t_Conn_SeedIndexes_Plane2);
   OutputTree->Branch("ConnOutputIndexes_Plane2",&t_Conn_OutputIndexes_Plane2);
   OutputTree->Branch("ConnOutputSizes_Plane2",&t_Conn_OutputSizes_Plane2);
   OutputTree->Branch("ConnSeedChannels_Plane2",&t_Conn_SeedChannels_Plane2);
   OutputTree->Branch("ConnSeedTicks_Plane2",&t_Conn_SeedTicks_Plane2);

   OutputTree->Branch("SysDials",&t_SysDials);
   OutputTree->Branch("SysWeights",&t_SysWeights);

   //////////////////////////////////////////
   //             Metadata Tree	           //
   //////////////////////////////////////////

   fNEvents=0;
   fNHyperons=0;
   fNSignal=0;
   fNGoodReco=0;

   fPOT=0;

   MetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");

   MetaTree->Branch("NEvents",&fNEvents);
   MetaTree->Branch("NHyperons",&fNHyperons);
   MetaTree->Branch("NSignal",&fNSignal);
   MetaTree->Branch("NGoodReco",&fNGoodReco);

   MetaTree->Branch("POT",&fPOT);

   if(fDebug) std::cout << "Finished begin job" << std::endl;

}



void hyperon::HyperonNtuples::endJob()
{
   MetaTree->Fill();
}

void hyperon::HyperonNtuples::beginSubRun(const art::SubRun& sr)
{
   if(fDebug) std::cout << "Getting Subrun POT Info" << std::endl;

   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(fPOTSummaryLabel,POTHandle)) fPOT += POTHandle->totpot;	
}

void hyperon::HyperonNtuples::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::HyperonNtuples)
