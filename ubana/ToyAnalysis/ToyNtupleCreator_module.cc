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

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      // Output trees
      TTree* OutputTree;
      TTree* MetaTree;

      // Basic event info
      unsigned int t_EventID;
      int t_run,t_subrun,t_event;

      double t_Weight;

      /////////////////////////
      // Metadata for sample //
      /////////////////////////

      int m_NEvents;
      double m_POT = 0; //total POT of the sample

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      std::string f_GeneratorLabel;
      std::string f_TrackLabel;
      std::string f_POTSummaryLabel;
      bool f_Debug = false;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::HyperonNtuples::HyperonNtuples(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_GeneratorLabel(p.get<std::string>("GeneratorLabel")),
   f_TrackLabel(p.get<std::string>("TrackLabel")),
   f_POTSummaryLabel(p.get<std::string>("POTSummaryLabel")),
   f_Debug(p.get<bool>("Debug",false))
{

}

void hyperon::HyperonNtuples::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "New Event" << std::endl;

   // General Event Info
   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   //begin by resetting everything
   t_Weight = 1.0;


   OutputTree->Fill();
}

///////////////////////////////////////////////////////////////	

void hyperon::HyperonNtuples::beginJob(){

   if(f_Debug) std::cout << "Creating TFileService and Setting Up Trees" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //             Output Tree	           //
   //////////////////////////////////////////

   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");
   OutputTree->Branch("EventID",&t_EventID);
   OutputTree->Branch("run",&t_run);
   OutputTree->Branch("subrun",&t_subrun);
   OutputTree->Branch("event",&t_event);
   OutputTree->Branch("Weight",&t_Weight);

   //////////////////////////////////////////
   //             Metadata Tree	           //
   //////////////////////////////////////////

   m_NEvents=0;
   m_POT=0;

   MetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");
   MetaTree->Branch("NEvents",&m_NEvents);
   MetaTree->Branch("POT",&m_POT);
}

void hyperon::HyperonNtuples::endJob()
{
   MetaTree->Fill();
}

void hyperon::HyperonNtuples::beginSubRun(const art::SubRun& sr)
{
   if(f_Debug) std::cout << "Getting Subrun POT Info" << std::endl;
   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(f_POTSummaryLabel,POTHandle)) m_POT += POTHandle->totpot;	
}

void hyperon::HyperonNtuples::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::HyperonNtuples)
