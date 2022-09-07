////////////////////////////////////////////////////////////////////////
// Class:       ToyNtupleCreator
// Plugin Type: analyzer (art v3_03_01)
// File:        ToyNtupleCreator_module.cc
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

//root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//local includes

namespace tutorial {
   class ToyNtupleCreator;
}


class tutorial::ToyNtupleCreator : public art::EDAnalyzer {
   public:
      explicit ToyNtupleCreator(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      ToyNtupleCreator(ToyNtupleCreator const&) = delete;
      ToyNtupleCreator(ToyNtupleCreator&&) = delete;
      ToyNtupleCreator& operator=(ToyNtupleCreator const&) = delete;
      ToyNtupleCreator& operator=(ToyNtupleCreator&&) = delete;

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

      std::string t_Mode;
      std::string t_CCNC;
      int t_NReconstructedTracks;

      /////////////////////////
      // Metadata for sample //
      /////////////////////////

      int m_NEvents;
      double m_POT = 0; //total POT of the sample

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      std::string f_GeneratorLabel;
      std::string f_PFParticleLabel;
      std::string f_TrackLabel;
      std::string f_POTSummaryLabel;
      bool f_Debug = false;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

tutorial::ToyNtupleCreator::ToyNtupleCreator(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_GeneratorLabel(p.get<std::string>("GeneratorLabel")),
   f_PFParticleLabel(p.get<std::string>("PFParticleLabel")),
   f_TrackLabel(p.get<std::string>("TrackLabel")),
   f_POTSummaryLabel(p.get<std::string>("POTSummaryLabel")),
   f_Debug(p.get<bool>("Debug",false))
{

}

void tutorial::ToyNtupleCreator::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "New Event" << std::endl;

   // General Event Info
   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   //begin by resetting everything
   t_Weight = 1.0;
   t_Mode = "NONE";
   t_CCNC = "NONE";
   t_NReconstructedTracks=0;

   // Get the generator truth information
   art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
   std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;

   if(!e.getByLabel(f_GeneratorLabel,Handle_MCTruth))  
      throw cet::exception("ToyNtupleCreator") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  

   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){

      simb::MCNeutrino Nu = theMCTruth->GetNeutrino();

      int mode = Nu.Mode();
      int ccnc = Nu.CCNC();

      if(ccnc == 0) t_CCNC = "CC";
      else t_CCNC = "NC";

      if(mode == 0) t_Mode = "QEL";
      else if(mode == 1) t_Mode = "RES";
      else if(mode == 2) t_Mode = "DIS";
      else if(mode == 3) t_Mode = "COH";
      else if(mode == 5) t_Mode = "ElectronScattering";
      else if(mode == 10) t_Mode = "MEC";
      else if(mode == 11) t_Mode = "Diffractive";
      else t_Mode = "Other";
   }

   // Load the reco'd tracks from the event 
   art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
   std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
   art::Handle<std::vector<recob::Track>> Handle_Track;
   std::vector<art::Ptr<recob::Track>> Vect_Track;

   if(!e.getByLabel(f_PFParticleLabel,Handle_PFParticle)) 
      throw cet::exception("ToyNtupleCreator") << "No PFParticle Data Products Found!" << std::endl;

   if(!e.getByLabel(f_TrackLabel,Handle_Track)) 
      throw cet::exception("ToyNtupleCreator") << "No Track Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);
   art::fill_ptr_vector(Vect_Track,Handle_Track);
   art::FindManyP<recob::Track> Assoc_PFParticleTrack = art::FindManyP<recob::Track>(Vect_PFParticle,e,f_TrackLabel);

   // Get the ID of the neutrino candidate
   size_t neutrinoid = -9999;
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle)
      if(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))
         neutrinoid = pfp->Self();

   // Get the number of tracks in the neutrino slice     
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
      if(pfp->Parent() != neutrinoid) continue;
      std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack.at(pfp.key());
      if(pfpTracks.size() != 1) continue;
      t_NReconstructedTracks++;
   }

   OutputTree->Fill();
}

///////////////////////////////////////////////////////////////	

void tutorial::ToyNtupleCreator::beginJob(){

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
   OutputTree->Branch("CCNC",&t_CCNC);
   OutputTree->Branch("Mode",&t_Mode);
   OutputTree->Branch("Weight",&t_Weight);
   OutputTree->Branch("NReconstructedTracks",&t_NReconstructedTracks);

   //////////////////////////////////////////
   //             Metadata Tree	           //
   //////////////////////////////////////////

   m_NEvents=0;
   m_POT=0;

   MetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");
   MetaTree->Branch("NEvents",&m_NEvents);
   MetaTree->Branch("POT",&m_POT);
}

void tutorial::ToyNtupleCreator::endJob()
{
   MetaTree->Fill();
}

void tutorial::ToyNtupleCreator::beginSubRun(const art::SubRun& sr)
{
   if(f_Debug) std::cout << "Getting Subrun POT Info" << std::endl;
   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(f_POTSummaryLabel,POTHandle)) m_POT += POTHandle->totpot;	
}

void tutorial::ToyNtupleCreator::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(tutorial::ToyNtupleCreator)
