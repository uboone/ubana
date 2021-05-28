////////////////////////////////////////////////////////////////////////
// Class:       WireTreeMaker
// Plugin Type: analyzer (art v3_03_01)
// File:        WireTreeMaker_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <string>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// Root includes
#include "TTree.h"
#include "TVector3.h"

// Local includes
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"


namespace hyperon {
   class WireTreeMaker;
}


class hyperon::WireTreeMaker : public art::EDAnalyzer {
   public:
      explicit WireTreeMaker(fhicl::ParameterSet const& p);

      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      WireTreeMaker(WireTreeMaker const&) = delete;
      WireTreeMaker(WireTreeMaker&&) = delete;
      WireTreeMaker& operator=(WireTreeMaker const&) = delete;
      WireTreeMaker& operator=(WireTreeMaker&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;



      // Selected optional functions.
      void beginJob() override;
      void endJob() override;


      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);


   private:

      // General event info
      unsigned int t_EventID;
      int t_run,t_subrun,t_event;

      ///////////////////////////////
      //     Wire Signal Info      //
      ///////////////////////////////

      TTree *t_WireTree;

      // Track Start (in xyz and tick/channel space)
      std::vector<int> t_TrackStart_Channel_Plane0;
      std::vector<int> t_TrackStart_Time_Plane0;

      std::vector<int> t_TrackStart_Channel_Plane1;
      std::vector<int> t_TrackStart_Time_Plane1;

      std::vector<int> t_TrackStart_Channel_Plane2;
      std::vector<int> t_TrackStart_Time_Plane2;

      std::vector<double> t_TrackStart_X;
      std::vector<double> t_TrackStart_Y;
      std::vector<double> t_TrackStart_Z;

      std::vector<double> t_TrackDir_X;
      std::vector<double> t_TrackDir_Y;
      std::vector<double> t_TrackDir_Z;


      // Calo Start (in xyz and tick/channel space)
      std::vector<int> t_CaloStart_Channel_Plane0;
      std::vector<int> t_CaloStart_Time_Plane0;
      std::vector<double> t_CaloStart_X_Plane0;
      std::vector<double> t_CaloStart_Y_Plane0;
      std::vector<double> t_CaloStart_Z_Plane0;

      std::vector<int> t_CaloStart_Channel_Plane1;
      std::vector<int> t_CaloStart_Time_Plane1;
      std::vector<double> t_CaloStart_X_Plane1;
      std::vector<double> t_CaloStart_Y_Plane1;
      std::vector<double> t_CaloStart_Z_Plane1;

      std::vector<int> t_CaloStart_Channel_Plane2;
      std::vector<int> t_CaloStart_Time_Plane2;
      std::vector<double> t_CaloStart_X_Plane2;
      std::vector<double> t_CaloStart_Y_Plane2;
      std::vector<double> t_CaloStart_Z_Plane2;


      // Wire Signals
      std::vector<int> t_Wire_Channel_Plane0;
      std::vector<int> t_Wire_Tick_Plane0;
      std::vector<double> t_Wire_Signal_Plane0;

      std::vector<int> t_Wire_Channel_Plane1;
      std::vector<int> t_Wire_Tick_Plane1;
      std::vector<double> t_Wire_Signal_Plane1;

      std::vector<int> t_Wire_Channel_Plane2;
      std::vector<int> t_Wire_Tick_Plane2;
      std::vector<double> t_Wire_Signal_Plane2;


      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      bool fDebug;

      // Producer module labels
      std::string fTrackLabel;
      std::string fPFParticleLabel;
      std::string fCaloLabel;
      std::string fWireLabel;


};



////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////



hyperon::WireTreeMaker::WireTreeMaker(fhicl::ParameterSet const& p)
   : EDAnalyzer{p}   // ,
   // More initializers here.
{

   fDebug = p.get<bool>("Debug","false");

   // Module labels
   fPFParticleLabel = p.get<std::string>("PFParticleLabel");
   fTrackLabel = p.get<std::string>("TrackLabel");
   fCaloLabel = p.get<std::string>("CaloLabel");
   fWireLabel = p.get<std::string>("WireLabel");

}

void hyperon::WireTreeMaker::analyze(art::Event const& e)
{

   if(fDebug) std::cout << "New Event" << std::endl;

   ////////////////////////////////
   //  Reset All of the Vectors  //
   ////////////////////////////////

   t_TrackStart_Channel_Plane0.clear();
   t_TrackStart_Time_Plane0.clear();

   t_TrackStart_Channel_Plane1.clear();
   t_TrackStart_Time_Plane1.clear();

   t_TrackStart_Channel_Plane2.clear();
   t_TrackStart_Time_Plane2.clear();

   t_TrackStart_X.clear();
   t_TrackStart_Y.clear();
   t_TrackStart_Z.clear();

   t_TrackDir_X.clear();
   t_TrackDir_Y.clear();
   t_TrackDir_Z.clear();

   t_CaloStart_Channel_Plane0.clear();
   t_CaloStart_Time_Plane0.clear();
   t_CaloStart_X_Plane0.clear();
   t_CaloStart_Y_Plane0.clear();
   t_CaloStart_Z_Plane0.clear();

   t_CaloStart_Channel_Plane1.clear();
   t_CaloStart_Time_Plane1.clear();
   t_CaloStart_X_Plane1.clear();
   t_CaloStart_Y_Plane1.clear();
   t_CaloStart_Z_Plane1.clear();

   t_CaloStart_Channel_Plane2.clear();
   t_CaloStart_Time_Plane2.clear();
   t_CaloStart_X_Plane2.clear();
   t_CaloStart_Y_Plane2.clear();
   t_CaloStart_Z_Plane2.clear();


   t_Wire_Channel_Plane0.clear();
   t_Wire_Tick_Plane0.clear();
   t_Wire_Signal_Plane0.clear();

   t_Wire_Channel_Plane1.clear();
   t_Wire_Tick_Plane1.clear();
   t_Wire_Signal_Plane1.clear();

   t_Wire_Channel_Plane2.clear();
   t_Wire_Tick_Plane2.clear();
   t_Wire_Signal_Plane2.clear();

   ////////////////////////////
   //  Event ID Information  //
   ////////////////////////////

   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   /////////////////////////////////////////////////////////
   //  Obtain Start Positions of Tracks in the Hierarchy  //
   /////////////////////////////////////////////////////////

   if(fDebug) std::cout << "Getting Reco'd Particles" << std::endl;

   // Setup handles
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
   art::Handle<std::vector<recob::Track>> trackHandle;
   art::Handle<std::vector<anab::Calorimetry>> caloHandle;

   std::vector<art::Ptr<recob::Track>> trackVect;
   std::vector<art::Ptr<recob::PFParticle>> pfparticleVect;
   std::vector<art::Ptr<anab::Calorimetry>> caloVect;

   // Fill PFP vector
   if(e.getByLabel(fPFParticleLabel,pfparticleHandle)){
      art::fill_ptr_vector(pfparticleVect,pfparticleHandle);
   }

   // If PFP vector is empty, add blank entry to tree and go to next event
   if(!pfparticleVect.size()) {
      t_WireTree->Fill();
      return;
   }


   // Fill track vector
   if(e.getByLabel(fTrackLabel,trackHandle)) art::fill_ptr_vector(trackVect,trackHandle);
   else
      std::cout << "Track handle not setup" << std::endl;


   // Setup Assns

   // Tracks Assoc with PFPs
   art::FindManyP<recob::Track> trackAssoc(pfparticleVect,e,fTrackLabel);

   // Calo Assoc to Tracks
   art::FindManyP<anab::Calorimetry> caloTrackAssoc(trackVect,e,fCaloLabel);


   // Go through the list of pandora PFP's, find the reconstructed neutrino
   size_t neutrinoID = 99999;

   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect)
      if((pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12 )))
         neutrinoID = pfp->Self();


   // Collect start points of tracks
   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){

      // Only store information from particles in the hierarchy
      if(pfp->Parent() != neutrinoID) continue;

      std::vector<art::Ptr<recob::Track>> pfpTracks = trackAssoc.at(pfp.key());

      if(pfpTracks.size() == 1){

         art::Ptr<recob::Track> trk = pfpTracks.at(0);

         // Setup Calo Assn for next couple of steps	
         std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());


         // Get Track Start point
         TVector3 TrackStart(trk->Start().X(),trk->Start().Y(),trk->Start().Z());
         TVector3 TrackDir(trk->StartDirection().X(),trk->StartDirection().Y(),trk->StartDirection().Z());


         t_TrackStart_Channel_Plane0.push_back(U_wire(TrackStart)); 
         t_TrackStart_Time_Plane0.push_back(tick(TrackStart)); 

         t_TrackStart_Channel_Plane1.push_back(V_wire(TrackStart)); 
         t_TrackStart_Time_Plane1.push_back(tick(TrackStart)); 

         t_TrackStart_Channel_Plane2.push_back(Y_wire(TrackStart)); 
         t_TrackStart_Time_Plane2.push_back(tick(TrackStart)); 

         t_TrackStart_X.push_back(TrackStart.X());
         t_TrackStart_Y.push_back(TrackStart.Y());
         t_TrackStart_Z.push_back(TrackStart.Z());

         t_TrackDir_X.push_back(TrackDir.X());
         t_TrackDir_Y.push_back(TrackDir.Y());
         t_TrackDir_Z.push_back(TrackDir.Z());


         // Get Calo Start Point

         int U = -1000;
         int V = -1000;
         int Y = -1000;
         int tick_Plane0 = -1000;                                        
         int tick_Plane1 = -1000;                                        
         int tick_Plane2 = -1000;                                        

         TVector3 CaloStart_Plane0(-1000,-1000,-1000);
         TVector3 CaloStart_Plane1(-1000,-1000,-1000);
         TVector3 CaloStart_Plane2(-1000,-1000,-1000);

         for(size_t i_plane=0;i_plane<caloFromTrack.size();i_plane++){

            auto thisPlaneCalo = caloFromTrack.at(i_plane);
            int planeno = thisPlaneCalo->PlaneID().Plane;

            // Get the last point in the calo object
            if(thisPlaneCalo->XYZ().size()){

               // Start of the track is the last point in the Calo point list
               TVector3 CaloStart(thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).X(),thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).Y(),thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).Z());


               if(planeno == 0){
                  U = U_wire(CaloStart);
                  tick_Plane0 = tick(CaloStart);
                  CaloStart_Plane0 = CaloStart;
               } 
               if(planeno == 1){
                  V = V_wire(CaloStart);
                  tick_Plane1 = tick(CaloStart);
                  CaloStart_Plane1 = CaloStart;
               }
               if(planeno == 2){
                  Y = Y_wire(CaloStart);
                  tick_Plane2 = tick(CaloStart);
                  CaloStart_Plane2 = CaloStart;
               } 
            }
         }              

         t_CaloStart_Channel_Plane0.push_back(U);
         t_CaloStart_Time_Plane0.push_back(tick_Plane0);
         t_CaloStart_X_Plane0.push_back(CaloStart_Plane0.X());
         t_CaloStart_Y_Plane0.push_back(CaloStart_Plane0.Y());
         t_CaloStart_Z_Plane0.push_back(CaloStart_Plane0.Z());

         t_CaloStart_Channel_Plane1.push_back(V);
         t_CaloStart_Time_Plane1.push_back(tick_Plane1);
         t_CaloStart_X_Plane1.push_back(CaloStart_Plane1.X());
         t_CaloStart_Y_Plane1.push_back(CaloStart_Plane1.Y());
         t_CaloStart_Z_Plane1.push_back(CaloStart_Plane1.Z());

         t_CaloStart_Channel_Plane2.push_back(Y);
         t_CaloStart_Time_Plane2.push_back(tick_Plane2);
         t_CaloStart_X_Plane2.push_back(CaloStart_Plane2.X());
         t_CaloStart_Y_Plane2.push_back(CaloStart_Plane2.Y());
         t_CaloStart_Z_Plane2.push_back(CaloStart_Plane2.Z());

      } //pfpTracks.size() == 1 

   } //end of PFP loop


   /////////////////////////////
   //   Obtain Wire Signals   //
   /////////////////////////////


   // Setup handles
   art::Handle<std::vector<recob::Wire>> wireHandle;
   std::vector<art::Ptr<recob::Wire>> wireVect;

   // Fill Wire vector
   if(e.getByLabel(fWireLabel,wireHandle)) art::fill_ptr_vector(wireVect,wireHandle);
   else
      std::cout << "Wire handle not setup" << std::endl;

   // Iterate through all of the wires, record signal at every tick with nonzero signal
   for(const art::Ptr<recob::Wire> &wire : wireVect){

      // Get regions of interest        
      unsigned int NROI = wire->SignalROI().n_ranges();
      for(size_t i_roi=0; i_roi<NROI; ++i_roi){

         // Region of tick space with nonzero activity
         recob::Wire::RegionsOfInterest_t::datarange_t const& range = wire->SignalROI().range(i_roi);

         // Iterate through the ticks in this ROI, record signal
         unsigned int thisTick = range.begin_index();

         while(thisTick < range.end_index()){

            if(wire->View() == 0){
               t_Wire_Channel_Plane0.push_back(wire->Channel());
               t_Wire_Tick_Plane0.push_back(thisTick);
               t_Wire_Signal_Plane0.push_back(wire->Signal().at(thisTick));
            }

            if(wire->View() == 1){
               t_Wire_Channel_Plane1.push_back(wire->Channel());
               t_Wire_Tick_Plane1.push_back(thisTick);
               t_Wire_Signal_Plane1.push_back(wire->Signal().at(thisTick));
            }

            if(wire->View() == 2){
               t_Wire_Channel_Plane2.push_back(wire->Channel());
               t_Wire_Tick_Plane2.push_back(thisTick);
               t_Wire_Signal_Plane2.push_back(wire->Signal().at(thisTick));
            }

            thisTick++;

         } // while(thisTick < range.end_index

      } // loop over ROI

   } // loop over wires


   // Fill tree! 
   t_WireTree->Fill();

}



//////////////////////////////////////////////////////////////////


void hyperon::WireTreeMaker::beginJob(){


   if(fDebug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //              Wire Tree		   //
   //////////////////////////////////////////


   t_WireTree=tfs->make<TTree>("WireTree","Wire Tree");

   t_WireTree->Branch("EventID",&t_EventID);
   t_WireTree->Branch("run",&t_run);
   t_WireTree->Branch("subrun",&t_subrun);
   t_WireTree->Branch("event",&t_event);

   t_WireTree->Branch("Wire_Channel_Plane0",&t_Wire_Channel_Plane0);
   t_WireTree->Branch("Wire_Tick_Plane0",&t_Wire_Tick_Plane0);
   t_WireTree->Branch("Wire_Signal_Plane0",&t_Wire_Signal_Plane0);

   t_WireTree->Branch("Wire_Channel_Plane1",&t_Wire_Channel_Plane1);
   t_WireTree->Branch("Wire_Tick_Plane1",&t_Wire_Tick_Plane1);
   t_WireTree->Branch("Wire_Signal_Plane1",&t_Wire_Signal_Plane1);

   t_WireTree->Branch("Wire_Channel_Plane2",&t_Wire_Channel_Plane2);
   t_WireTree->Branch("Wire_Tick_Plane2",&t_Wire_Tick_Plane2);
   t_WireTree->Branch("Wire_Signal_Plane2",&t_Wire_Signal_Plane2);

   t_WireTree->Branch("TrackStart_Channel_Plane0",&t_TrackStart_Channel_Plane0);
   t_WireTree->Branch("TrackStart_Time_Plane0",&t_TrackStart_Time_Plane0);

   t_WireTree->Branch("TrackStart_Channel_Plane1",&t_TrackStart_Channel_Plane1);
   t_WireTree->Branch("TrackStart_Time_Plane1",&t_TrackStart_Time_Plane1);

   t_WireTree->Branch("TrackStart_Channel_Plane2",&t_TrackStart_Channel_Plane2);
   t_WireTree->Branch("TrackStart_Time_Plane2",&t_TrackStart_Time_Plane2);

   t_WireTree->Branch("TrackStart_X",&t_TrackStart_X);
   t_WireTree->Branch("TrackStart_Y",&t_TrackStart_Y);
   t_WireTree->Branch("TrackStart_Z",&t_TrackStart_Z);

   t_WireTree->Branch("TrackDir_X",&t_TrackDir_X);
   t_WireTree->Branch("TrackDir_Y",&t_TrackDir_Y);
   t_WireTree->Branch("TrackDir_Z",&t_TrackDir_Z);


   t_WireTree->Branch("CaloStart_Channel_Plane0",&t_CaloStart_Channel_Plane0);
   t_WireTree->Branch("CaloStart_Time_Plane0",&t_CaloStart_Time_Plane0);        
   t_WireTree->Branch("CaloStart_X_Plane0",&t_CaloStart_X_Plane0);
   t_WireTree->Branch("CaloStart_Y_Plane0",&t_CaloStart_Y_Plane0);
   t_WireTree->Branch("CaloStart_Z_Plane0",&t_CaloStart_Z_Plane0);

   t_WireTree->Branch("CaloStart_Channel_Plane1",&t_CaloStart_Channel_Plane1);
   t_WireTree->Branch("CaloStart_Time_Plane1",&t_CaloStart_Time_Plane1);
   t_WireTree->Branch("CaloStart_X_Plane1",&t_CaloStart_X_Plane1);
   t_WireTree->Branch("CaloStart_Y_Plane1",&t_CaloStart_Y_Plane1);
   t_WireTree->Branch("CaloStart_Z_Plane1",&t_CaloStart_Z_Plane1);

   t_WireTree->Branch("CaloStart_Channel_Plane2",&t_CaloStart_Channel_Plane2);
   t_WireTree->Branch("CaloStart_Time_Plane2",&t_CaloStart_Time_Plane2);
   t_WireTree->Branch("CaloStart_X_Plane2",&t_CaloStart_X_Plane2);
   t_WireTree->Branch("CaloStart_Y_Plane2",&t_CaloStart_Y_Plane2);
   t_WireTree->Branch("CaloStart_Z_Plane2",&t_CaloStart_Z_Plane2);


   if(fDebug) std::cout << "Finished begin job" << std::endl;


}



void hyperon::WireTreeMaker::endJob()
{

}



void hyperon::WireTreeMaker::beginSubRun(const art::SubRun& sr)
{


}


void hyperon::WireTreeMaker::endSubRun(const art::SubRun& sr){}


DEFINE_ART_MODULE(hyperon::WireTreeMaker)
