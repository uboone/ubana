////////////////////////////////////////////////////////////////////////
// Class:       HyperonSelection
// Plugin Type: analyzer (art v3_03_01)
// File:        HyperonSelection_module.cc
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

#include "larcoreobj/SummaryData/POTSummary.h"


#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "TTree.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include <vector>
#include <string>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

//track momentum calculators
#include "larreco/RecoAlg/TrackMomentumCalculator.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

//root includes
#include "TVector3.h"
#include "TLorentzVector.h"

//local includes

//objects and helpers
#include "ubana/HyperonProduction/util/SimParticle.h"
#include "ubana/HyperonProduction/util/RecoParticle.h"
#include "ubana/HyperonProduction/util/Helpers.h"

//algorithms
#include "ubana/HyperonProduction/Alg/ParticleTypes.h"
#include "ubana/HyperonProduction/Alg/FV.h"
#include "ubana/HyperonProduction/Alg/MeandEdX.h"
#include "ubana/HyperonProduction/Alg/ThreePlaneMeandEdX.h"
#include "ubana/HyperonProduction/Alg/ConnectednessHelper.h"

//PID components
#include "ubana/HyperonProduction/Alg/LLR_PID.h"
#include "ubana/HyperonProduction/Alg/LLRPID_proton_muon_lookup.h"


namespace hyperon {
   class HyperonSelection;
}


class hyperon::HyperonSelection : public art::EDAnalyzer {
   public:
      explicit HyperonSelection(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HyperonSelection(HyperonSelection const&) = delete;
      HyperonSelection(HyperonSelection&&) = delete;
      HyperonSelection& operator=(HyperonSelection const&) = delete;
      HyperonSelection& operator=(HyperonSelection&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;



      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void PrintInfo();
      void FinishEvent();

      //lookup the origin type (primary, decay, other) of mc particle by
      //its TrackId
      int getOrigin(int idnum);

      //check if event contains a reco'd muon, proton and pion from Lambda decay
      //records their positions in track vector if they exist
      void StoreTrackTruth();


      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);


   private:

      //takes information retrieved from event in analyze and 
      //applies selection criteria

      bool PerformSelection();


      //basic event info
      unsigned int fEventID;
      //run/subrun/event information
      int run,subrun,event;

      int fileID;

      //output trees
      TTree * fOutputTree; //info for each event

      TTree * fMetaTree; //metadata


      ///////////////////////////
      //   Truth level info    //
      ///////////////////////////

      //generator truth info

      int fNMCTruths=0;	

      std::string fMode; //interaction mode
      std::string fCCNC; //charged current/neutral current

      bool fInActiveTPC=false; //is event in active TPC

      bool fIsHyperon=false; //true if event contains a hyperon
      bool fIsSigmaZero=false; //true if event hyperon was a sigma zero		
      bool fIsLambda=false;
      bool fIsLambdaCharged=false; //true if event contains Lambda decaying into p + pi-

      bool fIsAssociatedHyperon=false;

      bool fIsSignal=false; //true if numu event in fiducial vol producing Lambda decaying to p + pi- -->Signal to search for	
      bool fGoodReco=false; //true if is signal event with both decay products truth matching to reconstructed tracks

      double fWeight=1.0;

      //neutrino
      std::vector<SimParticle> fNeutrino;

      //lepton produced at primary vtx		
      std::vector<SimParticle> fLepton;

      //hyperon produced at primary vtx
      std::vector<SimParticle> fHyperon;

      //nucleons produced at primary vtx
      std::vector<SimParticle> fPrimaryNucleon;

      //pions produced at primary vtx
      std::vector<SimParticle> fPrimaryPion;

      //kaons produced at primary vtx
      std::vector<SimParticle> fPrimaryKaon; 

      //vertex information
      TVector3 fTruePrimaryVertex;

      //g4 truth info
      TVector3 fDecayVertex;

      std::vector<SimParticle> fDecay; //hyperon decay products

      std::vector<SimParticle> fSigmaZeroDecayPhoton; //photon produced by Sigma zero decay
      std::vector<SimParticle> fSigmaZeroDecayLambda;	//lambda produced by Sigma zero decay

      std::vector<SimParticle> fKaonDecay; //hyperon decay products

      double fDecayOpeningAngle; //opening angle between hyperon decay products
      double fLeptonPionAngle; //openining angle between lepton and pion
      double fLeptonNucleonAngle; //opening angle between lepton and nucleon

      int fDecayHyperonID;
      int fDecaySigmaID;

      //data storage (should not be written to trees)

      //used by G4 to track particles
      std::vector<int>daughter_IDs; //ids of semistable hyperon decay products
      std::vector<int>Sigma0_daughter_IDs; //ids of sigma0 decay products
      std::vector<int>primary_IDs; //ids of particles produced at primary vertex
      std::vector<int>Kaon_daughter_IDs; //ids of sigma0 decay products

      //create map between particles and their ID's
      std::map<int,art::Ptr<simb::MCParticle>> partByID;
      std::pair<int,art::Ptr<simb::MCParticle>>  part_and_id;

      int mode; //interaction mode for event generator
      int ccnc;


      ///////////////////////////
      //   Reco level info    //
      ///////////////////////////

      bool fSelectedEvent = false; //true if event passes some selection criteria	

      TVector3 fRecoPrimaryVertex;

      int fNPrimaryDaughters; //num of primary daughters
      int fNPrimaryTrackDaughters; //num of track like primary daughters
      int fNPrimaryShowerDaughters; //num of shower like primary daughters


      //Primary daughters
      std::vector<RecoParticle> fTrackPrimaryDaughters;
      std::vector<RecoParticle> fShowerPrimaryDaughters;
      int fMuonIndex=-1; //index of muon candidate in fTrackPrimaryDaughters , -1 if no muon candidate found

      ////////////////////////////
      //   Connectedness test   //
      ////////////////////////////

      std::vector<std::vector<int>> Conn_SeedIndexes_Plane0;
      std::vector<std::vector<int>> Conn_OutputIndexes_Plane0;
      std::vector<std::vector<int>> Conn_OutputSizes_Plane0;
      std::vector<std::vector<int>> Conn_SeedChannels_Plane0;
      std::vector<std::vector<int>> Conn_SeedTicks_Plane0;

      std::vector<std::vector<int>> Conn_SeedIndexes_Plane1;
      std::vector<std::vector<int>> Conn_OutputIndexes_Plane1;
      std::vector<std::vector<int>> Conn_OutputSizes_Plane1;
      std::vector<std::vector<int>> Conn_SeedChannels_Plane1;
      std::vector<std::vector<int>> Conn_SeedTicks_Plane1;

      std::vector<std::vector<int>> Conn_SeedIndexes_Plane2;
      std::vector<std::vector<int>> Conn_OutputIndexes_Plane2;
      std::vector<std::vector<int>> Conn_OutputSizes_Plane2;
      std::vector<std::vector<int>> Conn_SeedChannels_Plane2;
      std::vector<std::vector<int>> Conn_SeedTicks_Plane2;

      ///////////////////////////
      //  Truth Matching info  //
      ///////////////////////////

      //indices in the reco track vector of the true muon
      //and proton and pion from hyperon decay if the exist

      //will have values of -1 if they do not exist

      int fTrueMuonIndex=-1;
      int fTrueDecayProtonIndex=-1;
      int fTrueDecayPionIndex=-1;

      /////////////////////////
      // Metadata for sample //
      /////////////////////////

      //truth level metadata

      int fNEvents; //total events in sample
      int fNEventsInActiveVol; //total events in active volume (at truth level)

      int fNChargedCurrent; //number of cc events
      int fNNeutralCurrent; //number of nc events
      int fNnuMu; //number of numu events
      int fNnue; //number of nue events
      int fNnuMuBar; //number of numubar events
      int fNnueBar; //number of nuebar events

      int fNHyperons;

      int fNSignal; //number of signal events

      int fNGoodReco; //number of signal events with both pion and proton reco'd

      //reco level metadata

      int fNSelectedEvents; //total events passing selection

      int fNSelectedHyperons; //total true hyperon events passing selection

      int fNSelectedSignal; //number of signal events passing selection

      int fNSelectedGoodReco; //number of signal events passing selection

      //selection performance metrics

      double fHyperonEfficiency; //hyperon selection efficiency
      double fHyperonPurity;  //hyperon selection purity
      double fHyperonTruePurity; //hyperon selection efficiency x purity
      double	fHyperonEfficiencyTimesPurity; //hyperon selection efficiency x purity
      double	fHyperonEfficiencyTimesTruePurity; //hyperon selection efficiency x true purity

      double	fSignalEfficiency=0; //hyperon selection efficiency
      double	fSignalPurity=0;  //hyperon selection purity
      double	fSignalTruePurity=0; //hyperon purity after converting from enriched sample to real sample
      double	fSignalEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
      double	fSignalEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

      double	fGoodRecoEfficiency=0; //hyperon selection efficiency
      double	fGoodRecoPurity=0;  //hyperon selection purity
      double	fGoodRecoTruePurity=0; //hyperon purity after converting from enriched sample to real sample
      double	fGoodRecoEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
      double	fGoodRecoEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

      double fPOT=0; //total POT of the sample

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////


      bool fIsData;

      //producer module labels

      std::string fGenieGenModuleLabel;
      std::string fGeantModuleLabel;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fVertexLabel;
      std::string fPFParticleLabel;
      std::string fPIDLabel;
      std::string fCaloLabel;
      std::string fHitLabel;
      std::string fSpacePointLabel;
      std::string fTrackHitAssnLabel;
      std::string fHitTruthAssnLabel;
      std::string fMetadataLabel;
      std::string fShowerHitAssnLabel;
      std::string fPFPSpacePointAssnLabel;
      std::string fWireLabel;


      //POT module label
      std::string fPOTSummaryLabel;

      //misc
      bool fPrint;
      int fPrintPdg=-1; //-1 means print everything
      bool fDebug=false;

      double fHyperonEnrichment = 1.0; //hyperon production enrichment used

      ///////////////////////////////////
      // Hyperon selection parameters  //
      ///////////////////////////////////

      //Decay Proton/Pion Momentum thresholds (in GeV)
      double fDecayProtonThreshold;
      double fDecayPionThreshold;	

      int fMinDaughters;//minimum number of daughters
      int fMaxPrimaryShowerDaughters;
      int fMinDaughterTracks; //minimum number of track like primary daughters
      double fMuonPIDScore;
      double fMinimumMuonLength;
      double fSecondaryTrackLengthCut;
      double fTertiaryTrackLengthCut;

      double fTrackScore; //minimum track/shower score required for pfp to be classified as track

      ///////////////////////
      //      Objects      //
      ///////////////////////

      //LLR PID components
      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

      ConnectednessHelper Conn_Helper;


};


////////////////////////////////////////////////////
// Decides if an event passes selection criteria  //
////////////////////////////////////////////////////

bool hyperon::HyperonSelection::PerformSelection(){

   if(fDebug) std::cout << "Applying Selection" << std::endl;

   //Apply fiducial volume cut
   if(!inActiveTPC(fRecoPrimaryVertex)) return false;
   if(fNPrimaryDaughters < fMinDaughters) return false;

   if(fNPrimaryShowerDaughters > fMaxPrimaryShowerDaughters || fNPrimaryTrackDaughters < fMinDaughterTracks){

      return false;

   }

   return true;

}


////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////



hyperon::HyperonSelection::HyperonSelection(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   Conn_Helper(p.get<bool>("DrawConnectedness",false))   // ,
   // More initializers here.
{

   fIsData = p.get<bool>("IsData");

   //module labels
   fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
   fPFParticleLabel = p.get<std::string>("PFParticleLabel");
   fGeantModuleLabel = p.get<std::string>("GeantLabel");
   fVertexLabel = p.get<std::string>("VertexLabel");
   fTrackLabel = p.get<std::string>("TrackLabel");
   fShowerLabel = p.get<std::string>("ShowerLabel");
   fPIDLabel = p.get<std::string>("PIDLabel");
   fCaloLabel = p.get<std::string>("CaloLabel");
   fHitLabel  = p.get<std::string>("HitLabel");
   fSpacePointLabel = p.get<std::string>("SpacePointLabel");
   fTrackHitAssnLabel = p.get<std::string>("TrackHitAssnLabel");
   fHitTruthAssnLabel = p.get<std::string>("HitTruthAssnLabel");
   fShowerHitAssnLabel = p.get<std::string>("ShowerHitAssnLabel");
   fPFPSpacePointAssnLabel = p.get<std::string>("PFPSpacePointAssnLabel");
   fMetadataLabel = p.get<std::string>("MetadataLabel");
   fWireLabel = p.get<std::string>("WireLabel");

   //POT module label
   fPOTSummaryLabel = p.get<std::string>("POTSummaryLabel");

   //Decay Proton/Pion Momentum thresholds (in GeV)
   fDecayProtonThreshold = p.get<double>("DecayProtonThreshold",0.2);	
   fDecayPionThreshold = p.get<double>("DecayPionThreshold",0.05);	


   //selection parameters
   fMinDaughters = p.get<int>("MinDaughters",0);	
   fMinDaughterTracks = p.get<int>("MinDaughterTracks",0);
   fMaxPrimaryShowerDaughters = p.get<int>("MaxDaughterShowers",1000);
   fMuonPIDScore = p.get<double>("MuonPIDScore",0.6);
   fMinimumMuonLength = p.get<double>("MinimumMuonLength",0.0);
   fSecondaryTrackLengthCut = p.get<double>("SecondaryTrackLengthCut",65);
   fTertiaryTrackLengthCut = p.get<double>("TertiaryTrackLengthCut",35);

   //track/shower classification threshold
   fTrackScore = p.get<double>("TrackScore",0.5);

   //misc parameters
   fPrint = p.get<bool>("Print",false);
   fPrintPdg = p.get<int>("PrintPdg",-1);
   fDebug = p.get<bool>("Debug",false);

   //Hyperon enrichment used
   fHyperonEnrichment = p.get<double>("HyperonEnrichment",1.0);

   // set dedx pdf parameters for LLR PID
   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

}

void hyperon::HyperonSelection::analyze(art::Event const& e)
{

   if(fDebug) std::cout << "New Event" << std::endl;

   //begin by resetting everything

   /////////////
   // General //
   /////////////

   fNMCTruths=0;

   fInActiveTPC=false;
   fIsHyperon=false;
   fIsSigmaZero=false;
   fIsLambda=false;
   fIsLambdaCharged=false;
   fIsAssociatedHyperon=false;
   fIsSignal=false;	
   fGoodReco=false;

   fWeight=1.0;

   //default mode for real data
   fMode = "NONE";

   /////////////
   // G4 Info //
   /////////////

   //neutrino that interacted
   fNeutrino.clear();

   //lepton produced in interaction
   fLepton.clear();

   //hyperon produced
   fHyperon.clear();

   //nucleons, pions and kaons produced at primary vtx
   fPrimaryNucleon.clear();
   fPrimaryPion.clear();
   fPrimaryKaon.clear();

   //vertex information
   fTruePrimaryVertex.SetXYZ(-1000,-1000,-1000);

   //hyperon decay products
   fDecayVertex.SetXYZ(-1000,-1000,-1000);
   fDecay.clear();

   //sigma zero decay products
   fSigmaZeroDecayPhoton.clear();
   fSigmaZeroDecayLambda.clear();

   fDecayOpeningAngle=-1; //opening angle between hyperon decay products   
   fLeptonPionAngle=-1; //openining angle between lepton and pion
   fLeptonNucleonAngle=-1; //opening angle between lepton and nucleon

   fKaonDecay.clear();

   ///////////////
   // Reco Info //
   ///////////////

   fSelectedEvent = false; //true if event passes some selection criteria	
   fNPrimaryDaughters = 0; //number of primary daughters

   fNPrimaryTrackDaughters=0; //num of track like primary daughters
   fNPrimaryShowerDaughters=0; //num of shower like primary daughters

   fRecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

   fTrackPrimaryDaughters.clear();
   fShowerPrimaryDaughters.clear();
   fMuonIndex=-1;

   //Truth Matched info - indices of muon, proton and pion from decay
   //in fTrackPrimaryDaughters (if they exist)		
   fTrueMuonIndex=-1;	
   fTrueDecayProtonIndex=-1;
   fTrueDecayPionIndex=-1;

   ////////////////////////////
   //   Connectedness test   //
   ////////////////////////////

   Conn_SeedIndexes_Plane0.clear();
   Conn_OutputIndexes_Plane0.clear();
   Conn_OutputSizes_Plane0.clear();
   Conn_SeedChannels_Plane0.clear();
   Conn_SeedTicks_Plane0.clear();

   Conn_SeedIndexes_Plane1.clear();
   Conn_OutputIndexes_Plane1.clear();
   Conn_OutputSizes_Plane1.clear();
   Conn_SeedChannels_Plane1.clear();
   Conn_SeedTicks_Plane1.clear();

   Conn_SeedIndexes_Plane2.clear();
   Conn_OutputIndexes_Plane2.clear();
   Conn_OutputSizes_Plane2.clear();
   Conn_SeedChannels_Plane2.clear();
   Conn_SeedTicks_Plane2.clear();

   /////////////////////////
   // Event ID information /
   /////////////////////////

   fEventID = e.id().event();
   run = e.run();
   subrun = e.subRun();
   event = e.event();

   //////////////////////////////
   // GET EVENT GENERATOR INFO //
   //////////////////////////////

   if(!fIsData){

      if(fDebug) std::cout << "Getting Event Generator Info" << std::endl;

      art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
      std::vector<art::Ptr<simb::MCTruth> > mcTrVect;

      if(e.getByLabel(fGenieGenModuleLabel,mctruthListHandle)) art::fill_ptr_vector(mcTrVect,mctruthListHandle);    

      if(!mcTrVect.size()) {
         FinishEvent();
         return;
      }

      fNMCTruths = mcTrVect.size();

      art::Ptr<simb::MCTruth> MCtruth = mcTrVect.at(0);

      //get the neutrino info

      simb::MCNeutrino Nu = MCtruth->GetNeutrino();
      mode = Nu.Mode();
      ccnc = Nu.CCNC();

      if(ccnc == 0) fCCNC = "CC";
      else fCCNC = "NC";

      // NOTE: Hyperon events produced in GENIE use mode == 0
      // NuWro uses 1095 
      if(mode == 0) fMode = "QEL";
      else if(mode == 1) fMode = "RES";
      else if(mode == 2) fMode = "DIS";
      else if(mode == 3) fMode = "COH";
      else if(mode == 5) fMode = "ElectronScattering";
      else if(mode == 10) fMode = "MEC";
      else if(mode == 11) fMode = "Diffractive";
      else if(mode == 1095) { fMode = "HYP";  fWeight *= 1.0/fHyperonEnrichment; }
      else fMode = "Other";	

      for(int k_particles=0;k_particles<MCtruth->NParticles();k_particles++){

         simb::MCParticle Part = MCtruth->GetParticle(k_particles);

         // Get list of particles from true PV, if lepton or neutrino set PV
         if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1){ 
            fTruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());
            fInActiveTPC=inActiveTPC(fTruePrimaryVertex);
         }//if lepton

         // Record info about the neutrino
         if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
            SimParticle P;
            P.SetKinematics(Part.Momentum(),Part.Mass());				
            P.SetPositions(Part.Position(),Part.EndPosition());	
            P.PDG = Part.PdgCode();
            P.Origin = 0;
            fNeutrino.push_back( P );
         }

      }//k_particles (particles associated with MC truth)

   } // !fIsData

   ////////////////////////////////////////////
   // Get Geant Information
   ///////////////////////////////////////////

   if(!fIsData){

      if(fDebug) std::cout << "Getting Event G4 Info" << std::endl;

      //get list of geant4 particles
      art::Handle<std::vector<simb::MCParticle>>g4particleHandle;
      std::vector< art::Ptr<simb::MCParticle>>g4particleVect;
      g4particleVect.clear();

      if(e.getByLabel(fGeantModuleLabel,g4particleHandle)) art::fill_ptr_vector(g4particleVect,g4particleHandle);
      else
         std::cout << "Geant Particles Missing!" << std::endl;


      //id numbers of hyperon daughters
      daughter_IDs.clear(); 

      //id numbers of particles produced at primary vertex
      primary_IDs.clear();
      //id nubmers of Sigma Zero decay products
      Sigma0_daughter_IDs.clear();

      Kaon_daughter_IDs.clear();

      //map between particle ID numbers (g4p->TrackId()) and pointers to simb::MCParticle
      partByID.clear();

      for(const art::Ptr<simb::MCParticle> &g4p : g4particleVect){

         //if mother == 0 particle is produced at primary vertex
         if(g4p->Mother() == 0){ 

            //store vector of id's of primary particles	
            primary_IDs.push_back(g4p->TrackId());

            //if particle is a hyperon, record start position, pdg etc
            if(isHyperon(g4p->PdgCode())){

               //if it is a sigma 0, add its daughters to semistable decay products
               if(g4p->PdgCode() == 3212){

                  fIsSigmaZero = true;

                  if(g4p->EndProcess() == "Decay"){
                     for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

                        Sigma0_daughter_IDs.push_back( g4p->Daughter(i_d) );

                     }
                  }

               }
               //if not sigma 0, add daughters to list of decay products, record end point as decay vertex and distance travelled
               else
               {

                  if(g4p->PdgCode() == 3122) fIsLambda=true;

                  if(g4p->EndProcess() == "Decay"){

                     for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

                        daughter_IDs.push_back( g4p->Daughter(i_d) );

                     }

                     //if not Sigma0, record decay vertex

                     fDecayVertex.SetXYZ( g4p->EndPosition().X() , g4p->EndPosition().Y() , g4p->EndPosition().Z() );

                  }
               }		

            }
            else if(isKaon(g4p->PdgCode()) && g4p->EndProcess() == "Decay"){

               for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

                  Kaon_daughter_IDs.push_back( g4p->Daughter(i_d) );

               }


            }


         }

         part_and_id = std::make_pair(g4p->TrackId() , g4p);

         partByID.insert( part_and_id );

      }


      //you now have a map between track ids and particles
      //and a vector of id numbers of particles produced at primary vtx, 
      //id numbers of daughters of hyperon if it is Lambda/SigmaM/SigmaP
      //id numbers of daughters of Sigma0 if there is one


      //go through each of these lists, collect data as needed

      for(size_t i_p=0;i_p<primary_IDs.size();i_p++){

         //geant does not always keep everything it simulates, make sure particle is in list of IDs (crashes otherwise!)

         if(partByID.find(primary_IDs[i_p]) == partByID.end()) continue;

         art::Ptr<simb::MCParticle> part = partByID[primary_IDs.at(i_p)];

         if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

         SimParticle P = MakeSimParticle(*part);
         P.Origin = getOrigin(part->TrackId());

         //hyperon produced at primary vertex
         if( isHyperon(part->PdgCode()) ) {

            fHyperon.push_back( P );
            fIsHyperon = true;	
         }

         //lepton produced at primary vertex
         if( isLepton(part->PdgCode()) || isNeutrino(part->PdgCode()) ) fLepton.push_back( P );

         //nucleon produced at primary vertex
         if( isNucleon(part->PdgCode()) ) fPrimaryNucleon.push_back( P );

         //pion produced at primary vertex
         if( isPion(part->PdgCode()) ) fPrimaryPion.push_back( P );

         //kaon from primary vertex
         if( isKaon(part->PdgCode()) ) fPrimaryKaon.push_back( P );


      }





      //check all decay products are actually produced at decay vertex - sometimes you get some electrons thrown in

      if(fHyperon.size() == 1){

         //if using GENIE as evgen, hyperon events get labelled as QEL - change this to HYP
         if( fMode == "QEL" ) fMode = "HYP";

         std::vector<int> daughter_IDs_tmp;

         for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){

            //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
            if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

            art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

            if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

            if(part->Position().X() == fHyperon.at(0).EndX && part->Position().Y() == fHyperon.at(0).EndY && part->Position().Z() == fHyperon.at(0).EndZ){

               daughter_IDs_tmp.push_back(daughter_IDs.at(i_d));

            }
         }

         daughter_IDs = daughter_IDs_tmp;

      }


      if(fPrimaryKaon.size() == 1){

         std::vector<int> Kaon_daughter_IDs_tmp;

         for(size_t i_d=0;i_d<Kaon_daughter_IDs.size();i_d++){

            //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
            if(partByID.find(Kaon_daughter_IDs[i_d]) == partByID.end()) continue;

            art::Ptr<simb::MCParticle> part = partByID[Kaon_daughter_IDs[i_d]];

            if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

            if(part->Position().X() == fPrimaryKaon.at(0).EndX && part->Position().Y() == fPrimaryKaon.at(0).EndY && part->Position().Z() == fPrimaryKaon.at(0).EndZ){

               Kaon_daughter_IDs_tmp.push_back(Kaon_daughter_IDs.at(i_d));

            }
         }

         Kaon_daughter_IDs = Kaon_daughter_IDs_tmp;

      }




      //now go through list of hyperon daughters, get info about the decay
      for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){
         //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
         if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

         art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

         if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

         SimParticle Decay = MakeSimParticle(*part);
         Decay.Origin = getOrigin(part->TrackId());
         fDecay.push_back( Decay );



         if(isPion(part->PdgCode()) && fLepton.size() == 1 ) {

            SimParticle Lepton = fLepton.at(0);
            fLeptonPionAngle = (180/3.1416)*TMath::ACos( (Decay.Px*Lepton.Px + Decay.Py*Lepton.Py + Decay.Pz*Lepton.Pz )/(Lepton.ModMomentum*Decay.ModMomentum) );

         }

         if(isNucleon(part->PdgCode()) && fLepton.size() == 1){


            SimParticle Lepton = fLepton.at(0);
            fLeptonNucleonAngle = (180/3.1416)*TMath::ACos( (Decay.Px*Lepton.Px + Decay.Py*Lepton.Py + Decay.Pz*Lepton.Pz )/(Lepton.ModMomentum*Decay.ModMomentum) );

         }



      }




      if(fDecay.size() == 2){

         fDecayOpeningAngle = (180/3.1416)*TMath::ACos( (fDecay.at(0).Px*fDecay.at(1).Px + fDecay.at(0).Py*fDecay.at(1).Py + fDecay.at(0).Pz*fDecay.at(1).Pz)/(fDecay.at(0).ModMomentum*fDecay.at(1).ModMomentum));


         //if hyperon was a lambda, check the pdg codes of its decay products
         if(fIsLambda && ( (fDecay.at(0).PDG == 2212 && fDecay.at(1).PDG == -211) || (fDecay.at(1).PDG == 2212 && fDecay.at(0).PDG == -211) ) ) fIsLambdaCharged=true;


      }


      //now deal with sigma zeros

      if(fIsSigmaZero){
         for(size_t i_d=0; i_d<Sigma0_daughter_IDs.size();i_d++){

            //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
            if(partByID.find(Sigma0_daughter_IDs[i_d]) == partByID.end()) continue;

            art::Ptr<simb::MCParticle> part = partByID[Sigma0_daughter_IDs[i_d]];

            if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

            SimParticle P = MakeSimParticle(*part);			
            P.Origin = getOrigin(part->TrackId());

            //one particle should be a Lambda, one should be a photon	

            //if is photon
            if(part->PdgCode() == 22){

               fSigmaZeroDecayPhoton.push_back( P );

               //add its ID to list of primary ID's - produced at primary vtx so from reco point of view makes sense to put it in
               //list of primary IDs

               primary_IDs.push_back(part->TrackId());

            }
            //Lambda produced from sigma zero decay
            else if(part->PdgCode() == 3122){

               //now search for Lambda decay products, find its decay vertex etc			

               fSigmaZeroDecayLambda.push_back( P );

               if(part->EndProcess() == "Decay"){		

                  //get decay vertex of resulting lambda and its decay products

                  fDecayVertex.SetXYZ( part->EndPosition().X() , part->EndPosition().Y() , part->EndPosition().Z() );


                  for(int i_d2=0;i_d2<part->NumberDaughters();i_d2++){

                     //search the map for the daughter IDs
                     if(partByID.find(part->Daughter(i_d2)) == partByID.end()) continue;

                     art::Ptr<simb::MCParticle> part2 = partByID[part->Daughter(i_d2)];

                     if(part2->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

                     SimParticle P2 = MakeSimParticle(*part2);
                     P2.Origin = getOrigin(part2->TrackId());

                     fDecay.push_back( P2 );

                     daughter_IDs.push_back( part2->TrackId() );

                  }

               }



            }

            else std::cout << "Unrecognized Sigma0 daughter: " << part->PdgCode() << std::endl; 

         }

      }


      // Add Kaon daughters

      //now go through list of hyperon daughters, get info about the decay
      for(size_t i_d=0;i_d<Kaon_daughter_IDs.size();i_d++){
         //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
         if(partByID.find(Kaon_daughter_IDs[i_d]) == partByID.end()) continue;

         art::Ptr<simb::MCParticle> part = partByID[Kaon_daughter_IDs[i_d]];

         if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

         SimParticle KaonDecay = MakeSimParticle(*part);
         KaonDecay.Origin = getOrigin(part->TrackId());
         fKaonDecay.push_back( KaonDecay );

      }


   }//if !IsData


   ///////////////////////////////////////////////////////////////////////////
   //Get Reconstructed Info
   //////////////////////////////////////////////////////////////////////////

   //genie uses QEL for hyperon events, NuWro uses HYP
   if(fNeutrino.size() == 1 && fInActiveTPC && fIsLambdaCharged && fNeutrino.at(0).PDG == -14 && ( fMode == "QEL" || fMode == "HYP")){ 

      //add kinematic thresholds
      if(fDecay.at(0).PDG == 2212 && fDecay.at(0).ModMomentum > fDecayProtonThreshold && fDecay.at(1).PDG == -211 && fDecay.at(1).ModMomentum > fDecayPionThreshold) fIsSignal = true;
      else if(fDecay.at(1).PDG == 2212 && fDecay.at(1).ModMomentum > fDecayProtonThreshold && fDecay.at(0).PDG == -211 && fDecay.at(0).ModMomentum > fDecayPionThreshold) fIsSignal = true;
      else fIsSignal = false;

   } 

   if(fDebug) std::cout << "Getting Reco'd Particles" << std::endl;

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   //calculate momentum of tracks under various hypotheses
   trkf::TrackMomentumCalculator trkm{0};

   //setup handles
   art::Handle<std::vector<recob::PFParticle>>pfparticleHandle; //PFParticles
   art::Handle<std::vector<recob::Track>  > trackHandle; //reconstructed tracks
   art::Handle<std::vector<recob::Shower> > showerHandle; //reconstructed showers
   art::Handle<std::vector<anab::ParticleID>>pidHandle;	
   art::Handle<std::vector<anab::Calorimetry>>caloHandle;
   art::Handle<std::vector<recob::Hit> > hitHandle;
   art::Handle<larpandoraobj::PFParticleMetadata> particleMetadataHandle;
   art::Handle<std::vector<recob::SpacePoint>> spacepointHandle;


   std::vector <art::Ptr <recob::Track>  > trackVect; //vector of tracks
   std::vector <art::Ptr <recob::Shower> > showerVect;//vector of showers
   std::vector< art::Ptr<recob::PFParticle>>pfparticleVect; //vector of PFparticles
   std::vector < art::Ptr<anab::ParticleID>> pidVect; //pids
   std::vector<art::Ptr<anab::Calorimetry>>caloVect; //calorimetry
   std::vector< art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVect; //metadata
   std::vector <art::Ptr <recob::Hit>  > hitVect; //hits
   std::vector<art::Ptr<recob::SpacePoint>> spacepointVect; //PFP spacepoints

   //Fill PFP vector

   if(e.getByLabel(fPFParticleLabel,pfparticleHandle)){
      art::fill_ptr_vector(pfparticleVect,pfparticleHandle);
   }

   if(!pfparticleVect.size()) {
      FinishEvent();
      return;
   }



   //fill track vector
   if(e.getByLabel(fTrackLabel,trackHandle)) art::fill_ptr_vector(trackVect,trackHandle);
   else
      std::cout << "Track handle not setup" << std::endl;

   //fill shower vector
   if(e.getByLabel(fShowerLabel,showerHandle)) art::fill_ptr_vector(showerVect,showerHandle);
   else
      std::cout << "Shower handle not setup" << std::endl;



   e.getByLabel(fHitLabel,hitHandle);
   e.getByLabel(fSpacePointLabel,spacepointHandle);


   //fill hit vector
   if(e.getByLabel(fHitLabel,hitHandle)) art::fill_ptr_vector(hitVect,hitHandle);
   else
      std::cout << "Hit handle not setup" << std::endl;


   //vertices, tracks and showers assoc with PFPs
   art::FindManyP<recob::Vertex> vertexAssoc(pfparticleVect,e,fVertexLabel);
   art::FindManyP<recob::Track> trackAssoc(pfparticleVect,e,fTrackLabel);
   art::FindManyP<recob::Shower> showerAssoc(pfparticleVect,e,fShowerLabel);
   art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssoc(pfparticleVect,e,fMetadataLabel);

   art::FindManyP<anab::ParticleID> PIDAssoc(trackVect,e,fPIDLabel);
   art::FindManyP<anab::Calorimetry> caloTrackAssoc(trackVect,e,fCaloLabel);

   //get hits assoc with tracks
   art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackHitAssnLabel);
   //get hits assoc with showers
   art::FindManyP<recob::Hit> showerHitAssoc(showerHandle,e,fShowerHitAssnLabel);

   //backtracker
   art::FindMany< simb::MCParticle , anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle,e,fHitTruthAssnLabel);

   //spacepoints assoc with PFPs
   art::FindManyP<recob::SpacePoint> pfpSpacePointAssoc(pfparticleHandle,e,fPFPSpacePointAssnLabel);

   size_t neutrinoID = 99999;


   //go through the list of pandora PFP's
   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){

      std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());	

      //get reconstructed neutrino
      if((pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12 ))){

         neutrinoID = pfp->Self();

         fNPrimaryDaughters = pfp->NumDaughters();

         //get the reconstructed primary vertex
         for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

            //correct for space charge

            geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };                
            geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

            //no SC correction
            //fRecoPrimaryVertex.SetXYZ( vtx->position().X() , vtx->position().Y() , vtx->position().Z() );

            //w SC correction - forward
            fRecoPrimaryVertex.SetXYZ( vtx->position().X() + sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z() );



         }

      }

   }


   std::vector<TVector3> TrackStarts;


   //go through rest of particles, get lots of useful info!
   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){

      RecoParticle ThisPrimaryDaughter;

      //get data from every PFP, not just neutrino daughters
      if(pfp->Parent() != neutrinoID) continue;

      std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
      std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());

      std::vector< art::Ptr<recob::Shower> > pfpShowers = showerAssoc.at(pfp.key());
      std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> >pfpMeta = metadataAssoc.at(pfp.key());

      std::vector< art::Ptr<recob::SpacePoint> > pfpSpacePoints = pfpSpacePointAssoc.at(pfp.key());

      ThisPrimaryDaughter.PDG = pfp->PdgCode();

      for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){

         const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());

         if (!pfParticlePropertiesMap.empty()){
            for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){


               if(it->first == "TrackScore"){
                  ThisPrimaryDaughter.TrackShowerScore = it->second;
               }
            }
         }	
      }

      //////////////////////////////////////////////
      //            Get track info                //
      //////////////////////////////////////////////


      //if pfp is a track like object
      if(pfpTracks.size() == 1 && !pfpVertex.empty()){

         for(const art::Ptr<recob::Track> &trk : pfpTracks){			


            //sets track length/position related variables in ThisPrimaryDaughter
            SetTrackVariables(ThisPrimaryDaughter , trk);
            TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z())); 

            /////////////////////
            //   Truth match   //
            /////////////////////


            if(fIsData) ThisPrimaryDaughter.HasTruth = false;

            if(!fIsData){

               //loop through hits in track, find which MC particle
               //deposited most energy in track

               //get hits assoc with track
               std::vector< art::Ptr< recob::Hit> > hits = trackHitAssoc.at(trk.key());

               //use backtracker to find true pdg 

               std::unordered_map<int , double>  trkide;
               double maxe = -1;
               double tote = 0;

               int maxhits=-1;

               simb::MCParticle const* matchedParticle = NULL;

               std::vector< simb::MCParticle const*> particleVec;
               std::vector< anab::BackTrackerHitMatchingData const*> matchVec;

               //loop through hits in this track
               for (size_t i_hit = 0; i_hit < hits.size(); ++i_hit){

                  //clear vectors
                  particleVec.clear();
                  matchVec.clear();
                  particlesPerHit.get(hits[i_hit].key(), particleVec, matchVec);

                  //loop over particles that deposit energy in this hit
                  for (size_t i_particle = 0; i_particle < particleVec.size(); ++i_particle){

                     //	trkide[ particleVec[i_particle]->TrackId() ] += matchVec[i_particle]->energy;
                     trkide[ particleVec[i_particle]->TrackId() ] ++; //just increment the number of hits

                     tote += matchVec[i_particle]->energy;


                     /*
                     // old method - choose particle depositing most energy in track
                     if ( trkide[ particleVec[i_particle]->TrackId() ] > maxe ){
                     maxe = trkide[ particleVec[i_particle]->TrackId() ];
                     matchedParticle = particleVec[i_particle];

                     }
                     */


                     //new method - choose particle depositing energy in the most hits
                     if ( trkide[ particleVec[i_particle]->TrackId() ] > maxhits ){
                        maxhits = trkide[ particleVec[i_particle]->TrackId() ];
                        matchedParticle = particleVec[i_particle];

                     }

                  }

               }



               if(matchedParticle != NULL){ 

                  SimParticle P = MakeSimParticle(*matchedParticle);
                  P.Origin = getOrigin(matchedParticle->TrackId());

                  ThisPrimaryDaughter.HasTruth = true;
                  ThisPrimaryDaughter.TrackTruthPurity = maxe/tote;

                  ThisPrimaryDaughter.TrackTruePDG = P.PDG;
                  ThisPrimaryDaughter.TrackTrueE = P.E;
                  ThisPrimaryDaughter.TrackTruePx = P.Px;
                  ThisPrimaryDaughter.TrackTruePy = P.Py;
                  ThisPrimaryDaughter.TrackTruePz = P.Pz;

                  ThisPrimaryDaughter.TrackTrueModMomentum = P.ModMomentum;
                  ThisPrimaryDaughter.TrackTrueKE = P.KE;

                  ThisPrimaryDaughter.TrackTrueLength = P.Travel;

                  ThisPrimaryDaughter.TrackTrueOrigin = P.Origin;


               }
               else ThisPrimaryDaughter.HasTruth = false;

            }//if !IsData


            ///////////////////////
            // get PID for track //
            ///////////////////////

            // Setup Calo Assn
            std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());

            //Nicolo's PID

            double this_llr_pid=0;
            double this_llr_pid_score=0;
            for(auto const &calo : caloFromTrack){

               auto const &plane = calo->PlaneID().Plane;
               auto const &dedx_values = calo->dEdx();
               auto const &rr = calo->ResidualRange();
               auto const &pitch = calo->TrkPitchVec();
               std::vector<std::vector<float>> par_values;
               par_values.push_back(rr);
               par_values.push_back(pitch);

               if (calo->ResidualRange().size() == 0) continue;

               float calo_energy = 0;
               for (size_t i = 0; i < dedx_values.size(); i++)
               {
                  calo_energy += dedx_values[i] * pitch[i];
               }

               float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values, par_values, plane);
               this_llr_pid +=  llr_pid;

            }
            this_llr_pid_score = atan( this_llr_pid / 100.) * 2 / 3.14159266;

            ThisPrimaryDaughter.Track_LLR_PID = this_llr_pid_score;

            //Mean dE/dX
            std::vector<std::pair<int,double>> MeandEdXs = MeandEdX(caloFromTrack);

            double thisThreePlaneMeandEdX = ThreePlaneMeandEdX(trk,MeandEdXs);

            if( thisThreePlaneMeandEdX != thisThreePlaneMeandEdX ) std::cout << "NAN for three plane dedx!" << std::endl;

            ThisPrimaryDaughter.MeandEdX_ThreePlane = thisThreePlaneMeandEdX;

            //add single plane mean dEdX
            for(size_t i_plane=0;i_plane<MeandEdXs.size();i_plane++){

               if( MeandEdXs.at(i_plane).first == 0 ){
                  ThisPrimaryDaughter.MeandEdX_Plane0 = MeandEdXs.at(i_plane).second;
               }

               if( MeandEdXs.at(i_plane).first == 1 ){
                  ThisPrimaryDaughter.MeandEdX_Plane1 = MeandEdXs.at(i_plane).second;
               }

               if( MeandEdXs.at(i_plane).first == 2 ){
                  ThisPrimaryDaughter.MeandEdX_Plane0 = MeandEdXs.at(i_plane).second;
               }
            }

            //stock uboone PIDs
            std::vector<art::Ptr<anab::ParticleID>> trackPID = PIDAssoc.at(trk.key());

            std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

            for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){

               anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);


               //3 Plane Proton PID
               if(TMath::Abs(AlgScore.fAssumedPdg) == 2212  &&   AlgScore.fAlgName=="ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){	

                  //take log
                  ThisPrimaryDaughter.TrackPID = std::log(AlgScore.fValue);

               }

            }


         }


         ////////////////////////////
         // Get Vertex information //
         ////////////////////////////


         for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

            geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };                
            geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

            //no SC correction
            //TVector3 pos( vtx->position().X() , vtx->position().Y() , vtx->position().Z() );

            //w SC correction - forward
            TVector3 pos( vtx->position().X() + sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z() );

            ThisPrimaryDaughter.SetVertex( pos );

            ThisPrimaryDaughter.Displacement = (pos - fRecoPrimaryVertex).Mag();

         }

      }



      if(ThisPrimaryDaughter.PDG == 13 && pfpTracks.size() == 1){

         fTrackPrimaryDaughters.push_back( ThisPrimaryDaughter );
         fNPrimaryTrackDaughters++;
      }
      else if(ThisPrimaryDaughter.PDG == 11 && pfpShowers.size() == 1){ 
         fShowerPrimaryDaughters.push_back( ThisPrimaryDaughter );
         fNPrimaryShowerDaughters++;
      }





   }//end of PFP loop

   //set indices in particle vectors	
   for(size_t i_tr=0;i_tr<fTrackPrimaryDaughters.size();i_tr++) fTrackPrimaryDaughters[i_tr].Index = i_tr;
   for(size_t i_sh=0;i_sh<fShowerPrimaryDaughters.size();i_sh++) fShowerPrimaryDaughters[i_sh].Index = i_sh;

   // Run connectedness test

   art::Handle<std::vector<recob::Wire>> wireHandle;
   std::vector<art::Ptr<recob::Wire>> wireVect;

   // Fill Wire vector
   if(e.getByLabel(fWireLabel,wireHandle)) art::fill_ptr_vector(wireVect,wireHandle);
   else
      std::cout << "Wire handle not setup" << std::endl;

   if(fTrackPrimaryDaughters.size() >= 3){

      Conn_Helper.LoadWireActivity(wireVect);

      Conn_Helper.AddStartPositions(TrackStarts);

      std::vector<ConnectednessOutcome> Outcomes = Conn_Helper.RunTest();

      for(size_t i=0;i<Outcomes.size();i++){

         ConnectednessOutcome Outcome = Outcomes.at(i);

         std::vector<int> this_SeedIndexes = Outcome.SeedIndexes;
         std::vector<int> this_OutputIndexes = Outcome.OutputIndexes;
         std::vector<int> this_OutputSizes = Outcome.OutputSizes;
         std::vector<int> this_SeedChannels = Outcome.SeedChannels;
         std::vector<int> this_SeedTicks = Outcome.SeedTicks;

         if(Outcome.Plane == 0){
            Conn_SeedIndexes_Plane0.push_back(this_SeedIndexes);
            Conn_OutputIndexes_Plane0.push_back(this_OutputIndexes);
            Conn_OutputSizes_Plane0.push_back(this_OutputSizes);
            Conn_SeedChannels_Plane0.push_back(this_SeedChannels);
            Conn_SeedTicks_Plane0.push_back(this_SeedTicks);
         } 
         else if(Outcome.Plane == 1){
            Conn_SeedIndexes_Plane1.push_back(this_SeedIndexes);
            Conn_OutputIndexes_Plane1.push_back(this_OutputIndexes);
            Conn_OutputSizes_Plane1.push_back(this_OutputSizes);
            Conn_SeedChannels_Plane1.push_back(this_SeedChannels);
            Conn_SeedTicks_Plane1.push_back(this_SeedTicks);
         } 
         else if(Outcome.Plane == 2){
            Conn_SeedIndexes_Plane2.push_back(this_SeedIndexes);
            Conn_OutputIndexes_Plane2.push_back(this_OutputIndexes);
            Conn_OutputSizes_Plane2.push_back(this_OutputSizes);
            Conn_SeedChannels_Plane2.push_back(this_SeedChannels);
            Conn_SeedTicks_Plane2.push_back(this_SeedTicks);
         } 


      }


   }

   //store truth matching info for muon, decay proton and pion
   StoreTrackTruth();

   if(fPrint) PrintInfo();

   FinishEvent();

}



///////////////////////////////////////////
// Print some useful info from the event //
///////////////////////////////////////////

void hyperon::HyperonSelection::PrintInfo(){


   if( fHyperon.size() != 1 ) return;

   if( !(fHyperon.at(0).PDG == fPrintPdg || fPrintPdg == -1) || !fInActiveTPC) return;

   std::cout << std::endl;	
   std::cout << "EventID: " << fEventID-1 << std::endl;

   std::cout << "Truth info" << std::endl;

   std::cout << "Hyperon: " << std::endl;
   for(size_t i=0;i<fHyperon.size();i++) fHyperon.at(i).Print();

   std::cout << "Lepton: " << std::endl;
   for(size_t i=0;i<fLepton.size();i++) fLepton.at(i).Print();

   std::cout << "Hyperon decay products: " << std::endl;
   for(size_t i=0;i<fDecay.size();i++) fDecay.at(i).Print();

   if(fDecay.size() == 2){
      std::cout << "Decay Opening Angle:  " << fDecayOpeningAngle << std::endl;
      std::cout << "Angle between lepton and pion:  " << fLeptonPionAngle << std::endl; //openining angle between lepton and pion
      std::cout << "Angle between lepton and nucleon:  " <<   fLeptonNucleonAngle << std::endl; //opening angle between lepton and nucleon

   }



   std::cout << std::endl;
   std::cout << "Reco info" << std::endl;
   std::cout << "Num tracklike daughters: " << fNPrimaryTrackDaughters << "  Num showerlike daughters:  " << fNPrimaryShowerDaughters << std::endl;
   std::cout << "List tracklike of daughters: " << std::endl;

   for(size_t i=0;i<fTrackPrimaryDaughters.size();i++){

      fTrackPrimaryDaughters.at(i).Print();

   }



   std::cout << std::endl;
}


/////////////////////////////////////////////
// Check origin of particle by its TrackId //
/////////////////////////////////////////////

int hyperon::HyperonSelection::getOrigin(int idnum){

   //search list of primaries
   for(size_t i_p=0;i_p<primary_IDs.size();i_p++){
      if(primary_IDs.at(i_p) == idnum){ 
         return 1;
      }
   }

   //search list of hyperon decay products	
   for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){
      if(daughter_IDs.at(i_d) == idnum){ 
         return 2;
      }
   }

   //search list of kaon decay products	
   for(size_t i_d=0;i_d<Kaon_daughter_IDs.size();i_d++){
      if(Kaon_daughter_IDs.at(i_d) == idnum){ 
         return 4;
      }
   }

   return 3;

}





///////////////////////////////////////////////////////////////	
// Finished processing event - update Metadata and fill tree //
///////////////////////////////////////////////////////////////

void hyperon::HyperonSelection::FinishEvent(){

   if(fDebug) std::cout << "Finishing Event" << std::endl;

   if(!fIsData){

      if(fCCNC == "CC") fNChargedCurrent++;
      else fNNeutralCurrent++;

      if( fNeutrino.size() != 1 ) std::cout << "Number of simulated neutrinos in this event != 1 !!" << std::endl;
      else if(fNeutrino.at(0).PDG == 12) fNnue++;
      else if(fNeutrino.at(0).PDG == 14) fNnuMu++;
      else if(fNeutrino.at(0).PDG == -12) fNnueBar++;
      else if(fNeutrino.at(0).PDG == -14) fNnuMuBar++;


      //genie uses QEL for hyperon events, NuWro uses HYP
      if(fNeutrino.size() == 1 && fInActiveTPC && fIsLambdaCharged && fNeutrino.at(0).PDG == -14 && ( fMode == "QEL" || fMode == "HYP")){ 

         //add kinematic thresholds
         if(fDecay.at(0).PDG == 2212 && fDecay.at(0).ModMomentum > fDecayProtonThreshold && fDecay.at(1).PDG == -211 && fDecay.at(1).ModMomentum > fDecayPionThreshold) fIsSignal = true;
         else if(fDecay.at(1).PDG == 2212 && fDecay.at(1).ModMomentum > fDecayProtonThreshold && fDecay.at(0).PDG == -211 && fDecay.at(0).ModMomentum > fDecayPionThreshold) fIsSignal = true;
         else fIsSignal = false;

      } 

      //asscociated hyperon/kaon production tagger - hyperon and kaon in the final state
      if( fHyperon.size() && fPrimaryKaon.size() ) fIsAssociatedHyperon = true;

      //if is a signal event, check if decay products were reconstructed
      if(fIsSignal && fTrueDecayProtonIndex >= 0 && fTrueDecayPionIndex >= 0) fGoodReco = true;	

   }

   //Perform selection
   fSelectedEvent = PerformSelection();

   //store info for this event
   fOutputTree->Fill();

   //update metadata

   fNEvents++; //total events in sample
   if(fIsHyperon && fInActiveTPC) fNHyperons++; //total hyperon events in active vol
   if(fSelectedEvent) fNSelectedEvents++; //total events passing selection
   if(fSelectedEvent && fInActiveTPC && fIsHyperon) fNSelectedHyperons++; //total hyperons in active vol passing preselection

   //signal events
   if(fIsSignal) fNSignal++;
   if(fIsSignal && fSelectedEvent) fNSelectedSignal++;

   //signal events
   if(fGoodReco) fNGoodReco++;
   if(fGoodReco && fSelectedEvent) fNSelectedGoodReco++;

}


/////////////////////////////////////////////////////////////////
//   Find tracks truth matching to muon, proton and pion from  //
//   hyperon decay (if they exist), store positions in track   //
//   vector in tree                                            //
/////////////////////////////////////////////////////////////////

void hyperon::HyperonSelection::StoreTrackTruth(){

   //make sure index stores have been reset!
   fTrueMuonIndex=-1;	
   fTrueDecayProtonIndex=-1;
   fTrueDecayPionIndex=-1;

   //can be multiple tracks corresponding to the muon/proton/pion
   //first get all the tracks truth matching to primary muon, decay proton, decay pion

   std::vector<int> Muons;
   std::vector<int> DecayProtons;
   std::vector<int> DecayPions;

   for(size_t i_tr=0;i_tr<fTrackPrimaryDaughters.size();i_tr++){

      //if track does not have a truth matching, skip
      if(!fTrackPrimaryDaughters.at(i_tr).HasTruth) continue;

      //if muon produced at primary vertex
      if( abs(fTrackPrimaryDaughters.at(i_tr).TrackTruePDG) == 13 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 1 ) Muons.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

      //if proton produced from hyperon decay
      if( fTrackPrimaryDaughters.at(i_tr).TrackTruePDG == 2212 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 2 ) DecayProtons.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

      //if pi minus produced from hyperon decay
      if( fTrackPrimaryDaughters.at(i_tr).TrackTruePDG == -211 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 2 ) DecayPions.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

   }

   //if there are no muons, protons or pions, exit here
   if( !Muons.size() && !DecayProtons.size() && !DecayPions.size() ) return;

   //set muon information
   //if multiple muons found, choose the one closest to reco'd primary vertex
   double min_dist=10000;
   for(size_t i_m=0;i_m<Muons.size();i_m++){

      TVector3 MuonStart(fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartX,fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartY,fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartZ);

      double d = (MuonStart - fRecoPrimaryVertex).Mag();

      if(d < min_dist) { fTrueMuonIndex = Muons.at(i_m); min_dist = d;  }

   }

   //if there are no decay products, exit here
   if( !DecayProtons.size() && !DecayPions.size() ) return;


   //check you have both decay products, if you have both, choose the pair of tracks starting closest together

   if( DecayProtons.size() && DecayPions.size() ){

      //std::cout << "Got both decay products" << std::endl;	

      double min_sep = 10000;

      for(size_t i_pr=0;i_pr<DecayProtons.size();i_pr++){
         for(size_t i_pi=0;i_pi<DecayPions.size();i_pi++){

            TVector3 ProtonStart(fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartX,fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartY,fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartZ);

            TVector3 PionStart(fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartX,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartY,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartZ);

            double d = (ProtonStart - PionStart).Mag();
            if( d < min_sep ){ fTrueDecayProtonIndex = DecayProtons.at(i_pr); fTrueDecayPionIndex = DecayPions.at(i_pi); min_sep = d; }


         }
      }

   }

   //if missing one or both of the decay products, choose the tracks starting closest to PV
   else {

      //for protons
      min_dist = 10000;
      for(size_t i_pr=0;i_pr<DecayProtons.size();i_pr++){

         TVector3 ProtonStart(fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartX,fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartY,fTrackPrimaryDaughters.at(DecayProtons.at(i_pr)).TrackStartZ);

         double d = (ProtonStart - fRecoPrimaryVertex).Mag();

         if(d < min_dist) { fTrueDecayProtonIndex = DecayProtons.at(i_pr); min_dist = d;  }

      }

      //for pions
      min_dist = 10000;
      for(size_t i_pi=0;i_pi<DecayPions.size();i_pi++){

         TVector3 PionStart(fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartX,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartY,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartZ);

         double d = (PionStart - fRecoPrimaryVertex).Mag();

         if(d < min_dist) { fTrueDecayPionIndex = DecayPions.at(i_pi); min_dist = d;  }

      }


   }


}

//////////////////////////////////////////////////////////////////


void hyperon::HyperonSelection::beginJob(){


   if(fDebug) std::cout << "Begin job" << std::endl;



   // Implementation of optional member functions

   art::ServiceHandle<art::TFileService> tfs;

   fOutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");

   fOutputTree->Branch("IsData",&fIsData);

   //write a second tree to store decay product info
   fOutputTree->Branch("EventID",&fEventID);
   fOutputTree->Branch("run",&run);
   fOutputTree->Branch("subrun",&subrun);
   fOutputTree->Branch("event",&event);

   fOutputTree->Branch("Mode",&fMode);
   fOutputTree->Branch("CCNC",&fCCNC);

   //generator truth info
   fOutputTree->Branch("NMCTruths",&fNMCTruths);
   fOutputTree->Branch("InActiveTPC",&fInActiveTPC);
   fOutputTree->Branch("IsHyperon",&fIsHyperon);
   fOutputTree->Branch("IsSigmaZero",&fIsSigmaZero);
   fOutputTree->Branch("IsLambda",&fIsLambda);
   fOutputTree->Branch("IsLambdaCharged",&fIsLambdaCharged);
   fOutputTree->Branch("IsSignal",&fIsSignal);
   fOutputTree->Branch("GoodReco",&fGoodReco);
   fOutputTree->Branch("IsAssociatedHyperon",&fIsAssociatedHyperon);
   fOutputTree->Branch("SelectedEvent",&fSelectedEvent);

   fOutputTree->Branch("Weight",&fWeight);

   //neutrino information
   fOutputTree->Branch("Neutrino","vector<SimParticle>",&fNeutrino);
   fOutputTree->Branch("TruePrimaryVertex","TVector3",&fTruePrimaryVertex);

   //lepton
   fOutputTree->Branch("Lepton","vector<SimParticle>",&fLepton);
   //hyperon
   fOutputTree->Branch("Hyperon","vector<SimParticle>",&fHyperon);

   //nucleons, pions and kaons at the primary vtx
   fOutputTree->Branch("PrimaryNucleon","vector<SimParticle>",&fPrimaryNucleon);
   fOutputTree->Branch("PrimaryPion","vector<SimParticle>",&fPrimaryPion);
   fOutputTree->Branch("PrimaryKaon","vector<SimParticle>",&fPrimaryKaon);

   //hyperon decay info
   fOutputTree->Branch("DecayVertex","TVector3",&fDecayVertex); //position of hyperon decay vertex
   fOutputTree->Branch("Decay","vector<SimParticle>",&fDecay);

   //sigma zero decay products
   fOutputTree->Branch("SigmaZeroDecayPhoton","vector<SimParticle>",&fSigmaZeroDecayPhoton);
   fOutputTree->Branch("SigmaZeroDecayLambda","vector<SimParticle>",&fSigmaZeroDecayLambda);

   fOutputTree->Branch("KaonDecay","vector<SimParticle>",&fKaonDecay);

   //angles between particles involved

   fOutputTree->Branch("DecayOpeningAngle",&fDecayOpeningAngle); //opening angle between hyperon decay products     
   fOutputTree->Branch("LeptonPionAngle",&fLeptonPionAngle); //openining angle between lepton and pion
   fOutputTree->Branch("LeptonNucleonAngle",&fLeptonNucleonAngle); //opening angle between lepton and nucleon

   //////////////////////
   // Reco Information //
   //////////////////////

   fOutputTree->Branch("RecoPrimaryVertex","TVector3",&fRecoPrimaryVertex);
   fOutputTree->Branch("NPrimaryTrackDaughters",&fNPrimaryTrackDaughters);
   fOutputTree->Branch("NPrimaryShowerDaughters",&fNPrimaryShowerDaughters);

   fOutputTree->Branch("TracklikePrimaryDaughters","vector<RecoParticle>",&fTrackPrimaryDaughters);
   fOutputTree->Branch("ShowerlikePrimaryDaughters","vector<RecoParticle>",&fShowerPrimaryDaughters);
   fOutputTree->Branch("MuonIndex",&fMuonIndex);

   fOutputTree->Branch("TrueMuonIndex",&fTrueMuonIndex);
   fOutputTree->Branch("TrueDecayProtonIndex",&fTrueDecayProtonIndex);
   fOutputTree->Branch("TrueDecayPionIndex",&fTrueDecayPionIndex);

   ////////////////////////////
   //   Connectedness test   //
   ////////////////////////////

   fOutputTree->Branch("ConnSeedIndexes_Plane0",&Conn_SeedIndexes_Plane0);
   fOutputTree->Branch("ConnOutputIndexes_Plane0",&Conn_OutputIndexes_Plane0);
   fOutputTree->Branch("ConnOutputSizes_Plane0",&Conn_OutputSizes_Plane0);
   fOutputTree->Branch("ConnSeedChannels_Plane0",&Conn_SeedChannels_Plane0);
   fOutputTree->Branch("ConnSeedTicks_Plane0",&Conn_SeedTicks_Plane0);

   fOutputTree->Branch("ConnSeedIndexes_Plane1",&Conn_SeedIndexes_Plane1);
   fOutputTree->Branch("ConnOutputIndexes_Plane1",&Conn_OutputIndexes_Plane1);
   fOutputTree->Branch("ConnOutputSizes_Plane1",&Conn_OutputSizes_Plane1);
   fOutputTree->Branch("ConnSeedChannels_Plane1",&Conn_SeedChannels_Plane1);
   fOutputTree->Branch("ConnSeedTicks_Plane1",&Conn_SeedTicks_Plane1);

   fOutputTree->Branch("ConnSeedIndexes_Plane2",&Conn_SeedIndexes_Plane2);
   fOutputTree->Branch("ConnOutputIndexes_Plane2",&Conn_OutputIndexes_Plane2);
   fOutputTree->Branch("ConnOutputSizes_Plane2",&Conn_OutputSizes_Plane2);
   fOutputTree->Branch("ConnSeedChannels_Plane2",&Conn_SeedChannels_Plane2);
   fOutputTree->Branch("ConnSeedTicks_Plane2",&Conn_SeedTicks_Plane2);

   //////////////////////////////////////////
   //            Metadata Tree		//
   //////////////////////////////////////////

   fNEvents=0; //total events in sample
   fNHyperons=0;

   fNSelectedEvents=0; //total events passing selection
   fNSelectedHyperons=0; //total true hyperon events passing selection

   fHyperonEfficiency=0; //hyperon selection efficiency
   fHyperonPurity=0;  //hyperon selection purity
   fHyperonTruePurity=0; //hyperon purity after converting from enriched sample to real sample
   fHyperonEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
   fHyperonEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

   fNSignal=0;
   fNSelectedSignal=0;

   fNGoodReco=0;
   fNSelectedGoodReco=0;

   fNChargedCurrent=0; //number of cc events
   fNNeutralCurrent=0; //number of nc events
   fNnuMu=0; //number of numu events
   fNnue=0; //number of nue events
   fNnuMuBar=0; //number of numubar events
   fNnueBar=0; //number of nuebar events

   fPOT=0;

   fMetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");

   fMetaTree->Branch("TotalEvents",&fNEvents);
   fMetaTree->Branch("NSelectedEvents",&fNSelectedEvents);

   fMetaTree->Branch("NChargedCurrent",&fNChargedCurrent); //number of cc events
   fMetaTree->Branch("NNeutralCurrent",& fNNeutralCurrent); //number of nc events
   fMetaTree->Branch("Nnumu",&fNnuMu); //number of numu events
   fMetaTree->Branch("Nnue",&fNnue); //number of nue events
   fMetaTree->Branch("NnuMuBar",&fNnuMuBar); //number of numubar events
   fMetaTree->Branch("NnueBar",&fNnueBar); //number of nuebar events


   fMetaTree->Branch("NHyperons",&fNHyperons);
   fMetaTree->Branch("NSelectedHyperons",&fNSelectedHyperons);
   fMetaTree->Branch("HyperonEfficiency",&fHyperonEfficiency);
   fMetaTree->Branch("HyperonPurity",&fHyperonPurity);
   fMetaTree->Branch("HyperonEfficiencyTimesPurity",&fHyperonEfficiencyTimesPurity);
   fMetaTree->Branch("HyperonTruePurity",&fHyperonTruePurity);
   fMetaTree->Branch("HyperonEfficiencyTimesTruePurity",&fHyperonEfficiencyTimesTruePurity);


   fMetaTree->Branch("NSignal",&fNSignal);
   fMetaTree->Branch("NSelectedSignal",&fNSelectedSignal);	

   fMetaTree->Branch("SignalEfficiency",&fSignalEfficiency);
   fMetaTree->Branch("SignalPurity",&fSignalPurity);
   fMetaTree->Branch("SignalEfficiencyTimesPurity",&fSignalEfficiencyTimesPurity);
   fMetaTree->Branch("SignalTruePurity",&fSignalTruePurity);
   fMetaTree->Branch("SignalEfficiencyTimesTruePurity",&fSignalEfficiencyTimesTruePurity);

   fMetaTree->Branch("NGoodReco",&fNGoodReco);
   fMetaTree->Branch("NSelectedGoodReco",&fNSelectedGoodReco);	

   fMetaTree->Branch("GoodRecoEfficiency",&fGoodRecoEfficiency);
   fMetaTree->Branch("GoodRecoPurity",&fGoodRecoPurity);
   fMetaTree->Branch("GoodRecoEfficiencyTimesPurity",&fGoodRecoEfficiencyTimesPurity);
   fMetaTree->Branch("GoodRecoTruePurity",&fGoodRecoTruePurity);
   fMetaTree->Branch("GoodRecoEfficiencyTimesTruePurity",&fGoodRecoEfficiencyTimesTruePurity);

   fMetaTree->Branch("POT",&fPOT);

   if(fDebug) std::cout << "Finished begin job" << std::endl;


}



void hyperon::HyperonSelection::endJob()
{
   //calculate efficiency and purity

   if(fNSelectedEvents > 0){
      fHyperonPurity = (double)fNSelectedHyperons/fNSelectedEvents;
      fSignalPurity = (double)fNSelectedSignal/fNSelectedEvents;	
      fGoodRecoPurity = (double)fNSelectedGoodReco/fNSelectedEvents;
   }

   if(fNHyperons > 0) fHyperonEfficiency = (double)fNSelectedHyperons/fNHyperons;

   if(fNSelectedSignal > 0) fSignalEfficiency = (double)fNSelectedSignal/fNSignal;

   if(fNSelectedGoodReco > 0) fGoodRecoEfficiency = (double)fNSelectedGoodReco/fNGoodReco;

   fHyperonEfficiencyTimesPurity = fHyperonPurity * fHyperonEfficiency;
   fHyperonTruePurity = 1/( 1+(1/fHyperonPurity-1)*30 );
   fHyperonEfficiencyTimesTruePurity = fHyperonEfficiency * fHyperonTruePurity;

   fSignalEfficiencyTimesPurity = fSignalPurity * fSignalEfficiency;
   fSignalTruePurity = 1/( 1+(1/fSignalPurity-1)*30 );
   fSignalEfficiencyTimesTruePurity = fSignalEfficiency * fSignalTruePurity;

   fGoodRecoEfficiencyTimesPurity = fGoodRecoPurity * fGoodRecoEfficiency;
   fGoodRecoTruePurity = 1/( 1+(1/fGoodRecoPurity-1)*30 );
   fGoodRecoEfficiencyTimesTruePurity = fGoodRecoEfficiency * fGoodRecoTruePurity;

   //print useful metadata

   std::cout << std::endl << std::endl<< std::endl << "Some Metadata for this sample:" << std::endl;

   std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
   std::cout << "Selected Hyperons in active volume: " << fNSelectedHyperons << std::endl;
   std::cout << "Enriched Sample Hyperon Efficiency: " << fHyperonEfficiency*100 << std::endl;
   std::cout << "Enriched Sample Hyperon Purity: " << fHyperonPurity*100 << std::endl;	
   std::cout << "Enriched Sample Hyperon E x P: " << fHyperonEfficiencyTimesPurity*100 << std::endl;
   std::cout << "True Hyperon Purity assuming enrichment factor of 30: " << fHyperonTruePurity*100 <<  std::endl;

   std::cout << std::endl;
   std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
   std::cout << "Signal Events: " << fNSignal << std::endl;
   std::cout << "Selected Signal events: " << fNSelectedSignal  << std::endl;	
   std::cout << "Enriched Sample Efficiency: " << fSignalEfficiency*100 << std::endl;
   std::cout << "Enriched Sample Purity: " << fSignalPurity*100 << std::endl;	
   std::cout << "Enriched Sample E x P: " << fSignalEfficiencyTimesPurity*100 << std::endl;
   std::cout << "True Purity assuming enrichment factor of 30: " << fSignalTruePurity*100 <<  std::endl;



   std::cout << std::endl;
   std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
   std::cout << "GoodReco Events: " << fNGoodReco << std::endl;
   std::cout << "Selected GoodReco events: " << fNSelectedGoodReco  << std::endl;	
   std::cout << "Enriched Sample Efficiency: " << fGoodRecoEfficiency*100 << std::endl;
   std::cout << "Enriched Sample Purity: " << fGoodRecoPurity*100 << std::endl;	
   std::cout << "Enriched Sample E x P: " << fGoodRecoEfficiencyTimesPurity*100 << std::endl;
   std::cout << "True Purity assuming enrichment factor of 30: " << fGoodRecoTruePurity*100 <<  std::endl;


   std::cout << std::endl << std::endl << std::endl;

   fMetaTree->Fill();


}



void hyperon::HyperonSelection::beginSubRun(const art::SubRun& sr)
{

   if(fDebug) std::cout << "Getting Subrun POT Info" << std::endl;

   art::Handle<sumdata::POTSummary> POTHandle;

   if(sr.getByLabel(fPOTSummaryLabel,POTHandle)) fPOT += POTHandle->totpot;	

   //std::cout << "Adding subrun POT of " << POTHandle->totpot << std::endl;

}


void hyperon::HyperonSelection::endSubRun(const art::SubRun& sr){}


DEFINE_ART_MODULE(hyperon::HyperonSelection)
