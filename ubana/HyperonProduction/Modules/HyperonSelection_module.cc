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

//#include "art_root_io/TFileService.h"
//#include "art_root_io/TFileDirectory.h"
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

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"

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

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TVector3.h"
#include "TLorentzVector.h"

//local includes
#include "ubana/HyperonProduction/util/SimParticle.h"
#include "ubana/HyperonProduction/util/RecoParticle.h"

#include "ubana/HyperonProduction/Alg/ParticleTypes.h"
#include "ubana/HyperonProduction/Alg/FV.h"
#include "ubana/HyperonProduction/Alg/Muon_ID.h"
#include "ubana/HyperonProduction/Alg/Track_Length_Cut.h"
#include "ubana/HyperonProduction/Alg/Gap_Cut.h"

#include "ubana/HyperonProduction/util/Helpers.h"

#include "TRandom2.h"


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

		int getOrigin(int idnum);



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

		TTree * fPIDTree; //PID Scores

		TTree * fMetaTree; //metadata


		///////////////////////////
		//   Truth level info    //
		///////////////////////////

		//generator truth info

		std::string fMode; //interaction mode
		std::string fCCNC; //charged current/neutral current

		bool fInActiveTPC=false; //is event in active TPC

		bool fIsHyperon=false; //true if event contains a hyperon
		bool fIsSigmaZero=false; //true if event hyperon was a sigma zero		
		bool fIsLambda=false;
		bool fIsLambdaCharged; //true if event contains Lambda decaying into p + pi-
		bool fIsSignal=false; //true if numu event in fiducial vol producing Lambda decaying to p + pi- -->Signal to search for		

		//lepton produced at primary vtx		
		std::vector<SimParticle> fLepton;

		//hyperon produced at primary vtx
		std::vector<SimParticle> fHyperon;

		//nucleons produced at primary vtx
		std::vector<SimParticle> fPrimaryNucleon;

		//pions produced at primary vtx
		std::vector<SimParticle> fPrimaryPion;

		//vertex information
		TVector3 fTruePrimaryVertex;

		//neutrino information
		double fNuEnergy;
		int fNuPDG;

		//g4 truth info
		TVector3 fDecayVertex;

		std::vector<SimParticle> fDecay; //hyperon decay products

		std::vector<SimParticle> fSigmaZeroDecayPhoton; //photon produced by Sigma zero decay
		std::vector<SimParticle> fSigmaZeroDecayLambda;	//lambda produced by Sigma zero decay

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

		//data storage (not to be written to tree)		
		std::vector<int> fRecoTrackTruthMatchedID; //list of ID numbers of truth matched particles


		////////////////////////
		//     PID INFO       //
		////////////////////////

		std::vector<std::vector<int>> fAlgNumber; //pid scores for tracks 
		std::vector<std::vector<int>> fPlaneNo; //pid scores for tracks 
		std::vector<std::vector<int>> fAssumedPDG; //pid scores for tracks 
		std::vector<std::vector<double>> fScore; //pid scores for tracks 
		std::vector<double> fRecoLength; //pid scores for tracks 
		std::vector<double> fTruePDG; //pid scores for tracks 
		std::vector<double> fTrueLength; //pid scores for tracks 





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

		//reco level metadata

		int fNSelectedEvents; //total events passing selection

		int fNSelectedHyperons; //total true hyperon events passing selection

		int fNSelectedSignal; //number of signal events passing selection


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


		///////////////////////
		//   Module Labels   //
		///////////////////////

		std::string fGenieGenModuleLabel;
		std::string fGeantModuleLabel;
		std::string fTrackLabel;
		std::string fShowerLabel;
		std::string fVertexLabel;
		std::string fPFParticleLabel;
		std::string fPIDLabel;
		std::string fCaloLabel;
		std::string fHitLabel;
		std::string fTrackHitAssnLabel;
		std::string fHitTruthAssnLabel;
		std::string fMetadataLabel;
		std::string fShowerHitAssnLabel;


		//misc
		bool fPrint;
		int fPrintPdg=-1; //-1 means print everything


		///////////////////////////////////
		// Hyperon selection parameters  //
		///////////////////////////////////

		int fMinDaughters;//minimum number of daughters
		int fMaxPrimaryShowerDaughters;
		int fMinDaughterTracks; //minimum number of track like primary daughters
		double fMuonPIDScore;
		double fSecondaryTrackLengthCut;
		double fTertiaryTrackLengthCut;

		double fTrackScore; //minimum track/shower score required for pfp to be classified as track


};


////////////////////////////////////////////////////
// Decides if an event passes selection criteria  //
////////////////////////////////////////////////////

//TODO: Make proper system for storing cut history in trees/metadata

bool hyperon::HyperonSelection::PerformSelection(){


	//if(fRecoPrimaryVertex.X() == -1000 || fRecoPrimaryVertex.Y() == -1000 || fRecoPrimaryVertex.Y() == -1000) return false; //no reco'd neutrino 

	//Apply fiducial volume cut

	if(!inActiveTPC(fRecoPrimaryVertex)) return false;

	if(fNPrimaryDaughters < fMinDaughters) return false;
	//	std::cout << "Passed Min Daughter cut" << std::endl;

	if(fNPrimaryShowerDaughters > fMaxPrimaryShowerDaughters) return false;
	//	std::cout << "Pass Max Shower cut" << std::endl;

	if(fNPrimaryTrackDaughters < fMinDaughterTracks) return false;
	//	std::cout << "Passed Min Track cut" << std::endl;


	//Muon PID
	int i_muon = Muon_ID(fTrackPrimaryDaughters,fMuonPIDScore);
	if(i_muon == -1) return false;
	fMuonIndex = i_muon;

	//Secondary and Tertiary Track Length Cuts
	if(!Track_Length_Cut(fTrackPrimaryDaughters,i_muon,fSecondaryTrackLengthCut,fTertiaryTrackLengthCut)) return false;


	//useful if you want to get ED of events passing selection

	/*
	   std::cout << "Event Num: " << event << "   ";
	   if(fIsSignal) std::cout << "__ACCEPTED_SIGNAL__"  << "  ";
	   else std::cout << "__ACCEPTED_BACKGROUND__" << "  ";   	
	   std::cout << fMode << std::endl; 
	   */

	return true;


}


////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////



hyperon::HyperonSelection::HyperonSelection(fhicl::ParameterSet const& p)
	: EDAnalyzer{p}   // ,
	// More initializers here.
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
	fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
	fPFParticleLabel = p.get<std::string>("PFParticleLabel");
	fGeantModuleLabel = p.get<std::string>("GeantLabel");
	fVertexLabel = p.get<std::string>("VertexLabel");
	fTrackLabel = p.get<std::string>("TrackLabel");
	fShowerLabel = p.get<std::string>("ShowerLabel");
	fPIDLabel = p.get<std::string>("PIDLabel");
	fCaloLabel = p.get<std::string>("CaloLabel");
	fHitLabel  = p.get<std::string>("HitLabel");
	fTrackHitAssnLabel = p.get<std::string>("TrackHitAssnLabel");
	fHitTruthAssnLabel = p.get<std::string>("HitTruthAssnLabel");
	fShowerHitAssnLabel = p.get<std::string>("ShowerHitAssnLabel");
	fMetadataLabel = p.get<std::string>("MetadataLabel");

	//selection parameters
	fMinDaughters = p.get<int>("MinDaughters",0);	
	fMinDaughterTracks = p.get<int>("MinDaughterTracks",0);
	fMaxPrimaryShowerDaughters = p.get<int>("MaxDaughterShowers",1000);
	fMuonPIDScore = p.get<double>("MuonPIDScore",0.6);
	fSecondaryTrackLengthCut = p.get<double>("SecondaryTrackLengthCut",65);
	fTertiaryTrackLengthCut = p.get<double>("TertiaryTrackLengthCut",35);

	//misc
	fTrackScore = p.get<double>("TrackScore",0.5);




	//misc parameters
	fPrint = p.get<bool>("Print",false);
	fPrintPdg = p.get<int>("PrintPdg",-1);

}

void hyperon::HyperonSelection::analyze(art::Event const& e)
{



	//begin by resetting everything

	//Generator Info

	fInActiveTPC=false;
	fIsHyperon=false; //true if event contains a hyperon
	fIsSigmaZero=false; //true if event hyperon was a sigma zero		
	fIsSignal=false;	

	//G4 Info

	//lepton
	fLepton.clear();

	//hyperon
	fHyperon.clear();

	//nucleons and pions produced at primary vtx
	fPrimaryNucleon.clear();
	fPrimaryPion.clear();

	//vertex information
	fTruePrimaryVertex.SetXYZ(-1000,-1000,-1000);

	//neutrino information
	fNuEnergy=-1; //neutrino energy
	fNuPDG=0; //neutrino PDG code

	//g4 truth info
	fDecayVertex.SetXYZ(-1000,-1000,-1000);

	//hyperon decay products
	fDecay.clear();

	//sigma zero decay products
	fSigmaZeroDecayPhoton.clear();
	fSigmaZeroDecayLambda.clear();

	fDecayOpeningAngle=-1; //opening angle between hyperon decay products   
	fLeptonPionAngle=-1; //openining angle between lepton and pion
	fLeptonNucleonAngle=-1; //opening angle between lepton and nucleon


	//Reco Info

	fSelectedEvent = false; //true if event passes some selection criteria	
	fNPrimaryDaughters = 0; //number of primary daughters

	fNPrimaryTrackDaughters=0; //num of track like primary daughters
	fNPrimaryShowerDaughters=0; //num of shower like primary daughters

	fRecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

	fTrackPrimaryDaughters.clear();
	fShowerPrimaryDaughters.clear();
	fMuonIndex=-1;

	//PID Info
	fAlgNumber.clear(); //pid scores for tracks 
	fPlaneNo.clear(); //pid scores for tracks 
	fAssumedPDG.clear(); //pid scores for tracks 
	fScore.clear(); //pid scores for tracks 
	fRecoLength.clear(); //pid scores for tracks 
	fTruePDG.clear(); //pid scores for tracks 
	fTrueLength.clear(); //pid scores for tracks 





	//Event ID information

	fEventID = e.id().event();
	run = e.run();
	subrun = e.subRun();
	event = e.event();


	//////////////////////////////
	// GET EVENT GENERATOR INFO //
	//////////////////////////////


	art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
	std::vector<art::Ptr<simb::MCTruth> > mcTrVect;

	if(e.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){

		art::fill_ptr_vector(mcTrVect,mctruthListHandle);

	}

	if(!mcTrVect.size()) {
		std::cout << "EMPTY EVENT" << std::endl;
		FinishEvent();
		return;
	}



	for(const art::Ptr<simb::MCTruth> &MCtruth : mcTrVect){

		//get the neutrino info

		simb::MCNeutrino Nu = MCtruth->GetNeutrino();
		mode = Nu.Mode();
		ccnc = Nu.CCNC();

		//TODO: Find out why this is always 0
		//std::cout << ccnc << std::endl;

		if(mode == 0) fMode = "QEL";
		else if(mode == 1) fMode = "RES";
		else if(mode == 2) fMode = "DIS";
		else if(mode == 3) fMode = "COH";
		else if(mode == 10) fMode = "MEC";
		else if(mode == 1095) fMode = "HYP";
		else fMode = "Other";	



		for(int k_particles=0;k_particles<MCtruth->NParticles();k_particles++){


			simb::MCParticle Part = MCtruth->GetParticle(k_particles);

			//Get list of particles from true PV, if lepton or neutrino set PV
			if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1){ 


				fTruePrimaryVertex.SetXYZ( Part.Vx() , Part.Vy() , Part.Vz() );

				fInActiveTPC=inActiveTPC(fTruePrimaryVertex);



			}//if lepton


			//neutrino in the interaction
			if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){

				fNuEnergy = Part.E();
				fNuPDG =  Part.PdgCode();

			}



		}//k_particles (particles associated with MC truth)




	}//i (num of MC truths in event)



	////////////////////////////////////////////
	//Get Geant Information
	///////////////////////////////////////////



	//if event contains a hyperon, find decay vertex and decay products


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

		SimParticle P = MakeSimParticle(part);
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



	}



	//check all decay products are actually produced at decay vertex - sometimes you get some electrons thrown in

	if(fHyperon.size() == 1){

		//	std::cout << "Checking hyperon vector for fake daughters" << std::endl;

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




	//now go through list of hyperon daughters, get info about the decay
	for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){
		//geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
		if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

		art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

		SimParticle Decay = MakeSimParticle(part);
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

			SimParticle P = MakeSimParticle(part);			
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

						SimParticle P2 = MakeSimParticle(part2);
						P2.Origin = getOrigin(part2->TrackId());

						fDecay.push_back( P );

						daughter_IDs.push_back( part2->TrackId() );

						//add its ID to list od decay products for searching later



					}

				}



			}

			else std::cout << "Unrecognized Sigma0 daughter: " << part->PdgCode() << std::endl; 


		}

	}




	///////////////////////////////////////////////////////////////////////////
	//Get Reconstructed Info
	//////////////////////////////////////////////////////////////////////////



	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();


	fRecoTrackTruthMatchedID.clear(); //ids of truth matched MC particles

	//setup handles
	art::Handle<std::vector<recob::PFParticle>>pfparticleHandle; //PFParticles
	art::Handle<std::vector<recob::Track>  > trackHandle; //reconstructed tracks
	art::Handle<std::vector<recob::Shower> > showerHandle; //reconstructed showers
	art::Handle<std::vector<anab::ParticleID>>pidHandle;	
	art::Handle<std::vector<anab::Calorimetry>>caloHandle;
	art::Handle<std::vector<recob::Hit> > hitHandle;
	art::Handle<larpandoraobj::PFParticleMetadata> particleMetadataHandle;

	std::vector <art::Ptr <recob::Track>  > trackVect; //vector of tracks
	std::vector <art::Ptr <recob::Shower> > showerVect;//vector of showers
	std::vector< art::Ptr<recob::PFParticle>>pfparticleVect; //vector of PFparticles
	std::vector < art::Ptr<anab::ParticleID>> pidVect; //pids
	std::vector<art::Ptr<anab::Calorimetry>>caloVect; //calorimetry
	std::vector< art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVect; //metadata
	std::vector <art::Ptr <recob::Hit>  > hitVect; //hits

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
				//		fRecoPrimaryVertex.SetXYZ( vtx->position().X() , vtx->position().Y() , vtx->position().Z() );

				//w SC correction - forward
				fRecoPrimaryVertex.SetXYZ( vtx->position().X() - sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z() );


			}

		}

	}



	//go through rest of particles, get lots of useful info!
	for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){


		RecoParticle ThisPrimaryDaughter;



		//create PID store
		std::vector<int> thisTrackAlgNumbers;
		std::vector<int> thisTrackPlaneNos;
		std::vector<int> thisTrackAssumedPDGs;
		std::vector<double> thisTrackScores;
		double thisTrackTruePDG;
		double thisTrackRecoLength;
		double thisTrackTrueLength;

		thisTrackAlgNumbers.clear();
		thisTrackPlaneNos.clear();
		thisTrackAssumedPDGs.clear();
		thisTrackScores.clear();




		//get data from every PFP, not just neutrino daughters
		if(pfp->Parent() != neutrinoID) continue;

		std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
		std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());

		std::vector< art::Ptr<recob::Shower> > pfpShowers = showerAssoc.at(pfp.key());
		std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> >pfpMeta = metadataAssoc.at(pfp.key());

		ThisPrimaryDaughter.PDG = pfp->PdgCode();

		if(pfp->PdgCode() == 11) fNPrimaryShowerDaughters++;
		if(pfp->PdgCode() == 13) fNPrimaryTrackDaughters++;

		//	std::cout << "Tracks assoc to this pfp: " << pfpTracks.size() << "  Showers assoc to this pfp: " << pfpShowers.size() << std::endl;



		for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){

			const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());

			if (!pfParticlePropertiesMap.empty()){
				//			std::cout << " Found PFParticle " << pfp->Self() << " with: " << std::endl;
				for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){

					//                    std::cout << "  - " << it->first << " = " << it->second << std::endl;

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
		if(!pfpTracks.empty() && !pfpVertex.empty()){

			for(const art::Ptr<recob::Track> &trk : pfpTracks){

				//reconstructed track length
				ThisPrimaryDaughter.TrackLength = trk->Length();			
				thisTrackRecoLength = trk->Length();

				/////////////////////
				//   Truth match   //
				/////////////////////

			
				//loop through hits in track, find which MC particle
				//deposited most energy in track

				//get hits assoc with track
				std::vector< art::Ptr< recob::Hit> > hits = trackHitAssoc.at(trk.key());

				//use backtracker to find true pdg 

				std::unordered_map<int , double>  trkide;
				double maxe = -1;
				double tote = 0;

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

						trkide[ particleVec[i_particle]->TrackId() ] += matchVec[i_particle]->energy;

						tote += matchVec[i_particle]->energy;

						if ( trkide[ particleVec[i_particle]->TrackId() ] > maxe ){

							maxe = trkide[ particleVec[i_particle]->TrackId() ];
							matchedParticle = particleVec[i_particle];

						}

					}

				}



				if(matchedParticle != NULL){ 

					ThisPrimaryDaughter.HasTruth = true;
					ThisPrimaryDaughter.TrackEdepPurity = maxe/tote;

					ThisPrimaryDaughter.TrackTruePDG = matchedParticle->PdgCode();
					ThisPrimaryDaughter.TrackTrueE = matchedParticle->Momentum().E();
					ThisPrimaryDaughter.TrackTruePx = matchedParticle->Momentum().X();
					ThisPrimaryDaughter.TrackTruePy = matchedParticle->Momentum().Y();
					ThisPrimaryDaughter.TrackTruePz = matchedParticle->Momentum().Z();


					ThisPrimaryDaughter.TrackTrueModMomentum = sqrt( matchedParticle->E()*matchedParticle->E() - matchedParticle->Mass()*matchedParticle->Mass() );
					ThisPrimaryDaughter.TrackTrueKE = matchedParticle->E() - matchedParticle->Mass();

					ThisPrimaryDaughter.TrackTrueLength = sqrt( (matchedParticle->Vx() - matchedParticle->EndX())*(matchedParticle->Vx() - matchedParticle->EndX())
							+ (matchedParticle->Vy() - matchedParticle->EndY())*(matchedParticle->Vy() - matchedParticle->EndY())
							+ (matchedParticle->Vz() - matchedParticle->EndZ())*(matchedParticle->Vz() - matchedParticle->EndZ()) ); 

					ThisPrimaryDaughter.TrackTrueOrigin = getOrigin(matchedParticle->TrackId());

					thisTrackTruePDG = matchedParticle->PdgCode();
					thisTrackTrueLength = ThisPrimaryDaughter.TrackTrueLength;


				}
				else ThisPrimaryDaughter.HasTruth = false;





				///////////////////////
				// get PID for track //
				///////////////////////


				std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());
				std::vector<art::Ptr<anab::ParticleID>> trackPID = PIDAssoc.at(trk.key());

				std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

				for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){

					anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

					if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward) continue;

					int algno=-1;
					if(AlgScore.fAlgName == "Chi2") algno = 1;
					else if(AlgScore.fAlgName == "BraggPeakLLH") algno = 2;
					else if(AlgScore.fAlgName == "BraggPeakLLH_shift") algno = 3;
					else if(AlgScore.fAlgName == "PIDA_median") algno = 4;
					else if(AlgScore.fAlgName == "PIDA_mean") algno = 5;
					else if(AlgScore.fAlgName == "TruncatedMean") algno = 6;
					else if(AlgScore.fAlgName == "DepEvsRangeE") algno = 7;
					else if(AlgScore.fAlgName == "ThreePlaneProtonPID") algno = 8;
					else std::cout << "Unrecognized algorithm name: " << AlgScore.fAlgName << std::endl;

					thisTrackAlgNumbers.push_back(algno);


					thisTrackPlaneNos.push_back(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask));
					thisTrackAssumedPDGs.push_back(TMath::Abs(AlgScore.fAssumedPdg));
					thisTrackScores.push_back(AlgScore.fValue);



					//Just use 3 plane proton PID for the time being
					if(  TMath::Abs(AlgScore.fAssumedPdg) == 2212  &&   AlgScore.fAlgName=="ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){	

						ThisPrimaryDaughter.TrackPID = AlgScore.fValue;

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
				//TVector3 pos(vtx->position().X(),vtx->position().Y(),vtx->position().Z());
				//w SC correction - forward
				TVector3 pos( vtx->position().X() - sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z() );

				ThisPrimaryDaughter.SetVertex( pos );


				double d = sqrt(  (pos.X() - fRecoPrimaryVertex.X())*(pos.X() - fRecoPrimaryVertex.X())
						+ (pos.Y() - fRecoPrimaryVertex.Y())*(pos.Y() - fRecoPrimaryVertex.Y())
						+ (pos.Z() - fRecoPrimaryVertex.Z())*(pos.Z() - fRecoPrimaryVertex.Z()) );

				ThisPrimaryDaughter.Displacement = d;

			}

		}


		if(ThisPrimaryDaughter.PDG == 13){

			fTrackPrimaryDaughters.push_back( ThisPrimaryDaughter );

			//fill PID info

			//for(size_t j=0;j<thisTrackAlgNames.size();j++) std::cout << thisTrackAlgNames.at(j) << std::endl;

			fAlgNumber.push_back(thisTrackAlgNumbers); //pid scores for tracks 
			fPlaneNo.push_back(thisTrackPlaneNos); //pid scores for tracks 
			fAssumedPDG.push_back(thisTrackAssumedPDGs); //pid scores for tracks 
			fScore.push_back(thisTrackScores); //pid scores for tracks 
			fRecoLength.push_back(thisTrackRecoLength); //pid scores for tracks 
			fTruePDG.push_back(thisTrackTruePDG); //pid scores for tracks 
			fTrueLength.push_back(thisTrackTrueLength); //pid scores for tracks 


		}
		else fShowerPrimaryDaughters.push_back( ThisPrimaryDaughter );

	}//end of PFP loop
	

	if(fPrint) PrintInfo();


	fSelectedEvent = PerformSelection();

	FinishEvent();


}



///////////////////////////////////////////
// Print some useful info from the event //
///////////////////////////////////////////

void hyperon::HyperonSelection::PrintInfo(){


	if(!fHyperon.size()) return;

	if(!(fHyperon.at(0).PDG == fPrintPdg || fPrintPdg == -1) || !fInActiveTPC) return;

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

		if(primary_IDs.at(i_p)	== idnum){ 
			return 1;
		}

	}

	//search list of hyperon decay products	
	for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){

		if(daughter_IDs.at(i_d)	== idnum){ 
			return 2;
		}
	}


	return 3;

}





///////////////////////////////////////////////////////////////	
// Finished processing event - update Metadata and fill tree //
///////////////////////////////////////////////////////////////

void hyperon::HyperonSelection::FinishEvent(){

	//	if(ccnc == 0) {fNChargedCurrent++; fCCNC="CC"; }
	//	else { fNNeutralCurrent++; fCCNC="NC"; }

	
	
	if( fLepton.size() != 1 || isLepton(fLepton.at(0).PDG) ) { fNChargedCurrent++; fCCNC="CC"; }
	else { fNNeutralCurrent++; fCCNC="NC"; }

//	if( isNeutrino(fLepton.at(0).PDG) ) { fNNeutralCurrent++; fCCNC="NC"; }
//	else { fNChargedCurrent++; fCCNC="CC"; }

	if(fNuPDG == 12) fNnue++;
	else if(fNuPDG == 14) fNnuMu++;
	else if(fNuPDG == -12) fNnueBar++;
	else fNnuMuBar++;


	if(fInActiveTPC && fIsLambdaCharged && fNuPDG == -14) fIsSignal = true;

	//store info for this event
	fOutputTree->Fill();
	
	//store PID info
	fPIDTree->Fill();
	
	//update metadata

	fNEvents++; //total events in sample
	if(fIsHyperon && fInActiveTPC) fNHyperons++; //total hyperon events in active vol
	if(fSelectedEvent) fNSelectedEvents++; //total events passing selection
	if(fSelectedEvent && fInActiveTPC && fIsHyperon) fNSelectedHyperons++; //total hyperons in active vol passing preselection
	if(fSelectedEvent && fIsHyperon) fNSelectedHyperons++; //total true hyperon events passing selection


	//signal events
	if(fIsSignal) fNSignal++;
	if(fIsSignal && fSelectedEvent) fNSelectedSignal++;


	//efficiency = fNSelectedHyperonsInActiveVol/fNHyperonsInActiv;
	//purity = fNSelectedHyperonsInActiveVol/fSelectedEvents;


}




void hyperon::HyperonSelection::beginJob(){

	fileID=0;	


	// Implementation of optional member functions

	art::ServiceHandle<art::TFileService> tfs;

	fOutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");


	//write a second tree to store decay product info
	fOutputTree->Branch("EventID",&fEventID);
	fOutputTree->Branch("run",&run);
	fOutputTree->Branch("subrun",&subrun);
	fOutputTree->Branch("event",&event);
	fOutputTree->Branch("fileID",&fileID);

	fOutputTree->Branch("Mode",&fMode);
	fOutputTree->Branch("CCNC",&fCCNC);

	//generator truth info

	fOutputTree->Branch("InActiveTPC",&fInActiveTPC);
	fOutputTree->Branch("IsHyperon",&fIsHyperon);
	fOutputTree->Branch("IsSigmaZero",&fIsSigmaZero);
	fOutputTree->Branch("IsLambda",&fIsLambda);
	fOutputTree->Branch("IsLambdaCharged",&fIsLambdaCharged);
	fOutputTree->Branch("IsSignal",&fIsSignal);

	fOutputTree->Branch("SelectedEvent",&fSelectedEvent);

	//neutrino information
	fOutputTree->Branch("NuEnergy",&fNuEnergy); //neutrino energy
	fOutputTree->Branch("NuPDG",&fNuPDG); //neutrino PDG code
	fOutputTree->Branch("TruePrimaryVertex","TVector3",&fTruePrimaryVertex);

	//lepton
	fOutputTree->Branch("Lepton","vector<SimParticle>",&fLepton);
	//hyperon
	fOutputTree->Branch("Hyperon","vector<SimParticle>",&fHyperon);

	//nucleons and pions at the primary vtx
	fOutputTree->Branch("PrimaryNucleon","vector<SimParticle>",&fPrimaryNucleon);
	fOutputTree->Branch("PrimaryPion","vector<SimParticle>",&fPrimaryPion);

	//g4 truth info

	fOutputTree->Branch("DecayVertex","TVector3",&fDecayVertex); //position of hyperon decay vertex
	fOutputTree->Branch("Decay","vector<SimParticle>",&fDecay);

	fOutputTree->Branch("SigmaZeroDecayPhoton","vector<SimParticle>",&fSigmaZeroDecayPhoton);
	fOutputTree->Branch("SigmaZeroDecayLambda","vector<SimParticle>",&fSigmaZeroDecayLambda);


	//angles between particles involved

	fOutputTree->Branch("DecayOpeningAngle",&fDecayOpeningAngle); //opening angle between hyperon decay products     
	fOutputTree->Branch("LeptonPionAngle",&fLeptonPionAngle); //openining angle between lepton and pion
	fOutputTree->Branch("LeptonNucleonAngle",&fLeptonNucleonAngle); //opening angle between lepton and nucleon

	//////////////////////
	// Reco Information //
	//////////////////////

	fOutputTree->Branch("RecoPrimaryVertex","TVector3",&fRecoPrimaryVertex);
	fOutputTree->Branch("NPrimaryTrackDaughters",&fNPrimaryTrackDaughters); //num ofOutputTree->Branch("",f track like primary daughters
	fOutputTree->Branch("NPrimaryShowerDaughters",&fNPrimaryShowerDaughters); //num ofOutputTree->Branch("",f shower like primary daughters

	fOutputTree->Branch("TracklikePrimaryDaughters","vector<RecoParticle>",&fTrackPrimaryDaughters);
	fOutputTree->Branch("ShowerlikePrimaryDaughters","vector<RecoParticle>",&fShowerPrimaryDaughters);
	fOutputTree->Branch("MuonIndex",&fMuonIndex);


	///////////////////////
        //     PID Scores    //
        ///////////////////////     

        fPIDTree=tfs->make<TTree>("PIDTree","PID Scores Tree");
        fPIDTree->Branch("SelectedEvent",&fSelectedEvent);
        fPIDTree->Branch("IsSignal",&fIsSignal);
        fPIDTree->Branch("AlgNumber",&fAlgNumber); //pid scores for tracks 
        fPIDTree->Branch("PlaneNo",&fPlaneNo); //pid scores for tracks 
        fPIDTree->Branch("AssumedPDG",&fAssumedPDG); //pid scores for tracks 
        fPIDTree->Branch("Score",&fScore); //pid scores for tracks 
        fPIDTree->Branch("RecoLength",&fRecoLength); //pid scores for tracks 
        fPIDTree->Branch("TruePDG",&fTruePDG); //pid scores for tracks 
        fPIDTree->Branch("trueLength",&fRecoLength); //pid scores for tracks 





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


	fNChargedCurrent=0; //number of cc events
	fNNeutralCurrent=0; //number of nc events
	fNnuMu=0; //number of numu events
	fNnue=0; //number of nue events
	fNnuMuBar=0; //number of numubar events
	fNnueBar=0; //number of nuebar events


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

}



void hyperon::HyperonSelection::endJob()
{
	//calculate efficiency and purity

	if(fNSelectedEvents > 0){
		fHyperonPurity = (double)fNSelectedHyperons/fNSelectedEvents;
		fSignalPurity = (double)fNSelectedSignal/fNSelectedEvents;	
	}

	if(fNHyperons > 0) fHyperonEfficiency = (double)fNSelectedHyperons/fNHyperons;

	if(fNSelectedSignal > 0) fSignalEfficiency = (double)fNSelectedSignal/fNSignal;


	fHyperonEfficiencyTimesPurity = fHyperonPurity * fHyperonEfficiency;
	fHyperonTruePurity = 1/( 1+(1/fHyperonPurity-1)*30 );
	fHyperonEfficiencyTimesTruePurity = fHyperonEfficiency * fHyperonTruePurity;

	fSignalEfficiencyTimesPurity = fSignalPurity * fSignalEfficiency;
	fSignalTruePurity = 1/( 1+(1/fSignalPurity-1)*30 );
	fSignalEfficiencyTimesTruePurity = fSignalEfficiency * fSignalTruePurity;


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

	std::cout << std::endl << std::endl << std::endl;

	fMetaTree->Fill();

}



void hyperon::HyperonSelection::beginSubRun(const art::SubRun& sr)
{

	std::cout << "New file, ID: " << fileID << std::endl;
	fileID++;

}


void hyperon::HyperonSelection::endSubRun(const art::SubRun& sr){}


DEFINE_ART_MODULE(hyperon::HyperonSelection)
