////////////////////////////////////////////////////////////////////////
// Class:       HyperonPID
// Plugin Type: analyzer (art v3_03_01)
// File:        HyperonPID_module.cc
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
#include "TString.h"

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

//#include "larevt/SpaceCharge/SpaceChargeStandard.h"
//#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

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
#include "ubana/HyperonProduction/util/PIDStore.h"

#include "ubana/HyperonProduction/Alg/ParticleTypes.h"
#include "ubana/HyperonProduction/Alg/FV.h"
#include "ubana/HyperonProduction/Alg/Muon_ID.h"
#include "ubana/HyperonProduction/Alg/Track_Length_Cut.h"
#include "ubana/HyperonProduction/Alg/Gap_Cut.h"

#include "ubana/HyperonProduction/util/Helpers.h"

#include "TRandom2.h"




namespace hyperon {
	class HyperonPID;
}


class hyperon::HyperonPID : public art::EDAnalyzer {
	public:
		explicit HyperonPID(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		HyperonPID(HyperonPID const&) = delete;
		HyperonPID(HyperonPID&&) = delete;
		HyperonPID& operator=(HyperonPID const&) = delete;
		HyperonPID& operator=(HyperonPID&&) = delete;

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



		TRandom2 *R = new TRandom2();


		//takes information retrieved from event in analyze and 
		//applies selection criteria

		bool PerformSelection();


		//basic event info
		unsigned int fEventID;
		//run/subrun/event information
		int run,subrun,event;


		int fileID;

		//output trees
		TTree * fOutputTree;

		TTree * fPIDTree;

		TTree * fMetaTree;


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

		//lepton information				
		std::vector<SimParticle> fLepton; //Lepton produced
		
		//hyperon information
		std::vector<SimParticle> fHyperon;
		
		//nucleons produced at primary vtx
		std::vector<SimParticle> fPrimaryNucleon;

		//pions produced at primary vtx
		std::vector<SimParticle> fPrimaryPion;

		//vertex information
		TVector3 fTruePrimaryVertex;

		//neutrino information
		double fNuEnergy; //neutrino energy
		int fNuPDG; //neutrino PDG code

		//g4 truth info
		TVector3 fDecayVertex;

		std::vector<SimParticle> fDecay; //hyperon decay products


		std::vector<SimParticle> fSigmaZeroDecayPhoton;
		std::vector<SimParticle> fSigmaZeroDecayLambda;	

                double fDecayOpeningAngle; //opening angle between hyperon decay products
                double fLeptonPionAngle; //openining angle between lepton and pion
                double fLeptonNucleonAngle; //opening angle between lepton and nucleon

		int fDecayHyperonID;
		int fDecaySigmaID;



	//data storage (should not be written to output)

	//used by G4 to track particles
	std::vector<int>daughter_IDs; //ids of semistable hyperon decay products
	std::vector<int>Sigma0_daughter_IDs; //ids of sigma0 decay products
	std::vector<int>primary_IDs; //ids of particles produced at primary vertex

	//create map between particles and their ID's
	std::map< int , art::Ptr<simb::MCParticle> > partByID;
	std::pair< int ,art::Ptr<simb::MCParticle>  >  part_and_id;

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
		int fMuonIndex=-1; //index of muon candidate in fTrackPrimaryDaughters 

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


		int fTotalEvents; //total events in sample
		int fTruthEventsInActiveVol; //total events in active volume

		int fNChargedCurrent; //number of cc events
		int fNNeutralCurrent; //number of nc events
		int fNnuMu; //number of numu events
		int fNnue; //number of nue events
		int fNnuMuBar; //number of numubar events
		int fNnueBar; //number of nuebar events

		int fTotalHyperons;
		int fTotalHyperonsInActiveVol;

		int fTotalLambdas;
		int fTotalLambdasInActiveVol;

		int fTotalChargedLambdas;
		int fTotalChargedLambdasInActiveVol;		

		int fNSelectedEvents; //total events passing selection
		int fNSelectedHyperons; //total true hyperon events passing selection
		int fNSelectedHyperonsInActiveVol; //total true hyperon events with true primary vtx in active tpc that pass selection 
	
		int fNSelectedLambdas;
		int fNSelectedLambdasInActiveVol;

		int fNSelectedChargedLambdas;
		int fNSelectedChargedLambdasInActiveVol;

		int fTotalSignal;
		int fSelectedSignal;


		double fHyperonEfficiency; //hyperon selection efficiency
		double fHyperonPurity;  //hyperon selection purity
		double fHyperonTruePurity; //hyperon selection efficiency x purity
		double	fHyperonEfficiencyTimesPurity; //hyperon selection efficiency x purity
		double	fHyperonEfficiencyTimesTruePurity; //hyperon selection efficiency x true purity

		double fLambdaEfficiency; //hyperon selection efficiency
		double fLambdaPurity;  //hyperon selection purity
		double fLambdaTruePurity; //hyperon selection efficiency x purity
		double	fLambdaEfficiencyTimesPurity; //hyperon selection efficiency x purity
		double	fLambdaEfficiencyTimesTruePurity; //hyperon selection efficiency x true purity

		double	fChargedLambdaEfficiency=0; //hyperon selection efficiency
		double	fChargedLambdaPurity=0;  //hyperon selection purity
		double	fChargedLambdaTruePurity=0; //hyperon purity after converting from enriched sample to real sample
		double	fChargedLambdaEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
		double	fChargedLambdaEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

		double	fSignalEfficiency=0; //hyperon selection efficiency
		double	fSignalPurity=0;  //hyperon selection purity
		double	fSignalTruePurity=0; //hyperon purity after converting from enriched sample to real sample
		double	fSignalEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
		double	fSignalEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

		//data storage (not to be written to output)		

		std::vector<int> fRecoTrackTruthMatchedID; //list of ID numbers of truth matched particles

		///////////////////////
		// Module Labels     //
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
		bool fPrintChargedLambda;


		///////////////////////////////////
		// Hyperon selection parameters  //
		///////////////////////////////////

		int fMinDaughters;//minimum number of daughters
		int fMaxPrimaryShowerDaughters;
		int fMinDaughterTracks; //minimum number of track like primary daughters
		double fMuonPIDScore;
		double fLeadingTrackLengthCut;
		double fSubleadingTrackLengthCut;
		double fMaxGap;	
		double fMinGap;	

		double fTrackScore; //minimum track/shower score required for pfp to be classified as track

			

};


////////////////////////////////////////////////////
// Decides if an event passes selection criteria  //
////////////////////////////////////////////////////

bool hyperon::HyperonPID::PerformSelection(){


if(fRecoPrimaryVertex.X() == -1000 || fRecoPrimaryVertex.Y() == -1000 || fRecoPrimaryVertex.Y() == -1000) return false; //no reco'd neutrino 

/*
if(fIsLambdaCharged && fInActiveTPC) {
std::cout << std::endl;
std::cout << "EventID: " << fEventID << std::endl;
std::cout << "Num Shower Daughters: " <<  fNPrimaryShowerDaughters << std::endl;
std::cout << "Num Track Daughters: " << fNPrimaryTrackDaughters << std::endl; 
std::cout << std::endl;
}
*/

//Apply fiducial volume cut

	if(!inActiveTPC(fRecoPrimaryVertex)) return false;

	if(fNPrimaryDaughters < fMinDaughters) return false;
//	std::cout << "Passed Min Daughter cut" << std::endl;
	if(fNPrimaryShowerDaughters > fMaxPrimaryShowerDaughters) return false;
//	std::cout << "Pass Max Shower cut" << std::endl;
	if(fNPrimaryTrackDaughters < fMinDaughterTracks) return false;
//	std::cout << "Passed Min Track cut" << std::endl;



//add Muon ID
int i_muon = Muon_ID(fTrackPrimaryDaughters,fMuonPIDScore);
fMuonIndex = i_muon;
if(i_muon == -1) return false;

//if rejected by track length cut, reject event
if(!Track_Length_Cut(fTrackPrimaryDaughters,i_muon,fLeadingTrackLengthCut,fSubleadingTrackLengthCut)) return false;

//add gap cut
//if(!Gap_Cut(fTrackPrimaryDaughters,i_muon,fMinGap,fMaxGap)) return false;

std::cout << "Event Num: " << event << "   ";


if(fIsLambdaCharged && fInActiveTPC) std::cout << "__ACCEPTED_SIGNAL__"  << "  ";
else std::cout << "__ACCEPTED_BACKGROUND__" << "  ";   	

std::cout << fMode << std::endl; 

		return true;


}





hyperon::HyperonPID::HyperonPID(fhicl::ParameterSet const& p)
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
	fLeadingTrackLengthCut = p.get<double>("LeadingTrackLengthCut",65);
	fSubleadingTrackLengthCut = p.get<double>("SubleadingTrackLengthCut",35);

	fMaxGap = p.get<double>("MaxGap",80);	
	fMinGap = p.get<double>("MinGap",0);	

	
	fTrackScore = p.get<double>("TrackScore",0.5);




	//misc parameters
	fPrint = p.get<bool>("Print",false);

	fPrintPdg = p.get<int>("PrintPdg",-1);
	fPrintChargedLambda = p.get<bool>("PrintChargedLambda",false);

}

void hyperon::HyperonPID::analyze(art::Event const& e)
{


	//RESET EVERYTHING

	//generator truth info

	fInActiveTPC=false;
	fIsHyperon=false; //true if event contains a hyperon
	fIsSigmaZero=false; //true if event hyperon was a sigma zero		
	fIsLambda=false;
	fIsLambdaCharged=false;
	fIsSignal=false;	

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


	//RECO LEVEL INFO

	fSelectedEvent = false; //true if event passes some selection criteria	
	fNPrimaryDaughters = 0; //number of primary daughters

	fNPrimaryTrackDaughters=0; //num of track like primary daughters
	fNPrimaryShowerDaughters=0; //num of shower like primary daughters

	fRecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

	fTrackPrimaryDaughters.clear();
	fShowerPrimaryDaughters.clear();
        fMuonIndex=-1;
	


	//TRACK PID INFO
	fAlgNumber.clear(); //pid scores for tracks 
	fPlaneNo.clear(); //pid scores for tracks 
	fAssumedPDG.clear(); //pid scores for tracks 
	fScore.clear(); //pid scores for tracks 
	fRecoLength.clear(); //pid scores for tracks 
	fTruePDG.clear(); //pid scores for tracks 
	fTrueLength.clear(); //pid scores for tracks 





	fEventID = e.id().event();

	run = e.run();
	subrun = e.subRun();
	event = e.event();

 
	// GET EVENT GENERATOR INFO

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

			//get antimuons and positrons in the sample
			if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1){ 


				fTruePrimaryVertex.SetXYZ( Part.Vx() , Part.Vy() , Part.Vz() );

				fInActiveTPC=inActiveTPC(fTruePrimaryVertex);



			}//if lepton


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


	//	std::cout << "List of G4 Particles in this event: " << std::endl;

	daughter_IDs.clear();
	primary_IDs.clear();
	Sigma0_daughter_IDs.clear();
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

		//geant does not always keep everything it simulates, make sure particle is in list of IDs

		if(partByID.find(primary_IDs[i_p]) == partByID.end()) continue;

		art::Ptr<simb::MCParticle> part = partByID[primary_IDs[i_p]];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these


		//hyperon produced at primary vertex
		if(isHyperon(part->PdgCode())) {

			fHyperon.push_back( MakeSimParticle(part) );
			fHyperon.back().Origin = getOrigin(part->TrackId());	
			fIsHyperon = true;	
					}

		//lepton produced at primary vertex
		if(isLepton(part->PdgCode()) || isNeutrino(part->PdgCode())){

			SimParticle P = MakeSimParticle(part);
			P.Origin = getOrigin(part->TrackId());	
			
			fLepton.push_back( P );
//			fLepton.back().Origin = getOrigin(part->TrackId());	
			
			}
		

		if(isNucleon(part->PdgCode()) ){
	
                        fPrimaryNucleon.push_back( MakeSimParticle(part) );
			fPrimaryNucleon.back().Origin = getOrigin(part->TrackId());	
		}

		if(isPion(part->PdgCode())){
                        		
                        fPrimaryPion.push_back( MakeSimParticle(part) );
			fPrimaryPion.back().Origin = getOrigin(part->TrackId());	
	
		}

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


			std::vector<double> p = {part->E() , part->Px(),part->Py(),part->Pz()  };

			//one particle should be a Lambda, one should be a photon

			//if is photon
			if(part->PdgCode() == 22){

                        fSigmaZeroDecayPhoton.push_back( MakeSimParticle(part) );

				//add its ID to list of primary ID's - produced at primary vtx so from reco point of view makes sense to put it in
				//list of primary IDs

				primary_IDs.push_back(part->TrackId());

			}
			//Lambda produced from sigma zero decay
			else if(part->PdgCode() == 3122){

				//now search for Lambda decay products, find its decay vertex etc			

                        fSigmaZeroDecayLambda.push_back( MakeSimParticle(part) );



				if(part->EndProcess() == "Decay"){		

					//get Lambdas decay vertex

					fDecayVertex.SetXYZ( part->EndPosition().X() , part->EndPosition().Y() , part->EndPosition().Z() );


					for(int i_d2=0;i_d2<part->NumberDaughters();i_d2++){

						//search the map for the daughter IDs
						if(partByID.find(part->Daughter(i_d2)) == partByID.end()) continue;

						art::Ptr<simb::MCParticle> part2 = partByID[part->Daughter(i_d2)];

						if(part2->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
	
				                fDecay.push_back( MakeSimParticle(part2) );


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

//try to get SC corrections for some random points




fRecoTrackTruthMatchedID.clear();
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
std::vector < art::Ptr<anab::ParticleID>> pidVect;
std::vector<art::Ptr<anab::Calorimetry>>caloVect;
std::vector< art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVect;
std::vector <art::Ptr <recob::Hit>  > hitVect; 
//Fill PFP and Track vectors


if(e.getByLabel(fPFParticleLabel,pfparticleHandle)){
	art::fill_ptr_vector(pfparticleVect,pfparticleHandle);
}

if(!pfparticleVect.size()) {
	FinishEvent();
	return;
}



if(e.getByLabel(fTrackLabel,trackHandle)) art::fill_ptr_vector(trackVect,trackHandle);
else
std::cout << "Track handle not setup" << std::endl;

if(e.getByLabel(fShowerLabel,showerHandle)) art::fill_ptr_vector(showerVect,showerHandle);
else
std::cout << "Shower handle not setup" << std::endl;




e.getByLabel(fHitLabel,hitHandle);

//	if(e.getByLabel(fHitLabel,hitHandle)) art::fill_ptr_vector(hitVect,hitHandle);
//	else
//		std::cout << "Hit handle not setup" << std::endl;


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
		//w SC correction - backward
		//fRecoPrimaryVertex.SetXYZ( vtx->position().X() + sce_corr.X() , vtx->position().Y() + sce_corr.Y() , vtx->position().Z() + sce_corr.Z() );


		}




	}

}



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

//std::cout << pfp->PdgCode() << std::endl;
	
	ThisPrimaryDaughter.PDG = pfp->PdgCode();

	if(pfp->PdgCode() == 11) fNPrimaryShowerDaughters++;
	if(pfp->PdgCode() == 13) fNPrimaryTrackDaughters++;

//	std::cout << "Tracks assoc to this pfp: " << pfpTracks.size() << "  Showers assoc to this pfp: " << pfpShowers.size() << std::endl;

	//print out the track/shower scores


//////////////////////////////////////////////
// Get metadata - record track/shower score //
//////////////////////////////////////////////


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



	//if pfp is a track like object
	if(!pfpTracks.empty() && !pfpVertex.empty()){

		for(const art::Ptr<recob::Track> &trk : pfpTracks){

			ThisPrimaryDaughter.TrackLength = trk->Length();			
			thisTrackRecoLength = trk->Length();
			//start of truth match

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



			//get PID for track

		
			std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());
			std::vector<art::Ptr<anab::ParticleID>> trackPID = PIDAssoc.at(trk.key());

			std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

			//get Chi2 scores for each particle and each plane
			std::vector<double> Proton_Chi2_Vec;
			std::vector<double> Muon_Chi2_Vec;
			std::vector<double> Kaon_Chi2_Vec;
			std::vector<double> Pion_Chi2_Vec;

			for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){

			//setup empty PID Store

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
				
				//get muon, proton and pion PID scores for each particle			

				if(  TMath::Abs(AlgScore.fAssumedPdg) == 2212  &&   AlgScore.fAlgName=="ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){	

					ThisPrimaryDaughter.TrackPID = AlgScore.fValue;

				}





if( AlgScore.fAlgName=="Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward ){
//	int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);

//	if(TMath::Abs(AlgScore.fAssumedPdg) == 2212) Proton_Chi2_Vec.push_back(AlgScore.fValue);
//	if(TMath::Abs(AlgScore.fAssumedPdg) == 13) Muon_Chi2_Vec.push_back(AlgScore.fValue);
//	if(TMath::Abs(AlgScore.fAssumedPdg) == 321) Kaon_Chi2_Vec.push_back(AlgScore.fValue);
//	if(TMath::Abs(AlgScore.fAssumedPdg) == 211) Pion_Chi2_Vec.push_back(AlgScore.fValue);

if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 2){

//std::cout << "Plane=" << planeid << " Assumed PDG=" << TMath::Abs(AlgScore.fAssumedPdg) << " Score=" << AlgScore.fValue << std::endl; 

if(TMath::Abs(AlgScore.fAssumedPdg) == 2212) ThisPrimaryDaughter.TrackProtonChi2 = AlgScore.fValue;
if(TMath::Abs(AlgScore.fAssumedPdg) == 13) ThisPrimaryDaughter.TrackMuonChi2 = AlgScore.fValue;
if(TMath::Abs(AlgScore.fAssumedPdg) == 321) ThisPrimaryDaughter.TrackKaonChi2 = AlgScore.fValue;
if(TMath::Abs(AlgScore.fAssumedPdg) == 211) ThisPrimaryDaughter.TrackPionChi2 = AlgScore.fValue;
	
}

}





			}

			





		}



		for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

		geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };                
		geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

		//no SC correction
//		TVector3 pos(vtx->position().X(),vtx->position().Y(),vtx->position().Z());
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


}



if(fPrint) PrintInfo();

fSelectedEvent = PerformSelection();

FinishEvent();


}

// Print some useful info from the event


void hyperon::HyperonPID::PrintInfo(){


if(fPrintChargedLambda && !fIsLambdaCharged && fNuPDG != -14) return;

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


//check origin of particle by its TrackId

int hyperon::HyperonPID::getOrigin(int idnum){


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


//updates metadata information
void hyperon::HyperonPID::FinishEvent(){

	if(ccnc == 0) {fNChargedCurrent++; fCCNC="CC"; }
	else { fNNeutralCurrent++; fCCNC="NC"; }

	if(fNuPDG == 12) fNnue++;
	else if(fNuPDG == 14) fNnuMu++;
	else if(fNuPDG == -12) fNnueBar++;
	else fNnuMuBar++;


	if(fInActiveTPC && fIsLambdaCharged && fNuPDG == -14) fIsSignal = true;


	//store info for this event
	fOutputTree->Fill();
	fPIDTree->Fill();

	//update metadata

	fTotalEvents++; //total events in sample
	if(fInActiveTPC) fTruthEventsInActiveVol++; //total events in active volume
	if(fIsHyperon) fTotalHyperons++; //total hyperon events
	if(fIsHyperon && fInActiveTPC) fTotalHyperonsInActiveVol++; //total hyperons in active volume
	if(fSelectedEvent) fNSelectedEvents++; //total events passing selection
	if(fSelectedEvent && fInActiveTPC && fIsHyperon) fNSelectedHyperonsInActiveVol++;   //total events that pass selection that are true hyperons and true primary vtx in active vol
	if(fSelectedEvent && fIsHyperon) fNSelectedHyperons++; //total true hyperon events passing selection

	if(fIsLambda) fTotalLambdas++;
	if(fIsLambda && fInActiveTPC) fTotalLambdasInActiveVol++;
	if(fSelectedEvent && fIsLambda) fNSelectedLambdas++;
	if(fSelectedEvent && fIsLambda && fInActiveTPC) fNSelectedLambdasInActiveVol++;

	if(fIsLambdaCharged) fTotalChargedLambdas++;
	if(fIsLambdaCharged && fInActiveTPC) fTotalChargedLambdasInActiveVol++;
	if(fSelectedEvent && fIsLambdaCharged) fNSelectedChargedLambdas++;
	if(fSelectedEvent && fIsLambdaCharged && fInActiveTPC) fNSelectedChargedLambdasInActiveVol++;

	if(fIsSignal) fTotalSignal++;
	if(fIsSignal && fSelectedEvent) fSelectedSignal++;


	//efficiency = fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiv;
	//purity = fNSelectedHyperonsInActiveVol/fSelectedEvents;


}




void hyperon::HyperonPID::beginJob(){

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
	fOutputTree->Branch("IsHyperon",&fIsHyperon); //true if event contains a hyperon
	fOutputTree->Branch("IsSigmaZero",&fIsSigmaZero); //true if event hyperon was a sigma zero		
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


	//information about plane containing decay products and 
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
	//            metadata tree		//
	//////////////////////////////////////////

	fTotalEvents=0; //total events in sample
	fTruthEventsInActiveVol=0; //total events in active volume
	fTotalHyperons=0;
	fTotalHyperonsInActiveVol=0;

	fTotalLambdas=0;
	fTotalLambdasInActiveVol=0;

	 fTotalChargedLambdas=0;
	 fTotalChargedLambdasInActiveVol=0;



	fNSelectedEvents=0; //total events passing selection
	fNSelectedHyperons=0; //total true hyperon events passing selection
	fNSelectedHyperonsInActiveVol=0; //total events that pass selection that are true hyperons and true primary vtx in active vol
	fNSelectedLambdas=0; //total true hyperon events passing selection
	fNSelectedLambdasInActiveVol=0; //total events that pass selection that are true hyperons and true primary vtx in active vol


	fNSelectedChargedLambdas=0;
	fNSelectedChargedLambdasInActiveVol=0;
	
	fHyperonEfficiency=0; //hyperon selection efficiency
	fHyperonPurity=0;  //hyperon selection purity
	fHyperonTruePurity=0; //hyperon purity after converting from enriched sample to real sample
	fHyperonEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
	fHyperonEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity
	fLambdaEfficiency=0; //hyperon selection efficiency
	fLambdaPurity=0;  //hyperon selection purity
	fLambdaTruePurity=0; //hyperon purity after converting from enriched sample to real sample
	fLambdaEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
	fLambdaEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

	fChargedLambdaEfficiency=0; //hyperon selection efficiency
	fChargedLambdaPurity=0;  //hyperon selection purity
	fChargedLambdaTruePurity=0; //hyperon purity after converting from enriched sample to real sample
	fChargedLambdaEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
	fChargedLambdaEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

	
	fTotalSignal=0;
	fSelectedSignal=0;


	fNChargedCurrent=0; //number of cc events
	fNNeutralCurrent=0; //number of nc events
	fNnuMu=0; //number of numu events
	fNnue=0; //number of nue events
	fNnuMuBar=0; //number of numubar events
	fNnueBar=0; //number of nuebar events

	

	fMetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");

	fMetaTree->Branch("TotalEvents",&fTotalEvents);
	fMetaTree->Branch("TruthEventsInActiveVol",&fTruthEventsInActiveVol);

	fMetaTree->Branch("NChargedCurrent",&fNChargedCurrent); //number of cc events
	fMetaTree->Branch("NNeutralCurrent",& fNNeutralCurrent); //number of nc events
	fMetaTree->Branch("Nnumu",&fNnuMu); //number of numu events
	fMetaTree->Branch("Nnue",&fNnue); //number of nue events
	fMetaTree->Branch("NnuMuBar",&fNnuMuBar); //number of numubar events
	fMetaTree->Branch("NnueBar",&fNnueBar); //number of nuebar events




	fMetaTree->Branch("TotalHyperons",&fTotalHyperons);
	fMetaTree->Branch("TruthHyperonsInActiveVol",&fTotalHyperonsInActiveVol);
	fMetaTree->Branch("NSelectedHyperons",&fNSelectedHyperons);
	fMetaTree->Branch("NSelectedHyperonsInActiveVol",&fNSelectedHyperonsInActiveVol);
	fMetaTree->Branch("HyperonEfficiency",&fHyperonEfficiency);
	fMetaTree->Branch("HyperonPurity",&fHyperonPurity);
	fMetaTree->Branch("HyperonEfficiencyTimesPurity",&fHyperonEfficiencyTimesPurity);
	fMetaTree->Branch("HyperonTruePurity",&fHyperonTruePurity);
	fMetaTree->Branch("HyperonEfficiencyTimesTruePurity",&fHyperonEfficiencyTimesTruePurity);

	fMetaTree->Branch("TotalLambdas",&fTotalLambdas);
	fMetaTree->Branch("TruthLambdasInActiveVol",&fTotalLambdasInActiveVol);
	fMetaTree->Branch("NSelectedLambdas",&fNSelectedLambdas);
	fMetaTree->Branch("NSelectedLambdasInActiveVol",&fNSelectedLambdasInActiveVol);
	fMetaTree->Branch("LambdaEfficiency",&fLambdaEfficiency);
	fMetaTree->Branch("LambdaPurity",&fLambdaPurity);
	fMetaTree->Branch("LambdaEfficiencyTimesPurity",&fLambdaEfficiencyTimesPurity);
	fMetaTree->Branch("LambdaTruePurity",&fLambdaTruePurity);
	fMetaTree->Branch("LambdaEfficiencyTimesTruePurity",&fLambdaEfficiencyTimesTruePurity);

	fMetaTree->Branch("TotalChargedLambdas",&fTotalChargedLambdas);
	fMetaTree->Branch("TruthChargedLambdasInActiveVol",&fTotalChargedLambdasInActiveVol);
	fMetaTree->Branch("NSelectedChargedLambdas",&fNSelectedChargedLambdas);
	fMetaTree->Branch("NSelectedChargedLambdasInActiveVol",&fNSelectedChargedLambdasInActiveVol);
	fMetaTree->Branch("ChargedLambdaEfficiency",&fChargedLambdaEfficiency);
	fMetaTree->Branch("ChargedLambdaPurity",&fChargedLambdaPurity);
	fMetaTree->Branch("ChargedLambdaEfficiencyTimesPurity",&fChargedLambdaEfficiencyTimesPurity);
	fMetaTree->Branch("ChargedLambdaTruePurity",&fChargedLambdaTruePurity);
	fMetaTree->Branch("ChargedLambdaEfficiencyTimesTruePurity",&fChargedLambdaEfficiencyTimesTruePurity);
	
	fMetaTree->Branch("NSelectedEvents",&fNSelectedEvents);
	fMetaTree->Branch("TotalSignal",&fTotalSignal);
	fMetaTree->Branch("NSelectedSignal",&fSelectedSignal);	
	fMetaTree->Branch("SignalEfficiency",&fSignalEfficiency);
	fMetaTree->Branch("SignalPurity",&fSignalPurity);
	fMetaTree->Branch("SignalEfficiencyTimesPurity",&fSignalEfficiencyTimesPurity);
	fMetaTree->Branch("SignalTruePurity",&fSignalTruePurity);
	fMetaTree->Branch("SignalEfficiencyTimesTruePurity",&fSignalEfficiencyTimesTruePurity);

}



void hyperon::HyperonPID::endJob()
{
	//calculate efficiency and purity

	//efficiency = fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiveVol;
	//purity = fNSelectedHyperonsInActiveVol/fSelectedEvents;

	if(fNSelectedEvents > 0){
 	fHyperonPurity = (double)fNSelectedHyperonsInActiveVol/fNSelectedEvents;
	fLambdaPurity = (double)fNSelectedLambdasInActiveVol/fNSelectedEvents;
	fChargedLambdaPurity = (double)fNSelectedChargedLambdasInActiveVol/fNSelectedEvents;
	fSignalPurity = (double)fSelectedSignal/fNSelectedEvents;	
	}

	if(fTotalHyperonsInActiveVol > 0) fHyperonEfficiency = (double)fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiveVol;

	if(fTotalLambdasInActiveVol > 0) fLambdaEfficiency = (double)fNSelectedLambdasInActiveVol/fTotalLambdasInActiveVol;

	if(fTotalChargedLambdasInActiveVol > 0) fChargedLambdaEfficiency = (double)fNSelectedChargedLambdasInActiveVol/fTotalChargedLambdasInActiveVol;

	if(fSelectedSignal > 0) fSignalEfficiency = (double)fSelectedSignal/fTotalSignal;




	fHyperonEfficiencyTimesPurity = fHyperonPurity * fHyperonEfficiency;
	fHyperonTruePurity = 1/( 1+(1/fHyperonPurity-1)*50 );
	fHyperonEfficiencyTimesTruePurity = fHyperonEfficiency * fHyperonTruePurity;

		
	fLambdaEfficiencyTimesPurity = fLambdaPurity * fLambdaEfficiency;
	fLambdaTruePurity = 1/( 1+(1/fLambdaPurity-1)*50 );
	fLambdaEfficiencyTimesTruePurity = fLambdaEfficiency * fLambdaTruePurity;


	fChargedLambdaEfficiencyTimesPurity = fChargedLambdaPurity * fChargedLambdaEfficiency;
	fChargedLambdaTruePurity = 1/( 1+(1/fChargedLambdaPurity-1)*50 );
	fChargedLambdaEfficiencyTimesTruePurity = fChargedLambdaEfficiency * fChargedLambdaTruePurity;

	fSignalEfficiencyTimesPurity = fSignalPurity * fSignalEfficiency;
	fSignalTruePurity = 1/( 1+(1/fSignalPurity-1)*50 );
	fSignalEfficiencyTimesTruePurity = fSignalEfficiency * fSignalTruePurity;





	std::cout << std::endl << std::endl<< std::endl << "Some Metadata for this sample:" << std::endl;

	std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
	std::cout << "Hyperons in active volume: " << fTotalHyperonsInActiveVol << std::endl;
	std::cout << "Selected Hyperons in active volume: " << fNSelectedHyperonsInActiveVol << std::endl;
	std::cout << "Enriched Sample Hyperon Efficiency: " << fHyperonEfficiency*100 << std::endl;
	std::cout << "Enriched Sample Hyperon Purity: " << fHyperonPurity*100 << std::endl;	
	std::cout << "Enriched Sample Hyperon E x P: " << fHyperonEfficiencyTimesPurity*100 << std::endl;
	std::cout << "True Hyperon Purity assuming enrichment factor of 50: " << fHyperonTruePurity*100 <<  std::endl;

	std::cout << std::endl;
	std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
	std::cout << "Lambdas in active volume: " << fTotalLambdasInActiveVol << std::endl;
	std::cout << "Selected Lambdas in active volume: " << fNSelectedLambdasInActiveVol << std::endl;	
	std::cout << "Enriched Sample Lambda Efficiency: " << fLambdaEfficiency*100 << std::endl;
	std::cout << "Enriched Sample Lambda Purity: " << fLambdaPurity*100 << std::endl;	
	std::cout << "Enriched Sample Lambda E x P: " << fLambdaEfficiencyTimesPurity*100 << std::endl;
	std::cout << "True Lambda Purity assuming enrichment factor of 50: " << fLambdaTruePurity*100 <<  std::endl;


	std::cout << std::endl;
	std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
	std::cout << "ChargedLambdas in active volume: " << fTotalChargedLambdasInActiveVol << std::endl;
	std::cout << "Selected ChargedLambdas in active volume: " << fNSelectedChargedLambdasInActiveVol << std::endl;	
	std::cout << "Enriched Sample ChargedLambda Efficiency: " << fChargedLambdaEfficiency*100 << std::endl;
	std::cout << "Enriched Sample ChargedLambda Purity: " << fChargedLambdaPurity*100 << std::endl;	
	std::cout << "Enriched Sample ChargedLambda E x P: " << fChargedLambdaEfficiencyTimesPurity*100 << std::endl;
	std::cout << "True ChargedLambda Purity assuming enrichment factor of 50: " << fChargedLambdaTruePurity*100 <<  std::endl;

	std::cout << std::endl;
	std::cout << "Selected Events: " <<  fNSelectedEvents << std::endl;
	std::cout << "Signal Events: " << fTotalSignal << std::endl;
	std::cout << "Selected Signal events in active volume: " << fSelectedSignal  << std::endl;	
	std::cout << "Enriched Sample Efficiency: " << fSignalEfficiency*100 << std::endl;
	std::cout << "Enriched Sample Purity: " << fSignalPurity*100 << std::endl;	
	std::cout << "Enriched Sample E x P: " << fSignalEfficiencyTimesPurity*100 << std::endl;
	std::cout << "True Purity assuming enrichment factor of 50: " << fSignalTruePurity*100 <<  std::endl;

	std::cout << std::endl << std::endl << std::endl;




	fMetaTree->Fill();

}



void hyperon::HyperonPID::beginSubRun(const art::SubRun& sr)
{

	std::cout << "New file, ID: " << fileID << std::endl;
	fileID++;

}


void hyperon::HyperonPID::endSubRun(const art::SubRun& sr){}


DEFINE_ART_MODULE(hyperon::HyperonPID)
