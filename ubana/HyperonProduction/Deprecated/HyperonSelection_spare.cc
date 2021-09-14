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

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TVector3.h"
#include "TLorentzVector.h"

//local includes
#include "Alg/ParticleTypes.h"
#include "Alg/Functions.h"


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

		std::string getOrigin(int idnum);
	


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
		TTree * fTruthTree;
		TTree * fRecoTree;
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


		//lepton information
		std::vector<double> fLeptonMomentum; //lepton component wise momentum			
		double fLeptonModMomentum; //magnitude of lepton momentum
		double fLeptonKE;
		double fLeptonTheta;
		double fLeptonPhi;
		int fLeptonPDG; //outgoing lepton PDG code				
		


		//hyperon inforamtion
		std::vector<double>fHyperonMomentum; //hyperon momentum component wise
		double fHyperonModMomentum; //magnitude of hyperon momentum
		double fHyperonKE; //hyperon kinetic energy
		double fHyperonTheta;
		double fHyperonPhi;
		int fHyperonPDG; //hyperon pdg code

		//info for other nucleons/pions produced at primary vertex
		std::vector<std::vector<double>>fPrimaryNucleonMomentum; //componentwise momentum of nucleons
		std::vector<double> fPrimaryNucleonModMomentum; //magnitude of nucleon momentum
		std::vector<double> fPrimaryNucleonKE;
		std::vector<int> fPrimaryNucleonPDG; //PDG codes of nucleons

		std::vector<std::vector<double>>fPrimaryPionMomentum; //componentwise momentum of nucleons
		std::vector<double> fPrimaryPionModMomentum; //magnitude of nucleon momentum
		std::vector<double> fPrimaryPionKE;
		std::vector<int> fPrimaryPionPDG; //PDG codes of nucleons


		//vertex information
		std::vector<double>fTruePrimaryVertex;

		//neutrino information
		double fNuEnergy; //neutrino energy
		int fNuPDG; //neutrino PDG code

		//g4 truth info
		double fLeptonTravel; //distance travelled by lepton produced at primary vertex
		std::vector<double>fDecayVertex; //position of hyperon decay vertex
		std::vector<double>fDecayPDG; //pdg codes of decay products
		std::vector<std::vector<double>>fDecayMomenta; //component wise momenta of decay products
		std::vector<double> fDecayModMomenta; //magnitude of momenta of hyperon decay products
		std::vector<double> fDecayKE;
		std::vector<double> fDecayTheta;
		std::vector<double> fDecayPhi;
		double fHyperonTravel; //distance travlled by hyperon
		std::vector<double>fDecayTravel; //distance travelled by hyperon decay products

		std::vector<double>fSigmaZeroDecayPhotonMomentum; //component wise momentum of photon produced during Sigma0 decay
		double fSigmaZeroDecayPhotonModMomentum; //magnitude of momentum of photon produced during Sigma0 decay
		std::vector<double>fSigmaZeroDecayLambdaMomentum; //
		double fSigmaZeroDecayLambdaModMomentum;
		double fSigmaZeroDecayLambdaKE;



                double fDecayOpeningAngle; //opening angle between hyperon decay products
                std::vector<double> fDecayPlaneVector; //vector normal to plane containing hyperon (cross product of pion and nucleon momenta)
                double fLeptonPionAngle; //openining angle between lepton and pion
                double fLeptonNucleonAngle; //opening angle between lepton and nucleon
		//angles of decay plane
		double fDecayPlaneTheta;
		double fDecayPlanePhi;



		int fDecayHyperonID;
		int fDecaySigmaID;

		bool fLambdaCharged; //true if event contains Lambda decaying into p + pi-


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

		std::vector<double>fRecoPrimaryVertex;

		int fNPrimaryDaughters; //num of primary daughters
		int fNPrimaryTrackDaughters; //num of track like primary daughters
		int fNPrimaryShowerDaughters; //num of shower like primary daughters

		// Track like primary daughters	

		std::vector<double>fDaughterTrackLengths;
		std::vector<double>fDaughterTrackPIDs; //3 plane proton PID scores for daughter tracks
		std::vector<std::vector<double>>fDaughterTrackVertices; //vertex positions of daughters
		std::vector<double>fDaughterVertexDisplacement;
		std::vector<double>fTrackTrackScores; //pandora track/shower score : <0.5 is shower, >0.5 is track


		// Shower like primary daughters
		std::vector<double>fShowerTrackScores; //pandora track/shower score : <0.5 is shower, >0.5 is track
		std::vector<std::vector<double>>fShowerDirections; //directions of shower objects - compare to momenta of mc particles

		//reco truth matched data

		//reconstructed tracks
		std::vector<int>fRecoTrackTruthMatchedPDG; //truth matched pdg codes of track objects
		std::vector<std::vector<double>>fRecoTrackTruthMatchedVertex;
		std::vector<std::vector<double>>fRecoTrackTruthMatchedMomentum; //momentum of matched particle
		std::vector<double>fRecoTrackTruthMatchedModMomentum; //magnitude of momentum of matched particle
		std::vector<double>fRecoTrackTruthMatchedKE; //KE of matched particle
		std::vector<double>fRecoTrackTruthMatchedLength; //distance travelled by matched particle
		std::vector<std::string>fRecoTrackTruthMatchedOrigin; //possible values: "Primary" , "Decay" , "Other"  -> origin of mc particle matched to this track
		
		
		std::vector<std::string> fRecoShowerTruthMatchedOrigin;
		std::vector<int> fRecoShowerTruthMatchedPDG;

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


		//other fhicl variables
		double fVertexSep=0; //maximum displacement between vertices of two PFPs that can be considered the same vertex
		double fVertexLowCut=0;     //minimum distance from reco primary vertex to be used when searching for decay vertex
		double fVertexHighCut=0;   //maximum distance from reco primary vertex to be used when searching for decay vertex

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
		double fTrackScore; //minimum track/shower score required for pfp to be classified as track



};


////////////////////////////////////////////////////
// Decides if an event passes selection criteria  //
////////////////////////////////////////////////////

bool hyperon::HyperonSelection::PerformSelection(){


//std::cout << "NPrimary Daughters: " << fNPrimaryDaughters << " NTracklikeDaughters:  " << fNPrimaryTrackDaughters << "  NShowerlikeDaughters:  " << fNPrimaryShowerDaughters << std::endl;


if(fRecoPrimaryVertex.size() != 3) return false; //if PV vector in uninitialized there is no reco'd neutrino


if(fLambdaCharged && fInActiveTPC) {
std::cout << std::endl;
std::cout << "EventID: " << fEventID << std::endl;
std::cout << "Num Shower Daughters: " <<  fNPrimaryShowerDaughters << std::endl;
std::cout << "Num Track Daughters: " << fNPrimaryTrackDaughters << std::endl; 
std::cout << std::endl;

}

//if(!inActiveTPC(fRecoPrimaryVertex)) return false;



	if(fNPrimaryDaughters < fMinDaughters) return false;
//	std::cout << "Passed Min Daughter cut" << std::endl;
	if(fNPrimaryShowerDaughters > fMaxPrimaryShowerDaughters) return false;
//	std::cout << "Pass Max Shower cut" << std::endl;
	if(fNPrimaryTrackDaughters < fMinDaughterTracks) return false;
//	std::cout << "Passed Min Track cut" << std::endl;
	
		return true;


}



/////////////////////////////////////////////////
// Other member functions                      //
/////////////////////////////////////////////////

 





hyperon::HyperonSelection::HyperonSelection(fhicl::ParameterSet const& p)
	: EDAnalyzer{p}   // ,
	// More initializers here.
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
	fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
	//std::cout << "BEFORE OPENING PARAMS" << std::endl;
	fPFParticleLabel = p.get<std::string>("PFParticleLabel");
	//std::cout<<fPFParticleLabel <<std::endl;
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
	fVertexSep = p.get<double>("VertexSeparation");
	fVertexLowCut = p.get<double>("VertexLowCut");
	fVertexHighCut = p.get<double>("VertexHighCut");


	fMinDaughterTracks = p.get<int>("MinDaughterTracks",0);
	fMaxPrimaryShowerDaughters = p.get<int>("MaxDaughterShowers",1000);
	fMinDaughters = p.get<int>("MinDaughters",0);	
	
	fTrackScore = p.get<double>("TrackScore",0.5);

	//misc parameters
	fPrint = p.get<bool>("Print",false);

	fPrintPdg = p.get<int>("PrintPdg",-1);
	fPrintChargedLambda = p.get<bool>("PrintChargedLambda",false);

}

void hyperon::HyperonSelection::analyze(art::Event const& e)
{



	//RESET EVERYTHING

	//generator truth info

	fInActiveTPC=false;
	fIsHyperon=false; //true if event contains a hyperon
	fIsSigmaZero=false; //true if event hyperon was a sigma zero		
	fIsLambda=false;

	//lepton information
	fLeptonMomentum.clear(); //lepton component wise momentum			
	fLeptonModMomentum=-1; //magnitude of lepton momentum
	fLeptonKE=-1;
	fLeptonTheta=-1;
	fLeptonPhi=-1;
	fLeptonPDG=-1; //outgoing lepton PDG code				

	//hyperon inforamtion
	fHyperonMomentum.clear(); //hyperon momentum component wise
	fHyperonModMomentum=-1; //magnitude of hyperon momentum
	fHyperonKE=-1;
	fHyperonTheta=-1;
	fHyperonPhi=-1;
	fHyperonPDG=-1; //hyperon pdg code

	//info for other nucleons/pions produced at primary vertex
	fPrimaryNucleonMomentum.clear(); //componentwise momentum of nucleons
	fPrimaryNucleonModMomentum.clear(); //magnitude of nucleon momentum
	fPrimaryNucleonKE.clear();
	fPrimaryNucleonPDG.clear(); //PDG codes of nucleons

	fPrimaryPionMomentum.clear(); //componentwise momentum of nucleons
	fPrimaryPionKE.clear();
	fPrimaryPionModMomentum.clear(); //magnitude of nucleon momentum
	fPrimaryPionPDG.clear(); //PDG codes of nucleons


	//vertex information
	fTruePrimaryVertex.clear();

	//neutrino information
	fNuEnergy=-1; //neutrino energy
	fNuPDG=-1; //neutrino PDG code

	//g4 truth info

	fLeptonTravel=-1;
	fDecayVertex.clear(); //position of hyperon decay vertex
	fDecayPDG.clear(); //pdg codes of decay products
	fDecayMomenta.clear(); //component wise momenta of decay products
	fDecayModMomenta.clear(); //magnitude of momenta of hyperon decay products
	fDecayKE.clear(); //KE of decay products
	fDecayTheta.clear();
	fDecayPhi.clear();
	fHyperonTravel=-1; //distance travelled by hyperon before decaying	
	fDecayTravel.clear(); //distance travelled by hyperon decay products

	fSigmaZeroDecayPhotonMomentum.clear(); //component wise momentum of photon produced during Sigma0 decay
	fSigmaZeroDecayPhotonModMomentum=-1; //magnitude of momentum of photon produced during Sigma0 decay
	fSigmaZeroDecayLambdaMomentum.clear(); //component wise momentum of Lambda produced during Sigma0 decay
	fSigmaZeroDecayLambdaModMomentum=-1; //magnitude of momentum of Lambda produced during Sigma0 decay
	fSigmaZeroDecayLambdaKE=-1;

 
        fDecayOpeningAngle=-1; //opening angle between hyperon decay products   
        fDecayPlaneVector.clear(); //vector normal to plane containing hyperon (cross product of pion and nucleon momenta)
        fLeptonPionAngle=-1; //openining angle between lepton and pion
        fLeptonNucleonAngle=-1; //opening angle between lepton and nucleon
	fDecayPlaneTheta=-1;
	fDecayPlanePhi=-1;

	fLambdaCharged=false;

	//RECO LEVEL INFO

	fSelectedEvent = false; //true if event passes some selection criteria	
	fNPrimaryDaughters = 0; //number of primary daughters

	fNPrimaryTrackDaughters=0; //num of track like primary daughters
	fNPrimaryShowerDaughters=0; //num of shower like primary daughters



	fRecoPrimaryVertex.clear(); //position of reco'd primary vertex
	fDaughterTrackLengths.clear(); //lengths of track like primary daughters
	fDaughterTrackPIDs.clear(); //proton PID scores of track like primary daughters
	fDaughterTrackVertices.clear(); //vertices of track like primary daughters
	fDaughterVertexDisplacement.clear(); //distance between reco'd primary vertex and vertices of track like primary daughters

	//reco Truth matched info
	fRecoTrackTruthMatchedPDG.clear();
	fRecoTrackTruthMatchedVertex.clear();
	fRecoTrackTruthMatchedMomentum.clear(); //momentum of matched particle
	fRecoTrackTruthMatchedModMomentum.clear(); //magnitude of momentum of matched particle
	fRecoTrackTruthMatchedKE.clear(); //KE of matched particle
	fRecoTrackTruthMatchedLength.clear(); //distance travelled by matched particle
	fRecoTrackTruthMatchedOrigin.clear(); // origin of mc particle matched to this track

	fTrackTrackScores.clear(); //pandora track/shower score : <0.5 is shower, >0.5 is track
	fShowerTrackScores.clear(); //pandora track/shower score : <0.5 is shower, >0.5 is track
	fShowerDirections.clear(); //shower directions

	fRecoShowerTruthMatchedOrigin.clear();
	fRecoShowerTruthMatchedPDG.clear();

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
	ccnc = 	Nu.CCNC();

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


				fTruePrimaryVertex = {Part.Vx() , Part.Vy() , Part.Vz()};

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

		//if(abs(g4p->PdgCode()) != 11 && g4p->PdgCode()  != 22 && g4p->PdgCode() < 10000) 
		//	std::cout << "ID No: " << g4p->TrackId()  << " Mother ID: " << g4p->Mother() << " PDG:  " << g4p->PdgCode() << std::endl;


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

						fDecayVertex = {g4p->EndPosition().X() , g4p->EndPosition().Y() , g4p->EndPosition().Z() };
						fHyperonTravel = sqrt( (g4p->EndPosition().X() - fTruePrimaryVertex.at(0))*(g4p->EndPosition().X() - fTruePrimaryVertex.at(0))
								+ (g4p->EndPosition().Y() - fTruePrimaryVertex.at(1))*(g4p->EndPosition().Y() - fTruePrimaryVertex.at(1))
								+ (g4p->EndPosition().Z() - fTruePrimaryVertex.at(2))*(g4p->EndPosition().Z() - fTruePrimaryVertex.at(2)) );



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

			fIsHyperon = true;
			fHyperonMomentum = {part->E() , part->Px() , part->Py() , part->Pz()};
			fHyperonModMomentum = sqrt(part->Px()*part->Px() + part->Py()*part->Py() + part->Pz()*part->Pz());
			fHyperonKE = part->E() - part->Mass();	
			fHyperonPDG = part->PdgCode();				
			fHyperonTheta = (180/3.1416)*TMath::ACos(part->Pz()/fHyperonModMomentum);
	
			fHyperonPhi = (180/3.1416)*TMath::ASin(part->Py()/fHyperonModMomentum);
		

		
	//		std::cout << "Hyperon momentum: " << fHyperonMomentum.at(1) << "   " << fHyperonMomentum.at(2) << "   " << fHyperonMomentum.at(3) << "    Theta=" << fHyperonTheta << "   Phi=" << fHyperonPhi <<   std::endl; 

		}

		//lepton produced at primary vertex
		if(isLepton(part->PdgCode()) || isNeutrino(part->PdgCode())){
		


			fLeptonPDG = part->PdgCode();
			fLeptonMomentum = {part->E() , part->Px() , part->Py() , part->Pz()};
			fLeptonModMomentum = sqrt(part->Px()*part->Px() + part->Py()*part->Py() + part->Pz()*part->Pz());
			fLeptonKE = part->E() - part->Mass();

			fLeptonTheta = (180/3.1416)*TMath::ACos(part->Pz()/fLeptonModMomentum);
			fLeptonPhi = (180/3.1416)*TMath::ASin(part->Py()/fLeptonModMomentum);
		
					
	//		std::cout << "Lepton momentum: " << fLeptonMomentum.at(1) << "   " << fLeptonMomentum.at(2) << "   " << fLeptonMomentum.at(3) << "    Theta=" << fLeptonTheta << "   Phi=" << fLeptonPhi <<   std::endl; 


			//record distance travelled by lepton	
			fLeptonTravel = sqrt( (part->EndPosition().X() - fTruePrimaryVertex.at(0))*(part->EndPosition().X() - fTruePrimaryVertex.at(0))
					+ (part->EndPosition().Y() - fTruePrimaryVertex.at(1))*(part->EndPosition().Y() - fTruePrimaryVertex.at(1))
					+ (part->EndPosition().Z() - fTruePrimaryVertex.at(2))*(part->EndPosition().Z() - fTruePrimaryVertex.at(2)) );



		}

		if(isNucleon(part->PdgCode()) ){

			std::vector<double> v = {part->E() , part->Px() , part->Py() , part->Pz()};

			fPrimaryNucleonMomentum.push_back(v);
			fPrimaryNucleonModMomentum.push_back(sqrt(part->Px()*part->Px() + part->Py()*part->Py() + part->Pz()*part->Pz()));
			fPrimaryNucleonKE.push_back(part->E() - part->Mass());	

			fPrimaryNucleonPDG.push_back( part->PdgCode());		

		}

		if(isPion(part->PdgCode())){

			std::vector<double> v = {part->E() , part->Px() , part->Py() , part->Pz()};

			fPrimaryPionMomentum.push_back(v);
			fPrimaryPionModMomentum.push_back(sqrt(part->Px()*part->Px() + part->Py()*part->Py() + part->Pz()*part->Pz()));
			fPrimaryPionKE.push_back(part->E() - part->Mass());	

			fPrimaryPionPDG.push_back( part->PdgCode());		

		}


	}



	//now go through list of hyperon daughters, get info about the decay
	for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){
		//geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
		if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

		art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these


		std::vector<double> p = {part->E() , part->Px(),part->Py(),part->Pz()  };

		double P =  sqrt(p.at(1)*p.at(1) + p.at(2)*p.at(2) + p.at(3)*p.at(3));

		fDecayTravel.push_back( sqrt( (part->EndPosition().X() - fDecayVertex.at(0))*(part->EndPosition().X() - fDecayVertex.at(0))
					+ (part->EndPosition().Y() - fDecayVertex.at(1))*(part->EndPosition().Y() - fDecayVertex.at(1))
					+ (part->EndPosition().Z() - fDecayVertex.at(2))*(part->EndPosition().Z() - fDecayVertex.at(2)) ) );


		fDecayPDG.push_back( part->PdgCode() );

		fDecayMomenta.push_back(p);
                fDecayModMomenta.push_back(P);
		fDecayKE.push_back(part->E() - part->Mass());	
		
			fDecayTheta.push_back((180/3.1416)*TMath::ACos(part->Pz()/P));
			fDecayPhi.push_back((180/3.1416)*TMath::ASin(part->Py()/P));


 if(isPion(part->PdgCode()) && fLeptonMomentum.size()) {


fLeptonPionAngle = (180/3.1416)*TMath::ACos( (p.at(1)*fLeptonMomentum.at(1) + p.at(2)*fLeptonMomentum.at(2) + p.at(3)*fLeptonMomentum.at(3) )/(fLeptonModMomentum*P) );

        }

        if(isNucleon(part->PdgCode()) && fLeptonMomentum.size()){

fLeptonNucleonAngle = (180/3.1416)*TMath::ACos( (p.at(1)*fLeptonMomentum.at(1) + p.at(2)*fLeptonMomentum.at(2) + p.at(3)*fLeptonMomentum.at(3) )/(fLeptonModMomentum*P) );

        }


	}


if(fDecayMomenta.size() == 2){
fDecayOpeningAngle = (180/3.1416)*TMath::ACos( (fDecayMomenta.at(0).at(3)*fDecayMomenta.at(1).at(3) + fDecayMomenta.at(0).at(1)*fDecayMomenta.at(1).at(1) + fDecayMomenta.at(0).at(2)*fDecayMomenta.at(1).at(2) )/(fDecayModMomenta.at(0)*fDecayModMomenta.at(1))  );


//calculate cross product of nucleon and pion momenta
double x = fDecayMomenta.at(0).at(2)*fDecayMomenta.at(1).at(3) - fDecayMomenta.at(1).at(2)*fDecayMomenta.at(0).at(3);
double y = (-1)*( fDecayMomenta.at(0).at(1)*fDecayMomenta.at(1).at(3) - fDecayMomenta.at(0).at(3)*fDecayMomenta.at(1).at(1) );
double z = fDecayMomenta.at(0).at(1)*fDecayMomenta.at(1).at(2) - fDecayMomenta.at(1).at(1)*fDecayMomenta.at(0).at(2);

double r = sqrt(x*x + y*y + z*z);

//unit vector in direction of cross product of decay product momenta
fDecayPlaneVector = { x/r , y/r , z/r };

double phi = 90 - (180/3.1416)*TMath::ACos( fDecayPlaneVector.at(0) );
if(phi<0) phi += 180;

fDecayPlanePhi = phi;

double theta = 90 - (180/3.1416)*TMath::ACos( fDecayPlaneVector.at(2) );
if(theta<0) theta += 180;

fDecayPlaneTheta = theta;

}


	//if hyperon was a lambda, check the pdg codes of its decay products


	if(fIsLambda && fDecayPDG.size() == 2 && ( (fDecayPDG.at(0) == 2212 && fDecayPDG.at(1) == -211) || (fDecayPDG.at(1) == 2212 && fDecayPDG.at(0) == -211) ) ) fLambdaCharged=true;
		



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
				fSigmaZeroDecayPhotonMomentum = p;
				fSigmaZeroDecayPhotonModMomentum = sqrt(p.at(1)*p.at(1) + p.at(2)*p.at(2) + p.at(3)*p.at(3));

				//add its ID to list of primary ID's - produced at primary vtx so from reco point of view makes sense to put it in
				//list of primary IDs
				//

				primary_IDs.push_back(part->TrackId());



			}
			//Lambda produced from sigma zero decay
			else if(part->PdgCode() == 3122){


				fSigmaZeroDecayLambdaMomentum = p;
				fSigmaZeroDecayLambdaModMomentum = sqrt(p.at(1)*p.at(1) + p.at(2)*p.at(2) + p.at(3)*p.at(3));

				fSigmaZeroDecayLambdaKE = part->E() - part->Mass();	

				//now search for Lambda decay products, find its decay vertex etc			




				if(part->EndProcess() == "Decay"){		

					//get Lambdas decay vertex

					fDecayVertex = {part->EndPosition().X() , part->EndPosition().Y() , part->EndPosition().Z() };

					fHyperonTravel = sqrt( (part->EndPosition().X() - fTruePrimaryVertex.at(0))*(part->EndPosition().X() - fTruePrimaryVertex.at(0))
							+ (part->EndPosition().Y() - fTruePrimaryVertex.at(1))*(part->EndPosition().Y() - fTruePrimaryVertex.at(1))
							+ (part->EndPosition().Z() - fTruePrimaryVertex.at(2))*(part->EndPosition().Z() - fTruePrimaryVertex.at(2)) );





					for(int i_d2=0;i_d2<part->NumberDaughters();i_d2++){

						//search the map for the daughter IDs
						if(partByID.find(part->Daughter(i_d2)) == partByID.end()) continue;

						art::Ptr<simb::MCParticle> part2 = partByID[part->Daughter(i_d2)];

						if(part2->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these


						std::vector<double> p2 = {part2->E() , part2->Px(),part2->Py(),part2->Pz()  };


						fDecayTravel.push_back( sqrt( (part2->EndPosition().X() - fDecayVertex.at(0))*(part2->EndPosition().X() - fDecayVertex.at(0))
									+ (part2->EndPosition().Y() - fDecayVertex.at(1))*(part2->EndPosition().Y() - fDecayVertex.at(1))
									+ (part2->EndPosition().Z() - fDecayVertex.at(2))*(part2->EndPosition().Z() - fDecayVertex.at(2)) ) );


						fDecayPDG.push_back( part2->PdgCode() );

						fDecayMomenta.push_back(p2);
						fDecayModMomenta.push_back( sqrt(p2.at(1)*p2.at(1) + p2.at(2)*p2.at(2) + p2.at(3)*p2.at(3)) );
						fDecayKE.push_back(part2->E() - part2->Mass());	


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

//space charge correction

//auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();




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
		
		fRecoPrimaryVertex = {vtx->position().X() , vtx->position().Y() , vtx->position().Z()};


//		geo::Point_t point = {120.0 , 0.0 , 500.0};


//		geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
	//	std::cout << sce_corr.X() << "  " << sce_corr.Y() << "  " << sce_corr.Z() << std::endl;


		}




	}

}



for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){




	//get data from every PFP, not just neutrino daughters

	if(pfp->Parent() != neutrinoID) continue;



	std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
	std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());
	
	std::vector< art::Ptr<recob::Shower> > pfpShowers = showerAssoc.at(pfp.key());
	std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> >pfpMeta = metadataAssoc.at(pfp.key());

//std::cout << pfp->PdgCode() << std::endl;

	if(pfp->PdgCode() == 11) fNPrimaryShowerDaughters++;
	if(pfp->PdgCode() == 13) fNPrimaryTrackDaughters++;



	//track/shower score
	double thisParticleTrackScore;


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
			thisParticleTrackScore = it->second;
		//	if(it->second > fTrackScore) fNPrimaryTrackDaughters++;
	//		else fNPrimaryShowerDaughters++;
			}
		}
	}	
}



	//if pfp is a track like object
	if(!pfpTracks.empty() && !pfpVertex.empty()){

		for(const art::Ptr<recob::Track> &trk : pfpTracks){

			fTrackTrackScores.push_back(thisParticleTrackScore);

			fDaughterTrackLengths.push_back(trk->Length());

			

			//start of truth match


			//get hits assoc with track
			std::vector< art::Ptr< recob::Hit> > hits = trackHitAssoc.at(trk.key());
			//std::cout << "Hits assoc with this track: " << hits.size() << std::endl;

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

			//pdg of mc particle that contributed the edep to this track
			fRecoTrackTruthMatchedPDG.push_back(matchedParticle->PdgCode());

			fRecoTrackTruthMatchedVertex.push_back( { matchedParticle->Vx() , matchedParticle->Vy() , matchedParticle->Vz() } );

			fRecoTrackTruthMatchedID.push_back(matchedParticle->TrackId());


			std::vector<double> p = { matchedParticle->E() , matchedParticle->Px() , matchedParticle->Py() , matchedParticle->Pz() };
			fRecoTrackTruthMatchedMomentum.push_back(p);
			fRecoTrackTruthMatchedModMomentum.push_back( sqrt( p.at(1)*p.at(1) + p.at(2)*p.at(2) + p.at(3)*p.at(3) ) ); 		
			fRecoTrackTruthMatchedKE.push_back(matchedParticle->E() - matchedParticle->Mass());




			//calculate distance travelled by particle	
			fRecoTrackTruthMatchedLength.push_back( sqrt( (matchedParticle->Vx() - matchedParticle->EndX())*(matchedParticle->Vx() - matchedParticle->EndX())
						+ (matchedParticle->Vy() - matchedParticle->EndY())*(matchedParticle->Vy() - matchedParticle->EndY())	
						+ (matchedParticle->Vz() - matchedParticle->EndZ())*(matchedParticle->Vz() - matchedParticle->EndZ()) ) );
	

			//origin of truth matched particle
			fRecoTrackTruthMatchedOrigin.push_back(getOrigin(matchedParticle->TrackId()));

}



//std::cout << "Matched Track PDG: " <<  matchedParticle->PdgCode() <<  "  Origin: " << getOrigin(matchedParticle->TrackId()) <<  std::endl;

	


			//get PID for track


			std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());
			std::vector<art::Ptr<anab::ParticleID>> trackPID = PIDAssoc.at(trk.key());

			std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

			for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){

				anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);



				if(  TMath::Abs(AlgScore.fAssumedPdg) == 2212  &&   AlgScore.fAlgName=="ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){	
					//	std::cout << "Found thee plane PID!" << std::endl;			

					fDaughterTrackPIDs.push_back(AlgScore.fValue);


				}

			}
		}



		for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){



			double d = sqrt( (vtx->position().X() - fRecoPrimaryVertex.at(0))*(vtx->position().X() - fRecoPrimaryVertex.at(0))
					+ (vtx->position().Y() - fRecoPrimaryVertex.at(1))*(vtx->position().Y() - fRecoPrimaryVertex.at(1))
					+ (vtx->position().Z() - fRecoPrimaryVertex.at(2))*(vtx->position().Z() - fRecoPrimaryVertex.at(2)) );

			fDaughterTrackVertices.push_back(std::vector<double>{vtx->position().X() , vtx->position().Y() , vtx->position().Z()});
			fDaughterVertexDisplacement.push_back(d);



		}





	}




}


if(fPrint) PrintInfo();

fSelectedEvent = PerformSelection();

FinishEvent();



}

// Print some useful info from the event


void hyperon::HyperonSelection::PrintInfo(){

if(fPrintChargedLambda && !fLambdaCharged) return;

if(!(fHyperonPDG == fPrintPdg || fPrintPdg == -1) || !fInActiveTPC) return;

std::cout << std::endl;	
std::cout << "EventID: " << fEventID-1 << std::endl;


std::cout << "Truth info" << std::endl;
std::cout << "Hyperon PDG: " << fHyperonPDG << std::endl;
std::cout << "Lepton PDG: " << fLeptonPDG << std::endl;
std::cout << "Dist travelled by hyperon: " << fHyperonTravel << std::endl;
std::cout << "Hyperon decay products: "<<std::endl; 

for(size_t i=0;i<fDecayPDG.size();i++){

std::cout <<"pdg: " <<  fDecayPDG.at(i) << " length: " << fDecayTravel.at(i) << "  KE:  " << fDecayKE.at(i) <<  std::endl;
std::cout << "Direction:  " <<  fDecayMomenta.at(i).at(1)/fDecayModMomenta.at(i) << "   " <<   fDecayMomenta.at(i).at(2)/fDecayModMomenta.at(i) << "   " <<   fDecayMomenta.at(i).at(3)/fDecayModMomenta.at(i) << " Theta:  "  << fDecayTheta.at(i) << "  Phi:  " << fDecayPhi.at(i) <<  std::endl; 
}


if(fDecayPDG.size() == 2){
std::cout << "Decay Opening Angle:  " << fDecayOpeningAngle << std::endl;
std::cout << "Angle between lepton and pion:  " << fLeptonPionAngle << std::endl; //openining angle between lepton and pion
std::cout << "Angle between lepton and nucleon:  " <<   fLeptonNucleonAngle << std::endl; //opening angle between lepton and nucleon


//get angle of decay plane to x axis

std::cout << "Decay Plane theta: " << fDecayPlaneTheta << std::endl;
std::cout << "DecayPlane phi: " << fDecayPlanePhi << std::endl;

}

std::cout << std::endl;

std::cout << "Reco info" << std::endl;

std::cout << "Num tracklike daughters: " << fNPrimaryTrackDaughters << "  Num showerlike daughters:  " << fNPrimaryShowerDaughters << std::endl;

std::cout << "Tracklike Daughters: " << std::endl;
for(size_t i=0;i<fDaughterTrackVertices.size();i++){
std::cout << std::endl;
std::cout << "Track Score: " << fTrackTrackScores.at(i) << std::endl;
std::cout << "Dist from PV: " <<  fDaughterVertexDisplacement.at(i)  << "  Proton PID Score:  " << fDaughterTrackPIDs.at(i)  << "  Reco Track Len:  " <<  fDaughterTrackLengths.at(i) << " Matched PDG: " << fRecoTrackTruthMatchedPDG.at(i) << "  Origin:  " << getOrigin(fRecoTrackTruthMatchedID.at(i)) <<    std::endl;


}



std::cout << std::endl;
std::cout << "Showerlike daughters: " << std::endl;
for(size_t i_sh=0;i_sh<fShowerTrackScores.size();i_sh++){
std::cout << std::endl;
std::cout << "Track Score: " << fShowerTrackScores.at(i_sh) << std::endl;
std::cout << "Pdg: " <<  fRecoShowerTruthMatchedPDG.at(i_sh) <<  "  Origin: " << 	fRecoShowerTruthMatchedOrigin.at(i_sh) << std::endl;

}


std::cout << std::endl;
}


//check origin of particle by its TrackId

std::string hyperon::HyperonSelection::getOrigin(int idnum){

std::string origin;
bool found=false;

			for(size_t i_p=0;i_p<primary_IDs.size();i_p++){

				if(primary_IDs.at(i_p)	== idnum){ 
					origin = "Primary"; 
					found = true;
					break;
				}

			}

			if(!found){
				for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){

					if(daughter_IDs.at(i_d)	== idnum){ 
						origin = "HyperonDecay";
						found = true;
						break;
					}
				}
			}


			if(!found) origin = "Other";

return origin;

}


//updates metadata information
void hyperon::HyperonSelection::FinishEvent(){


	if(!isNeutrino(fLeptonPDG)) { fNChargedCurrent++; fCCNC="CC"; }
	else { fNNeutralCurrent++; fCCNC="NC"; }

	if(fNuPDG == 12) fNnue++;
	else if(fNuPDG == 14) fNnuMu++;
	else if(fNuPDG == -12) fNnueBar++;
	else fNnuMuBar++;


	//store reco'd info
	fTruthTree->Fill();
	fRecoTree->Fill();

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

	if(fLambdaCharged) fTotalChargedLambdas++;
	if(fLambdaCharged && fInActiveTPC) fTotalChargedLambdasInActiveVol++;
	if(fSelectedEvent && fLambdaCharged) fNSelectedChargedLambdas++;
	if(fSelectedEvent && fLambdaCharged && fInActiveTPC) fNSelectedChargedLambdasInActiveVol++;

	//efficiency = fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiv;
	//purity = fNSelectedHyperonsInActiveVol/fSelectedEvents;


}




void hyperon::HyperonSelection::beginJob(){

	fileID=0;	


	// Implementation of optional member functions

	art::ServiceHandle<art::TFileService> tfs;

	fTruthTree=tfs->make<TTree>("TruthTree","Truth Info Tree");



	//write a second tree to store decay product info
	fTruthTree->Branch("EventID",&fEventID);
	fTruthTree->Branch("run",&run);
	fTruthTree->Branch("subrun",&subrun);
	fTruthTree->Branch("event",&event);
	fTruthTree->Branch("fileID",&fileID);

	fTruthTree->Branch("Mode",&fMode);
	fTruthTree->Branch("CCNC",&fCCNC);

	//generator truth info

	fTruthTree->Branch("InActiveTPC",&fInActiveTPC);
	fTruthTree->Branch("IsHyperon",&fIsHyperon); //true if event contains a hyperon
	fTruthTree->Branch("IsSigmaZero",&fIsSigmaZero); //true if event hyperon was a sigma zero		
	fTruthTree->Branch("IsLambda",&fIsLambda);
	fTruthTree->Branch("LambdaCharged",&fLambdaCharged);

	fTruthTree->Branch("SelectedEvent",&fSelectedEvent);

	//lepton information
	fTruthTree->Branch("LeptonMomentum",&fLeptonMomentum); //lepton component wise momentum			
	fTruthTree->Branch("LeptonModMomentum",&fLeptonModMomentum); //magnitude of lepton momentum
	fTruthTree->Branch("LeptonKE",&fLeptonKE); //magnitude of lepton momentum
	fTruthTree->Branch("LeptonTheta",&fLeptonTheta);
	fTruthTree->Branch("LeptonPhi",&fLeptonPhi);

	fTruthTree->Branch("LeptonPDG",&fLeptonPDG); //outgoing lepton PDG code				

	//hyperon inforamtion
	fTruthTree->Branch("HyperonMomentum",&fHyperonMomentum); //hyperon momentum component wise
	fTruthTree->Branch("HyperonModMomentum",&fHyperonModMomentum); //magnitude of hyperon momentum
	fTruthTree->Branch("HyperonKE",&fHyperonKE); //magnitude of hyperon momentum
	fTruthTree->Branch("HyperonPDG",&fHyperonPDG); //hyperon pdg code
	fTruthTree->Branch("HyperonTheta",&fHyperonTheta);
	fTruthTree->Branch("HyperonPhi",&fHyperonPhi);

	//info for other nucleons/pions produced at primary vertex
	fTruthTree->Branch("PrimaryNucleonMomentum",&fPrimaryNucleonMomentum); //componentwise momentum of nucleons
	fTruthTree->Branch("PrimaryNucleonModMomentum",&fPrimaryNucleonModMomentum); //magnitude of nucleon momentum
	fTruthTree->Branch("PrimaryNucleonKE",&fPrimaryNucleonKE); //magnitude of nucleon momentum
	fTruthTree->Branch("PrimaryNucleonPDG",&fPrimaryNucleonPDG); //PDG codes of nucleons

	fTruthTree->Branch("PrimaryPionMomentum",&fPrimaryPionMomentum); //componentwise momentum of nucleons
	fTruthTree->Branch("PrimaryPionModMomentum",&fPrimaryPionModMomentum); //magnitude of nucleon momentum
	fTruthTree->Branch("PrimaryPionKE",&fPrimaryPionKE); //magnitude of nucleon momentum
	fTruthTree->Branch("PrimaryPionPDG",&fPrimaryPionPDG); //PDG codes of nucleons


	//vertex information
	fTruthTree->Branch("TruePrimaryVertex",&fTruePrimaryVertex);

	//neutrino information
	fTruthTree->Branch("NuEnergy",&fNuEnergy); //neutrino energy
	fTruthTree->Branch("NuPDG",&fNuPDG); //neutrino PDG code

	//g4 truth info

	fTruthTree->Branch("LeptonTravel",&fLeptonTravel);
	fTruthTree->Branch("DecayVertex",&fDecayVertex); //position of hyperon decay vertex
	fTruthTree->Branch("DecayPDG",&fDecayPDG); //pdg codes of decay products
	fTruthTree->Branch("DecayMomenta",&fDecayMomenta); //component wise momenta of decay products
	fTruthTree->Branch("DecayModMomenta",&fDecayModMomenta); //magnitude of momenta of hyperon decay products
	fTruthTree->Branch("DecayKE",&fDecayKE); //magnitude of momenta of hyperon decay products
	fTruthTree->Branch("DecayTheta",&fDecayTheta);
	fTruthTree->Branch("DecayPhi",&fDecayPhi);
	fTruthTree->Branch("HyperonTravel",&fHyperonTravel); //distance travlled by hyperon
	fTruthTree->Branch("DecayTravel",&fDecayTravel); //distance travelled by hyperon decay products


	fTruthTree->Branch("SigmaZeroDecayPhotonMomentum",&fSigmaZeroDecayPhotonMomentum); //component wise momentum of photon produced during Sigma0 decay
	fTruthTree->Branch("SigmaZeroDecayPhotonModMomentum",&fSigmaZeroDecayPhotonModMomentum); //magnitude of momentum of photon produced during Sigma0 decay
	fTruthTree->Branch("SigmaZeroDecayLambdaMomentum",&fSigmaZeroDecayLambdaMomentum); //component wise momentum of Lambda produced during Sigma0 decay
	fTruthTree->Branch("SigmaZeroDecayLambdaModMomentum",&fSigmaZeroDecayLambdaModMomentum); //magnitude of momentum of Lambda produced during Sigma0 decay
	fTruthTree->Branch("SigmaZeroDecayLambdaKE",&fSigmaZeroDecayLambdaKE); //magnitude of momentum of Lambda produced during Sigma0 decay


	//information about plane containing decay products and 
	//angles between particles involved

	fTruthTree->Branch("DecayOpeningAngle",&fDecayOpeningAngle); //opening angle between hyperon decay products     
        fTruthTree->Branch("DecayPlaneVector",&fDecayPlaneVector); //vector normal to plane containing hyperon (cross product of pion and nucleon momenta)
        fTruthTree->Branch("LeptonPionAngle",&fLeptonPionAngle); //openining angle between lepton and pion
        fTruthTree->Branch("LeptonNucleonAngle",&fLeptonNucleonAngle); //opening angle between lepton and nucleon
	fTruthTree->Branch("DecayPlanePhi",&fDecayPlanePhi);
	fTruthTree->Branch("DecayPlaneTheta",&fDecayPlaneTheta);

	std::cout << "Setup Reco Tree" << std::endl;

	fRecoTree=tfs->make<TTree>("RecoTree","Reconstructed Info Tree");

	fRecoTree->Branch("fileID",&fileID);
	fRecoTree->Branch("EventID",&fEventID);
	fRecoTree->Branch("run",&run);
	fRecoTree->Branch("subrun",&subrun);
	fRecoTree->Branch("event",&event);

	
	fRecoTree->Branch("Mode",&fMode);
	fRecoTree->Branch("CCNC",&fCCNC);

	fRecoTree->Branch("InActiveTPC",&fInActiveTPC);
	fRecoTree->Branch("IsHyperon",&fIsHyperon); //true if event contains a hyperon
	fRecoTree->Branch("IsLambda",&fIsLambda);	
	fRecoTree->Branch("LambdaCharged",&fLambdaCharged);

	fRecoTree->Branch("SelectedEvent",&fSelectedEvent);
	fRecoTree->Branch("RecoPrimaryVertex",&fRecoPrimaryVertex);
	fRecoTree->Branch("NPrimaryDaughters",&fNPrimaryDaughters);
	fRecoTree->Branch("NPrimaryTrackDaughters",&fNPrimaryTrackDaughters); //num ofRecoTree->Branch("",f track like primary daughters
	fRecoTree->Branch("NPrimaryShowerDaughters",&fNPrimaryShowerDaughters); //num ofRecoTree->Branch("",f shower like primary daughters

	fRecoTree->Branch("DaughterTrackLengths",&fDaughterTrackLengths);
	fRecoTree->Branch("DaughterTrackPIDs",&fDaughterTrackPIDs);
	fRecoTree->Branch("DaughterTrackVertices",&fDaughterTrackVertices);
	fRecoTree->Branch("DaughterVertexDisplacement",&fDaughterVertexDisplacement);
	fRecoTree->Branch("DaughterTrackScores",&fTrackTrackScores);

	fRecoTree->Branch("RecoTrackTruthMatchedPDG",&fRecoTrackTruthMatchedPDG);
	fRecoTree->Branch("RecoTrackTruthMatchedVertex",&fRecoTrackTruthMatchedVertex); //possible values: "Primary" , "Decay" , "Other"  -> origin of mc particle matched to this track
	fRecoTree->Branch("RecoTrackTruthMatchedMomentum",&fRecoTrackTruthMatchedMomentum);
	fRecoTree->Branch("RecoTrackTruthMatcheModMomentum",&fRecoTrackTruthMatchedModMomentum);
	fRecoTree->Branch("RecoTrackTruthMatchedKE",&fRecoTrackTruthMatchedKE);
	fRecoTree->Branch("RecoTrackTruthMatchedLength",&fRecoTrackTruthMatchedLength);
	fRecoTree->Branch("RecoTrackTruthMatchedOrigin",&fRecoTrackTruthMatchedOrigin); //possible values: "Primary" , "Decay" , "Other"  -> origin of mc particle matched to this track


	fRecoTree->Branch("ShowerTrackScores",&fShowerTrackScores);
	fRecoTree->Branch("RecoShowerTruthMatchedPDG",&fRecoShowerTruthMatchedPDG);
	fRecoTree->Branch("RecoShowerTruthMatchedOrigin",&fRecoShowerTruthMatchedOrigin);

	//metadata tree


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
	fMetaTree->Branch("NSelectedEvents",&fNSelectedEvents);
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


}



void hyperon::HyperonSelection::endJob()
{
	//calculate efficiency and purity

	//efficiency = fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiveVol;
	//purity = fNSelectedHyperonsInActiveVol/fSelectedEvents;

	if(fNSelectedEvents > 0){
 	fHyperonPurity = (double)fNSelectedHyperonsInActiveVol/fNSelectedEvents;
	fLambdaPurity = (double)fNSelectedLambdasInActiveVol/fNSelectedEvents;
	fChargedLambdaPurity = (double)fNSelectedChargedLambdasInActiveVol/fNSelectedEvents;
	}

	if(fTotalHyperonsInActiveVol > 0) fHyperonEfficiency = (double)fNSelectedHyperonsInActiveVol/fTotalHyperonsInActiveVol;

	if(fTotalLambdasInActiveVol > 0) fLambdaEfficiency = (double)fNSelectedLambdasInActiveVol/fTotalLambdasInActiveVol;

	if(fTotalChargedLambdasInActiveVol > 0) fChargedLambdaEfficiency = (double)fNSelectedChargedLambdasInActiveVol/fTotalChargedLambdasInActiveVol;

	fHyperonEfficiencyTimesPurity = fHyperonPurity * fHyperonEfficiency;
	fHyperonTruePurity = 1/( 1+(1/fHyperonPurity-1)*50 );
	fHyperonEfficiencyTimesTruePurity = fHyperonEfficiency * fHyperonTruePurity;

		
	fLambdaEfficiencyTimesPurity = fLambdaPurity * fLambdaEfficiency;
	fLambdaTruePurity = 1/( 1+(1/fLambdaPurity-1)*50 );
	fLambdaEfficiencyTimesTruePurity = fLambdaEfficiency * fLambdaTruePurity;


	fChargedLambdaEfficiencyTimesPurity = fChargedLambdaPurity * fChargedLambdaEfficiency;
	fChargedLambdaTruePurity = 1/( 1+(1/fChargedLambdaPurity-1)*50 );
	fChargedLambdaEfficiencyTimesTruePurity = fChargedLambdaEfficiency * fChargedLambdaTruePurity;


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
