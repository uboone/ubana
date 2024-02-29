////////////////////////////////////////////////////////////////////////
// Class:       BlipVertexAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        BlipVertexAnalysis_module.cc
//
// Generated at Tue Feb  6 00:18:35 2024 by Keng Lin using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Pandora/PdgTable.h"

// extra includes, trying to fix
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"

// for MCParticle.E() and such
#include "nusimdata/SimulationBase/MCTrajectory.h"
// for associations
#include "canvas/Persistency/Common/FindManyP.h" // find associations as pointers
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "lardataobj/RecoBase/Track.h"
// this line good for "new" larsoft:
#include "art/Persistency/Common/PtrMaker.h"
// for outdated versions (i.e MCC8) use this line:
//#include "lardata/Utilities/PtrMaker.h"

// additional framework includes
#include "art/Framework/Services/Optional/TFileService.h"

// pot summary
// from /uboone/app/users/afurmans/muonScatteringAnalysis/
//       larsoft/srcs/uboonecode/uboone/AnalysisTree/
//       AnalysisTree_module.cc
#include "larcoreobj/SummaryData/POTSummary.h"


#include "Libraries/init_branches.h"
#include "Libraries/variables.h"
#include "Libraries/BVA_onlyblips.h"
#include "Libraries/BVA_onlypandora.h"
#include "Libraries/BVA_blipsnearvertex.h"

namespace BVA_ana
{
//	class BlipVertexAnalysis;


	class BlipVertexAnalysis : public art::EDAnalyzer {
		public:
			explicit BlipVertexAnalysis(fhicl::ParameterSet const& p);
			// The compiler-generated destructor is fine for non-base
			// classes without bare pointers or other resource use.

			// Plugins should not be copied or assigned.
			BlipVertexAnalysis(BlipVertexAnalysis const&) = delete;
			BlipVertexAnalysis(BlipVertexAnalysis&&) = delete;
			BlipVertexAnalysis& operator=(BlipVertexAnalysis const&) = delete;
			BlipVertexAnalysis& operator=(BlipVertexAnalysis&&) = delete;

			// Required functions.
			void analyze(art::Event const& e) override;

			// Selected optional functions.
			void beginJob() override;
			void endJob() override;
			void endSubRun(const art::SubRun& sr);
			
			var_all vars;
			var_evt vars_evt;

		private:

			// Declare member data here.
			// blip
			blip::BlipRecoAlg* fBlipAlg;

			std::string fPOTModuleLabel;
	};


	BlipVertexAnalysis::BlipVertexAnalysis(fhicl::ParameterSet const& p)
		: EDAnalyzer{p}  // ,
	// More initializers here.
	{

		// Call appropriate consumes<>() for any products to be retrieved by this module.
		fhicl::ParameterSet pset_blipalg = p.get<fhicl::ParameterSet>("BlipAlg");
		fBlipAlg = new blip::BlipRecoAlg( pset_blipalg );
	}

	void BlipVertexAnalysis::analyze(art::Event const& e)
	{//Main code
		ClearVars(vars);
	
		// Implementation of required member function here.
		vars.evt = e.id().event();
		vars.run = e.id().run();

		//CHECK hardcode this
		std::string m_pandoraLabel = "pandora";
		std::string m_trackLabel = "pandora";
		std::string m_showerLabel = "pandora";

		//Some verticies.········
		//Collect the PFParticles from the event. 
		art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = e.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
		std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;//Prepare this vector
		art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

		//load vertices with PFParticle key()
		art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, e, m_pandoraLabel);

		//take neutrino's vertex
		int NeutrinoIndex = -1;
		int NuCounts = 0;
		for(size_t i=0; i< pfParticleVector.size(); ++i){
			auto pfp = pfParticleVector[i];
			// Check if this particle is identified as the neutrino
			const int pdg(pfp->PdgCode());
			const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
			if(isNeutrino){
			NeutrinoIndex = i;
			NuCounts++;
			}
		}
		
		if( NuCounts ==1 ){
			const lar_pandora::VertexVector &vertexVector = vertices_per_pfparticle.at(pfParticleVector[NeutrinoIndex].key());
			GrabPandoraVertex(vertexVector, vars);
		} else {

		std::cout <<"Warning, number of neutrino candidates in this event is "<<NuCounts<<std::endl;
		}

		//CHECK Hard code a vertex
		NeutrinoIndex = 0;
		const lar_pandora::VertexVector &vertexVector = vertices_per_pfparticle.at(pfParticleVector[NeutrinoIndex].key());
		GrabPandoraVertex(vertexVector, vars);
		//CHECK hardcode finished
		
		//Now Track and Shower
		std::vector<art::Ptr<recob::PFParticle>> pfNuDaughtersVector;
		std::vector<art::Ptr<recob::PFParticle>> pfpTrackVector;
		std::vector<art::Ptr<recob::PFParticle>> pfpShowerVector;
		std::vector<art::Ptr<recob::PFParticle>> pfpNIDVector;

		for (const size_t daughterId : pfParticleVector[NeutrinoIndex]->Daughters()){
			pfNuDaughtersVector.push_back( findpfpDaughter(pfParticleVector, daughterId));
		}


		art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, e, m_trackLabel);
		std::vector< art::Ptr<recob::Track> > TrackVector;

		art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, e, m_trackLabel);
		std::vector< art::Ptr<recob::Shower> > ShowerVector;

		for(art::Ptr<recob::PFParticle> tmp_pfp : pfNuDaughtersVector){
			std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(tmp_pfp.key()));
			if(associatedTracks.size() == 1){
				pfpTrackVector.push_back(tmp_pfp);
				TrackVector.push_back( associatedTracks.front() );
				continue;
			}

			std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(tmp_pfp.key()));
			if(associatedShowers.size() == 1){
				pfpShowerVector.push_back(tmp_pfp);
				ShowerVector.push_back( associatedShowers.front() );
				continue;
			}
			pfpNIDVector.push_back(tmp_pfp);
		}

		vars.m_num_tracks = TrackVector.size();
		vars.m_num_showers = ShowerVector.size();
		std::cout<<"CHECK Found pfp-track:"<<vars.m_num_tracks <<std::endl;
		std::cout<<"CHECK Found pfp-shower:"<<vars.m_num_showers <<std::endl;
		std::cout<<"CHECK Found pfp-empty:"<<pfpNIDVector.size()<<std::endl;

		for(size_t index = 0; index < ShowerVector.size(); index ++){
			AnalyzeShower(pfpShowerVector[index], ShowerVector[index]);
		}

		//BLIP: grab blip varaibles;
		fBlipAlg->RunBlipReco(e);
		AnalyzeBlips( fBlipAlg, vars);

		AnalyzeBlipsNearVertex( vars );



		vars.f_output_tree->Fill();
	}

	void BlipVertexAnalysis::endSubRun(const art::SubRun& sr)
	{//wrap up, fill up event info. here;

		ClearVarsEvt(vars_evt);
		vars_evt._totpot_run    = sr.run();
		vars_evt._totpot_subrun = sr.subRun();
		vars_evt._begintime     = sr.beginTime().value();
		vars_evt._endtime       = sr.endTime().value();
		art::Handle< sumdata::POTSummary > potListHandle;
		if(sr.getByLabel(fPOTModuleLabel,potListHandle)) {
			vars_evt._totpot = potListHandle->totpot;
		} else {
			vars_evt._totpot = 0.;
		}


		vars_evt.fPOT->Fill();

	}

	void BlipVertexAnalysis::beginJob()
	{
		// Implementation of optional member function here.
		art::ServiceHandle<art::TFileService> tfs;
		vars.f_output_tree = tfs->make<TTree>("BlipVertexAnalysis","BlipVertexAnalysis");
		vars_evt.fPOT = tfs->make<TTree>("pottree","pot tree");

		std::cout<<"CHECK running customized package"<<std::endl;
		CreateVarsBranches(vars);
		CreateVarsEvtBranches(vars_evt);
	}

	void BlipVertexAnalysis::endJob()
	{
		// Implementation of optional member function here.
	}

	DEFINE_ART_MODULE(BlipVertexAnalysis)
}
