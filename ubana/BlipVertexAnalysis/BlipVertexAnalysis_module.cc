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
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

// extra includes, trying to fix
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"

// blip reco includes
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// for MCParticle.E() and such
#include "nusimdata/SimulationBase/MCTrajectory.h"
// for associations
#include "canvas/Persistency/Common/FindManyP.h" // find associations as pointers
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "lardataobj/RecoBase/Track.h"
// access hits and c++ stuff
#include "lardataobj/RecoBase/Hit.h"
#include <vector> // to store multiple things per event
#include <string> // for product names
#include <set>
#include <string>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"

// this line good for "new" larsoft:
#include "art/Persistency/Common/PtrMaker.h"
// for outdated versions (i.e MCC8) use this line:
//#include "lardata/Utilities/PtrMaker.h"

// additional framework includes
//#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Optional/TFileService.h"

// pot summary
// from /uboone/app/users/afurmans/muonScatteringAnalysis/
//       larsoft/srcs/uboonecode/uboone/AnalysisTree/
//       AnalysisTree_module.cc
#include "larcoreobj/SummaryData/POTSummary.h"


namespace BVA_ana
{
	class BlipVertexAnalysis;


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


		private:

			// Declare member data here.
			// blip
			blip::BlipRecoAlg* fBlipAlg;

			std::string fPOTModuleLabel;
			TTree * fPOT;
			Double_t _totpot;
			Double_t _totpot_run;
			Double_t _totpot_subrun;
			Double_t _begintime;
			Double_t _endtime;
			TTree * f_output_tree; // output tree
			int run;
			int evt;
			int nsps;
			double max_plane2_integral;
			double max_plane2_true_integral;
			std::vector<double> *plane2_integrals = new std::vector<double>;
			std::vector<double> *plane2_true_integrals = new std::vector<double>;
			std::vector<int> *nhits_pl2 = new std::vector<int>;
			std::vector<int> *nhits_oth_max = new std::vector<int>;
			std::vector<int> *nhits_oth_min = new std::vector<int>;
			std::vector<double> *plane_oth_max_integrals = new std::vector<double>;
			std::vector<double> *plane_oth_min_integrals = new std::vector<double>;
			std::vector<double> *sps_x = new std::vector<double>;
			std::vector<double> *sps_y = new std::vector<double>;
			std::vector<double> *sps_z = new std::vector<double>;
			std::vector<double> *elec_x = new std::vector<double>;
			std::vector<double> *elec_y = new std::vector<double>;
			std::vector<double> *elec_z = new std::vector<double>;
			std::vector<double> *elec_E = new std::vector<double>;
			std::vector<int> *dep_electrons = new std::vector<int>;
			std::vector<double> *g4_energy = new std::vector<double>;
			std::vector<double> *g4_x = new std::vector<double>;
			std::vector<double> *g4_y = new std::vector<double>;
			std::vector<double> *g4_z = new std::vector<double>;
			std::vector<int> *matched_plane = new std::vector<int>;
			int ntrueelec;

			//add reco vertex
			m_reco_vertex_size = 0;
			m_vertex_pos_x=-99999;
			m_vertex_pos_y=-99999;
			m_vertex_pos_z=-99999;

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
	{
		run = 0;
		evt=0;
		nsps = 0;
		max_plane2_integral = 0;
		max_plane2_true_integral = 0;
		plane2_integrals->clear();
		plane2_true_integrals->clear();
		//true_pdg = 0.;
		//true_pdg_v->clear();
		sps_x->clear();
		sps_y->clear();
		sps_z->clear();
		elec_x->clear();
		elec_y->clear();
		elec_z->clear();
		elec_E->clear();
		nhits_pl2->clear();
		nhits_oth_max->clear();
		nhits_oth_min->clear();
		plane_oth_max_integrals->clear();
		plane_oth_min_integrals->clear();
		dep_electrons->clear();
		g4_energy->clear();
		g4_x->clear();
		g4_y->clear();
		g4_z->clear();
		matched_plane->clear();
		ntrueelec = 0;

		// Implementation of required member function here.
		evt = e.id().event();
		run = e.id().run();

		// run blip reconstruction
		fBlipAlg->RunBlipReco(e);
		std::vector<blip::Blip> blip_v      = fBlipAlg->blips;     // reco blip vector
		std::vector<blip::TrueBlip> true_blip_v = fBlipAlg->trueblips; // true blip vector
		std::vector<blip::ParticleInfo> &pinfo = fBlipAlg->pinfo;

		// truth primary electrons
		// this is what the generator made
		for ( auto p: pinfo ) {
			double true_x = -999.;
			double true_y = -999.;
			double true_z = -999.;
			double true_E = -999.;
			if ( p.isPrimary == 1 ) {
				ntrueelec++;
				true_x = p.position.X();
				true_y = p.position.Y();
				true_z = p.position.Z();
				true_E = p.KE/1000.;
				elec_x->push_back(true_x);
				elec_y->push_back(true_y);
				elec_z->push_back(true_z);
				elec_E->push_back(true_E);
			}
		}

		// loop over reco blips
		for ( auto iblip: blip_v ) {
			// blip position
			sps_x->push_back(iblip.X());
			sps_y->push_back(iblip.Y());
			sps_z->push_back(iblip.Z());

			/**
			 * truth matched reco blips go to "g4" branches
			 * branch in sync with regular blip branch
			 */
			auto tb = iblip.truth;   // truth info of current blip
			auto c = iblip.clusters; // clusters of the current blip 
			double true_x = -999.;
			double true_y = -999.;
			double true_z = -999.;
			double true_E = -999.;
			double true_integral[3] = { -999., -999., -999.0};
			int dep_e = -999;
			if ( tb.ID != -9 ) {
				true_x = tb.Position.X();
				true_y = tb.Position.Y();
				true_z = tb.Position.Z();
				true_E = tb.Energy/1000.; // blipreco Energy is kinetic energy
				dep_e = tb.DepElectrons;
				for ( int i = 0; i < 3; i++ ) {
					true_integral[i] = c[i].ADCs;
				}
			}
			g4_x->push_back(true_x);
			g4_y->push_back(true_y);
			g4_z->push_back(true_z);
			g4_energy->push_back(true_E);
			dep_electrons->push_back(dep_e); // deposited electrons
			plane2_true_integrals->push_back(true_integral[2]);
			if ( true_integral[2] > max_plane2_true_integral ) {
				max_plane2_true_integral = true_integral[2];
			}

			/**
			 * loop over clusters

			 * 3-plane matched blip is made of three clusters, one in each
			 * wire plane
			 */
			int nhits[3] = {0};         // nhits in each plane
			double integrals[3] = {0.}; // integral in each plane (ADCs)
			for ( int plane = 0; plane < 3; plane++ ) {
				nhits[plane]     = c[plane].NHits;
				//integrals[plane] = c[plane].Charge;
				integrals[plane] = c[plane].ADCs;
				if ( plane == 2 ) {
					if ( integrals[2] > max_plane2_integral ) {
						max_plane2_integral = integrals[2];
					}
					plane2_integrals->push_back(integrals[2]);
					nhits_pl2->push_back(nhits[2]);
				}
			}

			/**

			 * classify induction plane according to max/min integrals.
			 * save nhits and integrals
			 */
			if ( integrals[0] > integrals[1] ) {
				nhits_oth_max->push_back(nhits[0]);
				plane_oth_max_integrals->push_back(integrals[0]);
				nhits_oth_min->push_back(nhits[1]);
				plane_oth_min_integrals->push_back(integrals[1]);
				matched_plane->push_back(0);
			} else {
				nhits_oth_max->push_back(nhits[1]);
				plane_oth_max_integrals->push_back(integrals[1]);
				nhits_oth_min->push_back(nhits[0]);
				plane_oth_min_integrals->push_back(integrals[0]);
				matched_plane->push_back(1);
			}
		}

		nsps = blip_v.size();
		// store outputs in the tree
		f_output_tree->Fill();
	}

	void BlipVertexAnalysis::endSubRun(const art::SubRun& sr)
	{
		_totpot = 0;
		_totpot_run = 0;
		_totpot_subrun = 0;
		_begintime = 0;
		_endtime =0 ;

		_totpot_run    = sr.run();
		_totpot_subrun = sr.subRun();
		_begintime     = sr.beginTime().value();
		_endtime       = sr.endTime().value();

		art::Handle< sumdata::POTSummary > potListHandle;
		if(sr.getByLabel(fPOTModuleLabel,potListHandle)) {
			_totpot = potListHandle->totpot;
		} else {
			_totpot = 0.;
		}
		fPOT->Fill();
	}

	void BlipVertexAnalysis::beginJob()
	{
		// Implementation of optional member function here.
		art::ServiceHandle<art::TFileService> tfs;
		std::cout<<"CHECK running customized package"<<std::endl;
		fPOT = tfs->make<TTree>("pottree","pot tree");
		fPOT->Branch("totpot",&_totpot);
		fPOT->Branch("run",&_totpot_run);
		fPOT->Branch("subrun",&_totpot_subrun);
		fPOT->Branch("begintime",&_begintime);
		fPOT->Branch("endtime",&_endtime);

		f_output_tree = tfs->make<TTree>("BlipVertexAnalysis","BlipVertexAnalysis");
		f_output_tree->Branch("run",&run);
		f_output_tree->Branch("evt",&evt);
		f_output_tree->Branch("nsps",&nsps);
		f_output_tree->Branch("sps_x",&sps_x);
		f_output_tree->Branch("sps_y",&sps_y);
		f_output_tree->Branch("sps_z",&sps_z);
		f_output_tree->Branch("g4_E",&g4_energy);
		f_output_tree->Branch("g4_x",&g4_x);
		f_output_tree->Branch("g4_y",&g4_y);
		f_output_tree->Branch("g4_z",&g4_z);
		f_output_tree->Branch("max_pl2_integ",&max_plane2_integral);
		f_output_tree->Branch("max_pl2_true_integ",&max_plane2_true_integral);
		f_output_tree->Branch("pl2_integs",&plane2_integrals);
		f_output_tree->Branch("pl2_true_integs",&plane2_true_integrals);
		f_output_tree->Branch("pl2_nhits",&nhits_pl2);
		f_output_tree->Branch("pl_othmax_nhits",&nhits_oth_max);
		f_output_tree->Branch("pl2_othmin_nhits",&nhits_oth_min);
		f_output_tree->Branch("pl2_othmax_integs",&plane_oth_max_integrals);
		f_output_tree->Branch("pl2_othmin_integs",&plane_oth_min_integrals);
		f_output_tree->Branch("dep_electrons",&dep_electrons);
		f_output_tree->Branch("induction_plane_max_integral",&matched_plane);
		f_output_tree->Branch("nelec",&ntrueelec);
		f_output_tree->Branch("elec_x",&elec_x);
		f_output_tree->Branch("elec_y",&elec_y);
		f_output_tree->Branch("elec_z",&elec_z);
		f_output_tree->Branch("elec_E",&elec_E);

		//add_reco_vertex
		f_output_tree->Branch("reco_vertex_size", &m_reco_vertex_size);
		f_output_tree->Branch("reco_vertex_x", &m_vertex_pos_x);
		f_output_tree->Branch("reco_vertex_y", &m_vertex_pos_y);
		f_output_tree->Branch("reco_vertex_z", &m_vertex_pos_z);

	}

	void BlipVertexAnalysis::endJob()
	{
		// Implementation of optional member function here.
	}
}

DEFINE_ART_MODULE(BlipVertexAnalysis)
