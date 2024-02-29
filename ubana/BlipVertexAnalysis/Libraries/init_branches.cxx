#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_INIT_BRANCHES_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_INIT_BRANCHES_CXX

#include "init_branches.h"

namespace BVA_ana
{
	void ClearVars(var_all& vars){

		vars.run = 0;
		vars.evt=0;
		vars.nsps = 0;
		vars.max_plane2_integral = 0;
		vars.max_plane2_true_integral = 0;
		vars.plane2_integrals->clear();
		vars.plane2_true_integrals->clear();
		//true_pdg = 0.;
		//true_pdg_v->clear();
		vars.sps_x->clear();
		vars.sps_y->clear();
		vars.sps_z->clear();
		vars.elec_x->clear();
		vars.elec_y->clear();
		vars.elec_z->clear();
		vars.elec_E->clear();
		vars.nhits_pl2->clear();
		vars.nhits_oth_max->clear();
		vars.nhits_oth_min->clear();
		vars.plane_oth_max_integrals->clear();
		vars.plane_oth_min_integrals->clear();
		vars.dep_electrons->clear();
		vars.g4_energy->clear();
		vars.g4_x->clear();
		vars.g4_y->clear();
		vars.g4_z->clear();
		vars.matched_plane->clear();
		vars.ntrueelec = 0;

		vars.m_reco_vertex_size = 0;
		vars.m_vertex_pos_x=-99999;
		vars.m_vertex_pos_y=-99999;
		vars.m_vertex_pos_z=-99999;
		vars.m_num_tracks = 0;
		vars.m_num_showers = 0;

		vars.sps_dist->clear();

	}

	void CreateVarsBranches(var_all& vars){

		vars.f_output_tree->Branch("run",vars.run);
		vars.f_output_tree->Branch("evt",vars.evt);
		vars.f_output_tree->Branch("nsps",vars.nsps);
		vars.f_output_tree->Branch("sps_x",vars.sps_x);
		vars.f_output_tree->Branch("sps_y",vars.sps_y);
		vars.f_output_tree->Branch("sps_z",vars.sps_z);
		vars.f_output_tree->Branch("g4_E",vars.g4_energy);
		vars.f_output_tree->Branch("g4_x",vars.g4_x);
		vars.f_output_tree->Branch("g4_y",vars.g4_y);
		vars.f_output_tree->Branch("g4_z",vars.g4_z);
		vars.f_output_tree->Branch("max_pl2_integ",vars.max_plane2_integral);
		vars.f_output_tree->Branch("max_pl2_true_integ",vars.max_plane2_true_integral);
		vars.f_output_tree->Branch("pl2_integs",vars.plane2_integrals);
		vars.f_output_tree->Branch("pl2_true_integs",vars.plane2_true_integrals);
		vars.f_output_tree->Branch("pl2_nhits",vars.nhits_pl2);
		vars.f_output_tree->Branch("pl_othmax_nhits",vars.nhits_oth_max);
		vars.f_output_tree->Branch("pl2_othmin_nhits",vars.nhits_oth_min);
		vars.f_output_tree->Branch("pl2_othmax_integs",vars.plane_oth_max_integrals);
		vars.f_output_tree->Branch("pl2_othmin_integs",vars.plane_oth_min_integrals);
		vars.f_output_tree->Branch("dep_electrons",vars.dep_electrons);
		vars.f_output_tree->Branch("induction_plane_max_integral",vars.matched_plane);
		vars.f_output_tree->Branch("nelec",vars.ntrueelec);
		vars.f_output_tree->Branch("elec_x",vars.elec_x);
		vars.f_output_tree->Branch("elec_y",vars.elec_y);
		vars.f_output_tree->Branch("elec_z",vars.elec_z);
		vars.f_output_tree->Branch("elec_E",vars.elec_E);

		//add_reco_vertex
		vars.f_output_tree->Branch("reco_vertex_size", &vars.m_reco_vertex_size, "reco_vertex_size/I");
		vars.f_output_tree->Branch("reco_vertex_x", &vars.m_vertex_pos_x, "reco_vertex_x/D");
		vars.f_output_tree->Branch("reco_vertex_y", &vars.m_vertex_pos_y, "reco_vertex_y/D");
		vars.f_output_tree->Branch("reco_vertex_z", &vars.m_vertex_pos_z, "reco_vertex_z/D");
		vars.f_output_tree->Branch("reco_asso_tracks",  &vars.m_num_tracks ,"reco_asso_tracks/I");
		vars.f_output_tree->Branch("reco_asso_showers", &vars.m_num_showers,"reco_asso_showers/I");
		
		//add vertex-blips
		vars.f_output_tree->Branch("sps_dist",vars.sps_dist);
	}

	void ClearVarsEvt(var_evt& vars){
		vars._totpot = 0;
		vars._totpot_run = 0;
		vars._totpot_subrun = 0;
		vars._begintime = 0;
		vars._endtime =0 ;

	}

	void CreateVarsEvtBranches(var_evt& vars){
		vars.fPOT->Branch("totpot",vars._totpot);
		vars.fPOT->Branch("run",vars._totpot_run);
		vars.fPOT->Branch("subrun",vars._totpot_subrun);
		vars.fPOT->Branch("begintime",vars._begintime);
		vars.fPOT->Branch("endtime",vars._endtime);

	}



}

#endif
