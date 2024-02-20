#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_VARIABLES_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_VARIABLES_H

// ROOT includes
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

// access hits and c++ stuff
#include "lardataobj/RecoBase/Hit.h"
#include <vector> // to store multiple things per event
#include <string> // for product names
#include <set>
#include <string>
#include <iostream>
#include "TVector3.h"



namespace BVA_ana
{

	struct var_evt{
			TTree * fPOT;
			Double_t _totpot;
			Double_t _totpot_run;
			Double_t _totpot_subrun;
			Double_t _begintime;
			Double_t _endtime;
	};

	struct var_all{
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
			int m_reco_vertex_size;
			double m_vertex_pos_x;
			double m_vertex_pos_y;
			double m_vertex_pos_z;

			double true_x;
			double true_y;
			double true_z;
			double true_E;
			double dep_e;

	};



}




#endif
