#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYBLIPS_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYBLIPS_CXX

#include "BVA_onlyblips.h"

namespace BVA_ana
{
	void AnalyzeBlips(blip::BlipRecoAlg* fBlipAlg, var_all& vars){
		// run blip reconstruction
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
				vars.ntrueelec++;
				vars.true_x = p.position.X();
				vars.true_y = p.position.Y();
				vars.true_z = p.position.Z();
				vars.true_E = p.KE/1000.;
				vars.elec_x->push_back(true_x);
				vars.elec_y->push_back(true_y);
				vars.elec_z->push_back(true_z);
				vars.elec_E->push_back(true_E);
			}
		}

		// loop over reco blips
		for ( auto iblip: blip_v ) {
			// blip position
			vars.sps_x->push_back(iblip.X());
			vars.sps_y->push_back(iblip.Y());
			vars.sps_z->push_back(iblip.Z());

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
				vars.true_x = tb.Position.X();
				vars.true_y = tb.Position.Y();
				vars.true_z = tb.Position.Z();
				vars.true_E = tb.Energy/1000.; // blipreco Energy is kinetic energy
				vars.dep_e = tb.DepElectrons;
				for ( int i = 0; i < 3; i++ ) {
					true_integral[i] = c[i].ADCs;
				}
			}
			vars.g4_x->push_back(true_x);
			vars.g4_y->push_back(true_y);
			vars.g4_z->push_back(true_z);
			vars.g4_energy->push_back(true_E);
			vars.dep_electrons->push_back(dep_e); // deposited electrons
			vars.plane2_true_integrals->push_back(true_integral[2]);
			if ( true_integral[2] > vars.max_plane2_true_integral ) {
				vars.max_plane2_true_integral = true_integral[2];
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
					if ( integrals[2] > vars.max_plane2_integral ) {
						vars.max_plane2_integral = integrals[2];
					}
					vars.plane2_integrals->push_back(integrals[2]);
					vars.nhits_pl2->push_back(nhits[2]);
				}
			}

			/**

			 * classify induction plane according to max/min integrals.
			 * save nhits and integrals
			 */
			if ( integrals[0] > integrals[1] ) {
				vars.nhits_oth_max->push_back(nhits[0]);
				vars.plane_oth_max_integrals->push_back(integrals[0]);
				vars.nhits_oth_min->push_back(nhits[1]);
				vars.plane_oth_min_integrals->push_back(integrals[1]);
				vars.matched_plane->push_back(0);
			} else {
				vars.nhits_oth_max->push_back(nhits[1]);
				vars.plane_oth_max_integrals->push_back(integrals[1]);
				vars.nhits_oth_min->push_back(nhits[0]);
				vars.plane_oth_min_integrals->push_back(integrals[0]);
				vars.matched_plane->push_back(1);
			}
		}

		vars.nsps = blip_v.size();
		// store outputs in the tree


	}
}


#endif
