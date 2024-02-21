#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX

#include "BVA_blipsnearvertex.h"

namespace BVA_ana
{
	void AnalyzeBlipsNearVertex(var_all& vars){
	
		TVector3 f_vertex(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z);

		for( size_t index = 0 ; index < (vars.sps_x)->size(); index++ ){
			TVector3 f_bliplocation( vars.sps_x->at(index),vars.sps_y->at(index),vars.sps_z->at(index));
			TVector3 f_diff;
			f_diff = f_vertex - f_bliplocation;

			vars.sps_dist->push_back( f_diff.Mag());
			std::cout<<" Dist "<<f_diff.Mag()<<std::endl;
		}


	}
}


#endif
