#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYPANDORA_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYPANDORA_CXX
#include "BVA_onlypandora.h"

namespace BVA_ana
{
    void AnalyzePandora(lar_pandora::VertexVector vertexVector, var_all& vars){


		vars.m_reco_vertex_size = vertexVector.size();
		std::cout<<"CHECK number of vertex "<<vars.m_reco_vertex_size<<std::endl;

		if (!vertexVector.empty())
		{
                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);
				vars.m_vertex_pos_x = xyz[0];
				vars.m_vertex_pos_y = xyz[1];
				vars.m_vertex_pos_z = xyz[2];

				std::cout<<"CHECK vertex location ("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")"<<std::endl;

		}
    }
}








#endif
