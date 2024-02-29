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
			double tmp_dist = f_diff.Mag();

			vars.sps_dist.push_back( tmp_dist);
			std::cout<<" Dist "<<f_diff.Mag()<<std::endl;
		}

		//look at the distribution;
		vars.sps_dist_sorted = vars.sps_dist;
		std::sort(vars.sps_dist_sorted.begin(), vars.sps_dist_sorted.end()); //ascending order
		
		double ringRadius = 10.0;
		double max_val = 1050;
		int numRings = max_val / ringRadius;
		
		std::vector<int> ringCounts(numRings, 0);
		std::vector<double> ringSums(numRings, 0.0);

		for(double& tmp_dist : vars.sps_dist_sorted){
			int ringIndex = static_cast<int>(tmp_dist / ringRadius);

			ringCounts[ringIndex]++;
			ringSums[ringIndex] += tmp_dist;

		}
		for (int i = 0; i < numRings; ++i) {
			//double ringStart = static_cast<double>(i);
			double ringEnd = static_cast<double>(i + 1);

			if (ringEnd <= max_val) {
				double mean = ringCounts[i] > 0 ? ringSums[i] / ringCounts[i] : 0.0;

//				std::cout << "Bin " << i << " (" << ringStart << " to " << ringEnd << "): ";
//				std::cout << "Count: " << ringCounts[i] << ", Mean: " << mean << "\n";
				vars.sps_dist_mean10cmrings.push_back(mean);
			}
		}
		vars.sps_counts_10cmrings = ringCounts;




	}//Finish AnalyzeBlipsNearVertex
}


#endif
