#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX

#include "BVA_blipsnearvertex.h"

namespace BVA_ana
{
	void AnalyzeBlipsNearVertex(var_all& vars){

		TVector3 f_vertex(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z);
	
		//Use an unordered_map to store blips according to the distance in a histogram.
		unordered_map<int, vector<int>> groups_dist; //map(distance index , sps_index)
		double ringRadius = 30.0;
		double max_val = 1050;
		int numRings = max_val / ringRadius;//The last bin will be overflow


		//Loop over blips, and evaluate their distance to vertex
		for( size_t index = 0 ; index < (vars.sps_x)->size(); index++ ){
			TVector3 f_bliplocation( vars.sps_x->at(index),vars.sps_y->at(index),vars.sps_z->at(index));
			TVector3 f_diff;
			
			//geometric distance
			f_diff = f_vertex - f_bliplocation;
			double tmp_dist = f_diff.Mag();

			vars.sps_dist.push_back( tmp_dist);
			std::cout<<" Dist "<<f_diff.Mag()<<std::endl;

			
			//Build the groups_dist map here based on blip distance
			int groupKey = static_cast<int>(tmp_dist / ringRadius);
			groups_dist[groupKey].push_back(index);
		}

		// Iterate the groups_dist and evaluate other quantities;
		for (const auto& pair : groups_dist) {
			cout << "Group " << pair.first ;
			for (int index : pair.second) {
				cout << index << " ";
			}
			cout << endl;
		}

		SummarizeDist( vars, groups_dist);

		std::vector<int> ringCounts(numRings, 0);
		std::vector<double> ringSums(numRings, 0.0);

		//sort out the distance from small to large
		vars.sps_dist_sorted = vars.sps_dist;
		std::sort(vars.sps_dist_sorted.begin(), vars.sps_dist_sorted.end()); //ascending order
		
		


	}//End of AnalyzeBlipsNearVertex()


	//modify vars.sps_dist_mean30cmrings and vars.sps_counts_30cmrings;
	void SummarizeDist(var_all& vars, unordered_map<int, vector<int>>groups_dist){

		bool verbose=true;
		for (const auto& pair : groups_dist) {
			if(verbose) cout << "Group " << pair.first ;

			int tmp_counts = (pair.second).size();
			vars.sps_counts_30cmrings = tmp_counts;

			double tmp_mean = 0;

			for (int index : pair.second) {
				if(verbose)	cout << index << " " << " value: "<<sps_dist[index];
				tmp_mean += sps_dist[index];
			}
			if(verbose) cout << endl;

			tmp_mean /= tmp_counts;
			sps_dist_mean30cmrings.push_back(tmp_mean);
		}

	}//End of SummarizeDist

	void SummarizeAngInRad(var_all& vars, TVector3 vertex, unordered_map<int, vector<int>>groups_dist){

		bool verbose=true;

		sps_ang_unitthetaYZ.resize( (vars.sps_x)->size() );
		sps_ang_unitphiYZ.resize( (vars.sps_x)->size() );

		for (const auto& pair : groups_dist) {
			if(verbose) cout << "Group " << pair.first ;

			int tmp_counts = vars.sps_counts_30cmrings[pair.first];

			if(vars.sps_counts_30cmrings.size() < 1){
				tmp_counts = (pair.second).size();
				vars.sps_counts_30cmrings = tmp_counts;
			}

			double tmp_mean = 0;

			for (int index : pair.second) {
				if(verbose)	std::cout << index << " " << " value: "<<sps_dist[index];

				TVector3 f_bliplocation( vars.sps_x->at(index),vars.sps_y->at(index),vars.sps_z->at(index));
				TVector3 direction = f_bliplocation - vertex;
				direction /= direction.Mag();
				sps_ang_unitthetaYZ[index]=( acos(direction.Z() / direction.Y()) );
				sps_ang_unitphiXZ[index]=( acos(direction.Z() / direction.X()) );

				if(verbose) std::cout<<" YZ " <<sps_ang_unitthetaYZ[index] <<" XZ "<< sps_ang_unitphiXZ[index] <<"; ";
			}

			if(verbose) cout << endl;

		}

	}//End of SummarizeAngInRad



}


#endif
