#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_CXX

#include "BVA_blipsnearvertex.h"

namespace BVA_ana
{
	void AnalyzeBlipsNearVertex(var_all& vars){

		TVector3 f_vertex(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z);
	
		//Use an unordered_map to store blips according to the distance in a histogram.
		std::unordered_map<int, std::vector<int>> groups_dist; //map(distance index , sps_index)
		double shellRadius = 30.0;
		double max_val = 1050;
		int numRings = max_val / shellRadius;//The last bin will be overflow
		bool verbose = false;

//		TVector3 beamZ;
//		beamZ = (0,0,1); //on-axis experiment

		//Loop over blips, and evaluate their distance to vertex
		for( size_t index = 0 ; index < (vars.sps_x)->size(); index++ ){
			TVector3 f_bliplocation( vars.sps_x->at(index),vars.sps_y->at(index),vars.sps_z->at(index));
			TVector3 f_diff;
			
			//geometric distance
			f_diff = f_vertex - f_bliplocation;
			double tmp_dist = f_diff.Mag();

			vars.sps_dist.push_back( tmp_dist);

			
			//Build the groups_dist map here based on blip distance
			int groupKey = static_cast<int>(tmp_dist / shellRadius);
			groups_dist[groupKey].push_back(index);


			//geometric angular relation
			double tmpTheta = GetThetaYZDegree	(f_vertex, f_bliplocation);
			double tmpPhi = GetPhiXZDegree		(f_vertex, f_bliplocation);
			vars.sps_ang_thetaYZdegrees.push_back( tmpTheta);
			vars.sps_ang_phiXZdegrees.push_back( tmpPhi);

			if(verbose){
			std::cout<<" Dist "<<f_diff.Mag()<<std::endl;
			std::cout<<" Theta "<<tmpTheta<<std::endl;
			std::cout<<" Phi "<<tmpPhi<<std::endl;
			}
		}

		// Iterate the groups_dist and evaluate other quantities;
		for (const auto& pair : groups_dist) {
			std::cout << "Group " << pair.first ;
			for (int index : pair.second) {
				std::cout << index << " ";
			}
			std::cout << std::endl;
		}

		SummarizeDist( vars, groups_dist);

		std::vector<int> shellCounts(numRings, 0);
		std::vector<double> shellSums(numRings, 0.0);

		//sort out the distance from small to large
		vars.sps_dist_sorted = vars.sps_dist;
		std::sort(vars.sps_dist_sorted.begin(), vars.sps_dist_sorted.end()); //ascending order
		
		


	}//End of AnalyzeBlipsNearVertex()


	//modify vars.sps_dist_mean30cmEachShells and vars.sps_counts_30cmEachShells;
	//Each group contains blips within concentric spheres with steps of shellRadius
	void SummarizeDist(var_all& vars, std::unordered_map<int, std::vector<int>>groups_dist){

		bool verbose=true;
		for (const auto& pair : groups_dist) {
			if(verbose) std::cout << "Group " << pair.first ;

			int tmp_counts = (pair.second).size();
			vars.sps_counts_30cmEachShells = tmp_counts;

			double tmp_mean = 0;

			for (int index : pair.second) {
				if(verbose)	std::cout << index << " " << " value: "<<vars.sps_dist[index];
				tmp_mean += vars.sps_dist[index];
			}
			if(verbose)std::cout <<std::endl;

			tmp_mean /= tmp_counts;
			vars.sps_dist_mean30cmEachShells.push_back(tmp_mean);
		}

	}//End of SummarizeDist

	double GetThetaYZDegree(TVector3 vertex, TVector3 blip){
				TVector3 direction = (blip - vertex).Unit();
				double tan = direction.Z() / direction.Y();
//				std::cout<<"Direction "<<direction.Z()<<" and "<<direction.Y()<<std::endl;
//				std::cout<<"CHECK tan"<<tan<<std::endl;
				double degrees = atan(tan) /3.14159*180;
				return  degrees;
	}

	double GetPhiXZDegree(TVector3 vertex, TVector3 blip){
				TVector3 direction = (blip - vertex).Unit();
				double tan = direction.X() / direction.Y();
				double degrees = atan(tan) /3.14159*180;
				return  degrees;

	}

//	void SummarizeAngInRad(var_all& vars, TVector3 vertex, std::vector<double>groups_dist){
//
//		bool verbose=true;
//
//		vars.sps_ang_unitthetaYZ.resize( (vars.sps_x)->size() );
//		vars.sps_ang_unitphiXZ.resize( (vars.sps_x)->size() );
//
//		for (const auto& pair : groups_dist) {
//			if(verbose)std::cout << "Group " << pair.first ;
//
//			int tmp_counts = vars.sps_counts_30cmEachShells[pair.first];
//
//			if(vars.sps_counts_30cmEachShells.size() < 1){
//				tmp_counts = (pair.second).size();
//				vars.sps_counts_30cmEachShells.resize(tmp_counts);
//			}
//
//			for (int index : pair.second) {
//				if(verbose)	std::cout << index << " " << " value: "<<vars.sps_dist[index];
//
//				TVector3 f_bliplocation( vars.sps_x->at(index),vars.sps_y->at(index),vars.sps_z->at(index));
//				TVector3 direction = (f_bliplocation - vertex).Unit();
////				direction /= direction.Mag();
//				vars.sps_ang_unitthetaYZ[index]=( acos(direction.Z() / direction.Y()) );
//				vars.sps_ang_unitphiXZ[index]=( acos(direction.Z() / direction.X()) );
//
//				if(verbose) std::cout<<" YZ " <<vars.sps_ang_unitthetaYZ[index] <<" XZ "<< vars.sps_ang_unitphiXZ[index] <<"; ";
//			}
//
//			if(verbose)std::cout <<std::endl;
//
//		}
//
//	}//End of SummarizeAngInRad



}


#endif
