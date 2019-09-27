

#ifndef __BOBBYVERTEXBUILDER_H__
#define __BOBBYVERTEXBUILDER_H__

#include "SinglePhoton_module.h"
#include "Atlas.h"

//---VertexBuilder---
#include "VertexBuilder/VertexBuilder.h"
#include "VertexBuilder/ParticleAssociations.h" 
//#include "VertexBuilder/DetectorObjects.h"
//
//#include "VertexQuality.h"
//#include "FillTreeVariables.h"
//#include "RecoMCMatching.h"
//-------------------

//#include "TCanvas.h"
//#include "TGraph.h"
//#include "TFile.h"
//#include "TAxis.h"
//#include "TLine.h"
//#include "TLatex.h"
//#include "TLegend.h"
//#include "TPrincipal.h"
//#include "TVectorD.h"
//#include "TMatrixD.h"
//#include "TF1.h"
//#include "TEllipse.h"

//the following is not working yet
//#include <curses.h> //Use LINES and COLS to get the dimensions of the screen

using namespace std;

namespace single_photon
{

	/*****************************
	 * CollectTracksAndShowers_v2() - Fill in package.selected_showers and
	 *	package.selected_tracks as well as other nu-like objects (more_tracks/showers)
	 *
	 *  packahe (to be modified) - contains all vectors and maps that are
	 *		needed to access pandora_objects; it will be modified as
	 *		the code runs.
	 *  //not using run_count here now;
	 *	run_count (input) - has values 0,1,2; each number indicates differnt collections of slices: 
	 *		0 - selected nu_slice; 
	 *		1 - not-selected nu_slices; 
	 *		2 - cosmic slices (non-nu_slices).  NOT DOING THIS
	 *
	 * **************************/

	void SinglePhoton::CollectTracksAndShowers_v2(
			const art::Event &evt,
			class Atlas &package){
		//			int run_count){
		//		if(run_count>1){//Um.. something goes wrong if this is true
		//			mf::LogDebug("SinglePhoton") << "  It seems that we loop the VertexBuilder too much.\n";
		//		}

		std::vector< art::Ptr<recob::PFParticle> > candidate_particles(package.particles);//CHECK, always vertex the most nu-like slice.

		if(m_run_all_pfps){
			cout<<"d(O.O)b Looking at all slices. ";
		}else{
			cout<<"d(O.O)b Not looking at all slices."<<endl;
		}
		//we dont need the run_all_particle boolean, because it is considered under package.particles;

		if(candidate_particles.size() == 0){
			mf::LogDebug("SinglePhoton") << "  No PFParticles can be considered for vertexing.\n";
		}


		for(const art::Ptr<recob::PFParticle> &pParticle : candidate_particles){ //candidate_particles are determined in above if-else statement.
			const std::vector< art::Ptr<recob::Track> > ToBeAddedTracks((package.PFParticleAsATrack)->at(pParticle.key()));
			const std::vector< art::Ptr<recob::Shower> > ToBeAddedShowers((package.PFParticleAsAShower)->at(pParticle.key()));

			const unsigned int nTracks(ToBeAddedTracks.size());
			const unsigned int nShowers(ToBeAddedShowers.size());

			//Check if we can add a track/a shower that is identified from a PFParticle.
			//		cout<<"CHECK! # of tracks/showers are ready, over!"<<endl;
			//		cout<<nTracks<<endl;
			//		cout<<nShowers<<endl;
			if( nTracks + nShowers == 0 ){//Um... the PFParticle is not identified as track or shower;
				mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
			}else if( nTracks + nShowers > 1 ){
				//Check, do I need throw??
				throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
				//Ok, if going through the below if-elses, it means we have a shower/ a track!

			}else if( nTracks == 1 ){ //Add a Track
				if(package.PFPToNuSliceMap[pParticle]){//add the package.selected nu_slice particle;
					(package.selected_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(package.selected_tracks).back()]= pParticle;

				}else{//add other nu-like PFParticle;
					(package.more_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(package.more_tracks).back()]= pParticle;
				}

				if(m_is_verbose)std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  ToBeAddedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

			}else if( nShowers == 1 ){ //Add a Shower
				if(package.PFPToNuSliceMap[pParticle]) {//add the package.selected nu_slice particle;
					(package.selected_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(package.selected_showers).back()] = pParticle;
				}else{//add other nu-like PFParticle
					(package.more_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(package.more_showers).back()]= pParticle;
				}

				if(m_is_verbose)std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  ToBeAddedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
			} //repeat this loop for another PFParticle in the candidate_particles.
		}
	}
//		package.selected_showers.clear();//no need these, but they dont hurt;
//		package.selected_tracks.clear();
//
//		switch (run_count)
//		{
//			case 0://actually only case 0 allowed here;
//				cout<<"\n\n Load objects form the selected nu slice!"<<endl;
//				package.selected_showers = wanted_showers;
//				package.selected_tracks = wanted_tracks;
//				break;
//			case 1:
//				cout<<"\n\n Load objects form other nu slices!"<<endl;
//				(package.selected_showers).reserve( more_showers.size() + wanted_tracks.size());
//				(package.selected_showers).insert((package.selected_showers).end(), more_showers.begin(), more_showers.end());
//				(package.selected_showers).insert((package.selected_showers).end(), wanted_showers.begin(), wanted_showers.end());
//				//ok, same for tracks
//				(package.selected_tracks).reserve(more_tracks.size()+wanted_tracks.size());
//				(package.selected_tracks).insert((package.selected_tracks).end(), more_tracks.begin(), more_tracks.end());
//				(package.selected_tracks).insert((package.selected_tracks).end(), wanted_tracks.begin(), wanted_tracks.end());
//				break;
//
//			case 3://CHECK, this is for testing purpose.
//				cout<<"\n\n Only more Nuslices! CHECK IAM MESSING THIS UP"<<endl;
//				package.selected_showers = more_showers;
//				package.selected_tracks = more_tracks;
//		}
//				cout<<"Showers: "<<package.selected_showers.size()<<endl;
//				cout<<"Tracks: "<<package.selected_tracks.size()<<endl;
//	}


	/***************************
	 *
	 * ReconsiderMoreCandidates() - add more showers and tracks candidates into current vertex.
	 *
	 * ***************************

	void ReconsiderMoreCandidates(ParticleAssociations_all &candidates, class Atlas &package){

	// Ingredients:
		package.more_showers;
		package.more_tracks;
		vector<double> current_vertex = {candidates.GetRecoVertex().at(0), candidates.GetRecoVertex().at(1), candidates.GetRecoVertex().at(2)};//get the vertex position;

		//GetRecoVertex()
		//fvertex, geoalgo::Point_t in class ParticleAssociation
		//
		// fassociations a vector of ParticleAssociation class, run under AddAssociation()
		//

	//I think I need to do this:
		pas.GetDetectorObjects().AddShowers(package.more_showers);//load tracks
		pas.GetDetectorObjects().AddTracks(package.more_tracks);//load showers
		AssociateTracks(pas);//this is the code for associating the tracks, see 1027 at VertexBuilder.h
		AssociateShowers(pas);//this is the code for associating the showers
	}*/


	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	ParticleAssociations_all SinglePhoton::BobbyVertexBuilder_ext(class Atlas &package, bool more_objects){
		bool fverbose = m_is_verbose;

//		initscr();//initialize the COLS and LINES for the screen size
		int screen_width = 86;//ncurses::COLS;//get the width of the screen!
		
		//Initialize criteria for vertexing,  unit: cm.
		double start_prox = 4;//Max. track proximity threshold (t_max)
		double shower_prox = 10;//Max. shower proximity threshold (s_max)
		double cpoa_vert_prox = 10;//Max. distance btw shower start & cloest approach (dp_max), when no vertex, look at shower bkw projection and track distance
		double cpoa_trackend_prox = 5;//Max. distance btw midway point of impact parameter to a potential vertex (a_max)

		if(fverbose){//Over view of the inputs
			std::cout << "\n\nRun vertex builder with: ";
			cout<<(package.selected_showers).size()<<" shower candidates and "<<endl;
			cout<<(package.selected_tracks).size()<<" track candidates."<<endl;

			cout<<endl;//print out 
			for(int i = 0 ; i < screen_width/2 ;i++)
				cout<<"/*";
			cout<<"\nFinish loading tracks and showers! Start Revertexing."<<endl;
			for(int i = 0 ; i < screen_width/2 ;i++)
				cout<<"/*";
			cout<<"\n"<<endl;
		}
		
		//Declear two classes to kick off vertexing
		VertexBuilder vbuilder;//it was named vb
		ParticleAssociations_all candidates;// it was named pas

		//Initialize creteria for reconstruction.
		vbuilder.SetVerbose(fverbose);
		vbuilder.SetParameters({start_prox, shower_prox, cpoa_vert_prox, cpoa_trackend_prox});

		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);


		candidates.SetVerbose(fverbose);
		candidates.GetDetectorObjects().AddShowers(package.selected_showers);//load tracks
		candidates.GetDetectorObjects().AddTracks(package.selected_tracks);//load showers
		if(more_objects){//CHECK, now just load it without checking bobby result
		candidates.GetDetectorObjects().AddShowers(package.more_showers);//load tracks
		candidates.GetDetectorObjects().AddTracks(package.more_tracks);//load showers
		}



		vbuilder.Run(candidates);//here deals with the candidates and find the vertex.


		cout<<endl;//printout 
		for(int i = 0 ; i < screen_width/2 ;i++)
			cout<<"/*";
		cout<<"\nBobby Revertexing is finished. Now start to fill in the TTree."<<endl;
		for(int i = 0 ; i < screen_width/2 ;i++)
			cout<<"/*";
		cout<<endl;
		cout<<" Reports from the ParticleAssociations_all"<<endl;	

		if(fverbose)candidates.PrintAssociations_all();

		cout<<endl;//printout 
		for(int i = 0 ; i < screen_width/2 ;i++)
			cout<<"/*";
		cout<<"\n"<<endl;
		if(fverbose)
		candidates.PrintNodes();

		cout<<endl;//printout 
		for(int i = 0 ; i < screen_width/2 ;i++)
			cout<<"/*";
		cout<<"\n"<<endl;

		if(vbuilder.f_dist_tt[0]<999){
		m_dist_tt = vbuilder.f_dist_tt;
		}
		if(vbuilder.f_dist_sx[0]<999){
		m_dist_sx = vbuilder.f_dist_sx;
		}

		if(vbuilder.f_dist_st[0]<999){
		m_dist_st = vbuilder.f_dist_st;
		}

		if(vbuilder.f_dist_sst[0]<999){
		m_dist_sst = vbuilder.f_dist_sst;
		}



			return candidates;//and fill in the vertexed file (tree) in the SinglePhoton_module.cc

	}
}
#endif
