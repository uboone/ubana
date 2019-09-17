

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
	 * CollectTracksAndShowers_v2() - Fill in package.collected_showers and
	 *	package.collected_tracks based on the run_count.
	 *
	 *  packahe (to be modified) - contains all vectors and maps that are
	 *		needed to access pandora_objects; it will be modified as
	 *		the code runs.
	 *
	 *	run_count (input) - has values 0,1,2; each number indicates differnt collections of slices: 
	 *		0 - selected nu_slice; 
	 *		1 - not-selected nu_slices; 
	 *		2 - cosmic slices (non-nu_slices).  NOT DOING THIS
	 *
	 * **************************/

	void SinglePhoton::CollectTracksAndShowers_v2(
			const art::Event &evt,
			class Atlas &package,
			int run_count){
		if(run_count>1){//Um.. something goes wrong if this is true
			mf::LogDebug("SinglePhoton") << "  It seems that we loop the VertexBuilder too much.\n";
		}

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

		std::vector< art::Ptr<recob::Track> >	wanted_tracks;
		std::vector< art::Ptr<recob::Shower> >	wanted_showers;
		std::vector< art::Ptr<recob::Track> >	backup_tracks;
		std::vector< art::Ptr<recob::Shower> >	backup_showers;

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
				if(!(package.PFPToNuSliceMap[pParticle])){//add other nu-like PFParticle;
					(backup_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(backup_tracks).back()]= pParticle;

				}else{//add selected nu-slice particle.
					(wanted_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(wanted_tracks).back()]= pParticle;
				}

				if(m_is_verbose)std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  ToBeAddedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

			}else if( nShowers == 1 ){ //Add a Shower
				if(!(package.PFPToNuSliceMap[pParticle]) ){//add other nu-like PFParticle
					(backup_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(backup_showers).back()]= pParticle;

				}else{//add selected nu-slice particle
					(wanted_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(wanted_showers).back()] = pParticle;
				}

				if(m_is_verbose)std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  ToBeAddedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
			} //repeat this loop for another PFParticle in the candidate_particles.
		}

		package.collected_showers.clear();//no need these, but they dont hurt;
		package.collected_tracks.clear();

		switch (run_count)
		{
			case 0://actually only case 0 allowed here;
				cout<<"\n\n Load objects form the selected nu slice!"<<endl;
				package.collected_showers = wanted_showers;
				package.collected_tracks = wanted_tracks;
				break;
			case 1:
				cout<<"\n\n Load objects form other nu slices!"<<endl;
				(package.collected_showers).reserve( backup_showers.size() + wanted_tracks.size());
				(package.collected_showers).insert((package.collected_showers).end(), backup_showers.begin(), backup_showers.end());
				(package.collected_showers).insert((package.collected_showers).end(), wanted_showers.begin(), wanted_showers.end());
				//ok, same for tracks
				(package.collected_tracks).reserve(backup_tracks.size()+wanted_tracks.size());
				(package.collected_tracks).insert((package.collected_tracks).end(), backup_tracks.begin(), backup_tracks.end());
				(package.collected_tracks).insert((package.collected_tracks).end(), wanted_tracks.begin(), wanted_tracks.end());
				break;

			case 3://CHECK, this is for testing purpose.
				cout<<"\n\n Only backup Nuslices! CHECK IAM MESSING THIS UP"<<endl;
				package.collected_showers = backup_showers;
				package.collected_tracks = backup_tracks;
		}
				cout<<"Showers: "<<package.collected_showers.size()<<endl;
				cout<<"Tracks: "<<package.collected_tracks.size()<<endl;
	}


	/***************************
	 *
	 * ReconsiderMoreCandidates() - add more showers and tracks candidates into current vertex.
	 *
	 * ***************************

	void ReconsiderMoreCandidates(ParticleAssociations_all &candidates, class Atlas &package){

	// Ingredients:
		package.backup_showers;
		package.backup_tracks;
		vector<double> current_vertex = {candidates.GetRecoVertex().at(0), candidates.GetRecoVertex().at(1), candidates.GetRecoVertex().at(2)};//get the vertex position;

		//GetRecoVertex()
		//fvertex, geoalgo::Point_t in class ParticleAssociation
		//
		// fassociations a vector of ParticleAssociation class, run under AddAssociation()
		//

	//I think I need to do this:
		pas.GetDetectorObjects().AddShowers(package.backup_showers);//load tracks
		pas.GetDetectorObjects().AddTracks(package.backup_tracks);//load showers
		AssociateTracks(pas);//this is the code for associating the tracks, see 1027 at VertexBuilder.h
		AssociateShowers(pas);//this is the code for associating the showers
	}*/


	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	ParticleAssociations_all SinglePhoton::BobbyVertexBuilder_ext(class Atlas &package){
		bool fverbose = m_is_verbose;

//		initscr();//initialize the COLS and LINES for the screen size
		int screen_width = 86;//ncurses::COLS;//get the width of the screen!
		
		//Initialize criteria for vertexing,  unit: cm.
		double start_prox = 4;//Max. track proximity threshold (t_max)
		double shower_prox = 10;//Max. shower proximity threshold (s_max)
		double cpoa_vert_prox = 10;//Max. distance btw shower start & cloest approach (dp_max)
		double cpoa_trackend_prox = 5;//Max. distance btw midway point of impact parameter to a potential vertex (a_max)

		if(fverbose){//Over view of the inputs
			std::cout << "\n\nRun vertex builder with: \n";
			cout<<"Number of shower candidates: "<<(package.collected_showers).size()<<endl;
			cout<<"Number of track candidates : "<<(package.collected_tracks).size()<<endl;

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
		candidates.GetDetectorObjects().AddShowers(package.collected_showers);//load tracks
		candidates.GetDetectorObjects().AddTracks(package.collected_tracks);//load showers


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

		return candidates;//and fill in the vertexed file (tree) in the SinglePhoton_module.cc

	}
}
#endif