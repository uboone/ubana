

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


	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	void SinglePhoton::BobbyVertexBuilder_ext(class Atlas &package, bool more_objects){
		bool fverbose = m_is_verbose;

//		initscr();//initialize the COLS and LINES for the screen size
		int screen_width = 86;//ncurses::COLS;//get the width of the screen!
		
		//Initialize criteria for vertexing,  unit: cm.
		double start_prox = 8;//Max. track proximity threshold (t_max)
		double shower_prox = 15;//Max. shower proximity threshold (s_max)
		double cpoa_vert_prox = 13;//Max. distance btw shower start & cloest approach (dp_max), when no vertex, look at shower bkw projection and track distance
		double cpoa_trackend_prox = 18;//Max. distance btw midway point of impact parameter to a potential vertex (a_max)

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


		candidates.SetVerbose(fverbose);
		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);


//regroup tracks and showers info.
		std::vector< art::Ptr<recob::Track> > use_tracks;
		std::vector< art::Ptr<recob::Shower> > use_showers;

		if(more_objects){//CHECK, now just load it without checking bobby result
			//MERGE vector first; dont Addshowers/Tracks twice! otherwise will mess up the indices
			use_tracks.reserve( package.selected_tracks.size() + package.more_tracks.size() );
			use_tracks.insert( use_tracks.end(), package.selected_tracks.begin(), package.selected_tracks.end() );
			use_tracks.insert( use_tracks.end(), package.more_tracks.begin(), package.more_tracks.end() );

			use_showers.reserve( package.selected_showers.size() + package.more_showers.size() );
			use_showers.insert( use_showers.end(), package.selected_showers.begin(), package.selected_showers.end() );
			use_showers.insert( use_showers.end(), package.more_showers.begin(), package.more_showers.end() );
		} else{
			use_tracks = package.selected_tracks;
			use_showers = package.selected_showers;

		}

		candidates.GetDetectorObjects().AddShowers(use_showers);//load tracks
		candidates.GetDetectorObjects().AddTracks( use_tracks);//load showers
//finish loading tracks and showers.

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


		
std::cout<<"Filling in Bobby's Vertex info. with "<<candidates.GetSelectedAssociations().size()<<" Vertex candidates."<<std::endl;
		bool reset_bobbyvertex = true;
		if(candidates.GetSelectedAssociations().size()==0){
		cout<<"No vertex is reconstructed."<<endl;
		}

		for(size_t const nth_associations : candidates.GetSelectedAssociations()) {//Loop over all associations, which is a vector
			ParticleAssociation const & particle_associated = candidates.GetAssociations().at(nth_associations);//grab the "pn"th association;
			geoalgo::Point_t const & reco_vertex = particle_associated.GetRecoVertex();//Grab the vertec of the "pn"th association.
			if(reset_bobbyvertex){
				m_bobbyvertex_pos_xv.clear();
				m_bobbyvertex_pos_yv.clear();
				m_bobbyvertex_pos_zv.clear();
				m_bobbytracksv.clear();
				m_bobbyshowersv.clear();

				m_bobbyprotontrack.clear();
				m_bobbyphotonshower.clear();
				m_bobbypi0daughter.clear();
				reset_bobbyvertex = false;
			}

			m_bobbyvertex_pos_xv.push_back( reco_vertex.at(0));
			m_bobbyvertex_pos_yv.push_back( reco_vertex.at(1));
			m_bobbyvertex_pos_zv.push_back( reco_vertex.at(2));

			cout<<"Vertex Coordinates found by Bobby Vertex Builder: ";
			cout<<reco_vertex.at(0)<<", ";
			cout<<reco_vertex.at(1)<<", ";
			cout<<reco_vertex.at(2)<<endl;

			//calculate the # of tracks/showers;
			DetectorObjects_all const & detos = candidates.GetDetectorObjects();
			int temp_num_tracks = 0;
			int temp_num_showers = 0;
			
			int get_a_proton = 0;
			int get_a_photon = 0;
			int get_a_pi0daughter = 0;
			for(size_t const n : particle_associated.GetObjectIndices()) {
				if(detos.GetRecoType(n) == detos.ftrack_reco_type) {
					++temp_num_tracks;
					
					int trackindex = detos.GetTrackIndexFromObjectIndex(n);
					art::Ptr<simb::MCParticle> temp_MCtrack = package.trackToMCParticleMap.find(use_tracks[trackindex])->second;
					if(temp_MCtrack->PdgCode()==2212){
						get_a_proton++;
					}
				}
				if(detos.GetRecoType(n) == detos.fshower_reco_type) {

					++temp_num_showers;
					int showerindex = detos.GetShowerIndexFromObjectIndex(n);//CHECK
					art::Ptr<simb::MCParticle> temp_MCshower = package.showerToMCParticleMap.find(use_showers[showerindex])->second;
					if(temp_MCshower->PdgCode()==22){
						get_a_photon++;
					}
					art::Ptr<simb::MCParticle> amother = package.MCParticleToTrackIdMap[temp_MCshower->Mother()];
					if(amother){//sometime Mother is unknown..
						if(amother->PdgCode() == 111){
							get_a_pi0daughter++;
						}
					}
				}
			}
			cout<<"# of showers: "<<temp_num_showers<<endl;
			cout<<"# of tracks : "<<temp_num_tracks<<endl;

			m_bobbytracksv.push_back(temp_num_tracks);
			m_bobbyshowersv.push_back(temp_num_showers);

			m_bobbyprotontrack.push_back(get_a_proton);
			m_bobbyphotonshower.push_back(get_a_photon);
			m_bobbypi0daughter.push_back(get_a_pi0daughter);
		}

//			double best_vertex_dist = SIZE_MAX;
			for( size_t index = 0; index< m_bobbyvertex_pos_xv.size() ; index++){
//				double temp_dist =  pow(m_vertex_pos_x - m_bobbyvertex_pos_xv[index],2) + pow(m_vertex_pos_y - m_bobbyvertex_pos_yv[index],2)+ pow(m_vertex_pos_z - m_bobbyvertex_pos_zv[index],2);

//				if(temp_dist < best_vertex_dist){}//update when find a closer vertex to the pandora vertex.
				if((m_bobbyshowersv[index] > 0 && m_bobbytracksv[index] > 0)|| index == 0){
				//	best_vertex_dist = temp_dist;
					
					m_bobbyvertex_pos_x = m_bobbyvertex_pos_xv[index];
					m_bobbyvertex_pos_y = m_bobbyvertex_pos_yv[index];
					m_bobbyvertex_pos_z = m_bobbyvertex_pos_zv[index];
					m_bobbytracks =  m_bobbytracksv[index];
					m_bobbyshowers = m_bobbyshowersv[index];
					break;
				}
			}

//			if( m_run_all_pfps && m_bobbyvertexing_more &&m_bobbytracks + m_bobbyshowers < 2&& run_count < 2 ){//repeat with 0,1
//				run_count++;
//				cout<<"Need more objects. Now run the "<<run_count<<" time the vertexing."<<endl;
//				goto redo_event;
//			}

			std::cout<<"Got Bobby's info.!\n"<<std::endl;

//			return candidates;//and fill in the vertexed file (tree) in the SinglePhoton_module.cc

	}
}
#endif
