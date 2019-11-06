

#ifndef __BOBBYVERTEXBUILDER_H__
#define __BOBBYVERTEXBUILDER_H__

#include "SinglePhoton_module.h"
#include "Atlas.h"

//---VertexBuilder---
#include "VertexBuilder/VertexBuilder.h"
#include "VertexBuilder/ParticleAssociations.h" 
//#include "larsim/MCCheater/ParticleInventoryService.h"

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

//		if(m_run_all_pfps){
//			cout<<"d(O.O)b Looking at all slices. ";
//		}else{
//			cout<<"d(O.O)b Not looking at all slices."<<endl;
//		}
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


    void SinglePhoton::CollectMCParticles_v2(
		const art::Event &evt, 
		class Atlas & package){
		bool debug_message = false;

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(m_geantModuleLabel, theParticles);

        if (!theParticles.isValid()){

            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        } else {

            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, m_geantModuleLabel);
		package.MCParticleToTrackIdMap.clear();

		int num_of_mothers = 0;
		for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i){//theParticles are pointers to all MCParticles.
			//theTruthAssns are MCTruth with geantModuleLabel
			const art::Ptr<simb::MCParticle> particle(theParticles, i);
			const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));

			package.MCTruthToMCParticlesMap[truth].push_back(particle);
			package.MCParticleToMCTruthMap[particle] = truth;//truth ->GetParticles() will give a bunch of particles that has not been linked.
			package.MCParticleToTrackIdMap[particle->TrackId()] = particle;

			//---------- Find out MCParticles that will lead to the ancestors, i.e. particles produced by Genie
			art::Ptr<simb::MCParticle> temp_mc = particle;
			//			int temindex=0;
			//			if(debug_message){
			//				std::cout<<"\nNew MCParticle "<<i+1<<"th, with size"<<package.matchedMCParticleVector.size();
			//				std::cout<<"Looking at Track ID:"<<temp_mc->TrackId()<<std::endl; //gives the id of the MCParticle
			//			}
			//			//Id is not continuous OMG..
			while(true){//search for ancestor; initial state - 0; other state - 1; https://internal.dunescience.org/doxygen/namespacegenie.html#a05cd2ccc34b3e3a9e88bdd335f990118
				int mothersId = temp_mc->Mother();
				//				if(debug_message && mothersId == 0 && temp_mc == particle){
				//				}
				if(mothersId == 0){
					package.AncestorToPdgMap[particle] = temp_mc->PdgCode();//default ancestor pdg is -999, which means nothing;
					if(temp_mc==particle){//this is the first generation of particle from geant
						num_of_mothers++;
						if(debug_message){ 
							std::cout<<"SinglePhoton::CollectMCParticles_v2() \t||\t Track ID:"<<setw(5)<<temp_mc->TrackId(); //gives the id of the MCParticle
							std::cout<<" Mother ID:"<<setw(5)<<mothersId; //gives the id of the MCParticle
							std::cout<<" Pdg Code:"<<setw(5)<<temp_mc->PdgCode(); //gives the id of the MCParticle
							cout<<" Identify an ancestor particle."<<endl;
						}
					}
					break;
				}

				if(package.MCParticleToTrackIdMap.find(mothersId)==package.MCParticleToTrackIdMap.end()){//not found,

					if(debug_message)cout<<"MCParticle's mother is not found"<<endl;
					break;
				} else {//go up one generation to find mother's mother.

					temp_mc = package.MCParticleToTrackIdMap.find(mothersId)->second;
				}
			}

				package.MCParticleToAncestorMap.emplace(particle,temp_mc);

		}
				//				if(mothersId == 0 && temp_mc == particle){
				//					num_of_mothers++;
				//					if(debug_message) cout<<"Count mother"<<endl;
				//					break;
				//				}else{
				//				//	std::cout<<" Update MC to ID:"<<mothersId<<std::endl; //gives the id of the MCParticle
				//					if(package.MCParticleToTrackIdMap.find(mothersId)==package.MCParticleToTrackIdMap.end()){//not found,
				//						if(debug_message&& temp_mc == particle)cout<<"MCParticle's mother is not found"<<endl;
				//						break;
				//					}else{
				//						temp_mc = package.MCParticleToTrackIdMap.find(mothersId)->second;
				//					}
				//				}
				//----------- num_of_mothers gives # of MCParticles that has mother() as 0. -------------


		//determine whether mothers are delta radiative products by
		//	comparing energy sum (photon & proton/neutron).
		std::vector< art::Ptr< simb::MCParticle>> temp_MCPmother;
		std::vector< art::Ptr< simb::MCTruth>> temp_MCTruth;
		std::map< pair<size_t, size_t>, double >  ParticlesEnergySum;//pair up two MCParticles and sum up their energy

		for( int i = 0; i < num_of_mothers; i++ ){//collect MCParticle candidates
			art::Ptr<simb::MCParticle> check_mcp0 = package.matchedMCParticleVector[i];
			switch(check_mcp0->PdgCode()){
				case 22://gamma
				case 2112://neutron
				case 2212://proton
					temp_MCPmother.push_back(check_mcp0);//new vector for mother candidates!
					//From this line on, DONT USE package.matchedMCParticleVector!
					break;
				default://not adding others
				{}
			}
		}
		//if temp_MCPmother.size()==0, that means no radiative;

		for(size_t i = 0; i < temp_MCPmother.size(); i++){//fillin ParticlesEnergySum. from the selected simb::MCParticle, temp_MCPmother
			for(size_t j = 1; j < temp_MCPmother.size() - i; j++){
			double total_energy = temp_MCPmother[i]->E() + temp_MCPmother[i+j]->E();
			pair<size_t,size_t> mcparticles(i,j);
			if(debug_message){
				cout<<"Pair up "<<temp_MCPmother[i]->TrackId()<<" with "<<temp_MCPmother[i]->PdgCode();
				cout<<" and "<<temp_MCPmother[j]->TrackId()<<" with "<<temp_MCPmother[j]->PdgCode()<<endl;
			}

			ParticlesEnergySum[mcparticles] = total_energy;
			}
		}

		double delta_energy;
		int temp_pdgcode = -999;//record delta pdg code;
		int expected_pdg = 0;
		int check_counter = 0;
		pair<size_t, size_t> want_this_pair; 
		if(temp_MCPmother.size()>0){//no need to check energy if no deltaradiative;
			for(size_t i = 0; i<package.mcTruthVector.size(); i++){
				for(int j = 0; j<package.mcTruthVector[i]->NParticles(); j++){
					//collect MCTruth candidate;
					simb::MCParticle check_mct = package.mcTruthVector[i]->GetParticle(j);
					switch(check_mct.PdgCode()){//find MCTruth that is delta, hopefully only one delta;
						case 2224://delta++ (not likely)
						case 1114://delta- (not likely)
						case 2214://delta+
							expected_pdg = 2212 + 22;//proton & gamma;

						case 2114://delta0
							{
								if(expected_pdg ==0){
									expected_pdg = 2112 + 22;//proton & gamma;
								}

								delta_energy = check_mct.E();
								temp_pdgcode = check_mct.PdgCode();
								double resolution = 0.0005;
								while(resolution > 0.00005){	
									check_counter = 0;
									for(auto const & [this_pair, this_energy] : ParticlesEnergySum){
										art::Ptr<simb::MCParticle> mcp1 = temp_MCPmother[this_pair.first];
										art::Ptr<simb::MCParticle> mcp2 = temp_MCPmother[this_pair.second];
										if(debug_message){
											cout<<"delta Energy "<<delta_energy<<" and combined energy "<<this_energy;
											cout<<" Code "<< mcp1->PdgCode()<< " and Code "<<mcp2->PdgCode()<<endl;
										}
										if( abs(this_energy - delta_energy)< resolution ){
											if(debug_message)	cout<<"Energy matched!"<<endl;
											int check_pdgcode = mcp1->PdgCode()+mcp2->PdgCode();
											if(check_pdgcode == expected_pdg){
												if(debug_message) cout<<"PdgCode matched!"<<endl;
												want_this_pair = this_pair;
												check_counter++;
												break;
											}else{
												if(debug_message) cout<<"PdgCode not matched!"<<check_pdgcode<<" with "<<expected_pdg<<endl;
											}
										}
									}
									resolution = resolution/10;
									if(check_counter==1) break;
								}
								break;
							}
						default:{}//not adding others
					}
				}
			}
			package.AncestorToPdgMap[temp_MCPmother[want_this_pair.first]] = temp_pdgcode;//map the delta pdg to the AncestorToPdgMap;
			package.AncestorToPdgMap[temp_MCPmother[want_this_pair.second]] = temp_pdgcode;//map the delta pdg to the AncestorToPdgMap;
		}

		if(check_counter>1){
			cout<<"More than 1 delta or 1 gamma in genie? Check ";
			cout<<__FILE__<<" at "<<__LINE__<<endl;
			exit(1);
		}
			//----------------

        std::cout<<"SinglePhoton::CollectMCParticles_v2() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
		if(debug_message)cout<<"Take a look at the AncestorMap"<<endl;
		for(auto const & this_iterator : package.AncestorToPdgMap){//map the ancestor pdg to each MCParticle
			auto this_mcp = this_iterator.first;
			auto mother_mcp = package.MCParticleToAncestorMap.find(this_mcp)->second;
			int ancestor_pdg = package.AncestorToPdgMap.find(mother_mcp)->second;

			package.AncestorToPdgMap[this_mcp] = ancestor_pdg;
			if(debug_message){
			cout<<"Ancestor pdg is "<<ancestor_pdg;
			cout<<" Track Id: "<< this_mcp->TrackId()<<" has pdg: "<< this_mcp->PdgCode();
			cout<<" with ancestor: "<< mother_mcp->TrackId();
			cout<<" whose Pdg is: "<<this_iterator.second<<endl;
			}
		}
		if(debug_message)cout<<"Finish taking a look at the AncestorMap"<<endl;

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
//		for (auto it:package.trackToDistMap){
//			if(it.second > temp_max_value) temp_max_value = it.second;
//		}
//		for (auto it:package.showerToDistMap){
//			if(it.second > temp_max_value) temp_max_value = it.second;
//		}//ok, now find the furthest distance;


		//Declear two classes to kick off vertexing
		VertexBuilder vbuilder;//it was named vb
		ParticleAssociations_all candidates;// it was named pas

		//Initialize creteria for reconstruction.
		vbuilder.SetParameters({start_prox, shower_prox, cpoa_vert_prox, cpoa_trackend_prox});

		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);
		//prepare parameters for pre-check # of tracks and showers;
		std::vector< art::Ptr<recob::Track> > use_tracks(package.selected_tracks);//CHECK all tracks
		std::vector< art::Ptr<recob::Shower> > use_showers(package.selected_showers);
		auto trackmap = package.trackToDistMap;
		auto showermap = package.showerToDistMap;
//		cout<<"CHECK SIZE "<<package.trackToDistMap.size()<<endl;
		//-- Get a preview of # of tracks and showers first	
		std::vector<int>	tem_bobbyshowersv;
		std::vector<int> 	tem_bobbytracksv;
		double nearby_vertex_dist = 0;
		double temp_max_value = 200;//objects that is too far..
		bool successful_vertex = true;//true, do them all together; false, do them by step;
		bool one_for_all = true;//true, means load more_tracks/showers all at once; sucessful_vertex has to be true to make it work;
		std::vector< int > target_topo = {1,1};//{1,2}<->1s2t signal;
		bool reset_bobbyvertex = true;
		int loop_tracker = 1;
		int max_loops = 4;//consider all objects within temp_max_value distance after the 4th loop.
	
		//sort out values of both maps first
		vector<double> all_dist;
		for (auto it:trackmap){
			if(it.second > temp_max_value){
				trackmap.erase(it.first);
			} else{
				all_dist.push_back(it.second);
			}
		}
		for (auto it:showermap){
			if(it.second > temp_max_value){
				showermap.erase(it.first);
			} else{
				all_dist.push_back(it.second);
			}
		}//ok, now find the furthest distance;
		std::sort(all_dist.begin(), all_dist.end());// distance from low to high;
//		for(auto x:all_dist){
//			cout<<x<<" ";
//		}

		if(fverbose){
		cout<<"\nPreview Bobby Reco. Preformance"<<endl;
		cout<<setw(9)<<right<<"Attempts ";
		cout<<setw(14)<<right<<"|Input shower ";
		cout<<setw(14)<<left<< "|Input track";
		cout<<setw(14)<<right<<"|Reco. showers ";
		cout<<setw(14)<<left<< "|Reco. of tracks "<<endl;
		}
		while(!successful_vertex&&all_dist.size()>0){
			ParticleAssociations_all candidates_copy(candidates);
			candidates_copy.SetVerbose(false);
			vbuilder.SetVerbose(false);
			//CHECK control what tracks and showers to look at
			if(more_objects){
				nearby_vertex_dist = all_dist[(all_dist.size()-1)*loop_tracker/max_loops] + 0.01;
				cout<<"within the distance: "<<nearby_vertex_dist<<endl;
//				cout<<"Search for additional tracks"<<endl;
				for (auto it:trackmap){
//					cout<<it.second<<endl;
					//if(it.second >9998){//remove objects in nu slice
					//	trackmap.erase(it.first);
					//	continue;
					//}
					if(it.second < nearby_vertex_dist){
						use_tracks.push_back(it.first);
						trackmap.erase(it.first);//eliminate for next round of consideration
					}
				}
//				cout<<"Search for additional showers"<<endl;
				for (auto it:showermap){
//				cout<<it.second<<endl;
					//if(it.second >9998){
					//	showermap.erase(it.first);
					//	continue;
					//}
					if(it.second < nearby_vertex_dist){
						use_showers.push_back(it.first);
						showermap.erase(it.first);
					}
				}
			} 

			candidates_copy.GetDetectorObjects().AddShowers(use_showers);//load tracks
			candidates_copy.GetDetectorObjects().AddTracks( use_tracks);//load showers
			vbuilder.Run(candidates_copy);//here deals with the candidates_copy and find the vertex.

			//Coming update:
			//1. Label slice with Pandora reco. position
			//2. Add 4 slices at m.
			//3. Improve algo. later;


			//Fill in Bobby varaibles;
//			std::cout<<"Preview Topological Reco. result status: "<<candidates_copy.GetSelectedAssociations().size()<<" vertex candidates_copy."<<std::endl;


			if(candidates_copy.GetSelectedAssociations().size() == 0){
				cout<<"No vertex is reconstructed."<<endl;
			}


			for(size_t const nth_associations : candidates_copy.GetSelectedAssociations()) {//Loop over all associations, which is a vector
				ParticleAssociation const & particle_associated = candidates_copy.GetAssociations().at(nth_associations);//grab the "pn"th association;

				//calculate the # of tracks/showers;
				DetectorObjects_all const & detos = candidates_copy.GetDetectorObjects();
				int temp_num_tracks = 0;
				int temp_num_showers = 0;

				for(size_t const n : particle_associated.GetObjectIndices()) {
					if(detos.GetRecoType(n) == detos.ftrack_reco_type) {
						++temp_num_tracks;
					}
					if(detos.GetRecoType(n) == detos.fshower_reco_type) {
						++temp_num_showers;
					}
				}
				if(( temp_num_showers >= target_topo[0]  && temp_num_tracks >= target_topo[1])||all_dist.size()<1) successful_vertex = true;

				if(fverbose){//Over view of the inputs
					cout<<setw(9)<<left<<loop_tracker<<"|";
					cout<<setw(13)<<use_showers.size()<<"|";
					cout<<setw(13)<<use_tracks.size()<<"|";
					cout<<setw(13)<<left<<temp_num_showers<<"|";
					cout<<setw(13)<<left<<temp_num_tracks<<endl;
				}
			}
			if(fverbose){
				for(int i = 0; i<screen_width; i++) 
					cout<<"-";
				cout<<endl;
			}

			loop_tracker++;
			if(loop_tracker > max_loops){
				break;
			}
		}//end loop of testing diff. inputs for vertexing.
//		cout<<"CHECKCHECK"<<endl;
		candidates.SetVerbose(fverbose);
		vbuilder.SetVerbose(fverbose);
		//Looks good, then proceed to really fill in trees;

		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\nBobby Revertexing is finished. Now start to fill in the TTree."<<endl;
		////			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<endl;
		//			cout<<" Reports from the ParticleAssociations_all"<<endl;	

		//			if(fverbose)candidates.PrintAssociations_all();
		//
		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\n"<<endl;
		//			if(fverbose)
		//				candidates.PrintNodes();
		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\n"<<endl;

		candidates.GetDetectorObjects().AddShowers(use_showers);//load showers
		if(more_objects && one_for_all) candidates.GetDetectorObjects().AddShowers(package.more_showers);//load more showers
		candidates.GetDetectorObjects().AddTracks( use_tracks);//load tracks
		cout<<"ADD TRACKS:"<<use_tracks.size()<<endl;
		if(more_objects && one_for_all) candidates.GetDetectorObjects().AddTracks( package.more_tracks);//load more tracks
		cout<<"ADD TRACKS:"<<package.more_tracks.size()<<endl;
		vbuilder.Run(candidates);//here deals with the candidates_copy and find the vertex.

		m_dist_tt = vbuilder.f_dist_tt;

		m_dist_sx = vbuilder.f_dist_sx;

		m_dist_st = vbuilder.f_dist_st;

		m_dist_sst = vbuilder.f_dist_sst;


		//temporary varaibles
		std::vector<double> tem_bobbyvertex_pos_xv;
		std::vector<double> tem_bobbyvertex_pos_yv;
		std::vector<double> tem_bobbyvertex_pos_zv;
		std::vector<int> 	tem_bobbyphotonshowerv;
		std::vector<int> 	tem_bobbypi0daughterv;
		std::vector<int> 	tem_bobbydeltaraddaughterv;
		std::vector<int> 	tem_bobbyotherdaughterv;
		std::vector<int> 	tem_bobbyprotontrackv;

		std::vector<std::vector <int>> tem_track_daughter_pdgv;
		std::vector<std::vector <int>> tem_shower_daughter_pdgv;
		double min_bobbyvertexradius = 999;
		size_t min_index = 0;
		size_t temp_counter = 0;

		for(size_t const nth_associations : candidates.GetSelectedAssociations()) {//Loop over all associations, which is a vector
			ParticleAssociation const & particle_associated = candidates.GetAssociations().at(nth_associations);//grab the "pn"th association;
			geoalgo::Point_t const & reco_vertex = particle_associated.GetRecoVertex();//Grab the vertec of the "pn"th association.
			if(reset_bobbyvertex){
				m_bobbyvertex_pos_xv.clear();
				m_bobbyvertex_pos_yv.clear();
				m_bobbyvertex_pos_zv.clear();
				m_bobbytracksv.clear();
				m_bobbyshowersv.clear();

				m_bobbyprotontrackv.clear();
				m_bobbyphotonshowerv.clear();
				m_bobbypi0daughterv.clear();
				m_bobbydeltaraddaughterv.clear();
				m_bobbyotherdaughterv.clear();
				m_bobbyvertexradiusv.clear();
				m_bobbytrackdaughter_pdg.clear();
				m_bobbyshowerdaughter_pdg.clear();
				reset_bobbyvertex = false;
			}

			tem_bobbyvertex_pos_xv.push_back( reco_vertex.at(0));
			tem_bobbyvertex_pos_yv.push_back( reco_vertex.at(1));
			tem_bobbyvertex_pos_zv.push_back( reco_vertex.at(2));

			cout<<"BB Vertex: (";
			cout<<reco_vertex.at(0)<<", ";
			cout<<reco_vertex.at(1)<<", ";
			cout<<reco_vertex.at(2)<<") ";

			//calculate the # of tracks/showers;
			DetectorObjects_all const & detos = candidates.GetDetectorObjects();
			int temp_num_tracks = 0;
			int temp_num_showers = 0;

			int get_a_proton = 0;
			int get_a_photon = 0;
			int get_a_pi0daughter = 0;
			int get_a_deltaraddaughter = 0;
			int get_a_otherdaughter = 0;
			std::vector <int> get_a_track_daughter_pdg;
			std::vector <int> get_a_shower_daughter_pdg;

			for(auto const & [a,b]:package.trackToMCParticleMap){
				cout<<a->StartMomentum()<<" "<<b->TrackId()<<endl;
			}
			cout<<"CHECK GAPPP- ---"<<endl;
			for(auto const & [a,b]:package.showerToMCParticleMap){
				cout<<a->ID()<<" "<<b->TrackId()<<endl;
			}

			for(size_t const n : particle_associated.GetObjectIndices()) {
				int index;
				art::Ptr<simb::MCParticle> temp_mcp;

				if(detos.GetRecoType(n) == detos.ftrack_reco_type) {//it is a track
					cout<<"A track"<<endl;
					++temp_num_tracks;
					index = detos.GetTrackIndexFromObjectIndex(n);

				cout<<"CHECK 0 "<<__LINE__<<"Do "<<use_tracks[index]->StartMomentum()<<endl;
					temp_mcp = package.trackToMCParticleMap.find(use_tracks[index])->second;//CHECK THIS
				cout<<"CHECK 0 "<<__LINE__<<endl;
					get_a_track_daughter_pdg.push_back(temp_mcp->PdgCode());
				}

				if(detos.GetRecoType(n) == detos.fshower_reco_type) {//it is a shower
					cout<<"A shower"<<endl;
					++temp_num_showers;
					index = detos.GetShowerIndexFromObjectIndex(n);
					temp_mcp = package.showerToMCParticleMap.find(use_showers[index])->second;
					get_a_shower_daughter_pdg.push_back(temp_mcp->PdgCode());
				}
				//identify shower/track MCTruth info.
				if(m_is_verbose)cout<<"CHECK PARTICLE PdgCode "<<temp_mcp->PdgCode()<<endl;
				switch (temp_mcp->PdgCode()){
					case 2212:
						get_a_proton++;
						break;
					case 22:
						get_a_photon++;
						break;
					default:
					{}
				}
				//identify the ancestor!
				switch( package.AncestorToPdgMap.find(temp_mcp)->second){
					case 2214:
					case 2114:
						cout<<"Parent is delta"<<endl;
						get_a_deltaraddaughter++;
						break;
					case 111:
						cout<<"Parent is pi0"<<endl;
						get_a_pi0daughter++;
						break;
					default:
						get_a_otherdaughter++;
				}

				m_bobbyvertexradius = particle_associated.GetGoodness();
				//mark down the smallest radius;
				if(m_bobbyvertexradius < min_bobbyvertexradius){
					min_bobbyvertexradius = m_bobbyvertexradius;
					min_index = temp_counter;
				}
				temp_counter++;
			}
			cout<<"Reconstruct a vertex with "<<endl;
			cout<<temp_num_showers<<" showers, ";
			cout<<temp_num_tracks<<" tracks."<<endl;

			tem_bobbytracksv.push_back(temp_num_tracks);
			tem_bobbyshowersv.push_back(temp_num_showers);

			tem_bobbyprotontrackv.push_back(get_a_proton);
			tem_bobbyphotonshowerv.push_back(get_a_photon);
			tem_bobbypi0daughterv.push_back(get_a_pi0daughter);
			tem_bobbydeltaraddaughterv.push_back(get_a_deltaraddaughter);
			tem_bobbyotherdaughterv.push_back(get_a_otherdaughter);
//	CHECK, need to figure out how to push vector into a vector<vector>
//			tem_track_daughter_pdgv.push_back(get_a_track_daughter_pdg);
			tem_shower_daughter_pdgv.push_back(get_a_shower_daughter_pdg);

			m_bobbyvertexradiusv.push_back(m_bobbyvertexradius);
		}


		//Finalize Bobby info. and pick a vertex set.

		m_bobbyvertex_pos_xv = tem_bobbyvertex_pos_xv;
		m_bobbyvertex_pos_yv = tem_bobbyvertex_pos_yv;
		m_bobbyvertex_pos_zv = tem_bobbyvertex_pos_zv;
		m_bobbyprotontrackv  = tem_bobbyprotontrackv;
		m_bobbyphotonshowerv = tem_bobbyphotonshowerv;
		m_bobbypi0daughterv  = tem_bobbypi0daughterv;
		m_bobbydeltaraddaughterv  = tem_bobbydeltaraddaughterv;
		m_bobbyotherdaughterv  = tem_bobbyotherdaughterv;
		m_bobbytracksv = tem_bobbytracksv;
        m_bobbyshowersv= tem_bobbyshowersv;

//		for( size_t index = 0; index< m_bobbyvertex_pos_xv.size() ; index++){//loop over bobbyvertex
			//metric for vertex:
			//Topology:
			//VertexSize/doogness:

			//				if(temp_dist < best_vertex_dist){}//update when find a closer vertex to the pandora vertex.
			if(true){//maybe pick the smallest vertex radius (for now),
				//	best_vertex_dist = temp_dist;
				m_bobbyvertex_pos_x = m_bobbyvertex_pos_xv[min_index];
				m_bobbyvertex_pos_y = m_bobbyvertex_pos_yv[min_index];
				m_bobbyvertex_pos_z = m_bobbyvertex_pos_zv[min_index];
				m_bobbytracks =  m_bobbytracksv[min_index];
				m_bobbyshowers = m_bobbyshowersv[min_index];
				m_bobbyphotonshower = m_bobbyphotonshowerv[min_index];
				m_bobbypi0daughter = m_bobbypi0daughterv[min_index];
				m_bobbydeltaraddaughter = m_bobbydeltaraddaughterv[min_index];
				m_bobbyotherdaughter = m_bobbyotherdaughterv[min_index];
				m_bobbyprotontrack = m_bobbyprotontrackv[min_index];
//				m_bobbytrackdaughter_pdg = tem_track_daughter_pdgv[min_index];
//				m_bobbyshowerdaughter_pdg = tem_shower_daughter_pdgv[min_index];
				//CHECK add shower/track daughter here;
			}
//		}

		//			if( m_run_all_pfps && m_bobbyvertexing_more &&m_bobbytracks + m_bobbyshowers < 2&& run_count < 2 ){//repeat with 0,1
		//				run_count++;
		//				cout<<"Need more objects. Now run the "<<run_count<<" time the vertexing."<<endl;
		//				goto redo_event;
		//			}

		std::cout<<"Got Bobby's info.!\n"<<std::endl;

	}
}
#endif
