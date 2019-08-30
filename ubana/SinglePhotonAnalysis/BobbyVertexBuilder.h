

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


using namespace std;

namespace single_photon
{

	/*****************************
	 * CollectTracksAndShowers () - this associates tracks and showers to one event.
	 *		Tracks and showers come from pfParticles.
	 *	particles(input) - a std::vector< art::Ptr<recob::PFParticle> > for all pfparticle address of the nuslice
	 *	pfParticleMap(input) - a std::map< size_t, art::Ptr<recob::PFParticle>> for all pfparticle address (I think the address is stored in the *.second?).
	 *	pfParticleHandle (input) - ???
	 *
	 *
	 *	UPGRADE! Everything but the event
	 *		are packed to the Atlas in v2;
	 * **************************/

	void SinglePhoton::CollectTracksAndShowers_v2(
			const art::Event &evt,
			class Atlas &package,
			int run_count){
		if(run_count>2){//Um.. something goes wrong if this is true
			mf::LogDebug("SinglePhoton") << "  It seems that we loop the VertexBuilder too much.\n";
		}

//		cout<<"CHECK! Get in CollectTracksAndShowers_v2 over!\n\n\n"<<endl;
		//the following function labels objects with sliceID and tell
		//	if the slice is nu_slice or not.
		this->AnalyzeSlices(
				package.PFParticleToMetadataMap, 
				package.PFParticleMap,  
				package.primaryPFPSliceIdVec, 
				package.sliceIdToNuScoreMap, 
				package.PFPToClearCosmicMap, 
				package.PFPToSliceIdMap, 
				package.PFPToNuSliceMap, 
				package.PFPToTrackScoreMap);
		
//		cout<<"CHECK! AnalyzeSlices pass over!\n\n\n"<<endl;
		//set up what particles to look at
		std::vector< art::Ptr<recob::PFParticle> > candidate_particles;

		if(m_run_all_pfps){
			cout<<"d(O.O)b Looking at all slices. ";
		}else{
			cout<<"d(O.O)b Looking at only the most neutrino-like slice."<<endl;
		}
		//we dont need the run_all_particle boolean, because it is considered under package.particles;
		candidate_particles = package.particles;//CHECK, always vertex the most nu-like slice.

		if(candidate_particles.size() == 0){
			mf::LogDebug("SinglePhoton") << "  No PFParticles can be considered for vertexing.\n";
		}

		std::vector< art::Ptr<recob::Track> >		cosmic_tracks;//yea,, I think all PFParticles should include them;
		std::vector< art::Ptr<recob::Shower> >		cosmic_showers;
		std::vector< art::Ptr<recob::Track> >		wanted_tracks;//yea,, I think all PFParticles should include them;
		std::vector< art::Ptr<recob::Shower> >		wanted_showers;

		for(const art::Ptr<recob::PFParticle> &pParticle : candidate_particles){ //candidate_particles are determined in above if-else statement.

//		cout<<"CHECK! Inside the pfParticles, over!"<<endl;
			
			//pfParticleHandle is from recob::PFParticle on m_pandoraLabel;
			//pfParticleToTrackAssoc is a FindManyP<reco::Track> from <pfParticleHandle, evt,m_trackLabel)
			//associatedTracks is a copy of pfPartToTrackAssoc.at(pParticle.key()));
			//nTracks = associatedTracks.size()

			//Make a copy of vectors we need, lazy as I am~
			const std::vector< art::Ptr<recob::Track> > ToBeAddedTracks((package.PFParticleAsATrack)->at(pParticle.key()));
			const std::vector< art::Ptr<recob::Shower> > ToBeAddedShowers((package.PFParticleAsAShower)->at(pParticle.key()));

			const unsigned int nTracks(ToBeAddedTracks.size());
			const unsigned int nShowers(ToBeAddedShowers.size());

			//Check if we can add a track/a shower that is identified from a PFParticle.
//		cout<<"CHECK! # of tracks/showers are ready, over!"<<endl;
			if( nTracks + nShowers == 0 ){//Um... the PFParticle is not identified as track or shower;
//		cout<<"00"<<endl;
				mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
			}else if( nTracks + nShowers > 1 ){
//		cout<<"11"<<endl;
				//Check, do I need throw??
				throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
				
			//Ok, if going through the below if-elses, it means we have a shower/ a track!
			}else if( nTracks == 1 ){ //Add a Track

//				cout<<"10"<<endl;
				if( package.PFPToClearCosmicMap[pParticle] == 1 ){//add cosmic PFParticle;
					(cosmic_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(cosmic_tracks).back()]= pParticle;

				}else if(m_bobbyvertexing_more && !(package.PFPToNuSliceMap[pParticle])){//we want more particles that is not in the most nu-like slice;
					(package.backup_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(package.backup_tracks).back()]= pParticle;

				}else{//not cosmic; in most nu-like slice; record it;
					(wanted_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(wanted_tracks).back()]= pParticle;
				}
				std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  ToBeAddedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

			}else if( nShowers == 1 ){ //Add a Shower
//				cout<<"01"<<endl;
				if( package.PFPToClearCosmicMap[pParticle] == 1 ){//add cosmic PFParticle;
					(cosmic_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(cosmic_showers).back()]= pParticle;

				}else if( m_bobbyvertexing_more && !(package.PFPToNuSliceMap[pParticle]) ){//we want more particles that is not in the most nu-like slice;
					(package.backup_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(package.backup_showers).back()]= pParticle;

				}else{//not cosmic; in most nu-like slice; record it;
					(wanted_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(wanted_showers).back()] = pParticle;
				}
				std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  ToBeAddedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
			}

			//repeat this loop for another PFParticle in the candidate_particles.
		}

		package.collected_showers.clear();//dont know if these help;
		package.collected_tracks.clear();
		switch (run_count)
		{
			case 1:
				cout<<"\n\n Load objects form other nu slices!"<<endl;
				(package.collected_showers).reserve(package.backup_showers.size()+wanted_tracks.size());
				(package.collected_showers).insert((package.collected_showers).end(), package.backup_showers.begin(), package.backup_showers.end());
				(package.collected_showers).insert((package.collected_showers).end(), wanted_showers.begin(), wanted_showers.end());
				//ok, same for tracks
				(package.collected_tracks).reserve(package.backup_tracks.size()+wanted_tracks.size());
				(package.collected_tracks).insert((package.collected_tracks).end(), package.backup_tracks.begin(), package.backup_tracks.end());
				(package.collected_tracks).insert((package.collected_tracks).end(), wanted_tracks.begin(), wanted_tracks.end());

				break;

			case 2:
				cout<<"\n\n Load objects form other nu slices and cosmic slices!"<<endl;

				(package.collected_showers).reserve(cosmic_showers.size()+package.backup_showers.size()+wanted_tracks.size());
				(package.collected_showers).insert((package.collected_showers).end(), cosmic_showers.begin(), cosmic_showers.end());
				(package.collected_showers).insert((package.collected_showers).end(), package.backup_showers.begin(), package.backup_showers.end());
				(package.collected_showers).insert((package.collected_showers).end(), wanted_showers.begin(), wanted_showers.end());
				//ok, same for tracks
				(package.collected_tracks).reserve(cosmic_tracks.size()+package.backup_tracks.size()+wanted_tracks.size());
				(package.collected_tracks).insert((package.collected_tracks).end(), cosmic_tracks.begin(), cosmic_tracks.end());
				(package.collected_tracks).insert((package.collected_tracks).end(), package.backup_tracks.begin(), package.backup_tracks.end());
				(package.collected_tracks).insert((package.collected_tracks).end(), wanted_tracks.begin(), wanted_tracks.end());
				break;

			default://actually only case 0 allowed here;

				cout<<"\n\n Load objects form the selected nu slice!"<<endl;
				package.collected_showers = wanted_showers;
				package.collected_tracks = wanted_tracks;
		}

/*
		if(m_bobbyvertexing_more){//Include all PFParticles
			cout<<"Oh! You want to run on Nu particles, but with more objects from other possible slices."<<endl;
			cout<<"Note: here only record PFParticles and do nothing different to looking at the most neutrino-like slice."<<endl;
			package.collected_showers = wanted_showers;
			package.collected_tracks = wanted_tracks;
		} else{
			cout<<" Ok, actually you want all PFParticles; Here you go!"<<endl;
			(package.collected_showers).reserve(cosmic_showers.size()+package.backup_showers.size()+wanted_tracks.size());
			(package.collected_showers).insert((package.collected_showers).end(), cosmic_showers.begin(), cosmic_showers.end());
			(package.collected_showers).insert((package.collected_showers).end(), package.backup_showers.begin(), package.backup_showers.end());
			(package.collected_showers).insert((package.collected_showers).end(), wanted_showers.begin(), wanted_showers.end());
			//ok, same for tracks
			(package.collected_tracks).reserve(cosmic_tracks.size()+package.backup_tracks.size()+wanted_tracks.size());
			(package.collected_tracks).insert((package.collected_tracks).end(), cosmic_tracks.begin(), cosmic_tracks.end());
			(package.collected_tracks).insert((package.collected_tracks).end(), package.backup_tracks.begin(), package.backup_tracks.end());
			(package.collected_tracks).insert((package.collected_tracks).end(), wanted_tracks.begin(), wanted_tracks.end());
		}
*/

	}


	/***************************
	 *
	 * ReconsiderMoreCandidates() - add more showers and tracks candidates into current vertex.
	 *
	 * ***************************

	void ReconsiderMoreCandidates(ParticleAssociations &candidates, class Atlas &package){

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

	ParticleAssociations SinglePhoton::BobbyVertexBuilder_ext(class Atlas &package){
		bool fverbose = true;

		//PUT THIS FUNCTION INSIDE SINGLEPHOTON_MODULE.h
		VertexBuilder vbuilder;// definition. it was named vb
		ParticleAssociations candidates;// definition. it was named pas


		vbuilder.SetVerbose(fverbose);
		//to get info. of the following variables, see VertexBuilder.h
		//for associating tracks
		if(fverbose){	
			cout<<"Transfering creteria from SinglePhoton class to VertexBuilder class:"<<endl;
			cout<<right<<setw(82)<<" Max. track proximity threshold (t_max)= "<<fstart_prox<<endl;
			cout<<right<<setw(82)<<" Max. shower proximity threshold (s_max)= "<<fshower_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw shower start & cloest approach (dp_max)= "<<fcpoa_vert_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw midway point of impact parameter to a potential vertex (a_max)= "<<fcpoa_trackend_prox<<endl;
		}
		vbuilder.SetMaximumTrackEndProximity(fstart_prox);//Set the maximum track proximity threshold (in cm)		
		//for associating showers (this include connection to tracks)
		vbuilder.SetMaximumShowerIP(fshower_prox);
		vbuilder.CPOAToVert(fcpoa_vert_prox);
		vbuilder.SetMaximumTrackEndProx(fcpoa_trackend_prox);

		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);

		candidates.SetVerbose(fverbose);

		if(fverbose) std::cout << "\n\nRun vertex builder with: \n";
		//Object inside an object cause the problem, i.e. 
		//ParticleAssociations cannot link to DetectorOBjects;
		//CHekced: Members are identified; but the reference is not defined
		cout<<"Number of Showers: "<<(package.collected_showers).size()<<endl;
		cout<<"Number of Tracks: "<<(package.collected_tracks).size()<<endl;
		candidates.GetDetectorObjects().AddShowers(package.collected_showers);//load tracks
		candidates.GetDetectorObjects().AddTracks(package.collected_tracks);//load showers

		cout<<"\n/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*"<<endl;
		cout<<"Finish loading tracks and showers! Start Revertexing."<<endl;
		cout<<"/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*\n"<<endl;


		vbuilder.Run(candidates);//here deals with the candidates and find the vertex.


		cout<<"\n/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*"<<endl;
		cout<<"Bobby Revertexing is finished. Now start to fill in the TTree."<<endl;
		cout<<"/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*\n"<<endl;
		
		return candidates;//and fill in the vertexed file (tree) in the SinglePhoton_module.cc

	}
}
#endif
