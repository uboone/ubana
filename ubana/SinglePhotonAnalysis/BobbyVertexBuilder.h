

#ifndef __BOBBYVERTEXBUILDER_H__
#define __BOBBYVERTEXBUILDER_H__

#include "SinglePhoton_module.h"

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
		/***************
	 * A class that designed for storing addresses for all associated (to an event) tracks, showers, 
	 *		and their cooresponding PFParticles.
	 *	evt (input) - the event that we currently look at.
	 *	tracks (modified) - a vector contains associated track pointers. 
	 *	showers (modified) - a vector contains associated shower pointers.
	 *	trackToNuPFParticleMap (modified) - the map btw the track and the PFParticle.
	 *	showerToNuPFParticleMap (modified) - the map btw the shower and the PFParticle.
	 * *************/

//	typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
//	typedef std::vector< art::Ptr<recob::Track> > TrackVector;
//	typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
//	typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;
	class ObjectCandidates{//Initialize this with event address;
//		art::Event event;//dont care what event we are looking at.
		friend class SinglePhoton;//private as default

//		protected:
//		art::Event evt;
//		art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;

		public:
		//Constructor
		ObjectCandidates();
//		ObjectCandidates(const art::Event &evt);

//		~ObjectCandidates(){};
//		//main stuffs that we feed into the vertex builder.
		std::vector< art::Ptr<recob::PFParticle> >						particles;
		std::vector< art::Ptr<recob::Track> >							all_tracks;
		std::vector< art::Ptr<recob::Shower> >							all_showers;

		//Collections of tracks/showers that to be fed into the vertex builder, 
		//	it is different from the track/shower map.
		std::vector< art::Ptr<recob::Track> >							collected_tracks;
		std::vector< art::Ptr<recob::Shower> >							collected_showers;

		//Pairs that connect PFParticle to sliceID.
		std::vector<std::pair<art::Ptr<recob::PFParticle>,int>>			primaryPFPSliceIdVec;
		
		//Maps for more pandora objects.
		std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>	trackToNuPFParticleMap;
		std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;
		std::map< size_t, art::Ptr<recob::PFParticle>>					PFParticleMap;
		std::map<int, double>											sliceIdToNuScoreMap;
		std::map<art::Ptr<recob::PFParticle>,bool>						PFPToClearCosmicMap;
		std::map<art::Ptr<recob::PFParticle>, int> 						PFPToSliceIdMap;
		std::map<art::Ptr<recob::PFParticle>,bool> 						PFPToNuSliceMap;
		std::map<art::Ptr<recob::PFParticle>,double>					PFPToTrackScoreMap;
		std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > PFParticleToMetadataMap;

		//FindManyP's!
		//specially for the number of coorresponding recob (pandora_objects) to a PFParticle;
		// exp: PFParticleIsATrack[particles.key()] gives the vector that containts all 
		//		cooresponding tracks;
		art::FindManyP< recob::Track	>*								PFParticleAsATrack;
		art::FindManyP< recob::Shower	>*								PFParticleAsAShower;
//		art::FindManyP< recob::Track	> PFParticleAsATrack(PFParticleHandle, evt, m_trackLabel);
//		art::FindManyP< recob::Shower	> PFParticleAsAShower(PFParticleHandle, evt, m_showerLabel);
	};
		//Constructor
		ObjectCandidates::ObjectCandidates(){}
		//Overloaded Constructor
//		ObjectCandidates::ObjectCandidates(const art::Event &evt){
//			const PFParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);//This is useful for FindManyP< reco::Track/Shower>
//		}

	/*****************************
	 * CollectTracksAndShowers () - this associates tracks and showers to one event.
	 *		Tracks and showers come from pfParticles.
	 *	particles(input) - a std::vector< art::Ptr<recob::PFParticle> > for all pfparticle address of the nuslice
	 *	pfParticleMap(input) - a std::map< size_t, art::Ptr<recob::PFParticle>> for all pfparticle address (I think the address is stored in the *.second?).
	 *	pfParticleHandle (input) - ???
	 *
	 *
	 *	UPGRADE! Everything but the event
	 *		are packed to the ObjectCandidates in v2;
	 * **************************/

	void SinglePhoton::CollectTracksAndShowers_v2(
			const art::Event &evt,
			class ObjectCandidates &package){
		cout<<"CHECK! Get in CollectTracksAndShowers_v2 over!\n\n\n"<<endl;
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
		
		cout<<"CHECK! AnalyzeSlices pass over!\n\n\n"<<endl;
		//set up what particles to look at
		std::vector< art::Ptr<recob::PFParticle> > candidate_particles;

		if(m_run_all_pfps){//Include all PFParticles

			cout<<"d(O.O)b Looking at all slices."<<endl;
			//determine what tracks or showers to be rescued from other slices;
			for(std::pair<const long unsigned int, art::Ptr<recob::PFParticle> > pair : package.PFParticleMap){
				
                const art::Ptr<recob::PFParticle> &pParticle = pair.second;

				if(!m_bobbyvertexing_more || (package.PFPToClearCosmicMap[pParticle] < 1/*This means most nu-like slice&& package.PFPToNuSliceMap[pParticle]*/) ){
				//consider the pfParticle if it is not clear cosmic or we dont want more objects (for vertexing).
				candidate_particles.push_back(pParticle);
				}
			}

		}else{//Include PFParticles in the selected nu_slice
			cout<<"d(O.O)b Looking at only the most neutrino-like slice."<<endl;
			candidate_particles = package.particles;
		}

		if(candidate_particles.size() == 0){
			mf::LogDebug("SinglePhoton") << "  No PFParticles can be considered for vertexing.\n";
		}

		//look into the candidate PFParticles and add them to tracks or showers;
		for(const art::Ptr<recob::PFParticle> &pParticle : candidate_particles){ //candidate_particles are determined in above if-else statement.

		cout<<"CHECK! Inside the pfParticles, over!"<<endl;
			
			//pfParticleHandle is from recob::PFParticle on m_pandoraLabel;
			//pfParticleToTrackAssoc is a FindManyP<reco::Track> from <pfParticleHandle, evt,m_trackLabel)
			//associatedTracks is a copy of pfPartToTrackAssoc.at(pParticle.key()));
			//nTracks = associatedTracks.size()
			
			//Make a copy of vectors we need, lazy as I am~
			const std::vector< art::Ptr<recob::Track> > ToBeAddedTracks((package.PFParticleAsATrack)->at(pParticle.key()));
			const std::vector< art::Ptr<recob::Shower> > ToBeAddedShowers((package.PFParticleAsAShower)->at(pParticle.key()));
				
		//CHECK Idk why this works..
//			art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);//This is useful for FindManyP< reco::Track/Shower>
//			FillTracksAndShowers(ToBeAddedTracks,ToBeAddedShowers, pParticle,  pfParticleHandle, evt, package.collected_tracks, package.collected_showers, package.trackToNuPFParticleMap, package.showerToNuPFParticleMap);
//			FillTracksAndShowers_v2(ToBeAddedTracks,ToBeAddedShowers, pParticle, evt, package);
		cout<<"??? FIllTracksAndShowers works!!"<<endl;

			const unsigned int nTracks(ToBeAddedTracks.size());
			const unsigned int nShowers(ToBeAddedShowers.size());

			//Check if we can add a track/a shower that is identified from a PFParticle.
		cout<<"CHECK! # of tracks/showers are ready, over!"<<endl;
			if(nTracks + nShowers == 0){//Um... the PFParticle is not identified as track or shower;
		cout<<"00"<<endl;
				mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
			} else if(nTracks + nShowers > 1){
		cout<<"11"<<endl;
				//Check, do I need throw??
				throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
				
			//Ok, if going through the below if-elses, it means we have a shower/ a track!
			} else if( nTracks == 1){ //Add a Track
				
		cout<<"10"<<endl;
				(package.collected_tracks).push_back(ToBeAddedTracks.front());
				(package.trackToNuPFParticleMap)[(package.collected_tracks).back()]= pParticle;
				std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  ToBeAddedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

			} else if( nShowers == 1){ //Add a Shower
		cout<<"01"<<endl;

				(package.collected_showers).push_back(ToBeAddedShowers.front());
				(package.showerToNuPFParticleMap)[(package.collected_showers).back()] = pParticle;
				std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  ToBeAddedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
			}

		//repeat this loop another PFParticle in the candidate_particles.
		}
	}



	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	ParticleAssociations SinglePhoton::BobbyVertexBuilder_ext(std::vector<art::Ptr<recob::Track>> & tracks,  std::vector<art::Ptr<recob::Shower>> & showers ){
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
		//
		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);
		//
		candidates.SetVerbose(fverbose);

		if(fverbose) std::cout << "\n\nRun vertex builder with: \n";
		//Object inside an object cause the problem, i.e. 
		//ParticleAssociations cannot link to DetectorOBjects;
		//CHekced: Members are identified; but the reference is not defined
		cout<<"Number of Showers: "<<showers.size()<<endl;
		cout<<"Number of Tracks: "<<tracks.size()<<endl;
		candidates.GetDetectorObjects().AddShowers(showers);//load tracks
		candidates.GetDetectorObjects().AddTracks(tracks);//load showers

		cout<<"\n/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*"<<endl;
		cout<<"Finish loading tracks and showers! Start Revertexing."<<endl;
		cout<<"/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*\n"<<endl;

		vbuilder.Run(candidates);


		cout<<"\n/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*"<<endl;
		cout<<"Bobby Revertexing is finished. Now start to fill in the TTree."<<endl;
		cout<<"/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*\n"<<endl;
		
		return candidates;//and fill in the vertexed file (tree) in the SinglePhoton_module.cc

	}
}
#endif
