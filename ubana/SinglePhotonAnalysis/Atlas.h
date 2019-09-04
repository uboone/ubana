#ifndef __ATLAS_H__
#define __ATLAS_H__

#include "SinglePhoton_module.h"

namespace single_photon
{		/***************
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
	class Atlas{//Initialize this with event address;
		friend class SinglePhoton;//private as default

//		protected:
//		art::Event evt;
//		art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
		public:
		//Constructor
		Atlas();//this might be useless;
		//Overflow constructor 1;
		Atlas( const art::Event &evt,
			std::vector< std::string > labels);//this initialize vectors and maps below

//		~Atlas(){};
//		//main stuffs that we feed into the vertex builder.
		std::vector< art::Ptr<recob::PFParticle> >						particles;//this is loaded depends on the option m_run_all_pfps, configurable in the .fcl file.


		//the following recob objects are created through the overflow constructor 1 (see above);
		std::vector< art::Ptr<recob::Track> >							all_tracks;
		std::vector< art::Ptr<recob::Shower> >							all_showers;
		std::vector< art::Ptr<recob::Hit> >							all_hits;
		std::vector< art::Ptr<recob::OpFlash> >							all_opflashes;
		std::vector< art::Ptr<recob::Cluster> >							all_clusters;


		//Collections of tracks/showers that to be fed into the vertex builder, initially empty and it is filled from the CollectTracksAndShowers_v2() in BobbyVertexBuilder.h
		//	it is different from the track/shower map.
		std::vector< art::Ptr<recob::Track> >							collected_tracks;
		std::vector< art::Ptr<recob::Shower> >							collected_showers;

//		std::vector< art::Ptr<recob::PFParticle> >						backup_PFParticles;
		std::vector< art::Ptr<recob::Track> >							backup_tracks; //These two backup_objects are used to search for additional objects.
		std::vector< art::Ptr<recob::Shower> >							backup_showers; 


		//Maps constructed in the constructor of this class.
		std::map< size_t, art::Ptr<recob::PFParticle>>					IDToPFParticleMap;
		std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > PFParticlesToVerticesMap;//old name pfParticlesToVerticesMap;
		std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > PFParticleToMetadataMap;//old name pfParticleToMetadataMap;
		std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > PFParticleToSpacePointsMap;


		//Maps for more pandora objects.
		std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>	trackToNuPFParticleMap;
		std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;
		std::map<int, double>											sliceIdToNuScoreMap;
		//Pairs that connect PFParticle to sliceID.
		std::vector<std::pair<art::Ptr<recob::PFParticle>,int>>			primaryPFPSliceIdVec;
		std::map<art::Ptr<recob::PFParticle>,bool>						PFPToClearCosmicMap;
		std::map<art::Ptr<recob::PFParticle>, int> 						PFPToSliceIdMap;
		std::map<art::Ptr<recob::PFParticle>,bool> 						PFPToNuSliceMap;
		std::map<art::Ptr<recob::PFParticle>,double>					PFPToTrackScoreMap;

		//FindManyP's!
		//specially for the number of coorresponding recob (pandora_objects) to a PFParticle;
		// example: PFParticleIsATrack[particles.key()] gives the vector that containts all 
		//		cooresponding tracks;
		art::FindManyP< recob::Track	>*								PFParticleAsATrack;
		art::FindManyP< recob::Shower	>*								PFParticleAsAShower;

//		art::FindManyP< recob::Track	> PFParticleAsATrack(PFParticleHandle, evt, m_trackLabel);
//		art::FindManyP< recob::Shower	> PFParticleAsAShower(PFParticleHandle, evt, m_showerLabel);
		
		private:
		/***********************
		 *
		 * Wigets for constructing the different types of varaibels.in this class
		 * 1. HandleToVector;
		 *
		 * **********************/
		//1. Scratch from Handle, and return a equivalently useful vector.
		//sample usage:
		//		art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); //yea, this should be gone;
		//		recob::Hit dummy_hit;//This is to specify the template;
		//		std::vector<art::Ptr<recob::Hit>> hitVector = HandleToVector(dummy_hit, evt, m_hitfinderLabel);
		template <typename recob_object>//A helper template that allows you to make compliated types.
			struct temporary_types{ 
				using type1 = std::vector<art::Ptr<recob_object>>;
				using type2 = art::ValidHandle<std::vector<recob_object>>;
				using type3 = std::vector<recob_object>;
			};
		template <class recob_object>//ref_type is only used to identify the temporary_types late.
			typename temporary_types<recob_object>::type1 HandleToVector(recob_object ref_type, const art::Event &evt, std::string &label){

				typename temporary_types<recob_object>::type2 const & Handle = evt.getValidHandle<typename temporary_types<recob_object>::type3>(label);
				typename temporary_types<recob_object>::type1 Vector;
				art::fill_ptr_vector(Vector,Handle);
				return Vector;
			}

		//2. Wiget 2 comes here..

	};
	

//sth similar to class DetectorObjects

//create struct for Track and Shower to add more detail to them, i.e. the geoalgo class;

//struct Track
//struct Shower








//THE HEADER FILE IS THE ABOVE



	//Constructor
	Atlas::Atlas (){}
	//Overloaded Constructor, initialize the essential variables
	Atlas::Atlas ( const art::Event &evt,
				std::vector<std::string > labels){

		//PREPARE some recob objects;
		//vector<string> labels = {m_trackLabel, m_showerLabel, m_hitfinderLabel, m_flashLabel, m_pandoraLabel}
		recob::Track dummy_track;//This is to specify the template;
		all_tracks = HandleToVector(dummy_track, evt, labels[0]);//m_trackLabel

		recob::Shower dummy_shower;//This is to specify the template;
		all_showers = HandleToVector(dummy_shower, evt, labels[1]);//m_showerLabel
		recob::Hit	dummy_hit;
		all_hits = HandleToVector(dummy_hit, evt, labels[2]);//m_hitfinderLabel

		recob::OpFlash	dummy_opflash;
		all_opflashes = HandleToVector(dummy_opflash, evt, labels[3]);//m_flashLabel

		recob::Cluster	dummy_cluster;
		all_clusters = HandleToVector(dummy_cluster, evt, labels[4]);//m_pandoraLabel



		//CREATE maps!
		//Ingredient 1: all PFParticles; I temporary define it here for mapping purpose;
		art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(labels[4]);//This is useful for FindManyP< reco::Track/Shower>
		recob::PFParticle dummy_PFParticle;
		std::vector<art::Ptr<recob::PFParticle>> all_pfparticles = HandleToVector(dummy_PFParticle, evt, labels[4]);
		
		//Ingredient 2: FindManyPs; these will be gone when construction finished
		art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  labels[4]);
		art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, labels[4]);
		art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, labels[4]);

/* CANT do these here, because they will be gone when the construction is finished.
 *      art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, labels[0]);
        art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, labels[1]);
		PFParticleAsATrack = &pfPartToTrackAssoc;//record this to the atlas
		PFParticleAsAShower = &pfPartToShowerAssoc;//record this to the atlas
*/		
		//make maps here;
		for(size_t i=0; i< all_pfparticles.size(); ++i){//old name pfParticleVector, 0~52 in the test sample
			const art::Ptr<recob::PFParticle> pfp = all_pfparticles[i];
			PFParticlesToVerticesMap[pfp] = vertices_per_pfparticle.at(pfp.key());//old name: pfParticlesToVerticesMap;

			PFParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
			IDToPFParticleMap[pfp->Self()] = pfp;
			PFParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
		}
	}


}




#endif
