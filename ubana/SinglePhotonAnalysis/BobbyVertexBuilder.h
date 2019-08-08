

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
	struct SinglePhoton::ObjectCandidates{//Initialize this with event address;
//		art::Event event;
		std::vector< art::Ptr<recob::Track> > tracks;
		std::vector< art::Ptr<recob::Shower> > showers;
		std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>> trackToNuPFParticleMap;
		std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;
		
	};


	/*****************************
	 * CollectTracksAndShowers () - this associates tracks and showers to one event.
	 *		Tracks and showers come from pfParticles.
	 *	particles(input) - a std::vector< art::Ptr<recob::PFParticle> > for all pfparticle address of the nuslice
	 *	pfParticleMap(input) - a std::map< size_t, art::Ptr<recob::PFParticle>> for all pfparticle address (I think the address is stored in the *.second?).
	 *	pfParticleHandle (input) - ???
	 * **************************/

	void SinglePhoton::CollectTracksAndShowers_v2(const std::vector< art::Ptr<recob::PFParticle> > &particles,
			const std::map< size_t, art::Ptr<recob::PFParticle>> pfParticleMap,
			const art::Event &input_event,
//			const PFParticleHandle &pfParticleHandle,
			struct ObjectCandidates TracksAndShowers){

	}


	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	void SinglePhoton::BobbyVertexBuilder_ext(std::vector<art::Ptr<recob::Track>> & tracks,  std::vector<art::Ptr<recob::Shower>> & showers ){
		bool fverbose = true;

		//PUT THIS FUNCTION INSIDE SINGLEPHOTON_MODULE.h
		VertexBuilder vbuilder;// definition. it was named vb
		ParticleAssociations candidates;// definition. it was named pas



		vbuilder.SetVerbose(fverbose);
		//to get info. of the following variables, see VertexBuilder.h
		//for associating tracks
		if(fverbose){	
			cout<<"Transfering creteria from SinglePhoton class to VertexBuilder class:"<<endl;
			cout<<right<<setw(82)<<" Max. track proximity threshold (t_max)="<<fstart_prox<<endl;
			cout<<right<<setw(82)<<" Max. shower proximity threshold (s_max)="<<fshower_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw shower start & cloest approach (dp_max)="<<fcpoa_vert_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw midway point of impact parameter to a potential vertex (a_max)="<<fcpoa_trackend_prox<<endl;
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

		if(fverbose) std::cout << "\n\nRun vertex builder\n";
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
		//Vertex position are named: (reco_bobbyvertex_x, reco_bobbyvertex_y, reco_bobbyvertex_z);
		//Number of showers/tracks are named: reco_bobbyshowers and reco_bobbytracks.

		//m_bvertex_pos_x = 0;
		//m_bvertex_pos_y = 0;
		//m_bvertex_pos_z = 0;

		cout<<"\n/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*"<<endl;
		cout<<"Bobby Revertexing is finished."<<endl;
		cout<<"/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*\n"<<endl;


	}
}
#endif
