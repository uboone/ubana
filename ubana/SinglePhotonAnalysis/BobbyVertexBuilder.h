

#ifndef BOBBYVERTEXBUILDER_H
#define BOBBYVERTEXBUILDER_H

#include "SinglePhoton_module.h"

//---VertexBuilder---
#include "VertexBuilder/VertexBuilder.h"
#include "VertexBuilder/ParticleAssociations.h" 
#include "VertexBuilder/DetectorObjects.h"
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

		vbuilder.Run(candidates);

		m_bvertex_pos_x = 0;
		m_bvertex_pos_y = 0;
		m_bvertex_pos_z = 0;

		cout<<"\n\n Run Run Run Away, Run Away Baby"<<endl;

 }
 }
#endif
