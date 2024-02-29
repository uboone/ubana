#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_JUSTPANDORA_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_JUSTPANDORA_H

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "variables.h"

namespace BVA_ana
{
//Utility return daughter pfp, given the daughte id.
	art::Ptr<recob::PFParticle> findpfpDaughter(const std::vector<art::Ptr<recob::PFParticle>>& pfParticleVector, size_t daughter_id);


//Main content
	void GrabPandoraVertex(lar_pandora::VertexVector vertexVector, var_all& vars);

	void AnalyzeShower(art::Ptr<recob::PFParticle> pfpShowerVector , art::Ptr<recob::Shower> shower); 

//	//T = art::Ptr<recob::Track> or art::Ptr<recob::Shower>
//	template<typename T>
//	T GrabPandoraObjects( std::vector< T >& associatedTracks, var_all& vars){
//		const unsigned int nTracks(associatedTracks.size());
//	if (nTracks == 1 )
//	{
//		std::cout<<"CHECK 1 track "<<std::endl;
//		return associatedTracks.front();
//	}
//	std::cout<<"CHECK # of tracks: "<<nTracks<<std::endl;
//	return {};
//
//	}


}

#endif
