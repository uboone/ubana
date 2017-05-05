/**
 *  @file   ChargedTrackMultiplicityAlg.cxx
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices.
 *          It is based of original NuMuCCInclusiveAlg, geared towards ChargeTrack Multiplicty
 *          analysis.
 *
 *          Date of last update: 5 May 2017
 *          Authors: Algorithm/Analysis -- Aleena Rafique (aleena@ksu.edu)
 *                   Module/LArSoft technicals -- Tracy Usher (usher@slac.standford.edu) and Wesley Ketchum (wketchum@fnal.gov)
 * 
 */

// The main include
#include "uboone/TPCNeutrinoIDFilter/Algorithms/ChargedTrackMultiplicityAlg.h"

// Framework Includes
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

ChargedTrackMultiplicityAlg::ChargedTrackMultiplicityAlg(fhicl::ParameterSet const &pset) :
    fMyProducerModule(0),
    fGeometry(lar::providerFrom<geo::Geometry>()),
    fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ChargedTrackMultiplicityAlg::~ChargedTrackMultiplicityAlg()
{
}
    
void ChargedTrackMultiplicityAlg::reconfigure(fhicl::ParameterSet const &inputPset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& pset = inputPset.get<fhicl::ParameterSet>("ChargedTrackMultiplicityAlg");
    
    fTrackModuleLabel         = pset.get<std::string> ("TrackModuleLabel");
    fVertexModuleLabel        = pset.get<std::string> ("VertexModuleLabel");
    fOpFlashModuleLabel       = pset.get<std::string> ("OpFlashModuleLabel");
    
    fDistToEdgeX              = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",   6.);
    fDistToEdgeY              = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   10.);
    fDistToEdgeZ              = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   6.);
    
    fFlashWidth               = pset.get<double>      ("FlashWidth",                            100.);
    fBeamMin                  = pset.get<double>      ("BeamMin",                              3.2);
    fBeamMax                  = pset.get<double>      ("BeamMax",                              4.8);
    fPEThresh                 = pset.get<double>      ("PEThresh",                              50.);
    fMinTrk2VtxDist           = pset.get<double>      ("MinTrk2VtxDist",                         5.);
    fMinTrackLen              = pset.get<double>      ("MinTrackLen",                           60.);
    
    fDoHists                  = pset.get<bool>        ("FillHistograms",                      false);
    fCreateAnalysisCollection = pset.get<bool>        ("CreateAnalysisCollection",            false);
    if(fCreateAnalysisCollection)
      fAnaMaxTrackDistance = pset.get<double>         ("AnaMaxTrackDistance", 5.);
}
    
void ChargedTrackMultiplicityAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fNFlashPerEvent   = tfs->make<TH1D>("NFlashEvent", ";Flash/Event",     200,   0.,  200.);
        fFlashPE          = tfs->make<TH1D>("FlashPE",     ";PE",              100,   0.,  100.);
        fFlashTime        = tfs->make<TH1D>("FlashTime",   ";Flash Time(us)",  100, -10.,   30.);
    }
    
    return;
}
    
void ChargedTrackMultiplicityAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;

    fMyProducerModule->produces< art::Assns<recob::Track, recob::Vertex> >();
    //fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();

    if(fCreateAnalysisCollection){
      fMyProducerModule->produces< std::vector<recob::Vertex> >();
      fMyProducerModule->produces< std::vector<recob::Track> >();
      fMyProducerModule->produces< std::vector<recob::Hit> >();
      fMyProducerModule->produces< art::Assns<recob::Track, recob::Hit> >();
    }
}

    
bool ChargedTrackMultiplicityAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
  std::unique_ptr<art::Assns<recob::Track,recob::Vertex>>      vertexTrackAssociations(new art::Assns<recob::Track,recob::Vertex>);
  //std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

    //these will only be included in event if fCreateAnalysisCollection is set
    std::unique_ptr< std::vector<recob::Vertex> >          anaVertexCollection(new std::vector<recob::Vertex> );
    std::unique_ptr< std::vector<recob::Track> >           anaTrackCollection(new std::vector<recob::Track> );
    std::unique_ptr< std::vector<recob::Hit> >             anaHitCollection(new std::vector<recob::Hit> );
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::OpFlash>> flashListHandle;
    
    event.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,     trackVecHandle);
    
    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);

    //final vertex to be used
    int    FinalVertexCandidate=-1;
    int    FinalTrackCandidate=-1;
    
    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        // Recover associations to PFParticles...
        art::FindManyP<recob::PFParticle> trackToPFPartAssns(trackVecHandle,  event, fTrackModuleLabel);
        
        //----loop over all the flashes and check if there are flashes within the beam
        //window and above the PE threshold
        const recob::OpFlash* flashPtr(0);
        double                flashmax(0);
        bool                  flashtag(false);
        
        for(const auto& opFlash : flashlist)
        {
            if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
            {
                flashtag = true;
                
                // Keep track of the largest flash
                if (opFlash->TotalPE() > flashmax)
                {
                    flashPtr = opFlash.get();
                    flashmax = opFlash->TotalPE();
                }
            }
            
            if (fDoHists)
            {
                fFlashPE->Fill(opFlash->TotalPE(), 1.);
                fFlashTime->Fill(opFlash->Time(), 1.);
            }
        }  //end of loop over all the flashes
        
        if (fDoHists) fNFlashPerEvent->Fill(flashlist.size(), 1.);

	
        if(flashtag)
        {
            // We need to keep track of the best combination
            // Can we assign art ptrs? I don't think so...
            int    VertexCandidate=-1;
            int    TrackCandidate=-1;
            double TrackCandLength=0;
            double trackstartzcandidate=0;
            double trackstartxcandidate=0;
            double trackstartycandidate=0;
            double trackendzcandidate=0;
            double trackendxcandidate=0;
            double trackendycandidate=0;
            
            //-----------------------------------------------------------
            for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
            {
                // Recover art ptr to vertex
                art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
                
                // Get the position of the vertex
                // Ultimately we really want the vertex position in a TVector3 object...
                double vertexXYZ[3];
                
                vertex->XYZ(vertexXYZ);
                
                TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
            	   
                if(inFV(vertexPos.X(),vertexPos.Y(),vertexPos.Z()))
                {
                    // For each vertex we loop over all tracks looking for matching pairs
                    // The outer loop here, then is over one less than all tracks
                    for(size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
                    {
                        // Work with an art Ptr here
                        art::Ptr<recob::Track> track(trackVecHandle,trackIdx);
                        
                        // so we need to get the track direction sorted out.
                        TVector3 trackPos = track->Vertex();
                        TVector3 trackEnd = track->End();
                        
                        // Take the closer end---------------------------------
                        double trackToVertexDist = (trackPos - vertexPos).Mag();
                        
                        if ((trackEnd - vertexPos).Mag() < trackToVertexDist)
                        {
                            trackPos          = track->End();
                            trackEnd          = track->Vertex();
                            trackToVertexDist = (trackPos - vertexPos).Mag();
                        }
                        
                        //--------------------------------------------------------------------------
                        if(trackToVertexDist<fMinTrk2VtxDist)
                        {
                            if((trackEnd-trackPos).Mag()>TrackCandLength)
                            {
                                TrackCandLength = (trackEnd-trackPos).Mag();
                                TrackCandidate=trackIdx;
                                VertexCandidate=vertexIdx;
                                trackstartzcandidate=trackPos.z();
                                trackstartxcandidate=trackPos.x();
                                trackstartycandidate=trackPos.y();
                                trackendzcandidate=trackEnd.z();
                                trackendxcandidate=trackEnd.x();
                                trackendycandidate=trackEnd.y();
                            }
                        } //end of if track distance is within 5cm
                    }  //end of loop over the tracks
                }  //end of if the vertex is contained
            } //end of loop over all the vertex
            
            if(TrackCandidate > -1)
            {
                bool trackContainedFlag = FlashTrackDist(flashPtr->ZCenter(), trackstartzcandidate, trackendzcandidate) < fFlashWidth;
                
                // Check to see if we think we have a candidate
                if(TrackCandLength>fMinTrackLen && trackContainedFlag && inFV(trackstartxcandidate, trackstartycandidate, trackstartzcandidate) && inFV(trackendxcandidate, trackendycandidate, trackendzcandidate) )
                {
		  FinalVertexCandidate = VertexCandidate;
		  FinalTrackCandidate = TrackCandidate;
		  
		  /*		  		  
                    // Make an association between the best vertex and the matching tracks
                    art::Ptr<recob::Vertex> vertex(vertexVecHandle,VertexCandidate);
                    art::Ptr<recob::Track>  track(trackVecHandle,TrackCandidate);
                    
                    //util::CreateAssn(*fMyProducerModule, event, track, vertex, *vertexTrackAssociations);
                    util::CreateAssn(*fMyProducerModule, event, vertex, track, *vertexTrackAssociations);
		  */
		  /*  
                    // Find the associated PFParticle
                    std::vector<art::Ptr<recob::PFParticle>> pfParticleVec = trackToPFPartAssns.at(track.key());
                    
                    if (!pfParticleVec.empty())
                    {
                        util::CreateAssn(*fMyProducerModule, event, pfParticleVec[0], vertex, *vertexPFParticleAssociations);
                    }
		  */
                }
            }
        }  //end of if flag
    }

    // Now, if we want to create a separate collection for analysis objects, we do that here.
    if(fCreateAnalysisCollection && FinalVertexCandidate>-1 && FinalTrackCandidate>-1){

      recob::Vertex const& myVertex((*vertexVecHandle).at(FinalVertexCandidate));

      //put the vertex onto the analysis collection
      anaVertexCollection->push_back(myVertex);

      double vtx_xyz[3]; myVertex.XYZ(vtx_xyz);
      bool keep_track=false;
      
      for(size_t i_trk=0; i_trk<trackVecHandle->size(); ++i_trk){

	keep_track=false;
	auto const& trk( (*trackVecHandle)[i_trk] );
	
	for(size_t i_pt=0; i_pt<trk.NumberTrajectoryPoints(); ++i_pt){
	  auto const& pt = trk.TrajectoryPoint(i_pt);
	  double dist = std::sqrt( (pt.position.x()-vtx_xyz[0])*(pt.position.x()-vtx_xyz[0]) +
				   (pt.position.y()-vtx_xyz[1])*(pt.position.y()-vtx_xyz[1]) +
				   (pt.position.z()-vtx_xyz[2])*(pt.position.z()-vtx_xyz[2]) );
	  if(dist < fAnaMaxTrackDistance){
	    keep_track = true;
	    break;
	  }
	}//end loop over trajectory points

	if(keep_track || (int)i_trk==FinalTrackCandidate){

	  //	  std::cout << "Adding track " << i_trk << std::endl;

	  //push back track
	  anaTrackCollection->push_back(trk);
	  
	  if((int)i_trk==FinalTrackCandidate){
	    //std::cout << "Adding final track " << i_trk << std::endl;
	    util::CreateAssn(*fMyProducerModule,event,*anaTrackCollection,*anaVertexCollection,*vertexTrackAssociations,
			     anaVertexCollection->size()-1,anaVertexCollection->size());
	  }
	  
	  //get associated hits to track
	  //note, there is a chance here that the same hit, if associated to multiple tracks, could be added to
	  //the event multiple times, as I'm not checking for that here.
	  art::FindMany<recob::Hit> trackToHitAssns(trackVecHandle,  event, fTrackModuleLabel);
	  auto hitPtrVec = trackToHitAssns.at(i_trk);

	  //push back the hits
	  size_t prev_size = anaHitCollection->size();
	  for(auto hitptr : hitPtrVec)
	    anaHitCollection->push_back(*hitptr);

	  //create track<-->hit associations
	  util::CreateAssn(*fMyProducerModule,event,*anaTrackCollection,*anaHitCollection,*anaTrackHitAssociations,prev_size,anaHitCollection->size());
	  
	}//end keep track condition

      }//end loop over tracks
      
    }

    //std::cout << "There are " << vertexTrackAssociations->size() << " vtx-trk associations." << std::endl;

    // Add associations to event.
    event.put(std::move(vertexTrackAssociations));
    //event.put(std::move(vertexPFParticleAssociations));

    if(fCreateAnalysisCollection){

      event.put(std::move(anaVertexCollection));
      event.put(std::move(anaTrackCollection));
      event.put(std::move(anaHitCollection));
      event.put(std::move(anaTrackHitAssociations));
    }
    
    return true;
}
    
bool ChargedTrackMultiplicityAlg::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * fGeometry->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}

//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double ChargedTrackMultiplicityAlg::FlashTrackDist(double flash, double start, double end) const
{
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}

} // namespace
