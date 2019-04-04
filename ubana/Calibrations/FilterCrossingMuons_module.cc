#ifndef FILTER_FILTERCROSSINGMUONS_H
#define FILTER_FILTERCROSSINGMUONS_H

/// Frameworks includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

/// Frameworks includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

using namespace std;

namespace crossingmufilter {  
  
  class FilterCrossingMuons : public art::EDFilter
  {  
    // explicit EDFilter(ParameterSet const&)  
  public:
    
    explicit FilterCrossingMuons(fhicl::ParameterSet const &pset);
    virtual ~FilterCrossingMuons();
    
    bool filter(art::Event&) ;
    virtual void reconfigure(fhicl::ParameterSet const&)  ;
    
    virtual void beginJob()  ;
    
  private:
    
    std::string fTrackModuleLabel;
    double fXmax;
    double fXmin;
  };
  
} // namespace crossingmufilter


namespace crossingmufilter{
  
  //-----------------------------------------------------------------------
  // Constructor
  FilterCrossingMuons::FilterCrossingMuons(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }
  
  //-----------------------------------------------------------------------
  // Destructor
  FilterCrossingMuons::~FilterCrossingMuons()
  {
  }
  
  //-----------------------------------------------------------------------
  void FilterCrossingMuons::beginJob()
  {
    //art::ServiceHandle<art::TFileService> tfs;
    //art::ServiceHandle<geo::Geometry> geo;
  }
  
  //-----------------------------------------------------------------------
  void FilterCrossingMuons::reconfigure(fhicl::ParameterSet const& p)
  {
    fXmax            = p.get< double >("Xmax");
    fXmin            = p.get< double >("Xmin");
    fTrackModuleLabel = p.get< std::string >("TrackModuleLabel","");

    return; 
  }
  
  //========================================================================
  bool FilterCrossingMuons::filter(art::Event& evt){
    
    bool tracksWanted(false);
    
    //get the list of tracks from this event 
    art::Handle< std::vector<recob::Track> > trackListHandle; 
    std::vector<art::Ptr<recob::Track> > tracklist;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;
      auto& pos = track.Vertex();
      auto& end = track.End();
      float X = std::abs(pos.X()-end.X());
      if(X>fXmin && X<fXmax){
        std::cout << "I found a crossing track!" << std::endl;
	tracksWanted = true;
      }// crossing trks...
    } // loop over trks..     
    
    return tracksWanted;
    
  } // end FilterCrossingMuons() function
  
}// namespace crossingmufilter


namespace crossingmufilter {
  
  DEFINE_ART_MODULE(FilterCrossingMuons)
  
} // namespace crossingmufilter

#endif // FILTER_FILTERCROSSINGMUONS_H
