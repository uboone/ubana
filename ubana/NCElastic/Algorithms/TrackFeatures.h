/**
 * \file TrackFeatures.h
 *
 * @author Katherine Woodruff
 */

#ifndef TRACKFEATURES_H
#define TRACKFEATURES_H

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include <iostream>

namespace dtfeatures{
  
  /**
   \class TrackFeatures
   User defined class TrackFeatures ... these comments are used to generate
   doxygen documentation!
 */

  class TrackFeatures {
    
  public:
    
    /// Default constructor
    TrackFeatures();
    
    /// Default destructor
    ~TrackFeatures(){}

    void Configure(fhicl::ParameterSet const& p);

    /// Create all of the track features used for decision tree training and predicting
    std::vector<std::vector<float>> CreateFeatures( art::Event const&,
                                   art::Handle<std::vector<recob::Track>>& trackVecHandle );

    /// Create class label of the track
    std::vector<anab::CosmicTagID_t> CosmicTagLabel();

   protected:

    unsigned int _nendhits;
    float _beamwinstartt;
    float _beamwinendt;    

    std::string _pidassoclabel;
    std::string _caloassoclabel;
    std::string _cosmictagassoclabel;
    std::string _containtagassoclabel;

  };
}

#endif
/** @} */ // end of doxygen group 

