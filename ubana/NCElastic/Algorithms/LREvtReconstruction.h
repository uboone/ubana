/**
 * \file LREvtReconstruction.h
 *
 * @author Katherine Woodruff
 */

#ifndef LREVTRECONSTRUCTION_H
#define LREVTRECONSTRUCTION_H

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include <iostream>

namespace qeselection{
  
  /**
   \class LREvtReconstruction
   User defined class LREvtReconstruction ... these comments are used to generate
   doxygen documentation!
 */

  class LREvtReconstruction {
    
  public:
    
    /// Default constructor
    LREvtReconstruction();
    
    /// Default destructor
    ~LREvtReconstruction(){}

    void Configure(fhicl::ParameterSet const& p);

    /// Determine if event is NCEp, CCQEn, or neither
    std::vector<float> ReconstructEvent( art::Event const&);

   protected:

    unsigned int _nendhits;

    float _beamwinstartt;
    float _beamwinendt;    
    float _minflashpe;

    float _minplength;
    float _minpscore;

    float _fid_xmin;
    float _fid_xmax;
    float _fid_ymin;
    float _fid_ymax;
    float _fid_zmin;
    float _fid_zmax;

    float _minmuscore;
    float _minpiscore;

    std::string _trackmodulelabel;
    std::string _caloassoclabel;
    std::string _flashmodulelabel;
    std::string _dtassoclabel;

    std::vector<float> _lrcoefficients;

    /// Find minimum z and y flash distance per track
    std::vector<float> MinFlashDistance(art::Handle< std::vector<recob::OpFlash> >,
                                                             art::Ptr<recob::Track>);

    const anab::CosmicTagID_t TAGID_P  = anab::CosmicTagID_t::kGeometry_YY;
    const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
    const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
    const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
    const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;

  };
}

#endif
/** @} */ // end of doxygen group 

