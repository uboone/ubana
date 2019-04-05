////////////////////////////////////////////////////////////////////////
// Class:       DecisionTreeID
// Plugin Type: producer (art v2_05_00)
// File:        DecisionTreeID_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TString.h"
#include "TTree.h"

#include "ubana/XGBoost/xgboost/c_api.h"
#include "ubana/NCElastic/Algorithms/TrackFeatures.h"

#include <memory>
#include <cmath>
#include <stdexcept>

class DecisionTreeID;

namespace treeid {
      class DecisionTreeParticleID;
}

class DecisionTreeID : public art::EDProducer {
public:
  explicit DecisionTreeID(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DecisionTreeID(DecisionTreeID const &) = delete;
  DecisionTreeID(DecisionTreeID &&) = delete;
  DecisionTreeID & operator = (DecisionTreeID const &) = delete;
  DecisionTreeID & operator = (DecisionTreeID &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // Pointer to the feature algorithm
  ::dtfeatures::TrackFeatures fTrackFeatures;

  std::string ftrackmodulelabel;

  std::string fmodelname;

  const std::vector<float> endPt1 = {-9999., -9999., -9999.};
  const std::vector<float> endPt2 = {-9999., -9999., -9999.};

  std::vector<anab::CosmicTagID_t> ctlabelvec;

  BoosterHandle h_booster;
};


DecisionTreeID::DecisionTreeID(fhicl::ParameterSet const & p) : EDProducer{p}
{
  ftrackmodulelabel     = p.get<std::string>("TrackModuleLabel");
  fmodelname            = p.get<std::string>("ModelFile");

  fTrackFeatures.Configure(p.get<fhicl::ParameterSet>("FeaturesConfig"));
  ctlabelvec = fTrackFeatures.CosmicTagLabel();
   
  produces< std::vector<anab::CosmicTag> >();  
  produces< art::Assns<recob::Track,anab::CosmicTag> >();  

  // create booster and load in model
  //std::string _env = std::getenv("MRB_INSTALL");
  //std::string _filename = _env+"/"+fmodelname;
  //std::string _filename = fmodelname;
  std::string _filename;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fmodelname,_filename);
  std::cout << "Loading model from " << _filename.c_str() << std::endl;
  DMatrixHandle h_tmp;
  XGBoosterCreate( &h_tmp, 0, &h_booster );
  XGBoosterSetParam( h_booster, "seed", "0" );
  int rc = XGBoosterLoadModel( h_booster, _filename.c_str() );
  if(rc == -1)
  {
    throw std::runtime_error("Failed to load XGBoost model.");
  }
}

void DecisionTreeID::produce(art::Event & e)
{

  // Make the ouput vector and assoc.
  std::unique_ptr< std::vector<anab::CosmicTag> > cosmicTagTrackVector(new std::vector<anab::CosmicTag>);
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag> > assnOutCosmicTagTrack(new art::Assns<recob::Track,anab::CosmicTag>);

  // recover handle for tracks that we want to analyze
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(ftrackmodulelabel, trackVecHandle);

  if(trackVecHandle.isValid())
  {
    std::vector< std::vector<float> > evtdata = fTrackFeatures.CreateFeatures(e, trackVecHandle);

    // Loop over input tracks and create associations
    for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
    {

      std::vector<float> trkdata = evtdata.at(trkIdx);
      int ncols = trkdata.size();

      // put data in DMatrix format
      DMatrixHandle h_data;
      XGDMatrixCreateFromMat(&trkdata[0], 1, ncols, NAN, &h_data);

      // predict classes
      bst_ulong out_len;
      const float *preds;
      
      XGBoosterPredict( h_booster, h_data, 0,0, &out_len, &preds);

      art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);
      for(size_t itag=0; itag < out_len; ++itag)
      {
        cosmicTagTrackVector->emplace_back(endPt1, endPt2, preds[itag], ctlabelvec[itag]);
        util::CreateAssn(*this, e, *cosmicTagTrackVector, trk, *assnOutCosmicTagTrack);
      }

      XGDMatrixFree(h_data);

    }
  }
  e.put(std::move(cosmicTagTrackVector));
  e.put(std::move(assnOutCosmicTagTrack));
}



DEFINE_ART_MODULE(DecisionTreeID)
