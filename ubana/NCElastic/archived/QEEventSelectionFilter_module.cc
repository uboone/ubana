////////////////////////////////////////////////////////////////////////
// Class:       QEEventSelectionFilter
// Plugin Type: filter (art v2_05_00)
// File:        QEEventSelectionFilter_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "TString.h"
#include "TTree.h"

#include "ubana/NCElastic/Algorithms/LREvtReconstruction.h"

#include <memory>
#include <cmath>
#include <stdexcept>

class QEEventSelectionFilter;

class QEEventSelectionFilter : public art::EDFilter {
public:
  explicit QEEventSelectionFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  QEEventSelectionFilter(QEEventSelectionFilter const &) = delete;
  QEEventSelectionFilter(QEEventSelectionFilter &&) = delete;
  QEEventSelectionFilter & operator = (QEEventSelectionFilter const &) = delete;
  QEEventSelectionFilter & operator = (QEEventSelectionFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Pointer to the feature algorithm
  ::qeselection::LREvtReconstruction fEvtReconstruction;

  float       _minLR;

  float _reco_q2,_reco_Tp,_reco_Rp,_reco_LR;
  int _reco_TID;
  float _reco_TH,_reco_PH;
  float _reco_nux,_reco_nuy,_reco_nuz;

};


QEEventSelectionFilter::QEEventSelectionFilter(fhicl::ParameterSet const & p)
{
  _minLR            = p.get<float>("MinimumLRScore",0.9);

  fEvtReconstruction.Configure(p.get<fhicl::ParameterSet>("SelectionConfig"));

}

bool QEEventSelectionFilter::filter(art::Event & e)
{

  std::vector<float> recovec = fEvtReconstruction.ReconstructEvent(e);

  _reco_q2   = recovec.at(0);
  _reco_Tp   = recovec.at(1);
  _reco_Rp   = recovec.at(2);
  _reco_LR   = recovec.at(3);
  _reco_TID  = recovec.at(4);
  _reco_TH   = recovec.at(5);
  _reco_PH   = recovec.at(6);
  _reco_nux  = recovec.at(7);
  _reco_nuy  = recovec.at(8);
  _reco_nuz  = recovec.at(9);


  if(_reco_LR < _minLR)
  {
    return false;
  }
  return true;
}

DEFINE_ART_MODULE(QEEventSelectionFilter)
