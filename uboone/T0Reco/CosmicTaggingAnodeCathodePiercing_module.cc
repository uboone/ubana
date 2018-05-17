////////////////////////////////////////////////////////////////////////
// Class:       CosmicTaggingAnodeCathodePiercing
// Module Type: producer
// File:        CosmicTaggingAnodeCathodePiercing_module.cc
//
// David Caratelli - davidc1@fnal.gov   - July 13 2016
// Chris Barnes    - barnchri@umich.edu 
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

// services etc...
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT
#include "TVector3.h"

// C++
#include <memory>
#include <iostream>
#include <utility>

class CosmicTaggingAnodeCathodePiercing;

class CosmicTaggingAnodeCathodePiercing : public art::EDProducer {
public:
  explicit CosmicTaggingAnodeCathodePiercing(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicTaggingAnodeCathodePiercing(CosmicTaggingAnodeCathodePiercing const &) = delete;
  CosmicTaggingAnodeCathodePiercing(CosmicTaggingAnodeCathodePiercing &&) = delete;
  CosmicTaggingAnodeCathodePiercing & operator = (CosmicTaggingAnodeCathodePiercing const &) = delete;
  CosmicTaggingAnodeCathodePiercing & operator = (CosmicTaggingAnodeCathodePiercing &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  // producer of 3D reconstructed track to be used
  std::string fTrackProducer;

  // producer of reconstructed optical flashes
  std::string fFlashProducer;
  
  // pfparticle producer
  std::string fPFPartProducer;

  // set "resolution". How far away from the detector bounds
  // do we want to be to make a claim.
  double fTPCResolution; // [cm]

  // drift velocity // cm / us
  double fDriftVelocity;

  // tag which types of tracks to reconstruct
  bool top2side;
  bool front2side;
  bool back2side;
  bool side2bottom;
  bool side2front;
  bool side2back;

  // debug (verbose) mode?
  bool _debug;

  // define top, bottom, front and back boundaries of TPC
  double _TOP, _BOTTOM, _FRONT, _BACK;
  
  // vector to hold flash-times for the event
  std::vector<double> _flash_times;
  std::vector<size_t> _flash_idx_v;

  // detector width [drift-coord]
  double _det_width; // [cm]

  // time resolution for flashes to match [us] (separate values for Anode and Cathode)
  double fTimeResA, fTimeResC;
  
  // minimum PE threshold for flash to make it into match
  double fPEmin;

  // time-adjustment (us) to align reconstructed T0 for Anode and Cathode-crossing tracks to reconstructed flash-times
  double fRecoT0TimeOffsetA, fRecoT0TimeOffsetC;

  // time-bounds for ADC truncation peak removal (negative)
  double fT0negMin, fT0negMax;
  // (positive pek)
  double fT0posMin, fT0posMax;

  // functions to be used throughout module
  bool   TrackEntersTop     (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersFront   (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersBack    (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersAnode   (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersSide    (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsBottom   (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsFront    (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsBack     (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsAnode    (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsSide     (const std::vector<TVector3>& sorted_trk);

  // functions to be used for organization in the module
  void   SortTrackPoints      (const recob::Track& track,
			       std::vector<TVector3>& sorted_trk);
  double GetEnteringTimeCoord (const std::vector<TVector3>& sorted_trk);
  double GetExitingTimeCoord  (const std::vector<TVector3>& sorted_trk);

  // validate flash matching by requiring PMT flash
  std::pair<double,size_t> FlashMatch(const double reco_time);

};


CosmicTaggingAnodeCathodePiercing::CosmicTaggingAnodeCathodePiercing(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  //produces< std::vector< anab::T0 > >();
  produces< std::vector< anab::CosmicTag > >();
  //produces< art::Assns <recob::PFParticle, recob::OpFlash> >();
  //produces< art::Assns <recob::PFParticle, anab::T0> >();
  produces< art::Assns <recob::PFParticle, anab::CosmicTag> >();

  fPFPartProducer    = p.get<std::string>("PFPartProducer");
  fTrackProducer     = p.get<std::string>("TrackProducer");
  fFlashProducer     = p.get<std::string>("FlashProducer");
  fTPCResolution     = p.get<double>     ("Resolution");
  fTimeResA          = p.get<double>     ("TimeResA");
  fTimeResC          = p.get<double>     ("TimeResC");
  fRecoT0TimeOffsetA = p.get<double>     ("RecoT0TimeOffsetA");
  fRecoT0TimeOffsetC = p.get<double>     ("RecoT0TimeOffsetC");
  fT0negMin          = p.get<double>     ("T0negMin");
  fT0negMax          = p.get<double>     ("T0negMax");
  fT0posMin          = p.get<double>     ("T0posMin");
  fT0posMax          = p.get<double>     ("T0posMax");
  fPEmin             = p.get<double>     ("PEmin");
  top2side           = p.get<bool>       ("top2side");
  front2side         = p.get<bool>       ("front2side");
  back2side          = p.get<bool>       ("back2side");
  side2bottom        = p.get<bool>       ("side2bottom");
  side2front         = p.get<bool>       ("side2front");
  side2back          = p.get<bool>       ("side2back");
  _debug             = p.get<bool>       ("debug");
  
  // get boundaries based on detector bounds
  auto const* geom = lar::providerFrom<geo::Geometry>();

  _TOP    =   geom->DetHalfHeight() - fTPCResolution;
  _BOTTOM = - geom->DetHalfHeight() + fTPCResolution;
  _FRONT  =   fTPCResolution;
  _BACK   =   geom->DetLength() - fTPCResolution;
  
  _det_width = geom->DetHalfWidth() * 2;

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  fDriftVelocity = _detp -> DriftVelocity(efield,temp);

}

void CosmicTaggingAnodeCathodePiercing::produce(art::Event & e)
{
  if (_debug) { std::cout << "NEW EVENT" << std::endl; }

  _flash_times.clear();
  _flash_idx_v.clear();

  // produce data-products and associations
  //std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< std::vector<anab::CosmicTag> > CosmicTag_v(new std::vector<anab::CosmicTag>);
  //std::unique_ptr< art::Assns <recob::PFParticle, recob::OpFlash> > pfp_flash_assn_v( new art::Assns<recob::PFParticle, recob::OpFlash> );
  //std::unique_ptr< art::Assns <recob::PFParticle, anab::T0> >       pfp_t0_assn_v   ( new art::Assns<recob::PFParticle, anab::T0>       );
  std::unique_ptr< art::Assns <recob::PFParticle, anab::CosmicTag> > pfp_cosmictag_assn_v ( new art::Assns<recob::PFParticle, anab::CosmicTag>  );


  // load Flash
  if (_debug) { std::cout << "loading flash from producer " << fFlashProducer << std::endl; }
  art::Handle<std::vector<recob::OpFlash> > flash_h;
  e.getByLabel(fFlashProducer,flash_h);

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    throw std::exception();
  }

  // load PFParticles for which T0 reconstruction should occur
  if (_debug) { std::cout << "loading PFParticles from producer " << fPFPartProducer << std::endl; }
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(fPFPartProducer,pfp_h);

  // make sure pfparticles look good
  if(!pfp_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate PFParticles!"<<std::endl;
    throw std::exception();
  }

  // grab tracks associated with PFParticles
  art::FindMany<recob::Track> pfp_track_assn_v(pfp_h, e, fTrackProducer);
  if (_debug)
    std::cout << "There are " << pfp_track_assn_v.size() << " pfpart -> track associations" << std::endl;

  std::vector<art::Ptr<recob::PFParticle> > PFPVec;
  art::fill_ptr_vector(PFPVec, pfp_h);

  // prepare a vector of optical flash times, if flash above some PE cut value

  size_t flash_ctr = 0;
  for (auto const& flash : *flash_h){
    if (flash.TotalPE() > fPEmin){
      _flash_times.push_back( flash.Time() );
      _flash_idx_v.push_back(flash_ctr);
      if (_debug) { std::cout << "\t flash time : " << flash.Time() << ", PE : " << flash.TotalPE() << std::endl; }
    }
    flash_ctr += 1;
  }// for all flashes

  if (_debug) { std::cout << "Selected a total of " << _flash_times.size() << " OpFlashes" << std::endl; }

  for (size_t i=0; i < PFPVec.size(); i++) {

    auto pfp = PFPVec.at(i);

    // grab associated tracks
    const std::vector<const recob::Track*>& track_v = pfp_track_assn_v.at(i);

    if (_debug) {
      std::cout << "Looping through pfpart number " << i << std::endl;
      std::cout << "PFPart has " << track_v.size() << " tracks associated" << std::endl;
    }

    // Declare the variable 'trkT' up here so that I can continue and not fill the t0 object if trkT is still equal to 0
    double trkT = 0.;
    // dT w.r.t. closest flash
    //double deltaT = 0.;
    // matched flash index
    //size_t flashidx = 0;

    // how many times has a time been reconstructed for this PFPart?
    size_t nt0s = 0;
    
    // loop through all tracks in PFParticle
    // if ANY of them matched -> tag entire PFParticle 
    for (auto const *track : track_v) {
      
      // get sorted points for the track object [assuming downwards going]
      std::vector<TVector3> sorted_trk;
      SortTrackPoints(*track,sorted_trk);

      // keep track of whether it goes thorugh the anode or cathode
      bool anode = 0;
      
      // 1st category: tracks which ENTER SIDE
      if ( TrackEntersSide(sorted_trk) == true ) {
	
	if (_debug) std::cout << "\t track enters side" << std::endl;
	
	// we are not done. We need to check that the track either: 1) Exits the bottom. 2) exits the front or 3) exits the back of the TPC.
	bool tagged = false;
	
	// tracks that exit the bottom
	if ( (TrackExitsBottom(sorted_trk) == true) and (side2bottom == true) ) {
	  tagged = true;
	  if (_debug) std::cout << "\t track exits bottom" << std::endl;
	}
	// tracks that exit the front
	if ( (TrackExitsFront(sorted_trk) == true) and (TrackEntersFront(sorted_trk) == false) and (side2front == true) ) {
	  tagged = true;
	  if (_debug) std::cout << "\t track exits front" << std::endl;
	}
	// tracks that exit the back
	if ( (TrackExitsBack(sorted_trk) == true) and (TrackExitsBack(sorted_trk) == false) and (side2back == true) ) {
	  tagged = true;
	  if (_debug) std::cout << "\t track exits back" << std::endl;
	}
	
	// has either of these 3 conditions been met? if no, skip this track
	if (tagged == false) continue;
	
	// figure out if it enters the anode or cathode                                                                                                              
	bool enters_anode = TrackEntersAnode(sorted_trk);
	
	// get the X coordinate of the point piercing the anode/cathode (upon ENTERING) 
	double trkX = GetEnteringTimeCoord(sorted_trk);
	
	// reconstruct track T0 w.r.t. trigger time                                                 	
	// The 'trkX' enters on the anode, the side of the TPC with a lower x value than the cathode
	if (enters_anode){
	  trkT = trkX / fDriftVelocity + fRecoT0TimeOffsetA;
	  anode = 1;
	}
	// This will also give a small T0 value, because the cathode is a distance of _det_width from the anode
	else{
	  trkT = (trkX - _det_width) / fDriftVelocity + fRecoT0TimeOffsetC; 
	  anode = 0;
	}
	
      }// if the track enters the side
      
      // case in which the track exits the side
      if (TrackExitsSide(sorted_trk) == true) {
	
	if (_debug) std::cout << "\t track exits side" << std::endl;
	
	// we are not done. We need to check that the track either: 1) Enters the bottom. 2) enters the front or 3) enters the back of the TPC.
	bool tagged = false;
	
	// track enters the top
	if ( (TrackEntersTop(sorted_trk) == true) and (top2side == true) ) {       
	  tagged = true;
	  if (_debug) std::cout << "\t track enters the top" << std::endl;
	}
	
	if ( (TrackEntersFront(sorted_trk) == true) and (TrackExitsFront(sorted_trk) == false) and (front2side == true) ) {
	  tagged = true;
	  if (_debug) std::cout << "\t track enters front" << std::endl;
	}
	
	if ( (TrackEntersBack(sorted_trk) == true) and (TrackExitsBack(sorted_trk) == false) and (back2side == true) ) {
	  tagged = true;
	  if (_debug) std::cout << "\t track enters back" << std::endl;
	}
	
	// has either of these 3 conditions been met? if no, skip this track
	if (tagged == false) continue;
	
	// figure out if it enters the anode or cathode                                                                                                              
	bool enters_anode = TrackEntersAnode(sorted_trk);
	
	// get the X coordinate of the point piercing the anode/cathode (upon ENTERING) 
	double trkX = GetEnteringTimeCoord(sorted_trk);
	
	// reconstruct track T0 w.r.t. trigger time                                                                                                                              
	
	// The 'trkX' enters on the anode, the side of the TPC with a lower x value than the cathode
	if (enters_anode){
	  trkT = trkX / fDriftVelocity + fRecoT0TimeOffsetA;
	  anode = 1;
	}
	// This will also give a small T0 value, because the cathode is a distance of _det_width from the anode
	else{
	  trkT = (trkX - _det_width) / fDriftVelocity + fRecoT0TimeOffsetC; 
	  anode = 0;
	}
	
      }// if the track exits the side
      
      
      if (_debug)
	std::cout << "\t this track has a reconstructed time = " << trkT << std::endl;
      
      // if the time does not match one from optical flashes -> don't reconstruct
      auto const& flash_match_result = FlashMatch(trkT);
      // flash_match_result is std::pair
      // 1st element is dt w.r.t. closest flash of light in PMTs
      // 2nd element is index of PMT flash matched to
      if ( (flash_match_result.first > fTimeResA) && (anode == 1) )
	continue;
      if ( (flash_match_result.first > fTimeResC) && (anode == 0) )
	continue;

      if (_debug) {std::cout << "\t matched to flash with dt = " << flash_match_result.first << std::endl; }
      
      // some T0 reconstructed values mean that the track hits were truncated due
      // to ADC waveform truncation. They can be identified by the distribution
      // of reconstructed T0s
      if ( (trkT > fT0negMin) && (trkT < fT0negMax) ) continue;
      if ( (trkT > fT0posMin) && (trkT < fT0posMax) ) continue;

      if (_debug) { std::cout << "\t track passes all cuts. T0 successfully reconstructed!" << std::endl; }

      nt0s += 1;
      
      // DON'T CREATE the t0 object unless the reconstructed t0 is some value other than 0
      
      //deltaT   = flash_match_result.first;
      //flashidx = flash_match_result.second;
      
    }// for all tracks in this PFParticle
    
    // create T0 object with this information!
    //anab::T0 t0(trkT, 0, deltaT);
    //T0_v->emplace_back(t0);
    //util::CreateAssn(*this, e, *T0_v, pfp, *pfp_t0_assn_v);
    
    // create CosmicTag object
    if ( (nt0s == 1) && (trkT != 0) ) {
      anab::CosmicTag ct(trkT);
      CosmicTag_v->emplace_back(ct);
      util::CreateAssn(*this, e, *CosmicTag_v, pfp, *pfp_cosmictag_assn_v);
      if (_debug) { std::cout << "\t found T0 for PFParticle. Time = " << trkT << std::endl; }
    }
    // get pointer to individual track
    // TMP const art::Ptr<recob::Track>   trk_ptr(track_h,trk_ctr-1);
    //const art::Ptr<recob::OpFlash> flash_ptr(flash_h, flashidx );
    //if (_debug)
    // std::cout << "\t mathed to flash w/ index " << flashidx << " w/ PE " << flash_ptr->TotalPE() << " and time " << flash_ptr->Time() << " vs reco time " << trkT << std::endl;
    //pfp_flash_assn_v->addSingle( pfp, flash_ptr );
    
  }// for all PFParticles
  
  //e.put(std::move(T0_v));
  //e.put(std::move(pfp_t0_assn_v));
  e.put(std::move(CosmicTag_v));
  e.put(std::move(pfp_cosmictag_assn_v));
  if (_debug)
    std::cout << "create track flash association " << std::endl;
  //e.put(std::move(pfp_flash_assn_v));


}

std::pair<double,size_t> CosmicTaggingAnodeCathodePiercing::FlashMatch(const double reco_time){
  
  // loop through all reco'd flash times and see if one matches
  // the reco time from the track
  double dt_min = 4000.; // us
  size_t idx_min = _flash_times.size();

  for (size_t i=0; i < _flash_times.size(); i++){
    auto const& time = _flash_times[i];
    double dt = fabs(time - reco_time);
    if (dt < dt_min){
      dt_min  = dt;
      idx_min = _flash_idx_v[i];
    }
  }

  std::pair<double,size_t> ret(dt_min,idx_min);
  return ret;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackEntersTop(const std::vector<TVector3>& sorted_trk)
{
  // check that the first point in the track
  // pierces the top boundary of the TPC
  // This track either will pierce the top of the TPC or is just about to (the '_TOP' variable is just below the actual coordinate position of the top in Y)

  if (sorted_trk.at(0).Y() > _TOP)
    return true;

  return false;
}


bool CosmicTaggingAnodeCathodePiercing::TrackEntersFront(const std::vector<TVector3>& sorted_trk)
{

  // Determine if the track enters the                                                                                                                                       
  // front of the TPC based on if the position                                                                                                                                
  // of its initial Z-coordinate is less than 
  // the location of the front of the TPC in Z                                                                                                                               
  
  // First define 'top_pt' to mean the point at the start of the track                                                                                            
  auto const& top_pt = sorted_trk.at(0);

  if (top_pt.Z() < _FRONT)
    return true;

  // I may include the case in which I check                                                                                                                       
  // the y-coordinates as well, but I will not                                                                                                 
  // implement that at this time                                                                                                                                 
  
  // If this condition is not satisfied, then return 'false' (the track was not determined                                                                  
  // within resolution to enter the front of the TPC)                                                                                                            
  return false;
}


bool CosmicTaggingAnodeCathodePiercing::TrackEntersBack(const std::vector<TVector3>& sorted_trk)
{

  // Determines if the track enters the                                                                                              
  // back of the TPC based on if the position                                                                                       
  // of its initial Z-coordinate is greater                                                                                          
  // than the location of the back of the                                                                                                                          
  // TPC in Z                                                                                                                                                           
  
  // First define 'top_pt' to mean the point at the start of the track                                                                     
  auto const& top_pt = sorted_trk.at(0);

  if (top_pt.Z() > _BACK)
    return true;

  // If this condition is not satisfied, then return 'false' (the track was not determined                                                       
  // within resolution to enter the back of the TPC)                                                                                                               
  return false;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackEntersAnode(const std::vector<TVector3>& sorted_trk)
{

  // we know the track enters either the                                                                                                                               
  // anode or cathode                                                                                                                                                        
  // at this point figure out                                                                                                                                                
  // if it ENTERS the ANODE or CATHODE                                                                                         
  // ANODE: top point must be at lower X-coord                                                                                                
  // than bottom point                                                                                              
  // CATHODE: top point must be at larger X-coord                                                                                                                              
  // than bottom point
  // assume track has already been sorted                                                                                                                                     
  // such that the 1st point is the most elevated in Y coord.                                                                                                                 
  // return TRUE if passes the ANODE                                                                                                                                                
  
  auto const& top    = sorted_trk.at(0);
  auto const& bottom = sorted_trk.at( sorted_trk.size() - 1 );

  if (top.X() < bottom.X())
    return true;
  
  return false;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackEntersSide(const std::vector<TVector3>& sorted_trk)
{
  
  // check that the top-most point                                                                                                                                            
  // is not on the top of the TPC                                                                                                                                              
  // nor on the front & back of the TPC                                                                                                                                           
  
  auto const& top_pt = sorted_trk.at(0);

  // if highest point above the TOP -> false                                                                                                                                   
  if (top_pt.Y() > _TOP)
    return false;

  // if highest point in Z close to front or back                                                                                                                              
  // -> FALSE                                                                                                                                                                 
  if ( (top_pt.Z() < _FRONT) or (top_pt.Z() > _BACK) )
    return false;


  // If the function makes it this far, then it will enter through one of the sides of the TPC
  return true;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackExitsBottom(const std::vector<TVector3>& sorted_trk)
{

  // check that the last point in the track                                                                                                                                    
  // pierces the bottom boundary of the TPC                                                                                                                                   
  if ( sorted_trk.at( sorted_trk.size() - 1).Y() < _BOTTOM )
    return true;

  return false;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackExitsFront(const std::vector<TVector3>& sorted_trk)
{

  // Determine if the track exits the                                                                                                                                      
  // front of the TPC based on if the position                                                                                                                             
  // of its final Z-coordinate is less than                                                                                                                      
  // the location of the front of the TPC in Z                                                                                                                           
  
  // First define 'bottom_pt' to mean the point at the end of the track                                                                                             
  auto const& bottom_pt = sorted_trk.at(sorted_trk.size() - 1);

  if (bottom_pt.Z() < _FRONT)
    return true;
  
  return false;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackExitsBack(const std::vector<TVector3>& sorted_trk)
{

  // Determine if the track exits the                                                                                            
  // front of the TPC based on if the position                                                                                      
  // of its final Z-coordinate is less than                                                                                    
  // the location of the front of the TPC in Z                                                                                                                          
  
  // First define 'bottom_pt' to mean the point at the end of the track                                              
  auto const& bottom_pt = sorted_trk.at(sorted_trk.size() - 1);

  if (bottom_pt.Z() > _BACK)
    return true;

  return false;
}


bool   CosmicTaggingAnodeCathodePiercing::TrackExitsAnode(const std::vector<TVector3>& sorted_trk)
{

  // Check, once it's known that the track doesn't exit out of the bottom, whether it's the anode or                                                                           
  // the cathode that it exits out of                                                                                                                                         
  // This can be done by direct analogy with the 'Anode' function (shown in this file as the 'TrackEntersAnode') function written by D. Caratelli                             
  // Define 'top' as the point at the start of the track, and 'bottom' as the point at the end of the track                                                                     

  auto const& top    = sorted_trk.at(0);
  auto const& bottom = sorted_trk.at(sorted_trk.size() - 1);

  // Check to see which point has a lower x coordinate                                                                                                                         
  // If the bottom does, then it exits out of the anode                                                                                                                       
  // If the top does, then it exits out of the cathode                                                                                                                          
  if (bottom.X() < top.X()) 
    return true;

  return false; // Otherwise, the top is less than the bottom, so the track ended closer to the cathode and exited there                                                          
}


bool   CosmicTaggingAnodeCathodePiercing::TrackExitsSide(const std::vector<TVector3>& sorted_trk)
{

  // check that the bottom-most point                                                                                                                                           
  // is not on the bottom of the TPC                                                                                                                                            
  // nor on the front & back of the TPC                                                                                                                                              

  auto const& bottom_pt = sorted_trk.at(sorted_trk.size() - 1);

  // if lowest point below the BOTTOM -> false                                                                                                            
  // Within this resolution, this means that it's likely that the track exited out of the bottom (at a point earlier on in the process than the last point) OR is just about to

  if (bottom_pt.Y() <  _BOTTOM)
    return false;

  // if lowest point in Z close to front or back                                                                                                         
  // -> FALSE                                                                                                                                              
  // If the the bottom point is less than the front, then the track has already pierced the front of the TPC and exited that way OR is likely just about to
  // If the bottom point is greater than the back, then the track has already pierced the back of the TPC and exited that way OR is likely just about to
  if ( (bottom_pt.Z() < _FRONT) or (bottom_pt.Z() > _BACK) )
    return false;

  return true;
}

void   CosmicTaggingAnodeCathodePiercing::SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_trk)
{

  // vector to store 3D coordinates of                                                                                                                                           
  // ordered track                              
  sorted_trk.clear();

  // take the reconstructed 3D track                                                                                                                                           
  // and assuming it is downwards                                                                                                    
  // going, sort points so that                                                                                                              
  // the track starts at the top                                                                                                     
  // which point is further up in Y coord?                                                                                                                  
  // start or end?                                                                                                                 
  auto const&N = track.NumberTrajectoryPoints();
  auto const&start = track.LocationAtPoint(0);
  auto const&end   = track.LocationAtPoint( N - 1 );

  // if points are ordered correctly                                                                                                                                       
  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint(i) );
  }
  
  // otherwise flip order                                                                                                                                                 
  else {
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint( N - i - 1) );
  }
}


double CosmicTaggingAnodeCathodePiercing::GetEnteringTimeCoord(const std::vector<TVector3>& sorted_trk)
{

  // get the drift-coordinate value                                                                                                             
  // associated with the point                                                                                                                                      
  // along the track piercing the anode / cathode                                                                                                 
  // ** WHEN the track enters the anode / cathode
  return sorted_trk.at(0).X();
}


double CosmicTaggingAnodeCathodePiercing::GetExitingTimeCoord(const std::vector<TVector3>& sorted_trk) 
{
  // get the drift-coordinate value                                                                                                                                           
  // associated with the point                                                                                                                                                
  // along the track piercing the anode / cathode                                                                                                                            
  // ** WHEN the track exits the anode / cathode
  return sorted_trk.at(sorted_trk.size() - 1).X();
}

DEFINE_ART_MODULE(CosmicTaggingAnodeCathodePiercing)
