////////////////////////////////////////////////////////////////////////
// Class:       UsingAssnsAna
// Plugin Type: analyzer (art v2_11_03)
// File:        UsingAssnsAna_module.cc
//
// Generated at Wed Sep  5 15:05:25 2018 by Adam Lister using cetskelgen
// from cetlib version v3_03_01.
// 
// You can generate the skeleton of a LArSoft module by using the 
// cetskelgen command. This module was generated by using the command:
// 
// cetskelgen analyzer UsingAssnsAna
//
// This generates an _analyzer_ module with the name 
// UsingAssnsAna_module.cc. Analyzer modules access information in 
// the artroot file and can print to screen or create ROOT trees. 
// There are also _producer_ modules, which add data products to the
// artroot file, and _filter_ modules, which filter events based on
// criteria defined by the user.
//
// This module will simply access the recob::Track data 
// product in a MCC8 simulated event, and find the hits which are 
// associated to the track, meaning the hits used in construction of
// the track. It'll then save those variables to a simple ROOT tree.
//
// Tested against v07_04_00.
//
////////////////////////////////////////////////////////////////////////

// includes generated by default with cetskelgen
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Tracking down includes can be a pain, the way I usually do it is
// by going to the LArSoft doxygen, looking for the relevant .h
// file, and looking at the directory structure on the left.
//
// It's not perfect, but it works most of the time, unless there's
// been any re-organisation of header files.

// ART includes
#include "art_root_io/TFileService.h" // used for ROOT file
#include "canvas/Persistency/Common/FindManyP.h" // used for assns

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include "TTree.h"

class UsingAssnsAna;


class UsingAssnsAna : public art::EDAnalyzer {

  // By convention, LArSoft modules are completely self contained,
  // meaning that there's no .h file for a _module.cc
  // instead, there's this section at the top for variables and 
  // function definitions.
  //
  // You can force LArSoft to create a separate header files
  // (see cetskelgen --help), but it's ill-advised.

  public:
    explicit UsingAssnsAna(fhicl::ParameterSet const & p);

    UsingAssnsAna(UsingAssnsAna const &) = delete;
    UsingAssnsAna(UsingAssnsAna &&) = delete;
    UsingAssnsAna & operator = (UsingAssnsAna const &) = delete;
    UsingAssnsAna & operator = (UsingAssnsAna &&) = delete;

    void analyze(art::Event const & e) override;

    // this is an optional function, which can be generated automatically
    // using the '-e beginJob' flag with cetskelgen
    void beginJob() override;

    // function to clear out vectors and set dummy variables every event
    void resetVariables();

  private:

    // initialise service handles
    // services basically contain a lot of functions which
    // make your life easier. 
    // There are a few Services you might want to use as an analyzer:
    // * art::TFileService
    //   -- used for saving information to root files (Trees, histograms, etc.)
    // * cheat::BackTracker
    //   -- used to do reco-true matching
    // * geo::Geometry
    //   -- used to access detector geometry information (length, height, 
    //      convert ticks to X position, etc.)
    // * sim::LArG4Parameters
    //   -- used to access information about LAr
    art::ServiceHandle< art::TFileService > tfs;

    // defining private member variables here
    TTree* out_tree;

    int run;
    int sub_run;
    int event;
    // going to set the tree up such that one entry
    // is one event, so anything with multiple
    // objects per event (i.e. tracks) needs a vector
    std::vector<double> track_length_fromptrvec_v;
    std::vector<double> track_theta_fromptrvec_v;
    std::vector<double> track_phi_fromptrvec_v;

    // each track has a set of hits associated to it, so
    // define these to be a vector of vectors
    std::vector< std::vector< double > > hit_peak_time_v;
    std::vector< std::vector< double > > hit_integral_v;

    // fhicl parameters
    std::string fTrackLabel;
    std::string fTrackHitAssnLabel;

    // other variables
    double track_length;
    double track_theta;
    double track_phi;
    double hit_peak_time;
    double hit_integral;
};


// This function is the constructor for the module. It's going to be used
// to read in fhicl parameters from the associated fhicl file. 
UsingAssnsAna::UsingAssnsAna(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{

  // this is how you access variables from fhicl files. The first argument 
  // says "get the fhicl parameter that is called TrackLabel and treat it 
  // as a std::string", the second argument is the default argument in case
  // the fhicl parameter isn't set.
  fTrackLabel = p.get<std::string>("TrackLabel", "pandoraNu::McRecoStage2");
  fTrackHitAssnLabel = p.get<std::string>("TrackHitAssnLabel", "pandoraNu::McRecoStage2");

}

// this function is run once at the beginning of a job, not once per event
// we'll use this to setup the ROOT file
void UsingAssnsAna::beginJob()
{

  // the general usage for making objects with tfs is 
  // tfs->make<OBJ>("op1", "op2", "op3" ...), 
  // where the options are the variables needed for the 
  // constructor of the OBJ
  out_tree = tfs->make<TTree>("out_tree"      , "out_tree");
  out_tree->Branch("track_length_fromptrvec_v", "std::vector<double>", &track_length_fromptrvec_v);
  out_tree->Branch("track_theta_fromptrvec_v" , "std::vector<double>", &track_theta_fromptrvec_v);
  out_tree->Branch("track_phi_fromptrvec_v"   , "std::vector<double>", &track_phi_fromptrvec_v);
  out_tree->Branch("hit_peak_time_v"     , "std::vector< std::vector<double> >", &hit_peak_time_v);
  out_tree->Branch("hit_integral_v"      , "std::vector< std::vector<double> >" , &hit_integral_v);

}

// this acts as an event loop, everything inside these brackets are 
// applied to each event
void UsingAssnsAna::analyze(art::Event const & e)
{

  resetVariables();

  // the event itself contains some auxilliary information which we can make use of
  run = e.run();
  sub_run = e.subRun();
  event = e.event();

  std::cout << "Processing event " << run << "." << sub_run << "." << event << std::endl;

  // The standard way to access a data product from an artroot file.
  // How do you know what data products you can access? Check
  // whats in the file! You can do this by running
  // lar -c eventdump.fcl my_input_artroot_file.root
  art::Handle< std::vector<recob::Track> > trackHandle;

  // getByLabel will fill the track handle with a
  // std::vector<recob::Track> with the label defined in
  // fTrackLabel (which comes from the fhicl file)
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackPtrVector;
  art::fill_ptr_vector(trackPtrVector, trackHandle);
 
  // we're now going to build a "Map" between all of the tracks and all of the hits
  // in the event. We're using art::FindManyP which returns a vector of art::Ptrs to 
  // objects. There is also FindMany, FindOne and FindOneP, but this is a good 
  // catch-all example
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fTrackHitAssnLabel);
  
  for (size_t i = 0; i < trackPtrVector.size(); i++){
    
    art::Ptr<recob::Track> thisTrack = trackPtrVector.at(i);

    // we have an art::Ptr to an object. The accessor for this is ->, just
    // like a regular pointer, i.e:

    track_length = thisTrack->Length();
    track_theta  = thisTrack->Theta();
    track_phi    = thisTrack->Phi();

    track_length_fromptrvec_v.push_back(track_length);
    track_theta_fromptrvec_v.push_back(track_theta);
    track_phi_fromptrvec_v.push_back(track_phi);

    // using the map we build earlier using art::FindManyP, we can get all of the 
    // hits associatd withh this track.
    // We do this by using the key of the art::Ptr (note this is a property of the ptr,
    // not the object). If you don't have an art::Ptr to an object, this becomes a little
    // more complicated.
    //
    // IN THE CASE YOU DO NOT HAVE AN ART PTR then you can create a ptr vector to the
    // track collection, as above, and then loop over the tracks and find the 
    // art::Ptr<recob::Track> which has the same track->ID(), and then get the ptr key
    // of that track.
    //
    // IF USING PANDORA then rather than use the key of the Ptr you _can_ use the 
    // thisTrack->ID() method as they are the same by convention. This is not garunteed
    // for other producers becuase the ID is user-set.
    std::vector< art::Ptr< recob::Hit > > hitCollection = hitsFromTrack.at(thisTrack.key());

    // we now have a simple vector we can loop over!
    std::vector<double> hit_peak_time_tmpvec = {};
    std::vector<double> hit_integral_tmpvec = {};
    for (size_t i_hit = 0; i_hit < hitCollection.size(); i_hit++){

      art::Ptr<recob::Hit> thisHit = hitCollection.at(i_hit);

      hit_peak_time = (double)(thisHit->PeakTime()); 
      hit_integral = (double)(thisHit->Integral());

      hit_peak_time_tmpvec.push_back(hit_peak_time);
      hit_integral_tmpvec.push_back(hit_integral);

    }

    hit_peak_time_v.push_back(hit_peak_time_tmpvec);
    hit_integral_v.push_back(hit_integral_tmpvec);

  }

  out_tree->Fill();

}


void UsingAssnsAna::resetVariables()
{

  run = -1;
  sub_run = -1;
  event = -1;

  track_length_fromptrvec_v.resize(0);
  track_theta_fromptrvec_v.resize(0);
  track_phi_fromptrvec_v.resize(0);
  hit_peak_time_v.resize(0);
  hit_integral_v.resize(0);

}


DEFINE_ART_MODULE(UsingAssnsAna)
