/**
 * @file    AnalysisTree_module.cc
 * @brief   Module to create a TTree for analysis
 * @authors tjyang@fnal.gov, sowjanyag@phys.ksu.edu
 *
 *
 * Current implementation
 * =======================
 * 
 * There is one tree only, with one set of branches for each tracking algorithm.
 * The data structure which hosts the addresses of the tree branches is
 * dynamically allocated on demand, and it can be optionally destroyed at the
 * end of each event.
 * The data structure (AnalysisTreeDataStruct) directly contains the truth and
 * simulation information as C arrays. The data from tracking algorithms is the
 * largest, and it is contained in a C++ vector of structures (TrackDataStruct),
 * one per algorithm. These structures can also be allocated on demand.
 * Each of these structures is connected to a set of branches, one branch per
 * data member. Data members are vectors of numbers or vectors of fixed-size
 * C arrays. The vector index represents the tracks reconstructed by the
 * algorithm, and each has a fixed size pool for hits (do ROOT trees support
 * branches with more than one dimension with variable size?).
 * The data structures can assign default values to their data, connect to a
 * ROOT tree (creating the branches they need) and resize.
 * The AnalysisTreeDataStruct is constructed with as many tracking algorithms as
 * there are named in the module configuration (even if they are not backed by
 * any available tracking data).
 * By default construction, TrackDataStruct is initialized in a state which does
 * not allow any track (maximum tracks number is zero), and in such state trying
 * to connect to a tree has no effect. This is done so that the
 * AnalysisTreeDataStruct can be initialized first (and with unusable track data
 * structures), and then the TrackDataStruct instances are initialized one by
 * one when the number of tracks needed is known.
 * A similar mechanism is implemented for the truth information.
 * 
 * The "UseBuffers: false" mode assumes that on each event a new
 * AnalysisTreeDataStruct is created with unusable tracker data, connected to
 * the ROOT tree (the addresses of the available branches are assigned), then
 * each of the tracking algorithm data is resized to host the correct number
 * of reconstructed tracks and connected to the tree. Then the normal process of
 * filling the event data and then the tree take place. Finally, the whole
 * data structure is freed and the tree is left in a invalid state (branch
 * addresses are invalid). It could be possible to make the tree in a valid
 * state by resetting the addresses, but there is no advantage in that.
 * 
 * The "UseBuffers: true" mode assumes that on the first event a new
 * AnalysisTreeDataStruct is created and used just as in the other mode
 * described above. At the end of the first event, the data structure is left
 * around (and the tree is in a valid state). On the next event, all the
 * addresses are checked, then for each tracker the data is resized to
 * accommodate the right number of tracks for tis event. If the memory is
 * increased, the address will be changed. All the branches are reconnected to
 * the data structure, and the procedure goes on as normal.
 * 
 * Note that reducing the maximum number of tracks in a TrackDataStruct does not
 * necessarily make memory available, because of how std::vector::resize()
 * works; that feature can be implemented, but it currently has not been.
 * 
 * The BoxedArray<> class is a wrapper around a normal C array; it is needed
 * to be able to include such structure in a std::vector. This container
 * requires its objects to be default-constructible and copy-constructible,
 * and a C array is neither. BoxedArray<> is: the default construction leaves it
 * uninitialized (for speed reasons) while the copy construction is performed
 * as in a Plain Old Data structure (memcpy; really!).
 * 
 * 
 * Steps to reduce memory usage
 * -----------------------------
 * 
 * [x] create the data structure connected to the tree only when needed
 * [x] reduce the size of the elemental items (Double_t => Float_t could damage precision)
 * [x] create a different structure for each tracker, allocate only what needed
 * [x] use variable size array buffers for each tracker datum instead of [kMaxTrack]
 * [x] turn the truth/GEANT information into vectors
 * [ ] move hit_trkid into the track information, remove kMaxTrackers
 * [ ] turn the hit information into vectors (~1 MB worth), remove kMaxHits
 * [ ] fill the tree branch by branch
 * 
 * 
 * Quick guide to addition of branches
 * ====================================
 * 
 * A number of steps need to be taken to add an entry to the tree:
 * -# declare a variable for the datum
 * -# create a branch in the tree
 * -# assign the address of that branch to the datum
 * -# on each event: update all the data
 * -# on each event: fill the tree with the updated data
 * 
 * Depending on the data type, the following "reference guide" describes the
 * concrete implementation of each of the concepts explained here.
 * 
 * 
 * Declare a variable
 * -------------------
 * 
 * If the variable is a track variable, add a data member with the proper
 * "Data structure" in AnalysisTreeDataStruct::TrackDataStruct; if it is a
 * shower variable, add it to AnalysisTreeDataStruct::ShowerDataStruct; if
 * it is a PFParticle variable, add it to 
 * AnalysisTreeDataStruct::PFParticleDataStruct. Otherwise, add it directly 
 * to AnalysisTreeDataStruct.
 * 
 * | category of variable | where to add the variable                   |
 * | -------------------- | ------------------------------------------- |
 * | track                | AnalysisTreeDataStruct::TrackDataStruct     |
 * | shower               | AnalysisTreeDataStruct::ShowerDataStruct    |
 * | PFParticle           | AnalysisTreeDataStruct::PFParticleDataStruct|
 * | other                | AnalysisTreeDataStruct                      |
 * 
 * 
 * Create a branch/assign the branch address
 * ------------------------------------------
 * 
 * Assigning a branch an address is done by adding in the location described
 * in this table the "Create branch" command shown below.
 * The BranchCreator creates the branch if it does not exist, and then reassigns
 * the address.
 * 
 * | category of variable | where to add the variable                                   |
 * | -------------------- | ----------------------------------------------------------- |
 * | track                | AnalysisTreeDataStruct::TrackDataStruct::SetAddresses()     |
 * | shower               | AnalysisTreeDataStruct::ShowerDataStruct::SetAddresses()    |
 * | PFParticle           | AnalysisTreeDataStruct::PFParticleDataStruct::SetAddresses()|
 * | other                | AnalysisTreeDataStruct::SetAddresses()                      |
 * 
 * 
 * Prepare to update the content
 * ------------------------------
 * 
 * If the variable is a vector (like most of them are), add a proper resize()
 * call for it. This is needed for the track and shower variables.
 * The variables of simple `std::vector<>` types can be resized anywhere.
 * 
 * | category  | where to add the call                                 | call to be added                 |
 * | --------- | ----------------------------------------------------- | -------------------------------- |
 * | track     | AnalysisTreeDataStruct::TrackDataStruct::Resize()     | `VarName.resize(MaxTracks);`     |
 * | shower    | AnalysisTreeDataStruct::ShowerDataStruct::Resize()    | `VarName.resize(MaxShowers);`    |
 * | PFParticle| AnalysisTreeDataStruct::PFParticleDataStruct::Resize()| `VarName.resize(MaxPFParticles);`|
 * | other     | AnalysisTreeDataStruct::ClearLocalData()              | specific to the type             |
 * 
 * 
 * Then, assign the variable with an "invalid" value (e.g., -9999.); using a
 * FillWith() call (see below).
 * 
 * | category  | where to add the call                                |
 * | --------- | ---------------------------------------------------- |
 * | track     | AnalysisTreeDataStruct::TrackDataStruct::Clear()     |
 * | shower    | AnalysisTreeDataStruct::ShowerDataStruct::Clear()    |
 * | PFParticle| AnalysisTreeDataStruct::PFParticleDataStruct::Clear()|
 * | other     | AnalysisTreeDataStruct::ClearLocalData()             |
 * 
 * Finally, do what you need to fill the data, in the AnalysisTree::analyze()
 * method.
 * 
 * 
 * Fill the tree with the data event
 * ----------------------------------
 * 
 * This is already done. No need to be explicit about that.
 * 
 * 
 * 
 * Reference guide to addition of branches
 * ========================================
 * 
 * In template expressions, the type `T` is a type that can properly represent
 * the single value being stored in the branch (e.g., UInt_t, Double_t, etc.).
 * 
 * The `CreateBranch` in the "create branch" row is an instance of
 * BranchCreator constructed with a pointer to the ROOT tree; in its arguments,
 * "/T" must be replaced with the correct ROOT branch code.
 * 
 * 
 * Simple data types
 * ------------------
 * 
 * | Example        | number of generated particles                           |
 * | Data structure | `T`                                                     |
 * | Indices        | none                                                    |
 * | Fill call      | `VarName = Value;`                                      |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + "/T")` |
 * 
 * 
 * One-dimension vectors
 * ----------------------
 * 
 * | Example        | energy of the generated particles            |
 * | Data structure | `std::vector<T>`                             |
 * | Indices        | first: index of the item in the list         |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + "/T")` |
 * 
 * 
 * Track/shower information
 * -------------------------
 * 
 * Multiple track and shower collections (from different algorithms or algorithm
 * settings) live in different branches.
 * 
 * ### track/shower collection information
 * 
 * | Example        | Track ID                                                                  |
 * | Data structure | `TrackData_t<T>` / `ShowerData<T>`                                        |
 * | Index          | first: index of the track                                                 |
 * | Fill call      | `FillWith(VarName, Value);`                                               |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + NTracksIndexStr + "/T")` |
 * |                | `CreateBranch(BranchName, &VarName, BranchName + NShowerIndexStr + "/T")` |
 * 
 * 
 * ### plane-wide information for a track/shower
 * 
 * | Example        | Number of hits in the plane                                                  |
 * | Data structure | `PlaneData_t<T>`                                                             |
 * | Indices (2)    | first: index of the track/shower                                             |
 * |                | second: index of the plane (fixed dimension kNplanes)                        |
 * | Fill call      | `for (auto& elem: VarName) FillWith(elem, Value);`                           |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + NTracksIndexStr + "[3]/T")` |
 * |                | `CreateBranch(BranchName, &VarName, BranchName + NShowerIndexStr + "[3]/T")` |
 * 
 * 
 * ### information about hits in a track/shower
 * 
 * | Example        | dE/dX of a single hit                                                                                     |
 * | Data structure | `HitData_t<T>`                                                                                            |
 * | Indices (3)    | first: index of the track/shower                                                                          |
 * |                | second: index of the plane (fixed dimension kNplanes)                                                     |
 * |                | third: index of the hit (fixed dimension kMaxTrackHits/kMaxShowerHits)                                    |
 * | Fill call      | `for (auto& elem: VarName) FillWith(elem, Value);`                                                        |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/T")`  |
 * |                | `CreateBranch(BranchName, &VarName, BranchName + NShowerIndexStr + "[3]" + MaxShowerHitsIndexStr + "/T")` |
 * 
 * 
 * ### 3D coordinates of hits in a track/shower
 * 
 * | Example        | position coordinates of a hit                                                                                |
 * | Data structure | `HitCoordData_t<T>`                                                                                          |
 * | Indices (4)    | first: index of the track/shower                                                                             |
 * |                | second: index of the plane (fixed dimension kNplanes)                                                        |
 * |                | third: index of the hit (fixed dimension kMaxTrackHits/kMaxShowerHits)                                       |
 * |                | fourth: coordinate index (0 = x, 1 = y, 2 = z)                                                               |
 * | Fill call      | `for (auto& elem: VarName) FillWith(elem, Value);`                                                           |
 * | Create branch  | `CreateBranch(BranchName, &VarName, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "[3]/T")`  |
 * |                | `CreateBranch(BranchName, &VarName, BranchName + NShowerIndexStr + "[3]" + MaxShowerHitsIndexStr + "[3]/T")` |
 * 
 * 
 * Auxiliary detector information
 * -------------------------------
 * 
 * ### per-true-particle data
 * 
 * | Example        | number of auxiliary detectors a particle crossed                         |
 * | Data structure | `std::vector<T>`                                                         |
 * | Index          | first: index of the particle                                             |
 * | Fill call      | `FillWith(VarName, Value);`                                              |
 * | Create branch  | `CreateBranch(BranchName, VarName, BranchName + "[geant_list_size]/T");` |
 * 
 * ### detector data for a given particle
 * 
 * | Example        | energy observed by the detector                                                                   |
 * | Data structure | `AuxDetMCData_t<T>`                                                                               |
 * | Indices (2)    | first: index of the particle                                                                      |
 * |                | second: index of the detector                                                                     |
 * | Fill call      | `for (auto& elem: VarName) FillWith(elem, Value);`                                                |
 * | Create branch  | `CreateBranch(BranchName, VarName, BranchName + "[geant_list_size]" + MaxAuxDetIndexStr + "/T");` |
 * 
 * 
 ******************************************************************************/

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"

constexpr int kNplanes       = 3;     //number of wire planes
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTrackHits  = 2000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 15;    //number of trackers passed into fTrackModuleLabel
constexpr int kMaxVertices   = 500;    //max number of 3D vertices
constexpr int kMaxVertexAlgos = 10;    //max number of vertex algorithms
constexpr int kMaxFlashAlgos = 10;    //max number of flash algorithms
constexpr int kNOpDets = 32;          ///< number of optical detectors (PMTs)
constexpr unsigned short kMaxAuxDets = 4; ///< max number of auxiliary detector cells per MC particle
constexpr int kMaxFlashes    = 1000;   //maximum number of flashes
constexpr int kMaxShowerHits = 10000;  //maximum number of hits on a shower
constexpr int kMaxTruth      = 10;     //maximum number of neutrino truth interactions
constexpr int kMaxClusters   = 1000;   //maximum number of clusters;
constexpr int kMaxTicks   = 9600;   //maximum number of ticks (time samples)

constexpr int kMaxNDaughtersPerPFP = 100; //maximum number of daughters per PFParticle
constexpr int kMaxNClustersPerPFP  = 100; //maximum number of clusters per PFParticle
constexpr int kMaxNPFPNeutrinos    = 10;  //maximum number of reconstructed neutrino PFParticles

//constexpr int kMaxSysts = 1000;
//constexpr int kMaxWeights = 1000;

/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
  = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace microboone {

  /// Data structure with all the tree information.
  /// 
  /// Can connect to a tree, clear its fields and resize its data.
  class AnalysisTreeDataStruct {
  public:
    
    /// A wrapper to a C array (needed to embed an array into a vector)
    template <typename Array_t>
    class BoxedArray {
    protected:
      Array_t array; // actual data
      
    public:
      using This_t = BoxedArray<Array_t>;
      typedef typename std::remove_all_extents<Array_t>::type Data_t;
      
      BoxedArray() {} // no initialization
      BoxedArray(const This_t& from)
      { std::memcpy((char*) &(data()), (char*) &(from.data()), sizeof(Array_t)); }
      
      Array_t& data() { return array; }
      const Array_t& data() const { return array; }
      
      //@{
      /// begin/end interface
      static constexpr size_t size() { return total_extent<Array_t>::value; }
      Data_t* begin() { return reinterpret_cast<Data_t*>(&array); }
      const Data_t* begin() const { return reinterpret_cast<const Data_t*>(&array); }
      Data_t* end() { return begin() + size(); }
      const Data_t* end() const { return begin() + size(); }
      //@}
      
      //@{
      /// Array interface
      auto operator[] (size_t index) -> decltype(*array) { return array[index]; }
      auto operator[] (size_t index) const -> decltype(*array) { return array[index]; }
      auto operator+ (ptrdiff_t index) -> decltype(&*array) { return array + index; }
      auto operator+ (ptrdiff_t index) const -> decltype(&*array) { return array + index; }
      auto operator- (ptrdiff_t index) -> decltype(&*array) { return array - index; }
      auto operator- (ptrdiff_t index) const -> decltype(&*array) { return array - index; }
      auto operator* () -> decltype(*array) { return *array; }
      auto operator* () const -> decltype(*array) { return *array; }
      
      operator decltype(&array[0]) () { return &array[0]; }
      operator decltype(&array[0]) () const { return &array[0]; }
      //@}
      
    }; // BoxedArray
    
    /// Raw Waveform results
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class RawDataStruct {
    public:
      template <typename T>
      using RawData_t = std::vector<T>;
      template <typename T>
      using RawWaveFormData_t = std::vector<BoxedArray<T[kMaxTicks]>>;
      
      size_t MaxRawChannels; ///< maximum number of storable 
      
      Short_t  raw_nchannels;             //number of raw channels
      RawData_t<Short_t>            raw_channelId;
      RawWaveFormData_t<Short_t>    raw_wf; 
      
      /// Creates a RawData data structure allowing up to maxRawChannels Raw channels
      RawDataStruct(size_t maxRawChannels = 0):
        MaxRawChannels(maxRawChannels) { Clear(); }
      
      void Clear();
      void SetMaxRawChannels(size_t maxRawChannels)
      { MaxRawChannels = maxRawChannels; Resize(MaxRawChannels); }
      void Resize(size_t raw_nChannels);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxRawChannels() const { return MaxRawChannels; }
      size_t GetMaxTicksPerRawChannel(int /* irawChannel */ = 0) const
      { return (size_t) kMaxTicks; }    
    }; // class RawDataStruct
    
    
    /// Calibration Waveform results
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class CalibDataStruct {
    public:
      template <typename T>
      using CalibData_t = std::vector<T>;
      template <typename T>
      using CalibWaveFormData_t = std::vector<BoxedArray<T[kMaxTicks]>>;
      
      size_t MaxCalibChannels; ///< maximum number of storable 
      
      Short_t  calib_nchannels;        //number of calibration channels
      CalibData_t<Short_t>            calib_channelId;
      CalibWaveFormData_t<Short_t>    calib_wf; 
      
      /// Creates a CalibData data structure allowing up to maxCalibChannels Calib channels
      CalibDataStruct(size_t maxCalibChannels = 0):
        MaxCalibChannels(maxCalibChannels) { Clear(); }
      
      void Clear();
      void SetMaxCalibChannels(size_t maxCalibChannels)
      { MaxCalibChannels = maxCalibChannels; Resize(MaxCalibChannels); }
      void Resize(size_t calib_nChannels);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxCalibChannels() const { return MaxCalibChannels; }
      size_t GetMaxTicksPerCalibChannel(int /* icalibChannel */ = 0) const
      { return (size_t) kMaxTicks; }    
    }; // class CalibDataStruct
    
    /// Tracker algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class TrackDataStruct {
    public:
      /* Data structure size:
       *
       * TrackData_t<Short_t>                    :  2  bytes/track
       * TrackData_t<Float_t>                    :  4  bytes/track
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/track
       * HitData_t<Float_t>                      : 24k bytes/track
       * HitCoordData_t<Float_t>                 : 72k bytes/track
       */
      template <typename T>
      using TrackData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits][3]>>;
      
      size_t MaxTracks; ///< maximum number of storable tracks
      
      Short_t  ntracks;             //number of reconstructed tracks
      PlaneData_t<Float_t>    trkke;
      PlaneData_t<Float_t>    trkrange;
      PlaneData_t<Int_t>      trkidtruth;  //true geant trackid
      PlaneData_t<Short_t>    trkorigin;   //_ev_origin 0: unknown, 1: neutrino, 2: cosmic, 3: supernova, 4: singles
      PlaneData_t<Int_t>      trkpdgtruth; //true pdg code
      PlaneData_t<Float_t>    trkefftruth; //completeness
      PlaneData_t<Float_t>    trkpurtruth; //purity of track
      PlaneData_t<Float_t>    trkpitchc;
      PlaneData_t<Short_t>    ntrkhits;
      HitData_t<Float_t>      trkdedx;
      HitData_t<Float_t>      trkdqdx;
      HitData_t<Float_t>      trkresrg;
      HitCoordData_t<Float_t> trkxyz;

      // more track info
      TrackData_t<Short_t> trkId;
      TrackData_t<Short_t> trkncosmictags_tagger;
      TrackData_t<Float_t> trkcosmicscore_tagger;
      TrackData_t<Short_t> trkcosmictype_tagger;
      TrackData_t<Short_t> trkncosmictags_containmenttagger;
      TrackData_t<Float_t> trkcosmicscore_containmenttagger;
      TrackData_t<Short_t> trkcosmictype_containmenttagger;
      TrackData_t<Short_t> trkncosmictags_flashmatch;
      TrackData_t<Float_t> trkcosmicscore_flashmatch;
      TrackData_t<Short_t> trkcosmictype_flashmatch;
      TrackData_t<Float_t> trkstartx;     // starting x position.
      TrackData_t<Float_t> trkstarty;     // starting y position.
      TrackData_t<Float_t> trkstartz;     // starting z position.
      TrackData_t<Float_t> trkstartd;     // starting distance to boundary.
      TrackData_t<Float_t> trkendx;       // ending x position.
      TrackData_t<Float_t> trkendy;       // ending y position.
      TrackData_t<Float_t> trkendz;       // ending z position.
      TrackData_t<Float_t> trkendd;       // ending distance to boundary.
      TrackData_t<Float_t> trkACpierceT0;   // t0 per track from anode or cathode piercing tracks (in ns)     
      TrackData_t<Float_t> trkflashT0;   // t0 per track from matching tracks to flashes (in ns)
      TrackData_t<Float_t> trktrueT0;    // t0 per track from truth information (in ns)
      TrackData_t<Float_t> trkpurity;    // track purity based on hit information
      TrackData_t<Float_t> trkcompleteness; //track completeness based on hit information
      TrackData_t<int> trkg4id;        //true g4 track id for the reconstructed track
      TrackData_t<int> trkorig;        //origin of the track 
      TrackData_t<Float_t> trktheta;      // theta.
      TrackData_t<Float_t> trkphi;        // phi.
      TrackData_t<Float_t> trkstartdcosx;
      TrackData_t<Float_t> trkstartdcosy;
      TrackData_t<Float_t> trkstartdcosz;
      TrackData_t<Float_t> trkenddcosx;
      TrackData_t<Float_t> trkenddcosy;
      TrackData_t<Float_t> trkenddcosz;
      TrackData_t<Float_t> trkthetaxz;    // theta_xz.
      TrackData_t<Float_t> trkthetayz;    // theta_yz.
      TrackData_t<Float_t> trkmom;        // momentum.
      TrackData_t<Float_t> trklen;        // length.
      TrackData_t<Float_t> trkmomrange;    // track momentum from range using CSDA tables
      TrackData_t<Float_t> trkmommschi2;   // track momentum from multiple scattering Chi2 method
      TrackData_t<Float_t> trkmommsllhd;   // track momentum from multiple scattering LLHD method
      TrackData_t<Float_t> trkmcsfwdmom;   // track momentum assuming forward direction
      TrackData_t<Float_t> trkmcsfwdll;    // MCS likelihood assuming forward direction
      TrackData_t<Float_t> trkmcsfwderr;   // MCS uncertainty assuming forward direction
      TrackData_t<Float_t> trkmcsbwdmom;   // track momentum assuming backward direction
      TrackData_t<Float_t> trkmcsbwdll;    // MCS likelihood assuming backward direction
      TrackData_t<Float_t> trkmcsbwderr;   // MCS uncertainty assuming backward direction
      TrackData_t<Short_t> trksvtxid;     // Vertex ID associated with the track start
      TrackData_t<Short_t> trkevtxid;     // Vertex ID associated with the track end
      TrackData_t<int>     trkntraj; // Number of trajectory points on the track (there is a one-to-one correspondence between trajectory points and hits, so this should also tell you the number of hits)
      PlaneData_t<Int_t> trkpidpdg;       // particle PID pdg code
      PlaneData_t<Float_t> trkpidchi;
      PlaneData_t<Float_t> trkpidchipr;   // particle PID chisq for proton
      PlaneData_t<Float_t> trkpidchika;   // particle PID chisq for kaon
      PlaneData_t<Float_t> trkpidchipi;   // particle PID chisq for pion
      PlaneData_t<Float_t> trkpidchimu;   // particle PID chisq for muon
      PlaneData_t<Float_t> trkpidpida;    // particle PIDA
      TrackData_t<Short_t> trkpidbestplane; // this is defined as the plane with most hits
	
	  TrackData_t<Short_t> trkhasPFParticle; // whether this belongs to a PFParticle 
	  TrackData_t<Short_t> trkPFParticleID;  // if hasPFParticle, its ID
       
      /// Creates an empty tracker data structure
      TrackDataStruct(): MaxTracks(0) { Clear(); }
      /// Creates a tracker data structure allowing up to maxTracks tracks
      TrackDataStruct(size_t maxTracks): MaxTracks(maxTracks) { Clear(); }
      void Clear();
      void SetMaxTracks(size_t maxTracks)
      { MaxTracks = maxTracks; Resize(MaxTracks); }
      void Resize(size_t nTracks);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);
      
      size_t GetMaxTracks() const { return MaxTracks; }
      size_t GetMaxPlanesPerTrack(int /* iTrack */ = 0) const
      { return (size_t) kNplanes; }
      size_t GetMaxHitsPerTrack(int /* iTrack */ = 0, int /* ipl */ = 0) const
      { return (size_t) kMaxTrackHits; }
      
    }; // class TrackDataStruct

    //Vertex data struct
    class VertexDataStruct {
    public:
      template <typename T>
      using VertexData_t = std::vector<T>;

      size_t MaxVertices; ///< maximum number of storable vertices

      Short_t  nvtx;             //number of reconstructed vertices
      VertexData_t<Short_t> vtxId;    // the vertex ID.
      VertexData_t<Float_t> vtxx;     // x position.
      VertexData_t<Float_t> vtxy;     // y position.
      VertexData_t<Float_t> vtxz;     // z position.
	  
	  VertexData_t<Short_t> vtxhasPFParticle; // whether this belongs to a PFParticle 
	  VertexData_t<Short_t> vtxPFParticleID;  // if hasPFParticle, its ID

      VertexDataStruct(): MaxVertices(0) { Clear(); }
      VertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
      void Clear();
      void SetMaxVertices(size_t maxVertices)
      { MaxVertices = maxVertices; Resize(MaxVertices); }
      void Resize(size_t nVertices);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxVertices() const { return MaxVertices; }
    }; // class VertexDataStruct 

    //Neutrino vertex data struct
    class NeutrinoVertexDataStruct {
    public:
      template <typename T>
      using VertexData_t = std::vector<T>;

      size_t MaxVertices;             ///< maximum number of storable vertices

      Short_t  nvtx;                  ///< number of neutrino reconstructed vertices
      VertexData_t<Short_t> vtxId;    ///< the vertex ID.
      VertexData_t<Float_t> vtxx;     ///< x position.
      VertexData_t<Float_t> vtxy;     ///< y position.
      VertexData_t<Float_t> vtxz;     ///< z position.
      VertexData_t<Int_t>   vtxpdg;   ///< pdg of pfp

      VertexData_t<Short_t> vtxhasPFParticle; // whether this belongs to a PFParticle 
      VertexData_t<Short_t> vtxPFParticleID;  // if hasPFParticle, its ID

      NeutrinoVertexDataStruct(): MaxVertices(0) { Clear(); }
      NeutrinoVertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
      void Clear();
      void SetMaxVertices(size_t maxVertices)
      { MaxVertices = maxVertices; Resize(MaxVertices); }
      void Resize(size_t nVertices);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxVertices() const { return MaxVertices; }
    }; // class NeutrinoVertexDataStruct 

    // Flash data struct
    class FlashDataStruct {
    public:
      template <typename T>
      using FlashData_t     = std::vector<T>;
      template <typename T>
      using FlashDataSpec_t = std::vector<BoxedArray<T[kNOpDets]>>;

      size_t MaxFlashes;                                ///< Maximum number of storable flashes

      Short_t nfls;                                     ///< Number of reconstructed flashes
      FlashData_t<Float_t> flsTime;                     ///< Flash time (us)
      FlashData_t<Float_t> flsPe;                       ///< Flash total PE
      FlashDataSpec_t<Double_t> flsPePerOpDet;          ///< Flash PE per optical detector
      FlashData_t<Float_t> flsXcenter;                  ///< Flash X center (cm)
      FlashData_t<Float_t> flsYcenter;                  ///< Flash Y center (cm)
      FlashData_t<Float_t> flsZcenter;                  ///< Flash Z center (cm)
      FlashData_t<Float_t> flsYwidth;                   ///< Flash Y width (cm)
      FlashData_t<Float_t> flsZwidth;                   ///< Flash Z width (cm)
      FlashData_t<Float_t> flsTwidth;                   ///< Flash time width (us)

      void Clear();
      void SetMaxFlashes(size_t maxFlashes)
      { MaxFlashes = maxFlashes; Resize(MaxFlashes); }
      void Resize(size_t nFlashes);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxFlashes() const { return MaxFlashes; }

      size_t GetNOpDet() const { auto const* geom = lar::providerFrom<geo::Geometry>(); return geom->NOpDets(); }
    }; // class FlashDataStruct

    /// Shower algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class ShowerDataStruct {
    public:
      /* Data structure size:
       *
       * ShowerData_t<Short_t>                   :  2  bytes/shower
       * ShowerData_t<Float_t>                   :  4  bytes/shower
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/shower
       * HitData_t<Float_t>                      : 24k bytes/shower
       * HitCoordData_t<Float_t>                 : 72k bytes/shower
       */
      template <typename T>
      using ShowerData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits][3]>>;
      
      std::string name; ///< name of the shower algorithm (for branch names)
      
      size_t MaxShowers; ///< maximum number of storable showers
      
      /// @{
      /// @name Branch data structures
      Short_t  nshowers;                      ///< number of showers
      ShowerData_t<Short_t>  showerID;        ///< Shower ID
      ShowerData_t<Short_t>  shwr_bestplane;  ///< Shower best plane
      ShowerData_t<Float_t>  shwr_length;     ///< Shower length
      ShowerData_t<Float_t>  shwr_theta;      ///< Shower theta
      ShowerData_t<Float_t>  shwr_thetaxz;    ///< Shower theta xz projection
      ShowerData_t<Float_t>  shwr_thetayz;    ///< Shower theta yz projection
      ShowerData_t<Float_t>  shwr_phi;        ///< Shower phi
      ShowerData_t<Float_t>  shwr_startdcosx; ///< X directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startdcosy; ///< Y directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startdcosz; ///< Z directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startx;     ///< startx of shower
      ShowerData_t<Float_t>  shwr_starty;     ///< starty of shower
      ShowerData_t<Float_t>  shwr_startz;     ///< startz of shower
      PlaneData_t<Float_t>   shwr_totEng;     ///< Total energy of the shower per plane
      PlaneData_t<Float_t>   shwr_dedx;       ///< dE/dx of the shower per plane
      PlaneData_t<Float_t>   shwr_mipEng;     ///< Total MIP energy of the shower per plane
	  
	  ShowerData_t<Short_t>  shwr_hasPFParticle; // whether this belongs to a PFParticle 
	  ShowerData_t<Short_t>  shwr_PFParticleID;  // if hasPFParticle, its ID
      /// @}
      
      /// Creates a shower data structure allowing up to maxShowers showers
      ShowerDataStruct(std::string new_name = "", size_t maxShowers = 0):
        name(new_name), MaxShowers(maxShowers) { Clear(); }
      
      std::string Name() const { return name; }
      
      void Clear();
      
      /// Applies a special prescription to mark shower information as missing
      void MarkMissing(TTree* pTree);
      void SetName(std::string new_name) { name = new_name; }
      void SetMaxShowers(size_t maxShowers)
      { MaxShowers = maxShowers; Resize(MaxShowers); }
      void Resize(size_t nShowers);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxShowers() const { return MaxShowers; }
      size_t GetMaxPlanesPerShower(int /* iShower */ = 0) const
      { return (size_t) kNplanes; }
      size_t GetMaxHitsPerShower(int /* iShower */ = 0, int /* ipl */ = 0) const
      { return (size_t) kMaxShowerHits; }
      
    }; // class ShowerDataStruct
    
    class PFParticleDataStruct {
    public:
      /* Data structure size:
       *
       * PFParticleData_t<Short_t>   :  2  bytes/PFParticle
       * PFParticleData_t<Int_t>     :  4  bytes/PFParticle
       * DaughterData_t<Short_t>     :  20 bytes/PFParticle
       * ClusterData_t<Short_t>      :  20 bytes/PFParticle
       * Short_t [kMaxNPFPNeutrinos] :  10 bytes in total
       */
      template <typename T>
      using PFParticleData_t = std::vector<T>;
      template <typename T>
      using DaughterData_t = std::vector<BoxedArray<T[kMaxNDaughtersPerPFP]>>;
      template <typename T>
      using ClusterData_t = std::vector<BoxedArray<T[kMaxNClustersPerPFP]>>;
      
      size_t MaxPFParticles; ///< maximum number of storable PFParticles
      
      /// @{
      /// @name Branch data structures
      Short_t                   nPFParticles;     ///< the total number of PFParticles
      PFParticleData_t<Short_t> pfp_selfID;       ///< the PFParticles' own IDs
      PFParticleData_t<Short_t> pfp_isPrimary;    ///< whether the PFParticle is a primary particle
      
      PFParticleData_t<Short_t> pfp_numDaughters; ///< the number of daughters belonging to this PFParticle
      DaughterData_t<Short_t>   pfp_daughterIDs;  ///< the IDs of the daughter PFParticles
      PFParticleData_t<Short_t> pfp_parentID;     ///< the ID of this PFParticle's immediate parent
      
      PFParticleData_t<Short_t> pfp_vertexID;     ///< the ID of the vertex belonging to this PFParticle
      PFParticleData_t<Short_t> pfp_isShower;     ///< whether this PFParticle corresponds to a shower
      PFParticleData_t<Short_t> pfp_isTrack;      ///< whether this PFParticle corresponds to a track
      PFParticleData_t<Short_t> pfp_trackID;      ///< the ID of the track object corresponding to this PFParticle, if !isShower
      PFParticleData_t<Short_t> pfp_showerID;     ///< the ID of the shower object corresponding to this PFParticle, if isShower
      
      PFParticleData_t<Short_t> pfp_isNeutrino;   ///< whether this PFParticle is a neutrino
      PFParticleData_t<Int_t>   pfp_pdgCode;      ///< the preliminary estimate of the PFParticle type using the PDG code
      
      PFParticleData_t<Short_t> pfp_numClusters;  ///< the number of associated clusters
      ClusterData_t<Short_t>    pfp_clusterIDs;   ///< the IDs of any associated clusters
      
      Short_t                   pfp_numNeutrinos; ///< the number of reconstructed neutrinos
      Short_t pfp_neutrinoIDs[kMaxNPFPNeutrinos]; ///< the PFParticle IDs of the neutrinos
      /// @}
      
      /// Creates a PFParticle data structure allowing up to maxPFParticles PFParticles
      PFParticleDataStruct(size_t maxPFParticles = 0):
        MaxPFParticles(maxPFParticles) { Clear(); }
      
      void Clear();
      void SetMaxPFParticles(size_t maxPFParticles)
        { MaxPFParticles = maxPFParticles; Resize(MaxPFParticles); }
      void Resize(size_t numPFParticles);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxPFParticles() const { return MaxPFParticles; }
      size_t GetMaxDaughtersPerPFParticle(int /* iPFParticle */ = 0) const
        { return (size_t) kMaxNDaughtersPerPFP; }
      size_t GetMaxClustersPerPFParticle(int /* iPFParticle */ = 0) const
        { return (size_t) kMaxNClustersPerPFP; }
      
    }; // class PFParticleDataStruct

    enum DataBits_t: unsigned int {
      tdAuxDet = 0x01,
	tdCry = 0x02,
	tdGenie = 0x04,
	tdGeant = 0x08,
	tdHit = 0x10,
	tdTrack = 0x20,
	tdVertex = 0x40,
	tdFlash = 0x80,
	tdShower = 0x100,
	tdMCshwr = 0x200,
	tdMCtrk  = 0x400,
	tdCluster = 0x800,
	tdRawDigit = 0x1000,
	tdCalWire  = 0x2000,
	tdSimChannel = 0x4000,
	tdRawWaveForm = 0x8000,
	tdCalibWaveForm = 0x10000,
        tdPandoraNuVertex = 0x20000,
        tdPFParticle = 0x40000,
        tdSWTrigger = 0x80000,
        tdOpticalFilter = 0x100000,
	tdDefault = 0
	}; // DataBits_t
    
    /*    /// information from the run
	  struct RunData_t {
	  public:
	  RunData_t() { Clear(); }
	  void Clear() {}
	  }; // struct RunData_t
    */
    /// information from the subrun
    struct SubRunData_t {
      SubRunData_t() { Clear(); }
      void Clear() { 
        pot = -99999.; 
        potbnbETOR860 = -99999.;
        potbnbETOR875 = -99999.;
        potnumiETORTGT = -99999.;
        run = -99999;
        subrun = -99999;
        begintime = -99999;
        endtime = -99999;
      }
      Double_t pot; //protons on target
      Double_t potbnbETOR860;
      Double_t potbnbETOR875;
      Double_t potnumiETORTGT;
      Int_t    run;
      Int_t    subrun;
      Double_t begintime;
      Double_t endtime;
    }; // struct SubRunData_t

    //    RunData_t    RunData; ///< run data collected at begin of run
    SubRunData_t SubRunData; ///< subrun data collected at begin of subrun

    //run information
    Int_t      run;                  //run number
    Int_t      subrun;               //subrun number
    Int_t      event;                //event number
    Double_t   evttime;              //event time in sec
    Double_t   beamtime;             //beam time
    //  Double_t   pot;                  //protons on target moved in subrun data
    Double_t   taulife;              //electron lifetime
    Char_t     isdata;               //flag, 0=MC 1=data
    unsigned int triggernumber;      //trigger counter
    Double_t     triggertime;        //trigger time w.r.t. electronics clock T0
    Double_t     beamgatetime;       //beamgate time w.r.t. electronics clock T0
    unsigned int triggerbits;        //trigger bits
    Double_t     potbnb;             //pot per event (BNB E:TOR860)
    Double_t     potnumitgt;         //pot per event (NuMI E:TORTGT)
    Double_t     potnumi101;         //pot per event (NuMI E:TOR101)

    std::vector<std::string> evtwgt_funcname;          // the name of the functions used
    std::vector<std::vector<double>> evtwgt_weight;    // the weights (a vector for each function used)
    std::vector<int> evtwgt_nweight;                   // number of weights for each function
    Int_t evtwgt_nfunc;                                // number of functions used

    //Optical Filter Information
    Float_t opfilter_pebeam;      // pe in the beam window period
    Float_t opfilter_peveto;      // pe in the veto window period
    Float_t opfilter_pmtmaxfrac;  // max fraction of pe in single pmt

    // Software trigger information
    std::vector<std::string> swtrigger_name;           // name of the software trigger used 
    std::vector<bool>        swtrigger_triggered;      // if this event was triggerd or not (based on the relative swtrigger_name logic)

    // hit information (non-resizeable, 45x kMaxHits = 900k bytes worth)
    Int_t    no_hits;                  //number of hits
    Int_t    no_hits_stored;                  //number of hits actually stored in the tree    
    Short_t  hit_plane[kMaxHits];      //plane number
    Short_t  hit_wire[kMaxHits];       //wire number
    Short_t  hit_channel[kMaxHits];    //channel ID
    Float_t  hit_peakT[kMaxHits];      //peak time
    Float_t  hit_charge[kMaxHits];     //charge (area)
    Float_t  hit_ph[kMaxHits];         //amplitude
    Float_t  hit_startT[kMaxHits];     //hit start time
    Float_t  hit_endT[kMaxHits];       //hit end time
    Float_t  hit_rms[kMaxHits];       //hit rms from the hit object
    Float_t  hit_goodnessOfFit[kMaxHits]; //chi2/dof goodness of fit 
    Short_t  hit_multiplicity[kMaxHits];  //multiplicity of the given hit					 
    Float_t  hit_trueX[kMaxHits];      // hit true X (cm)
    Float_t  hit_nelec[kMaxHits];     //hit number of electrons
    Float_t  hit_energy[kMaxHits];       //hit energy
    Short_t  hit_trkid[kMaxHits];      //is this hit associated with a reco track?
    Short_t  hit_trkKey[kMaxHits];      //is this hit associated with a reco track,  if so associate a unique track key ID?
    Short_t  hit_clusterid[kMaxHits];  //is this hit associated with a reco cluster?
    Short_t  hit_clusterKey[kMaxHits];  //is this hit associated with a reco cluster, if so associate a unique cluster key ID?

    //RawDigit information
    Float_t  rawD_ph[kMaxHits];  
    Int_t    rawD_peakT[kMaxHits];  
    Float_t  rawD_charge[kMaxHits];  
    Float_t  rawD_fwhh[kMaxHits];  
    Double_t rawD_rms[kMaxHits]; 
    
    //CalWire information
    Float_t  calwire_ph[kMaxHits];  
    Int_t    calwire_peakT[kMaxHits];  
    Float_t  calwire_charge[kMaxHits];  
    Float_t  calwire_fwhh[kMaxHits];  
    Double_t calwire_rms[kMaxHits]; 
    
    //SimChannel information
    Float_t  sim_ph[kMaxHits];  
    Int_t    sim_tdc[kMaxHits];  
    Float_t  sim_charge[kMaxHits];  
    Float_t  sim_fwhh[kMaxHits];  
    Double_t sim_rms[kMaxHits]; 

    //Cluster Information
    Short_t nclusters;				      //number of clusters in a given event
    Short_t clusterId[kMaxClusters];		      //ID of this cluster	 
    Short_t clusterView[kMaxClusters];		      //which plane this cluster belongs to	      
    Int_t   cluster_isValid[kMaxClusters];	      //is this cluster valid? will have a value of -1 if it is not valid
    Float_t cluster_StartCharge[kMaxClusters];	       //charge on the first wire of the cluster in ADC
    Float_t cluster_StartAngle[kMaxClusters];	      //starting angle of the cluster
    Float_t cluster_EndCharge[kMaxClusters];	      //charge on the last wire of the cluster in ADC
    Float_t cluster_EndAngle[kMaxClusters];	      //ending angle of the cluster
    Float_t cluster_Integral[kMaxClusters];	      //returns the total charge of the cluster from hit shape in ADC
    Float_t cluster_IntegralAverage[kMaxClusters];    //average charge of the cluster hits in ADC
    Float_t cluster_SummedADC[kMaxClusters];	      //total charge of the cluster from signal ADC counts
    Float_t cluster_SummedADCaverage[kMaxClusters];   //average signal ADC counts of the cluster hits.
    Float_t cluster_MultipleHitDensity[kMaxClusters]; //Density of wires in the cluster with more than one hit. 
    Float_t cluster_Width[kMaxClusters];	      //cluster width in ? units
    Short_t cluster_NHits[kMaxClusters];	      //Number of hits in the cluster
    Short_t cluster_StartWire[kMaxClusters];	      //wire coordinate of the start of the cluster 
    Short_t cluster_StartTick[kMaxClusters];	      //tick coordinate of the start of the cluster in time ticks
    Short_t cluster_EndWire[kMaxClusters];	      //wire coordinate of the end of the cluster
    Short_t cluster_EndTick[kMaxClusters];            //tick coordinate of the end of the cluster in time ticks    
    //Cluster cosmic tagging information
    Short_t cluncosmictags_tagger[kMaxClusters];      //No. of cosmic tags associated to this cluster    
    Float_t clucosmicscore_tagger[kMaxClusters];      //Cosmic score associated to this cluster. In the case of more than one tag, the first one is associated.    
    Short_t clucosmictype_tagger[kMaxClusters];       //Cosmic tag type for this cluster.    
    
    // Flash information
    Char_t kNFlashAlgos;
    std::vector<FlashDataStruct> FlashData;

    //track information
    Char_t   kNTracker;
    std::vector<TrackDataStruct> TrackData;

    //vertex information
    Char_t   kNVertexAlgos;
    std::vector<VertexDataStruct> VertexData;
   
    //neutrino vertex information
    Char_t   kNNeutrinoVertexAlgos;
    std::vector<NeutrinoVertexDataStruct> NeutrinoVertexData;

    // shower information
    Char_t   kNShowerAlgos;
    std::vector<ShowerDataStruct> ShowerData;

    // PFParticle information
    PFParticleDataStruct PFParticleData;
    
    // Raw Waveform information
    RawDataStruct RawData;
    
    // Calibration Waveform information
    CalibDataStruct CalibData;

    //mctruth information
    Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
    Int_t     nuPDG_truth[kMaxTruth];     //neutrino PDG code
    Int_t     ccnc_truth[kMaxTruth];      //0=CC 1=NC
    Int_t     mode_truth[kMaxTruth];      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    Float_t  enu_truth[kMaxTruth];       //true neutrino energy
    Float_t  Q2_truth[kMaxTruth];        //Momentum transfer squared
    Float_t  W_truth[kMaxTruth];         //hadronic invariant mass
    Float_t  X_truth[kMaxTruth];
    Float_t  Y_truth[kMaxTruth];
    Int_t     hitnuc_truth[kMaxTruth];    //hit nucleon
    Float_t  nuvtxx_truth[kMaxTruth];    //neutrino vertex x
    Float_t  nuvtxy_truth[kMaxTruth];    //neutrino vertex y
    Float_t  nuvtxz_truth[kMaxTruth];    //neutrino vertex z
    Float_t  sp_charge_corrected_nuvtxx_truth[kMaxTruth]; 
    Float_t  sp_charge_corrected_nuvtxy_truth[kMaxTruth]; 
    Float_t  sp_charge_corrected_nuvtxz_truth[kMaxTruth]; 
    Float_t  nu_dcosx_truth[kMaxTruth];  //neutrino dcos x
    Float_t  nu_dcosy_truth[kMaxTruth];  //neutrino dcos y
    Float_t  nu_dcosz_truth[kMaxTruth];  //neutrino dcos z
    Float_t  lep_mom_truth[kMaxTruth];   //lepton momentum
    Float_t  lep_dcosx_truth[kMaxTruth]; //lepton dcos x
    Float_t  lep_dcosy_truth[kMaxTruth]; //lepton dcos y
    Float_t  lep_dcosz_truth[kMaxTruth]; //lepton dcos z

    //flux information
    Float_t  vx_flux[kMaxTruth];          //X position of hadron/muon decay (cm)
    Float_t  vy_flux[kMaxTruth];          //Y position of hadron/muon decay (cm)
    Float_t  vz_flux[kMaxTruth];          //Z position of hadron/muon decay (cm)
    Float_t  pdpx_flux[kMaxTruth];        //Parent X momentum at decay point (GeV)
    Float_t  pdpy_flux[kMaxTruth];        //Parent Y momentum at decay point (GeV)
    Float_t  pdpz_flux[kMaxTruth];        //Parent Z momentum at decay point (GeV)
    Float_t  ppdxdz_flux[kMaxTruth];      //Parent dxdz direction at production
    Float_t  ppdydz_flux[kMaxTruth];      //Parent dydz direction at production
    Float_t  pppz_flux[kMaxTruth];        //Parent Z momentum at production (GeV)
    
    Int_t    ptype_flux[kMaxTruth];        //Parent GEANT code particle ID
    Float_t  ppvx_flux[kMaxTruth];        //Parent production vertex X (cm)
    Float_t  ppvy_flux[kMaxTruth];        //Parent production vertex Y (cm)
    Float_t  ppvz_flux[kMaxTruth];        //Parent production vertex Z (cm)
    Float_t  muparpx_flux[kMaxTruth];     //Muon neutrino parent production vertex X (cm)
    Float_t  muparpy_flux[kMaxTruth];     //Muon neutrino parent production vertex Y (cm)
    Float_t  muparpz_flux[kMaxTruth];     //Muon neutrino parent production vertex Z (cm)
    Float_t  mupare_flux[kMaxTruth];      //Muon neutrino parent energy (GeV)
    
    Int_t    tgen_flux[kMaxTruth];        //Parent generation in cascade. (1 = primary proton, 
    					  //2=particles produced by proton interaction, 3 = particles 
					  //produced by interactions of the 2's, ...
    Int_t    tgptype_flux[kMaxTruth];     //Type of particle that created a particle flying of the target  
    Float_t  tgppx_flux[kMaxTruth];       //X Momentum of a particle, that created a particle that flies 
    					  //off the target, at the interaction point. (GeV)
    Float_t  tgppy_flux[kMaxTruth];       //Y Momentum of a particle, that created a particle that flies  
					  //off the target, at the interaction point. (GeV)
    Float_t  tgppz_flux[kMaxTruth];       //Z Momentum of a particle, that created a particle that flies
					  //off the target, at the interaction point. (GeV)
    Float_t  tprivx_flux[kMaxTruth];      //Primary particle interaction vertex X (cm)  
    Float_t  tprivy_flux[kMaxTruth];      //Primary particle interaction vertex Y (cm)  
    Float_t  tprivz_flux[kMaxTruth];      //Primary particle interaction vertex Z (cm) 
    Float_t  dk2gen_flux[kMaxTruth];      //distance from decay to ray origin (cm) 
    Float_t  gen2vtx_flux[kMaxTruth];     //distance from ray origin to event vtx (cm)
 
    Float_t  tpx_flux[kMaxTruth];        //Px of parent particle leaving BNB/NuMI target (GeV)
    Float_t  tpy_flux[kMaxTruth];        //Py of parent particle leaving BNB/NuMI target (GeV)
    Float_t  tpz_flux[kMaxTruth];        //Pz of parent particle leaving BNB/NuMI target (GeV)
    Int_t    tptype_flux[kMaxTruth];     //Type of parent particle leaving BNB target
     
    //genie information
    size_t MaxGeniePrimaries = 0;
    Int_t     genie_no_primaries;
    std::vector<Int_t>    genie_primaries_pdg;
    std::vector<Float_t>  genie_Eng;
    std::vector<Float_t>  genie_Px;
    std::vector<Float_t>  genie_Py;
    std::vector<Float_t>  genie_Pz;
    std::vector<Float_t>  genie_P;
    std::vector<Int_t>    genie_status_code;
    std::vector<Float_t>  genie_mass;
    std::vector<Int_t>    genie_trackID;
    std::vector<Int_t>    genie_ND;
    std::vector<Int_t>    genie_mother;
    
    //cosmic cry information
    Int_t     mcevts_truthcry;    //number of neutrino Int_teractions in the spill
    Int_t     cry_no_primaries;
    std::vector<Int_t>    cry_primaries_pdg;
    std::vector<Float_t>  cry_Eng;
    std::vector<Float_t>  cry_Px;
    std::vector<Float_t>  cry_Py;
    std::vector<Float_t>  cry_Pz;
    std::vector<Float_t>  cry_P;
    std::vector<Float_t>  cry_StartPointx;
    std::vector<Float_t>  cry_StartPointy;
    std::vector<Float_t>  cry_StartPointz;
    std::vector<Float_t>  cry_StartPointt;
    std::vector<Int_t>    cry_status_code;
    std::vector<Float_t>  cry_mass;
    std::vector<Int_t>    cry_trackID;
    std::vector<Int_t>    cry_ND;
    std::vector<Int_t>    cry_mother;
    
    //G4 MC Particle information
    size_t MaxGEANTparticles = 0; ///! how many particles there is currently room for
    Int_t     no_primaries;      //number of primary geant particles
    Int_t     geant_list_size;  //number of all geant particles
    Int_t     geant_list_size_in_tpcAV;
    std::vector<Int_t>    pdg;
    std::vector<Int_t>    status;    
    std::vector<Float_t>  Eng;
    std::vector<Float_t>  EndE;
    std::vector<Float_t>  Mass;
    std::vector<Float_t>  Px;
    std::vector<Float_t>  Py;
    std::vector<Float_t>  Pz;
    std::vector<Float_t>  P;
    std::vector<Float_t>  StartPointx;
    std::vector<Float_t>  sp_charge_corrected_StartPointx; 
    std::vector<Float_t>  sp_charge_corrected_StartPointy; 
    std::vector<Float_t>  sp_charge_corrected_StartPointz; 
    std::vector<Float_t>  StartPointy;
    std::vector<Float_t>  StartPointz;
    std::vector<Float_t>  StartT;  
    std::vector<Float_t>  EndT;          
    std::vector<Float_t>  EndPointx;
    std::vector<Float_t>  EndPointy;
    std::vector<Float_t>  EndPointz;
    std::vector<Float_t>  sp_charge_corrected_EndPointx; 
    std::vector<Float_t>  sp_charge_corrected_EndPointy;
    std::vector<Float_t>  sp_charge_corrected_EndPointz; 
    std::vector<Float_t>  theta;    
    std::vector<Float_t>  phi;    
    std::vector<Float_t>  theta_xz;    
    std::vector<Float_t>  theta_yz;    
    std::vector<Float_t>  pathlen;    
    std::vector<Int_t>    inTPCActive;    
    std::vector<Float_t>  StartPointx_tpcAV;
    std::vector<Float_t>  StartPointy_tpcAV;
    std::vector<Float_t>  StartPointz_tpcAV;
    std::vector<Float_t>  sp_charge_corrected_StartPointx_tpcAV; 
    std::vector<Float_t>  sp_charge_corrected_StartPointy_tpcAV; 
    std::vector<Float_t>  sp_charge_corrected_StartPointz_tpcAV; 
    std::vector<Float_t>  StartT_tpcAV;
    std::vector<Float_t>  StartE_tpcAV;
    std::vector<Float_t>  StartP_tpcAV;
    std::vector<Float_t>  StartPx_tpcAV;
    std::vector<Float_t>  StartPy_tpcAV;
    std::vector<Float_t>  StartPz_tpcAV;
    std::vector<Float_t>  EndPointx_tpcAV;
    std::vector<Float_t>  EndPointy_tpcAV;
    std::vector<Float_t>  EndPointz_tpcAV;
    std::vector<Float_t>  sp_charge_corrected_EndPointx_tpcAV; 
    std::vector<Float_t>  sp_charge_corrected_EndPointy_tpcAV; 
    std::vector<Float_t>  sp_charge_corrected_EndPointz_tpcAV;
    std::vector<Float_t>  EndT_tpcAV;
    std::vector<Float_t>  EndE_tpcAV;
    std::vector<Float_t>  EndP_tpcAV;
    std::vector<Float_t>  EndPx_tpcAV;
    std::vector<Float_t>  EndPy_tpcAV;
    std::vector<Float_t>  EndPz_tpcAV;
    std::vector<Float_t>  pathlen_drifted;    
    std::vector<Int_t>    inTPCDrifted;    
    std::vector<Float_t>  StartPointx_drifted;
    std::vector<Float_t>  StartPointy_drifted;
    std::vector<Float_t>  StartPointz_drifted;
    std::vector<Float_t>  sp_charge_corrected_StartPointx_drifted; 
    std::vector<Float_t>  sp_charge_corrected_StartPointy_drifted;
    std::vector<Float_t>  sp_charge_corrected_StartPointz_drifted;
    std::vector<Float_t>  StartT_drifted;
    std::vector<Float_t>  StartE_drifted;
    std::vector<Float_t>  StartP_drifted;
    std::vector<Float_t>  StartPx_drifted;
    std::vector<Float_t>  StartPy_drifted;
    std::vector<Float_t>  StartPz_drifted;
    std::vector<Float_t>  EndPointx_drifted;
    std::vector<Float_t>  EndPointy_drifted;
    std::vector<Float_t>  EndPointz_drifted;
    std::vector<Float_t>  sp_charge_corrected_EndPointx_drifted;
    std::vector<Float_t>  sp_charge_corrected_EndPointy_drifted;
    std::vector<Float_t>  sp_charge_corrected_EndPointz_drifted;
    std::vector<Float_t>  EndT_drifted;
    std::vector<Float_t>  EndE_drifted;
    std::vector<Float_t>  EndP_drifted;
    std::vector<Float_t>  EndPx_drifted;
    std::vector<Float_t>  EndPy_drifted;
    std::vector<Float_t>  EndPz_drifted;
    std::vector<Int_t>    NumberDaughters;
    std::vector<Int_t>    TrackId;
    std::vector<Int_t>    Mother;
    std::vector<Int_t>    process_primary;
    std::vector<std::string> processname;
    std::vector<Int_t>    MergedId; //geant track segments, which belong to the same particle, get the same
    std::vector<Int_t>    origin;   ////0: unknown, 1: neutrino, 2: cosmic, 3: supernova, 4: singles 
    std::vector<Int_t>    MCTruthIndex; //this geant particle comes from the neutrino interaction of the _truth variables with this index

    //MC Shower information
    Int_t     no_mcshowers;                         //number of MC Showers in this event.
    //MC Shower particle information
    std::vector<Int_t>       mcshwr_origin;	    //MC Shower origin information.  
    std::vector<Int_t>       mcshwr_pdg;	    //MC Shower particle PDG code.   
    std::vector<Int_t>       mcshwr_TrackId;        //MC Shower particle G4 track ID.
    std::vector<std::string> mcshwr_Process;	    //MC Shower particle's creation process. 
    std::vector<Float_t>     mcshwr_startX;	    //MC Shower particle G4 startX 
    std::vector<Float_t>     mcshwr_startY;	    //MC Shower particle G4 startY 
    std::vector<Float_t>     mcshwr_startZ;	    //MC Shower particle G4 startZ 
    std::vector<Float_t>     mcshwr_endX;	    //MC Shower particle G4 endX 
    std::vector<Float_t>     mcshwr_endY;	    //MC Shower particle G4 endY 
    std::vector<Float_t>     mcshwr_endZ;	    //MC Shower particle G4 endZ 
    std::vector<Float_t>    mcshwr_CombEngX;	    //MC Shower Combined energy deposition information, Start Point X Position. 
    std::vector<Float_t>    mcshwr_CombEngY;	    //MC Shower Combined energy deposition information, Start Point Y Position.
    std::vector<Float_t>    mcshwr_CombEngZ;	    //MC Shower Combined energy deposition information, Start Point Z Position.
    std::vector<Float_t>     mcshwr_CombEngPx;	    //MC Shower Combined energy deposition information, Momentum X direction.
    std::vector<Float_t>     mcshwr_CombEngPy;	    //MC Shower Combined energy deposition information, Momentum X direction.
    std::vector<Float_t>     mcshwr_CombEngPz;	    //MC Shower Combined energy deposition information, Momentum X direction.
    std::vector<Float_t>     mcshwr_CombEngE;	    //MC Shower Combined energy deposition information, Energy
    std::vector<Float_t>     mcshwr_dEdx;           //MC Shower dEdx, MeV/cm
    std::vector<Float_t>     mcshwr_StartDirX;      //MC Shower Direction of begining of shower, X direction 
    std::vector<Float_t>     mcshwr_StartDirY;      //MC Shower Direction of begining of shower, Y direction 
    std::vector<Float_t>     mcshwr_StartDirZ;      //MC Shower Direction of begining of shower, Z direction 
    std::vector<Int_t>       mcshwr_isEngDeposited;  //tells whether if this shower deposited energy in the detector or not.
    //yes = 1; no =0;	
    //MC Shower mother information
    std::vector<Int_t>       mcshwr_Motherpdg;       //MC Shower's mother PDG code. 
    std::vector<Int_t>       mcshwr_MotherTrkId;     //MC Shower's mother G4 track ID.
    std::vector<std::string> mcshwr_MotherProcess;   //MC Shower's mother creation process. 
    std::vector<Float_t>     mcshwr_MotherstartX;    //MC Shower's mother  G4 startX .
    std::vector<Float_t>     mcshwr_MotherstartY;    //MC Shower's mother  G4 startY .
    std::vector<Float_t>     mcshwr_MotherstartZ;    //MC Shower's mother  G4 startZ .
    std::vector<Float_t>     mcshwr_MotherendX;	     //MC Shower's mother  G4 endX   .
    std::vector<Float_t>     mcshwr_MotherendY;	     //MC Shower's mother  G4 endY   .
    std::vector<Float_t>     mcshwr_MotherendZ;	     //MC Shower's mother  G4 endZ   .
    //MC Shower ancestor information
    std::vector<Int_t>       mcshwr_Ancestorpdg;       //MC Shower's ancestor PDG code. 
    std::vector<Int_t>       mcshwr_AncestorTrkId;     //MC Shower's ancestor G4 track ID.
    std::vector<std::string> mcshwr_AncestorProcess;   //MC Shower's ancestor creation process. 
    std::vector<Float_t>     mcshwr_AncestorstartX;    //MC Shower's ancestor  G4 startX
    std::vector<Float_t>     mcshwr_AncestorstartY;    //MC Shower's ancestor  G4 startY
    std::vector<Float_t>     mcshwr_AncestorstartZ;    //MC Shower's ancestor  G4 startZ
    std::vector<Float_t>     mcshwr_AncestorendX;      //MC Shower's ancestor  G4 endX  
    std::vector<Float_t>     mcshwr_AncestorendY;      //MC Shower's ancestor  G4 endY  
    std::vector<Float_t>     mcshwr_AncestorendZ;      //MC Shower's ancestor  G4 endZ  
    
    //MC track information
    Int_t     no_mctracks;                         //number of MC tracks in this event.
    //MC track particle information
    std::vector<Int_t>       mctrk_origin;	    //MC track origin information.  
    std::vector<Int_t>       mctrk_pdg;	    //MC track particle PDG code.   
    std::vector<Int_t>       mctrk_TrackId;        //MC track particle G4 track ID.
    std::vector<std::string> mctrk_Process;	    //MC track particle's creation process. 
    std::vector<Float_t>     mctrk_startX;	    //MC track particle G4 startX 
    std::vector<Float_t>     mctrk_startY;	    //MC track particle G4 startY 
    std::vector<Float_t>     mctrk_startZ;	    //MC track particle G4 startZ 
    std::vector<Float_t>     mctrk_endX;		//MC track particle G4 endX 
    std::vector<Float_t>     mctrk_endY;		//MC track particle G4 endY 
    std::vector<Float_t>     mctrk_endZ;		//MC track particle G4 endZ
    std::vector<Float_t>     mctrk_startX_drifted;        //MC track particle first step in TPC x
    std::vector<Float_t>     mctrk_startY_drifted;        //MC track particle first step in TPC y
    std::vector<Float_t>     mctrk_startZ_drifted;        //MC track particle first step in TPC z
    std::vector<Float_t>     mctrk_endX_drifted;          //MC track particle last step in TPC x
    std::vector<Float_t>     mctrk_endY_drifted;          //MC track particle last step in TPC y
    std::vector<Float_t>     mctrk_endZ_drifted;          //MC track particle last step in TPC z
    std::vector<Float_t>     mctrk_len_drifted;           //MC track length within TPC
    std::vector<Float_t>     mctrk_p_drifted;             //MC track momentum at start point in TPC
    std::vector<Float_t>     mctrk_px_drifted;            //MC track x momentum at start point in TPC
    std::vector<Float_t>     mctrk_py_drifted;            //MC track y momentum at start point in TPC
    std::vector<Float_t>     mctrk_pz_drifted;            //MC track z momentum at start point in TPC
    //MC Track mother information
    std::vector<Int_t>       mctrk_Motherpdg;       //MC Track's mother PDG code. 
    std::vector<Int_t>       mctrk_MotherTrkId;     //MC Track's mother G4 track ID.
    std::vector<std::string> mctrk_MotherProcess;   //MC Track's mother creation process. 
    std::vector<Float_t>     mctrk_MotherstartX;    //MC Track's mother  G4 startX .
    std::vector<Float_t>     mctrk_MotherstartY;    //MC Track's mother  G4 startY .
    std::vector<Float_t>     mctrk_MotherstartZ;    //MC Track's mother  G4 startZ .
    std::vector<Float_t>     mctrk_MotherendX;	     //MC Track's mother  G4 endX   .
    std::vector<Float_t>     mctrk_MotherendY;	     //MC Track's mother  G4 endY   .
    std::vector<Float_t>     mctrk_MotherendZ;	     //MC Track's mother  G4 endZ   .
    //MC Track ancestor information
    std::vector<Int_t>       mctrk_Ancestorpdg;       //MC Track's ancestor PDG code. 
    std::vector<Int_t>       mctrk_AncestorTrkId;     //MC Track's ancestor G4 track ID.
    std::vector<std::string> mctrk_AncestorProcess;   //MC Track's ancestor creation process. 
    std::vector<Float_t>     mctrk_AncestorstartX;    //MC Track's ancestor  G4 startX
    std::vector<Float_t>     mctrk_AncestorstartY;    //MC Track's ancestor  G4 startY
    std::vector<Float_t>     mctrk_AncestorstartZ;    //MC Track's ancestor  G4 startZ
    std::vector<Float_t>     mctrk_AncestorendX;      //MC Track's ancestor  G4 endX  
    std::vector<Float_t>     mctrk_AncestorendY;      //MC Track's ancestor  G4 endY  
    std::vector<Float_t>     mctrk_AncestorendZ;      //MC Track's ancestor  G4 endZ    
						      
    // Auxiliary detector variables saved for each geant track
    // This data is saved as a vector (one item per GEANT particle) of C arrays
    // (wrapped in a BoxedArray for technical reasons), one item for each
    // affected detector cell (which one is saved in AuxDetID
    template <typename T>
    using AuxDetMCData_t = std::vector<BoxedArray<T[kMaxAuxDets]>>;
    
    std::vector<UShort_t> NAuxDets;         ///< Number of AuxDets crossed by this particle
    AuxDetMCData_t<Short_t> AuxDetID;       ///< Which AuxDet this particle went through
    AuxDetMCData_t<Float_t> entryX;         ///< Entry X position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryY;         ///< Entry Y position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryZ;         ///< Entry Z position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryT;         ///< Entry T position of particle into AuxDet
    AuxDetMCData_t<Float_t> exitX;          ///< Exit X position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitY;          ///< Exit Y position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitZ;          ///< Exit Z position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitT;          ///< Exit T position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPx;         ///< Exit x momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPy;         ///< Exit y momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPz;         ///< Exit z momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> CombinedEnergyDep; ///< Sum energy of all particles with this trackID (+ID or -ID) in AuxDet
    
    unsigned int bits; ///< complementary information
 
    /// Returns whether we have auxiliary detector data
    bool hasAuxDetector() const { return bits & tdAuxDet; }
    
    /// Returns whether we have Cry data
    bool hasCryInfo() const { return bits & tdCry; }
    
    /// Returns whether we have Genie data
    bool hasGenieInfo() const { return bits & tdGenie; }
    
    /// Returns whether we have MCShower data
    bool hasMCShowerInfo() const { return bits & tdMCshwr; }
    
    /// Returns whether we have MCTrack data
    bool hasMCTrackInfo() const { return bits & tdMCtrk; }
    
    /// Returns whether we have Hit data
    bool hasHitInfo() const { return bits & tdHit; }

    /// Returns whether we have Hit data
    bool hasSWTriggerInfo() const { return bits & tdSWTrigger; }

    /// Returns whether we have OpFilter data
    bool hasOpticalFilterInfo() const { return bits & tdOpticalFilter; }

    /// Returns whether we have Hit data
    bool hasRawDigitInfo() const { return bits & tdRawDigit; }
    
    /// Returns whether we have Sim Channel data
    bool hasSimChannelInfo() const { return bits & tdSimChannel; }
    
    /// Returns whether we have CalWire data
    bool hasCalWireInfo() const { return bits & tdCalWire; }
    
    /// Returns whether we have Raw WaveForm data
    bool hasRawWaveFormInfo() const { return bits & tdRawWaveForm; }

    /// Returns whether we have Calibration WaveForm data
    bool hasCalibWaveFormInfo() const { return bits & tdCalibWaveForm; }

    /// Returns whether we have Track data
    bool hasTrackInfo() const { return bits & tdTrack; }
    
    /// Returns whether we have Shower data
    bool hasShowerInfo() const { return bits & tdShower; }
    
    /// Returns whether we have Vertex data
    bool hasVertexInfo() const { return bits & tdVertex; }

    /// Returns whether we have PFParticle data
    bool hasPFParticleInfo() const { return bits & tdPFParticle; }
    
    /// Returns whether we have Cluster data
    bool hasClusterInfo() const { return bits & tdCluster; }

    /// Returns whether we have Pandora Nu Vertex data
    bool hasPandoraNuVertexInfo() const { return bits & tdPandoraNuVertex; }
    
    /// Returns whether we have Geant data
    bool hasGeantInfo() const { return bits & tdGeant; }

    /// Returns whether we have Flash data
    bool hasFlashInfo() const { return bits & tdFlash; }
    
    /// Sets the specified bits
    void SetBits(unsigned int setbits, bool unset = false)
    { if (unset) bits &= ~setbits; else bits |= setbits; }
      
    /// Constructor; clears all fields
    AnalysisTreeDataStruct(size_t nTrackers = 0, size_t nVertexAlgos = 0, size_t nNuVertexAlgos = 0, size_t nFlashAlgos = 0,
			   std::vector<std::string> const& ShowerAlgos = {}):
      bits(tdDefault) 
    { SetTrackers(nTrackers); SetVertexAlgos(nVertexAlgos); SetNeutrinoVertexAlgos(nNuVertexAlgos); SetFlashAlgos(nFlashAlgos); SetShowerAlgos(ShowerAlgos); Clear(); }

    TrackDataStruct& GetTrackerData(size_t iTracker)
    { return TrackData.at(iTracker); }
    const TrackDataStruct& GetTrackerData(size_t iTracker) const
    { return TrackData.at(iTracker); }
   
    FlashDataStruct& GetFlashData(size_t iTracker)
    { return FlashData.at(iTracker); }
    const FlashDataStruct& GetFlashData(size_t iTracker) const
    { return FlashData.at(iTracker); }
 
    ShowerDataStruct& GetShowerData(size_t iShower)
    { return ShowerData.at(iShower); }
    ShowerDataStruct const& GetShowerData(size_t iShower) const
    { return ShowerData.at(iShower); }
    
    VertexDataStruct& GetVertexData(size_t iVertex)
    { return VertexData.at(iVertex); }
    const VertexDataStruct& GetVertexData(size_t iVertex) const
    { return VertexData.at(iVertex); }

    NeutrinoVertexDataStruct& GetNeutrinoVertexData(size_t iNuVertex)
    { return NeutrinoVertexData.at(iNuVertex); }
    const NeutrinoVertexDataStruct& GetNeutrinoVertexData(size_t iNuVertex) const
    { return NeutrinoVertexData.at(iNuVertex); }

    PFParticleDataStruct& GetPFParticleData()
      { return PFParticleData; }
    const PFParticleDataStruct& GetPFParticleData() const
      { return PFParticleData; }
      
    RawDataStruct& GetRawData()
      { return RawData; }
    const RawDataStruct& GetRawData() const
      { return RawData; }  
      
    CalibDataStruct& GetCalibData()
      { return CalibData; }
    const CalibDataStruct& GetCalibData() const
      { return CalibData; }    
    
    /// Clear all fields if this object (not the tracker algorithm data)
    void ClearLocalData();
    
    /// Clear all fields
    void Clear();
    
    
    /// Allocates data structures for the given number of trackers (no Clear())
    void SetTrackers(size_t nTrackers) { TrackData.resize(nTrackers); }

    /// Allocates data structures for the given number of vertex algos (no Clear())
    void SetVertexAlgos(size_t nVertexAlgos) { VertexData.resize(nVertexAlgos); }

    /// Allocates data structures for the given number of neutrino vertex algos (no Clear())
    void SetNeutrinoVertexAlgos(size_t nNuVertexAlgos) { NeutrinoVertexData.resize(nNuVertexAlgos); }

    ///Allocates data structures for the given number of flash algos (no Clear())
    void SetFlashAlgos(size_t nFlashAlgos) { FlashData.resize(nFlashAlgos); }

    /// Allocates data structures for the given number of trackers (no Clear())
    void SetShowerAlgos(std::vector<std::string> const& ShowerAlgos);
    
    /// Resize the data strutcure for GEANT particles
    void ResizeGEANT(int nParticles);
    
    /// Resize the data strutcure for Genie primaries
    void ResizeGenie(int nPrimaries);
    
    /// Resize the data strutcure for Cry primaries
    void ResizeCry(int nPrimaries);
    
    /// Resize the data strutcure for  MC Showers
    void ResizeMCShower(int nMCShowers);
    
    /// Resize the data strutcure for  MC Tracks
    void ResizeMCTrack(int nMCTracks);
    
    /// Connect this object with a tree
    void SetAddresses(
		      TTree* pTree,
		      std::vector<std::string> const& trackers,
                      std::vector<std::string> const& vertexalgos,
                      std::vector<std::string> const& nuvertexalgos,
                      std::vector<std::string> const& flashalgos,
		      std::vector<std::string> const& showeralgos,
		      bool isCosmics
		      );
    
    
    /// Returns the number of trackers for which data structures are allocated
    size_t GetNTrackers() const { return TrackData.size(); }
   
    /// Returns the number of Vertex algos for which data structures are allocated
    size_t GetNVertexAlgos() const { return VertexData.size(); }

    /// Returns the number of Neutrinoo Vertex algos for which data structures are allocated
    size_t GetNNuVertexAlgos() const { return NeutrinoVertexData.size(); }

    /// Returns the number of Flash algos for which data structures are allocated
    size_t GetNFlashAlgos() const { return FlashData.size(); }

    /// Returns the number of trackers for which data structures are allocated
    size_t GetNShowerAlgos() const { return ShowerData.size(); }
    
    /// Returns the number of hits for which memory is allocated
    size_t GetMaxHits() const { return kMaxHits; }
    
    /// Returns the number of trackers for which memory is allocated
    size_t GetMaxTrackers() const { return TrackData.capacity(); }

    /// Returns the number of vertex algos for which memory is allocated
    size_t GetMaxVertexAlgos() const { return VertexData.capacity(); }

    /// Returns the number of flash algos for which memory is allocated
    size_t GetMaxFlashAlgos() const { return FlashData.capacity(); }
 
    /// Returns the number of trackers for which memory is allocated
    size_t GetMaxShowers() const { return ShowerData.capacity(); }
    
    /// Returns the number of GEANT particles for which memory is allocated
    size_t GetMaxGEANTparticles() const { return MaxGEANTparticles; }
        
    /// Returns the number of GENIE primaries for which memory is allocated
    size_t GetMaxGeniePrimaries() const { return MaxGeniePrimaries; }
    
    
  private:
    /// Little helper functor class to create or reset branches in a tree
    class BranchCreator {
    public:
      TTree* pTree; ///< the tree to be worked on
      BranchCreator(TTree* tree): pTree(tree) {}

      //@{
      /// Create a branch if it does not exist, and set its address
      void operator()
      (std::string name, void* address, std::string leaflist /*, int bufsize = 32000 */)
      {
	if (!pTree) return;
	TBranch* pBranch = pTree->GetBranch(name.c_str());
	if (!pBranch) {
	  pTree->Branch(name.c_str(), address, leaflist.c_str() /*, bufsize */);
	  MF_LOG_DEBUG("AnalysisTreeStructure")
	    << "Creating branch '" << name << " with leaf '" << leaflist << "'";
	}
	else if (pBranch->GetAddress() != address) {
	  pBranch->SetAddress(address);
	  MF_LOG_DEBUG("AnalysisTreeStructure")
	    << "Reassigning address to branch '" << name << "'";
	}
	else {
	  MF_LOG_DEBUG("AnalysisTreeStructure")
	    << "Branch '" << name << "' is fine";
	}
      } // operator()
      void operator()
      (std::string name, void* address, const std::stringstream& leaflist /*, int bufsize = 32000 */)
      { return this->operator() (name, address, leaflist.str() /*, int bufsize = 32000 */); }
      template <typename T>
      void operator()
      (std::string name, std::vector<T>& data, std::string leaflist /*, int bufsize = 32000 */)
      { return this->operator() (name, (void*) data.data(), leaflist /*, int bufsize = 32000 */); }

      template <typename T>
      void operator() (std::string name, std::vector<T>& data)
      {
	// overload for a generic object expressed directly by reference
	// (as opposed to a generic object expressed by a pointer or
	// to a simple leaf sequence specification);
	// TTree::Branch(name, T* obj, Int_t bufsize, splitlevel) and
	// TTree::SetObject() are used.
	if (!pTree) return;
	TBranch* pBranch = pTree->GetBranch(name.c_str());
	if (!pBranch) {
	  pTree->Branch(name.c_str(), &data);
	  // ROOT needs a TClass definition for T in order to create a branch,
	  // se we are sure that at this point the TClass exists;
	  // well, except for when it does not.
	  MF_LOG_DEBUG("AnalysisTreeStructure")
	    << "Creating object branch '" << name
	    << "' with "
	    << (TClass::GetClass(typeid(T))? TClass::GetClass(typeid(T))->ClassName(): "some")
	    << " data type";
	}
	else if
	  (*(reinterpret_cast<std::vector<T>**>(pBranch->GetAddress())) != &data)
          {
            // when an object is provided directly, the address of the object
            // is assigned in TBranchElement::fObject (via TObject::SetObject())
            // and the address itself is set to the address of the fObject
            // member. Here we check that the address of the object in fObject
            // is the same as the address of our current data type
            pBranch->SetObject(&data);
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Reassigning object to branch '" << name << "'";
          }
	else {
	  MF_LOG_DEBUG("AnalysisTreeStructure")
	    << "Branch '" << name << "' is fine";
	}
      } // operator()
      //@}
    }; // class BranchCreator

  }; // class AnalysisTreeDataStruct
  
  
  /// Contains ROOTTreeCode<>::code, ROOT tree character for branch of type T
  template <typename T> struct ROOTTreeCode; // generally undefined
  
  template<> struct ROOTTreeCode<Short_t>  { static constexpr char code = 'S'; };
  template<> struct ROOTTreeCode<Int_t>    { static constexpr char code = 'I'; };
  template<> struct ROOTTreeCode<Double_t> { static constexpr char code = 'D'; };
  
  
  /// Class whose "type" contains the base data type of the container
  template <typename C> struct ContainerValueType; // generally undefined
  
  template <typename A>
  struct ContainerValueType<std::vector<AnalysisTreeDataStruct::BoxedArray<A>>>
  { using type = typename AnalysisTreeDataStruct::BoxedArray<A>::value_type; };
  
  template <typename T>
  struct ContainerValueType<std::vector<T>>
  { using type = typename std::vector<T>::value_type; };
  
  
  /**
   * @brief Creates a simple ROOT tree with tracking and calorimetry information
   * 
   * <h2>Configuration parameters</h2>
   * - <b>UseBuffers</b> (default: false): if enabled, memory is allocated for
   *   tree data for all the run; otherwise, it's allocated on each event, used
   *   and freed; use "true" for speed, "false" to save memory
   * - <b>SaveAuxDetInfo</b> (default: false): if enabled, auxiliary detector
   *   data will be extracted and included in the tree
   */
  class AnalysisTree : public art::EDAnalyzer {

  public:

    explicit AnalysisTree(fhicl::ParameterSet const& pset);
    virtual ~AnalysisTree();

    /// read access to event
    void analyze(const art::Event& evt);
    //  void beginJob() {}
    void beginSubRun(const art::SubRun& sr);
    void endSubRun(const art::SubRun& sr);

  private:

    void   HitsPurity(detinfo::DetectorClocksData const& clockData,
                      std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe);
    double length(const recob::Track& track);
    double driftedLength(detinfo::DetectorPropertiesData const& detProp,
                         const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
    double driftedLength(detinfo::DetectorPropertiesData const& detProp,
                         const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom);
    double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
    double bdist(const TVector3& pos);

    TTree* fTree;
    TTree* fPOT;
      
    // For keeping track of the replacement backtracker
    std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;

    // event information is huge and dynamic;
    // run information is much smaller and we still store it statically
    // in the event
    std::unique_ptr<AnalysisTreeDataStruct> fData;
    //    AnalysisTreeDataStruct::RunData_t RunData;
    AnalysisTreeDataStruct::SubRunData_t SubRunData;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fCryGenModuleLabel;
    std::string fG4ModuleLabel;
    std::string fClusterModuleLabel;
    std::vector<std::string> fPandoraNuVertexModuleLabel;
    std::vector<std::string> fOpFlashModuleLabel;
    std::string fMCShowerModuleLabel;
    std::string fMCTrackModuleLabel;
    std::vector<std::string> fTrackModuleLabel;
    std::string fPFParticleModuleLabel;
    std::vector<std::string> fVertexModuleLabel;
    std::vector<std::string> fShowerModuleLabel;
    std::vector<std::string> fCalorimetryModuleLabel;
    std::vector<std::string> fParticleIDModuleLabel;
    std::vector<std::string> fT0FinderLabel;   
    std::vector<std::string> fT0RecoAnodeCathodePiercingLabel; 
    std::vector<std::string> fFlashT0FinderLabel;
    std::vector<std::string> fMCT0FinderLabel;
    std::vector<std::string> fTrackMCSFitLabel;
    std::string fPOTModuleLabel;
    std::string fCosmicClusterTaggerAssocLabel;
    std::string fSWTriggerLabel;
    std::string fOpticalFilterLabel;
    bool fUseBuffer; ///< whether to use a permanent buffer (faster, huge memory)    
    bool fSaveAuxDetInfo; ///< whether to extract and save auxiliary detector data
    bool fSaveCryInfo; ///whether to extract and save CRY particle data
    bool fSaveGenieInfo; ///whether to extract and save Genie information
    bool fSaveGeantInfo; ///whether to extract and save Geant information
    bool fSaveMCShowerInfo; ///whether to extract and save MC Shower information
    bool fSaveMCTrackInfo; ///whether to extract and save MC Track information
    bool fSaveHitInfo; ///whether to extract and save Hit information
    bool fSaveRawDigitInfo; ///whether to extract and save Raw Digit information
    bool fSaveCalWireInfo; ///whether to extract and save CalWire information
    bool fSaveSimChannelInfo; ///whether to extract and save SimChannel information
    bool fSaveRawWaveFormInfo; ///whether to extract and save Raw WaveForm information    
    bool fSaveCalibWaveFormInfo; ///whether to extract and save Calib WaveForm information    
    bool fSaveTrackInfo; ///whether to extract and save Track information
    bool fSaveVertexInfo; ///whether to extract and save Vertex information
    bool fSaveClusterInfo;  ///whether to extract and save Cluster information
    bool fSavePandoraNuVertexInfo; ///whether to extract and save nu vertex information from Pandora
    bool fSaveFlashInfo;  ///whether to extract and save Flash information
    bool fSaveShowerInfo;  ///whether to extract and save Shower information
    bool fSavePFParticleInfo; ///whether to extract and save PFParticle information
    bool fSaveSWTriggerInfo; ///whether to extract and save software trigger information
    bool fSaveOpticalFilterInfo; ///whether to extract and save optical filter information
    bool fIsOverlay; //whether we are looking at the overlay sample instead of regular MC/data

    std::vector<std::string> fCosmicTaggerAssocLabel;
    std::vector<std::string> fContainmentTaggerAssocLabel;
    std::vector<std::string> fFlashMatchAssocLabel;

    bool bIgnoreMissingShowers; ///< whether to ignore missing shower information
    
    bool isCosmics;      ///< if it contains cosmics
    bool fSaveCaloCosmics; ///< save calorimetry information for cosmics
    float fG4minE;         ///< Energy threshold to save g4 particle info
    float fRawDigitROI;    ///< Define ROI +/-sigma for raw digit range
    float fCalWireROI;     ///< Define ROI +/-sigma for calibrated signal range
    float fRawDigitThresh; ///< Define Raw digit threshold in terms of % of pulse height for raw digit acceptance 
    float fCalWireThresh;  ///< Define cal wire threshold in terms of % of pulse height for cal wire acceptance 
    /// Returns the number of trackers configured
    size_t GetNTrackers() const { return fTrackModuleLabel.size(); }

    size_t GetNVertexAlgos() const { return fVertexModuleLabel.size(); }

    size_t GetNNuVertexAlgos() const { return fPandoraNuVertexModuleLabel.size(); }     
 
    size_t GetNFlashAlgos() const { return fOpFlashModuleLabel.size(); }
 
    /// Returns the number of shower algorithms configured
    size_t GetNShowerAlgos() const { return fShowerModuleLabel.size(); }
    
    /// Returns the name of configured shower algorithms (converted to string)
    std::vector<std::string> GetShowerAlgos() const
    { return { fShowerModuleLabel.begin(), fShowerModuleLabel.end() }; }
    
    /// Creates the structure for the tree data; optionally initializes it
    void CreateData(bool bClearData = false)
    {
      if (!fData) {
	fData.reset
	  (new AnalysisTreeDataStruct(GetNTrackers(), GetNVertexAlgos(), GetNNuVertexAlgos(), GetNFlashAlgos(), GetShowerAlgos()));
	fData->SetBits(AnalysisTreeDataStruct::tdCry,    !fSaveCryInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdGenie,  !fSaveGenieInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdGeant,  !fSaveGeantInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdMCshwr, !fSaveMCShowerInfo); 
	fData->SetBits(AnalysisTreeDataStruct::tdMCtrk,  !fSaveMCTrackInfo); 
	fData->SetBits(AnalysisTreeDataStruct::tdHit,    !fSaveHitInfo);
        fData->SetBits(AnalysisTreeDataStruct::tdSWTrigger,    !fSaveSWTriggerInfo);
        fData->SetBits(AnalysisTreeDataStruct::tdOpticalFilter,    !fSaveOpticalFilterInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdRawDigit,    !fSaveRawDigitInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdCalWire,    !fSaveCalWireInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdSimChannel,    !fSaveSimChannelInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdRawWaveForm,    !fSaveRawWaveFormInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdCalibWaveForm,    !fSaveCalibWaveFormInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdFlash,  !fSaveFlashInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdShower, !fSaveShowerInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdCluster,!fSaveClusterInfo);
        fData->SetBits(AnalysisTreeDataStruct::tdPandoraNuVertex,!fSavePandoraNuVertexInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdTrack,  !fSaveTrackInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdVertex, !fSaveVertexInfo);
	fData->SetBits(AnalysisTreeDataStruct::tdAuxDet, !fSaveAuxDetInfo);
        fData->SetBits(AnalysisTreeDataStruct::tdPFParticle, !fSavePFParticleInfo);
      }
      else {
	fData->SetTrackers(GetNTrackers());
	fData->SetVertexAlgos(GetNVertexAlgos());
        fData->SetNeutrinoVertexAlgos(GetNNuVertexAlgos());
        fData->SetFlashAlgos(GetNFlashAlgos());
	fData->SetShowerAlgos(GetShowerAlgos());

	if (bClearData) fData->Clear();
      }
    } // CreateData()
    
    /// Sets the addresses of all the tree branches, creating the missing ones
    void SetAddresses()
    {
      CheckData(__func__); CheckTree(__func__);
      fData->SetAddresses
	(fTree, fTrackModuleLabel, fVertexModuleLabel, fPandoraNuVertexModuleLabel, fOpFlashModuleLabel, fShowerModuleLabel, isCosmics);
    } // SetAddresses()
    
    /// Sets the addresses of all the tree branches of the specified tracking algo,
    /// creating the missing ones
    void SetTrackerAddresses(size_t iTracker)
    {
      CheckData(__func__); CheckTree(__func__);
      if (iTracker >= fData->GetNTrackers()) {
	throw art::Exception(art::errors::LogicError)
	  << "AnalysisTree::SetTrackerAddresses(): no tracker #" << iTracker
	  << " (" << fData->GetNTrackers() << " available)";
      }
      fData->GetTrackerData(iTracker)
	.SetAddresses(fTree, fTrackModuleLabel[iTracker], isCosmics);
    } // SetTrackerAddresses()


    void SetVertexAddresses(size_t iVertexAlg)
    {
      CheckData(__func__); CheckTree(__func__);
      if (iVertexAlg >= fData->GetNVertexAlgos()) {
	throw art::Exception(art::errors::LogicError)
	  << "AnalysisTree::SetVertexAddresses(): no vertex alg #" << iVertexAlg
	  << " (" << fData->GetNVertexAlgos() << " available)";
      }
      fData->GetVertexData(iVertexAlg)
	.SetAddresses(fTree, fVertexModuleLabel[iVertexAlg], isCosmics);
    } // SetVertexAddresses()
  
    void SetNeutrinoVertexAddresses(size_t iNuVertexAlg)
    {
      CheckData(__func__); CheckTree(__func__);
      if (iNuVertexAlg >= fData->GetNNuVertexAlgos()) {
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::SetNeutrinoVertexAddresses(): no vertex alg #" << iNuVertexAlg
          << " (" << fData->GetNNuVertexAlgos() << " available)";
      }
      fData->GetNeutrinoVertexData(iNuVertexAlg)
        .SetAddresses(fTree, fPandoraNuVertexModuleLabel[iNuVertexAlg], isCosmics);
    } // SetVertexAddresses()
 
    void SetFlashAddresses(size_t iFlashAlg)
    {
      CheckData(__func__); CheckTree(__func__);
      if (iFlashAlg >= fData->GetNFlashAlgos()) {
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::SetFlashAddresses(): no flash alg #" << iFlashAlg
          << " (" << fData->GetNFlashAlgos() << " available)";
      }
      fData->GetFlashData(iFlashAlg)
        .SetAddresses(fTree, fOpFlashModuleLabel[iFlashAlg], isCosmics);
    } // SetVertexAddresses()

 
    /// Sets the addresses of all the tree branches of the specified shower algo,
    /// creating the missing ones
    void SetShowerAddresses(size_t iShower)
    {
      CheckData(__func__); CheckTree(__func__);
      if (iShower >= fData->GetNShowerAlgos()) {
	throw art::Exception(art::errors::LogicError)
	  << "AnalysisTree::SetShowerAddresses(): no shower algo #" << iShower
	  << " (" << fData->GetNShowerAlgos() << " available)";
      }
      fData->GetShowerData(iShower).SetAddresses(fTree);
    } // SetShowerAddresses()

    /// Sets the addresses of the tree branch of the PFParticle,
    /// creating it if missing
    void SetPFParticleAddress()
    {
      CheckData(__func__); CheckTree(__func__);
      fData->GetPFParticleData().SetAddresses(fTree);
    } // SetPFParticleAddress()
    
    /// Sets the addresses of the tree branch of the RawData,
    /// creating it if missing
    void SetRawAddress()
    {
      CheckData(__func__); CheckTree(__func__);
      fData->GetRawData().SetAddresses(fTree);
    } // SetRawAddress()
    
     /// Sets the addresses of the tree branch of the CalibData,
    /// creating it if missing
    void SetCalibAddress()
    {
      CheckData(__func__); CheckTree(__func__);
      fData->GetCalibData().SetAddresses(fTree);
    } // SetCalibAddress()
    
    /// Create the output tree and the data structures, if needed
    void CreateTree(bool bClearData = false);
    
    /// Destroy the local buffers (existing branches will point to invalid address!)
    void DestroyData() { fData.reset(); }
    
    /// Helper function: throws if no data structure is available
    void CheckData(std::string caller) const
    {
      if (fData) return;
      throw art::Exception(art::errors::LogicError)
	<< "AnalysisTree::" << caller << ": no data";
    } // CheckData()
    /// Helper function: throws if no tree is available
    void CheckTree(std::string caller) const
    {
      if (fTree) return;
      throw art::Exception(art::errors::LogicError)
	<< "AnalysisTree::" << caller << ": no tree";
    } // CheckData()
    
    /// Stores the information of shower in slot iShower of showerData
    void FillShower(
		    AnalysisTreeDataStruct::ShowerDataStruct& showerData,
		    size_t iShower, recob::Shower const& showers, const bool fSavePFParticleInfo, 
            const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
		    ) const;
    
    /// Stores the information of all showers into showerData
    void FillShowers(
		     AnalysisTreeDataStruct::ShowerDataStruct& showerData,
		     std::vector<recob::Shower> const& showers, const bool fSavePFParticleInfo, 
             const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
		     ) const;
    
  }; // class microboone::AnalysisTree
} // namespace microboone


namespace { // local namespace
  /// Simple stringstream which empties its buffer on operator() call
  class AutoResettingStringSteam: public std::ostringstream {
  public:
    AutoResettingStringSteam& operator() () { str(""); return *this; }
  }; // class AutoResettingStringSteam

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, ITER to, TYPE value)
  { std::fill(from, to, value); }

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, size_t n, TYPE value)
  { std::fill(from, from + n, value); }

  /// Fills a container with begin()/end() interface
  template <typename CONT, typename V>
  inline void FillWith(CONT& data, const V& value)
  { FillWith(std::begin(data), std::end(data), value); }

} // local namespace

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::RawDataStruct
//---

void microboone::AnalysisTreeDataStruct::RawDataStruct::Resize(size_t raw_nChannels)
{
  MaxRawChannels = raw_nChannels;
  
  raw_channelId.resize(MaxRawChannels);
  raw_wf.resize(MaxRawChannels);
} // microboone::AnalysisTreeDataStruct::RawDataStruct::Resize()

void microboone::AnalysisTreeDataStruct::RawDataStruct::Clear() {
  Resize(MaxRawChannels);
  raw_nchannels = 0;
  
  FillWith(raw_channelId        , -9999  );
  
  for (size_t irawch = 0; irawch < MaxRawChannels; ++irawch){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(raw_wf[irawch]      , -9999);
  } // for Raw Channel
  
} // microboone::AnalysisTreeDataStruct::RawDataStruct::Clear()

void microboone::AnalysisTreeDataStruct::RawDataStruct::SetAddresses(
								       TTree* pTree
								       ) {
  if (MaxRawChannels == 0) return; // no raw channels, no tree!
  
  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);
  AutoResettingStringSteam sstr;
  sstr() << kMaxTicks;
  std::string MaxRawTicksIndexStr("[" + sstr.str() + "]");
  
  std::string BranchName;

  BranchName = "raw_nchannels";
  CreateBranch(BranchName, &raw_nchannels, BranchName + "/S");
  std::string rawnChannelsIndexStr = "[" + BranchName + "]";

  BranchName = "raw_channelId";
  CreateBranch(BranchName, raw_channelId, BranchName + rawnChannelsIndexStr + "/S"); 
  
  BranchName = "raw_wf";
  CreateBranch(BranchName, raw_wf, BranchName + rawnChannelsIndexStr + MaxRawTicksIndexStr + "/S");  
  
} // microboone::AnalysisTreeDataStruct::RawDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::CalibDataStruct
//---

void microboone::AnalysisTreeDataStruct::CalibDataStruct::Resize(size_t calib_nChannels)
{
  MaxCalibChannels = calib_nChannels;
  
  calib_channelId.resize(MaxCalibChannels);
  calib_wf.resize(MaxCalibChannels);
} // microboone::AnalysisTreeDataStruct::CalibDataStruct::Resize()

void microboone::AnalysisTreeDataStruct::CalibDataStruct::Clear() {
  Resize(MaxCalibChannels);
  calib_nchannels = 0;
  
  FillWith(calib_channelId        , -9999  );
  
  for (size_t icalch = 0; icalch < MaxCalibChannels; ++icalch){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(calib_wf[icalch]      , -9999);
  } // for Calib Channel
  
} // microboone::AnalysisTreeDataStruct::CalibDataStruct::Clear()

void microboone::AnalysisTreeDataStruct::CalibDataStruct::SetAddresses(
								       TTree* pTree
								       ) {
  if (MaxCalibChannels == 0) return; // no calib channels, no tree!
  
  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);
  AutoResettingStringSteam sstr;
  sstr() << kMaxTicks;
  std::string MaxCalibTicksIndexStr("[" + sstr.str() + "]");
  
  std::string BranchName;

  BranchName = "calib_nchannels";
  CreateBranch(BranchName, &calib_nchannels, BranchName + "/S");
  std::string calibnChannelsIndexStr = "[" + BranchName + "]";

  BranchName = "calib_channelId";
  CreateBranch(BranchName, calib_channelId, BranchName + calibnChannelsIndexStr + "/S"); 
  
  BranchName = "calib_wf";
  CreateBranch(BranchName, calib_wf, BranchName + calibnChannelsIndexStr + MaxCalibTicksIndexStr + "/S");  
  
} // microboone::AnalysisTreeDataStruct::CalibDataStruct::SetAddresses()


//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::TrackDataStruct
//---

void microboone::AnalysisTreeDataStruct::TrackDataStruct::Resize(size_t nTracks)
{
  MaxTracks = nTracks;
  
  trkId.resize(MaxTracks);
  trkncosmictags_tagger.resize(MaxTracks);
  trkcosmicscore_tagger.resize(MaxTracks);
  trkcosmictype_tagger.resize(MaxTracks);
  trkncosmictags_containmenttagger.resize(MaxTracks);
  trkcosmicscore_containmenttagger.resize(MaxTracks);
  trkcosmictype_containmenttagger.resize(MaxTracks);
  trkncosmictags_flashmatch.resize(MaxTracks);
  trkcosmicscore_flashmatch.resize(MaxTracks);
  trkcosmictype_flashmatch.resize(MaxTracks);
  trkstartx.resize(MaxTracks);
  trkstarty.resize(MaxTracks);
  trkstartz.resize(MaxTracks);
  trkstartd.resize(MaxTracks);
  trkendx.resize(MaxTracks);
  trkendy.resize(MaxTracks);
  trkendz.resize(MaxTracks);
  trkendd.resize(MaxTracks);
  trkACpierceT0.resize(MaxTracks);    
  trkflashT0.resize(MaxTracks);  
  trktrueT0.resize(MaxTracks);  
  trktheta.resize(MaxTracks);
  trkphi.resize(MaxTracks);
  trkstartdcosx.resize(MaxTracks);
  trkstartdcosy.resize(MaxTracks);
  trkstartdcosz.resize(MaxTracks);
  trkenddcosx.resize(MaxTracks);
  trkenddcosy.resize(MaxTracks);
  trkenddcosz.resize(MaxTracks);
  trkthetaxz.resize(MaxTracks);
  trkthetayz.resize(MaxTracks);
  trkmom.resize(MaxTracks);
  trkmomrange.resize(MaxTracks);
  trkmommschi2.resize(MaxTracks);
  trkmommsllhd.resize(MaxTracks);  
  trkmcsfwdmom.resize(MaxTracks);
  trkmcsfwdll.resize(MaxTracks);
  trkmcsfwderr.resize(MaxTracks);
  trkmcsbwdmom.resize(MaxTracks);
  trkmcsbwdll.resize(MaxTracks);
  trkmcsbwderr.resize(MaxTracks);
  trklen.resize(MaxTracks);
  trksvtxid.resize(MaxTracks);
  trkevtxid.resize(MaxTracks);
  trkntraj.resize(MaxTracks);
  // PID variables
  trkpidpdg.resize(MaxTracks);
  trkpidchi.resize(MaxTracks);
  trkpidchipr.resize(MaxTracks);
  trkpidchika.resize(MaxTracks);
  trkpidchipi.resize(MaxTracks);
  trkpidchimu.resize(MaxTracks);
  trkpidpida.resize(MaxTracks);
  trkpidbestplane.resize(MaxTracks);
  
  trkke.resize(MaxTracks);
  trkrange.resize(MaxTracks);
  trkidtruth.resize(MaxTracks);
  trkorigin.resize(MaxTracks);
  trkpdgtruth.resize(MaxTracks);
  trkefftruth.resize(MaxTracks);
  trkpurtruth.resize(MaxTracks);
  trkpurity.resize(MaxTracks);
  trkcompleteness.resize(MaxTracks);
  trkg4id.resize(MaxTracks);
  trkorig.resize(MaxTracks);
  trkpitchc.resize(MaxTracks);
  ntrkhits.resize(MaxTracks);
  
  trkdedx.resize(MaxTracks);
  trkdqdx.resize(MaxTracks);
  trkresrg.resize(MaxTracks);
  trkxyz.resize(MaxTracks);
  
  trkhasPFParticle.resize(MaxTracks);
  trkPFParticleID.resize(MaxTracks);
  
} // microboone::AnalysisTreeDataStruct::TrackDataStruct::Resize()

void microboone::AnalysisTreeDataStruct::TrackDataStruct::Clear() {
  Resize(MaxTracks);
  ntracks = 0;
  
  FillWith(trkId        , -9999  );
  FillWith(trkncosmictags_tagger, -9999  );
  FillWith(trkcosmicscore_tagger, -99999.);
  FillWith(trkcosmictype_tagger, -9999  );
  FillWith(trkncosmictags_containmenttagger, -9999  );
  FillWith(trkcosmicscore_containmenttagger, -99999.);
  FillWith(trkcosmictype_containmenttagger, -9999  );
  FillWith(trkncosmictags_flashmatch, -9999  );
  FillWith(trkcosmicscore_flashmatch, -99999.);
  FillWith(trkcosmictype_flashmatch, -9999  );
  FillWith(trkstartx    , -99999.);
  FillWith(trkstarty    , -99999.);
  FillWith(trkstartz    , -99999.);
  FillWith(trkstartd    , -99999.);
  FillWith(trkendx      , -99999.);
  FillWith(trkendy      , -99999.);
  FillWith(trkendz      , -99999.);
  FillWith(trkendd      , -99999.);
  FillWith(trkACpierceT0, -99999.);    
  FillWith(trkflashT0   , -99999.);  
  FillWith(trktrueT0    , -99999.);  
  FillWith(trkg4id      , -99999 );
  FillWith(trkpurity    , -99999.);
  FillWith(trkcompleteness, -99999.);
  FillWith(trkorig      , -99999 );
  FillWith(trktheta     , -99999.);
  FillWith(trkphi       , -99999.);
  FillWith(trkstartdcosx, -99999.);
  FillWith(trkstartdcosy, -99999.);
  FillWith(trkstartdcosz, -99999.);
  FillWith(trkenddcosx  , -99999.);
  FillWith(trkenddcosy  , -99999.);
  FillWith(trkenddcosz  , -99999.);
  FillWith(trkthetaxz   , -99999.);
  FillWith(trkthetayz   , -99999.);
  FillWith(trkmom       , -99999.);
  FillWith(trkmomrange  , -99999.);  
  FillWith(trkmommschi2 , -99999.);  
  FillWith(trkmommsllhd , -99999.);  
  FillWith(trkmcsfwdmom , -99999.);
  FillWith(trkmcsfwdll  , -99999.);
  FillWith(trkmcsfwderr , -99999.);
  FillWith(trkmcsbwdmom , -99999.);
  FillWith(trkmcsbwdll  , -99999.);
  FillWith(trkmcsbwderr , -99999.);
  FillWith(trklen       , -99999.);
  FillWith(trksvtxid    , -1);
  FillWith(trkevtxid    , -1);
  FillWith(trkntraj     , -99999 );
  FillWith(trkpidbestplane, -1); 
  
  FillWith(trkhasPFParticle, -1);
  FillWith(trkPFParticleID , -1);
 
  for (size_t iTrk = 0; iTrk < MaxTracks; ++iTrk){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(trkke[iTrk]      , -99999.);
    FillWith(trkrange[iTrk]   , -99999.);
    FillWith(trkidtruth[iTrk] , -99999 );
    FillWith(trkorigin[iTrk]  , -1 );
    FillWith(trkpdgtruth[iTrk], -99999 );
    FillWith(trkefftruth[iTrk], -99999.);
    FillWith(trkpurtruth[iTrk], -99999.);
    FillWith(trkpitchc[iTrk]  , -99999.);
    FillWith(ntrkhits[iTrk]   ,  -9999 );
    
    FillWith(trkdedx[iTrk], 0.);
    FillWith(trkdqdx[iTrk], 0.);
    FillWith(trkresrg[iTrk], 0.);
    
    FillWith(trkxyz[iTrk], 0.);
 
    FillWith(trkpidpdg[iTrk]    , -1);
    FillWith(trkpidchi[iTrk]    , -99999.);
    FillWith(trkpidchipr[iTrk]  , -99999.);
    FillWith(trkpidchika[iTrk]  , -99999.);
    FillWith(trkpidchipi[iTrk]  , -99999.);
    FillWith(trkpidchimu[iTrk]  , -99999.);
    FillWith(trkpidpida[iTrk]   , -99999.);
  } // for track
  
} // microboone::AnalysisTreeDataStruct::TrackDataStruct::Clear()


void microboone::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses(
								       TTree* pTree, std::string tracker, bool isCosmics
								       ) {
  if (MaxTracks == 0) return; // no tracks, no tree!
  
  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxTrackHits;
  std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");
  
  std::string TrackLabel = tracker;
  std::string BranchName;

  BranchName = "ntracks_" + TrackLabel;
  CreateBranch(BranchName, &ntracks, BranchName + "/S");
  std::string NTracksIndexStr = "[" + BranchName + "]";
  
  BranchName = "trkId_" + TrackLabel;
  CreateBranch(BranchName, trkId, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkncosmictags_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_tagger, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_tagger, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_tagger, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkncosmictags_containmenttagger_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_containmenttagger, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkcosmicscore_containmenttagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_containmenttagger, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkcosmictype_containmenttagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_containmenttagger, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkncosmictags_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_flashmatch, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkke_" + TrackLabel;
  CreateBranch(BranchName, trkke, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkrange_" + TrackLabel;
  CreateBranch(BranchName, trkrange, BranchName + NTracksIndexStr + "[3]/F");
   
  BranchName = "trkidtruth_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkorigin_" + TrackLabel;
  CreateBranch(BranchName, trkorigin, BranchName + NTracksIndexStr + "[3]/S");
  
  BranchName = "trkpdgtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpdgtruth, BranchName + NTracksIndexStr + "[3]/I");
  
  BranchName = "trkefftruth_" + TrackLabel;
  CreateBranch(BranchName, trkefftruth, BranchName + NTracksIndexStr + "[3]/F");
 
  BranchName = "trkpurtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpurtruth, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkpitchc_" + TrackLabel;
  CreateBranch(BranchName, trkpitchc, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "ntrkhits_" + TrackLabel;
  CreateBranch(BranchName, ntrkhits, BranchName + NTracksIndexStr + "[3]/S");
  
  if (!isCosmics){
    BranchName = "trkdedx_" + TrackLabel;
    CreateBranch(BranchName, trkdedx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
    BranchName = "trkdqdx_" + TrackLabel;
    CreateBranch(BranchName, trkdqdx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    BranchName = "trkresrg_" + TrackLabel;
    CreateBranch(BranchName, trkresrg, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    BranchName = "trkxyz_" + TrackLabel;
    CreateBranch(BranchName, trkxyz, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "[3]" + "/F");
  }

  BranchName = "trkstartx_" + TrackLabel;
  CreateBranch(BranchName, trkstartx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstarty_" + TrackLabel;
  CreateBranch(BranchName, trkstarty, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartz_" + TrackLabel;
  CreateBranch(BranchName, trkstartz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartd_" + TrackLabel;
  CreateBranch(BranchName, trkstartd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendx_" + TrackLabel;
  CreateBranch(BranchName, trkendx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendy_" + TrackLabel;
  CreateBranch(BranchName, trkendy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendz_" + TrackLabel;
  CreateBranch(BranchName, trkendz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendd_" + TrackLabel;
  CreateBranch(BranchName, trkendd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkACpierceT0_" + TrackLabel;
  CreateBranch(BranchName, trkACpierceT0, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkflashT0_" + TrackLabel;
  CreateBranch(BranchName, trkflashT0, BranchName + NTracksIndexStr + "/F");

  BranchName = "trktrueT0_" + TrackLabel;
  CreateBranch(BranchName, trktrueT0, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkg4id_" + TrackLabel;
  CreateBranch(BranchName, trkg4id, BranchName + NTracksIndexStr + "/I");

  BranchName = "trkorig_" + TrackLabel;
  CreateBranch(BranchName, trkorig, BranchName + NTracksIndexStr + "/I");
  
  BranchName = "trkpurity_" + TrackLabel;
  CreateBranch(BranchName, trkpurity, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcompleteness_" + TrackLabel;
  CreateBranch(BranchName, trkcompleteness, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trktheta_" + TrackLabel;
  CreateBranch(BranchName, trktheta, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkphi_" + TrackLabel;
  CreateBranch(BranchName, trkphi, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosx_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosy_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosz_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosx_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosy_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosz_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetaxz_" + TrackLabel;
  CreateBranch(BranchName, trkthetaxz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetayz_" + TrackLabel;
  CreateBranch(BranchName, trkthetayz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmom_" + TrackLabel;
  CreateBranch(BranchName, trkmom, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmomrange_" + TrackLabel;
  CreateBranch(BranchName, trkmomrange, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommschi2_" + TrackLabel;
  CreateBranch(BranchName, trkmommschi2, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommsllhd_" + TrackLabel;
  CreateBranch(BranchName, trkmommsllhd, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsfwdmom_" + TrackLabel;
  CreateBranch(BranchName, trkmcsfwdmom, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsfwdll_" + TrackLabel;
  CreateBranch(BranchName, trkmcsfwdll, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsfwderr_" + TrackLabel;
  CreateBranch(BranchName, trkmcsfwderr, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsbwdmom_" + TrackLabel;
  CreateBranch(BranchName, trkmcsbwdmom, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsbwdll_" + TrackLabel;
  CreateBranch(BranchName, trkmcsbwdll, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmcsbwderr_" + TrackLabel;
  CreateBranch(BranchName, trkmcsbwderr, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trklen_" + TrackLabel;
  CreateBranch(BranchName, trklen, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trksvtxid_" + TrackLabel;
  CreateBranch(BranchName, trksvtxid, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkevtxid_" + TrackLabel;
  CreateBranch(BranchName, trkevtxid, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkntraj_" + TrackLabel;
  CreateBranch(BranchName, trkntraj, BranchName + NTracksIndexStr + "/I");

  BranchName = "trkpidpdg_" + TrackLabel;
  CreateBranch(BranchName, trkpidpdg, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkpidchi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipr_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipr, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchika_" + TrackLabel;
  CreateBranch(BranchName, trkpidchika, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchimu_" + TrackLabel;
  CreateBranch(BranchName, trkpidchimu, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidpida_" + TrackLabel;
  CreateBranch(BranchName, trkpidpida, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidbestplane_" + TrackLabel;
  CreateBranch(BranchName, trkpidbestplane, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkhasPFParticle_" + TrackLabel;
  CreateBranch(BranchName, trkhasPFParticle, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkPFParticleID_" + TrackLabel;
  CreateBranch(BranchName, trkPFParticleID, BranchName + NTracksIndexStr + "/S");

} // microboone::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses()



void microboone::AnalysisTreeDataStruct::VertexDataStruct::Resize(size_t nVertices)
{
  MaxVertices = nVertices;
  vtxId.resize(MaxVertices);
  vtxx.resize(MaxVertices);
  vtxy.resize(MaxVertices);
  vtxz.resize(MaxVertices);
  
  vtxhasPFParticle.resize(MaxVertices);
  vtxPFParticleID.resize(MaxVertices);
}

void microboone::AnalysisTreeDataStruct::VertexDataStruct::Clear() {
  Resize(MaxVertices);
  nvtx = -9999;

  FillWith(vtxId       , -9999  );
  FillWith(vtxx        , -9999  );
  FillWith(vtxy        , -9999  );
  FillWith(vtxz        , -9999  );
  FillWith(vtxhasPFParticle, -1  );
  FillWith(vtxPFParticleID , -1  );
}

void microboone::AnalysisTreeDataStruct::VertexDataStruct::SetAddresses(
									TTree* pTree, std::string alg, bool isCosmics
									) {
  if (MaxVertices == 0) return; // no tracks, no tree!

  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;

  std::string VertexLabel = alg;
  std::string BranchName;

  BranchName = "nvtx_" + VertexLabel;
  CreateBranch(BranchName, &nvtx, BranchName + "/S");
  std::string NVertexIndexStr = "[" + BranchName + "]";

  BranchName = "vtxId_" + VertexLabel;
  CreateBranch(BranchName, vtxId, BranchName + NVertexIndexStr + "/S");

  BranchName = "vtxx_" + VertexLabel;
  CreateBranch(BranchName, vtxx, BranchName + NVertexIndexStr + "/F");

  BranchName = "vtxy_" + VertexLabel;
  CreateBranch(BranchName, vtxy, BranchName + NVertexIndexStr + "/F");

  BranchName = "vtxz_" + VertexLabel;
  CreateBranch(BranchName, vtxz, BranchName + NVertexIndexStr + "/F");
  
  BranchName = "vtxhasPFParticle_" + VertexLabel;
  CreateBranch(BranchName, vtxhasPFParticle, BranchName + NVertexIndexStr + "/S");
  
  BranchName = "vtxPFParticleID_" + VertexLabel;
  CreateBranch(BranchName, vtxPFParticleID, BranchName + NVertexIndexStr + "/S");
}


//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::NeutrinoVertexDataStruct
//---

void microboone::AnalysisTreeDataStruct::NeutrinoVertexDataStruct::Resize(size_t nVertices)
{ 
  MaxVertices = nVertices;
  vtxId.resize(MaxVertices);
  vtxx.resize(MaxVertices);
  vtxy.resize(MaxVertices);
  vtxz.resize(MaxVertices);
  vtxpdg.resize(MaxVertices);

  vtxhasPFParticle.resize(MaxVertices);
  vtxPFParticleID.resize(MaxVertices);
}

void microboone::AnalysisTreeDataStruct::NeutrinoVertexDataStruct::Clear() {
  Resize(MaxVertices);
  nvtx = -9999;

  FillWith(vtxId       , -9999  );
  FillWith(vtxx        , -9999  );
  FillWith(vtxy        , -9999  );
  FillWith(vtxz        , -9999  );
  FillWith(vtxpdg      , -9999  );
  FillWith(vtxhasPFParticle, -1  );
  FillWith(vtxPFParticleID , -1  );
}

void microboone::AnalysisTreeDataStruct::NeutrinoVertexDataStruct::SetAddresses(
                                                                        TTree* pTree, std::string alg, bool isCosmics
                                                                        ) {
  if (MaxVertices == 0) return; // no tracks, no tree!

  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;

  std::string VertexLabel = alg;
  std::string BranchName;

  BranchName = "nnuvtx_" + VertexLabel;
  CreateBranch(BranchName, &nvtx, BranchName + "/S");
  std::string NVertexIndexStr = "[" + BranchName + "]";

  BranchName = "nuvtxId_" + VertexLabel;
  CreateBranch(BranchName, vtxId, BranchName + NVertexIndexStr + "/S");

  BranchName = "nuvtxx_" + VertexLabel;
  CreateBranch(BranchName, vtxx, BranchName + NVertexIndexStr + "/F");

  BranchName = "nuvtxy_" + VertexLabel;
  CreateBranch(BranchName, vtxy, BranchName + NVertexIndexStr + "/F");

  BranchName = "nuvtxz_" + VertexLabel;
  CreateBranch(BranchName, vtxz, BranchName + NVertexIndexStr + "/F");

  BranchName = "nuvtxpdg_" + VertexLabel;
  CreateBranch(BranchName, vtxpdg, BranchName + NVertexIndexStr + "/I");

  BranchName = "nuvtxhasPFParticle_" + VertexLabel;
  CreateBranch(BranchName, vtxhasPFParticle, BranchName + NVertexIndexStr + "/S");

  BranchName = "nuvtxPFParticleID_" + VertexLabel;
  CreateBranch(BranchName, vtxPFParticleID, BranchName + NVertexIndexStr + "/S");
}

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::FlashDataStruct
//---

void microboone::AnalysisTreeDataStruct::FlashDataStruct::Resize(size_t nFlashes)
{
  MaxFlashes = nFlashes;

  flsTime.resize(MaxFlashes);
  flsPe.resize(MaxFlashes);
  flsPePerOpDet.resize(MaxFlashes);
  flsXcenter.resize(MaxFlashes);
  flsYcenter.resize(MaxFlashes);
  flsZcenter.resize(MaxFlashes);
  flsYwidth.resize(MaxFlashes);
  flsZwidth.resize(MaxFlashes);
  flsTwidth.resize(MaxFlashes);

  size_t nOpDets = this->GetNOpDet();
  if (kNOpDets != nOpDets) {
    mf::LogError("AnalysisTree") 
        << "Number of optical detectors from geometry services is" << nOpDets 
        << ", which is different from the one expected of " << kNOpDets 
        << ". Check variale kNOpDets in analysis tree module.";
  }

  for (size_t opfls = 0; opfls < MaxFlashes; opfls++) {
    FillWith(flsPePerOpDet[opfls]   , -99999.);
  }
}

void microboone::AnalysisTreeDataStruct::FlashDataStruct::Clear() {
  Resize(MaxFlashes);
  nfls = -9999;

  FillWith(flsTime        , -9999  );
  FillWith(flsPe          , -9999  );
  FillWith(flsXcenter     , -9999  );
  FillWith(flsYcenter     , -9999  );
  FillWith(flsZcenter     , -9999  );
  FillWith(flsYwidth      , -9999  );
  FillWith(flsZwidth      , -9999  );
  FillWith(flsTwidth      , -9999  );

  for (size_t opfls = 0; opfls < MaxFlashes; opfls++) {
    FillWith(flsPePerOpDet[opfls]  , -9999  );
  }
}

void microboone::AnalysisTreeDataStruct::FlashDataStruct::SetAddresses(
                                                                        TTree* pTree, std::string alg, bool isCosmics
                                                                        ) {
  if (MaxFlashes == 0) return; // no tracks, no tree!

  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;

  std::string FlashLabel = alg;
  std::string BranchName;

  std::ostringstream oss;
  oss << "[" << kNOpDets << "]";
  std::string NOpDetIndexString = oss.str();

  BranchName = "nfls_" + FlashLabel;
  CreateBranch(BranchName, &nfls, BranchName + "/S");
  std::string NFlashIndexStr = "[" + BranchName + "]";

  BranchName = "flsTime_" + FlashLabel;
  CreateBranch(BranchName, flsTime, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsPe_" + FlashLabel;
  CreateBranch(BranchName, flsPe, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsPePerOpDet_" + FlashLabel;
  CreateBranch(BranchName, flsPePerOpDet, BranchName + NFlashIndexStr + NOpDetIndexString + "/D");

  BranchName = "flsXcenter_" + FlashLabel;
  CreateBranch(BranchName, flsXcenter, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsYcenter_" + FlashLabel;
  CreateBranch(BranchName, flsYcenter, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsZcenter_" + FlashLabel;
  CreateBranch(BranchName, flsZcenter, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsYwidth_" + FlashLabel;
  CreateBranch(BranchName, flsYwidth, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsZwidth_" + FlashLabel;
  CreateBranch(BranchName, flsZwidth, BranchName + NFlashIndexStr + "/F");

  BranchName = "flsTwidth_" + FlashLabel;
  CreateBranch(BranchName, flsTwidth, BranchName + NFlashIndexStr + "/F");
}



//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::PFParticleDataStruct
//---

void microboone::AnalysisTreeDataStruct::PFParticleDataStruct::Resize(size_t nPFParticles)
{
  MaxPFParticles = nPFParticles;
  
  pfp_selfID.resize(MaxPFParticles);
  pfp_isPrimary.resize(MaxPFParticles);
  pfp_numDaughters.resize(MaxPFParticles);
  pfp_daughterIDs.resize(MaxPFParticles);
  pfp_parentID.resize(MaxPFParticles);
  pfp_vertexID.resize(MaxPFParticles);
  pfp_isShower.resize(MaxPFParticles);
  pfp_isTrack.resize(MaxPFParticles);
  pfp_trackID.resize(MaxPFParticles);
  pfp_showerID.resize(MaxPFParticles);
  pfp_pdgCode.resize(MaxPFParticles);
  pfp_numClusters.resize(MaxPFParticles);
  pfp_clusterIDs.resize(MaxPFParticles);
  pfp_isNeutrino.resize(MaxPFParticles);
}

void microboone::AnalysisTreeDataStruct::PFParticleDataStruct::Clear() {
  Resize(MaxPFParticles);
  
  nPFParticles = -9999;
  FillWith(pfp_selfID, -9999);
  FillWith(pfp_isPrimary, -9999);
  FillWith(pfp_numDaughters, -9999);
  FillWith(pfp_parentID, -9999);
  FillWith(pfp_vertexID, -9999);
  FillWith(pfp_isShower, -9999);
  FillWith(pfp_isTrack, -9999);
  FillWith(pfp_trackID, -9999);
  FillWith(pfp_showerID, -9999);
  FillWith(pfp_pdgCode, -9999);
  FillWith(pfp_isNeutrino, -9999);
  pfp_numNeutrinos = -9999;
  FillWith(pfp_neutrinoIDs, -9999);

  for (size_t iPFParticle = 0; iPFParticle < MaxPFParticles; ++iPFParticle){
    // the following are BoxedArrays;
    // their iterators traverse all the array dimensions
    FillWith(pfp_daughterIDs[iPFParticle], -9999);
    FillWith(pfp_clusterIDs[iPFParticle], -9999);
  }
}

void microboone::AnalysisTreeDataStruct::PFParticleDataStruct::SetAddresses(
  TTree* pTree
) {
    
  if (MaxPFParticles == 0) { return; } // no PFParticles, no tree

  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxNDaughtersPerPFP;
  std::string MaxNDaughtersIndexStr("[" + sstr.str() + "]");
  
  sstr.str("");
  sstr() << kMaxNClustersPerPFP;
  std::string MaxNClustersIndexStr("[" + sstr.str() + "]");
  
  sstr.str("");
  sstr() << kMaxNPFPNeutrinos;
  std::string MaxNNeutrinosIndexStr("[" + sstr.str() + "]");
  
  std::string BranchName;

  BranchName = "nPFParticles";
  CreateBranch(BranchName, &nPFParticles, BranchName + "/S");
  std::string NPFParticleIndexStr = "[" + BranchName + "]";

  BranchName = "pfp_selfID";
  CreateBranch(BranchName, pfp_selfID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_isPrimary";
  CreateBranch(BranchName, pfp_isPrimary, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_numDaughters";
  CreateBranch(BranchName, pfp_numDaughters, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_daughterIDs";
  CreateBranch(BranchName, pfp_daughterIDs, BranchName + NPFParticleIndexStr + MaxNDaughtersIndexStr + "/S");
  
  BranchName = "pfp_parentID";
  CreateBranch(BranchName, pfp_parentID, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_vertexID";
  CreateBranch(BranchName, pfp_vertexID, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_isShower";
  CreateBranch(BranchName, pfp_isShower, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_isTrack";
  CreateBranch(BranchName, pfp_isTrack, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_trackID";
  CreateBranch(BranchName, pfp_trackID, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_showerID";
  CreateBranch(BranchName, pfp_showerID, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_pdgCode";
  CreateBranch(BranchName, pfp_pdgCode, BranchName + NPFParticleIndexStr + "/I");
  
  BranchName = "pfp_numClusters";
  CreateBranch(BranchName, pfp_numClusters, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_clusterIDs";
  CreateBranch(BranchName, pfp_clusterIDs, BranchName + NPFParticleIndexStr + MaxNClustersIndexStr + "/S");
  
  BranchName = "pfp_isNeutrino";
  CreateBranch(BranchName, pfp_isNeutrino, BranchName + NPFParticleIndexStr + "/S");
  
  BranchName = "pfp_numNeutrinos";
  CreateBranch(BranchName, &pfp_numNeutrinos, BranchName + "/S");
  
  BranchName = "pfp_neutrinoIDs";
  CreateBranch(BranchName, pfp_neutrinoIDs, BranchName + MaxNNeutrinosIndexStr + "/S");
}

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::ShowerDataStruct
//---

void microboone::AnalysisTreeDataStruct::ShowerDataStruct::Resize
(size_t nShowers)
{
  MaxShowers = nShowers;
  
  showerID.resize(MaxShowers);
  shwr_bestplane.resize(MaxShowers);
  shwr_length.resize(MaxShowers);
  shwr_theta.resize(MaxShowers);
  shwr_thetaxz.resize(MaxShowers);
  shwr_thetayz.resize(MaxShowers);
  shwr_phi.resize(MaxShowers);
  shwr_startdcosx.resize(MaxShowers);
  shwr_startdcosy.resize(MaxShowers);
  shwr_startdcosz.resize(MaxShowers);
  shwr_startx.resize(MaxShowers);
  shwr_starty.resize(MaxShowers);
  shwr_startz.resize(MaxShowers);
  shwr_totEng.resize(MaxShowers);
  shwr_dedx.resize(MaxShowers);
  shwr_mipEng.resize(MaxShowers);
  
  shwr_hasPFParticle.resize(MaxShowers);
  shwr_PFParticleID.resize(MaxShowers);
  
} // microboone::AnalysisTreeDataStruct::ShowerDataStruct::Resize()

void microboone::AnalysisTreeDataStruct::ShowerDataStruct::Clear() {
  Resize(MaxShowers);
  nshowers = 0;
  
  FillWith(showerID,         -9999 );
  FillWith(shwr_bestplane,   -9999 );
  FillWith(shwr_length,     -99999.);
  FillWith(shwr_theta,      -99999.);
  FillWith(shwr_thetaxz,    -99999.);
  FillWith(shwr_thetayz,    -99999.);
  FillWith(shwr_phi,        -99999.);
  FillWith(shwr_startdcosx, -99999.);
  FillWith(shwr_startdcosy, -99999.);
  FillWith(shwr_startdcosz, -99999.);
  FillWith(shwr_startx,     -99999.);
  FillWith(shwr_starty,     -99999.);
  FillWith(shwr_startz,     -99999.);
  
  FillWith(shwr_hasPFParticle, -1);
  FillWith(shwr_PFParticleID,  -1);
 
  for (size_t iShw = 0; iShw < MaxShowers; ++iShw){
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(shwr_totEng[iShw], -99999.);
    FillWith(shwr_dedx[iShw],   -99999.);
    FillWith(shwr_mipEng[iShw], -99999.);
  } // for shower
  
} // microboone::AnalysisTreeDataStruct::ShowerDataStruct::Clear()


void microboone::AnalysisTreeDataStruct::ShowerDataStruct::MarkMissing
(TTree* pTree)
{
  // here we implement the policy prescription for a missing set of showers;
  // this means that no shower data product was found in the event,
  // yet the user has accepted to go on.
  // We now need to mark this product in a unmistakably clear way, so that it
  // is not confused with a valid collection of an event where no showers
  // were reconstructed, not as a list of valid showers.
  // The prescription currently implemented is:
  // - have only one shower in the list;
  // - set the ID of that shower as -9999
  //
  
  // first set the data structures to contain one invalid shower:
  SetMaxShowers(1); // includes resize to a set of one shower
  Clear(); // initializes all the showers in the set (one) as invalid
  // now set the tree addresses to the newly allocated memory;
  // this creates the tree branches in case they are not there yet
  SetAddresses(pTree);
  
  // then, set the variables so that ROOT tree knows there is one shower only
  nshowers = 1;
  
} // microboone::AnalysisTreeDataStruct::ShowerDataStruct::MarkMissing()


void microboone::AnalysisTreeDataStruct::ShowerDataStruct::SetAddresses
(TTree* pTree)
{
  if (MaxShowers == 0) return; // no showers, no tree!
  
  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxShowerHits;
  std::string MaxShowerHitsIndexStr("[" + sstr.str() + "]");
  
  std::string ShowerLabel = Name();
  std::string BranchName;

  BranchName = "nshowers_" + ShowerLabel;
  CreateBranch(BranchName, &nshowers, BranchName + "/S");
  std::string NShowerIndexStr = "[" + BranchName + "]";
  
  BranchName = "showerID_" + ShowerLabel;
  CreateBranch(BranchName, showerID, BranchName + NShowerIndexStr + "/S");
  
  BranchName = "shwr_bestplane_" + ShowerLabel;
  CreateBranch(BranchName, shwr_bestplane, BranchName + NShowerIndexStr + "/S");
  
  BranchName = "shwr_length_" + ShowerLabel;
  CreateBranch(BranchName, shwr_length, BranchName + NShowerIndexStr + "/F");
 
  BranchName = "shwr_theta_" + ShowerLabel;
  CreateBranch(BranchName, shwr_theta, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_thetaxz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_thetaxz, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_thetayz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_thetayz, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_phi_" + ShowerLabel;
  CreateBranch(BranchName, shwr_phi, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startdcosx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosx, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_startdcosy_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosy, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_startdcosz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosz, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_startx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startx, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_starty_" + ShowerLabel;
  CreateBranch(BranchName, shwr_starty, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_startz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startz, BranchName + NShowerIndexStr + "/F");
  
  BranchName = "shwr_totEng_" + ShowerLabel;
  CreateBranch(BranchName, shwr_totEng, BranchName + NShowerIndexStr + "[3]/F");
  
  BranchName = "shwr_dedx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_dedx, BranchName + NShowerIndexStr + "[3]/F");
  
  BranchName = "shwr_mipEng_" + ShowerLabel;
  CreateBranch(BranchName, shwr_mipEng, BranchName + NShowerIndexStr + "[3]/F");
  
  BranchName = "shwr_hasPFParticle_" + ShowerLabel;
  CreateBranch(BranchName, shwr_hasPFParticle, BranchName + NShowerIndexStr + "/S");
  
  BranchName = "shwr_PFParticleID_" + ShowerLabel;
  CreateBranch(BranchName, shwr_PFParticleID, BranchName + NShowerIndexStr + "/S");
  
} // microboone::AnalysisTreeDataStruct::ShowerDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct
//---

void microboone::AnalysisTreeDataStruct::ClearLocalData() {

  //  RunData.Clear();
  SubRunData.Clear();

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  beamtime = -99999;
  isdata = -99;
  taulife = -99999;
  triggernumber = 0;
  triggertime = -99999;
  beamgatetime = -99999;
  triggerbits = 0;
  potbnb = 0;
  potnumitgt = 0;
  potnumi101 = 0;

  opfilter_pebeam = -999;
  opfilter_peveto = -999;
  opfilter_pmtmaxfrac = -999;

  evtwgt_nfunc = 0;
  FillWith(evtwgt_funcname, "noname");
  FillWith(evtwgt_nweight, 0);

  no_hits = 0;
  no_hits_stored = 0;  
  
  std::fill(hit_plane, hit_plane + sizeof(hit_plane)/sizeof(hit_plane[0]), -9999);
  std::fill(hit_wire, hit_wire + sizeof(hit_wire)/sizeof(hit_wire[0]), -9999);
  std::fill(hit_channel, hit_channel + sizeof(hit_channel)/sizeof(hit_channel[0]), -9999);
  std::fill(hit_peakT, hit_peakT + sizeof(hit_peakT)/sizeof(hit_peakT[0]), -99999.);
  std::fill(hit_charge, hit_charge + sizeof(hit_charge)/sizeof(hit_charge[0]), -99999.);
  std::fill(hit_ph, hit_ph + sizeof(hit_ph)/sizeof(hit_ph[0]), -99999.);
  std::fill(hit_startT, hit_startT + sizeof(hit_startT)/sizeof(hit_startT[0]), -99999.);
  std::fill(hit_endT, hit_endT + sizeof(hit_endT)/sizeof(hit_endT[0]), -99999.);
  std::fill(hit_rms, hit_rms + sizeof(hit_rms)/sizeof(hit_rms[0]), -99999.);
  std::fill(hit_trueX, hit_trueX + sizeof(hit_trueX)/sizeof(hit_trueX[0]), -99999.);
  std::fill(hit_goodnessOfFit, hit_goodnessOfFit + sizeof(hit_goodnessOfFit)/sizeof(hit_goodnessOfFit[0]), -99999.);
  std::fill(hit_multiplicity, hit_multiplicity + sizeof(hit_multiplicity)/sizeof(hit_multiplicity[0]), -99999.);
  std::fill(hit_trkid, hit_trkid + sizeof(hit_trkid)/sizeof(hit_trkid[0]), -9999);
  std::fill(hit_trkKey, hit_trkKey + sizeof(hit_trkKey)/sizeof(hit_trkKey[0]), -9999);
  std::fill(hit_clusterid, hit_clusterid + sizeof(hit_clusterid)/sizeof(hit_clusterid[0]), -99999);
  std::fill(hit_clusterKey, hit_clusterKey + sizeof(hit_clusterKey)/sizeof(hit_clusterKey[0]), -9999);
  std::fill(hit_nelec, hit_nelec + sizeof(hit_nelec)/sizeof(hit_nelec[0]), -99999.);
  std::fill(hit_energy, hit_energy + sizeof(hit_energy)/sizeof(hit_energy[0]), -99999.);
  //raw digit information
  std::fill(rawD_ph, rawD_ph + sizeof(rawD_ph)/sizeof(rawD_ph[0]), -99999.);
  std::fill(rawD_peakT, rawD_peakT + sizeof(rawD_peakT)/sizeof(rawD_peakT[0]), -99999.);
  std::fill(rawD_charge, rawD_charge + sizeof(rawD_charge)/sizeof(rawD_charge[0]), -99999.);
  std::fill(rawD_fwhh, rawD_fwhh + sizeof(rawD_fwhh)/sizeof(rawD_fwhh[0]), -99999.);
  std::fill(rawD_rms, rawD_rms + sizeof(rawD_rms)/sizeof(rawD_rms[0]), -99999.);
  //cal wire information
  std::fill(calwire_ph, calwire_ph + sizeof(calwire_ph)/sizeof(calwire_ph[0]), -99999.);
  std::fill(calwire_peakT, calwire_peakT + sizeof(calwire_peakT)/sizeof(calwire_peakT[0]), -9999);
  std::fill(calwire_charge, calwire_charge + sizeof(calwire_charge)/sizeof(calwire_charge[0]), -99999.);
  std::fill(calwire_fwhh, calwire_fwhh + sizeof(calwire_fwhh)/sizeof(calwire_fwhh[0]), -99999.);
  std::fill(calwire_rms, calwire_rms + sizeof(calwire_rms)/sizeof(calwire_rms[0]), -99999.); 
  //Sim Channel information
  std::fill(sim_ph, sim_ph + sizeof(sim_ph)/sizeof(sim_ph[0]), -99999.);
  std::fill(sim_tdc, sim_tdc + sizeof(sim_tdc)/sizeof(sim_tdc[0]), -9999);
  std::fill(sim_charge, sim_charge + sizeof(sim_charge)/sizeof(sim_charge[0]), -99999.);
  std::fill(sim_fwhh, sim_fwhh + sizeof(sim_fwhh)/sizeof(sim_fwhh[0]), -99999.);
  std::fill(sim_rms, sim_rms + sizeof(sim_rms)/sizeof(sim_rms[0]), -99999.); 

/* 
  no_flashes = 0;
  std::fill(flash_time, flash_time + sizeof(flash_time)/sizeof(flash_time[0]), -9999);
  std::fill(flash_pe, flash_pe + sizeof(flash_pe)/sizeof(flash_pe[0]), -9999);
  std::fill(flash_ycenter, flash_ycenter + sizeof(flash_ycenter)/sizeof(flash_ycenter[0]), -9999);
  std::fill(flash_zcenter, flash_zcenter + sizeof(flash_zcenter)/sizeof(flash_zcenter[0]), -9999);
  std::fill(flash_ywidth, flash_ywidth + sizeof(flash_ywidth)/sizeof(flash_ywidth[0]), -9999);
  std::fill(flash_zwidth, flash_zwidth + sizeof(flash_zwidth)/sizeof(flash_zwidth[0]), -9999);
  std::fill(flash_timewidth, flash_timewidth + sizeof(flash_timewidth)/sizeof(flash_timewidth[0]), -9999);
*/

  nclusters = 0;
  std::fill(clusterId, clusterId + sizeof(clusterId)/sizeof(clusterId[0]), -9999);
  std::fill(clusterView, clusterView + sizeof(clusterView)/sizeof(clusterView[0]), -9999);
  std::fill(cluster_isValid, cluster_isValid + sizeof(cluster_isValid)/sizeof(cluster_isValid[0]), -1);
  std::fill(cluster_StartCharge, cluster_StartCharge +  sizeof(cluster_StartCharge)/sizeof(cluster_StartCharge[0]), -99999.);
  std::fill(cluster_StartAngle, cluster_StartAngle + sizeof(cluster_StartAngle)/sizeof(cluster_StartAngle[0]), -99999.);
  std::fill(cluster_EndCharge, cluster_EndCharge + sizeof(cluster_EndCharge)/sizeof(cluster_EndCharge[0]), -99999.);
  std::fill(cluster_EndAngle , cluster_EndAngle + sizeof(cluster_EndAngle)/sizeof(cluster_EndAngle[0]), -99999.);
  std::fill(cluster_Integral , cluster_Integral + sizeof(cluster_Integral)/sizeof(cluster_Integral[0]), -99999.);
  std::fill(cluster_IntegralAverage, cluster_IntegralAverage + sizeof(cluster_IntegralAverage)/sizeof(cluster_IntegralAverage[0]), -99999.);
  std::fill(cluster_SummedADC, cluster_SummedADC + sizeof(cluster_SummedADC)/sizeof(cluster_SummedADC[0]), -99999.);
  std::fill(cluster_SummedADCaverage, cluster_SummedADCaverage + sizeof(cluster_SummedADCaverage)/sizeof(cluster_SummedADCaverage[0]), -99999.);
  std::fill(cluster_MultipleHitDensity, cluster_MultipleHitDensity + sizeof(cluster_MultipleHitDensity)/sizeof(cluster_MultipleHitDensity[0]), -99999.);
  std::fill(cluster_Width, cluster_Width + sizeof(cluster_Width)/sizeof(cluster_Width[0]), -99999.);
  std::fill(cluster_NHits, cluster_NHits + sizeof(cluster_NHits)/sizeof(cluster_NHits[0]), -9999);
  std::fill(cluster_StartWire, cluster_StartWire + sizeof(cluster_StartWire)/sizeof(cluster_StartWire[0]), -9999);
  std::fill(cluster_StartTick, cluster_StartTick + sizeof(cluster_StartTick)/sizeof(cluster_StartTick[0]), -9999);
  std::fill(cluster_EndWire, cluster_EndWire + sizeof(cluster_EndWire)/sizeof(cluster_EndWire[0]), -9999);
  std::fill(cluster_EndTick, cluster_EndTick + sizeof(cluster_EndTick)/sizeof(cluster_EndTick[0]), -9999);
  std::fill(cluncosmictags_tagger, cluncosmictags_tagger + sizeof(cluncosmictags_tagger)/sizeof(cluncosmictags_tagger[0]), -9999);
  std::fill(clucosmicscore_tagger, clucosmicscore_tagger + sizeof(clucosmicscore_tagger)/sizeof(clucosmicscore_tagger[0]), -99999.);
  std::fill(clucosmictype_tagger , clucosmictype_tagger  + sizeof(clucosmictype_tagger )/sizeof(clucosmictype_tagger [0]), -9999);

  mcevts_truth = -99999;
  mcevts_truthcry = -99999;
  std::fill(nuPDG_truth, nuPDG_truth + sizeof(nuPDG_truth)/sizeof(nuPDG_truth[0]), -99999.);
  std::fill(ccnc_truth, ccnc_truth + sizeof(ccnc_truth)/sizeof(ccnc_truth[0]), -99999.);
  std::fill(mode_truth, mode_truth + sizeof(mode_truth)/sizeof(mode_truth[0]), -99999.);
  std::fill(enu_truth, enu_truth + sizeof(enu_truth)/sizeof(enu_truth[0]), -99999.);
  std::fill(Q2_truth, Q2_truth + sizeof(Q2_truth)/sizeof(Q2_truth[0]), -99999.);
  std::fill(W_truth, W_truth + sizeof(W_truth)/sizeof(W_truth[0]), -99999.);
  std::fill(X_truth, X_truth + sizeof(X_truth)/sizeof(X_truth[0]), -99999.);
  std::fill(Y_truth, Y_truth + sizeof(Y_truth)/sizeof(Y_truth[0]), -99999.);
  std::fill(hitnuc_truth, hitnuc_truth + sizeof(hitnuc_truth)/sizeof(hitnuc_truth[0]), -99999.);
  std::fill(nuvtxx_truth, nuvtxx_truth + sizeof(nuvtxx_truth)/sizeof(nuvtxx_truth[0]), -99999.);
  std::fill(nuvtxy_truth, nuvtxy_truth + sizeof(nuvtxy_truth)/sizeof(nuvtxy_truth[0]), -99999.);
  std::fill(nuvtxz_truth, nuvtxz_truth + sizeof(nuvtxz_truth)/sizeof(nuvtxz_truth[0]), -99999.);
  std::fill(sp_charge_corrected_nuvtxx_truth, sp_charge_corrected_nuvtxx_truth + sizeof(sp_charge_corrected_nuvtxx_truth)/sizeof(sp_charge_corrected_nuvtxx_truth[0]), -99999.); 
  std::fill(sp_charge_corrected_nuvtxy_truth, sp_charge_corrected_nuvtxy_truth + sizeof(sp_charge_corrected_nuvtxy_truth)/sizeof(sp_charge_corrected_nuvtxy_truth[0]), -99999.); 
  std::fill(sp_charge_corrected_nuvtxz_truth, sp_charge_corrected_nuvtxz_truth + sizeof(sp_charge_corrected_nuvtxz_truth)/sizeof(sp_charge_corrected_nuvtxz_truth[0]), -99999.); 
  std::fill(nu_dcosx_truth, nu_dcosx_truth + sizeof(nu_dcosx_truth)/sizeof(nu_dcosx_truth[0]), -99999.);
  std::fill(nu_dcosy_truth, nu_dcosy_truth + sizeof(nu_dcosy_truth)/sizeof(nu_dcosy_truth[0]), -99999.);
  std::fill(nu_dcosz_truth, nu_dcosz_truth + sizeof(nu_dcosz_truth)/sizeof(nu_dcosz_truth[0]), -99999.);
  std::fill(lep_mom_truth, lep_mom_truth + sizeof(lep_mom_truth)/sizeof(lep_mom_truth[0]), -99999.);
  std::fill(lep_dcosx_truth, lep_dcosx_truth + sizeof(lep_dcosx_truth)/sizeof(lep_dcosx_truth[0]), -99999.);
  std::fill(lep_dcosy_truth, lep_dcosy_truth + sizeof(lep_dcosy_truth)/sizeof(lep_dcosy_truth[0]), -99999.);
  std::fill(lep_dcosz_truth, lep_dcosz_truth + sizeof(lep_dcosz_truth)/sizeof(lep_dcosz_truth[0]), -99999.);

  //Flux information
  std::fill(vx_flux, vx_flux + sizeof(vx_flux)/sizeof(vx_flux[0]), -99999.);
  std::fill(vy_flux, vy_flux + sizeof(vy_flux)/sizeof(vy_flux[0]), -99999.);
  std::fill(vz_flux, vz_flux + sizeof(vz_flux)/sizeof(vz_flux[0]), -99999.);
  std::fill(pdpx_flux, pdpx_flux + sizeof(pdpx_flux)/sizeof(pdpx_flux[0]), -99999.);
  std::fill(pdpy_flux, pdpy_flux + sizeof(pdpy_flux)/sizeof(pdpy_flux[0]), -99999.);
  std::fill(pdpz_flux, pdpz_flux + sizeof(pdpz_flux)/sizeof(pdpz_flux[0]), -99999.);
  std::fill(ppdxdz_flux, ppdxdz_flux + sizeof(ppdxdz_flux)/sizeof(ppdxdz_flux[0]), -99999.);
  std::fill(ppdydz_flux, ppdydz_flux + sizeof(ppdydz_flux)/sizeof(ppdydz_flux[0]), -99999.);
  std::fill(pppz_flux, pppz_flux + sizeof(pppz_flux)/sizeof(pppz_flux[0]), -99999.);
  std::fill(ptype_flux, ptype_flux + sizeof(ptype_flux)/sizeof(ptype_flux[0]), -9999);
  std::fill(ppvx_flux, ppvx_flux + sizeof(ppvx_flux)/sizeof(ppvx_flux[0]), -99999.);
  std::fill(ppvy_flux, ppvy_flux + sizeof(ppvy_flux)/sizeof(ppvy_flux[0]), -99999.);
  std::fill(ppvz_flux, ppvz_flux + sizeof(ppvz_flux)/sizeof(ppvz_flux[0]), -99999.);
  std::fill(muparpx_flux, muparpx_flux + sizeof(muparpx_flux)/sizeof(muparpx_flux[0]), -99999.);
  std::fill(muparpy_flux, muparpy_flux + sizeof(muparpy_flux)/sizeof(muparpy_flux[0]), -99999.);
  std::fill(muparpz_flux, muparpz_flux + sizeof(muparpz_flux)/sizeof(muparpz_flux[0]), -99999.);
  std::fill(mupare_flux, mupare_flux + sizeof(mupare_flux)/sizeof(mupare_flux[0]), -99999.);
  std::fill(tgen_flux, tgen_flux + sizeof(tgen_flux)/sizeof(tgen_flux[0]), -9999);
  std::fill(tgptype_flux, tgptype_flux + sizeof(tgptype_flux)/sizeof(tgptype_flux[0]), -9999);
  std::fill(tgppx_flux, tgppx_flux + sizeof(tgppx_flux)/sizeof(tgppx_flux[0]), -99999.);
  std::fill(tgppy_flux, tgppy_flux + sizeof(tgppy_flux)/sizeof(tgppy_flux[0]), -99999.);
  std::fill(tgppz_flux, tgppz_flux + sizeof(tgppz_flux)/sizeof(tgppz_flux[0]), -99999.);
  std::fill(tprivx_flux, tprivx_flux + sizeof(tprivx_flux)/sizeof(tprivx_flux[0]), -99999.);
  std::fill(tprivy_flux, tprivy_flux + sizeof(tprivy_flux)/sizeof(tprivy_flux[0]), -99999.);
  std::fill(tprivz_flux, tprivz_flux + sizeof(tprivz_flux)/sizeof(tprivz_flux[0]), -99999.);
  std::fill(dk2gen_flux, dk2gen_flux + sizeof(dk2gen_flux)/sizeof(dk2gen_flux[0]), -99999.);
  std::fill(gen2vtx_flux, gen2vtx_flux + sizeof(gen2vtx_flux)/sizeof(gen2vtx_flux[0]), -99999.);
  std::fill(tpx_flux, tpx_flux + sizeof(tpx_flux)/sizeof(tpx_flux[0]), -99999.);
  std::fill(tpy_flux, tpy_flux + sizeof(tpy_flux)/sizeof(tpy_flux[0]), -99999.);
  std::fill(tpz_flux, tpz_flux + sizeof(tpz_flux)/sizeof(tpz_flux[0]), -99999.);
  std::fill(tptype_flux, tptype_flux + sizeof(tptype_flux)/sizeof(tptype_flux[0]), -99999.);

  genie_no_primaries = 0;
  cry_no_primaries = 0;
  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcAV = 0;
  no_mcshowers = 0;
  no_mctracks = 0;
  
  FillWith(pdg, -99999);
  FillWith(status, -99999);
  FillWith(Mass, -99999.);
  FillWith(Eng, -99999.);
  FillWith(EndE, -99999.);
  FillWith(Px, -99999.);
  FillWith(Py, -99999.);
  FillWith(Pz, -99999.);
  FillWith(P, -99999.);
  FillWith(StartPointx, -99999.); 
  FillWith(StartPointy, -99999.);
  FillWith(StartPointz, -99999.);
  FillWith(sp_charge_corrected_StartPointx,-99999.); 
  FillWith(sp_charge_corrected_StartPointy,-99999.); 
  FillWith(sp_charge_corrected_StartPointz,-99999.); 
  FillWith(StartT, -99999.);
  FillWith(EndT, -99999.);    
  FillWith(EndPointx, -99999.);
  FillWith(EndPointy, -99999.);
  FillWith(EndPointz, -99999.);
  FillWith(sp_charge_corrected_EndPointx,-99999.); 
  FillWith(sp_charge_corrected_EndPointy,-99999.); 
  FillWith(sp_charge_corrected_EndPointz,-99999.); 
  FillWith(EndT, -99999.);
  FillWith(theta, -99999.);
  FillWith(phi, -99999.);
  FillWith(theta_xz, -99999.);
  FillWith(theta_yz, -99999.);
  FillWith(pathlen, -99999.);
  FillWith(inTPCActive, -99999);
  FillWith(StartPointx_tpcAV, -99999.);
  FillWith(StartPointy_tpcAV, -99999.);
  FillWith(StartPointz_tpcAV, -99999.);
  FillWith(sp_charge_corrected_StartPointx_tpcAV, -99999.); 
  FillWith(sp_charge_corrected_StartPointy_tpcAV, -99999.); 
  FillWith(sp_charge_corrected_StartPointz_tpcAV, -99999.); 
  FillWith(StartT_tpcAV, -99999.);
  FillWith(StartE_tpcAV, -99999.);
  FillWith(StartP_tpcAV, -99999.);
  FillWith(StartPx_tpcAV, -99999.);
  FillWith(StartPy_tpcAV, -99999.);
  FillWith(StartPz_tpcAV, -99999.);
  FillWith(EndPointx_tpcAV, -99999.);
  FillWith(EndPointy_tpcAV, -99999.);
  FillWith(EndPointz_tpcAV, -99999.);
  FillWith(sp_charge_corrected_EndPointx_tpcAV, -99999.); 
  FillWith(sp_charge_corrected_EndPointy_tpcAV, -99999.); 
  FillWith(sp_charge_corrected_EndPointz_tpcAV, -99999.); 
  FillWith(EndT_tpcAV, -99999.);
  FillWith(EndE_tpcAV, -99999.);
  FillWith(EndP_tpcAV, -99999.);
  FillWith(EndPx_tpcAV, -99999.);
  FillWith(EndPy_tpcAV, -99999.);
  FillWith(EndPz_tpcAV, -99999.);
  FillWith(pathlen_drifted, -99999.);
  FillWith(inTPCDrifted, -99999);
  FillWith(StartPointx_drifted, -99999.);
  FillWith(StartPointy_drifted, -99999.);
  FillWith(StartPointz_drifted, -99999.);
  FillWith(sp_charge_corrected_StartPointx_drifted, -99999.);
  FillWith(sp_charge_corrected_StartPointy_drifted, -99999.);
  FillWith(sp_charge_corrected_StartPointz_drifted, -99999.);
  FillWith(StartT_drifted, -99999.);
  FillWith(StartE_drifted, -99999.);
  FillWith(StartP_drifted, -99999.);
  FillWith(StartPx_drifted, -99999.);
  FillWith(StartPy_drifted, -99999.);
  FillWith(StartPz_drifted, -99999.);
  FillWith(EndPointx_drifted, -99999.);
  FillWith(EndPointy_drifted, -99999.);
  FillWith(EndPointz_drifted, -99999.); 
  FillWith(sp_charge_corrected_EndPointx_drifted, -99999.);
  FillWith(sp_charge_corrected_EndPointy_drifted, -99999.);
  FillWith(sp_charge_corrected_EndPointz_drifted, -99999.);
  FillWith(EndT_drifted, -99999.);
  FillWith(EndE_drifted, -99999.);
  FillWith(EndP_drifted, -99999.); 
  FillWith(EndPx_drifted, -99999.); 
  FillWith(EndPy_drifted, -99999.); 
  FillWith(EndPz_drifted, -99999.); 
  FillWith(NumberDaughters, -99999);
  FillWith(Mother, -99999);
  FillWith(TrackId, -99999);
  FillWith(process_primary, -99999);
  FillWith(processname, "noname");
  FillWith(MergedId, -99999);
  FillWith(origin, -99999);
  FillWith(MCTruthIndex, -99999);
  FillWith(genie_primaries_pdg, -99999);
  FillWith(genie_Eng, -99999.);
  FillWith(genie_Px, -99999.);
  FillWith(genie_Py, -99999.);
  FillWith(genie_Pz, -99999.);
  FillWith(genie_P, -99999.);
  FillWith(genie_status_code, -99999);
  FillWith(genie_mass, -99999.);
  FillWith(genie_trackID, -99999);
  FillWith(genie_ND, -99999);
  FillWith(genie_mother, -99999);
  FillWith(cry_primaries_pdg, -99999);
  FillWith(cry_Eng, -99999.);
  FillWith(cry_Px, -99999.);
  FillWith(cry_Py, -99999.);
  FillWith(cry_Pz, -99999.);
  FillWith(cry_P, -99999.);
  FillWith(cry_StartPointx, -99999.);
  FillWith(cry_StartPointy, -99999.);
  FillWith(cry_StartPointz, -99999.);
  FillWith(cry_StartPointt, -99999.);  
  FillWith(cry_status_code, -99999);
  FillWith(cry_mass, -99999.);
  FillWith(cry_trackID, -99999);
  FillWith(cry_ND, -99999);
  FillWith(cry_mother, -99999); 
  FillWith(mcshwr_origin, -1);   
  FillWith(mcshwr_pdg, -99999);
  FillWith(mcshwr_TrackId, -99999);
  FillWith(mcshwr_Process, "noname");
  FillWith(mcshwr_startX, -99999.);
  FillWith(mcshwr_startY, -99999.);
  FillWith(mcshwr_startZ, -99999.);
  FillWith(mcshwr_endX, -99999.);
  FillWith(mcshwr_endY, -99999.);
  FillWith(mcshwr_endZ, -99999.);
  FillWith(mcshwr_CombEngX, -99999.);
  FillWith(mcshwr_CombEngY, -99999.);
  FillWith(mcshwr_CombEngZ, -99999.);
  FillWith(mcshwr_CombEngPx, -99999.);
  FillWith(mcshwr_CombEngPy, -99999.);
  FillWith(mcshwr_CombEngPz, -99999.);
  FillWith(mcshwr_CombEngE, -99999.);
  FillWith(mcshwr_dEdx, -99999.);
  FillWith(mcshwr_StartDirX, -99999.);
  FillWith(mcshwr_StartDirY, -99999.);
  FillWith(mcshwr_StartDirZ, -99999.);
  FillWith(mcshwr_isEngDeposited, -9999);  
  FillWith(mcshwr_Motherpdg, -99999);
  FillWith(mcshwr_MotherTrkId, -99999);
  FillWith(mcshwr_MotherProcess, "noname");
  FillWith(mcshwr_MotherstartX, -99999.);
  FillWith(mcshwr_MotherstartY, -99999.);
  FillWith(mcshwr_MotherstartZ, -99999.);
  FillWith(mcshwr_MotherendX, -99999.);
  FillWith(mcshwr_MotherendY, -99999.);
  FillWith(mcshwr_MotherendZ, -99999.);
  FillWith(mcshwr_Ancestorpdg, -99999);
  FillWith(mcshwr_AncestorTrkId, -99999);  
  FillWith(mcshwr_AncestorProcess, "noname");
  FillWith(mcshwr_AncestorstartX, -99999.);
  FillWith(mcshwr_AncestorstartY, -99999.);
  FillWith(mcshwr_AncestorstartZ, -99999.);
  FillWith(mcshwr_AncestorendX, -99999.);
  FillWith(mcshwr_AncestorendY, -99999.);
  FillWith(mcshwr_AncestorendZ, -99999.);
  
  // auxiliary detector information;
  FillWith(NAuxDets, 0);
  // - set to -9999 all the values of each of the arrays in AuxDetID;
  //   this auto is BoxedArray<Short_t>
  for (auto& partInfo: AuxDetID) FillWith(partInfo, -9999);
  // - pythonish C++: as the previous line, for each one in a list of containers
  //   of the same type (C++ is not python yet), using pointers to avoid copy;
  for (AuxDetMCData_t<Float_t>* cont: {
      &entryX, &entryY, &entryZ, &entryT,
	&exitX , &exitY , &exitZ, &exitT, &exitPx, &exitPy, &exitPz,
	&CombinedEnergyDep
	})
    {
      // this auto is BoxedArray<Float_t>
      for (auto& partInfo: *cont) FillWith(partInfo, -99999.);
    } // for container
  
} // microboone::AnalysisTreeDataStruct::ClearLocalData()


void microboone::AnalysisTreeDataStruct::Clear() {
  ClearLocalData();
  std::for_each
    (TrackData.begin(), TrackData.end(), std::mem_fn(&TrackDataStruct::Clear));
  std::for_each
    (VertexData.begin(), VertexData.end(), std::mem_fn(&VertexDataStruct::Clear));
  std::for_each
    (NeutrinoVertexData.begin(), NeutrinoVertexData.end(), std::mem_fn(&NeutrinoVertexDataStruct::Clear));
  std::for_each
    (FlashData.begin(), FlashData.end(), std::mem_fn(&FlashDataStruct::Clear));
  std::for_each
    (ShowerData.begin(), ShowerData.end(), std::mem_fn(&ShowerDataStruct::Clear));
} // microboone::AnalysisTreeDataStruct::Clear()


void microboone::AnalysisTreeDataStruct::SetShowerAlgos
(std::vector<std::string> const& ShowerAlgos)
{
  
  size_t const nShowerAlgos = ShowerAlgos.size();
  ShowerData.resize(nShowerAlgos);
  for (size_t iAlgo = 0; iAlgo < nShowerAlgos; ++iAlgo)
    ShowerData[iAlgo].SetName(ShowerAlgos[iAlgo]);
  
} // microboone::AnalysisTreeDataStruct::SetShowerAlgos()


void microboone::AnalysisTreeDataStruct::ResizeGEANT(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxGEANTparticles = (size_t) std::max(nParticles, 1);
  
  pdg.resize(MaxGEANTparticles);
  status.resize(MaxGEANTparticles);
  Mass.resize(MaxGEANTparticles);  
  Eng.resize(MaxGEANTparticles);
  EndE.resize(MaxGEANTparticles);
  Px.resize(MaxGEANTparticles);
  Py.resize(MaxGEANTparticles);
  Pz.resize(MaxGEANTparticles);
  P.resize(MaxGEANTparticles);
  StartPointx.resize(MaxGEANTparticles);
  StartPointy.resize(MaxGEANTparticles);
  StartPointz.resize(MaxGEANTparticles);
  sp_charge_corrected_StartPointx.resize(MaxGEANTparticles); 
  sp_charge_corrected_StartPointy.resize(MaxGEANTparticles); 
  sp_charge_corrected_StartPointz.resize(MaxGEANTparticles); 
  StartT.resize(MaxGEANTparticles); 
  EndT.resize(MaxGEANTparticles);    
  EndPointx.resize(MaxGEANTparticles);
  EndPointy.resize(MaxGEANTparticles);
  EndPointz.resize(MaxGEANTparticles);
  sp_charge_corrected_EndPointx.resize(MaxGEANTparticles); 
  sp_charge_corrected_EndPointy.resize(MaxGEANTparticles); 
  sp_charge_corrected_EndPointz.resize(MaxGEANTparticles); 
  EndT.resize(MaxGEANTparticles);  
  theta.resize(MaxGEANTparticles);
  phi.resize(MaxGEANTparticles);
  theta_xz.resize(MaxGEANTparticles);
  theta_yz.resize(MaxGEANTparticles);
  pathlen.resize(MaxGEANTparticles);
  inTPCActive.resize(MaxGEANTparticles);
  StartPointx_tpcAV.resize(MaxGEANTparticles);
  StartPointy_tpcAV.resize(MaxGEANTparticles);
  StartPointz_tpcAV.resize(MaxGEANTparticles);
  sp_charge_corrected_StartPointx_tpcAV.resize(MaxGEANTparticles); 
  sp_charge_corrected_StartPointy_tpcAV.resize(MaxGEANTparticles); 
  sp_charge_corrected_StartPointz_tpcAV.resize(MaxGEANTparticles); 
  StartT_tpcAV.resize(MaxGEANTparticles);
  StartE_tpcAV.resize(MaxGEANTparticles);
  StartP_tpcAV.resize(MaxGEANTparticles);
  StartPx_tpcAV.resize(MaxGEANTparticles);
  StartPy_tpcAV.resize(MaxGEANTparticles);
  StartPz_tpcAV.resize(MaxGEANTparticles);
  EndPointx_tpcAV.resize(MaxGEANTparticles);
  EndPointy_tpcAV.resize(MaxGEANTparticles);
  EndPointz_tpcAV.resize(MaxGEANTparticles); 
  sp_charge_corrected_EndPointx_tpcAV.resize(MaxGEANTparticles); 
  sp_charge_corrected_EndPointy_tpcAV.resize(MaxGEANTparticles); 
  sp_charge_corrected_EndPointz_tpcAV.resize(MaxGEANTparticles); 
  EndT_tpcAV.resize(MaxGEANTparticles);
  EndE_tpcAV.resize(MaxGEANTparticles);
  EndP_tpcAV.resize(MaxGEANTparticles);
  EndPx_tpcAV.resize(MaxGEANTparticles);
  EndPy_tpcAV.resize(MaxGEANTparticles);
  EndPz_tpcAV.resize(MaxGEANTparticles);
  pathlen_drifted.resize(MaxGEANTparticles);
  inTPCDrifted.resize(MaxGEANTparticles);
  StartPointx_drifted.resize(MaxGEANTparticles);
  StartPointy_drifted.resize(MaxGEANTparticles);
  StartPointz_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_StartPointx_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_StartPointy_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_StartPointz_drifted.resize(MaxGEANTparticles);
  StartT_drifted.resize(MaxGEANTparticles);
  StartE_drifted.resize(MaxGEANTparticles);
  StartP_drifted.resize(MaxGEANTparticles);
  StartPx_drifted.resize(MaxGEANTparticles);
  StartPy_drifted.resize(MaxGEANTparticles);
  StartPz_drifted.resize(MaxGEANTparticles);
  EndPointx_drifted.resize(MaxGEANTparticles);
  EndPointy_drifted.resize(MaxGEANTparticles);
  EndPointz_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_EndPointx_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_EndPointy_drifted.resize(MaxGEANTparticles);
  sp_charge_corrected_EndPointz_drifted.resize(MaxGEANTparticles);
  EndT_drifted.resize(MaxGEANTparticles);   
  EndE_drifted.resize(MaxGEANTparticles);   
  EndP_drifted.resize(MaxGEANTparticles);   
  EndPx_drifted.resize(MaxGEANTparticles);   
  EndPy_drifted.resize(MaxGEANTparticles);   
  EndPz_drifted.resize(MaxGEANTparticles);   
  NumberDaughters.resize(MaxGEANTparticles);
  Mother.resize(MaxGEANTparticles);
  TrackId.resize(MaxGEANTparticles);
  process_primary.resize(MaxGEANTparticles);
  processname.resize(MaxGEANTparticles);
  MergedId.resize(MaxGEANTparticles);
  origin.resize(MaxGEANTparticles);
  MCTruthIndex.resize(MaxGEANTparticles);  

  // auxiliary detector structure
  NAuxDets.resize(MaxGEANTparticles);
  AuxDetID.resize(MaxGEANTparticles);
  entryX.resize(MaxGEANTparticles);
  entryY.resize(MaxGEANTparticles);
  entryZ.resize(MaxGEANTparticles);
  entryT.resize(MaxGEANTparticles);
  exitX.resize(MaxGEANTparticles);
  exitY.resize(MaxGEANTparticles);
  exitZ.resize(MaxGEANTparticles);
  exitT.resize(MaxGEANTparticles);
  exitPx.resize(MaxGEANTparticles);
  exitPy.resize(MaxGEANTparticles);
  exitPz.resize(MaxGEANTparticles);
  CombinedEnergyDep.resize(MaxGEANTparticles);
  
} // microboone::AnalysisTreeDataStruct::ResizeGEANT()

void microboone::AnalysisTreeDataStruct::ResizeGenie(int nPrimaries) {
  
  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);
  genie_primaries_pdg.resize(MaxGeniePrimaries);
  genie_Eng.resize(MaxGeniePrimaries);
  genie_Px.resize(MaxGeniePrimaries);
  genie_Py.resize(MaxGeniePrimaries);
  genie_Pz.resize(MaxGeniePrimaries);
  genie_P.resize(MaxGeniePrimaries);
  genie_status_code.resize(MaxGeniePrimaries);
  genie_mass.resize(MaxGeniePrimaries);
  genie_trackID.resize(MaxGeniePrimaries);
  genie_ND.resize(MaxGeniePrimaries);
  genie_mother.resize(MaxGeniePrimaries);
} // microboone::AnalysisTreeDataStruct::ResizeGenie()

void microboone::AnalysisTreeDataStruct::ResizeCry(int nPrimaries) {

  cry_primaries_pdg.resize(nPrimaries);
  cry_Eng.resize(nPrimaries);
  cry_Px.resize(nPrimaries);
  cry_Py.resize(nPrimaries);
  cry_Pz.resize(nPrimaries);
  cry_P.resize(nPrimaries);
  cry_StartPointx.resize(nPrimaries);
  cry_StartPointy.resize(nPrimaries);
  cry_StartPointz.resize(nPrimaries); 
  cry_StartPointt.resize(nPrimaries);  
  cry_status_code.resize(nPrimaries);
  cry_mass.resize(nPrimaries);
  cry_trackID.resize(nPrimaries);
  cry_ND.resize(nPrimaries);
  cry_mother.resize(nPrimaries);

} // microboone::AnalysisTreeDataStruct::ResizeCry()

void microboone::AnalysisTreeDataStruct::ResizeMCShower(int nMCShowers) {
  mcshwr_origin.resize(nMCShowers);		
  mcshwr_pdg.resize(nMCShowers);		
  mcshwr_TrackId.resize(nMCShowers); 	
  mcshwr_Process.resize(nMCShowers); 
  mcshwr_startX.resize(nMCShowers); 	
  mcshwr_startY.resize(nMCShowers);	 
  mcshwr_startZ.resize(nMCShowers);	 
  mcshwr_endX.resize(nMCShowers);	 
  mcshwr_endY.resize(nMCShowers);	 
  mcshwr_endZ.resize(nMCShowers);	 
  mcshwr_CombEngX.resize(nMCShowers);	 
  mcshwr_CombEngY.resize(nMCShowers);	 
  mcshwr_CombEngZ.resize(nMCShowers);	 
  mcshwr_CombEngPx.resize(nMCShowers);	 
  mcshwr_CombEngPy.resize(nMCShowers);	 
  mcshwr_CombEngPz.resize(nMCShowers);	 
  mcshwr_CombEngE.resize(nMCShowers);
  mcshwr_dEdx.resize(nMCShowers);
  mcshwr_StartDirX.resize(nMCShowers);
  mcshwr_StartDirY.resize(nMCShowers);
  mcshwr_StartDirZ.resize(nMCShowers);
  mcshwr_isEngDeposited.resize(nMCShowers);	   	
  mcshwr_Motherpdg.resize(nMCShowers);	
  mcshwr_MotherTrkId.resize(nMCShowers);	
  mcshwr_MotherProcess.resize(nMCShowers);
  mcshwr_MotherstartX.resize(nMCShowers);	 
  mcshwr_MotherstartY.resize(nMCShowers);	 
  mcshwr_MotherstartZ.resize(nMCShowers);	 
  mcshwr_MotherendX.resize(nMCShowers);	 
  mcshwr_MotherendY.resize(nMCShowers);	 
  mcshwr_MotherendZ.resize(nMCShowers);	   	
  mcshwr_Ancestorpdg.resize(nMCShowers);	
  mcshwr_AncestorTrkId.resize(nMCShowers);	
  mcshwr_AncestorProcess.resize(nMCShowers); 
  mcshwr_AncestorstartX.resize(nMCShowers);	 
  mcshwr_AncestorstartY.resize(nMCShowers);	 
  mcshwr_AncestorstartZ.resize(nMCShowers);	 
  mcshwr_AncestorendX.resize(nMCShowers);	 
  mcshwr_AncestorendY.resize(nMCShowers);	 
  mcshwr_AncestorendZ.resize(nMCShowers);	 	

} // microboone::AnalysisTreeDataStruct::ResizeMCShower()

void microboone::AnalysisTreeDataStruct::ResizeMCTrack(int nMCTracks) {
  mctrk_origin.resize(nMCTracks);		
  mctrk_pdg.resize(nMCTracks);		
  mctrk_TrackId.resize(nMCTracks); 	
  mctrk_Process.resize(nMCTracks); 
  mctrk_startX.resize(nMCTracks); 	
  mctrk_startY.resize(nMCTracks);	 
  mctrk_startZ.resize(nMCTracks);	 
  mctrk_endX.resize(nMCTracks);	 
  mctrk_endY.resize(nMCTracks);	 
  mctrk_endZ.resize(nMCTracks);
  mctrk_startX_drifted.resize(nMCTracks);
  mctrk_startY_drifted.resize(nMCTracks);
  mctrk_startZ_drifted.resize(nMCTracks);
  mctrk_endX_drifted.resize(nMCTracks);  
  mctrk_endY_drifted.resize(nMCTracks);  
  mctrk_endZ_drifted.resize(nMCTracks);  
  mctrk_len_drifted.resize(nMCTracks);
  mctrk_p_drifted.resize(nMCTracks);
  mctrk_px_drifted.resize(nMCTracks);
  mctrk_py_drifted.resize(nMCTracks);
  mctrk_pz_drifted.resize(nMCTracks);
  mctrk_Motherpdg.resize(nMCTracks);	
  mctrk_MotherTrkId.resize(nMCTracks);	
  mctrk_MotherProcess.resize(nMCTracks);
  mctrk_MotherstartX.resize(nMCTracks);	 
  mctrk_MotherstartY.resize(nMCTracks);	 
  mctrk_MotherstartZ.resize(nMCTracks);	 
  mctrk_MotherendX.resize(nMCTracks);	 
  mctrk_MotherendY.resize(nMCTracks);	 
  mctrk_MotherendZ.resize(nMCTracks);	   	
  mctrk_Ancestorpdg.resize(nMCTracks);	
  mctrk_AncestorTrkId.resize(nMCTracks);	
  mctrk_AncestorProcess.resize(nMCTracks); 
  mctrk_AncestorstartX.resize(nMCTracks);	 
  mctrk_AncestorstartY.resize(nMCTracks);	 
  mctrk_AncestorstartZ.resize(nMCTracks);	 
  mctrk_AncestorendX.resize(nMCTracks);	 
  mctrk_AncestorendY.resize(nMCTracks);	 
  mctrk_AncestorendZ.resize(nMCTracks);	 
  
} // microboone::AnalysisTreeDataStruct::ResizeMCTrack()
  


void microboone::AnalysisTreeDataStruct::SetAddresses(
						      TTree* pTree,
						      const std::vector<std::string>& trackers,
						      const std::vector<std::string>& vertexalgos,
                                                      const std::vector<std::string>& nuvertexalgos,
                                                      const std::vector<std::string>& flashalgos,
						      const std::vector<std::string>& showeralgos,
						      bool isCosmics
						      ) {
  BranchCreator CreateBranch(pTree);

  CreateBranch("run",&run,"run/I");
  CreateBranch("subrun",&subrun,"subrun/I");
  CreateBranch("event",&event,"event/I");
  CreateBranch("evttime",&evttime,"evttime/D");
  CreateBranch("beamtime",&beamtime,"beamtime/D");
  CreateBranch("pot",&SubRunData.pot,"pot/D");
  CreateBranch("isdata",&isdata,"isdata/B");
  CreateBranch("taulife",&taulife,"taulife/D");
  CreateBranch("triggernumber",&triggernumber,"triggernumber/i");
  CreateBranch("triggertime",&triggertime,"triggertime/D");
  CreateBranch("beamgatetime",&beamgatetime,"beamgatetime/D");
  CreateBranch("triggerbits",&triggerbits,"triggerbits/i");
  CreateBranch("potbnb",&potbnb,"potbnb/D");
  CreateBranch("potnumitgt",&potnumitgt,"potnumitgt/D");
  CreateBranch("potnumi101",&potnumi101,"potnumi101/D");

  CreateBranch("evtwgt_funcname",evtwgt_funcname);
  CreateBranch("evtwgt_weight",evtwgt_weight);
  CreateBranch("evtwgt_nweight",evtwgt_nweight);
  CreateBranch("evtwgt_nfunc",&evtwgt_nfunc,"evtwgt_nfunc/I");

  if(hasOpticalFilterInfo()){
    CreateBranch("opfilter_pebeam",&opfilter_pebeam,"opfilter_pebeam/F");
    CreateBranch("opfilter_peveto",&opfilter_peveto,"opfilter_peveto/F");
    CreateBranch("opfilter_pmtmaxfrac",&opfilter_pmtmaxfrac,"opfilter_pmtmaxfrac/F");
  }

  if (hasSWTriggerInfo()){
    CreateBranch("swtrigger_name",      swtrigger_name);
    CreateBranch("swtrigger_triggered", swtrigger_triggered);
  }

  if (hasHitInfo()){    
    CreateBranch("no_hits",&no_hits,"no_hits/I");
    CreateBranch("no_hits_stored",&no_hits_stored,"no_hits_stored/I");    
    CreateBranch("hit_plane",hit_plane,"hit_plane[no_hits_stored]/S");
    CreateBranch("hit_wire",hit_wire,"hit_wire[no_hits_stored]/S");
    CreateBranch("hit_channel",hit_channel,"hit_channel[no_hits_stored]/S");
    CreateBranch("hit_peakT",hit_peakT,"hit_peakT[no_hits_stored]/F");
    CreateBranch("hit_charge",hit_charge,"hit_charge[no_hits_stored]/F");
    CreateBranch("hit_ph",hit_ph,"hit_ph[no_hits_stored]/F");
    CreateBranch("hit_startT",hit_startT,"hit_startT[no_hits_stored]/F");
    CreateBranch("hit_endT",hit_endT,"hit_endT[no_hits_stored]/F");
    CreateBranch("hit_rms",hit_rms,"hit_rms[no_hits_stored]/F"); 
    CreateBranch("hit_trueX",hit_trueX,"hit_trueX[no_hits_stored]/F");    
    CreateBranch("hit_goodnessOfFit",hit_goodnessOfFit,"hit_goodnessOfFit[no_hits_stored]/F");    
    CreateBranch("hit_multiplicity",hit_multiplicity,"hit_multiplicity[no_hits_stored]/S");    
    CreateBranch("hit_trkid",hit_trkid,"hit_trkid[no_hits_stored]/S");
    CreateBranch("hit_trkKey",hit_trkKey,"hit_trkKey[no_hits_stored]/S");   
    CreateBranch("hit_clusterid",hit_clusterid,"hit_clusterid[no_hits_stored]/S");
    CreateBranch("hit_clusterKey",hit_clusterKey,"hit_clusterKey[no_hits_stored]/S");
    if (!isCosmics){
      CreateBranch("hit_nelec",hit_nelec,"hit_nelec[no_hits_stored]/F");
      CreateBranch("hit_energy",hit_energy,"hit_energy[no_hits_stored]/F");
    }
    if (hasRawDigitInfo()){
      CreateBranch("rawD_ph",rawD_ph,"rawD_ph[no_hits_stored]/F");  
      CreateBranch("rawD_peakT",rawD_peakT,"rawD_peakT[no_hits_stored]/F");
      CreateBranch("rawD_charge",rawD_charge,"rawD_charge[no_hits_stored]/F");
      CreateBranch("rawD_fwhh",rawD_fwhh,"rawD_fwhh[no_hits_stored]/F");
      CreateBranch("rawD_rms",rawD_rms,"rawD_rms[no_hits_stored]/D"); 
    } 
    if (hasCalWireInfo()){
      CreateBranch("calwire_ph",calwire_ph,"calwire_ph[no_hits_stored]/F");  
      CreateBranch("calwire_peakT",calwire_peakT,"calwire_peakT[no_hits_stored]/I");
      CreateBranch("calwire_charge",calwire_charge,"calwire_charge[no_hits_stored]/F");
      CreateBranch("calwire_fwhh",calwire_fwhh,"calwire_fwhh[no_hits_stored]/F");
      CreateBranch("calwire_rms",calwire_rms,"calwire_rms[no_hits_stored]/D"); 
    }	
    if (hasSimChannelInfo()){
      CreateBranch("sim_ph",sim_ph,"sim_ph[no_hits_stored]/F");  
      CreateBranch("sim_tdc",sim_tdc,"sim_tdc[no_hits_stored]/I");
      CreateBranch("sim_charge",sim_charge,"sim_charge[no_hits_stored]/F");
      CreateBranch("sim_fwhh",sim_fwhh,"sim_fwhh[no_hits_stored]/F");
      CreateBranch("sim_rms",sim_rms,"sim_rms[no_hits_stored]/D"); 
    }	   
  }

  if (hasClusterInfo()){
    CreateBranch("nclusters",&nclusters,"nclusters/S");
    CreateBranch("clusterId", clusterId, "clusterId[nclusters]/S");
    CreateBranch("clusterView", clusterView, "clusterView[nclusters]/S");
    CreateBranch("cluster_StartCharge", cluster_StartCharge, "cluster_StartCharge[nclusters]/F");
    CreateBranch("cluster_StartAngle", cluster_StartAngle, "cluster_StartAngle[nclusters]/F");
    CreateBranch("cluster_EndCharge", cluster_EndCharge, "cluster_EndCharge[nclusters]/F");
    CreateBranch("cluster_EndAngle", cluster_EndAngle, "cluster_EndAngle[nclusters]/F");
    CreateBranch("cluster_Integral", cluster_Integral, "cluster_Integral[nclusters]/F");
    CreateBranch("cluster_IntegralAverage", cluster_IntegralAverage, "cluster_IntegralAverage[nclusters]/F");
    CreateBranch("cluster_SummedADC", cluster_SummedADC, "cluster_SummedADC[nclusters]/F");
    CreateBranch("cluster_SummedADCaverage", cluster_SummedADCaverage, "cluster_SummedADCaverage[nclusters]/F");
    CreateBranch("cluster_MultipleHitDensity", cluster_MultipleHitDensity, "cluster_MultipleHitDensity[nclusters]/F");
    CreateBranch("cluster_Width", cluster_Width, "cluster_Width[nclusters]/F");
    CreateBranch("cluster_NHits", cluster_NHits, "cluster_NHits[nclusters]/S");
    CreateBranch("cluster_StartWire", cluster_StartWire, "cluster_StartWire[nclusters]/S");
    CreateBranch("cluster_StartTick", cluster_StartTick, "cluster_StartTick[nclusters]/S");
    CreateBranch("cluster_EndWire", cluster_EndWire, "cluster_EndWire[nclusters]/S");
    CreateBranch("cluster_EndTick", cluster_EndTick, "cluster_EndTick[nclusters]/S");   
    CreateBranch("cluncosmictags_tagger", cluncosmictags_tagger, "cluncosmictags_tagger[nclusters]/S");   
    CreateBranch("clucosmicscore_tagger", clucosmicscore_tagger, "clucosmicscore_tagger[nclusters]/F");	 
    CreateBranch("clucosmictype_tagger", clucosmictype_tagger, "clucosmictype_tagger[nclusters]/S");	 
  }    

  if (hasFlashInfo()){
    kNFlashAlgos = flashalgos.size();
    CreateBranch("kNFlashAlgos",&kNFlashAlgos,"kNFlashAlgos/B");
    for(int i=0; i<kNFlashAlgos; i++) {
      std::string FlashLabel = flashalgos[i];

      // note that if the flashes data has maximum number of flashes 0,
      // nothing is initialized (branches are not even created)
      FlashData[i].SetAddresses(pTree, FlashLabel, isCosmics);
    }
  } // has flash info
  
  if (hasTrackInfo()){
    kNTracker = trackers.size();
    CreateBranch("kNTracker",&kNTracker,"kNTracker/B");
    for(int i=0; i<kNTracker; i++){
      std::string TrackLabel = trackers[i];

      // note that if the tracker data has maximum number of tracks 0,
      // nothing is initialized (branches are not even created)
      TrackData[i].SetAddresses(pTree, TrackLabel, isCosmics);
    } // for trackers
  } 

  if (hasVertexInfo()){
    kNVertexAlgos = vertexalgos.size();
    CreateBranch("kNVertexAlgos",&kNVertexAlgos,"kNVertexAlgos/B");
    for(int i=0; i<kNVertexAlgos; i++){
      std::string VertexLabel = vertexalgos[i];

      // note that if the vertes data has maximum number of vertices 0,
      // nothing is initialized (branches are not even created)
      VertexData[i].SetAddresses(pTree, VertexLabel, isCosmics);
    } // for trackers
  }

  if (hasPandoraNuVertexInfo()){
    kNNeutrinoVertexAlgos = nuvertexalgos.size();
    CreateBranch("kNNeutrinoVertexAlgos",&kNNeutrinoVertexAlgos,"kNNeutrinoVertexAlgos/B");
    for(int i=0; i<kNNeutrinoVertexAlgos; i++){
      std::string NuVertexLabel = nuvertexalgos[i];

      // note that if the vertes data has maximum number of vertices 0,
      // nothing is initialized (branches are not even created)
      NeutrinoVertexData[i].SetAddresses(pTree, NuVertexLabel, isCosmics);
    } // for trackers
  }

  if (hasShowerInfo()){
    kNShowerAlgos = showeralgos.size();
    CreateBranch("kNShowerAlgos",&kNShowerAlgos,"kNShowerAlgos/B");
    for(int i=0; i<kNShowerAlgos; i++){
      // note that if the shower data has maximum number of showers 0,
      // nothing is initialized (branches are not even created)
      ShowerData[i].SetAddresses(pTree);
    } // for showers
  } // if we have shower algos

  if (hasPFParticleInfo()){
    //CreateBranch("kNVertexAlgos",&kNVertexAlgos,"kNVertexAlgos/B"); // What would be the PFParticle equivalent of this? There's only 1 algo!
    PFParticleData.SetAddresses(pTree);
  }
  
  if (hasRawWaveFormInfo()){
    RawData.SetAddresses(pTree);        	  
  }
  
  if (hasCalibWaveFormInfo()){
    CalibData.SetAddresses(pTree);        	  
  }

  if (hasGenieInfo()){
    CreateBranch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    CreateBranch("nuPDG_truth",nuPDG_truth,"nuPDG_truth[mcevts_truth]/I");
    CreateBranch("ccnc_truth",ccnc_truth,"ccnc_truth[mcevts_truth]/I");
    CreateBranch("mode_truth",mode_truth,"mode_truth[mcevts_truth]/I");
    CreateBranch("enu_truth",enu_truth,"enu_truth[mcevts_truth]/F");
    CreateBranch("Q2_truth",Q2_truth,"Q2_truth[mcevts_truth]/F");
    CreateBranch("W_truth",W_truth,"W_truth[mcevts_truth]/F");
    CreateBranch("X_truth",X_truth,"X_truth[mcevts_truth]/F");
    CreateBranch("Y_truth",Y_truth,"Y_truth[mcevts_truth]/F");
    CreateBranch("hitnuc_truth",hitnuc_truth,"hitnuc_truth[mcevts_truth]/I");
    CreateBranch("nuvtxx_truth",nuvtxx_truth,"nuvtxx_truth[mcevts_truth]/F");
    CreateBranch("nuvtxy_truth",nuvtxy_truth,"nuvtxy_truth[mcevts_truth]/F");
    CreateBranch("nuvtxz_truth",nuvtxz_truth,"nuvtxz_truth[mcevts_truth]/F");
    CreateBranch("sp_charge_corrected_nuvtxx_truth",sp_charge_corrected_nuvtxx_truth,"sp_charge_corrected_nuvtxx_truth[mcevts_truth]/F"); 
    CreateBranch("sp_charge_corrected_nuvtxy_truth",sp_charge_corrected_nuvtxy_truth,"sp_charge_corrected_nuvtxy_truth[mcevts_truth]/F"); 
    CreateBranch("sp_charge_corrected_nuvtxz_truth",sp_charge_corrected_nuvtxz_truth,"sp_charge_corrected_nuvtxz_truth[mcevts_truth]/F"); 
    CreateBranch("nu_dcosx_truth",nu_dcosx_truth,"nu_dcosx_truth[mcevts_truth]/F");
    CreateBranch("nu_dcosy_truth",nu_dcosy_truth,"nu_dcosy_truth[mcevts_truth]/F");
    CreateBranch("nu_dcosz_truth",nu_dcosz_truth,"nu_dcosz_truth[mcevts_truth]/F");
    CreateBranch("lep_mom_truth",lep_mom_truth,"lep_mom_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosx_truth",lep_dcosx_truth,"lep_dcosx_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosy_truth",lep_dcosy_truth,"lep_dcosy_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosz_truth",lep_dcosz_truth,"lep_dcosz_truth[mcevts_truth]/F");

    CreateBranch("vx_flux",vx_flux,"vx_flux[mcevts_truth]/F");
    CreateBranch("vy_flux",vy_flux,"vy_flux[mcevts_truth]/F");
    CreateBranch("vz_flux",vz_flux,"vz_flux[mcevts_truth]/F");
    CreateBranch("pdpx_flux",pdpx_flux,"pdpx_flux[mcevts_truth]/F");
    CreateBranch("pdpy_flux",pdpy_flux,"pdpy_flux[mcevts_truth]/F");
    CreateBranch("pdpz_flux",pdpz_flux,"pdpz_flux[mcevts_truth]/F");
    CreateBranch("ppdxdz_flux",ppdxdz_flux,"ppdxdz_flux[mcevts_truth]/F");
    CreateBranch("ppdydz_flux",ppdydz_flux,"ppdydz_flux[mcevts_truth]/F");
    CreateBranch("pppz_flux",pppz_flux,"pppz_flux[mcevts_truth]/F");
    CreateBranch("ptype_flux",ptype_flux,"ptype_flux[mcevts_truth]/I");
    CreateBranch("ppvx_flux",ppvx_flux,"ppvx_flux[mcevts_truth]/F");
    CreateBranch("ppvy_flux",ppvy_flux,"ppvy_flux[mcevts_truth]/F");
    CreateBranch("ppvz_flux",ppvz_flux,"ppvz_flux[mcevts_truth]/F");
    CreateBranch("muparpx_flux",muparpx_flux,"muparpx_flux[mcevts_truth]/F");
    CreateBranch("muparpy_flux",muparpy_flux,"muparpy_flux[mcevts_truth]/F");
    CreateBranch("muparpz_flux",muparpz_flux,"muparpz_flux[mcevts_truth]/F");   
    CreateBranch("mupare_flux",mupare_flux,"mupare_flux[mcevts_truth]/F");
    CreateBranch("tgen_flux",tgen_flux,"tgen_flux[mcevts_truth]/I");
    CreateBranch("tgptype_flux",tgptype_flux,"tgptype_flux[mcevts_truth]/I");
    CreateBranch("tgppx_flux",tgppx_flux,"tgppx_flux[mcevts_truth]/F");
    CreateBranch("tgppy_flux",tgppy_flux,"tgppy_flux[mcevts_truth]/F");
    CreateBranch("tgppz_flux",tgppz_flux,"tgppz_flux[mcevts_truth]/F");
    CreateBranch("tprivx_flux",tprivx_flux,"tprivx_flux[mcevts_truth]/F");
    CreateBranch("tprivy_flux",tprivy_flux,"tprivy_flux[mcevts_truth]/F");
    CreateBranch("tprivz_flux",tprivz_flux,"tprivz_flux[mcevts_truth]/F");
    CreateBranch("dk2gen_flux",dk2gen_flux,"dk2gen_flux[mcevts_truth]/F");
    CreateBranch("gen2vtx_flux",gen2vtx_flux,"gen2vtx_flux[mcevts_truth]/F");
    CreateBranch("tpx_flux",tpx_flux,"tpx_flux[mcevts_truth]/F");
    CreateBranch("tpy_flux",tpy_flux,"tpy_flux[mcevts_truth]/F");
    CreateBranch("tpz_flux",tpz_flux,"tpz_flux[mcevts_truth]/F");
    CreateBranch("tptype_flux",tptype_flux,"tptype_flux[mcevts_truth]/I");

    CreateBranch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
    CreateBranch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
    CreateBranch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/F");
    CreateBranch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/F");
    CreateBranch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/F");
    CreateBranch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/F");
    CreateBranch("genie_P",genie_P,"genie_P[genie_no_primaries]/F");
    CreateBranch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
    CreateBranch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/F");
    CreateBranch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
    CreateBranch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
    CreateBranch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
  }

  if (hasCryInfo()){
    CreateBranch("mcevts_truthcry",&mcevts_truthcry,"mcevts_truthcry/I");  
    CreateBranch("cry_no_primaries",&cry_no_primaries,"cry_no_primaries/I");
    CreateBranch("cry_primaries_pdg",cry_primaries_pdg,"cry_primaries_pdg[cry_no_primaries]/I");
    CreateBranch("cry_Eng",cry_Eng,"cry_Eng[cry_no_primaries]/F");
    CreateBranch("cry_Px",cry_Px,"cry_Px[cry_no_primaries]/F");
    CreateBranch("cry_Py",cry_Py,"cry_Py[cry_no_primaries]/F");
    CreateBranch("cry_Pz",cry_Pz,"cry_Pz[cry_no_primaries]/F");
    CreateBranch("cry_P",cry_P,"cry_P[cry_no_primaries]/F");
    CreateBranch("cry_StartPointx",cry_StartPointx,"cry_StartPointx[cry_no_primaries]/F");
    CreateBranch("cry_StartPointy",cry_StartPointy,"cry_StartPointy[cry_no_primaries]/F");
    CreateBranch("cry_StartPointz",cry_StartPointz,"cry_StartPointz[cry_no_primaries]/F");
    CreateBranch("cry_StartPointt",cry_StartPointt,"cry_StartPointt[cry_no_primaries]/F");   
    CreateBranch("cry_status_code",cry_status_code,"cry_status_code[cry_no_primaries]/I");
    CreateBranch("cry_mass",cry_mass,"cry_mass[cry_no_primaries]/F");
    CreateBranch("cry_trackID",cry_trackID,"cry_trackID[cry_no_primaries]/I");
    CreateBranch("cry_ND",cry_ND,"cry_ND[cry_no_primaries]/I");
    CreateBranch("cry_mother",cry_mother,"cry_mother[cry_no_primaries]/I");
  }  

  if (hasGeantInfo()){  
    CreateBranch("no_primaries",&no_primaries,"no_primaries/I");
    CreateBranch("geant_list_size",&geant_list_size,"geant_list_size/I");
    CreateBranch("geant_list_size_in_tpcAV",&geant_list_size_in_tpcAV,"geant_list_size_in_tpcAV/I");  
    CreateBranch("pdg",pdg,"pdg[geant_list_size]/I");
    CreateBranch("status",status,"status[geant_list_size]/I");
    CreateBranch("Mass",Mass,"Mass[geant_list_size]/F");
    CreateBranch("Eng",Eng,"Eng[geant_list_size]/F");
    CreateBranch("EndE",EndE,"EndE[geant_list_size]/F");
    CreateBranch("Px",Px,"Px[geant_list_size]/F");
    CreateBranch("Py",Py,"Py[geant_list_size]/F");
    CreateBranch("Pz",Pz,"Pz[geant_list_size]/F");
    CreateBranch("P",P,"P[geant_list_size]/F");
    CreateBranch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
    CreateBranch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
    CreateBranch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_StartPointx",sp_charge_corrected_StartPointx,"sp_charge_corrected_StartPointx[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_StartPointy",sp_charge_corrected_StartPointy,"sp_charge_corrected_StartPointy[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_StartPointz",sp_charge_corrected_StartPointz,"sp_charge_corrected_StartPointz[geant_list_size]/F"); 
    CreateBranch("StartT",StartT,"StartT[geant_list_size]/F");
    CreateBranch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
    CreateBranch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
    CreateBranch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_EndPointx",sp_charge_corrected_EndPointx,"sp_charge_corrected_EndPointx[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_EndPointy",sp_charge_corrected_EndPointy,"sp_charge_corrected_EndPointy[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_EndPointz",sp_charge_corrected_EndPointz,"sp_charge_corrected_EndPointz[geant_list_size]/F"); 
    CreateBranch("EndT",EndT,"EndT[geant_list_size]/F");
    CreateBranch("theta",theta,"theta[geant_list_size]/F");
    CreateBranch("phi",phi,"phi[geant_list_size]/F");
    CreateBranch("theta_xz",theta_xz,"theta_xz[geant_list_size]/F");
    CreateBranch("theta_yz",theta_yz,"theta_yz[geant_list_size]/F");
    CreateBranch("pathlen",pathlen,"pathlen[geant_list_size]/F");
    CreateBranch("inTPCActive",inTPCActive,"inTPCActive[geant_list_size]/I");  
    CreateBranch("StartPointx_tpcAV",StartPointx_tpcAV,"StartPointx_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointy_tpcAV",StartPointy_tpcAV,"StartPointy_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointz_tpcAV",StartPointz_tpcAV,"StartPointz_tpcAV[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_StartPointx_tpcAV",sp_charge_corrected_StartPointx_tpcAV,"sp_charge_corrected_StartPointx_tpcAV[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_StartPointy_tpcAV",sp_charge_corrected_StartPointy_tpcAV,"sp_charge_corrected_StartPointy_tpcAV[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_StartPointz_tpcAV",sp_charge_corrected_StartPointz_tpcAV,"sp_charge_corrected_StartPointz_tpcAV[geant_list_size]/F"); 
    CreateBranch("StartT_tpcAV",StartT_tpcAV,"StartT_tpcAV[geant_list_size]/F");
    CreateBranch("StartE_tpcAV",StartE_tpcAV,"StartE_tpcAV[geant_list_size]/F");
    CreateBranch("StartP_tpcAV",StartP_tpcAV,"StartP_tpcAV[geant_list_size]/F");
    CreateBranch("StartPx_tpcAV",StartPx_tpcAV,"StartPx_tpcAV[geant_list_size]/F");
    CreateBranch("StartPy_tpcAV",StartPy_tpcAV,"StartPy_tpcAV[geant_list_size]/F");
    CreateBranch("StartPz_tpcAV",StartPz_tpcAV,"StartPz_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointx_tpcAV",EndPointx_tpcAV,"EndPointx_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointy_tpcAV",EndPointy_tpcAV,"EndPointy_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointz_tpcAV",EndPointz_tpcAV,"EndPointz_tpcAV[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_EndPointx_tpcAV",sp_charge_corrected_EndPointx_tpcAV,"sp_charge_corrected_EndPointx_tpcAV[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_EndPointy_tpcAV",sp_charge_corrected_EndPointy_tpcAV,"sp_charge_corrected_EndPointy_tpcAV[geant_list_size]/F"); 
    CreateBranch("sp_charge_corrected_EndPointz_tpcAV",sp_charge_corrected_EndPointz_tpcAV,"sp_charge_corrected_EndPointz_tpcAV[geant_list_size]/F"); 
    CreateBranch("EndT_tpcAV",EndT_tpcAV,"EndT_tpcAV[geant_list_size]/F");
    CreateBranch("EndE_tpcAV",EndE_tpcAV,"EndE_tpcAV[geant_list_size]/F");
    CreateBranch("EndP_tpcAV",EndP_tpcAV,"EndP_tpcAV[geant_list_size]/F");
    CreateBranch("EndPx_tpcAV",EndPx_tpcAV,"EndPx_tpcAV[geant_list_size]/F");
    CreateBranch("EndPy_tpcAV",EndPy_tpcAV,"EndPy_tpcAV[geant_list_size]/F");
    CreateBranch("EndPz_tpcAV",EndPz_tpcAV,"EndPz_tpcAV[geant_list_size]/F");
    CreateBranch("pathlen_drifted",pathlen_drifted,"pathlen_drifted[geant_list_size]/F");
    CreateBranch("inTPCDrifted",inTPCDrifted,"inTPCDrifted[geant_list_size]/I");  
    CreateBranch("StartPointx_drifted",StartPointx_drifted,"StartPointx_drifted[geant_list_size]/F");
    CreateBranch("StartPointy_drifted",StartPointy_drifted,"StartPointy_drifted[geant_list_size]/F");
    CreateBranch("StartPointz_drifted",StartPointz_drifted,"StartPointz_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_StartPointx_drifted",sp_charge_corrected_StartPointx_drifted,"sp_charge_corrected_StartPointx_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_StartPointy_drifted",sp_charge_corrected_StartPointy_drifted,"sp_charge_corrected_StartPointy_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_StartPointz_drifted",sp_charge_corrected_StartPointz_drifted,"sp_charge_corrected_StartPointz_drifted[geant_list_size]/F");
    CreateBranch("StartT_drifted",StartT_drifted,"StartT_drifted[geant_list_size]/F");
    CreateBranch("StartE_drifted",StartE_drifted,"StartE_drifted[geant_list_size]/F");
    CreateBranch("StartP_drifted",StartP_drifted,"StartP_drifted[geant_list_size]/F");
    CreateBranch("StartPx_drifted",StartPx_drifted,"StartPx_drifted[geant_list_size]/F");
    CreateBranch("StartPy_drifted",StartPy_drifted,"StartPy_drifted[geant_list_size]/F");
    CreateBranch("StartPz_drifted",StartPz_drifted,"StartPz_drifted[geant_list_size]/F");
    CreateBranch("EndPointx_drifted",EndPointx_drifted,"EndPointx_drifted[geant_list_size]/F");
    CreateBranch("EndPointy_drifted",EndPointy_drifted,"EndPointy_drifted[geant_list_size]/F");
    CreateBranch("EndPointz_drifted",EndPointz_drifted,"EndPointz_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_EndPointx_drifted",sp_charge_corrected_EndPointx_drifted,"sp_charge_corrected_EndPointx_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_EndPointy_drifted",sp_charge_corrected_EndPointy_drifted,"sp_charge_corrected_EndPointy_drifted[geant_list_size]/F");
    CreateBranch("sp_charge_corrected_EndPointz_drifted",sp_charge_corrected_EndPointz_drifted,"sp_charge_corrected_EndPointz_drifted[geant_list_size]/F");
    CreateBranch("EndT_drifted",EndT_drifted,"EndT_drifted[geant_list_size]/F");
    CreateBranch("EndE_drifted",EndE_drifted,"EndE_drifted[geant_list_size]/F");
    CreateBranch("EndP_drifted",EndP_drifted,"EndP_drifted[geant_list_size]/F");
    CreateBranch("EndPx_drifted",EndPx_drifted,"EndPx_drifted[geant_list_size]/F");
    CreateBranch("EndPy_drifted",EndPy_drifted,"EndPy_drifted[geant_list_size]/F");
    CreateBranch("EndPz_drifted",EndPz_drifted,"EndPz_drifted[geant_list_size]/F");
    CreateBranch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
    CreateBranch("Mother",Mother,"Mother[geant_list_size]/I");
    CreateBranch("TrackId",TrackId,"TrackId[geant_list_size]/I");
    CreateBranch("MergedId", MergedId, "MergedId[geant_list_size]/I");
    CreateBranch("origin", origin, "origin[geant_list_size]/I");
    CreateBranch("MCTruthIndex", MCTruthIndex, "MCTruthIndex[geant_list_size]/I");
    CreateBranch("process_primary",process_primary,"process_primary[geant_list_size]/I");
    CreateBranch("processname", processname);
  }

  if (hasMCShowerInfo()){
    CreateBranch("no_mcshowers",&no_mcshowers,"no_mcshowers/I");  
    CreateBranch("mcshwr_origin",mcshwr_origin,"mcshwr_origin[no_mcshowers]/I");
    CreateBranch("mcshwr_pdg",mcshwr_pdg,"mcshwr_pdg[no_mcshowers]/I");
    CreateBranch("mcshwr_TrackId",mcshwr_TrackId,"mcshwr_TrackId[no_mcshowers]/I");
    CreateBranch("mcshwr_Process",mcshwr_Process);
    CreateBranch("mcshwr_startX",mcshwr_startX,"mcshwr_startX[no_mcshowers]/F");
    CreateBranch("mcshwr_startY",mcshwr_startY,"mcshwr_startY[no_mcshowers]/F");
    CreateBranch("mcshwr_startZ",mcshwr_startZ,"mcshwr_startZ[no_mcshowers]/F");
    CreateBranch("mcshwr_endX",mcshwr_endX,"mcshwr_endX[no_mcshowers]/F");
    CreateBranch("mcshwr_endY",mcshwr_endY,"mcshwr_endY[no_mcshowers]/F");
    CreateBranch("mcshwr_endZ",mcshwr_endZ,"mcshwr_endZ[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngX",mcshwr_CombEngX,"mcshwr_CombEngX[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngY",mcshwr_CombEngY,"mcshwr_CombEngY[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngZ",mcshwr_CombEngZ,"mcshwr_CombEngZ[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngPx",mcshwr_CombEngPx,"mcshwr_CombEngPx[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngPy",mcshwr_CombEngPy,"mcshwr_CombEngPy[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngPz",mcshwr_CombEngPz,"mcshwr_CombEngPz[no_mcshowers]/F");
    CreateBranch("mcshwr_CombEngE",mcshwr_CombEngE,"mcshwr_CombEngE[no_mcshowers]/F");
    CreateBranch("mcshwr_dEdx",mcshwr_dEdx,"mcshwr_dEdx[no_mcshowers]/F");
    CreateBranch("mcshwr_StartDirX",mcshwr_StartDirX,"mcshwr_StartDirX[no_mcshowers]/F");
    CreateBranch("mcshwr_StartDirY",mcshwr_StartDirY,"mcshwr_StartDirY[no_mcshowers]/F");
    CreateBranch("mcshwr_StartDirZ",mcshwr_StartDirZ,"mcshwr_StartDirZ[no_mcshowers]/F");
    CreateBranch("mcshwr_isEngDeposited",mcshwr_isEngDeposited,"mcshwr_isEngDeposited[no_mcshowers]/I");
    CreateBranch("mcshwr_Motherpdg",mcshwr_Motherpdg,"mcshwr_Motherpdg[no_mcshowers]/I");
    CreateBranch("mcshwr_MotherTrkId",mcshwr_MotherTrkId,"mcshwr_MotherTrkId[no_mcshowers]/I");
    CreateBranch("mcshwr_MotherProcess",mcshwr_MotherProcess);
    CreateBranch("mcshwr_MotherstartX",mcshwr_MotherstartX,"mcshwr_MotherstartX[no_mcshowers]/F");
    CreateBranch("mcshwr_MotherstartY",mcshwr_MotherstartY,"mcshwr_MotherstartY[no_mcshowers]/F");
    CreateBranch("mcshwr_MotherstartZ",mcshwr_MotherstartZ,"mcshwr_MotherstartZ[no_mcshowers]/F");
    CreateBranch("mcshwr_MotherendX",mcshwr_MotherendX,"mcshwr_MotherendX[no_mcshowers]/F");
    CreateBranch("mcshwr_MotherendY",mcshwr_MotherendY,"mcshwr_MotherendY[no_mcshowers]/F");
    CreateBranch("mcshwr_MotherendZ",mcshwr_MotherendZ,"mcshwr_MotherendZ[no_mcshowers]/F");    
    CreateBranch("mcshwr_Ancestorpdg",mcshwr_Ancestorpdg,"mcshwr_Ancestorpdg[no_mcshowers]/I");
    CreateBranch("mcshwr_AncesotorTrkId",mcshwr_AncestorTrkId,"mcshwr_AncestorTrkId[no_mcshowers]/I");
    CreateBranch("mcshwr_AncesotorProcess",mcshwr_AncestorProcess);   
    CreateBranch("mcshwr_AncestorstartX",mcshwr_AncestorstartX,"mcshwr_AncestorstartX[no_mcshowers]/F");
    CreateBranch("mcshwr_AncestorstartY",mcshwr_AncestorstartY,"mcshwr_AncestorstartY[no_mcshowers]/F");
    CreateBranch("mcshwr_AncestorstartZ",mcshwr_AncestorstartZ,"mcshwr_AncestorstartZ[no_mcshowers]/F");
    CreateBranch("mcshwr_AncestorendX",mcshwr_AncestorendX,"mcshwr_AncestorendX[no_mcshowers]/F");
    CreateBranch("mcshwr_AncestorendY",mcshwr_AncestorendY,"mcshwr_AncestorendY[no_mcshowers]/F");
    CreateBranch("mcshwr_AncestorendZ",mcshwr_AncestorendZ,"mcshwr_AncestorendZ[no_mcshowers]/F");
  }   
  
  if (hasMCTrackInfo()){
    CreateBranch("no_mctracks",&no_mctracks,"no_mctracks/I");  
    CreateBranch("mctrk_origin",mctrk_origin,"mctrk_origin[no_mctracks]/I");
    CreateBranch("mctrk_pdg",mctrk_pdg,"mctrk_pdg[no_mctracks]/I");
    CreateBranch("mctrk_TrackId",mctrk_TrackId,"mctrk_TrackId[no_mctracks]/I");
    CreateBranch("mctrk_Process",mctrk_Process);
    CreateBranch("mctrk_startX",mctrk_startX,"mctrk_startX[no_mctracks]/F");
    CreateBranch("mctrk_startY",mctrk_startY,"mctrk_startY[no_mctracks]/F");
    CreateBranch("mctrk_startZ",mctrk_startZ,"mctrk_startZ[no_mctracks]/F");
    CreateBranch("mctrk_endX",mctrk_endX,"mctrk_endX[no_mctracks]/F");
    CreateBranch("mctrk_endY",mctrk_endY,"mctrk_endY[no_mctracks]/F");
    CreateBranch("mctrk_endZ",mctrk_endZ,"mctrk_endZ[no_mctracks]/F");
    CreateBranch("mctrk_startX_drifted",mctrk_startX_drifted,"mctrk_startX_drifted[no_mctracks]/F");
    CreateBranch("mctrk_startY_drifted",mctrk_startY_drifted,"mctrk_startY_drifted[no_mctracks]/F");
    CreateBranch("mctrk_startZ_drifted",mctrk_startZ_drifted,"mctrk_startZ_drifted[no_mctracks]/F");
    CreateBranch("mctrk_endX_drifted",mctrk_endX_drifted,"mctrk_endX_drifted[no_mctracks]/F");      
    CreateBranch("mctrk_endY_drifted",mctrk_endY_drifted,"mctrk_endY_drifted[no_mctracks]/F");      
    CreateBranch("mctrk_endZ_drifted",mctrk_endZ_drifted,"mctrk_endZ_drifted[no_mctracks]/F");      
    CreateBranch("mctrk_len_drifted",mctrk_len_drifted,"mctrk_len_drifted[no_mctracks]/F");
    CreateBranch("mctrk_p_drifted",mctrk_p_drifted,"mctrk_p_drifted[no_mctracks]/F");
    CreateBranch("mctrk_px_drifted",mctrk_px_drifted,"mctrk_px_drifted[no_mctracks]/F");
    CreateBranch("mctrk_py_drifted",mctrk_py_drifted,"mctrk_py_drifted[no_mctracks]/F");
    CreateBranch("mctrk_pz_drifted",mctrk_pz_drifted,"mctrk_pz_drifted[no_mctracks]/F");
    CreateBranch("mctrk_Motherpdg",mctrk_Motherpdg,"mctrk_Motherpdg[no_mctracks]/I");
    CreateBranch("mctrk_MotherTrkId",mctrk_MotherTrkId,"mctrk_MotherTrkId[no_mctracks]/I");
    CreateBranch("mctrk_MotherProcess",mctrk_MotherProcess);
    CreateBranch("mctrk_MotherstartX",mctrk_MotherstartX,"mctrk_MotherstartX[no_mctracks]/F");
    CreateBranch("mctrk_MotherstartY",mctrk_MotherstartY,"mctrk_MotherstartY[no_mctracks]/F");
    CreateBranch("mctrk_MotherstartZ",mctrk_MotherstartZ,"mctrk_MotherstartZ[no_mctracks]/F");
    CreateBranch("mctrk_MotherendX",mctrk_MotherendX,"mctrk_MotherendX[no_mctracks]/F");
    CreateBranch("mctrk_MotherendY",mctrk_MotherendY,"mctrk_MotherendY[no_mctracks]/F");
    CreateBranch("mctrk_MotherendZ",mctrk_MotherendZ,"mctrk_MotherendZ[no_mctracks]/F");    
    CreateBranch("mctrk_Ancestorpdg",mctrk_Ancestorpdg,"mctrk_Ancestorpdg[no_mctracks]/I");
    CreateBranch("mctrk_AncesotorTrkId",mctrk_AncestorTrkId,"mctrk_AncestorTrkId[no_mctracks]/I");
    CreateBranch("mctrk_AncesotorProcess",mctrk_AncestorProcess);   
    CreateBranch("mctrk_AncestorstartX",mctrk_AncestorstartX,"mctrk_AncestorstartX[no_mctracks]/F");
    CreateBranch("mctrk_AncestorstartY",mctrk_AncestorstartY,"mctrk_AncestorstartY[no_mctracks]/F");
    CreateBranch("mctrk_AncestorstartZ",mctrk_AncestorstartZ,"mctrk_AncestorstartZ[no_mctracks]/F");
    CreateBranch("mctrk_AncestorendX",mctrk_AncestorendX,"mctrk_AncestorendX[no_mctracks]/F");
    CreateBranch("mctrk_AncestorendY",mctrk_AncestorendY,"mctrk_AncestorendY[no_mctracks]/F");
    CreateBranch("mctrk_AncestorendZ",mctrk_AncestorendZ,"mctrk_AncestorendZ[no_mctracks]/F");
  }  

  if (hasAuxDetector()) {
    // Geant information is required to fill aux detector information.
    // if fSaveGeantInfo is not set to true, show an error message and quit!
    if (!hasGeantInfo()){
      throw art::Exception(art::errors::Configuration)
	<< "Saving Auxiliary detector information requies saving GEANT information, "
	<<"please set fSaveGeantInfo flag to true in your fhicl file and rerun.\n"; 
    }    
    std::ostringstream sstr;
    sstr << "[" << kMaxAuxDets << "]";
    std::string MaxAuxDetIndexStr = sstr.str();
    CreateBranch("NAuxDets",     NAuxDets, "NAuxDets[geant_list_size]/s");
    CreateBranch("AuxDetID",     AuxDetID, "AuxDetID[geant_list_size]" + MaxAuxDetIndexStr + "/S");
    CreateBranch("AuxDetEntryX", entryX,   "AuxDetEntryX[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryY", entryY,   "AuxDetEntryY[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryZ", entryZ,   "AuxDetEntryZ[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryT", entryT,   "AuxDetEntryT[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitX",  exitX,    "AuxDetExitX[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitY",  exitY,    "AuxDetExitY[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitZ",  exitZ,    "AuxDetExitZ[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitT",  exitT,    "AuxDetExitT[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPx", exitPx,   "AuxDetExitPx[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPy", exitPy,   "AuxDetExitPy[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPz", exitPz,   "AuxDetExitPz[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("CombinedEnergyDep", CombinedEnergyDep,
		 "CombinedEnergyDep[geant_list_size]" + MaxAuxDetIndexStr + "/F");
  } // if hasAuxDetector
  
} // microboone::AnalysisTreeDataStruct::SetAddresses()


//------------------------------------------------------------------------------
//---  AnalysisTree
//---

microboone::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr), fPOT(nullptr),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fCryGenModuleLabel        (pset.get< std::string >("CryGenModuleLabel")       ), 
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")           ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")     ),
  fPandoraNuVertexModuleLabel (pset.get< std::vector<std::string> >("PandoraNuVertexModuleLabel")     ),
  fOpFlashModuleLabel       (pset.get< std::vector<std::string> >("OpFlashModuleLabel")      ),
  fMCShowerModuleLabel      (pset.get< std::string >("MCShowerModuleLabel")     ),  
  fMCTrackModuleLabel      (pset.get< std::string >("MCTrackModuleLabel")     ),  
  fTrackModuleLabel         (pset.get< std::vector<std::string> >("TrackModuleLabel")),
  fVertexModuleLabel        (pset.get< std::vector<std::string> >("VertexModuleLabel")),
  fShowerModuleLabel        (pset.get< std::vector<std::string> >("ShowerModuleLabel")),
  fCalorimetryModuleLabel   (pset.get< std::vector<std::string> >("CalorimetryModuleLabel")),
  fParticleIDModuleLabel    (pset.get< std::vector<std::string> >("ParticleIDModuleLabel")   ),
  fT0RecoAnodeCathodePiercingLabel (pset.get< std::vector<std::string> >("T0RecoAnodeCathodePiercingLabel") ),
  fFlashT0FinderLabel       (pset.get< std::vector<std::string> >("FlashT0FinderLabel")   ),
  fMCT0FinderLabel          (pset.get< std::vector<std::string> >("MCT0FinderLabel")   ),
  fTrackMCSFitLabel         (pset.get< std::vector<std::string> >("TrackMCSFitLabel")),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")),   
  fCosmicClusterTaggerAssocLabel (pset.get< std::string >("CosmicClusterTaggerAssocLabel")), 
  fSWTriggerLabel           (pset.get< std::string >("SWTriggerModuleLabel")),
  fOpticalFilterLabel       (pset.get< std::string >("OpticalFilterLabel")),
  fUseBuffer                (pset.get< bool >("UseBuffers", false)),
  fSaveAuxDetInfo           (pset.get< bool >("SaveAuxDetInfo", false)),
  fSaveCryInfo              (pset.get< bool >("SaveCryInfo", false)),  
  fSaveGenieInfo	    (pset.get< bool >("SaveGenieInfo", false)), 
  fSaveGeantInfo	    (pset.get< bool >("SaveGeantInfo", false)), 
  fSaveMCShowerInfo	    (pset.get< bool >("SaveMCShowerInfo", false)), 
  fSaveMCTrackInfo	    (pset.get< bool >("SaveMCTrackInfo", false)), 
  fSaveHitInfo	            (pset.get< bool >("SaveHitInfo", false)), 
  fSaveRawDigitInfo	            (pset.get< bool >("SaveRawDigitInfo", false)), 
  fSaveCalWireInfo	    (pset.get< bool >("SaveRawDigitInfo", false)), 
  fSaveSimChannelInfo	    (pset.get< bool >("SaveSimChannelInfo", false)), 
  fSaveRawWaveFormInfo	    (pset.get< bool >("SaveRawWaveFormInfo", false)), 
  fSaveCalibWaveFormInfo    (pset.get< bool >("SaveCalibWaveFormInfo", false)),
  fSaveTrackInfo	    (pset.get< bool >("SaveTrackInfo", false)), 
  fSaveVertexInfo	    (pset.get< bool >("SaveVertexInfo", false)),
  fSaveClusterInfo	    (pset.get< bool >("SaveClusterInfo", false)),
  fSavePandoraNuVertexInfo  (pset.get< bool >("SavePandoraNuVertexInfo", false)),
  fSaveFlashInfo            (pset.get< bool >("SaveFlashInfo", false)),
  fSaveShowerInfo            (pset.get< bool >("SaveShowerInfo", false)),
  fSavePFParticleInfo	    (pset.get< bool >("SavePFParticleInfo", false)),
  fSaveSWTriggerInfo        (pset.get< bool >("SaveSWTriggerInfo", false)),
  fSaveOpticalFilterInfo    (pset.get< bool >("SaveOpticalFilterInfo", false)),
  fIsOverlay                (pset.get< bool >("IsOverlay", false)),
  fCosmicTaggerAssocLabel  (pset.get<std::vector< std::string > >("CosmicTaggerAssocLabel") ),
  fContainmentTaggerAssocLabel  (pset.get<std::vector< std::string > >("ContainmentTaggerAssocLabel") ),
  fFlashMatchAssocLabel (pset.get<std::vector< std::string > >("FlashMatchAssocLabel") ),
  bIgnoreMissingShowers     (pset.get< bool >("IgnoreMissingShowers", false)),
  isCosmics(false),
  fSaveCaloCosmics          (pset.get< bool >("SaveCaloCosmics",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01)),
  fRawDigitROI              (pset.get< float>("RawDigitROI",3.0)), 
  fCalWireROI               (pset.get< float>("CalWireROI",3.0)), 
  fRawDigitThresh           (pset.get< float>("RawDigitThresh",10)), 
  fCalWireThresh            (pset.get< float>("CalWireThresh",10))  
{

  if (fSavePFParticleInfo) fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");

  if (fSaveAuxDetInfo == true) fSaveGeantInfo = true;
  if (fSaveRawDigitInfo == true) fSaveHitInfo = true;
  if (fSaveCalWireInfo == true) fSaveHitInfo = true;
  if (fSaveSimChannelInfo == true) fSaveHitInfo = true;
  
  mf::LogInfo("AnalysisTree") << "Configuration:"
			      << "\n  UseBuffers: " << std::boolalpha << fUseBuffer
    ;
  if (GetNTrackers() > kMaxTrackers) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxTrackers
      << " tracking algorithms, but " << GetNTrackers() << " are specified."
      << "\nYou can increase kMaxTrackers and recompile.";
  } // if too many trackers
  if (fTrackModuleLabel.size() != fCalorimetryModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCalorimetryModuleLabel.size() = "<<fCalorimetryModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fParticleIDModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fParticleIDModuleLabel.size() = "<<fParticleIDModuleLabel.size();
  }
  
  if (fTrackModuleLabel.size() != fT0RecoAnodeCathodePiercingLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fT0RecoAnodeCathodePiercingLabel.size() = "<<fT0RecoAnodeCathodePiercingLabel.size();
  }
  
  if (fTrackModuleLabel.size() != fFlashT0FinderLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fFlashT0FinderLabel.size() = "<<fFlashT0FinderLabel.size();
  }
  if (fTrackModuleLabel.size() != fMCT0FinderLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fMCT0FinderLabel.size() = "<<fMCT0FinderLabel.size();
  }
  if (fTrackModuleLabel.size() != fCosmicTaggerAssocLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCosmicTaggerAssocLabel.size() = "<<fCosmicTaggerAssocLabel.size();
  }
  if (fTrackModuleLabel.size() != fFlashMatchAssocLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fFlashMatchAssocLabel.size() = "<<fFlashMatchAssocLabel.size();
  }
  if (fTrackModuleLabel.size() != fContainmentTaggerAssocLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fContainmentTaggerAssocLabel.size() = "<<fContainmentTaggerAssocLabel.size();
  }

  if (fTrackModuleLabel.size() != fTrackMCSFitLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fTrackMCSFitLabel.size() = "<<fTrackMCSFitLabel.size();
  }


  if (GetNVertexAlgos() > kMaxVertexAlgos) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxVertexAlgos
      << " vertex algorithms, but " << GetNVertexAlgos() << " are specified."
      << "\nYou can increase kMaxVertexAlgos and recompile.";
  } // if too many vertices

  if (GetNNuVertexAlgos() > kMaxVertexAlgos) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxVertexAlgos
      << " neutrino vertex algorithms, but " << GetNNuVertexAlgos() << " are specified."
      << "\nYou can increase kMaxVertexAlgos and recompile.";
  } // if too many vertices

  if (GetNFlashAlgos() > kMaxFlashAlgos) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxFlashAlgos
      << " flash algorithms, but " << GetNFlashAlgos() << " are specified."
      << "\nYou can increase kMaxFlashAlgos and recompile.";
  } // if too many flashes

  // Get the tool for MC Truth matching
  fMCTruthMatching = art::make_tool<truth::IMCTruthMatching>(pset.get<fhicl::ParameterSet>("MCTruthMatching"));
} // microboone::AnalysisTree::AnalysisTree()

//-------------------------------------------------
microboone::AnalysisTree::~AnalysisTree()
{
  DestroyData();
}

void microboone::AnalysisTree::CreateTree(bool bClearData /* = false */) {
  if (!fTree) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
  }
  if (!fPOT) {
    art::ServiceHandle<art::TFileService> tfs;
    fPOT = tfs->make<TTree>("pottree","pot tree");
    fPOT->Branch("pot",&SubRunData.pot,"pot/D");
    fPOT->Branch("potbnbETOR860",&SubRunData.potbnbETOR860,"potbnbETOR860/D");
    fPOT->Branch("potbnbETOR875",&SubRunData.potbnbETOR875,"potbnbETOR875/D");
    fPOT->Branch("potnumiETORTGT",&SubRunData.potnumiETORTGT,"potnumiETORTGT/D");
    fPOT->Branch("run",&SubRunData.run,"run/I");
    fPOT->Branch("subrun",&SubRunData.subrun,"subrun/I");
    fPOT->Branch("begintime",&SubRunData.begintime,"begintime/D");
    fPOT->Branch("endtime",&SubRunData.endtime,"endtime/D");
  }
  CreateData(bClearData);
  SetAddresses();
} // microboone::AnalysisTree::CreateTree()


void microboone::AnalysisTree::beginSubRun(const art::SubRun& sr)
{

//  art::Handle< sumdata::POTSummary > potListHandle;
//  //sr.getByLabel(fPOTModuleLabel,potListHandle);
//
//  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
//    SubRunData.pot=potListHandle->totpot;
//  else
//    SubRunData.pot=0.;

}

void microboone::AnalysisTree::endSubRun(const art::SubRun& sr)
{

  SubRunData.run = sr.run();
  SubRunData.subrun = sr.subRun();
  SubRunData.begintime = sr.beginTime().value();
  SubRunData.endtime = sr.endTime().value();

  art::Handle< sumdata::POTSummary > potListHandle;
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    SubRunData.pot=potListHandle->totpot;
  else
    SubRunData.pot=0.;
  
  art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR860;
  if (sr.getByLabel("beamdata","bnbETOR860",potSummaryHandlebnbETOR860)){
    SubRunData.potbnbETOR860 = potSummaryHandlebnbETOR860->totpot;
  }
  else
    SubRunData.potbnbETOR860 = 0;

  art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
  if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
    SubRunData.potbnbETOR875 = potSummaryHandlebnbETOR875->totpot;
  }
  else
    SubRunData.potbnbETOR875 = 0;

  art::Handle<sumdata::POTSummary> potSummaryHandlenumiETORTGT;
  if (sr.getByLabel("beamdata","numiETORTGT",potSummaryHandlenumiETORTGT)){
    SubRunData.potnumiETORTGT = potSummaryHandlenumiETORTGT->totpot;
  }
  else
    SubRunData.potnumiETORTGT = 0;

  if (fPOT) fPOT->Fill();

}

void microboone::AnalysisTree::analyze(const art::Event& evt)
{
  //services
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  // collect the sizes which might me needed to resize the tree data structure:
  bool isMC = !evt.isRealData();
  if (fIsOverlay) isMC = true;
    
  // If this is MC then we want to "rebuild"
  // For the BackTracker this call will be a noop (the interface intercepts) since it is a service
  // For the assocaitions version then it builds out the maps
  if(isMC)
    fMCTruthMatching->Rebuild(evt);

  // * hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // * clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (fSaveClusterInfo){
    if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);
  }        

  /* flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fOpFlashModuleLabel,flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);
  */

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (isMC){
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
  }    

  // *MC truth cosmic generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthcryListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistcry;
  if (isMC && fSaveCryInfo){
    if (evt.getByLabel(fCryGenModuleLabel,mctruthcryListHandle))
      art::fill_ptr_vector(mclistcry, mctruthcryListHandle);
  }       
  
  // *MC Shower information
  art::Handle< std::vector<sim::MCShower> > mcshowerh;
  if (isMC)
    evt.getByLabel(fMCShowerModuleLabel, mcshowerh);
  
  int nMCShowers = 0;
  if (fSaveMCShowerInfo && mcshowerh.isValid())
    nMCShowers = mcshowerh->size();
	
  // *MC Track information
  art::Handle< std::vector<sim::MCTrack> > mctrackh;
  if (isMC)
    evt.getByLabel(fMCTrackModuleLabel, mctrackh);
  
  int nMCTracks = 0;
  if (fSaveMCTrackInfo && mctrackh.isValid())
    nMCTracks = mctrackh->size();	
	
  art::Ptr<simb::MCTruth> mctruthcry;
  int nCryPrimaries = 0;
   
  if (fSaveCryInfo){
    mctruthcry = mclistcry[0];      
    nCryPrimaries = mctruthcry->NParticles();  
  } 
  
  int nGeniePrimaries = 0, nGEANTparticles = 0;
  
  art::Ptr<simb::MCTruth> mctruth;

  
  if (isMC) { //is MC

    //find origin
    auto allmclists = evt.getMany< std::vector<simb::MCTruth> >();
    for(size_t mcl = 0; mcl < allmclists.size(); ++mcl){
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
      for(size_t m = 0; m < mclistHandle->size(); ++m){
	art::Ptr<simb::MCTruth> mct(mclistHandle, m);
	if (mct->Origin() == simb::kCosmicRay) isCosmics = true;
      }
    }
    if (fSaveCaloCosmics) isCosmics = false; //override to save calo info

    // GENIE
    if (!mclist.empty()){//at least one mc record

      //        double maxenergy = -1;
      //        int imc0 = 0;
      //        for (std::map<art::Ptr<simb::MCTruth>,double>::iterator ii=mctruthemap.begin(); ii!=mctruthemap.end(); ++ii){
      //          if ((ii->second)>maxenergy){
      //            maxenergy = ii->second;
      //            mctruth = ii->first;
      //            imc = imc0;
      //          }
      //          imc0++;
      //        }

      //imc = 0; //set imc to 0 to solve a confusion for BNB+cosmic files where there are two MCTruth
      mctruth = mclist[0];

      if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
      //} //end (fSaveGenieInfo)
      
      const sim::ParticleList& plist = fMCTruthMatching->ParticleList();

      nGEANTparticles = plist.size();

      // to know the number of particles in AV would require
      // looking at all of them; so we waste some memory here
    } // if have MC truth
    MF_LOG_DEBUG("AnalysisTree") << "Expected "
			      << nGEANTparticles << " GEANT particles, "
			      << nGeniePrimaries << " GENIE particles";
  } // if MC
  
  CreateData(); // tracker data is created with default constructor
  if (fSaveGenieInfo)
    fData->ResizeGenie(nGeniePrimaries);
  if (fSaveCryInfo)
    fData->ResizeCry(nCryPrimaries);
  if (fSaveGeantInfo)    
    fData->ResizeGEANT(nGEANTparticles);  
  if (fSaveMCShowerInfo)
    fData->ResizeMCShower(nMCShowers); 
  if (fSaveMCTrackInfo)
    fData->ResizeMCTrack(nMCTracks);     
    
  fData->ClearLocalData(); // don't bother clearing tracker data yet
  
  //  const size_t Nplanes       = 3; // number of wire planes; pretty much constant...
  const size_t NTrackers = GetNTrackers(); // number of trackers passed into fTrackModuleLabel
  const size_t NShowerAlgos = GetNShowerAlgos(); // number of shower algorithms into fShowerModuleLabel
  const size_t NHits     = hitlist.size(); // number of hits
  const size_t NVertexAlgos = GetNVertexAlgos(); // number of vertex algos
  const size_t NNuVertexAlgos = GetNNuVertexAlgos(); // number of neutrino vertex algos
  const size_t NFlashAlgos = GetNFlashAlgos(); // number of flash algos
  const size_t NClusters = clusterlist.size(); //number of clusters
  // make sure there is the data, the tree and everything;
  CreateTree();

  /// transfer the run and subrun data to the tree data object
  //  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  fData->isdata = int(!isMC);

  // * raw trigger
  art::Handle< std::vector<raw::Trigger>> triggerListHandle;
  std::vector<art::Ptr<raw::Trigger>> triggerlist;
  if (evt.getByLabel(fDigitModuleLabel, triggerListHandle))
    art::fill_ptr_vector(triggerlist, triggerListHandle);

  if (triggerlist.size()){
    fData->triggernumber = triggerlist[0]->TriggerNumber();
    fData->triggertime   = triggerlist[0]->TriggerTime();
    fData->beamgatetime  = triggerlist[0]->BeamGateTime();
    fData->triggerbits   = triggerlist[0]->TriggerBits();
  }

  // * vertices
  std::vector< art::Handle< std::vector<recob::Vertex> > > vertexListHandle(NVertexAlgos);
  std::vector< std::vector<art::Ptr<recob::Vertex> > > vertexlist(NVertexAlgos);
  for (unsigned int it = 0; it < NVertexAlgos; ++it){
    if (evt.getByLabel(fVertexModuleLabel[it],vertexListHandle[it]))
      art::fill_ptr_vector(vertexlist[it], vertexListHandle[it]);
  }

  // * nu vertices
  std::vector< std::vector<art::Ptr<recob::Vertex> > >     nuvertexlist(NNuVertexAlgos);
  std::vector< std::vector<art::Ptr<recob::PFParticle> > > nuvertexlistToPfp(NNuVertexAlgos);
  for (unsigned int vtxLabel = 0; vtxLabel < NNuVertexAlgos; vtxLabel++){
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel[vtxLabel], vertexVector, particlesToVertices);
    
    for (lar_pandora::PFParticlesToVertices::iterator it = particlesToVertices.begin(); it != particlesToVertices.end(); ++it){

      art::Ptr<recob::PFParticle> pfParticle = it->first;
      lar_pandora::VertexVector   vertex_v   = it->second;
      if (vertex_v.empty()) continue;
      if (lar_pandora::LArPandoraHelper::IsNeutrino(pfParticle)) {
        if (vertex_v.size() == 1) { // require 1 vtx associated to the neutrino PFP
          nuvertexlist[vtxLabel].emplace_back(vertex_v[0]); 
          nuvertexlistToPfp[vtxLabel].emplace_back(pfParticle);
        }
      }
    }
  }

  // * flashes
  std::vector< art::Handle< std::vector<recob::OpFlash> > > flashListHandle(NFlashAlgos);
  std::vector< std::vector<art::Ptr<recob::OpFlash> > > flashlist(NFlashAlgos);
  for (unsigned int it = 0; it < NFlashAlgos; ++it){
    if (evt.getByLabel(fOpFlashModuleLabel[it],flashListHandle[it]))
      art::fill_ptr_vector(flashlist[it], flashListHandle[it]);
  }

  // * PFParticles
  lar_pandora::PFParticleVector pfparticlelist;
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPFParticleModuleLabel, pfparticlelist, pfParticleToClusterMap);
  
  // * tracks
  std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(NTrackers);
  std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(NTrackers);
  for (unsigned int it = 0; it < NTrackers; ++it){
    if (evt.getByLabel(fTrackModuleLabel[it],trackListHandle[it]))
      art::fill_ptr_vector(tracklist[it], trackListHandle[it]);
  }

  // * tracks MCS 
  std::vector< art::Handle< std::vector<recob::MCSFitResult> > > mcsfitListHandle(NTrackers);
  std::vector< std::vector<art::Ptr<recob::MCSFitResult> > > mcsfitlist(NTrackers);
  for (unsigned int it = 0; it < NTrackers; ++it){
    if (evt.getByLabel(fTrackMCSFitLabel[it],mcsfitListHandle[it]))
      art::fill_ptr_vector(mcsfitlist[it], mcsfitListHandle[it]);
  }  

  // * showers
  // It seems that sometimes the shower data product does not exist;
  // in that case, we store a nullptr in place of the pointer to the vector;
  // the data structure itself is owned by art::Event and we should not try
  // to manage its memory
  std::vector<std::vector<recob::Shower> const*> showerList;
  std::vector< art::Handle< std::vector<recob::Shower> > > showerListHandle(fShowerModuleLabel.size());
  showerList.reserve(fShowerModuleLabel.size());
  for (art::InputTag ShowerInputTag: fShowerModuleLabel) {
    art::Handle<std::vector<recob::Shower>> ShowerHandle;
    if (!evt.getByLabel(ShowerInputTag, ShowerHandle)) {
      showerList.push_back(nullptr);
      if (!bIgnoreMissingShowers) {
        throw art::Exception(art::errors::ProductNotFound)
          << "Showers with input tag '" << ShowerInputTag.encode()
          << "' were not found in the event."
	  " If you really know what you are doing,"
	  " set AnalysisTree's configuration parameter IgnoreMissingShowers"
	  " to \"true\"; the lack of any shower set will be tolerated,"
	  " and the shower list in the corresponding event will be set to"
	  " a list of one shower, with an invalid ID.\n";
      } // if bIgnoreMissingShowers
      else {
        // this message is more alarming than strictly necessary; by design.
        mf::LogError("AnalysisTree")
          << "No showers found for input tag '" << ShowerInputTag.encode()
          << "' ; FILLING WITH FAKE DATA AS FOR USER'S REQUEST";
      }
    }
    else showerList.push_back(ShowerHandle.product());
    
    showerListHandle.push_back(ShowerHandle); // either way, put it into the handle list
  } // for shower input tag
  
  art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > fluxlist;
  if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
    art::fill_ptr_vector(fluxlist, mcfluxListHandle);

  std::vector<const sim::AuxDetSimChannel*> fAuxDetSimChannels;
  if (fSaveAuxDetInfo){
    evt.getView(fLArG4ModuleLabel, fAuxDetSimChannels);
  }

  fData->run = evt.run();
  fData->subrun = evt.subRun();
  fData->event = evt.id().event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fData->evttime = tts.AsDouble();

  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beamdata",beam)){
    fData->beamtime = (double)beam->get_t_ms();
    fData->beamtime/=1000.; //in second
    std::map<std::string, std::vector<double>> datamap = beam->GetDataMap();
    if (datamap["E:TOR860"].size()){
      fData->potbnb = datamap["E:TOR860"][0];
    }
    if (datamap["E:TORTGT"].size()){
      fData->potnumitgt = datamap["E:TORTGT"][0];
    }
    if (datamap["E:TOR101"].size()){
      fData->potnumi101 = datamap["E:TOR101"][0];
    }
  }

  //*****************************
  //
  // EventWeight
  //
  //*****************************

  art::Handle< std::vector< evwgh::MCEventWeight > > evtWeights;

  if (evt.getByLabel("genieeventweight",evtWeights)) {
    const std::vector< evwgh::MCEventWeight > * evtwgt_vec = evtWeights.product();

    evwgh::MCEventWeight evtwgt = evtwgt_vec->at(0); // just for the first neutrino interaction
    std::map<std::string, std::vector<double>> evtwgt_map = evtwgt.fWeight;
    int countFunc = 0;
    // loop over the map and save the name of the function and the vector of weights for each function
    for(std::map<std::string, std::vector<double>>::iterator it = evtwgt_map.begin(); it != evtwgt_map.end(); ++it) {
      fData->evtwgt_funcname.push_back(it->first);      // filling the name of the function
      fData->evtwgt_weight.push_back(it->second);       // filling the vector with the weights
      std::vector<double> mytemp = it->second;          // getting the vector of weights

      fData->evtwgt_nweight.push_back(mytemp.size());   // filling the number of weights
      countFunc++;
    }
    fData->evtwgt_nweight.resize(countFunc);
    fData->evtwgt_funcname.resize(countFunc);
    fData->evtwgt_nfunc = countFunc;                    // saving the number of functions used
  }

  //*****************************
  //
  // Software Trigger
  //
  //***************************** 

  if (fSaveSWTriggerInfo) {
    art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    evt.getByLabel(fSWTriggerLabel, softwareTriggerHandle);

    if (!softwareTriggerHandle.isValid() || softwareTriggerHandle.failedToGet()){
      std::cerr << "Failed to get software trigget data product with label " << fSWTriggerLabel << std::endl;
    }

    int nAlgo = softwareTriggerHandle->getNumberOfAlgorithms();
    std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
    if ((unsigned int)nAlgo != algoNames.size()) {
      std::cerr << "Inconsistency. Check software trigger." << std::endl;
    }

    fData->swtrigger_name.resize(nAlgo);
    fData->swtrigger_triggered.resize(nAlgo);

    for (int trigger = 0; trigger < nAlgo; trigger++){
      fData->swtrigger_name[trigger]      = algoNames[trigger];
      fData->swtrigger_triggered[trigger] = softwareTriggerHandle->passedAlgo(algoNames[trigger]);
    }
  } // swtrigger


  //*****************************
  //
  // Common Optical Filter
  //
  //***************************** 
  if(fSaveOpticalFilterInfo){
    art::Handle<uboone::UbooneOpticalFilter> opticalFilterHandle;
    evt.getByLabel(fOpticalFilterLabel,opticalFilterHandle);

    if(!opticalFilterHandle.isValid() || opticalFilterHandle.failedToGet()){
      std::cerr << "Failed to get optical filter data product with label " << fOpticalFilterLabel << std::endl;
    }
    else{
      fData->opfilter_pebeam = opticalFilterHandle->PE_Beam();
      fData->opfilter_peveto = opticalFilterHandle->PE_Veto();
      fData->opfilter_pmtmaxfrac = opticalFilterHandle->PMT_MaxFraction();
    }
  }


  //  std::cout<<detProp.NumberTimeSamples()<<" "<<detProp.ReadOutWindowSize()<<std::endl;
  //  std::cout<<geom->DetHalfHeight()*2<<" "<<geom->DetHalfWidth()*2<<" "<<geom->DetLength()<<std::endl;
  //  std::cout<<geom->Nwires(0)<<" "<<geom->Nwires(1)<<" "<<geom->Nwires(2)<<std::endl;

  //hit information
  if (fSaveHitInfo){
    fData->no_hits = (int) NHits;
    fData->no_hits_stored = TMath::Min( (int) NHits, (int) kMaxHits);    
    if (NHits > kMaxHits) {
      // got this error? consider increasing kMaxHits
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NHits
					  << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
    }
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
      fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_charge[i]  = hitlist[i]->Integral();
      fData->hit_ph[i]  = hitlist[i]->PeakAmplitude();
      fData->hit_startT[i] = hitlist[i]->PeakTimeMinusRMS();
      fData->hit_endT[i] = hitlist[i]->PeakTimePlusRMS();
      fData->hit_rms[i] = hitlist[i]->RMS();
      fData->hit_goodnessOfFit[i] = hitlist[i]->GoodnessOfFit();
      fData->hit_multiplicity[i] = hitlist[i]->Multiplicity();

      //when the size of simIDEs is zero, the above function throws an exception
      //and crashes, so check that the simIDEs have non-zero size before 
      //extracting hit true XYZ from simIDEs
      
      /*
	for (unsigned int it=0; it<fTrackModuleLabel.size();++it){
        art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[it]);
        if (fmtk.at(i).size()!=0){
	hit_trkid[it][i] = fmtk.at(i)[0]->ID();
        }
        else
	hit_trkid[it][i] = 0;
	}
      */
      
      
      if (fSaveCalWireInfo){
      //Hit to CalWire information	    
      art::FindManyP<recob::Wire> fmwire(hitListHandle,evt,fHitsModuleLabel);
      
      if (hitlist[i]->WireID().Plane==2){
          int dataSize = fmwire.at(i)[0]->NSignal();
	  int t0 = hitlist[i]->PeakTime() - fCalWireROI*(hitlist[i]->RMS());
	  if (t0<0) t0 = 0;
	  int t1 = hitlist[i]->PeakTime() + fCalWireROI*(hitlist[i]->RMS());
	  if (t1>=dataSize) t1 = dataSize-1;
	  	  
          std::vector<float> signal(fmwire.at(i)[0]->Signal());
      
          fData->calwire_ph[i]    = -1.0;
          fData->calwire_peakT[i] = -1.0;
          for (int j = t0; j<=t1; ++j){
             if (signal[j]>fData->calwire_ph[i]){
               fData->calwire_ph[i]     = signal[j];
               fData->calwire_peakT[i]  = j;
             }   
          }
	  
	  fData->calwire_charge[i] = 0.0;  
	  fData->calwire_fwhh[i] = 0.0;  
      
          double mean_t = 0.0; 
          double mean_t2 = 0.0;
          for (int j = t0; j<=t1; ++j){
	    if (signal[j]>=0.5*fData->calwire_ph[i]){
		++fData->calwire_fwhh[i];
	    } 
            if (signal[j]>=(fCalWireThresh/100)*fData->calwire_ph[i]){
              fData->calwire_charge[i] += signal[j];
              mean_t  += double(j)*(signal[j]);
              mean_t2 += double(j)*double(j)*(signal[j]);
            }         
          }
          mean_t/=fData->calwire_charge[i];
          mean_t2/=fData->calwire_charge[i];
          fData->calwire_rms[i] = sqrt(mean_t2-mean_t*mean_t);
       }
     }    
     
      if (fSaveRawDigitInfo){
	//Hit to RawDigit information	       
	art::FindManyP<raw::RawDigit> fmrd(hitListHandle,evt,fHitsModuleLabel);
	if (hitlist[i]->WireID().Plane==2)
	  {
	    int dataSize = fmrd.at(i)[0]->Samples();
	    short ped = fmrd.at(i)[0]->GetPedestal();
     
	    std::vector<short> rawadc(dataSize);
	    raw::Uncompress(fmrd.at(i)[0]->ADCs(), rawadc, fmrd.at(i)[0]->Compression());
	    int t0 = hitlist[i]->PeakTime() - fRawDigitROI*(hitlist[i]->RMS());
	    if (t0<0) t0 = 0;
	    int t1 = hitlist[i]->PeakTime() + fRawDigitROI*(hitlist[i]->RMS());
	    if (t1>=dataSize) t1 = dataSize-1;
	    fData->rawD_ph[i] = -1;
	    fData->rawD_peakT[i] = -1;
	    for (int j = t0; j<=t1; ++j){
	      if (rawadc[j]-ped>fData->rawD_ph[i]){
		fData->rawD_ph[i] = rawadc[j]-ped;
		fData->rawD_peakT[i] = j;
	      }
	    }
	    fData->rawD_charge[i] = 0;
	    fData->rawD_fwhh[i] = 0;
	    double mean_t = 0.0;
	    double mean_t2 = 0.0;
	    for (int j = t0; j<=t1; ++j){
	      if (rawadc[j]-ped>=0.5*fData->rawD_ph[i]){
		++fData->rawD_fwhh[i];
	      }
	      if (rawadc[j]-ped>=(fRawDigitThresh/100)*fData->rawD_ph[i]){
		fData->rawD_charge[i] += rawadc[j]-ped;
		mean_t += (double)j*(rawadc[j]-ped);
		mean_t2 += (double)j*(double)j*(rawadc[j]-ped);
	      }
	    }
	    mean_t/=fData->rawD_charge[i];
	    mean_t2/=fData->rawD_charge[i];
	    fData->rawD_rms[i] = sqrt(mean_t2-mean_t*mean_t);
	  }   }
      
        //Hit to SimChannel information
        if (isMC && fSaveSimChannelInfo)
        {
            std::vector<const sim::SimChannel*> fSimChannels;
            evt.getView(fLArG4ModuleLabel, fSimChannels);
            
            const sim::SimChannel* chan = 0;
            for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
                if(fSimChannels[sc]->Channel() == hitlist[i]->Channel())
                    chan = fSimChannels[sc];
            }
            if (chan){
                auto const& tdcidemap = chan->TDCIDEMap();
                int k=-1;
                std::vector<double> elec(tdcidemap.size(),0.);
                std::vector<int> tdc(tdcidemap.size(),0.);
                for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
                    k++;
                    tdc[k]=(*mapitr).first;
                    const std::vector<sim::IDE> idevec = (*mapitr).second;
                    double nelec=0;
                    for(size_t iv = 0; iv < idevec.size(); ++iv){
                        nelec += idevec[iv].numElectrons;
                    }
                    elec[k] = nelec;
                }
                fData->sim_ph[i] = -1;
                fData->sim_tdc[i] = -1;
                for(unsigned int f=0;f<tdcidemap.size();f++){
                    if (elec[f]>fData->sim_ph[i]){
                        fData->sim_ph[i] = elec[f];
                        fData->sim_tdc[i] = tdc[f];
                    }
                }
                fData->sim_charge[i] = 0;
                fData->sim_fwhh[i] = 0;
                double mean_t = 0;
                double mean_t2 = 0;
                for (unsigned int f = 0; f<tdcidemap.size();f++){
                    if (elec[f]>=0.5*fData->sim_ph[i]){
                        ++fData->sim_fwhh[i];
                    }
                    if (elec[f]>=0.1*fData->sim_ph[i]){
                        fData->sim_charge[i]+= elec[f];
                        mean_t+= double(tdc[f])*elec[f];
                        mean_t2+= double(tdc[f])*double(tdc[f])*elec[f];
                    }
                }
                mean_t/=fData->sim_charge[i];
                mean_t2/=fData->sim_charge[i];
                fData->sim_rms[i] = sqrt(mean_t2-mean_t*mean_t);
            }
       } 
      
        if (isMC&&!isCosmics)
        {
            std::vector<sim::TrackIDE> trackIDEVec = fMCTruthMatching->HitToTrackID(clockData, hitlist[i]);
            fData -> hit_nelec[i] = 0;
            fData -> hit_energy[i] = 0;
            
            for(const auto& trackIDE : trackIDEVec)
            {
                fData -> hit_nelec[i]  += trackIDE.numElectrons;
                fData -> hit_energy[i] += trackIDE.energy;
            }
        }
    }

    if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      //Find tracks associated with hits
      art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[0]);
      for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
        if (fmtk.isValid()){
	  if (fmtk.at(i).size()!=0){
	    fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
	    fData->hit_trkKey[i] = fmtk.at(i)[0].key();

	  }
	  else
	    fData->hit_trkid[i] = -1;
        }
      }
    }
    
    //In the case of ClusterCrawler or linecluster, use "linecluster or clustercrawler" as HitModuleLabel.
    //using cchit will not make this association. In the case of gaushit, just use gaushit
    //Not initializing clusterID to -1 since some clustering algorithms assign negative IDs!
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      //Find clusters associated with hits
      art::FindManyP<recob::Cluster> fmcl(hitListHandle,evt,fClusterModuleLabel);
      for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
        if (fmcl.isValid()){
      	  if (fmcl.at(i).size()!=0){
	    fData->hit_clusterid[i] = fmcl.at(i)[0]->ID();
	    fData->hit_clusterKey[i] = fmcl.at(i)[0].key();
  	  }
        }
      }
    }
  }// end (fSaveHitInfo)
  
  if (fSaveRawWaveFormInfo){
    art::Handle< std::vector<raw::RawDigit> > rawdigit;
       if (! evt.getByLabel(fDigitModuleLabel, rawdigit)) {
           std::cout << "WARNING: no label " << fDigitModuleLabel << std::endl;
           return;
       }
       std::vector< art::Ptr<raw::RawDigit> >  wires;
       art::fill_ptr_vector(wires, rawdigit);

    AnalysisTreeDataStruct::RawDataStruct& RawData = fData->GetRawData();
    size_t NRawChannels = wires.size();
    
    RawData.SetMaxRawChannels(std::max(NRawChannels, (size_t) 1));
    RawData.Clear(); // clear all the data

    RawData.raw_nchannels = (short) NRawChannels;

    // now set the tree addresses to the newly allocated memory;
    // this creates the tree branches in case they are not there yet
    SetRawAddress();
      
    if (NRawChannels > RawData.GetMaxRawChannels()) {
      mf::LogError("AnalysisTree:limits") << "event has " << NRawChannels
                   << " raw channels, only "
                   << RawData.GetMaxRawChannels() << " stored in tree";
    }
    
    int irch=0;
    for (auto const& wire: wires) {
        RawData.raw_channelId[irch]    = wire->Channel();
    	size_t nSamples = wire->Samples();
    	std::vector<Short_t> uncompressed(nSamples);
	if (nSamples > RawData.GetMaxTicksPerRawChannel(irch)) {
	      // if you get this error, you'll have to increase kMaxTicks
	      mf::LogError("AnalysisTree:limits")
		<< "the " << fDigitModuleLabel << " RawChannel #" << irch
		<< " has " << nSamples << " in raw data channel"
		<<", only "
		<< RawData.GetMaxTicksPerRawChannel(irch) << " stored in tree";
	    }	
    	raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());
	for (size_t j = 0; j<nSamples && j < RawData.GetMaxTicksPerRawChannel(irch); ++j){	   
	    RawData.raw_wf[irch][j] = uncompressed[j];
	}    	
    	irch++;
    }
  }   
  
  if (fSaveCalibWaveFormInfo){
     art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! evt.getByLabel(fCalDataModuleLabel, wires_handle)) {
        std::cout << "WARNING: no label " << fCalDataModuleLabel << std::endl;
        return;
    }
    
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    AnalysisTreeDataStruct::CalibDataStruct& CalibData = fData->GetCalibData();
    size_t NCalibChannels = wires.size();
    
    CalibData.SetMaxCalibChannels(std::max(NCalibChannels, (size_t) 1));
    CalibData.Clear(); // clear all the data

    CalibData.calib_nchannels = (short) NCalibChannels;

    // now set the tree addresses to the newly allocated memory;
    // this creates the tree branches in case they are not there yet
    SetCalibAddress();
      
    if (NCalibChannels > CalibData.GetMaxCalibChannels()) {
      mf::LogError("AnalysisTree:limits") << "event has " << NCalibChannels
                   << " calib channels, only "
                   << CalibData.GetMaxCalibChannels() << " stored in tree";
    }
    
    int icch=0;
    for (auto const& wire: wires) {
        CalibData.calib_channelId[icch]    = wire->Channel();
	std::vector<float> calibwf = wire->Signal();
	
    	size_t nSamples = 9600;
	if (nSamples > CalibData.GetMaxTicksPerCalibChannel(icch)) {
	      // if you get this error, you'll have to increase kMaxTicks
	      mf::LogError("AnalysisTree:limits")
		<< "the " << fDigitModuleLabel << " CalibChannel #" << icch
		<< " has " << nSamples << " in calib data channel"
		<<", only "
		<< CalibData.GetMaxTicksPerCalibChannel(icch) << " stored in tree";
	    }	
	for (size_t j = 0; j<nSamples && j < CalibData.GetMaxTicksPerCalibChannel(icch); ++j){	   
	    CalibData.calib_wf[icch][j] = calibwf[j];
	}    	
    	icch++;
    }
  }   



  if (fSaveClusterInfo){
    fData->nclusters = (int) NClusters;
    if (NClusters > kMaxClusters){
      // got this error? consider increasing kMaxClusters
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NClusters
					  << " clusters, only kMaxClusters=" << kMaxClusters << " stored in tree";
    }
    for(unsigned int ic=0; ic<NClusters;++ic){//loop over clusters
      art::Ptr<recob::Cluster> clusterholder(clusterListHandle, ic);
      const recob::Cluster& cluster = *clusterholder;
      fData->clusterId[ic] = cluster.ID();
      fData->clusterView[ic] = cluster.View();
      fData->cluster_isValid[ic] = cluster.isValid();
      fData->cluster_StartCharge[ic] = cluster.StartCharge();
      fData->cluster_StartAngle[ic] = cluster.StartAngle();
      fData->cluster_EndCharge[ic] = cluster.EndCharge();
      fData->cluster_EndAngle[ic] = cluster.EndAngle();
      fData->cluster_Integral[ic] = cluster.Integral();
      fData->cluster_IntegralAverage[ic] = cluster.IntegralAverage();
      fData->cluster_SummedADC[ic] = cluster.SummedADC();
      fData->cluster_SummedADCaverage[ic] = cluster.SummedADCaverage();
      fData->cluster_MultipleHitDensity[ic] = cluster.MultipleHitDensity();
      fData->cluster_Width[ic] = cluster.Width();
      fData->cluster_NHits[ic] = cluster.NHits();
      fData->cluster_StartWire[ic] = cluster.StartWire();
      fData->cluster_StartTick[ic] = cluster.StartTick();
      fData->cluster_EndWire[ic] = cluster.EndWire();
      fData->cluster_EndTick[ic] = cluster.EndTick();
	 
      //Cosmic Tagger information for cluster
      art::FindManyP<anab::CosmicTag> fmcct(clusterListHandle,evt,fCosmicClusterTaggerAssocLabel);
      if (fmcct.isValid()){          
	fData->cluncosmictags_tagger[ic]     = fmcct.at(ic).size();
	if (fmcct.at(ic).size()>0){
	  if(fmcct.at(ic).size()>1)
	    std::cerr << "\n Warning : more than one cosmic tag per cluster in module! assigning the first tag to the cluster" << fCosmicClusterTaggerAssocLabel;
	  fData->clucosmicscore_tagger[ic] = fmcct.at(ic).at(0)->CosmicScore();
	  fData->clucosmictype_tagger[ic] = fmcct.at(ic).at(0)->CosmicType();
	}
      }
    }//end loop over clusters
  }//end fSaveClusterInfo
	 
  
  // Declare object-ID-to-PFParticleID maps so we can assign hasPFParticle and PFParticleID to the tracks, showers, vertices.
  std::map<Short_t, Short_t> trackIDtoPFParticleIDMap, vertexIDtoPFParticleIDMap, showerIDtoPFParticleIDMap;
  
  //Save PFParticle information
  if (fSavePFParticleInfo){
    AnalysisTreeDataStruct::PFParticleDataStruct& PFParticleData = fData->GetPFParticleData();
    size_t NPFParticles = pfparticlelist.size();
    
    PFParticleData.SetMaxPFParticles(std::max(NPFParticles, (size_t) 1));
    PFParticleData.Clear(); // clear all the data

    PFParticleData.nPFParticles = (short) NPFParticles;

    // now set the tree addresses to the newly allocated memory;
    // this creates the tree branches in case they are not there yet
    SetPFParticleAddress();
      
    if (NPFParticles > PFParticleData.GetMaxPFParticles()) {
      mf::LogError("AnalysisTree:limits") << "event has " << NPFParticles
                   << " PFParticles, only "
                   << PFParticleData.GetMaxPFParticles() << " stored in tree";
    }
    
    lar_pandora::PFParticleVector neutrinoPFParticles;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticlelist, neutrinoPFParticles);
    PFParticleData.pfp_numNeutrinos = neutrinoPFParticles.size();
    
    for (size_t i = 0; i < std::min(neutrinoPFParticles.size(), (size_t)kMaxNPFPNeutrinos); ++i) {
      PFParticleData.pfp_neutrinoIDs[i] = neutrinoPFParticles[i]->Self();
    }
    
    if (neutrinoPFParticles.size() > kMaxNPFPNeutrinos)
      std::cerr << "Warning: there were more than 5 reconstructed PFParticle neutrinos; only the first 5 being stored in tree" << std::endl;
    
    // Get a PFParticle-to-vertex map.
    lar_pandora::VertexVector allPfParticleVertices;
    lar_pandora::PFParticlesToVertices pfParticleToVertexMap;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPFParticleModuleLabel, allPfParticleVertices, pfParticleToVertexMap);
    
    // Get a PFParticle-to-track map.
    lar_pandora::TrackVector allPfParticleTracks;
    lar_pandora::PFParticlesToTracks pfParticleToTrackMap;
    lar_pandora::LArPandoraHelper::CollectTracks(evt, fPFParticleModuleLabel, allPfParticleTracks, pfParticleToTrackMap);
    
    // Get a PFParticle-to-shower map.
    lar_pandora::ShowerVector allPfParticleShowers;
    lar_pandora::PFParticlesToShowers pfParticleToShowerMap;
    lar_pandora::LArPandoraHelper::CollectShowers(evt, fPFParticleModuleLabel, allPfParticleShowers, pfParticleToShowerMap);
    
    for (size_t i = 0; i < NPFParticles && i < PFParticleData.GetMaxPFParticles() ; ++i){
      PFParticleData.pfp_selfID[i] = pfparticlelist[i]->Self();
      PFParticleData.pfp_isPrimary[i] = (Short_t)pfparticlelist[i]->IsPrimary();
      PFParticleData.pfp_numDaughters[i] = pfparticlelist[i]->NumDaughters();
      PFParticleData.pfp_parentID[i] = pfparticlelist[i]->Parent();
      PFParticleData.pfp_pdgCode[i] = pfparticlelist[i]->PdgCode();
      PFParticleData.pfp_isNeutrino[i] = lar_pandora::LArPandoraHelper::IsNeutrino(pfparticlelist[i]);
         
      // Set the daughter IDs.
      std::vector<size_t> daughterIDs = pfparticlelist[i]->Daughters();

      if (daughterIDs.size() > kMaxNDaughtersPerPFP)
        std::cerr << "Warning: the number of PFParticle daughters exceeds the maximum number. Only the first " << kMaxNDaughtersPerPFP << " will be used" << std::endl;

      for (size_t j = 0; j < daughterIDs.size() && j < kMaxNDaughtersPerPFP; ++j)
        PFParticleData.pfp_daughterIDs[i][j] = daughterIDs[j];
         
      // Set the vertex ID.
      auto vertexMapIter = pfParticleToVertexMap.find(pfparticlelist[i]);
      if (vertexMapIter != pfParticleToVertexMap.end()) {
          lar_pandora::VertexVector pfParticleVertices = vertexMapIter->second;
          
          if (pfParticleVertices.size() > 1)
            std::cerr << "Warning: there was more than one vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << ", storing only one" << std::endl;
          
          if (pfParticleVertices.size() > 0) {
            PFParticleData.pfp_vertexID[i] = pfParticleVertices.at(0)->ID();
            vertexIDtoPFParticleIDMap.insert(std::make_pair(pfParticleVertices.at(0)->ID(), pfparticlelist[i]->Self()));
          }
      }
      else
        std::cerr << "Warning: there was no vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
      
      if (lar_pandora::LArPandoraHelper::IsTrack(pfparticlelist[i])){
        PFParticleData.pfp_isTrack[i] = 1;
        
        // Set the track ID.
        auto trackMapIter = pfParticleToTrackMap.find(pfparticlelist[i]);
        if (trackMapIter != pfParticleToTrackMap.end()) {
            lar_pandora::TrackVector pfParticleTracks = trackMapIter->second;
            
            if (pfParticleTracks.size() > 1)
              std::cerr << "Warning: there was more than one track found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
            
            if (pfParticleTracks.size() > 0) {
              PFParticleData.pfp_trackID[i] = pfParticleTracks.at(0)->ID();
              trackIDtoPFParticleIDMap.insert(std::make_pair(pfParticleTracks.at(0)->ID(), pfparticlelist[i]->Self()));
            }
        }
        else
          std::cerr << "Warning: there was no track found for track-like PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
      }
      else
        PFParticleData.pfp_isTrack[i] = 0;
      
      if (lar_pandora::LArPandoraHelper::IsShower(pfparticlelist[i])) {
        PFParticleData.pfp_isShower[i] = 1;
        // Set the shower ID.
        auto showerMapIter = pfParticleToShowerMap.find(pfparticlelist[i]);
        if (showerMapIter != pfParticleToShowerMap.end()) {
          lar_pandora::ShowerVector pfParticleShowers = showerMapIter->second;
          
          if (pfParticleShowers.size() > 1)
            std::cerr << "Warning: there was more than one shower found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
          
          if (pfParticleShowers.size() > 0) {
            PFParticleData.pfp_showerID[i] = pfParticleShowers.at(0)->ID();
            showerIDtoPFParticleIDMap.insert(std::make_pair(pfParticleShowers.at(0)->ID(), pfparticlelist[i]->Self()));
          }
        }
        else
          std::cerr << "Warning: there was no shower found for shower-like PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
      }
      else
        PFParticleData.pfp_isShower[i] = 0;

      // Set the cluster IDs.
      auto clusterMapIter = pfParticleToClusterMap.find(pfparticlelist[i]);
      if (clusterMapIter != pfParticleToClusterMap.end()) {
          lar_pandora::ClusterVector pfParticleClusters = clusterMapIter->second;
          PFParticleData.pfp_numClusters[i] = pfParticleClusters.size();

          if (pfParticleClusters.size() > kMaxNClustersPerPFP)
            std::cerr << "Warning: the number of PFParticle clusters exceeds the maximum number. Only the first " << kMaxNClustersPerPFP << " will be used" << std::endl;           

          for (size_t j = 0; j < pfParticleClusters.size() && j < kMaxNClustersPerPFP; ++j) 
            PFParticleData.pfp_clusterIDs[i][j] = pfParticleClusters[j]->ID();
      }
      //else
      //  std::cerr << "Warning: there were no clusters found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
    }
  } // if fSavePFParticleInfo
  
  if (fSaveShowerInfo){
    
    // fill data from all the shower algorithms
    for (size_t iShowerAlgo = 0; iShowerAlgo < NShowerAlgos; ++iShowerAlgo) {
      AnalysisTreeDataStruct::ShowerDataStruct& ShowerData
        = fData->GetShowerData(iShowerAlgo);
      std::vector<recob::Shower> const* pShowers = showerList[iShowerAlgo];
      art::Handle< std::vector<recob::Shower> > showerHandle = showerListHandle[iShowerAlgo];
      
      if (pShowers) FillShowers(ShowerData, *pShowers, fSavePFParticleInfo, showerIDtoPFParticleIDMap);
      else ShowerData.MarkMissing(fTree); // tree should reflect lack of data
    } // for iShowerAlgo
    
  } // if fSaveShowerInfo

  //track information for multiple trackers
  if (fSaveTrackInfo) {
    for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
      AnalysisTreeDataStruct::TrackDataStruct& TrackerData = fData->GetTrackerData(iTracker);
    
      size_t NTracks = tracklist[iTracker].size();
      // allocate enough space for this number of tracks (but at least for one of them!)
      TrackerData.SetMaxTracks(std::max(NTracks, (size_t) 1));
      TrackerData.Clear(); // clear all the data
    
      TrackerData.ntracks = (int) NTracks;
    
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetTrackerAddresses(iTracker);
      if (NTracks > TrackerData.GetMaxTracks()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnalysisTree:limits") << "event has " << NTracks
					    << " " << fTrackModuleLabel[iTracker] << " tracks, only "
					    << TrackerData.GetMaxTracks() << " stored in tree";
      }
    
      //call the track momentum algorithm that gives you momentum based on track range
      // - change the minimal track length requirement to 50 cm
      trkf::TrackMomentumCalculator trkm{50};

      for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks
      
        //save t0 from anode/cathode piercing tracks
	art::FindManyP<anab::T0> fmacpt0(trackListHandle[iTracker],evt,fT0RecoAnodeCathodePiercingLabel[iTracker]);
	if (fmacpt0.isValid()){          
	  if(fmacpt0.at(iTrk).size()>0){
	    if(fmacpt0.at(iTrk).size()>1)
	      std::cerr << "\n Warning : more than one t0 tag per track in module! assigning the first tag to the track" << fT0RecoAnodeCathodePiercingLabel[iTracker];
	    TrackerData.trkACpierceT0[iTrk] = fmacpt0.at(iTrk).at(0)->Time();
	  }   
        }
	
	//save t0 from reconstructed flash track matching for every track
        art::FindManyP<anab::T0> fmt0(trackListHandle[iTracker],evt,fFlashT0FinderLabel[iTracker]);
        if (fmt0.isValid()){          
	  if(fmt0.at(iTrk).size()>0){
	    if(fmt0.at(iTrk).size()>1)
	      std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashT0FinderLabel[iTracker];
	    TrackerData.trkflashT0[iTrk] = fmt0.at(iTrk).at(0)->Time();
	  }   
        }
	
	//save t0 from reconstructed flash track matching for every track
        art::FindManyP<anab::T0> fmmct0(trackListHandle[iTracker],evt,fMCT0FinderLabel[iTracker]);
        if (fmmct0.isValid()){          
	  if(fmmct0.at(iTrk).size()>0){
	    if(fmmct0.at(iTrk).size()>1)
	      std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the cluster" << fMCT0FinderLabel[iTracker];
	    TrackerData.trktrueT0[iTrk] = fmmct0.at(iTrk).at(0)->Time();
	  }   
        }
        
        //Cosmic Tagger information
        art::FindManyP<anab::CosmicTag> fmct(trackListHandle[iTracker],evt,fCosmicTaggerAssocLabel[iTracker]);
        if (fmct.isValid()){          
          TrackerData.trkncosmictags_tagger[iTrk]     = fmct.at(iTrk).size();
          if (fmct.at(iTrk).size()>0){
            if(fmct.at(iTrk).size()>1)
              std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fCosmicTaggerAssocLabel[iTracker];
            TrackerData.trkcosmicscore_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicScore();
            TrackerData.trkcosmictype_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicType();
          }
        }

        //Containment Tagger information
        art::FindManyP<anab::CosmicTag> fmcnt(trackListHandle[iTracker],evt,fContainmentTaggerAssocLabel[iTracker]);
        if (fmcnt.isValid()){
          TrackerData.trkncosmictags_containmenttagger[iTrk]     = fmcnt.at(iTrk).size();
          if (fmcnt.at(iTrk).size()>0){
            if(fmcnt.at(iTrk).size()>1)
              std::cerr << "\n Warning : more than one containment tag per track in module! assigning the first tag to the track" << fContainmentTaggerAssocLabel[iTracker];
            TrackerData.trkcosmicscore_containmenttagger[iTrk] = fmcnt.at(iTrk).at(0)->CosmicScore();
            TrackerData.trkcosmictype_containmenttagger[iTrk] = fmcnt.at(iTrk).at(0)->CosmicType();
          }
        }

        //Flash match compatibility information
        //Unlike CosmicTagger, Flash match doesn't assign a cosmic tag for every track. For those tracks, AnalysisTree initializes them with -9999 or -99999
        art::FindManyP<anab::CosmicTag> fmbfm(trackListHandle[iTracker],evt,fFlashMatchAssocLabel[iTracker]);
        if (fmbfm.isValid()){  
          TrackerData.trkncosmictags_flashmatch[iTrk] = fmbfm.at(iTrk).size();
          if (fmbfm.at(iTrk).size()>0){
            if(fmbfm.at(iTrk).size()>1) 
              std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashMatchAssocLabel[iTracker];
	    TrackerData.trkcosmicscore_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicScore();
	    TrackerData.trkcosmictype_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicType();
	    //std::cout<<"\n"<<evt.event()<<"\t"<<iTrk<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicScore()<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicType();
          }
        }
     			 	   
        art::Ptr<recob::Track> ptrack(trackListHandle[iTracker], iTrk);
        const recob::Track& track = *ptrack;
      
        TVector3 pos, dir_start, dir_end, end;        

        double tlen = 0., mom = 0.;
        int TrackID = -1;
      
        int ntraj = track.NumberTrajectoryPoints();
          if (ntraj > 0) {
            pos       = track.Vertex<TVector3>();
            dir_start = track.VertexDirection<TVector3>();
            dir_end   = track.EndDirection<TVector3>();
            end       = track.End<TVector3>();

            tlen        = track.Length(); //length(track);
            if(track.NumberTrajectoryPoints() > 0)
              mom = track.VertexMomentum();
      
          TrackID = track.ID();
          double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
          double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
          double dpos = bdist(pos);
          double dend = bdist(end);
        
          TrackerData.trkId[iTrk]                 = TrackID;
          TrackerData.trkstartx[iTrk]             = pos.X();
          TrackerData.trkstarty[iTrk]             = pos.Y();
          TrackerData.trkstartz[iTrk]             = pos.Z();
          TrackerData.trkstartd[iTrk]		  = dpos;
          TrackerData.trkendx[iTrk]		  = end.X();
          TrackerData.trkendy[iTrk]		  = end.Y();
          TrackerData.trkendz[iTrk]		  = end.Z();
          TrackerData.trkendd[iTrk]		  = dend;
          TrackerData.trktheta[iTrk]		  = dir_start.Theta();
          TrackerData.trkphi[iTrk]		  = dir_start.Phi();
          TrackerData.trkstartdcosx[iTrk]	  = dir_start.X();
          TrackerData.trkstartdcosy[iTrk]	  = dir_start.Y();
          TrackerData.trkstartdcosz[iTrk]	  = dir_start.Z();
          TrackerData.trkenddcosx[iTrk] 	  = dir_end.X();
          TrackerData.trkenddcosy[iTrk] 	  = dir_end.Y();
          TrackerData.trkenddcosz[iTrk] 	  = dir_end.Z();
          TrackerData.trkthetaxz[iTrk]  	  = theta_xz;
          TrackerData.trkthetayz[iTrk]  	  = theta_yz;
          TrackerData.trkmom[iTrk]		  = mom;
          TrackerData.trklen[iTrk]		  = tlen;
          TrackerData.trkmomrange[iTrk] 	  = trkm.GetTrackMomentum(tlen,13);
          //TrackerData.trkmommschi2[iTrk]	  = trkm.GetMomentumMultiScatterChi2(ptrack);
          //TrackerData.trkmommsllhd[iTrk]	  = trkm.GetMomentumMultiScatterLLHD(ptrack);
	  TrackerData.trkntraj[iTrk]              = ntraj;
          if (tracklist[iTracker].size() == mcsfitlist[iTracker].size()){
            TrackerData.trkmcsfwdmom[iTrk]        = mcsfitlist[iTracker][iTrk]->fwdMomentum();
            TrackerData.trkmcsfwdll[iTrk]         = mcsfitlist[iTracker][iTrk]->fwdLogLikelihood();
            TrackerData.trkmcsfwderr[iTrk]        = mcsfitlist[iTracker][iTrk]->fwdMomUncertainty();
            TrackerData.trkmcsbwdmom[iTrk]        = mcsfitlist[iTracker][iTrk]->bwdMomentum();
            TrackerData.trkmcsbwdll[iTrk]         = mcsfitlist[iTracker][iTrk]->bwdLogLikelihood();
            TrackerData.trkmcsbwderr[iTrk]        = mcsfitlist[iTracker][iTrk]->bwdMomUncertainty();
          }
          if (fSavePFParticleInfo) {
            auto mapIter = trackIDtoPFParticleIDMap.find(TrackID);
            if (mapIter != trackIDtoPFParticleIDMap.end()) {
                // This track has a corresponding PFParticle.
                TrackerData.trkhasPFParticle[iTrk] = 1;
                TrackerData.trkPFParticleID[iTrk] = mapIter->second;
            }
            else 
                TrackerData.trkhasPFParticle[iTrk] = 0;
          }
          
        } // if we have trajectory

	// find vertices associated with this track
	/*
	  art::FindMany<recob::Vertex> fmvtx(trackListHandle[iTracker], evt, fVertexModuleLabel[iTracker]);
	  if(fmvtx.isValid()) {
	  std::vector<const recob::Vertex*> verts = fmvtx.at(iTrk);
	  // should have two at most
	  for(size_t ivx = 0; ivx < verts.size(); ++ivx) {
          verts[ivx]->XYZ(xyz);
          // find the vertex in TrackerData to get the index
          short theVtx = -1;
          for(short jvx = 0; jvx < TrackerData.nvtx; ++jvx) {
	  if(TrackerData.vtx[jvx][2] == xyz[2]) {
	  theVtx = jvx;
	  break;
	  }
          } // jvx
          // decide if it should be assigned to the track Start or End.
          // A simple dz test should suffice
          if(fabs(xyz[2] - TrackerData.trkstartz[iTrk]) < 
	  fabs(xyz[2] - TrackerData.trkendz[iTrk])) {
	  TrackerData.trksvtxid[iTrk] = theVtx;
          } else {
	  TrackerData.trkevtxid[iTrk] = theVtx;
          }
	  } // vertices
	  } // fmvtx.isValid()
	*/


	/* //commented out because now have several Vertices
	   Float_t minsdist = 10000;
	   Float_t minedist = 10000;
	   for (int ivx = 0; ivx < NVertices && ivx < kMaxVertices; ++ivx){
	   Float_t sdist = sqrt(pow(TrackerData.trkstartx[iTrk]-VertexData.vtxx[ivx],2)+
	   pow(TrackerData.trkstarty[iTrk]-VertexData.vtxy[ivx],2)+
	   pow(TrackerData.trkstartz[iTrk]-VertexData.vtxz[ivx],2));
	   Float_t edist = sqrt(pow(TrackerData.trkendx[iTrk]-VertexData.vtxx[ivx],2)+
	   pow(TrackerData.trkendy[iTrk]-VertexData.vtxy[ivx],2)+
	   pow(TrackerData.trkendz[iTrk]-VertexData.vtxz[ivx],2));
	   if (sdist<minsdist){
	   minsdist = sdist;
	   if (minsdist<10) TrackerData.trksvtxid[iTrk] = ivx;
	   }
	   if (edist<minedist){
	   minedist = edist;
	   if (minedist<10) TrackerData.trkevtxid[iTrk] = ivx;
	   }
	   }*/

	// find particle ID info
	art::FindMany<anab::ParticleID> fmpid(trackListHandle[iTracker], evt, fParticleIDModuleLabel[iTracker]);
	if(fmpid.isValid()) {
	  std::vector<const anab::ParticleID*> pids = fmpid.at(iTrk);
	  if (pids.size() == 0){
	    mf::LogError("AnalysisTree:limits")
	      << "No track-PID association found for " << fTrackModuleLabel[iTracker]
	      << " track " << iTrk << ". Not saving particleID information."; 
	  }
	  // Set dummy values
	  double pidpdg[3] = {-1,-1,-1};
	  double pidchi[3] = {99999.,99999.,99999.};
	  for (size_t ipid=0; ipid<pids.size(); ipid++){ 
	    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();

	    // Loop though AlgScoresVec and find the variables we want
	    for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
	      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

	      /*std::cout << "\n ParticleIDAlg " << AlgScore.fAlgName
			<< "\n -- Variable type: " << AlgScore.fVariableType
			<< "\n -- Track direction: " << AlgScore.fTrackDir
			<< "\n -- Assuming PDG: " << AlgScore.fAssumedPdg
			<< "\n -- Number of degrees of freedom: " << AlgScore.fNdf
			<< "\n -- Value: " << AlgScore.fValue
			<< "\n -- Using planeMask: " << AlgScore.fPlaneMask << " (plane " << UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) << ")" << std::endl;*/
	      
	      int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	      if (planenum<0 || planenum>2) continue;

	      if (AlgScore.fAlgName == "Chi2"){
		if (TMath::Abs(AlgScore.fAssumedPdg) == 13){ // chi2mu
		  TrackerData.trkpidchimu[iTrk][planenum] = AlgScore.fValue;
		  if (AlgScore.fValue<pidchi[planenum]){
		    pidchi[planenum] = AlgScore.fValue;
		    pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
		  }
		}
		else if (TMath::Abs(AlgScore.fAssumedPdg) == 2212){ // chi2pr
		  TrackerData.trkpidchipr[iTrk][planenum] = AlgScore.fValue;
		  if (AlgScore.fValue<pidchi[planenum]){
		    pidchi[planenum] = AlgScore.fValue;
		    pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
		  }
		}
		else if (TMath::Abs(AlgScore.fAssumedPdg) == 211){ // chi2pi
		  TrackerData.trkpidchipi[iTrk][planenum] = AlgScore.fValue;
		  if (AlgScore.fValue<pidchi[planenum]){
		    pidchi[planenum] = AlgScore.fValue;
		    pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
		  }
		}
		else if (TMath::Abs(AlgScore.fAssumedPdg) == 321){ // chi2ka
		  TrackerData.trkpidchika[iTrk][planenum] = AlgScore.fValue;
		  if (AlgScore.fValue<pidchi[planenum]){
		    pidchi[planenum] = AlgScore.fValue;
		    pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
		  }
		}
	      
	      }
	      else if (AlgScore.fVariableType==anab::kPIDA){
		TrackerData.trkpidpida[iTrk][planenum] = AlgScore.fValue;
	      }
	      
	    } // end loop though AlgScoresVec
	  } // end loop over pid[ipid]

	  // Finally, set min chi2
	  for (size_t planenum=0; planenum<3; planenum++){
	    TrackerData.trkpidchi[iTrk][planenum] = pidchi[planenum];
	    TrackerData.trkpidpdg[iTrk][planenum] = pidpdg[planenum];
	  }
	} // fmpid.isValid()
	
	art::FindMany<anab::Calorimetry> fmcal(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
	if (fmcal.isValid()){
	  std::vector<const anab::Calorimetry*> calos = fmcal.at(iTrk);
	  if (calos.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
	    // if you get this message, there is probably a bug somewhere since
	    // the calorimetry planes should be 3.
	    mf::LogError("AnalysisTree:limits")
	      << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
	      << " has " << calos.size() << " planes for calorimetry , only "
	      << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in tree";
	  }
	  for (size_t ical = 0; ical<calos.size(); ++ical){
	    if (!calos[ical]) continue;
	    if (!calos[ical]->PlaneID().isValid) continue;
	    int planenum = calos[ical]->PlaneID().Plane;
	    if (planenum<0||planenum>2) continue;
	    TrackerData.trkke[iTrk][planenum]    = calos[ical]->KineticEnergy();
	    TrackerData.trkrange[iTrk][planenum] = calos[ical]->Range();
	    //For now make the second argument as 13 for muons. 
	    TrackerData.trkpitchc[iTrk][planenum]= calos[ical] -> TrkPitchC();
	    const size_t NHits = calos[ical] -> dEdx().size();
	    TrackerData.ntrkhits[iTrk][planenum] = (int) NHits;
	    if (NHits > TrackerData.GetMaxHitsPerTrack(iTrk, planenum)) {
	      // if you get this error, you'll have to increase kMaxTrackHits
	      mf::LogError("AnalysisTree:limits")
		<< "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
		<< " has " << NHits << " hits on calorimetry plane #" << planenum
		<<", only "
		<< TrackerData.GetMaxHitsPerTrack(iTrk, planenum) << " stored in tree";
	    }
	    if (!isCosmics){
	      for(size_t iTrkHit = 0; iTrkHit < NHits && iTrkHit < TrackerData.GetMaxHitsPerTrack(iTrk, planenum); ++iTrkHit) {
		TrackerData.trkdedx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dEdx())[iTrkHit];
		TrackerData.trkdqdx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dQdx())[iTrkHit];
		TrackerData.trkresrg[iTrk][planenum][iTrkHit] = (calos[ical] -> ResidualRange())[iTrkHit];
		const auto& TrkPos = (calos[ical] -> XYZ())[iTrkHit];
		auto& TrkXYZ = TrackerData.trkxyz[iTrk][planenum][iTrkHit];
		TrkXYZ[0] = TrkPos.X();
		TrkXYZ[1] = TrkPos.Y();
		TrkXYZ[2] = TrkPos.Z();
	      } // for track hits
	    }
	  } // for calorimetry info
	  if(TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
	  else if(TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 1;
	  else if(TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
	  else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
	  else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0]) TrackerData.trkpidbestplane[iTrk] = 2;
	  else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
	  else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 2;
	} // if has calorimetry info

	//track truth information
	if (isMC){
	  //get the hits on each plane
	  art::FindManyP<recob::Hit>      fmht(trackListHandle[iTracker], evt, fTrackModuleLabel[iTracker]);
	  std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(iTrk);
	  std::vector< art::Ptr<recob::Hit> > hits[kNplanes];

	  for(size_t ah = 0; ah < allHits.size(); ++ah){
	    if (/* allHits[ah]->WireID().Plane >= 0 && */ // always true
		allHits[ah]->WireID().Plane <  3){
	      hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
	    }
	  }
	  for (size_t ipl = 0; ipl < 3; ++ipl){
	    double maxe = 0;
            HitsPurity(clockData, hits[ipl],TrackerData.trkidtruth[iTrk][ipl],TrackerData.trkpurtruth[iTrk][ipl],maxe);
	    //std::cout<<"\n"<<iTracker<<"\t"<<iTrk<<"\t"<<ipl<<"\t"<<TrackerData.trkidtruth[iTrk][ipl]<<"\t"<<TrackerData.trkpurtruth[iTrk][ipl]<<"\t"<<maxe;
	    if (TrackerData.trkidtruth[iTrk][ipl]>0){
	      const art::Ptr<simb::MCTruth> mc = fMCTruthMatching->TrackIDToMCTruth(TrackerData.trkidtruth[iTrk][ipl]);
	      TrackerData.trkorigin[iTrk][ipl] = mc->Origin();
	      const simb::MCParticle *particle = fMCTruthMatching->TrackIDToParticle(TrackerData.trkidtruth[iTrk][ipl]);

          double tote = 1.; //0;

	      TrackerData.trkpdgtruth[iTrk][ipl] = particle->PdgCode();
	      TrackerData.trkefftruth[iTrk][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
	      //std::cout<<"\n"<<trkpdgtruth[iTracker][iTrk][ipl]<<"\t"<<trkefftruth[iTracker][iTrk][ipl];
	    }
	  }

	  double maxe = 0;
          HitsPurity(clockData, allHits,TrackerData.trkg4id[iTrk],TrackerData.trkpurity[iTrk],maxe);
	  //std::cout<<"\n"<<fTrackModuleLabel[iTracker]<<"\t"<<iTrk<<"\t"<<"\t"<<TrackerData.trkg4id[iTrk]<<"\t"<<maxe<<std::endl;
	  if (TrackerData.trkg4id[iTrk]>0){
	    const art::Ptr<simb::MCTruth> mc = fMCTruthMatching->TrackIDToMCTruth(TrackerData.trkg4id[iTrk]);
	    TrackerData.trkorig[iTrk] = mc->Origin();
	  }
	  if (allHits.size()){
	    std::vector<art::Ptr<recob::Hit> > all_hits;
	    art::Handle<std::vector<recob::Hit> > hithandle;
	    float totenergy = 0;
	    if (evt.get(allHits[0].id(), hithandle)){
	      art::fill_ptr_vector(all_hits, hithandle);
	      for(size_t h = 0; h < all_hits.size(); ++h){

		art::Ptr<recob::Hit> hit = all_hits[h];
		std::vector<sim::IDE> ides;

		//fMCTruthMatching->HitToSimIDEs(hit,ides);
                std::vector<sim::TrackIDE> eveIDs = fMCTruthMatching->HitToEveID(clockData, hit);
		
		for(size_t e = 0; e < eveIDs.size(); ++e){
		  //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
		  if (eveIDs[e].trackID==TrackerData.trkg4id[iTrk]) totenergy += eveIDs[e].energy;
		}
	      }	      
	    }
	    if (totenergy) TrackerData.trkcompleteness[iTrk] = maxe/totenergy;
	  }
	}//end if (isMC)
      }//end loop over track
    }//end loop over track module labels
  }// end (fSaveTrackInfo) 
  
  /*trkf::TrackMomentumCalculator trkm;  
    std::cout<<"\t"<<trkm.GetTrackMomentum(200,2212)<<"\t"<<trkm.GetTrackMomentum(-10, 13)<<"\t"<<trkm.GetTrackMomentum(300,-19)<<"\n";
  */

  //Save Vertex information for multiple algorithms
  if (fSaveVertexInfo){
    for (unsigned int iVertexAlg=0; iVertexAlg < NVertexAlgos; ++iVertexAlg){
      AnalysisTreeDataStruct::VertexDataStruct& VertexData = fData->GetVertexData(iVertexAlg);

      size_t NVertices = vertexlist[iVertexAlg].size();

      VertexData.SetMaxVertices(std::max(NVertices, (size_t) 1));
      VertexData.Clear(); // clear all the data

      VertexData.nvtx = (short) NVertices;

      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetVertexAddresses(iVertexAlg);
      if (NVertices > VertexData.GetMaxVertices()) {
	// got this error? it might be a bug,
	// since we are supposed to have allocated enough space to fit all tracks
	mf::LogError("AnalysisTree:limits") << "event has " << NVertices
					    << " " << fVertexModuleLabel[iVertexAlg] << " vertices, only "
					    << VertexData.GetMaxVertices() << " stored in tree";
      }

      for (size_t i = 0; i < NVertices && i < kMaxVertices ; ++i){ //loop over vertices
	VertexData.vtxId[i] = vertexlist[iVertexAlg][i]->ID();
	Double_t xyz[3] = {};
	vertexlist[iVertexAlg][i] -> XYZ(xyz);
        VertexData.vtxx[i] = xyz[0];
        VertexData.vtxy[i] = xyz[1];
        VertexData.vtxz[i] = xyz[2];
        
        if (fSavePFParticleInfo) {
          auto mapIter = vertexIDtoPFParticleIDMap.find(vertexlist[iVertexAlg][i]->ID());
          if (mapIter != vertexIDtoPFParticleIDMap.end()) {
            // This vertex has a corresponding PFParticle.
            VertexData.vtxhasPFParticle[i] = 1;
            VertexData.vtxPFParticleID[i] = mapIter->second;
          }
          else 
            VertexData.vtxhasPFParticle[i] = 0;
        }
        
        // find PFParticle ID info
        art::FindMany<recob::PFParticle> fmPFParticle(vertexListHandle[iVertexAlg], evt, fPFParticleModuleLabel);
        if(fmPFParticle.isValid()) {
          std::vector<const recob::PFParticle*> pfparticles = fmPFParticle.at(i);
          if(pfparticles.size() > 1)
          std::cerr << "Warning: more than one associated PFParticle found for a vertex. Only one stored in tree." << std::endl;
          if (pfparticles.size() == 0)
          VertexData.vtxhasPFParticle[i] = 0;
          else {
            VertexData.vtxhasPFParticle[i] = 1;
            VertexData.vtxPFParticleID[i] = pfparticles.at(0)->Self();
          }
        } // fmPFParticle.isValid()
      }
    }
  }

  // Save neutrino vertex information for multiple algorithms
  if(fSavePandoraNuVertexInfo) {
    for (unsigned int iNuVertexAlg=0; iNuVertexAlg < NNuVertexAlgos; ++iNuVertexAlg){
      AnalysisTreeDataStruct::NeutrinoVertexDataStruct& NuVertexData = fData->GetNeutrinoVertexData(iNuVertexAlg);

      size_t NNuVertices = nuvertexlist[iNuVertexAlg].size();

      NuVertexData.SetMaxVertices(std::max(NNuVertices, (size_t) 1));
      NuVertexData.Clear(); // clear all the data

      NuVertexData.nvtx = (short) NNuVertices;

      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetNeutrinoVertexAddresses(iNuVertexAlg);
      if (NNuVertices > NuVertexData.GetMaxVertices()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all nu vertices
        mf::LogError("AnalysisTree:limits") << "event has " << NNuVertices
                                            << " " << fPandoraNuVertexModuleLabel[iNuVertexAlg] << " neutrino vertices, only "
                                            << NuVertexData.GetMaxVertices() << " stored in tree";
      }

      for (size_t i = 0; i < NNuVertices && i < kMaxVertices ; ++i){ //loop over vertices
        NuVertexData.vtxId[i] = nuvertexlist[iNuVertexAlg][i]->ID();
        Double_t xyz[3] = {};
        nuvertexlist[iNuVertexAlg][i] -> XYZ(xyz);
        NuVertexData.vtxx[i] = xyz[0];
        NuVertexData.vtxy[i] = xyz[1];
        NuVertexData.vtxz[i] = xyz[2];
        NuVertexData.vtxpdg[i] = nuvertexlistToPfp[iNuVertexAlg][i]->PdgCode();

        NuVertexData.vtxhasPFParticle[i] = 1;
        NuVertexData.vtxPFParticleID[i] = nuvertexlistToPfp[iNuVertexAlg][i]->Self();
      }
    }
  }


  // Save flash information for multiple algorithms
  if (fSaveFlashInfo) {

      auto const& channelMap = art::ServiceHandle<geo::WireReadout const>()->Get();

      for (unsigned int iFlashAlg=0; iFlashAlg < NFlashAlgos; ++iFlashAlg){
      AnalysisTreeDataStruct::FlashDataStruct& FlashData = fData->GetFlashData(iFlashAlg);

      size_t NFlashes = flashlist[iFlashAlg].size();
      
      FlashData.SetMaxFlashes(std::max(NFlashes, (size_t) 1));
      FlashData.Clear(); // clear all data
      FlashData.nfls = (short) NFlashes;
       
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetFlashAddresses(iFlashAlg);
      if (NFlashes > FlashData.GetMaxFlashes()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnalysisTree:limits") << "event has " << NFlashes
                                            << " " << fOpFlashModuleLabel[iFlashAlg] << " flashes, only "
                                            << FlashData.GetMaxFlashes() << " stored in tree";
      }
    
      for (size_t i = 0; i < NFlashes && i < kMaxFlashes ; ++i) { // loop over flashes
        FlashData.flsTime[i]       = flashlist[iFlashAlg][i]->Time();
        FlashData.flsPe[i]         = flashlist[iFlashAlg][i]->TotalPE();
        FlashData.flsYcenter[i]    = flashlist[iFlashAlg][i]->YCenter();
        FlashData.flsZcenter[i]    = flashlist[iFlashAlg][i]->ZCenter();
        FlashData.flsYwidth[i]     = flashlist[iFlashAlg][i]->YWidth();
        FlashData.flsZwidth[i]     = flashlist[iFlashAlg][i]->ZWidth();
        FlashData.flsTwidth[i]     = flashlist[iFlashAlg][i]->TimeWidth();

        for (int opch = 0; opch < kNOpDets; opch++) {
          Double_t pePerOpChannel = (Double_t) flashlist[iFlashAlg][i]->PE(opch);
          unsigned int opdet = channelMap.OpDetFromOpChannel(opch);
          FlashData.flsPePerOpDet[i][opdet] = pePerOpChannel;
        }

        // Now be dragons. simpleFlash always saves PEs per OpDet in a 32 position array
        // while opFlash in an n position one, depending on beam or cosmic opChannels
        // So in the case of opflashCosmic we need to look at cosmic discriminator channels
        // I know this is terrible, but I don't see a better way to do it now, but to go back
        // and have a better recob::OpFlash data product definition. --mdeltutt
        if (fOpFlashModuleLabel[iFlashAlg].find("opflashCosmic") != std::string::npos) {
          for (int opch = 200; opch < 200+kNOpDets; opch++) {
            Double_t pePerOpChannel = (Double_t) flashlist[iFlashAlg][i]->PE(opch);
            unsigned int opdet = channelMap.OpDetFromOpChannel(opch);
            FlashData.flsPePerOpDet[i][opdet] = pePerOpChannel;
          }
        }

        // Convert flash time to x coordinate (tjyang)
        geo::PlaneID pid(0, 0, 0);
        FlashData.flsXcenter[i] = detProp.ConvertTicksToX(flashlist[iFlashAlg][i]->Time()/sampling_rate(clockData)*1e3 + detProp.GetXTicksOffset(pid),pid);
      }
    }
  } // if fSaveFlashInfo

  //mc truth information
  if (isMC){
    if (fSaveCryInfo){ 
      //store cry (cosmic generator information) 
      fData->mcevts_truthcry = mclistcry.size();
      fData->cry_no_primaries = nCryPrimaries;
      //fData->cry_no_primaries;
      for(Int_t iPartc = 0; iPartc < mctruthcry->NParticles(); ++iPartc){
        const simb::MCParticle& partc(mctruthcry->GetParticle(iPartc));
        fData->cry_primaries_pdg[iPartc]=partc.PdgCode();
        fData->cry_Eng[iPartc]=partc.E();
        fData->cry_Px[iPartc]=partc.Px();
        fData->cry_Py[iPartc]=partc.Py();
        fData->cry_Pz[iPartc]=partc.Pz();
        fData->cry_P[iPartc]=partc.P();
	fData->cry_StartPointx[iPartc] = partc.Vx();
	fData->cry_StartPointy[iPartc] = partc.Vy();
	fData->cry_StartPointz[iPartc] = partc.Vz();
	fData->cry_StartPointt[iPartc] = partc.T();	
        fData->cry_status_code[iPartc]=partc.StatusCode();
        fData->cry_mass[iPartc]=partc.Mass();
        fData->cry_trackID[iPartc]=partc.TrackId();
        fData->cry_ND[iPartc]=partc.NumberDaughters();
        fData->cry_mother[iPartc]=partc.Mother();
      } // for cry particles  
    }// end fSaveCryInfo   
    //save neutrino interaction information
    fData->mcevts_truth = mclist.size();
    // Get space charge correction
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    // Get time offset for x space charge correction
    if (fData->mcevts_truth > 0){//at least one mc record
      if (fSaveGenieInfo){
	int neutrino_i = 0;
	for(unsigned int iList = 0; (iList < mclist.size()) && (neutrino_i < kMaxTruth) ; ++iList){
	  if (mclist[iList]->NeutrinoSet()){
	    // Get time offset for space charge correction
            double g4Ticks = clockData.TPCG4Time2Tick(mclist[iList]->GetNeutrino().Nu().T()) + detProp.GetXTicksOffset(0,0,0) - trigger_offset(clockData);
            double xtimeoffset = detProp.ConvertTicksToX(g4Ticks,0,0,0);
	    
	    fData->nuPDG_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().PdgCode();
	    fData->ccnc_truth[neutrino_i]   = mclist[iList]->GetNeutrino().CCNC();
	    fData->mode_truth[neutrino_i]   = mclist[iList]->GetNeutrino().Mode();
	    fData->Q2_truth[neutrino_i]     = mclist[iList]->GetNeutrino().QSqr();
	    fData->W_truth[neutrino_i]      = mclist[iList]->GetNeutrino().W();
	    fData->X_truth[neutrino_i]      = mclist[iList]->GetNeutrino().X();
	    fData->Y_truth[neutrino_i]      = mclist[iList]->GetNeutrino().Y();
	    fData->hitnuc_truth[neutrino_i] = mclist[iList]->GetNeutrino().HitNuc();
	    fData->enu_truth[neutrino_i]    = mclist[iList]->GetNeutrino().Nu().E();
	    fData->nuvtxx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vx();
	    fData->nuvtxy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vy();
	    fData->nuvtxz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vz();
	    fData->sp_charge_corrected_nuvtxx_truth[neutrino_i] = (mclist[iList]->GetNeutrino().Nu().Vx() - SCE->GetPosOffsets(geo::Point_t(mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz())).X() + xtimeoffset)*(1.114/1.098) + 0.6;
	    fData->sp_charge_corrected_nuvtxy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vy() + SCE->GetPosOffsets(geo::Point_t(mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz())).Y();
	    fData->sp_charge_corrected_nuvtxz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vz() + SCE->GetPosOffsets(geo::Point_t(mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz())).Z();
	    
	    if (mclist[iList]->GetNeutrino().Nu().P()){
	      fData->nu_dcosx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Px()/mclist[iList]->GetNeutrino().Nu().P();
	      fData->nu_dcosy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Py()/mclist[iList]->GetNeutrino().Nu().P();
	      fData->nu_dcosz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Pz()/mclist[iList]->GetNeutrino().Nu().P();
	    }
	    fData->lep_mom_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().P();
	    if (mclist[iList]->GetNeutrino().Lepton().P()){
	      fData->lep_dcosx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Px()/mclist[iList]->GetNeutrino().Lepton().P();
	      fData->lep_dcosy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Py()/mclist[iList]->GetNeutrino().Lepton().P();
	      fData->lep_dcosz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Pz()/mclist[iList]->GetNeutrino().Lepton().P();
	    }

	    //flux information
	    fData->vx_flux[neutrino_i]        = fluxlist[iList]->fvx;
	    fData->vy_flux[neutrino_i]        = fluxlist[iList]->fvy;
	    fData->vz_flux[neutrino_i]        = fluxlist[iList]->fvz;
	    fData->pdpx_flux[neutrino_i]      = fluxlist[iList]->fpdpx;
	    fData->pdpy_flux[neutrino_i]      = fluxlist[iList]->fpdpy;
	    fData->pdpz_flux[neutrino_i]      = fluxlist[iList]->fpdpz;
	    fData->ppdxdz_flux[neutrino_i]    = fluxlist[iList]->fppdxdz;
	    fData->ppdydz_flux[neutrino_i]    = fluxlist[iList]->fppdydz;
	    fData->pppz_flux[neutrino_i]      = fluxlist[iList]->fpppz;
	  
	    fData->ptype_flux[neutrino_i]      = fluxlist[iList]->fptype;
	    fData->ppvx_flux[neutrino_i]       = fluxlist[iList]->fppvx;
	    fData->ppvy_flux[neutrino_i]       = fluxlist[iList]->fppvy;
	    fData->ppvz_flux[neutrino_i]       = fluxlist[iList]->fppvz;
	    fData->muparpx_flux[neutrino_i]    = fluxlist[iList]->fmuparpx;
	    fData->muparpy_flux[neutrino_i]    = fluxlist[iList]->fmuparpy;
	    fData->muparpz_flux[neutrino_i]    = fluxlist[iList]->fmuparpz;
	    fData->mupare_flux[neutrino_i]     = fluxlist[iList]->fmupare;
	  
	    fData->tgen_flux[neutrino_i]     = fluxlist[iList]->ftgen;
	    fData->tgptype_flux[neutrino_i]  = fluxlist[iList]->ftgptype;
	    fData->tgppx_flux[neutrino_i]    = fluxlist[iList]->ftgppx;
	    fData->tgppy_flux[neutrino_i]    = fluxlist[iList]->ftgppy;
	    fData->tgppz_flux[neutrino_i]    = fluxlist[iList]->ftgppz;
	    fData->tprivx_flux[neutrino_i]   = fluxlist[iList]->ftprivx;
	    fData->tprivy_flux[neutrino_i]   = fluxlist[iList]->ftprivy;
	    fData->tprivz_flux[neutrino_i]   = fluxlist[iList]->ftprivz;

	    fData->dk2gen_flux[neutrino_i]   = fluxlist[iList]->fdk2gen;
	    fData->gen2vtx_flux[neutrino_i]   = fluxlist[iList]->fgen2vtx;

	    fData->tpx_flux[neutrino_i]    = fluxlist[iList]->ftpx;
	    fData->tpy_flux[neutrino_i]    = fluxlist[iList]->ftpy;
	    fData->tpz_flux[neutrino_i]    = fluxlist[iList]->ftpz;
	    fData->tptype_flux[neutrino_i] = fluxlist[iList]->ftptype;
	    neutrino_i++;
	  }//mclist is NeutrinoSet()
	}//loop over mclist

	if (mctruth->NeutrinoSet()){
	  //genie particles information
	  fData->genie_no_primaries = mctruth->NParticles();

	  size_t StoreParticles = std::min((size_t) fData->genie_no_primaries, fData->GetMaxGeniePrimaries());
	  if (fData->genie_no_primaries > (int) StoreParticles) {
	    // got this error? it might be a bug,
	    // since the structure should have enough room for everything
	    mf::LogError("AnalysisTree:limits") << "event has "
						<< fData->genie_no_primaries << " MC particles, only "
						<< StoreParticles << " stored in tree";
	  }
	  for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
	    const simb::MCParticle& part(mctruth->GetParticle(iPart));
	    fData->genie_primaries_pdg[iPart]=part.PdgCode();
	    fData->genie_Eng[iPart]=part.E();
	    fData->genie_Px[iPart]=part.Px();
	    fData->genie_Py[iPart]=part.Py();
	    fData->genie_Pz[iPart]=part.Pz();
	    fData->genie_P[iPart]=part.P();
	    fData->genie_status_code[iPart]=part.StatusCode();
	    fData->genie_mass[iPart]=part.Mass();
	    fData->genie_trackID[iPart]=part.TrackId();
	    fData->genie_ND[iPart]=part.NumberDaughters();
	    fData->genie_mother[iPart]=part.Mother();
	  } // for particle
	  //const simb::MCNeutrino& nu(mctruth->GetNeutrino());
	} //if neutrino set
      }// end (fSaveGenieInfo)  

      //Extract MC Shower information and fill the Shower branches
      if (fSaveMCShowerInfo && mcshowerh.isValid()){
	fData->no_mcshowers = nMCShowers;       
	size_t shwr = 0;
	for(std::vector<sim::MCShower>::const_iterator imcshwr = mcshowerh->begin();    
	    imcshwr != mcshowerh->end(); ++imcshwr) {
    	  const sim::MCShower& mcshwr = *imcshwr;
	  fData->mcshwr_origin[shwr]          = mcshwr.Origin();
    	  fData->mcshwr_pdg[shwr]	      = mcshwr.PdgCode();	   
    	  fData->mcshwr_TrackId[shwr]	      = mcshwr.TrackID();	   
    	  fData->mcshwr_Process[shwr]	      = mcshwr.Process();	   
    	  fData->mcshwr_startX[shwr]          = mcshwr.Start().X();	    
    	  fData->mcshwr_startY[shwr]          = mcshwr.Start().Y();      	 
    	  fData->mcshwr_startZ[shwr]          = mcshwr.Start().Z();      	 
    	  fData->mcshwr_endX[shwr]            = mcshwr.End().X();      	 
    	  fData->mcshwr_endY[shwr]            = mcshwr.End().Y();      	 
    	  fData->mcshwr_endZ[shwr]            = mcshwr.End().Z();
	  if (mcshwr.DetProfile().E()!= 0){   
	    fData->mcshwr_isEngDeposited[shwr] = 1;	
	    fData->mcshwr_CombEngX[shwr]        = mcshwr.DetProfile().X();      	 
	    fData->mcshwr_CombEngY[shwr]        = mcshwr.DetProfile().Y();       	 
	    fData->mcshwr_CombEngZ[shwr]        = mcshwr.DetProfile().Z();      	 
	    fData->mcshwr_CombEngPx[shwr]       = mcshwr.DetProfile().Px();     	 
	    fData->mcshwr_CombEngPy[shwr]       = mcshwr.DetProfile().Py();       	 
	    fData->mcshwr_CombEngPz[shwr]       = mcshwr.DetProfile().Pz();       	 
	    fData->mcshwr_CombEngE[shwr]        = mcshwr.DetProfile().E();
	    fData->mcshwr_dEdx[shwr]            = mcshwr.dEdx();
	    fData->mcshwr_StartDirX[shwr]       = mcshwr.StartDir().X();
	    fData->mcshwr_StartDirY[shwr]       = mcshwr.StartDir().Y();
	    fData->mcshwr_StartDirZ[shwr]       = mcshwr.StartDir().Z();
	  }
	  else
	    fData->mcshwr_isEngDeposited[shwr] = 0;
    	  fData->mcshwr_Motherpdg[shwr]       = mcshwr.MotherPdgCode(); 
    	  fData->mcshwr_MotherTrkId[shwr]     = mcshwr.MotherTrackID();    
    	  fData->mcshwr_MotherProcess[shwr]   = mcshwr.MotherProcess(); 
	  fData->mcshwr_MotherstartX[shwr]    = mcshwr.MotherStart().X();     
    	  fData->mcshwr_MotherstartY[shwr]    = mcshwr.MotherStart().Y();	  
    	  fData->mcshwr_MotherstartZ[shwr]    = mcshwr.MotherStart().Z();	  
    	  fData->mcshwr_MotherendX[shwr]      = mcshwr.MotherEnd().X(); 	  
    	  fData->mcshwr_MotherendY[shwr]      = mcshwr.MotherEnd().Y(); 	  
    	  fData->mcshwr_MotherendZ[shwr]      = mcshwr.MotherEnd().Z();  
    	  fData->mcshwr_Ancestorpdg[shwr]     = mcshwr.AncestorPdgCode();  
    	  fData->mcshwr_AncestorTrkId[shwr]   = mcshwr.AncestorTrackID();  
    	  fData->mcshwr_AncestorProcess[shwr] = mcshwr.AncestorProcess();
	  fData->mcshwr_AncestorstartX[shwr]  = mcshwr.AncestorStart().X();	  
    	  fData->mcshwr_AncestorstartY[shwr]  = mcshwr.AncestorStart().Y();		 
    	  fData->mcshwr_AncestorstartZ[shwr]  = mcshwr.AncestorStart().Z();		 
    	  fData->mcshwr_AncestorendX[shwr]    = mcshwr.AncestorEnd().X();	 
    	  fData->mcshwr_AncestorendY[shwr]    = mcshwr.AncestorEnd().Y();	 
    	  fData->mcshwr_AncestorendZ[shwr]    = mcshwr.AncestorEnd().Z();
    	  ++shwr; 
	}
	fData->mcshwr_Process.resize(shwr);
	fData->mcshwr_MotherProcess.resize(shwr);
	fData->mcshwr_AncestorProcess.resize(shwr);
      }//End if (fSaveMCShowerInfo){  
    
      //Extract MC Track information and fill the Shower branches
      if (fSaveMCTrackInfo && mctrackh.isValid()){
	fData->no_mctracks = nMCTracks;       
	size_t trk = 0;
	for(std::vector<sim::MCTrack>::const_iterator imctrk = mctrackh->begin();imctrk != mctrackh->end(); ++imctrk) {
    	  const sim::MCTrack& mctrk = *imctrk;
	  TLorentzVector tpcstart, tpcend, tpcmom;
          double plen = driftedLength(detProp, mctrk, tpcstart, tpcend, tpcmom);
	  fData->mctrk_origin[trk]          = mctrk.Origin();
    	  fData->mctrk_pdg[trk]	            = mctrk.PdgCode();	   
    	  fData->mctrk_TrackId[trk]	    = mctrk.TrackID();	   
    	  fData->mctrk_Process[trk]	    = mctrk.Process();	   
    	  fData->mctrk_startX[trk]          = mctrk.Start().X();	    
    	  fData->mctrk_startY[trk]          = mctrk.Start().Y();      	 
    	  fData->mctrk_startZ[trk]          = mctrk.Start().Z();      	 
    	  fData->mctrk_endX[trk]            = mctrk.End().X();      	 
    	  fData->mctrk_endY[trk]            = mctrk.End().Y();      	 
    	  fData->mctrk_endZ[trk]            = mctrk.End().Z(); 
	  fData->mctrk_Motherpdg[trk]       = mctrk.MotherPdgCode(); 
    	  fData->mctrk_MotherTrkId[trk]     = mctrk.MotherTrackID();    
    	  fData->mctrk_MotherProcess[trk]   = mctrk.MotherProcess(); 
	  fData->mctrk_MotherstartX[trk]    = mctrk.MotherStart().X();     
    	  fData->mctrk_MotherstartY[trk]    = mctrk.MotherStart().Y();	  
    	  fData->mctrk_MotherstartZ[trk]    = mctrk.MotherStart().Z();	  
    	  fData->mctrk_MotherendX[trk]      = mctrk.MotherEnd().X(); 	  
    	  fData->mctrk_MotherendY[trk]      = mctrk.MotherEnd().Y(); 	  
    	  fData->mctrk_MotherendZ[trk]      = mctrk.MotherEnd().Z();  
    	  fData->mctrk_Ancestorpdg[trk]     = mctrk.AncestorPdgCode();  
    	  fData->mctrk_AncestorTrkId[trk]   = mctrk.AncestorTrackID();  
    	  fData->mctrk_AncestorProcess[trk] = mctrk.AncestorProcess();
	  fData->mctrk_AncestorstartX[trk]  = mctrk.AncestorStart().X();	  
    	  fData->mctrk_AncestorstartY[trk]  = mctrk.AncestorStart().Y();		 
    	  fData->mctrk_AncestorstartZ[trk]  = mctrk.AncestorStart().Z();		 
    	  fData->mctrk_AncestorendX[trk]    = mctrk.AncestorEnd().X();	 
    	  fData->mctrk_AncestorendY[trk]    = mctrk.AncestorEnd().Y();	 
    	  fData->mctrk_AncestorendZ[trk]    = mctrk.AncestorEnd().Z();
        
	  fData->mctrk_len_drifted[trk]       = plen;

	  if (plen != 0){	  
	    fData->mctrk_startX_drifted[trk] = tpcstart.X();
	    fData->mctrk_startY_drifted[trk] = tpcstart.Y();
	    fData->mctrk_startZ_drifted[trk] = tpcstart.Z();
	    fData->mctrk_endX_drifted[trk]   = tpcend.X();
	    fData->mctrk_endY_drifted[trk]   = tpcend.Y();
	    fData->mctrk_endZ_drifted[trk]   = tpcend.Z();
	    fData->mctrk_p_drifted[trk]      = tpcmom.Vect().Mag();
	    fData->mctrk_px_drifted[trk]     = tpcmom.X();
	    fData->mctrk_py_drifted[trk]     = tpcmom.Y();
	    fData->mctrk_pz_drifted[trk]     = tpcmom.Z();
	  }
    	  ++trk; 
	}
       
	fData->mctrk_Process.resize(trk);
	fData->mctrk_MotherProcess.resize(trk);
	fData->mctrk_AncestorProcess.resize(trk);
      }//End if (fSaveMCTrackInfo){  
	  

      //GEANT particles information
      if (fSaveGeantInfo){ 

        const sim::ParticleList& plist = fMCTruthMatching->ParticleList();
        
        std::string pri("primary");
        int primary=0;
        int active = 0;
        size_t geant_particle=0;
        sim::ParticleList::const_iterator itPart = plist.begin(),
          pend = plist.end(); // iterator to pairs (track id, particle)
        	  
	// helper map track ID => index
	std::map<int, size_t> TrackIDtoIndex;
	std::vector<int> gpdg;
	std::vector<int> gmother;
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 
	    
        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart) {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }
          
          //++geant_particle;
	  bool isPrimary = pPart->Process() == pri;
	  int TrackID = pPart->TrackId();
	  TrackIDtoIndex.emplace(TrackID, iPart);
	  gpdg.push_back(pPart->PdgCode());
	  gmother.push_back(pPart->Mother());
          if (iPart < fData->GetMaxGEANTparticles()) {
	    if (pPart->E()<fG4minE&&(!isPrimary)) continue;
	    if (isPrimary) ++primary;
	    
	    art::Ptr<simb::MCTruth> const& mc_truth = fMCTruthMatching->ParticleToMCTruth(pPart);
	    // Get time offset for space charge correction
	    double xtimeoffset = 0;
	    if (mc_truth){
              double g4Ticks = detProp.GetXTicksOffset(0,0,0) - trigger_offset(clockData);
	      if (mc_truth->NeutrinoSet()){
                g4Ticks += clockData.TPCG4Time2Tick(mc_truth->GetNeutrino().Nu().T());
	      }
              xtimeoffset = detProp.ConvertTicksToX(g4Ticks,0,0,0);
	    }
	    
	    TLorentzVector mcstart, mcend, mcstartdrifted, mcenddrifted;
	    unsigned int pstarti, pendi, pstartdriftedi, penddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	    double plen = length(*pPart, mcstart, mcend, pstarti, pendi);
            double plendrifted = driftedLength(detProp, *pPart, mcstartdrifted, mcenddrifted, pstartdriftedi, penddriftedi);
	    
	    bool isActive = plen != 0;
	    bool isDrifted = plendrifted!= 0;
	    if (plen) ++active;

	    fData->process_primary[geant_particle] = int(isPrimary);
	    fData->processname[geant_particle]= pPart->Process();
	    fData->Mother[geant_particle]=pPart->Mother();
	    fData->TrackId[geant_particle]=TrackID;
	    fData->pdg[geant_particle]=pPart->PdgCode();
	    fData->status[geant_particle] = pPart->StatusCode();
	    fData->Eng[geant_particle]=pPart->E();
	    fData->EndE[geant_particle]=pPart->EndE();
	    fData->Mass[geant_particle]=pPart->Mass();
	    fData->Px[geant_particle]=pPart->Px();
	    fData->Py[geant_particle]=pPart->Py();
	    fData->Pz[geant_particle]=pPart->Pz();
	    fData->P[geant_particle]=pPart->Momentum().Vect().Mag();
	    fData->StartPointx[geant_particle]=pPart->Vx();
	    fData->StartPointy[geant_particle]=pPart->Vy();
	    fData->StartPointz[geant_particle]=pPart->Vz();
	    fData->sp_charge_corrected_StartPointx[geant_particle] = (pPart->Vx() - SCE->GetPosOffsets(geo::Point_t(pPart->Vx(),pPart->Vy(),pPart->Vz())).X() + xtimeoffset)*(1.114/1.098) + 0.6; 
	    fData->sp_charge_corrected_StartPointy[geant_particle] = pPart->Vy() + SCE->GetPosOffsets(geo::Point_t(pPart->Vx(),pPart->Vy(),pPart->Vz())).Y(); 
	    fData->sp_charge_corrected_StartPointz[geant_particle] = pPart->Vz() + SCE->GetPosOffsets(geo::Point_t(pPart->Vx(),pPart->Vy(),pPart->Vz())).Z(); 
	    fData->StartT[geant_particle] = pPart->T();
	    fData->EndPointx[geant_particle]=pPart->EndPosition()[0];
	    fData->EndPointy[geant_particle]=pPart->EndPosition()[1];
	    fData->EndPointz[geant_particle]=pPart->EndPosition()[2];
	    fData->sp_charge_corrected_EndPointx[geant_particle] = (pPart->EndPosition()[0] - SCE->GetPosOffsets(geo::Point_t(pPart->EndPosition()[0],pPart->EndPosition()[1],pPart->EndPosition()[2])).X() + xtimeoffset)*(1.114/1.098) + 0.6; 
	    fData->sp_charge_corrected_EndPointy[geant_particle] = pPart->EndPosition()[1] + SCE->GetPosOffsets(geo::Point_t(pPart->EndPosition()[0],pPart->EndPosition()[1],pPart->EndPosition()[2])).Y(); 
	    fData->sp_charge_corrected_EndPointz[geant_particle] = pPart->EndPosition()[2] + SCE->GetPosOffsets(geo::Point_t(pPart->EndPosition()[0],pPart->EndPosition()[1],pPart->EndPosition()[2])).Z(); 
	    fData->EndT[geant_particle] = pPart->EndT();
	    fData->theta[geant_particle] = pPart->Momentum().Theta();
	    fData->phi[geant_particle] = pPart->Momentum().Phi();
	    fData->theta_xz[geant_particle] = std::atan2(pPart->Px(), pPart->Pz());
	    fData->theta_yz[geant_particle] = std::atan2(pPart->Py(), pPart->Pz());
	    fData->pathlen[geant_particle]  = plen;
	    fData->pathlen_drifted[geant_particle]  = plendrifted;
	    fData->NumberDaughters[geant_particle]=pPart->NumberDaughters();
	    fData->inTPCActive[geant_particle] = int(isActive);
	    fData->inTPCDrifted[geant_particle] = int(isDrifted);
	    if (mc_truth){
	      fData->origin[geant_particle] = mc_truth->Origin();
	      fData->MCTruthIndex[geant_particle] = mc_truth.key();
	    }
	    if (isActive){	  
	      fData->StartPointx_tpcAV[geant_particle] = mcstart.X();
	      fData->StartPointy_tpcAV[geant_particle] = mcstart.Y();
	      fData->StartPointz_tpcAV[geant_particle] = mcstart.Z();
	      fData->sp_charge_corrected_StartPointx_tpcAV[geant_particle] = (mcstart.X() - SCE->GetPosOffsets(geo::Point_t(mcstart.X(),mcstart.Y(),mcstart.Z())).X() + xtimeoffset)*(1.114/1.098) + 0.6; 
	      fData->sp_charge_corrected_StartPointy_tpcAV[geant_particle] = mcstart.Y() + SCE->GetPosOffsets(geo::Point_t(mcstart.X(),mcstart.Y(),mcstart.Z())).Y(); 
	      fData->sp_charge_corrected_StartPointz_tpcAV[geant_particle] = mcstart.Z() + SCE->GetPosOffsets(geo::Point_t(mcstart.X(),mcstart.Y(),mcstart.Z())).Z(); 
	      fData->StartT_tpcAV[geant_particle] = mcstart.T();
	      fData->StartE_tpcAV[geant_particle] = pPart->E(pstarti);
	      fData->StartP_tpcAV[geant_particle] = pPart->P(pstarti);
	      fData->StartPx_tpcAV[geant_particle] = pPart->Px(pstarti);
	      fData->StartPy_tpcAV[geant_particle] = pPart->Py(pstarti);
	      fData->StartPz_tpcAV[geant_particle] = pPart->Pz(pstarti);
	      fData->EndPointx_tpcAV[geant_particle] = mcend.X();
	      fData->EndPointy_tpcAV[geant_particle] = mcend.Y();
	      fData->EndPointz_tpcAV[geant_particle] = mcend.Z();
	      fData->sp_charge_corrected_EndPointx_tpcAV[geant_particle] = (mcend.X() - SCE->GetPosOffsets(geo::Point_t(mcend.X(),mcend.Y(),mcend.Z())).X() + xtimeoffset)*(1.114/1.098) + 0.6; 
	      fData->sp_charge_corrected_EndPointy_tpcAV[geant_particle] = mcend.Y() + SCE->GetPosOffsets(geo::Point_t(mcend.X(),mcend.Y(),mcend.Z())).Y(); 
	      fData->sp_charge_corrected_EndPointz_tpcAV[geant_particle] = mcend.Z() + SCE->GetPosOffsets(geo::Point_t(mcend.X(),mcend.Y(),mcend.Z())).Z(); 
	      fData->EndT_tpcAV[geant_particle] = mcend.T();
	      fData->EndE_tpcAV[geant_particle] = pPart->E(pendi);
	      fData->EndP_tpcAV[geant_particle] = pPart->P(pendi);
	      fData->EndPx_tpcAV[geant_particle] = pPart->Px(pendi);
	      fData->EndPy_tpcAV[geant_particle] = pPart->Py(pendi);
	      fData->EndPz_tpcAV[geant_particle] = pPart->Pz(pendi);
	    }
	    if (isDrifted){	  
	      fData->StartPointx_drifted[geant_particle] = mcstartdrifted.X();
	      fData->StartPointy_drifted[geant_particle] = mcstartdrifted.Y();
	      fData->StartPointz_drifted[geant_particle] = mcstartdrifted.Z();
	      fData->sp_charge_corrected_StartPointx_drifted[geant_particle] = (mcstartdrifted.X() - SCE->GetPosOffsets(geo::Point_t(mcstartdrifted.X(),mcstartdrifted.Y(),mcstartdrifted.Z())).X() + xtimeoffset)*(1.114/1.098) + 0.6;
	      fData->sp_charge_corrected_StartPointy_drifted[geant_particle] = mcstartdrifted.Y() + SCE->GetPosOffsets(geo::Point_t(mcstartdrifted.X(),mcstartdrifted.Y(),mcstartdrifted.Z())).Y();
	      fData->sp_charge_corrected_StartPointz_drifted[geant_particle] = mcstartdrifted.Z() + SCE->GetPosOffsets(geo::Point_t(mcstartdrifted.X(),mcstartdrifted.Y(),mcstartdrifted.Z())).Z();
	      fData->StartT_drifted[geant_particle] = mcstartdrifted.T();
	      fData->StartE_drifted[geant_particle] = pPart->E(pstartdriftedi);
	      fData->StartP_drifted[geant_particle] = pPart->P(pstartdriftedi);
	      fData->StartPx_drifted[geant_particle] = pPart->Px(pstartdriftedi);
	      fData->StartPy_drifted[geant_particle] = pPart->Py(pstartdriftedi);
	      fData->StartPz_drifted[geant_particle] = pPart->Pz(pstartdriftedi);
	      fData->EndPointx_drifted[geant_particle] = mcenddrifted.X();
	      fData->EndPointy_drifted[geant_particle] = mcenddrifted.Y();
	      fData->EndPointz_drifted[geant_particle] = mcenddrifted.Z();
	      fData->sp_charge_corrected_EndPointx_drifted[geant_particle] = (mcenddrifted.X() - SCE->GetPosOffsets(geo::Point_t(mcenddrifted.X(),mcenddrifted.Y(),mcenddrifted.Z())).X() + xtimeoffset)*(1.114/1.098) + 0.6;
	      fData->sp_charge_corrected_EndPointy_drifted[geant_particle] = mcenddrifted.Y() + SCE->GetPosOffsets(geo::Point_t(mcenddrifted.X(),mcenddrifted.Y(),mcenddrifted.Z())).Y();
	      fData->sp_charge_corrected_EndPointz_drifted[geant_particle] = mcenddrifted.Z() + SCE->GetPosOffsets(geo::Point_t(mcenddrifted.X(),mcenddrifted.Y(),mcenddrifted.Z())).Z();
	      fData->EndT_drifted[geant_particle] = mcenddrifted.T();
	      fData->EndE_drifted[geant_particle] = pPart->E(penddriftedi);
	      fData->EndP_drifted[geant_particle] = pPart->P(penddriftedi);
	      fData->EndPx_drifted[geant_particle] = pPart->Px(penddriftedi);
	      fData->EndPy_drifted[geant_particle] = pPart->Py(penddriftedi);
	      fData->EndPz_drifted[geant_particle] = pPart->Pz(penddriftedi);
	    }	            
	    //access auxiliary detector parameters
            if (fSaveAuxDetInfo) {
              unsigned short nAD = 0; // number of cells that particle hit
              
              // find deposit of this particle in each of the detector cells
              for (const sim::AuxDetSimChannel* c: fAuxDetSimChannels) {
		
        	// find if this cell has a contribution (IDE) from this particle,
        	// and which one
        	const std::vector<sim::AuxDetIDE>& setOfIDEs = c->AuxDetIDEs();
        	// using a C++ "lambda" function here; this one:
        	// - sees only TrackID from the current scope
        	// - takes one parameter: the AuxDetIDE to be tested
        	// - returns if that IDE belongs to the track we are looking for
        	std::vector<sim::AuxDetIDE>::const_iterator iIDE
        	  = std::find_if(
				 setOfIDEs.begin(), setOfIDEs.end(),
				 [TrackID](const sim::AuxDetIDE& IDE){ return IDE.trackID == TrackID; }
				 );
        	if (iIDE == setOfIDEs.end()) continue;
        	
        	// now iIDE points to the energy released by the track #i (TrackID)
		
		// look for IDE with matching trackID
		// find trackIDs stored in setOfIDEs with the same trackID, but negative,
		// this is an untracked particle who's energy should be added as deposited by this original trackID
		float totalE = 0.; // total energy deposited around by the GEANT particle in this cell
		for(const auto& adtracks: setOfIDEs) {
		  if( fabs(adtracks.trackID) == TrackID )
		    totalE += adtracks.energyDeposited;
		} // for
		
		// fill the structure
		if (nAD < kMaxAuxDets) {
		  fData->AuxDetID[geant_particle][nAD] = c->AuxDetID();
		  fData->entryX[geant_particle][nAD]   = iIDE->entryX;
		  fData->entryY[geant_particle][nAD]   = iIDE->entryY;
		  fData->entryZ[geant_particle][nAD]   = iIDE->entryZ;
		  fData->entryT[geant_particle][nAD]   = iIDE->entryT;
		  fData->exitX[geant_particle][nAD]    = iIDE->exitX;
		  fData->exitY[geant_particle][nAD]    = iIDE->exitY;
		  fData->exitZ[geant_particle][nAD]    = iIDE->exitZ;
		  fData->exitT[geant_particle][nAD]    = iIDE->exitT;
		  fData->exitPx[geant_particle][nAD]   = iIDE->exitMomentumX;
		  fData->exitPy[geant_particle][nAD]   = iIDE->exitMomentumY;
		  fData->exitPz[geant_particle][nAD]   = iIDE->exitMomentumZ;
		  fData->CombinedEnergyDep[geant_particle][nAD] = totalE;
		}
		++nAD;
	      } // for aux det sim channels
	      fData->NAuxDets[geant_particle] = nAD; 
	      
	      if (nAD > kMaxAuxDets) {
		// got this error? consider increasing kMaxAuxDets
		mf::LogError("AnalysisTree:limits") 
		  << "particle #" << iPart
		  << " touches " << nAD << " auxiliary detector cells, only "
		  << kMaxAuxDets << " of them are saved in the tree";
	      } // if too many detector cells
	    } // if (fSaveAuxDetInfo) 
	    	    
	    ++geant_particle;
	  }
	  else if (iPart == fData->GetMaxGEANTparticles()) {
	    // got this error? it might be a bug,
	    // since the structure should have enough room for everything
	    mf::LogError("AnalysisTree:limits") << "event has "
						<< plist.size() << " MC particles, only "
						<< fData->GetMaxGEANTparticles() << " will be stored in tree";
	  }     
	} // for particles
	
	fData->geant_list_size_in_tpcAV = active;
	fData->no_primaries = primary;
	fData->geant_list_size = geant_particle;
	fData->processname.resize(geant_particle);
	MF_LOG_DEBUG("AnalysisTree") 
	  << "Counted "
	  << fData->geant_list_size << " GEANT particles ("
	  << fData->geant_list_size_in_tpcAV << " in AV), "
	  << fData->no_primaries << " primaries, "
	  << fData->genie_no_primaries << " GENIE particles";

	FillWith(fData->MergedId, 0);

	// for each particle, consider all the direct ancestors with the same
	// PDG ID, and mark them as belonging to the same "group"
	// (having the same MergedId)
	/* turn off for now
	   int currentMergedId = 1;
	   for(size_t iPart = 0; iPart < geant_particle; ++iPart){
	   // if the particle already belongs to a group, don't bother
	   if (fData->MergedId[iPart]) continue;
	   // the particle starts its own group
	   fData->MergedId[iPart] = currentMergedId;
	   int currentMotherTrackId = fData->Mother[iPart];
	   while (currentMotherTrackId > 0) {
	   if (TrackIDtoIndex.find(currentMotherTrackId)==TrackIDtoIndex.end()) break;
	   size_t gindex = TrackIDtoIndex[currentMotherTrackId];
	   if (gindex<0||gindex>=plist.size()) break;
	   // if the mother particle is of a different type,
	   // don't bother with iPart ancestry any further
	   if (gpdg[gindex]!=fData->pdg[iPart]) break;
	   if (TrackIDtoIndex.find(currentMotherTrackId)!=TrackIDtoIndex.end()){
	   size_t igeantMother = TrackIDtoIndex[currentMotherTrackId];
	   if (igeantMother>=0&&igeantMother<geant_particle){
	   fData->MergedId[igeantMother] = currentMergedId;
	   }
	   }
	   currentMotherTrackId = gmother[gindex];
	   }
	   ++currentMergedId;
	   }// for merging check
	*/
      } // if (fSaveGeantInfo) 
    }//if (mcevts_truth)
  }//if (isMC){
  fData->taulife = detProp.ElectronLifetime();
  fTree->Fill();
  
  if (mf::isDebugEnabled()) {
    // use mf::LogDebug instead of MF_LOG_DEBUG because we reuse it in many lines;
    // thus, we protect this part of the code with the line above
    mf::LogDebug logStream("AnalysisTreeStructure");
    logStream
      << "Tree data structure contains:"
      << "\n - " << fData->no_hits << " hits (" << fData->GetMaxHits() << ")"
      << "\n - " << fData->genie_no_primaries << " genie primaries (" << fData->GetMaxGeniePrimaries() << ")"
      << "\n - " << fData->geant_list_size << " GEANT particles (" << fData->GetMaxGEANTparticles() << "), "
      << fData->no_primaries << " primaries"
      << "\n - " << fData->geant_list_size_in_tpcAV << " GEANT particles in AV "
      << "\n - " << ((int) fData->kNTracker) << " trackers:"
      ;
    
    size_t iTracker = 0;
    for (auto tracker = fData->TrackData.cbegin();
	 tracker != fData->TrackData.cend(); ++tracker, ++iTracker
	 ) {
      logStream
	<< "\n -> " << tracker->ntracks << " " << fTrackModuleLabel[iTracker]
	<< " tracks (" << tracker->GetMaxTracks() << ")"
	;
      for (int iTrk = 0; iTrk < tracker->ntracks; ++iTrk) {
        logStream << "\n    [" << iTrk << "] "<< tracker->ntrkhits[iTrk][0];
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->ntrkhits[iTrk][ipl];
        logStream << " hits (" << tracker->GetMaxHitsPerTrack(iTrk, 0);
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->GetMaxHitsPerTrack(iTrk, ipl);
        logStream << ")";
      } // for tracks
    } // for trackers
  } // if logging enabled

  
  // if we don't use a permanent buffer (which can be huge),
  // delete the current buffer, and we'll create a new one on the next event
  if (!fUseBuffer) {
    MF_LOG_DEBUG("AnalysisTreeStructure") << "Freeing the tree data structure";
    DestroyData();
  }
} // microboone::AnalysisTree::analyze()


void microboone::AnalysisTree::FillShower(
					  AnalysisTreeDataStruct::ShowerDataStruct& showerData, size_t iShower,
					  recob::Shower const& shower, const bool fSavePFParticleInfo, 
                      const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
					  ) const {
  
  showerData.showerID[iShower]        = shower.ID();
  showerData.shwr_bestplane[iShower]  = shower.best_plane();
  showerData.shwr_length[iShower]     = shower.Length();
  
  TVector3 const& dir_start = shower.Direction();
  double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
  double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
  showerData.shwr_theta[iShower] = dir_start.Theta();
  showerData.shwr_phi[iShower] = dir_start.Phi();
  showerData.shwr_thetaxz[iShower] = theta_xz;
  showerData.shwr_thetayz[iShower] = theta_yz;
  showerData.shwr_startdcosx[iShower] = dir_start.X();
  showerData.shwr_startdcosy[iShower] = dir_start.Y();
  showerData.shwr_startdcosz[iShower] = dir_start.Z();

  TVector3 const& pos_start = shower.ShowerStart();
  showerData.shwr_startx[iShower]     = pos_start.X();
  showerData.shwr_starty[iShower]     = pos_start.Y();
  showerData.shwr_startz[iShower]     = pos_start.Z();
  
  if (fSavePFParticleInfo) {
    auto mapIter = showerIDtoPFParticleIDMap.find(shower.ID());
    if (mapIter != showerIDtoPFParticleIDMap.end()) {
      // This vertex has a corresponding PFParticle.
      showerData.shwr_hasPFParticle[iShower] = 1;
      showerData.shwr_PFParticleID[iShower] = mapIter->second;
    }
    else 
      showerData.shwr_hasPFParticle[iShower] = 0;
  }
  
  if (shower.Energy().size() == kNplanes)
    std::copy_n
      (shower.Energy().begin(),    kNplanes, &showerData.shwr_totEng[iShower][0]);
  if (shower.dEdx().size() == kNplanes)
    std::copy_n
      (shower.dEdx().begin(),      kNplanes, &showerData.shwr_dedx[iShower][0]);
  if (shower.MIPEnergy().size() == kNplanes)
    std::copy_n
      (shower.MIPEnergy().begin(), kNplanes, &showerData.shwr_mipEng[iShower][0]);
  
} // microboone::AnalysisTree::FillShower()


void microboone::AnalysisTree::FillShowers(
					   AnalysisTreeDataStruct::ShowerDataStruct& showerData,
					   std::vector<recob::Shower> const& showers, const bool fSavePFParticleInfo, 
                      const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
					   ) const {
  
  const size_t NShowers = showers.size();
  
  //
  // prepare the data structures, the tree and the connection between them
  //
  
  // allocate enough space for this number of showers
  // (but at least for one of them!)
  showerData.SetMaxShowers(std::max(NShowers, (size_t) 1));
  showerData.Clear(); // clear all the data
  
  // now set the tree addresses to the newly allocated memory;
  // this creates the tree branches in case they are not there yet
  showerData.SetAddresses(fTree);
  if (NShowers > showerData.GetMaxShowers()) {
    // got this error? it might be a bug,
    // since we are supposed to have allocated enough space to fit all showers
    mf::LogError("AnalysisTree:limits") << "event has " << NShowers
					<< " " << showerData.Name() << " showers, only "
					<< showerData.GetMaxShowers() << " stored in tree";
  }
  
  //
  // now set the data
  //
  // set the record of the number of showers
  // (data structures are already properly resized)
  showerData.nshowers = (Short_t) NShowers;
  
  // set all the showers one by one
  for (size_t i = 0; i < NShowers; ++i) FillShower(showerData, i, showers[i], fSavePFParticleInfo, showerIDtoPFParticleIDMap);
  
} // microboone::AnalysisTree::FillShowers()



void microboone::AnalysisTree::HitsPurity(detinfo::DetectorClocksData const& clockData,
                                          std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe){

  trackid = -1;
  purity = -1;

  std::map<int,double> trkide;

  for(size_t h = 0; h < hits.size(); ++h){

    art::Ptr<recob::Hit> hit = hits[h];
    std::vector<sim::IDE> ides;

    std::vector<sim::TrackIDE> eveIDs = fMCTruthMatching->HitToEveID(clockData, hit);

    for(size_t e = 0; e < eveIDs.size(); ++e){
      //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
      trkide[eveIDs[e].trackID] += eveIDs[e].energy;
    }
  }

  maxe = -1;
  double tote = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      trackid = ii->first;
    }
  }

  //std::cout << "the total energy of this reco track is: " << tote << std::endl;

  if (tote>0){
    purity = maxe/tote;
  }
  //std::cout << "tote = " << tote << ", maxe = " << maxe << ", purity = " << purity << ", trackid =" << trackid <<  std::endl;
}

// Calculate distance to boundary.
double microboone::AnalysisTree::bdist(const TVector3& pos)
{
  // Get geometry.
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC();

  double d1 = pos.X();                             // Distance to right side (wires).
  double d2 = 2.*tpc.HalfWidth() - pos.X();   // Distance to left side (cathode).
  double d3 = pos.Y() + tpc.HalfHeight();     // Distance to bottom.
  double d4 = tpc.HalfHeight() - pos.Y();     // Distance to top.
  double d5 = pos.Z();                             // Distance to front.
  double d6 = tpc.Length() - pos.Z();         // Distance to back.

  double result = std::min({d1, d2, d3, d4, d5, d6});
  return result;
}

// Length of reconstructed track, trajectory by trajectory.
double microboone::AnalysisTree::length(const recob::Track& track)
{
  return track.Length();
}


double microboone::AnalysisTree::driftedLength(detinfo::DetectorPropertiesData const& detProp,
                                               const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom){
  // Get geometry.
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC();

  //compute the drift x range
  double vDrift = detProp.DriftVelocity()*1e-3; //cm/ns
  double xrange[2] = {detProp.ConvertTicksToX(0,0,0,0),detProp.ConvertTicksToX(detProp.NumberTimeSamples(),0,0,0)};
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*tpc.HalfWidth(),-tpc.HalfHeight(),tpc.HalfHeight(),0.,tpc.Length()};

  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(auto step: mctrack) {
    // check if the particle is inside a TPC
    if (step.X() >= bnd[0] && step.X() <= bnd[1] && step.Y() >= bnd[2] && step.Y() <= bnd[3] && step.Z() >= bnd[4] && step.Z() <= bnd[5]){
      // Doing some manual shifting to account for
      // an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from
      // where the particle actually passed to to the time
      // being different from in-spill interactions
      double newX = step.X()+(step.T()*vDrift);
      if (newX < xrange[0] || newX > xrange[1]) continue;
     
      TLorentzVector pos(newX,step.Y(),step.Z(),step.T());
      if(first){
	tpcstart = pos;
	tpcmom = step.Momentum();
	first = false;
      }
      else {
	disp -= pos.Vect();
	result += disp.Mag();
      }
      disp = pos.Vect();
      tpcend = pos;
    }
  }
  return result;
}

// Length of MC particle, trajectory by trajectory (with the manual shifting for x correction)
double microboone::AnalysisTree::driftedLength(detinfo::DetectorPropertiesData const& detProp,
                                               const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  // Get geometry.
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC();
  
  //compute the drift x range
  double vDrift = detProp.DriftVelocity()*1e-3; //cm/ns
  double xrange[2] = {detProp.ConvertTicksToX(0,0,0,0),detProp.ConvertTicksToX(detProp.NumberTimeSamples(),0,0,0)};
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*tpc.HalfWidth(),-tpc.HalfHeight(),tpc.HalfHeight(),0.,tpc.Length()};

  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      // Doing some manual shifting to account for
      // an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from
      // where the particle actually passed to to the time
      // being different from in-spill interactions
      double newX = p.Vx(i)+(p.T(i)*vDrift);
      if (newX < xrange[0] || newX > xrange[1]) continue;
      TLorentzVector pos(newX,p.Vy(i),p.Vz(i),p.T());
      if(first){
	start = pos;
	starti=i;
	first = false;
      }
      else {
	disp -= pos.Vect();
	result += disp.Mag();
      }
      disp = pos.Vect();
      end = pos;
      endi = i;
    }
  }
  return result;
}

// Length of MC particle, trajectory by trajectory (with out the manual shifting for x correction)
double microboone::AnalysisTree::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  // Get geometry.
  auto const& tpc = art::ServiceHandle<geo::Geometry>()->TPC();
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*tpc.HalfWidth(),-tpc.HalfHeight(),tpc.HalfHeight(),0.,tpc.Length()};
  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      if(first){
	start = p.Position(i);
	first = false;
	starti = i;
      }else{
	disp -= p.Position(i).Vect();
	result += disp.Mag();
      }
      disp = p.Position(i).Vect();
      end = p.Position(i);
      endi = i;
    }
  }
  return result;
}



namespace microboone{

  DEFINE_ART_MODULE(AnalysisTree)

}

namespace lar_pandora{class LArPandoraHelper;}
