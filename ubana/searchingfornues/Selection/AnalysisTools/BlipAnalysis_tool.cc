#ifndef ANALYSIS_BLIPANALYSIS_CXX
#define ANALYSIS_BLIPANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "ubana/searchingfornues/Selection/CommonDefs/Typedefs.h"

// backtracking tools
#include "ubana/searchingfornues/Selection/CommonDefs/PIDFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/SCECorrections.h"
#include "ubana/searchingfornues/Selection/CommonDefs/Geometry.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "ubobj/Blip/DataTypes.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       BlipAnalysis
// File:        BlipAnalysis.cc
//
//  This analysis tool reads MeV-scale blip information from files
//  (blipobj::Blip data products) and saves them to the output ntuple.
//  To save space, only blips within some configurable distance from 
//  the identified neutrino vertex will be saved to the ntuple.
//
//
// Configuration parameters:
//
//  BlipProducer      : producer of blip collection saved to event to be read
//  NuBlipRadius      : radius of acceptance around neutrino vertex
//  SaveSCECorrLoc    : save the SCE-corrected blip location
//  SaveSCECorrEnergy : save the SCE- and lifetime-corrected charge and energy
//
// Original tool created by Afroditi Papadopoulou & Burke Irwin (burke.irwin7@gmail.com) on 01/24/2023.
// Updated for MCC10 by Will Foreman (wforeman.phys@gmail.com) on 01/28/2025.
//
////////////////////////////////////////////////////////////////////////

class BlipAnalysis : public AnalysisToolBase
{

public:
  
  /// @brief  Constructor
  BlipAnalysis(const fhicl::ParameterSet &pset);

  /// @brief  Destructor
  ~BlipAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /// @brief Analysis function
  void analyzeEvent(art::Event const &e, bool fData) override;

  /// @brief Analyze slice
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;
     
  /// @brief Method for adding a blip's info to the branch vectors
  void addTheBlip(art::Ptr<blipobj::Blip> blip);

  /// @brief set branches for TTree
  void setBranches(TTree *_tree) override;

  /// @brief reset ttree branches
  void resetTTree(TTree *_tree) override;
  
  /// @brief reset ttree branches
  void resetVariables();

  /// @brief truncate float variable to desired decimal precision
  float truncate(float input, double base);

private:

  // --- fcl parameters ---
  art::InputTag fBlipProducer;// blip collection producer
  float   fNuBlipRadius;      // cm (set <= 0 --> saves all blips in the event)
  bool    fSaveSCECorrLoc;    // option to save SCE-corrected location
  bool    fSaveSCECorrEnergy; // option to save SCE- and lifetime-corrected energy
  bool    fSaveOnlyNuEvts;    // only save blips for neutrino PFP events
  bool    fLiteMode;          // save bare minimum of variables

  // --- event identifiers --- 
  int _run, _sub, _evt;
  
  // --- 3D blip information ---
  int _nblips;                          // number of blips found in this event
  std::vector<int>    _blip_id;         // unique blip ID
  std::vector<int>    _blip_tpc;        // blip TPC
  std::vector<int>    _blip_nplanes;    // number of planes matched (2 or 3)
  std::vector<float>  _blip_proxtrkdist;// distance to closest track
  std::vector<int>    _blip_proxtrkid;  // track ID of closest track
  std::vector<bool>   _blip_touchtrk;   // is blip touching a track (ie, delta-like?)
  std::vector<int>    _blip_touchtrkid; // track ID of touched track
  std::vector<bool>   _blip_bydeadwire; // is blip adjacent to a dead channel on coll plane?
  std::vector<float>  _blip_badwirefrac;// frac of blip's wires (across all planes) deemed bad/noisy
  std::vector<float>  _blip_x;          // reconstructed X [cm]
  std::vector<float>  _blip_y;          // reconstructed Y [cm]
  std::vector<float>  _blip_z;          // reconstructed Z [cm]
  std::vector<float>  _blip_sigmayz;    // difference in central wire intersection points
  std::vector<float>  _blip_dx;         // length along drift direction [cm]
  std::vector<float>  _blip_dw;         // length projected onto axis perp. to wire orientation
  std::vector<float>  _blip_size;       // projected size in 3D (based on dx and dw)
  std::vector<float>  _blip_charge;     // charge reconstructed at calo plane (usually collection) [e-]
  std::vector<float>  _blip_energy;     // reconstructed energy [MeVee]
  std::vector<int>    _blip_true_g4id;  // truth-matched G4 track ID
  std::vector<int>    _blip_true_pdg;   // truth-matched PDG
  std::vector<float>  _blip_true_energy;// true energy deposited
  // plane-specific information
  std::vector<int>    _blip_pl0_nwires;       // number of wires on this plane
  std::vector<int>    _blip_pl1_nwires;       // number of wires on this plane
  std::vector<int>    _blip_pl2_nwires;       // number of wires on this plane
  std::vector<int>    _blip_pl0_nwiresbad;    // number of wires deemed noisy
  std::vector<int>    _blip_pl1_nwiresbad;    // number of wires deemed noisy
  std::vector<int>    _blip_pl2_nwiresbad;    // number of wires deemed noisy
  std::vector<bool>   _blip_pl0_bydeadwire;   // is adjacent to a dead channel?
  std::vector<bool>   _blip_pl1_bydeadwire;   // is adjacent to a dead channel?
  std::vector<bool>   _blip_pl2_bydeadwire;   // is adjacent to a dead channel?
  std::vector<int>    _blip_pl0_centerwire;   // central wire number
  std::vector<int>    _blip_pl1_centerwire;   // central wire number
  std::vector<int>    _blip_pl2_centerwire;   // central wire number

};


//----------------------------------------------------------------------------
BlipAnalysis::BlipAnalysis(const fhicl::ParameterSet &p)
{
  fBlipProducer       = p.get<art::InputTag>("BlipProducer",        "blipreco");
  fSaveOnlyNuEvts     = p.get<bool>         ("SaveOnlyNuEvts",      true);
  fNuBlipRadius       = p.get<float>        ("NuBlipRadius",        300.);
  fSaveSCECorrLoc     = p.get<bool>         ("SaveSCECorrLocation", true);
  fSaveSCECorrEnergy  = p.get<bool>         ("SaveSCECorrEnergy",   true);
  fLiteMode           = p.get<bool>         ("LiteMode",            false);
}


//----------------------------------------------------------------------------
void BlipAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//---------------------------------------------------------------------------
float BlipAnalysis::truncate(float input, double base){
  if( base > 0 ) return roundf( input / base ) * base;
  else return input;
}

//----------------------------------------------------------------------------
void BlipAnalysis::addTheBlip(art::Ptr<blipobj::Blip> blip){
    
    // get reconstructed position
    //TVector3 loc;
    TVector3 loc = (fSaveSCECorrLoc) ? blip->PositionSCE : blip->Position;
    
    // get reconstructed charge and energy
    float energy = -9; float charge = -9;
    if( fSaveSCECorrEnergy ) {
      energy = blip->Energy; 
      charge = blip->Charge;
    } else {
      energy = blip->EnergyCorr;  
      charge = blip->ChargeCorr;
    }

    float energyTrue = blip->truth.Energy;
    float size = sqrt( pow(blip->dX,2) + pow(blip->dYZ,2) );
    float x = loc.X();
    float y = loc.Y();
    float z = loc.Z();

    // bad wire frac
    int nwirestot = 0;
    int nwiresbad = 0;
    for(int j=0; j<3; j++){
      int nb = std::max(0,blip->clusters[j].NWiresBad);
      int nn = std::max(0,blip->clusters[j].NWiresNoisy);
      nwiresbad += (nb + nn);
      nwirestot = blip->clusters[j].NWires;
    }
    float badwirefrac = float(nwiresbad)/nwirestot;

    // truncate some floats for reduced data size
    x       = truncate(x,0.01);
    y       = truncate(y,0.01);
    z       = truncate(z,0.01);
    energy  = truncate(energy,0.001);
    charge  = truncate(charge,10);
    size    = truncate(size,0.01);
    energyTrue = truncate(energyTrue,0.001);


    bool isTouchTrk = (blip->TouchTrkID >= 0 );
    _blip_id          .push_back(blip->ID);
    _blip_tpc         .push_back(blip->TPC);
    _blip_nplanes     .push_back(blip->NPlanes);    
    _blip_proxtrkdist .push_back(blip->ProxTrkDist);
    _blip_proxtrkid   .push_back(blip->ProxTrkID);
    _blip_touchtrk    .push_back(isTouchTrk);
    _blip_touchtrkid  .push_back(blip->TouchTrkID);
    _blip_bydeadwire  .push_back(blip->clusters[2].DeadWireSep==0);
    _blip_badwirefrac .push_back(badwirefrac);
    _blip_x           .push_back(x);
    _blip_y           .push_back(y);
    _blip_z           .push_back(z);
    _blip_sigmayz     .push_back(blip->SigmaYZ);
    _blip_dx          .push_back(blip->dX);
    _blip_dw          .push_back(blip->dYZ);
    _blip_size        .push_back(size);
    _blip_charge      .push_back(charge);
    _blip_energy      .push_back(energy);
    _blip_true_g4id   .push_back(blip->truth.LeadG4ID);
    _blip_true_pdg    .push_back(blip->truth.LeadG4PDG);
    _blip_true_energy .push_back(energyTrue);
    _blip_pl0_nwires  .push_back(blip->clusters[0].NWires);
    _blip_pl1_nwires  .push_back(blip->clusters[1].NWires);
    _blip_pl2_nwires  .push_back(blip->clusters[2].NWires);
    _blip_pl0_nwiresbad  .push_back(blip->clusters[0].NWiresBad);
    _blip_pl1_nwiresbad  .push_back(blip->clusters[1].NWiresBad);
    _blip_pl2_nwiresbad  .push_back(blip->clusters[2].NWiresBad);
    _blip_pl0_bydeadwire  .push_back((blip->clusters[0].DeadWireSep==0));
    _blip_pl1_bydeadwire  .push_back((blip->clusters[1].DeadWireSep==0));
    _blip_pl2_bydeadwire  .push_back((blip->clusters[2].DeadWireSep==0));
    _blip_pl0_centerwire  .push_back(blip->clusters[0].CenterWire); 
    _blip_pl1_centerwire  .push_back(blip->clusters[1].CenterWire);
    _blip_pl2_centerwire  .push_back(blip->clusters[2].CenterWire);

}

//----------------------------------------------------------------------------
void BlipAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[BlipAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;

  //================================================
  // If a non-zero radius was configured, or we are
  // only saving neutrino events, then exit out of this 
  // function and wait for 'analyzeSlice'
  //================================================
  if( fSaveOnlyNuEvts   ) return;
  if( fNuBlipRadius > 0 ) return;
  
  //==========================================
  // ... otherwise, save all blips in the event
  //==========================================
  // reset branch vectors
  resetVariables(); 
  // loop over blips saved to the event
  art::Handle< std::vector<blipobj::Blip> > blipHandle;
  std::vector<art::Ptr<blipobj::Blip> > bliplist;
  if (e.getByLabel(fBlipProducer,blipHandle))
    art::fill_ptr_vector(bliplist, blipHandle);
  _nblips = bliplist.size();
  std::cout<<"Saving all "<<_nblips<<" blips in the event\n"; 
  for(size_t i=0; i<bliplist.size(); i++)
    addTheBlip(bliplist[i]);

}


//----------------------------------------------------------------------------
void BlipAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  
  //================================================
  // If NO radius was configured, it means the user 
  // wants to save all blips, which would have already 
  // been done in the call to 'analyzeEvent' unless 
  // the "SaveOnlyNuEvts" flag is enabled.
  //================================================
  if( !fSaveOnlyNuEvts   ) return;
  if( fNuBlipRadius <= 0 ) fNuBlipRadius = 99999.;

  std::cout
  <<"********** BlipAnalysis::analyzeSlice **********\n"
  <<" saving only blips < "<<fNuBlipRadius<<" cm from neutrino vtx\n";

  // reset branch vectors
  resetVariables();

  // get the neutrino vertex; skip if none found
  TVector3 nuvtx;
  bool foundVertex = false;
  for (auto pfp : slice_pfp_v){
    if (!pfp->IsPrimary()) continue;
    double xyz[3] = {};
    auto vtx = pfp.get<recob::Vertex>();
    if (vtx.size() != 1) break;
    vtx.at(0)->XYZ(xyz);
    nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
    foundVertex = true;
    break;
  }
  if( !foundVertex ) return;
  std::cout
  <<" nuvtx position: ("<<nuvtx.X()<<", "<<nuvtx.Y()<<", "<<nuvtx.Z()<<") cm\n";

  // add only blips around vertex
  art::Handle< std::vector<blipobj::Blip> > blipHandle;
  std::vector<art::Ptr<blipobj::Blip> > bliplist;
  if (e.getByLabel(fBlipProducer,blipHandle))
    art::fill_ptr_vector(bliplist, blipHandle);
  for(size_t i=0; i<bliplist.size(); i++){
    auto& blip = bliplist[i];
    //TVector3 loc;
    TVector3 loc = (fSaveSCECorrLoc) ? blip->PositionSCE : blip->Position;
    if( (loc-nuvtx).Mag() > fNuBlipRadius ) continue;
    addTheBlip(blip);
    _nblips++;
  }
  std::cout
  <<" total blips in the event: "<<bliplist.size()<<"\n"
  <<" blips saved: "<<_nblips<<"\n"
  <<"************************************************\n";

}


//----------------------------------------------------------------------------
void BlipAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("nblips_saved",           &_nblips,         "nblips_saved/I");
  //_tree->Branch("blip_id",        "std::vector< int >",  &_blip_id);
  _tree->Branch("blip_x",           "std::vector< float >",   &_blip_x);
  _tree->Branch("blip_y",           "std::vector< float >",   &_blip_y);
  _tree->Branch("blip_z",           "std::vector< float >",   &_blip_z);
  _tree->Branch("blip_size",        "std::vector< float >",   &_blip_size);
  _tree->Branch("blip_energy",      "std::vector< float >",   &_blip_energy);
  if( !fLiteMode ) {
  _tree->Branch("blip_charge",      "std::vector< float >",   &_blip_charge);
  _tree->Branch("blip_nplanes",     "std::vector< int >",     &_blip_nplanes);
  _tree->Branch("blip_proxtrkdist", "std::vector< float >",   &_blip_proxtrkdist);
  _tree->Branch("blip_proxtrkid",   "std::vector< int >",     &_blip_proxtrkid);
  _tree->Branch("blip_touchtrk",    "std::vector< bool >",    &_blip_touchtrk);
  _tree->Branch("blip_touchtrkid",  "std::vector< int >",     &_blip_touchtrkid);
  //_tree->Branch("blip_bydeadwire",  "std::vector< bool >",    &_blip_bydeadwire);
  _tree->Branch("blip_badwirefrac",  "std::vector< float >",    &_blip_badwirefrac);
  _tree->Branch("blip_pl0_nwires",  "std::vector< int >",     &_blip_pl0_nwires);
  _tree->Branch("blip_pl1_nwires",  "std::vector< int >",     &_blip_pl1_nwires);
  _tree->Branch("blip_pl2_nwires",  "std::vector< int >",     &_blip_pl2_nwires);
  //_tree->Branch("blip_pl0_nwiresbad", "std::vector< int >",   &_blip_pl0_nwiresbad);
  //_tree->Branch("blip_pl1_nwiresbad", "std::vector< int >",   &_blip_pl1_nwiresbad);
  //_tree->Branch("blip_pl2_nwiresbad", "std::vector< int >",   &_blip_pl2_nwiresbad);
  _tree->Branch("blip_pl0_bydeadwire","std::vector< bool >",  &_blip_pl0_bydeadwire);
  _tree->Branch("blip_pl1_bydeadwire","std::vector< bool >",  &_blip_pl1_bydeadwire);
  _tree->Branch("blip_pl2_bydeadwire","std::vector< bool >",  &_blip_pl2_bydeadwire);
  _tree->Branch("blip_pl0_centerwire","std::vector< int >",   &_blip_pl0_centerwire);
  _tree->Branch("blip_pl1_centerwire","std::vector< int >",   &_blip_pl1_centerwire);
  _tree->Branch("blip_pl2_centerwire","std::vector< int >",   &_blip_pl2_centerwire);
  }
  _tree->Branch("blip_true_g4id",   "std::vector< int >",     &_blip_true_g4id);
  _tree->Branch("blip_true_pdg",    "std::vector< int >",     &_blip_true_pdg);
  if( !fLiteMode ){ 
  _tree->Branch("blip_true_energy", "std::vector< float >",   &_blip_true_energy);
  }
}


//----------------------------------------------------------------------------
void BlipAnalysis::resetTTree(TTree *_tree)
{
  resetVariables();
}


//----------------------------------------------------------------------------
void BlipAnalysis::resetVariables()
{
    _nblips          = 0;
    _blip_id          .clear(); 
    _blip_tpc         .clear(); 
    _blip_nplanes     .clear(); 
    _blip_proxtrkdist .clear(); 
    _blip_proxtrkid   .clear(); 
    _blip_touchtrk    .clear(); 
    _blip_touchtrkid  .clear();
    _blip_bydeadwire  .clear();
    _blip_badwirefrac.clear();
    _blip_x           .clear(); 
    _blip_y           .clear(); 
    _blip_z           .clear(); 
    _blip_sigmayz     .clear(); 
    _blip_dx          .clear();  
    _blip_dw          .clear(); 
    _blip_size        .clear(); 
    _blip_charge      .clear(); 
    _blip_energy      .clear(); 
    _blip_true_g4id   .clear(); 
    _blip_true_pdg    .clear(); 
    _blip_true_energy .clear(); 
    _blip_pl0_nwires  .clear(); 
    _blip_pl1_nwires  .clear(); 
    _blip_pl2_nwires  .clear(); 
    _blip_pl0_nwiresbad .clear();
    _blip_pl1_nwiresbad .clear();
    _blip_pl2_nwiresbad .clear();
    _blip_pl0_bydeadwire.clear();
    _blip_pl1_bydeadwire.clear();
    _blip_pl2_bydeadwire.clear();
    _blip_pl0_centerwire.clear();
    _blip_pl1_centerwire.clear();
    _blip_pl2_centerwire.clear();
}


DEFINE_ART_CLASS_TOOL(BlipAnalysis)
} // namespace analysis

#endif
