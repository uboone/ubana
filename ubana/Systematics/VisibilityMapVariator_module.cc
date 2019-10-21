////////////////////////////////////////////////////////////////////////
// Class:       VisibilityMapVariator
// Plugin Type: analyze (art v2_05_00)
// File:        VisibilityMapVariator_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcore/Geometry/Geometry.h"

#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include <memory>
#include <cmath>
#include <stdexcept>

class VisibilityMapVariator;

namespace phot {
  class VisibilityMapVariator : public art::EDAnalyzer {
  public:
    explicit VisibilityMapVariator(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
  
    // Plugins should not be copied or assigned.
    VisibilityMapVariator(VisibilityMapVariator const &) = delete;
    VisibilityMapVariator(VisibilityMapVariator &&) = delete;
    VisibilityMapVariator & operator = (VisibilityMapVariator const &) = delete;
    VisibilityMapVariator & operator = (VisibilityMapVariator &&) = delete;
  
    // Required functions.
    void beginJob() override;
    void analyze(art::Event const & e) override;
    void init(TTree *tree); 

    Int_t           Voxel;
    Int_t           OpChannel;
    Float_t         Visibility;
    TBranch        *b_Voxel;   //!
    TBranch        *b_OpChannel;   //!
    TBranch        *b_Visibility;   //!

    float fXmin;  
    float fXmax;  
    float fYmin;  
    float fYmax;  
    float fZmin;  
    float fZmax;  
   
    TTree *ftree;

  private:
  
    // Pointer to the feature algorithm
  
    std::ofstream fout;
     
  };
}

namespace phot{

  VisibilityMapVariator::VisibilityMapVariator(fhicl::ParameterSet const & p)
    : EDAnalyzer(p) 
  {
  
  }
  
  void VisibilityMapVariator::analyze(art::Event const & e)
  {
  }
  
  
  void VisibilityMapVariator::beginJob()
  {
    art::ServiceHandle<geo::Geometry> geom;
    double CryoBounds[6];
    geom->CryostatBoundaries(CryoBounds);
    fXmin = CryoBounds[0];
    fXmax = CryoBounds[1];
    fYmin = CryoBounds[2];
    fYmax = CryoBounds[3];
    fZmin = CryoBounds[4];
    fZmax = CryoBounds[5];
   
    TFile *f = TFile::Open("/uboone/app/users/adi/Systematics/LY/variation/uboone_photon_library_v6_70kV.root");
    TDirectoryFile *d = (TDirectoryFile*)f->Get("pmtresponse");
    ftree = (TTree*)d->Get("PhotonLibraryData"); 
    init(ftree);

    TFile* fout = TFile::Open("uboone_photon_library_variation.root", "recreate");
    fout->mkdir("pmtresponse");
    fout->cd("pmtresponse"); ///PhotonLibraryData

    double xyz[3];
    float newVis = 0;
    TTree* t = (TTree*) ftree->CloneTree(0);
    t->Reset();
    t->SetBranchAddress("Visibility", &newVis);

    art::ServiceHandle<phot::PhotonVisibilityService const> pvs; 
    //for(size_t iopc = 0; iopc != pvs->NOpChannels(); iopc++)
    for(size_t iopc = 0; iopc != 32; iopc++)
    {
       for(int ivox = 0; ivox != pvs->GetVoxelDef().GetNVoxels(); ivox++)
       {
          std::vector<int> coor  = pvs->GetVoxelDef().GetVoxelCoords(ivox);
          TVector3    steps = pvs->GetVoxelDef().GetSteps();
          xyz[0] = fXmin+((fXmax - fXmin)/steps.X())*coor[0];
          xyz[1] = fYmin+((fYmax - fYmin)/steps.Y())*coor[1];
          xyz[2] = fZmin+((fZmax - fZmin)/steps.Z())*coor[2];

          float const vis = pvs->GetVisibility(xyz,ivox,true);
          newVis = vis*10;
          t->Fill();
       } 
    }
    fout->Write();
  }
  void VisibilityMapVariator::init(TTree *tree) 
  {
   if (!tree) return;
   ftree = tree;

   ftree->SetBranchAddress("Voxel", &Voxel, &b_Voxel);
   ftree->SetBranchAddress("OpChannel", &OpChannel, &b_OpChannel);
   ftree->SetBranchAddress("Visibility", &Visibility, &b_Visibility);
  }


}

namespace phot {
  DEFINE_ART_MODULE(VisibilityMapVariator)
}
