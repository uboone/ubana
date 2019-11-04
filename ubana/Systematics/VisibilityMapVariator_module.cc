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
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
//#include "larsim/LArG4/OpFastScintillation.hh"
//#include "larsim/LArG4/OpDetPhotonTable.h"

#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TMath.h"

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
    double Gaisser_Hillas(double x,double *par);
    double Gaisser_Hillas_Correction(std::vector<std::vector<double> > fGHvuvpars,std::vector<double> fborder_corr, double j, double distance_to_corner, double fReference_to_corner, double distance );

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
  
    std::vector<std::vector<double> > fGH_PARS_60;
    std::vector<std::vector<double> > fGH_PARS_120;
    std::vector<double> fBORDER_correction;
    bool fOverallSuppression;
    double fOverallSuppressionFactor;
    double fAttenuationLength;
    bool fAttenuationCorr;
    bool fRayleighCorr;
    bool fExtraInfo;
    //double fPMT_radius;

    TTree *ftree;

    /*vector<TH2F*> TwoDs_XY;
    vector<TH2F*> TwoDs_ZY;
    for (unsigned int i = 0; i<32 ; i++) {
       TwoDs_XY.push_back(new TH2F(Form("XY_%d",i),";X;Y", 75, x_min, x_max-((x_max-x_min)/75), 75, y_min, y_max-((y_max-y_min)/75)));
       TwoDs_ZY.push_back(new TH2F(Form("ZY_%d",i),";Z;Y", 400, z_min-((z_max-z_min)/400)/2, z_max-((z_max-z_min)/400)/2, 75, y_min, y_max-((y_max-y_min)/75)));
    }*/

  private:
  
    // Pointer to the feature algorithm
  
    std::ofstream fout;
     
  };
}

namespace phot{

  VisibilityMapVariator::VisibilityMapVariator(fhicl::ParameterSet const & p)
    : EDAnalyzer(p) 
  {
     fGH_PARS_60 = p.get<std::vector<std::vector<double> > >("GH_RS60cm_SBN");
     fGH_PARS_120 = p.get<std::vector<std::vector<double> > >("GH_RS120cm_SBN");
     std::vector<double> v0(2,0.0);
     fBORDER_correction = p.get<std::vector<double>>("BORDER_correction_SBND", v0);   
     fOverallSuppression = p.get<bool>("overallSuppression",true);
     fOverallSuppressionFactor = p.get<double>("overallSuppressionFactor",1.);
     fExtraInfo = p.get<bool>("extrainfo",false);
     fAttenuationCorr = p.get<bool>("AttenuationCorrection",true);
     fRayleighCorr = p.get<bool>("RayleighCorrection",true);
     fAttenuationLength = p.get<double>("AttenuationLength",1000);
  }
  
  void VisibilityMapVariator::analyze(art::Event const & e)
  {
  }
  
  double VisibilityMapVariator::Gaisser_Hillas(double x,double *par) {
    //This is the Gaisser-Hillas function
    double X_mu_0=par[3];
    double Normalization=par[0];
    double Diff=par[1]-X_mu_0;
    double Term=pow((x-X_mu_0)/Diff,Diff/par[2]);
    double Exponential=std::exp((par[1]-x)/par[2]);
    return (Normalization*Term*Exponential);
  }
 
  double VisibilityMapVariator::Gaisser_Hillas_Correction(std::vector<std::vector<double> > fGHvuvpars,std::vector<double> fborder_corr, double j, double distance_to_corner, double fReference_to_corner, double distance ) {
    double pars_ini_[4] = {fGHvuvpars[0][j] + fborder_corr[0] * (distance_to_corner - fReference_to_corner),
                           fGHvuvpars[1][j] + fborder_corr[1] * (distance_to_corner - fReference_to_corner),
                           fGHvuvpars[2][j],
                           fGHvuvpars[3][j]};
    return Gaisser_Hillas(distance, pars_ini_);   
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

    std::string variation;
    if (fAttenuationCorr && fRayleighCorr) variation = "All";
    else if (fAttenuationCorr) variation = "Attenuation";
    else if (fRayleighCorr)    variation = "Rayleigh";
    else if (fOverallSuppression) variation = Form("Suppression_%2F",fOverallSuppressionFactor);
    TFile* fout = TFile::Open(Form("uboone_photon_library_variation_%s.root",variation.c_str()), "recreate");
    fout->mkdir("pmtresponse");
    if (fExtraInfo) fout->mkdir("histos");

    float newVis = 0;
    int newVoxel = 0;
    int newOpChan = 0;
    float newX =0;
    float newY =0;
    float newZ =0;
    float GH_correction_120;
    float GH_correction_60;

    TTree* t = (TTree*) ftree->CloneTree(0);
    t->Reset();
    t->SetBranchAddress("Visibility", &newVis);
    t->SetBranchAddress("Voxel", &newVoxel);
    t->SetBranchAddress("OpChannel", &newOpChan);
    
    TList *l = new TList(); 
    std::vector<TH2F*> h2ds_XY;
    std::vector<TH2F*> h2ds_ZY;
    TBranch *bX = t->Branch("X",&newX,"X/F");
    TBranch *bY = t->Branch("Y",&newY,"Y/F");
    TBranch *bZ = t->Branch("Z",&newZ,"Z/F");
    TBranch *bGH_correction_120 = t->Branch("GH_correction_120",&GH_correction_120,"GH_correction_120/F");
    TBranch *bGH_correction_60 = t->Branch("GH_correction_60",&GH_correction_60,"GH_correction_60/F");

    if (fExtraInfo) { 
 
       for (unsigned int i = 0; i<32 ; i++) {
         h2ds_XY.push_back(new TH2F(Form("XY_%d",i),";X;Y", 75, fXmin, fXmax-((fXmax-fXmin)/75), 75, fYmin, fYmax-((fYmax-fYmin)/75)));
         h2ds_ZY.push_back(new TH2F(Form("ZY_%d",i),";Z;Y", 400, fZmin-((fZmax-fZmin)/400)/2, fZmax-((fZmax-fZmin)/400)/2, 75, fYmin, fYmax-((fYmax-fYmin)/75)));
	 l->Add(h2ds_XY[i]);
	 l->Add(h2ds_ZY[i]);
       }
    }

    art::ServiceHandle<PhotonVisibilityService> pvs;
    //auto const nOpChannels = pvs->NOpChannels();
    for(int iopc = 0; iopc < 32; iopc++)
    {
       unsigned int opdet = iopc; //TODO change to real mapping 
       newOpChan = iopc;
       for(int ivox = 0; ivox != pvs->GetVoxelDef().GetNVoxels(); ivox++)
       {
	  newVoxel = ivox;
          std::vector<int> coor  = pvs->GetVoxelDef().GetVoxelCoords(ivox);
          TVector3    steps = pvs->GetVoxelDef().GetSteps();
    	  double const xyz[3] = {
            fXmin+((fXmax - fXmin)/steps.X())*coor[0],
            fYmin+((fYmax - fYmin)/steps.Y())*coor[1],
            fZmin+((fZmax - fZmin)/steps.Z())*coor[2]
          };
	  
	  newX = xyz[0];
	  newY = xyz[1];
	  newZ = xyz[2];
	 
          float vis = pvs->GetVisibility(xyz,iopc,false);
          newVis = float(vis);
	  
	  //Overall Suppression
          if (fOverallSuppression) {
    	     newVis = newVis * fOverallSuppressionFactor;
	  }
	
          //Attenuation length correction
	  double distance = pvs->DistanceToOpDet(xyz,opdet);
	  if (fAttenuationCorr) {
	     newVis =  newVis * exp(-1.*distance / fAttenuationLength);
	  }

	  //Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
	  if (fRayleighCorr) {
	     double costheta = pvs->SolidAngleFactor(xyz,opdet);
	     double theta = acos(costheta)*180./CLHEP::pi;
	     double j = theta / 10; // angular bin width = 10
	     double fYactive_corner = (fYmax - fYmin)/2;
             double fZactive_corner = (fZmax - fZmin)/2;
             double fReference_to_corner = sqrt(pow(fYactive_corner,2) + pow(fZactive_corner,2));
	     double z_to_corner = abs(xyz[2] - fZactive_corner) - fZactive_corner;
    	     double y_to_corner = abs(xyz[1]) - fYactive_corner;
    	     double distance_to_corner = sqrt(y_to_corner*y_to_corner + z_to_corner*z_to_corner);// in the ph-cathode plane
	     std::vector<double> fborder_corr;
	     std::vector<std::vector<double> > fGHvuvpars;
             double pars_ini_60[4] = {fGH_PARS_60[0][j] + fBORDER_correction[0] * (distance_to_corner - fReference_to_corner),
                                    fGH_PARS_60[1][j] + fBORDER_correction[1] * (distance_to_corner - fReference_to_corner),
                                    fGH_PARS_60[2][j],
                                    fGH_PARS_60[3][j]};
             GH_correction_60 = Gaisser_Hillas(distance, pars_ini_60);
             double pars_ini_120[4] = {fGH_PARS_120[0][j] + fBORDER_correction[0] * (distance_to_corner - fReference_to_corner),
                                    fGH_PARS_120[1][j] + fBORDER_correction[1] * (distance_to_corner - fReference_to_corner),
                                    fGH_PARS_120[2][j],
                                    fGH_PARS_120[3][j]};
             GH_correction_120 = Gaisser_Hillas(distance, pars_ini_120);
	     if (GH_correction_60 > 0) newVis = newVis * GH_correction_120 / GH_correction_60; 
	  }

	  ///////////////////////////////////////////
	  if (fExtraInfo) {  
   	     bX->Fill();        
   	     bY->Fill();        
   	     bZ->Fill();        
    	     bGH_correction_120->Fill();
    	     bGH_correction_60->Fill();
	     h2ds_XY[iopc]->Fill(newX, newY, newVis); 
	     h2ds_ZY[iopc]->Fill(newZ, newY, newVis); 
	  }
	  ///////////////////////////////////////////
          t->Fill();

       } 
    }
    fout->cd("pmtresponse"); ///PhotonLibraryData
    t->Write();
    if (fExtraInfo) {
       fout->cd("histos");
       l->Write();
    }
    fout->Write();
    fout->Close();
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
