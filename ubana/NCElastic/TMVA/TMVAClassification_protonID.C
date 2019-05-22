#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification_protonID( TString myMethodList = "BDTG" )
{
  TMVA::Tools::Instance();

  std::map<std::string,int> Use;

  
  // Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  //
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification- BDTG" << std::endl;

  // Read training and test data
  // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
  TFile *input_pmtrack(0);
  TFile *input_pandora(0);
  TString fname_pandora = "/uboone/data/users/renlu23/book/April_13th/NCE_MCC9_v12/merged/overlay_training_pmtrack.root";
  //TString fname_pmtrack = "/uboone/data/users/renlu23/book/April_13th/NCE_MCC9_v12/merged/overlay_training_pandora.root";
  TString fname_pmtrack ="/uboone/data/users/renlu23/book/May_7th/NCE_MCC9_v12/overlay_training_trajcluster/overlay_training_trajcluster.root";
  if (!gSystem->AccessPathName( fname_pmtrack )) {
    input_pmtrack = TFile::Open( fname_pmtrack ); // check if file in local directory exists
    input_pandora = TFile::Open( fname_pandora ); // check if file in local directory exists
  }
   else {
     TFile::SetCacheFileDir("./");
     input_pmtrack = TFile::Open("trf.root", "CACHEREAD");
   }
  if (!input_pmtrack) {
    std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
  }

  std::cout << "--- TMVAClassification       : Using input file: " << input_pmtrack->GetName() << std::endl;
  
  // Register the training and test trees
  TTree *signalTree     = (TTree*)input_pmtrack->Get("decisiontreeidana/tree");
  TTree *background     = (TTree*)input_pmtrack->Get("decisiontreeidana/tree");
  
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "TMVA_trajcluster.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   dataloader->AddVariable( "nhits", 'F' );
   dataloader->AddVariable( "distlenratio", "straightness", "", 'F' );
   dataloader->AddVariable( "length",                "length", "cm", 'F' );
   dataloader->AddVariable( "theta",                "theta", "radian", 'F' );
   dataloader->AddVariable( "phi",                "phi", "radian", 'F' );
   dataloader->AddVariable( "startdedx",                "startdedx", "", 'F' );
   dataloader->AddVariable( "dedxratio",                "dedxratio", "", 'F' );
   dataloader->AddVariable( "traveragededx",                "traveragededx", "", 'F' );
   dataloader->AddVariable( "trtotaldedx",                "trtotaldedx", "", 'F' );
   dataloader->AddVariable( "starty",                "starty", "cm", 'F' );
   dataloader->AddVariable( "startz",                "startz", "cm", 'F' );
   dataloader->AddVariable( "endy",                "endy", "cm", 'F' );
   dataloader->AddVariable( "endz",                "endz", "cm", 'F' );
//   dataloader->AddVariable( "cosmicscore",                "cosmicscore", "", 'F' );
    
   dataloader->AddSpectator( "mcpdg",  "mcpdg", "", 'F' );
   dataloader->AddSpectator( "startx",  "startx", "cm", 'F' );
   dataloader->AddSpectator( "endx",  "endx", "cm", 'F' );
   dataloader->AddSpectator( "run",  "run", "", 'I' );
   dataloader->AddSpectator( "subrun",  "subrun", "", 'I' );
   dataloader->AddSpectator( "event",  "event", "", 'I' );
   dataloader->AddSpectator( "trackid",  "trackid", "", 'I' );
   dataloader->AddSpectator( "mc_ccnc","mc_ccnc","",'F');
   dataloader->AddSpectator( "mc_mode","mc_mode","",'F');
   dataloader->AddSpectator( "mc_hitnuc","mc_hitnuc","",'F');
   dataloader->AddSpectator( "mc_q2","mc_q2","",'F');
   dataloader->AddSpectator( "mc_nux","mc_nux","",'F');
   dataloader->AddSpectator( "mc_nuy","mc_nuy","",'F');
   dataloader->AddSpectator( "mc_nuz","mc_nuz","",'F');
   dataloader->AddSpectator( "mc_enu","mc_enu","",'F');
   dataloader->AddSpectator( "mcpdg","mcpdg","",'F');
   dataloader->AddSpectator( "mcprimary","mcprimary","",'F');
   dataloader->AddSpectator( "mcorigin","mcorigin","",'F');
   dataloader->AddSpectator( "mclength","mclength","",'F');
   dataloader->AddSpectator( "mcstartx","mcstartx","",'F');
   dataloader->AddSpectator( "mcstarty", "mcstarty","",'F');
   dataloader->AddSpectator( "mcstartz","mcstartz","",'F');
   dataloader->AddSpectator( "mcendx","mcendx","",'F');
   dataloader->AddSpectator( "mcendy", "mcendy","",'F');
   dataloader->AddSpectator( "mcendz","mcendz","",'F');
   dataloader->AddSpectator( "mctheta", "mctheta","",'F');
   dataloader->AddSpectator( "mcphi","mcphi","",'F');
   dataloader->AddSpectator( "mckinetic","mckinetic","",'F');
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

  
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

 
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "mcpdg==2212&&cosmicscore!=1&&nhits>0&&length>0&&distlenratio>0&&theta>-9999&&startdedx>-9999&&dedxratio>-9999&&traveragededx>-9999&&trtotaldedx>-9999&&starty>-100&&startz>0&&endy>-100&&endz>0&&starty<100&&startz<1000&&endy<100&&endz<1000&&startdedx<200&&dedxratio<10&&traveragededx<30&&trtotaldedx<1000&&startx>10&&endx<248"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = "mcpdg!=2212&&cosmicscore!=1&&nhits>0&&length>0&&distlenratio>0&&theta>-9999&&startdedx>-9999&&dedxratio>-9999&&traveragededx>-9999&&trtotaldedx>-9999&&starty>-100&&startz>0&&endy>-100&&endz>0&&starty<100&&startz<1000&&endy<100&&endz<1000&&startdedx<200&&dedxratio<10&&traveragededx<30&&trtotaldedx<1000&&startx>10&&endx<248";

   dataloader->PrepareTrainingAndTestTree( mycuts,mycutb,"SplitMode=Random:!V" );
      //   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );


   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=300:NegWeightTreatment=PairNegWeightsGlobal:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.75:nCuts=20:MaxDepth=6" );


     // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   // Save the output
   outputFile->Close();
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

