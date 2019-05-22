#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVAClassificationApplication_protonID( TString myMethodList = "BDTG" )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           trackid;
   Float_t         nhits;
   Float_t         length;
   Float_t         startx;
   Float_t         starty;
   Float_t         startz;
   Float_t         endx;
   Float_t         endy;
   Float_t         endz;
   Float_t         theta;
   Float_t         phi;
   Float_t         distlenratio;
   Float_t         startdedx;
   Float_t         dedxratio;
   Float_t         trtotaldedx;
   Float_t         traveragededx;
   Float_t         cosmicscore;
   Float_t         mc_ccnc;
   Float_t         mc_mode;
   Float_t         mc_hitnuc;
   Float_t         mc_q2;
   Float_t         mc_nux;
   Float_t         mc_nuy;
   Float_t         mc_nuz;
   Float_t         mc_enu;
   Float_t         mcpdg;
   Float_t         mcprimary;
   Float_t         mcorigin;
   Float_t         mclength;
   Float_t         mcstartx;
   Float_t         mcstarty;
   Float_t         mcstartz;
   Float_t         mcendx;
   Float_t         mcendy;
   Float_t         mcendz;
   Float_t         mctheta;
   Float_t         mcphi;
   Float_t         mckinetic;


   reader->AddVariable( "nhits", &nhits );
   reader->AddVariable( "distlenratio", &distlenratio );
   reader->AddVariable( "length",                &length );
   reader->AddVariable( "theta",                &theta);
   reader->AddVariable( "phi",                &phi );
   reader->AddVariable( "startdedx",                &startdedx);
   reader->AddVariable( "dedxratio",                &dedxratio );
   reader->AddVariable( "traveragededx",                &traveragededx);
   reader->AddVariable( "trtotaldedx",                &trtotaldedx);
   reader->AddVariable( "starty",                &starty );
   reader->AddVariable( "startz",                &startz );
   reader->AddVariable( "endy",                &endy );
   reader->AddVariable( "endz",                &endz);
   reader->AddVariable( "cosmicscore",                &cosmicscore );
   reader->AddSpectator( "mcpdg",  &mcpdg);
   reader->AddSpectator( "startx",  &startx);
   reader->AddSpectator( "endx",  &endx );
   reader->AddSpectator( "run",  &run );
   reader->AddSpectator( "subrun",  &subrun);
   reader->AddSpectator( "event",  &event);
   reader->AddSpectator( "trackid",  &trackid );
   reader->AddSpectator( "mc_ccnc",&mc_ccnc);
   reader->AddSpectator( "mc_mode",&mc_mode);
   reader->AddSpectator( "mc_hitnuc",&mc_hitnuc);
   reader->AddSpectator( "mc_q2",&mc_q2);
   reader->AddSpectator( "mc_nux",&mc_nux);
   reader->AddSpectator( "mc_nuy",&mc_nuy);
   reader->AddSpectator( "mc_nuz",&mc_nuz);
   reader->AddSpectator( "mc_enu",&mc_enu);
   reader->AddSpectator( "mcpdg",&mcpdg);
   reader->AddSpectator( "mcprimary",&mcprimary);
   reader->AddSpectator( "mcorigin",&mcorigin);
   reader->AddSpectator( "mclength",&mclength);
   reader->AddSpectator( "mcstartx",&mcstartx);
   reader->AddSpectator( "mcstarty", &mcstarty);
   reader->AddSpectator( "mcstartz",&mcstartz);
   reader->AddSpectator( "mcendx",&mcendx);
   reader->AddSpectator( "mcendy", &mcendy);
   reader->AddSpectator( "mcendz",&mcendz);
   reader->AddSpectator( "mctheta", &mctheta);
   reader->AddSpectator( "mcphi",&mcphi);
   reader->AddSpectator( "mckinetic",&mckinetic);


   // Book the MVA methods

   TString dir    = "output/";
   TString prefix = "TMVAClassification_BDTG_pmtrack";
 reader->BookMVA( "BDTG", "./output/TMVAClassification_BDTG_pmtrack.weights.xml");
  
   // Book output histograms
   UInt_t nbin = 100;
  
   TH1F *histBdtG(0);
  

   histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
 


   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   TString fname = "/uboone/data/users/renlu23/book/April_13th/NCE_MCC9_v12/merged/overlay_testing_pmtrack.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD"); // if not: download from ROOT server
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("decisiontreeidana/tree");
  
   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   //  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
  for (Long64_t ievt=0; ievt<theTree->500000;ievt++) {

      if (ievt%100000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      //   var1 = userVar1 + userVar2;
      //var2 = userVar1 - userVar2;

      // Return the MVA outputs and fill into histograms

   
     
      histBdtG   ->Fill( reader->EvaluateMVA( "BDTG"          ) );
    

   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

 

   // Write histograms

   TFile *target  = new TFile( "TMVA_pmtrack_testing.root","RECREATE" );
  
   histBdtG ->Write();
  

   std::cout << "--- Created root file: \"TMVA_bnb_data.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

