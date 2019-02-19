/*************************************************************
 *
 * This is a quick macro for counting how many daughter tracks a
 * reconstructed track has, using the number of tracks with start
 * or end points within a given distance of the end point of the
 * first track.
 *
 * Don't forget to set up gallery first:
 * setup gallery v1_03_08 -q e10:nu:prof
 *
 * Then compile:
 * make clean; make checkcalib
 *
 * ./CheckCalib "path/to/reco2/file.root or path/to/list/of/input/files.txt"
 *
 * Kirsty Duffy (kduffy@fnal.gov), Fermilab, Jan 28 2018
 *
 *************************************************************/


//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// Function to get root file directly *or* get root file names from a text file
#include <boost/algorithm/string/predicate.hpp>
std::vector<std::string> GetFileList(std::string input_file)
{
  std::vector<std::string> filenames;
  if(boost::algorithm::ends_with(input_file,".root")){
    filenames.emplace_back(input_file);
    std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
  }
  else{
    std::ifstream input(input_file);
    for(std::string line; std::getline(input,line);){
      filenames.emplace_back(line);
      std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
    }
  }
  return filenames;
}




// ---------------------- Main function ------------------------ //

int main(int argv, char** argc)
{
  // Get input file list from first argument
  // There should be no more arguments - if there are they will be ignored!
  std::string input_files(argc[1]);

  // Calibration file

  // For MC:
  TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_x_calibration_mcc_8_4_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_x_calibration_mcc_8_4_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_x_calibration_mcc_8_4_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_y_z_calibration_mcc_8_4_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_y_z_calibration_mcc_8_4_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/latest_y_z_calibration_mcc_8_4_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 1.0;
  double correction_t_plane1 = 1.0;
  double correction_t_plane2 = 1.0;

  unsigned int expected_run = 0;

  // Below are files for data (time-dependent so you need different calibration files depending on when the data was taken)

  // File: /pnfs/uboone/persistent/users/greenlee/devel/v06_26_01_11/reco/test_recal_extbnb/4352331_11/PhysicsRun-2016_3_25_11_31_14-0005588-00219_20160326T063056_ext_bnb_20160327T090909_merged_20171123T155120_reco1_20171123T155611_reco2_20171123T162803_merged_20180224T055449_cali.root
  // Run number 5588, on 3/25/2016 11:31:01 to 18:32:39. Calibration constants taken for this run
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_25_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_25_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_25_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.963905;
  double correction_t_plane1 = 1.00089;
  double correction_t_plane2 = 1.00376;

  unsigned int expected_run = 5588;*/

  // File: xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr/uboone/persistent/users/greenlee/devel/v06_26_01_11/reco/test_recal_bnb/4172999_18/PhysicsRun-2016_2_29_0_18_42-0005204-00014_20160301T163428_bnb_20160301T214949_merged_20171213T103225_reco1_20171213T103630_reco2_20171213T125018_merged_20180226T054515_cali.root
  // Run number 5762, on to . Calibration constants taken for this run
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_4_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_4_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_4_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 1.04434;
  double correction_t_plane1 = 1.0601;
  double correction_t_plane2 = 1.04943;

  unsigned int expected_run = 5762;*/

  // File: xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr/uboone/persistent/users/greenlee/devel/v06_26_01_11/reco/test_recal_bnb/4172988_14/PhysicsRun-2016_3_11_20_33_17-0005390-00063_20160326T104002_bnb_20160330T032048_merged_20171209T164500_reco1_20171209T165202_reco2_20171209T191845_merged_20180226T055144_cali.root
  // Run number 5924, on 4/15/2016 10:34:26 to 11:55:06. Calibration constants taken for this run
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_15_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_15_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_15_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.976321;
  double correction_t_plane1 = 0.987975;
  double correction_t_plane2 = 0.999762;

  unsigned int expected_run = 5924;*/

  // Test: look at run 5179, from 02/27/2016
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.976124;
  double correction_t_plane1 = 0.994605;
  double correction_t_plane2 = 1.0034;

  unsigned int expected_run = 5179;*/

  // File: /pnfs/uboone/scratch/users/tjyang/PhysicsRun-2016_2_27_6_25_1-0005179-00169_20160301T101208_bnb_20160301T142750_merged_20171213T113100_reco1_20171213T113559_reco2_20171213T151521_merged_20180305T165948_cali.root
  // Run number 5179, from 02/27/2016
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_2_27_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.976124;
  double correction_t_plane1 = 0.994605;
  double correction_t_plane2 = 1.0034;

  unsigned int expected_run = 5179;*/

  // Test from the same file: run 5424, on 3/14/2016
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_14_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_14_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_3_14_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.973103;
  double correction_t_plane1 = 0.992772;
  double correction_t_plane2 = 1.00424;

  unsigned int expected_run = 5424;*/


  // File /pnfs/uboone/persistent/users/greenlee/devel/v06_26_01_11/reco/test_recal_bnb3/4362454_99/PhysicsRun-2016_4_6_20_41_21-0005804-00107_20160407T062331_bnb_20160407T065730_merged_20171207T172701_reco1_20171207T172901_reco2_20171207T213112_merged_20180305T172845_cali.root
  // Run 5804, from 04/06/2016
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_6_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_6_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_6_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 1.05005;
  double correction_t_plane1 = 1.04505;
  double correction_t_plane2 = 1.03369;

  unsigned int expected_run = 5804;*/

  // File /pnfs/uboone/persistent/users/greenlee/devel/v06_26_01_11/reco/test_recal_bnb3/4685876_70/PhysicsRun-2016_4_14_13_15_24-0005914-00136_20160416T125120_bnb_20160503T035406_merged_20171206T224811_reco1_20171206T230458_reco2_20171207T010738_merged_20180305T173802_cali.root
  // Run 5914, from 04/14/2016
  /*TFile *fcalib_x_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_14_plane_0.root","read");
  TH1F *correction_x_plane0 = (TH1F*)fcalib_x_plane0->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_14_plane_1.root","read");
  TH1F *correction_x_plane1 = (TH1F*)fcalib_x_plane1->Get("dq_dx_x_error_hist");
  TFile *fcalib_x_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/X_correction_factors_2016_4_14_plane_2.root","read");
  TH1F *correction_x_plane2 = (TH1F*)fcalib_x_plane2->Get("dq_dx_x_error_hist");

  TFile *fcalib_yz_plane0 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_0.root","read");
  TH2F *correction_yz_plane0 = (TH2F*)fcalib_yz_plane0->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane1 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_1.root","read");
  TH2F *correction_yz_plane1 = (TH2F*)fcalib_yz_plane1->Get("error_dq_dx_z_vs_y_hist");
  TFile *fcalib_yz_plane2 = new TFile("/pnfs/uboone/persistent/users/vmeddage/final_calibration_root_files/YZ_correction_factors_2016_month_2_month_3_month_4_month_5_plane_2.root","read");
  TH2F *correction_yz_plane2 = (TH2F*)fcalib_yz_plane2->Get("error_dq_dx_z_vs_y_hist");

  double correction_t_plane0 = 0.984039;
  double correction_t_plane1 = 0.983625;
  double correction_t_plane2 = 1.00004;

  unsigned int expected_run = 5914;*/

  gStyle->SetOptStat(0);

  int n_total = 0;
  int n_correct = 0;

  // Format files list
  std::vector<std::string> filenames = GetFileList(input_files);

  // Loop through events
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //gallery::Event ev(filenames);

  std::cout << "Event timestamp: " << ev.eventAuxiliary().time().value() << std::endl;
  std::cout << "Run number: " << ev.eventAuxiliary().run() << std::endl;
  std::cout << "     Subrun number: " << ev.eventAuxiliary().subRun() << std::endl
	    << "     Event number:  " << ev.eventAuxiliary().event() << std::endl;

  if (expected_run != 0 && ev.eventAuxiliary().run() != expected_run){
    std::cout << "[ERROR] run number " << ev.eventAuxiliary().run() << " not " << expected_run << " as expected. Skipping event because calibration constants will be wrong!" << std::endl;
    continue;
    }

    // Get pandoraNu tracks
    auto const trk_handle = ev.getValidHandle<std::vector<recob::Track>>("pandoraNu");
    std::vector<recob::Track> trackVec(*trk_handle);
    art::FindManyP<anab::Calorimetry> uncalib_calos_from_tracks(trk_handle, ev, art::InputTag("pandoraNucalo"));
    art::FindManyP<anab::Calorimetry> calib_calos_from_tracks(trk_handle, ev, art::InputTag("pandoraNucali"));

    for (int i_track=0; i_track<(int)trackVec.size(); i_track++){
      recob::Track track = trk_handle->at(i_track);
      unsigned int trkid = track.ID();
      std::vector<art::Ptr<anab::Calorimetry>> uncalib_calos = uncalib_calos_from_tracks.at(trkid);
      std::vector<art::Ptr<anab::Calorimetry>> calib_calos = calib_calos_from_tracks.at(trkid);

      // Loop through the three planes
      for (int plane=0; plane < 3; plane++){

        std::cout << "------ Looking at calibration for plane " << plane << " -------- " << std::endl;

	std::vector<double> uncalib_dqdx_v;
	std::vector<TVector3> uncalib_xyz_v;
	std::vector<double> calib_dqdx_v;
	std::vector<TVector3> calib_xyz_v;

	// Loop through calorimetry objects and get dqdx and xyz: uncalib
	for (auto c : uncalib_calos) {
	  if (!c) continue; // avoid art errors if c doesn't exist
	  if (!c->PlaneID().isValid) continue; // check has valid plane
	  int planenum = c->PlaneID().Plane;
	  if (planenum != plane) continue; // only use calorimetry from plane 2

	  uncalib_dqdx_v = c->dQdx();
	  uncalib_xyz_v  = c->XYZ();
	} // close loop over uncalib calos

	// Loop through calorimetry objects and get dqdx and xyz: calib
	for (auto c : calib_calos) {
	  if (!c) continue; // avoid art errors if c doesn't exist
	  if (!c->PlaneID().isValid) continue; // check has valid plane
	  int planenum = c->PlaneID().Plane;
	  if (planenum != plane) continue; // only use calorimetry from plane 2

	  calib_dqdx_v = c->dQdx();
	  calib_xyz_v  = c->XYZ();
	} // close loop over calib calos


	for (unsigned int i_v=0; i_v < calib_dqdx_v.size(); i_v++){

	  double truecorr;
	  double corryz, corrx, corrt;

	  if (plane==0){
	    int binyz = correction_yz_plane0->FindBin(calib_xyz_v.at(i_v).Z(),calib_xyz_v.at(i_v).Y());
	    int binx = correction_x_plane0->FindBin(calib_xyz_v.at(i_v).X());
	    corryz = correction_yz_plane0->GetBinContent(binyz);
	    corrx = correction_x_plane0->GetBinContent(binx);
	    corrt = correction_t_plane0;

	    /*std::cout << "Correction_yz_plane0 = " << correction_yz_plane0->GetBinContent(binyz) << std::endl
		      << "Correction_x_plane0 = " << correction_x_plane0->GetBinContent(binx) << std::endl
		      << "Correction_t_plane0 = " << correction_t_plane0 << std::endl;*/
	  }
	  else if (plane==1){
	    int binyz = correction_yz_plane1->FindBin(calib_xyz_v.at(i_v).Z(),calib_xyz_v.at(i_v).Y());
	    int binx = correction_x_plane1->FindBin(calib_xyz_v.at(i_v).X());
	    corryz = correction_yz_plane1->GetBinContent(binyz);
	    corrx = correction_x_plane1->GetBinContent(binx);
	    corrt = correction_t_plane1;

	    /*std::cout << "binyz = " << binyz << ", Correction_yz_plane1 = " << correction_yz_plane1->GetBinContent(binyz) << std::endl
		      << "binx = " << binx << ", Correction_x_plane1 = " << correction_x_plane1->GetBinContent(binx) << std::endl
		      << "Correction_t_plane1 = " << correction_t_plane1 << std::endl;*/
	  }
	  else if (plane==2){
	    int binyz = correction_yz_plane2->FindBin(calib_xyz_v.at(i_v).Z(),calib_xyz_v.at(i_v).Y());
	    int binx = correction_x_plane2->FindBin(calib_xyz_v.at(i_v).X());
	    corryz = correction_yz_plane2->GetBinContent(binyz);
	    corrx = correction_x_plane2->GetBinContent(binx);
	    corrt = correction_t_plane2;

	    /*std::cout << "Correction_yz_plane2 = " << correction_yz_plane2->GetBinContent(binyz) << std::endl
		      << "Correction_x_plane2 = " << correction_x_plane2->GetBinContent(binx) << std::endl
		      << "Correction_t_plane2 = " << correction_t_plane2 << std::endl;*/
	  }

	  if (corryz == 0.) //{std::cout << "corrYZ=0, setting to 1" << std::endl; corryz = 1.;}
	    corryz = 1.;
	  if (corrx == 0.) //{std::cout << "corrX=0, setting to 1" << std::endl; corrx = 1.;}
	    corrx = 1.;

	  truecorr = corryz*corrx*corrt;

	  n_total++;
	  if (TMath::Abs(calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v) - truecorr) < 1e-4){ n_correct++; }

	  if (TMath::Abs(calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v) - truecorr) > 1e-4){

	    std::cout << "--- Entry " << i_v << " in calib vector (plane " << plane << ") ---" << std::endl;

	    std::cout << "X = " << calib_xyz_v.at(i_v).X() << ", Y = " << calib_xyz_v.at(i_v).Y() << ", Z = " << calib_xyz_v.at(i_v).Z() << std::endl
		      << "    calib dqdx = " << calib_dqdx_v.at(i_v) << std::endl;
	    std::cout << "uncalib X = " << uncalib_xyz_v.at(i_v).X() << ", uncalib Y = " << uncalib_xyz_v.at(i_v).Y() << ", uncalib Z = " << uncalib_xyz_v.at(i_v).Z() << std::endl
		      << "    uncalib dqdx = " << uncalib_dqdx_v.at(i_v) << std::endl;
	    std::cout << "Ratio calib/uncalib dqdx = " << calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v)
		      << ", should be " << truecorr << std::endl;
	    std::cout << "Difference is " << TMath::Abs(calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v) - truecorr) << std::endl;

	    std::cout << "Ratio calib/uncalib dqdx divided by true X/Y/Z corr = " << (calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v))/(corrx*corryz) << std::endl;
	    //std::cout << bool(TMath::Abs(calib_dqdx_v.at(i_v)/uncalib_dqdx_v.at(i_v) - truecorr) < 1e-4) << std::endl;
	  }

	}
      } // close loop over planes

    } // close loop over tracks

  } // loop through gallery events

    std::cout << std::endl;
    std::cout << "Looked at " << n_total << " hits, found " << n_correct << " correctly calibrated" << std::endl;

  return 0;
}
