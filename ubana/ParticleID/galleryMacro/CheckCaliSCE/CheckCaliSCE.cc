/*************************************************************
 *
 * This is a quick macro for counting how many daughter tracks a
 * reconstructed track has, using the number of tracks with start
 * or end points within a given distance of the end point of the
 * first track.
 *
 * Don't forget to set up gallery first:
 * setup gallery v1_11_02 -q e17:prof
 *
 * Then run with:
 * root
 * root[0] .L CheckCaliSCE.cc
 * root[1] CheckCaliSCE()
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
void CheckCaliSCE()
{
  // Get input file list from first argument
  // There should be no more arguments - if there are they will be ignored!
  // std::string input_files = "/pnfs/uboone/mc/uboone/reconstructed/prod_v08_00_00_07/prodgenie_bnb_nu_corsika_noSCE_genie2_new/noSCE_reco2/prodgenie_bnb_nu_cosmic_uboone_20190308T061338_gen_20190308T075946_g4_nospacecharge_20190309T010439_detsim_20190309T044244_reco1_postdlmctruth_postdlm_f8a3f245-cf0e-4fda-ab88-bee681838041.root"; // MC no SCE

  // std::string input_files = "/pnfs/uboone/mc/uboone/reconstructed/prod_v08_00_00_07/prodgenie_bnb_nu_corsika_SCE_genie2_newy/SCE_reco2/prodgenie_bnb_nu_cosmic_uboone_20190305T020333_gen_20190305T022737_g4_20190305T062922_detsim_20190305T075301_reco1_postdlmctruth_postdlmc_postwcct_f5b_d4c57bd1-75e9-4ad0-8597-4833f1684b19.root"; // MC with SCE

  // std::string input_files = "/pnfs/uboone/data/uboone/reconstructed/prod_v08_00_00_07/data_bnb_optfilter_C1_5e19/reco2/00/00/51/60/PhysicsRun-2016_2_26_7_58_47-0005160-00014_20160303T140250_bnb_20160304T103227_merged_20181106T062926_optfilter_20181221T145715_reco1_postwcct_postdl_20181221T180347_20190304T004958_reco2.root"; // BNB data

  // std::string input_files = "/pnfs/uboone/data/uboone/reconstructed/prod_v08_00_00_07/data_extbnb_mcc9.0_validate/run1_reco2/00/00/65/53/PhysicsRun-2016_6_12_2_6_32-0006553-00101_20160612T123323_ext_bnb_20160612T150236_merged_20181106T233733_optfilter_20190103T124224_reco1_postwcct_postdl_20190104T120349_20190310T101224_reco2.root"; // EXT data

  // std::string input_files = "/pnfs/uboone/overlay/uboone/reconstructed/prod_v08_00_00_08/prodgenie_bnb_nu_uboone_overlay_mcc9/run1_reco2/PhysicsRun-2016_3_24_0_23_41-0005568-00027_20160326T033950_ext_unbiased_20160326T082530_merged_gen_20190309T190903_g4_detsim_mix_r1a_r1b_postdlmctruth_172acaa6-c307-4015-b6f1-fa0600ddac66.root"; // Run 1 overlay

  // std::string input_files = "/uboone/app/users/kduffy/CCincMCC9/srcs/ubana/ubana/ParticleID/galleryMacro/CheckCaliSCE/datafiles_joseph_20files_20190312.list"; // New data for testing

  // std::string input_files = "/uboone/app/users/kduffy/CCincMCC9/srcs/ubana/ubana/ParticleID/galleryMacro/CheckCaliSCE/datafiles_joseph_onbeam_20files_20190313.list"; // New data for testing

  std::string input_files = "/uboone/app/users/kduffy/CCincMCC9/srcs/ubana/ubana/ParticleID/galleryMacro/CheckCaliSCE/datafiles_joseph_offbeam_20files_20190314.list"; // New data for testing

  // std::string input_files = "/uboone/app/users/kduffy/CCincMCC9/srcs/ubana/ubana/ParticleID/galleryMacro/CheckCaliSCE/overlay_20files_wes_20190312.list"; // New overlay for testing

  gStyle->SetOptStat(0);

  int n_total = 0;
  int n_same = 0;
  int n_different = 0;
  int n_inf = 0;

  // Format files list
  std::vector<std::string> filenames = GetFileList(input_files);

  // Loop through events
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //gallery::Event ev(filenames);
  //
  // std::cout << "Run number: " << ev.eventAuxiliary().run() << std::endl;
  // std::cout << "     Subrun number: " << ev.eventAuxiliary().subRun() << std::endl
	//     << "     Event number:  " << ev.eventAuxiliary().event() << std::endl;


    // Get pandora tracks
    auto const trk_handle = ev.getValidHandle<std::vector<recob::Track>>("pandora");
    std::vector<art::Ptr<recob::Track>> trackVec;
    art::fill_ptr_vector(trackVec,trk_handle);
    art::FindManyP<anab::Calorimetry> calis_from_tracks(trk_handle, ev, art::InputTag("pandoracali"));
    art::FindManyP<anab::Calorimetry> caliSCEs_from_tracks(trk_handle, ev, art::InputTag("pandoracaliSCE"));

    for (auto& track : trackVec){
      // recob::Track track = trk_handle->at(i_track);
      std::vector<art::Ptr<anab::Calorimetry>> calis = calis_from_tracks.at(track.key());
      std::vector<art::Ptr<anab::Calorimetry>> caliSCEs = caliSCEs_from_tracks.at(track.key());

      // Loop through the three planes
      for (int plane=0; plane < 3; plane++){

        // std::cout << "------ Looking at plane " << plane << " -------- " << std::endl;

      	std::vector<float> cali_dEdx_v;
      	std::vector<float> caliSCE_dEdx_v;

	// Loop through calorimetry objects and get dEdx
	for (auto c : calis) {

	  if (!c) continue; // avoid art errors if c doesn't exist
	  if (!c->PlaneID().isValid) continue; // check has valid plane
	  int planenum = c->PlaneID().Plane;
	  if (planenum != plane) continue; // only use calorimetry from plane 2

	  cali_dEdx_v = c->dEdx();
	} // close loop over uncalib calos

	// Loop through calorimetry objects and get dqdx and xyz: calib
	for (auto c2 : caliSCEs) {
	  if (!c2) continue; // avoid art errors if c doesn't exist
	  if (!c2->PlaneID().isValid) continue; // check has valid plane
	  int planenum = c2->PlaneID().Plane;
	  if (planenum != plane) continue; // only use calorimetry from plane 2

	  caliSCE_dEdx_v = c2->dEdx();
	} // close loop over calib calos


	for (unsigned int i_v=0; i_v < cali_dEdx_v.size(); i_v++){

	  n_total++;

	  if (TMath::Abs(cali_dEdx_v.at(i_v) - caliSCE_dEdx_v.at(i_v)) < 1e-10){

	    // std::cout << "--- Entry " << i_v << " in calib vector (plane " << plane << ") ---" << std::endl;
      //
      // std::cout << "Cali: " << cali_dEdx_v.at(i_v) << ", CaliSCE: " << caliSCE_dEdx_v.at(i_v) << std::endl;

      n_same++;
	  }
    else{

	    // std::cout << "--- Entry " << i_v << " in calib vector (plane " << plane << ") ---" << std::endl;
      //
      // std::cout << "Cali: " << cali_dEdx_v.at(i_v) << ", CaliSCE: " << caliSCE_dEdx_v.at(i_v) << std::endl;
      n_different++;
    }

    if (!(std::isnormal(cali_dEdx_v.at(i_v)) && std::isnormal(caliSCE_dEdx_v.at(i_v)))){
      n_inf++;

      std::cout << "Run number: " << ev.eventAuxiliary().run() << std::endl;
      std::cout << "     Subrun number: " << ev.eventAuxiliary().subRun() << std::endl
    	    << "     Event number:  " << ev.eventAuxiliary().event() << std::endl;

	    std::cout << "--- Entry " << i_v << " in calib vector (plane " << plane << ") ---" << std::endl;

      std::cout << "Cali: " << cali_dEdx_v.at(i_v) << ", CaliSCE: " << caliSCE_dEdx_v.at(i_v) << std::endl;
    }

	}
      } // close loop over planes

    } // close loop over tracks

  } // loop through gallery events

    std::cout << std::endl;
    std::cout << "Looked at " << n_total << " hits, found " << n_same << " the same for cali and caliSCE, and " << n_different << " different between the two" << std::endl;
    std::cout << "No. inf = " << n_inf << std::endl;

  return 0;
}
