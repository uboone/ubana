/*************************************************************
 * 
 * This is a quick macro to compare truth matching methods
 * (for now, on pandoraNu tracks only)
 *
 * Don't forget to set up gallery first:
 * setup gallery v1_03_08 -q e10:nu:prof
 *
 * Then compile:
 * make clean; make checktruth
 *
 * ./checktruth "path/to/reco2/file.root or path/to/list/of/input/files.txt"
 *
 * Kirsty Duffy (kduffy@fnal.gov), Fermilab, March 27 2018
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
#include "canvas/Persistency/Common/FindMany.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// Local includes: algorithm to find well-reconstructed tracks and return PDG code
#include "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/ParticleID/Algorithms/WellReconstructedTrackFinder.h"
//#include "BackTrackerTruthMatch.cxx"

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

  // Format files list
  std::vector<std::string> filenames = GetFileList(input_files);

  // Output file
  TFile *fout = new TFile("out_CheckTruthMatching.root","recreate");

  // Make histograms to store results
  TH1F *h_AssnsGaus_alltracks = new TH1F("h_AssnsGaus_alltracks", "Truth PDGs using associations (gaushit, all tracks);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_AssnsPCHR_alltracks = new TH1F("h_AssnsPCHR_alltracks", "Truth PDGs using associations (pandoraCosmicHitRemoval, all tracks);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_AssnsGaus_WRtracks = new TH1F("h_AssnsGaus_WRtracks", "Truth PDGs using associations (gaushit, well-reconstructed tracks only);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_AssnsPCHR_WRtracks = new TH1F("h_AssnsPCHR_WRtracks", "Truth PDGs using associations (pandoraCosmicHitRemoval, well-reconstructed tracks only);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_Afro_alltracks = new TH1F("h_Afro_alltracks", "Truth PDGs using Afro's backtracker (all tracks);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_Afro_WRtracks = new TH1F("h_Afro_WRtracks", "Truth PDGs using Afro's backtracker (well-reconstructed tracks only);Abs(PDG code);No. reco tracks",2500,0,2500);
  TH1F *h_WellRecod = new TH1F("h_WellRecod", "Truth PDGs using track start/end point (well-reconstructed tracks only);Abs(PDG code);No. reco tracks",2500,0,2500);

  int n_same = 0;
  int n_different = 0;
  int n_WR = 0;
  int n_notWR = 0;

  int n_Gaus_mu_all = 0;
  int n_Gaus_p_all  = 0;
  int n_Gaus_pi_all = 0;
  int n_Gaus_e_all  = 0;
  int n_Gaus_oth_all = 0;
  int n_Gaus_mu_WR = 0;
  int n_Gaus_p_WR  = 0;
  int n_Gaus_pi_WR = 0;
  int n_Gaus_e_WR  = 0;
  int n_Gaus_oth_WR = 0;

  int n_pCHR_mu_all = 0;
  int n_pCHR_p_all  = 0;
  int n_pCHR_pi_all = 0;
  int n_pCHR_e_all  = 0;
  int n_pCHR_oth_all = 0;
  int n_pCHR_mu_WR = 0;
  int n_pCHR_p_WR  = 0;
  int n_pCHR_pi_WR = 0;
  int n_pCHR_oth_WR = 0;
  int n_pCHR_e_WR  = 0;

  int n_Afro_mu_all = 0;
  int n_Afro_p_all  = 0;
  int n_Afro_pi_all = 0;
  int n_Afro_e_all  = 0;
  int n_Afro_oth_all = 0;
  int n_Afro_mu_WR = 0;
  int n_Afro_p_WR  = 0;
  int n_Afro_pi_WR = 0;
  int n_Afro_e_WR  = 0;
  int n_Afro_oth_WR = 0;

  int n_WR_mu = 0;
  int n_WR_p  = 0;
  int n_WR_pi = 0;
  int n_WR_e  = 0;
  int n_WR_oth = 0;
  
  // Loop through events
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
  //gallery::Event ev(filenames);{
  
  std::cout << "Run number: " << ev.eventAuxiliary().run() << std::endl;
  std::cout << "     Subrun number: " << ev.eventAuxiliary().subRun() << std::endl
	    << "     Event number:  " << ev.eventAuxiliary().event() << std::endl;
  
    // Get pandoraNu tracks
    auto const trk_handle = ev.getValidHandle<std::vector<recob::Track>>("pandoraNu::McRecoStage2");
    std::vector<recob::Track> trackVec(*trk_handle);
    art::FindManyP<recob::Hit> hits_per_tracks(trk_handle, ev, "pandoraNu::McRecoStage2");

   auto const hit_handle_gaus = ev.getValidHandle<std::vector<recob::Hit>>("gaushit::McRecoStage1");
    auto const hit_handle_pCHR = ev.getValidHandle<std::vector<recob::Hit>>("pandoraCosmicHitRemoval::McRecoStage2");

    // Get MCParticles
    auto const mcp_handle = ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");
    std::vector<simb::MCParticle> mcpVec(*mcp_handle);

    // Loop through tracks
    for (int i_track=0; i_track<(int)trackVec.size(); i_track++){
	
      std::cout << std::endl;
      std::cout << " -------------------------------------------------------------- " << std::endl;
      std::cout << std::endl;
      
      recob::Track track = trk_handle->at(i_track);
      unsigned int trkid = track.ID();

      std::vector<art::Ptr<recob::Hit>> hits_per_track = hits_per_tracks.at(trkid);

      
      int Assns_gaushit_PDG = -99999;
      int Assns_pandoraCHitR_PDG = -99999;
      int Afro_PDG = -99999;
      int WR_PDG = -99999;

      bool is_WR = false;
      
      // ----------------- Wes's truth matching: gaushit ----------------- //
      
      std::unordered_map<int,double> trkide_gaus;
      double maxe_gaus=-1, tote_gaus=0;
      simb::MCParticle const* maxp_me_gaus = NULL; //pointer for the particle match we will calculate

       art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit_gaus(hit_handle_gaus,ev,"gaushitTruthMatch::McRecoStage1");
      std::vector<simb::MCParticle const*> particle_vec_gaus;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec_gaus;
      
      //loop only over our hits
      for(size_t i_h=0; i_h<hits_per_track.size(); i_h++){
	
	particle_vec_gaus.clear(); match_vec_gaus.clear();
	particles_per_hit_gaus.get(hits_per_track[i_h].key(),particle_vec_gaus,match_vec_gaus);
	//the .key() gives us the index in the original collection

	//loop over particles
	for(size_t i_p=0; i_p<particle_vec_gaus.size(); ++i_p){
	  trkide_gaus[ particle_vec_gaus[i_p]->TrackId() ] += match_vec_gaus[i_p]->energy; //store energy per track id
	  tote_gaus += match_vec_gaus[i_p]->energy; //calculate total energy deposited
	  if( trkide_gaus[ particle_vec_gaus[i_p]->TrackId() ] > maxe_gaus ){ //keep track of maximum
	    maxe_gaus = trkide_gaus[ particle_vec_gaus[i_p]->TrackId() ];
	    maxp_me_gaus = particle_vec_gaus[i_p];
	  }
	}//end loop over particles per hit

      }

      //for (std::pair<int,double> element : trkide_gaus){
      // std::cout << "trkide_gaus: " << element.first << "::" << element.second << std::endl;
      //}

      Assns_gaushit_PDG = maxp_me_gaus->PdgCode();
      
      // ----------------- Try only considering number of shared hits ----------------- //

      // Loop over all MCParticles and look at all hits associated with them
      // If the hits associated to an MCParticle also appear in the collection of hits associated
      // to our track, count them up to find the MCP that has most hits in common with the track
      
      /*art::FindMany<recob::Hit, anab::BackTrackerHitMatchingData> hits_from_mcps(mcp_handle, ev, "gaushitTruthMatch::McRecoStage1");
      std::vector<recob::Hit const*> hits_from_mcp;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

      std::unordered_map<int,int> mcptrkid_nhits;
      std::unordered_map<int,int> mcptrkid_pdg;
      std::unordered_map<int,double> mcptrkid_energy;

      std::cout << "Total hits in reco track: " << hits_per_track.size() << std::endl;
      double total_energy = 0.0;
      
      for (size_t i_mcp=0; i_mcp < mcpVec.size(); i_mcp++){
	int mcptrackID = mcp_handle->at(i_mcp).TrackId();

	mcptrkid_pdg[ mcptrackID ] = mcp_handle->at(i_mcp).PdgCode();

	hits_from_mcp.clear();
	match_vec.clear();
	hits_from_mcps.get(i_mcp,hits_from_mcp,match_vec);
	
	//std::vector<art::Ptr<recob::Hit>> hits_from_mcp = hits_from_mcps.at(i_mcp);
	for (size_t i_h=0; i_h < hits_from_mcp.size(); i_h++){
	  total_energy += match_vec.at(i_h)->energy;
	  for (size_t i_htrack=0; i_htrack < hits_per_track.size(); i_htrack++){
	    if (hits_per_track[i_htrack].key() == i_h){
	      mcptrkid_nhits[ mcptrackID ]++;
	      mcptrkid_energy[ mcptrackID ] += match_vec.at(i_h)->energy;
	    }
	  }
	}
      }

      std::cout << "Total energy: " << total_energy << std::endl;
      
      for (std::pair<int,int> element : mcptrkid_nhits){
	std::cout << "mcptrkid_nhits: " << element.first << " :: " << element.second << ", PDG " << mcptrkid_pdg[ element.first ] << ", energy = " << mcptrkid_energy[ element.first ] << std::endl;
      }

      // Try the other way round (but rewriting it myself)
       art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit_k(hit_handle_gaus,ev,"gaushitTruthMatch::McRecoStage1");
      std::vector<simb::MCParticle const*> particle_vec_k;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec_k;

      std::unordered_map<int,double> mcptrkid_energy_k;
      std::unordered_map<int,int> mcptrkid_pdg_k;
      std::unordered_map<int,int> mcptrkid_nhits_k;
      double total_energy_k = 0.0;

      for (size_t i_h=0; i_h < hits_per_track.size(); i_h++){
	particle_vec_k.clear();
	match_vec_k.clear();
	particles_per_hit_k.get(hits_per_track.at(i_h).key(),particle_vec_k,match_vec_k);

	for (size_t i_mcp=0; i_mcp < particle_vec_k.size(); i_mcp++){
	  int mcptrackID = particle_vec_k.at(i_mcp)->TrackId();
	  if (match_vec_k.at(i_mcp)->energy != 0.0){
	    mcptrkid_pdg_k[ mcptrackID ] = particle_vec_k.at(i_mcp)->PdgCode();
	    mcptrkid_energy_k[ mcptrackID ] += match_vec_k.at(i_mcp)->energy;
	    mcptrkid_nhits_k[ mcptrackID ] ++;
	    total_energy_k += match_vec_k.at(i_mcp)->energy;
	  }
	}
      }

      std::cout << "Total energy (MCPs from hits): " << total_energy << std::endl;
      
      for (std::pair<int,double> element : mcptrkid_energy_k){
	//if (
	std::cout << "mcptrkid_energy_k: " << element.first << " :: " << element.second << ", PDG " << mcptrkid_pdg_k[ element.first ] << ", nhits: " << mcptrkid_nhits_k[ element.first ]  << std::endl;
	}*/

      // ----------------- Wes's truth matching: pandoraCosmicHitRemoval ----------------- //
      
      std::unordered_map<int,double> trkide_pCHR;
      double maxe_pCHR=-1, tote_pCHR=0;
      simb::MCParticle const* maxp_me_pCHR = NULL; //pointer for the particle match we will calculate

       art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit_pCHR(hit_handle_pCHR,ev,"crHitRemovalTruthMatch::McRecoStage2");
      std::vector<simb::MCParticle const*> particle_vec_pCHR;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec_pCHR;
      
      //loop only over our hits
      for(size_t i_h=0; i_h<hits_per_track.size(); ++i_h){
	
	particle_vec_pCHR.clear(); match_vec_pCHR.clear();
	particles_per_hit_pCHR.get(hits_per_track[i_h].key(),particle_vec_pCHR,match_vec_pCHR);
	//the .key() gives us the index in the original collection

	//loop over particles
	for(size_t i_p=0; i_p<particle_vec_pCHR.size(); ++i_p){
	  trkide_pCHR[ particle_vec_pCHR[i_p]->TrackId() ] += match_vec_pCHR[i_p]->energy; //store energy per track id
	  tote_pCHR += match_vec_pCHR[i_p]->energy; //calculate total energy deposited
	  if( trkide_pCHR[ particle_vec_pCHR[i_p]->TrackId() ] > maxe_pCHR ){ //keep track of maximum
	    maxe_pCHR = trkide_pCHR[ particle_vec_pCHR[i_p]->TrackId() ];
	    maxp_me_pCHR = particle_vec_pCHR[i_p];
	  }
	}//end loop over particles per hit

      }

      std::cout << std::endl;
      std::cout << "Final Match (Assns: pandoraCosmicHitRemoval) is pdg = " << maxp_me_pCHR->PdgCode() << " with energy " << maxe_pCHR << " over " << tote_pCHR << " (" << maxe_pCHR/tote_pCHR << ")"
		<< " trkid=" << maxp_me_pCHR->TrackId()
		<< " ke=" << maxp_me_pCHR->E()-maxp_me_pCHR->Mass()
		<< "\n\tstart (x,y,z)=(" << maxp_me_pCHR->Vx()
		<< "," << maxp_me_pCHR->Vy()
		<< "," << maxp_me_pCHR->Vz()
		<< ")\tend (x,y,z)=(" << maxp_me_pCHR->EndX()
		<< "," << maxp_me_pCHR->EndY()
		<< "," << maxp_me_pCHR->EndZ() << ")" << std::endl;

      Assns_pandoraCHitR_PDG = maxp_me_pCHR->PdgCode();
      
      // ----------------- Afro's method using new backtracker ----------------- //

      /*BackTrackerTruthMatch backtrackertruthmatch;
      backtrackertruthmatch.MatchToMCParticle(hit_handle_pCHR,ev,hits_per_track);
      art::Ptr< simb::MCParticle > MCP_Afro = backtrackertruthmatch.ReturnMCParticle();
      
      Afro_PDG = MCP_Afro->PdgCode();

      std::cout << std::endl;
      std::cout << "Final Match (Afro: pandoraCosmicHitRemoval) is pdg = " << Afro_PDG << std::endl
		<< "  MCTrackID = " << MCP_Afro->TrackId() << std::endl
		<< "  with purity " << backtrackertruthmatch.ReturnPurity() << std::endl;*/

      // ----------------- "Well reconstructed tracks" method ----------------- //

      int WR_trackID = -99999;
      simb::MCParticle WR_MCP;
      
      for(size_t i_mcp=0; i_mcp<mcpVec.size(); i_mcp++){
	simb::MCParticle thisMCP = mcpVec.at(i_mcp);
	if (thisMCP.StatusCode() != 1) continue;
	
	if (isWellReconstructed(track, thisMCP)){
	  is_WR = true;
	  WR_PDG = thisMCP.PdgCode();
	  WR_trackID = thisMCP.TrackId();
	  WR_MCP = thisMCP;
	  break;
	}
      }// end loop over MCParticles

      
      // ------------------------------------------------------------ //      
      // ----------------- Finally, fill histograms ----------------- //
      // ------------------------------------------------------------ //

      h_AssnsGaus_alltracks->Fill(TMath::Abs(Assns_gaushit_PDG));
      h_AssnsPCHR_alltracks->Fill(TMath::Abs(Assns_pandoraCHitR_PDG));
      h_Afro_alltracks->Fill(TMath::Abs(Afro_PDG));

      if (is_WR){
	h_AssnsGaus_WRtracks->Fill(TMath::Abs(Assns_gaushit_PDG));
	h_AssnsPCHR_WRtracks->Fill(TMath::Abs(Assns_pandoraCHitR_PDG));
	h_Afro_WRtracks->Fill(TMath::Abs(Afro_PDG));
	h_WellRecod->Fill(TMath::Abs(WR_PDG));

	n_WR++;

	if (TMath::Abs(Assns_gaushit_PDG) == 13) n_Gaus_mu_WR++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 2212) n_Gaus_p_WR++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 211) n_Gaus_pi_WR++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 11) n_Gaus_e_WR++;
	else n_Gaus_oth_WR++;
	
	if (TMath::Abs(Assns_pandoraCHitR_PDG) == 13) n_pCHR_mu_WR++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 2212) n_pCHR_p_WR++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 211) n_pCHR_pi_WR++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 11) n_pCHR_e_WR++;
	else n_pCHR_oth_WR++;
	
	if (TMath::Abs(Afro_PDG) == 13) n_Afro_mu_WR++;
	else if (TMath::Abs(Afro_PDG) == 2212) n_Afro_p_WR++;
	else if (TMath::Abs(Afro_PDG) == 211) n_Afro_pi_WR++;
	else if (TMath::Abs(Afro_PDG) == 11) n_Afro_e_WR++;
	else n_Afro_oth_WR++;
	
	if (TMath::Abs(WR_PDG) == 13) n_WR_mu++;
	else if (TMath::Abs(WR_PDG) == 2212) n_WR_p++;
	else if (TMath::Abs(WR_PDG) == 211) n_WR_pi++;
	else if (TMath::Abs(WR_PDG) == 11) n_WR_e++;
	else n_WR_oth++;
	
	if (Assns_pandoraCHitR_PDG == WR_PDG) n_same++;
	else {
	  n_different++;

	  std::cout << "Track match from associations != track match by track start/end." << std::endl;
	  std::cout << "Track start, end: " << "(" << track.Start().X() << ", " << track.Start().Y() << ", " << track.Start().Z() << "), (" << track.End().X() << ", " << track.End().Y() << ", " << track.End().Z() << ")" << std::endl;
	  std::cout << std::endl;
	  std::cout << "Final Match (Assns: pCHR) is pdg = " << maxp_me_pCHR->PdgCode() << " with energy " << maxe_pCHR << " over " << tote_pCHR << " (" << maxe_pCHR/tote_pCHR << ")"
		    << "\n\ttrkid=" << maxp_me_pCHR->TrackId()
		    << " ke=" << maxp_me_pCHR->E()-maxp_me_pCHR->Mass()
		    << "\n\tstart (x,y,z)=(" << maxp_me_pCHR->Vx()
		    << "," << maxp_me_pCHR->Vy()
		    << "," << maxp_me_pCHR->Vz()
		    << ")\tend (x,y,z)=(" << maxp_me_pCHR->EndX()
		    << "," << maxp_me_pCHR->EndY()
		    << "," << maxp_me_pCHR->EndZ() << ")" << std::endl;

	  std::cout << std::endl;
	  std::cout << "Final Match (Well Reconstructed track) is " << WR_trackID
		    << "\n\tpdg=" << WR_PDG 
		    << "\n\tstart(x,y,z) = (" << WR_MCP.Vx() << "," << WR_MCP.Vy() << "," << WR_MCP.Vz() << ")"
		    << "\n\tend(x,y,z) = (" << WR_MCP.EndX() << "," << WR_MCP.EndY() << "," << WR_MCP.EndZ() << ")" << std::endl;
		    }
      }
      else{
	n_notWR++;

	if (TMath::Abs(Assns_gaushit_PDG) == 13) n_Gaus_mu_all++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 2212) n_Gaus_p_all++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 211) n_Gaus_pi_all++;
	else if (TMath::Abs(Assns_gaushit_PDG) == 11) n_Gaus_e_all++;
	else n_Gaus_oth_all++;
	
	if (TMath::Abs(Assns_pandoraCHitR_PDG) == 13) n_pCHR_mu_all++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 2212) n_pCHR_p_all++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 211) n_pCHR_pi_all++;
	else if (TMath::Abs(Assns_pandoraCHitR_PDG) == 21) n_pCHR_e_all++;
	else n_pCHR_oth_all++;
	
	if (TMath::Abs(Afro_PDG) == 13) n_Afro_mu_all++;
	else if (TMath::Abs(Afro_PDG) == 2212) n_Afro_p_all++;
	else if (TMath::Abs(Afro_PDG) == 211) n_Afro_pi_all++;
	else if (TMath::Abs(Afro_PDG) == 11) n_Afro_e_all++;
	else n_Afro_oth_all++;
      }
	
      std::cout << std::endl;
      std::cout << " -------------------------------------------------------------- " << std::endl;
      std::cout << std::endl;

    } // close loop over tracks
    
  } // loop through gallery events

  fout->cd();
  
  h_AssnsGaus_alltracks->Write();
  h_AssnsPCHR_alltracks->Write();
  h_Afro_alltracks->Write();

  h_AssnsGaus_WRtracks->Write();
  h_AssnsPCHR_WRtracks->Write();
  h_Afro_WRtracks->Write();
  h_WellRecod->Write();
  
    std::cout << std::endl;
    std::cout << "Saw " << n_notWR << " not well reconstructed tracks." << std::endl;
    std::cout << "Saw " << n_WR << " well reconstructed tracks, of which: " << std::endl;
    
    /* std::cout << "   gaushittruthmatch: " << n_Gaus_mu_WR << " muons" << std::endl
	      << "                      " << n_Gaus_p_WR << " protons" << std::endl
	      << "                      " << n_Gaus_pi_WR << " pions" << std::endl
	      << "                      " << n_Gaus_e_WR << " electrons" << std::endl
	      << "                      " << n_Gaus_oth_WR << " other" << std::endl;*/
    std::cout << "   pCHRtruthmatch:    " << n_pCHR_mu_WR << " muons" << std::endl
	      << "                      " << n_pCHR_p_WR << " protons" << std::endl
	      << "                      " << n_pCHR_pi_WR << " pions" << std::endl
	      << "                      " << n_pCHR_e_WR << " electrons" << std::endl
	      << "                      " << n_pCHR_oth_WR << " other" << std::endl;
    /*std::cout << "   Afro's truthmatch: " << n_Afro_mu_WR << " muons" << std::endl
	      << "                      " << n_Afro_p_WR << " protons" << std::endl
	      << "                      " << n_Afro_pi_WR << " pions" << std::endl
	      << "                      " << n_Afro_e_WR << " electrons" << std::endl
	      << "                      " << n_Afro_oth_WR << " other" << std::endl;*/
    std::cout << "   WR truthmatch:     " << n_WR_mu << " muons" << std::endl
	      << "                      " << n_WR_p << " protons" << std::endl
	      << "                      " << n_WR_pi << " pions" << std::endl
	      << "                      " << n_WR_e << " electrons" << std::endl
	      << "                      " << n_WR_oth << " other" << std::endl;

    
    std::cout << "Saw " << n_notWR+n_WR << "total tracks, of which: " << std::endl;
    std::cout << "   gaushittruthmatch: " << n_Gaus_mu_all << " muons" << std::endl
	      << "                      " << n_Gaus_p_all << " protons" << std::endl
	      << "                      " << n_Gaus_pi_all << " pions" << std::endl
	      << "                      " << n_Gaus_e_all << " electrons" << std::endl
	      << "                      " << n_Gaus_oth_all << " other" << std::endl;
    std::cout << "   pCHRtruthmatch:    " << n_pCHR_mu_all << " muons" << std::endl
	      << "                      " << n_pCHR_p_all << " protons" << std::endl
	      << "                      " << n_pCHR_pi_all << " pions" << std::endl
	      << "                      " << n_pCHR_e_all << " electrons" << std::endl
	      << "                      " << n_pCHR_oth_all << " other" << std::endl;
    std::cout << "   Afro's truthmatch: " << n_Afro_mu_all << " muons" << std::endl
	      << "                      " << n_Afro_p_all << " protons" << std::endl
	      << "                      " << n_Afro_pi_all << " pions" << std::endl
	      << "                      " << n_Afro_e_all << " electrons" << std::endl
	      << "                      " << n_Afro_oth_all << " other" << std::endl;

  return 0;
}
