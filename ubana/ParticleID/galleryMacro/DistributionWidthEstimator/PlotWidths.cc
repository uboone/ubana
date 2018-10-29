/**************************************************
 * Macro to estimate uncertainty in dE/dx (for use in PID algorithms).
 * Idea is to look at all tracks in the range where the muon dE/dx distribution
 * should be flat (RR=100-150 cm), true MC muon tracks in that range, and true
 * MC proton tracks with residual range > 30 cm
 *
 * Kirsty Duffy, Fermilab, 1st March 2018
 **************************************************/


// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

// root includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"

// art includes
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// data object includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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

bool isInTPC(double x, double y, double z){

  double DET_X_HIGH = 256.0;
  double DET_X_LOW = 0.0;
  double DET_Y_HIGH = 116.5;
  double DET_Y_LOW = -116.5;
  double DET_Z_HIGH = 1040;
  double DET_Z_LOW = 0.0;

  if (x > DET_X_LOW && x < DET_X_HIGH &&
      y > DET_Y_LOW && y < DET_Y_HIGH &&
      z > DET_Z_LOW && z < DET_Z_HIGH)
    return true;
  else return false;

}

bool isWellReconstructed(recob::Track track, simb::MCParticle mcp){

  double tsy = track.Start().Y();
  double tey = track.End().Y();
  double tsz = track.Start().Z();
  double tez = track.End().Z();

  double msy = mcp.Vy();
  double mey = mcp.EndY();
  double msz = mcp.Vz();
  double mez = mcp.EndZ();

  double twoDStartRes = std::sqrt(std::pow(msy-tsy,2)+std::pow(msz-tsz,2));
  double twoDStartResFlip = std::sqrt(std::pow(msy-tey,2)+std::pow(msz-tez,2));
  double twoDEndRes = std::sqrt(std::pow(mey-tey,2)+std::pow(mez-tez,2));
  double twoDEndResFlip = std::sqrt(std::pow(mey-tsy,2)+std::pow(mez-tsz,2));

  if ((twoDStartRes < 2.0 && twoDEndRes < 2.0) ||
      (twoDStartResFlip < 2.0 && twoDEndResFlip < 2.0)){
    std::cout << "[MCP Matching] Found a match!" << std::endl;
    std::cout << "[MCP Matching] Start Res Fwd : " << twoDStartRes << std::endl;
    std::cout << "[MCP Matching] End Res Fwd   : " << twoDEndRes << std::endl;
    std::cout << "[MCP Matching] Start Res Bwd : " << twoDStartResFlip << std::endl;
    std::cout << "[MCP Matching] End Res Bwd   : " << twoDEndResFlip << std::endl;
    return true;
  }
  else {
    /*    std::cout << "[MCP Matching] Start Res Fwd : " << twoDStartRes << std::endl;
          std::cout << "[MCP Matching] End Res Fwd   : " << twoDEndRes << std::endl;
          std::cout << "[MCP Matching] Start Res Bwd : " << twoDStartResFlip << std::endl;
          std::cout << "[MCP Matching] End Res Bwd   : " << twoDEndResFlip << std::endl;
          */    return false;
  }
}

// ---------------------------------------------------------------------------- //

int main(int argv, char** argc){

  // extract input filename from input textfile
  std::string input_files(argc[1]);
  std::vector<std::string> filenames = GetFileList(input_files);

  TFile *fout = new TFile("output_PlotWidths.root","recreate");

  // Make histograms for plots
  TH2D *hdEdxResRange_all = new TH2D("hdEdxResRange_all","All particles;Residual Range (cm);dE/dx per hit (units?)",200,0,200,500,0,50);
  TH1D *hdEdx_all = new TH1D("hdEdx_all","All particles;dE/dx per hit (units?);No. hits",500,0,50);
  TH2D *hdEdxResRange_Muon = new TH2D("hdEdxResRange_Muon","True Muons Only;Residual Range (cm);dE/dx per hit (units?)",200,0,200,500,0,50);
  TH1D *hdEdx_Muon = new TH1D("hdEdx_Muon","True Muons Only;dE/dx per hit (units?);No. hits",500,0,50);
  TH2D *hdEdxResRange_Proton = new TH2D("hdEdxResRange_Proton","True Protons Only;Residual Range (cm);dE/dx per hit (units?)",200,0,200,500,0,50);
  TH1D *hdEdx_Proton = new TH1D("hdEdx_Proton","True Protons Only;dE/dx per hit (units?);No. hits",500,0,50);



  // Begin event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

    bool isData = ev.eventAuxiliary().isRealData();

    std::cout << "------ Processing "
	      << "Run " << ev.eventAuxiliary().run() << ", "
	      << "Event " << ev.eventAuxiliary().event() << ", "
	      << "Time " << ev.eventAuxiliary().time().timeHigh() << ", "
	      << "IsData " << isData << " ------" << std::endl;


    art::InputTag trackTag;
    art::InputTag caliTag;
    if (isData){
      trackTag=art::InputTag("pandoraNu::DataRecoStage2");
      caliTag=art::InputTag("pandoraNucali::DataRecoCali"); // Change this when data calibration available
    }
    else{
      trackTag=art::InputTag("pandoraNu::McRecoStage2");
      caliTag=art::InputTag("pandoraNucali");
    }

    //const auto& mcpHandle = ev.getValidHandle< std::vector< simb::MCParticle > >("largeant");
    //const auto& hitHandle = ev.getValidHandle< std::vector<recob::Hit> >(hitTag);

    // Get handles to needed information
    gallery::ValidHandle<std::vector<recob::Track>> trackHandle = ev.getValidHandle< std::vector<recob::Track> >(trackTag);
    std::vector<recob::Track> trackVec(*trackHandle);
    std::cout << "trackVec.size() = " << trackVec.size() << std::endl;

    //art::FindManyP<recob::Hit> trackHitAssn(trackHandle, ev, trackTag);
    art::FindManyP<anab::Calorimetry> trackCaliAssn(trackHandle, ev, caliTag);
    std::cout << "trackCaliAssn.size() = " << trackCaliAssn.size() << std::endl;
    bool nocaliobj = true;
    for (size_t i=0; i<trackCaliAssn.size(); i++){
      std::cout << "trackCaliAssn.at(" << i << ").size() = " << trackCaliAssn.at(i).size() << std::endl;
      if (trackCaliAssn.at(i).size() != 0) nocaliobj = false;
    }
    if (nocaliobj) continue;

    // Loop over tracks in event
    for (size_t i_tr = 0; i_tr < trackVec.size(); i_tr++){

      recob::Track track = trackVec.at(i_tr);
      std::cout << "track.ID() = " << track.ID() << std::endl;

      int pdgCode = 0;
      if (!isData){
	// Get truth information by looping through MCParticles and looking for MCParticles
	// that have start and end points within a given distance of the reconstructed start and
	// end. This has two effects: lets us find the true particle and rejects tracks that are
	// badly reconstructed. That means that these distributions will be valid for well-
	// -reconstructed tracks only
	const auto& mcpHandle = ev.getValidHandle< std::vector< simb::MCParticle > >("largeant");
	simb::MCParticle mcp;
	bool isDefined = false;
	std::cout << "Looping MCParticles..." << std::endl;
	for(auto const& thisMcp : (*mcpHandle)){

	  if (thisMcp.StatusCode() != 1 || !isInTPC(thisMcp.EndX(), thisMcp.EndY(), thisMcp.EndZ()) || !isInTPC(thisMcp.Vx(), thisMcp.Vy(), thisMcp.Vz())) continue;

	  if (isWellReconstructed(track, thisMcp)){
	    isDefined = true;
	    mcp = thisMcp;
	    break;
	  }

	}

	if (!isDefined) continue; // Skip not-well-reconstructed tracks
	std::cout << "Found good match" << std::endl;
	pdgCode = mcp.PdgCode();
	std::cout << "PDG code = " << pdgCode << std::endl;
      }

      // Get calorimetry object
      std::vector<art::Ptr<anab::Calorimetry>> caliFromTrack = trackCaliAssn.at(track.ID());
      art::Ptr<anab::Calorimetry> cali;
      for (auto c : caliFromTrack){
	int planenum = c->PlaneID().Plane;
	if (planenum != 2) continue; // only use calorimetry from collection plane (plane 2)
	cali = c;
      }
      if (!cali){
	std::cout << "Cali not set, skipping track" << std::endl;
	continue;
      }
      std::cout << "CALI USING PLANE: " << cali->PlaneID() << std::endl;

      std::vector<double> dEdx     = cali->dEdx();
      std::vector<double> resRange = cali->ResidualRange();


      // Now loop through hits in cali object
      for (size_t i_hit = 0; i_hit < resRange.size(); i_hit++){

	// All tracks: select hits with RR between 100 and 150 cm
	if (resRange.at(i_hit)>=100 && resRange.at(i_hit)<=150){
	  hdEdxResRange_all->Fill(resRange.at(i_hit), dEdx.at(i_hit));
	  hdEdx_all->Fill(dEdx.at(i_hit));
	} // end if RR 100-150cm

	if (isData) continue; // Truth information for MC only

	// Muons: select hits with RR between 100 and 150 cm
	if (abs(pdgCode) == 13 && resRange.at(i_hit)>=100 && resRange.at(i_hit)<=150){
	  hdEdxResRange_Muon->Fill(resRange.at(i_hit), dEdx.at(i_hit));
	  hdEdx_Muon->Fill(dEdx.at(i_hit));
	}

	// Protons: select hits with RR above 30 cm
	if (pdgCode == 2212 && resRange.at(i_hit)>30){
	  hdEdxResRange_Proton->Fill(resRange.at(i_hit), dEdx.at(i_hit));
	  hdEdx_Proton->Fill(dEdx.at(i_hit));
	}

      } // end loop over hits in cali object

    }// end track loop


  } // end event loop

  fout->cd();
  hdEdxResRange_all->Write();
  hdEdx_all->Write();
  hdEdxResRange_Muon->Write();
  hdEdx_Muon->Write();
  hdEdxResRange_Proton->Write();
  hdEdx_Proton->Write();
}
