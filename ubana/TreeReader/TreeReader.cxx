///////////////////////////////////////////////////////////////////////
/// \file  TreeReader.cxx
/// \brief Source to read beam flux files
/// \author zarko@fnal.gov
///
/// Adapted from FluxReader by Zarko Pavlovic, by Adam Lister &
/// Andy Mastbaum.
////////////////////////////////////////////////////////////////////////

//LArSoft 
#include "larcoreobj/SummaryData/POTSummary.h"

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//root
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "uboone/TreeReader/TreeReader.h"
#include "GSimpleInterface.h"
#include "NTupleInterface.h"
#include "cetlib/exception.h"

namespace uboone {

TreeReader::TreeReader(fhicl::ParameterSet const& ps, 
                       art::ProductRegistryHelper &helper,
                       art::SourceHelper const &pm)
    : fSourceHelper(pm), fCurrentSubRunID() {
  // General parameters
  fEventCounter = 0; 
  fEntry = ps.get<uint32_t>("skipEvents", 0);
  fMaxEvents = ps.get<int>("maxEvents", -1);
  fInputType = ps.get<std::string>("inputType");
  fDataProductName = ps.get<std::string>("dataProductName", "TreeReader");
  fVerbose = ps.get<bool>("verbose", false);
  fDLMode = ps.get<bool>("dlMode", false);

  // Define which data products this module produces
  if (fInputType == "gsimple") {
    helper.reconstitutes<std::vector<simb::MCFlux>, art::InEvent>(fDataProductName);
    helper.reconstitutes<std::vector<simb::MCTruth>, art::InEvent>(fDataProductName);
    helper.reconstitutes<sumdata::POTSummary, art::InSubRun>(fDataProductName);
  }
  else if (fInputType == "ntuple") {
    helper.reconstitutes<std::vector<simb::MCFlux>, art::InEvent>(fDataProductName);
    helper.reconstitutes<std::vector<simb::MCTruth>, art::InEvent>(fDataProductName);
    helper.reconstitutes<std::vector<simb::GTruth>, art::InEvent>(fDataProductName);
    fTreeName = ps.get<std::string>("treeName");
    fBranchDef = ps.get<fhicl::ParameterSet>("branches");
  }
  else {
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Ntuple format " << fInputType << " not supported" << std::endl;
  }
}


void TreeReader::closeCurrentFile() {    
  mf::LogInfo(__FUNCTION__) << "File boundary (processed " << fEventCounter
                            << " events)" << std::endl;
  fCurrentSubRunID.flushSubRun();
  fEventCounter = 0;
  fInputFile->Close();
  delete fInputFile;
  fSkipEvents = 0; // If crossing file boundary don't skip events in next file
  fEntry = 0;
}


void TreeReader::readFile(std::string const &name, art::FileBlock* &fb) {
  // Fill and return a new Fileblock.
  fb = new art::FileBlock(art::FileFormatVersion(1, "TreeReader"), name);

  fInputFile = new TFile(name.c_str());
  if (fInputFile->IsZombie()) {
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Failed to open " << fInputFile << std::endl;
  } 

  // Set up the ROOT file
  if (fInputType == "gsimple") {
    fInterface = new GSimpleInterface();
    ((GSimpleInterface*)fInterface)->SetRootFile(fInputFile);
  } 
  else if (fInputType == "ntuple") {
    NTupleInterface* iface = new NTupleInterface();
    iface->SetVerbose(fVerbose);
    iface->SetDLMode(fDLMode);
    iface->SetRootFile(fInputFile, fTreeName, fBranchDef);
    fInterface = iface;
  }
  else {
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Ntuple format " << fInputType << " not supported" << std::endl;
  }

  std::cout << "File has " << fInterface->GetEntries() << " entries" << std::endl;
  std::cout << "POT = " << fInterface->GetPOT() << std::endl;
  std::cout << "Run = " << fInterface->GetRun() << std::endl;
  fPOT += fInterface->GetPOT();
}


bool TreeReader::readNext(art::RunPrincipal* const&,
                          art::SubRunPrincipal* const&,
                          art::RunPrincipal* &outR,
                          art::SubRunPrincipal* &outSR,
                          art::EventPrincipal* &outE) {
  if (fMaxEvents > 0 && fEventCounter >= unsigned(fMaxEvents)) {
    return false;
  }

  if (fEventCounter % 10000 == 0) {
    mf::LogInfo(__FUNCTION__)
      << "Attempting to read event: " << fEventCounter << std::endl;
  }

  std::unique_ptr<std::vector<simb::MCFlux> > mcfluxvec(new std::vector<simb::MCFlux>);
  std::unique_ptr<std::vector<simb::MCTruth> > mctruthvec(new std::vector<simb::MCTruth>);
  std::unique_ptr<std::vector<simb::GTruth> > gtruthvec(new std::vector<simb::GTruth>);

  simb::MCFlux flux;
  if (!fInterface->FillMCFlux(fEntry, flux)) {
    return false;
  }

  simb::MCTruth mctruth;
  if (!fInterface->FillMCTruth(fEntry, mctruth)) {
    return false;
  }

  simb::GTruth gtruth;
  if (!fInterface->FillGTruth(fEntry, gtruth)) {
    return false;
  }

  // Fake mctruth product to cheat eventweight that gets neutrino energy from it
  if (fInputType == "gsimple") {
    simb::MCParticle mcpnu(0, flux.fntype, "Flux");
    mcpnu.AddTrajectoryPoint(fInterface->GetNuPosition(), fInterface->GetNuMomentum());
    mctruth.Add(mcpnu);
    mctruth.SetNeutrino(0,0,0,0,0,0,0,0,0,0);
    mctruthvec->push_back(mctruth);
  }

  if (fInputType == "gsimple") {
    mcfluxvec->push_back(flux);
    mctruthvec->push_back(mctruth);
  }
  else if (fInputType == "ntuple") {
    mcfluxvec->push_back(flux);
    mctruthvec->push_back(mctruth);
    gtruthvec->push_back(gtruth);
  }
  else {
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Ntuple format " << fInputType << " not supported" << std::endl;
  }

  fEventCounter++;
  fEntry++;

  art::RunNumber_t rn = fInterface->GetRun();
  if (rn==0) rn=999999;
  art::Timestamp tstamp(time(0));

  art::SubRunID newID(rn, 0);
  if (fCurrentSubRunID.runID() != newID.runID()) {
    outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
  }

  if (fCurrentSubRunID != newID) {
    outSR = fSourceHelper.makeSubRunPrincipal(rn,0,tstamp);
    std::unique_ptr<sumdata::POTSummary> pot(new sumdata::POTSummary);    
    pot->totpot = fPOT;
    pot->totgoodpot = fPOT;
    fPOT=0;

    if (fInputType == "gsimple") {
      art::put_product_in_principal(std::move(pot), *outSR, fDataProductName);
    }

    fCurrentSubRunID = newID;
  }

  outE = fSourceHelper.makeEventPrincipal(
    fCurrentSubRunID.run(), fCurrentSubRunID.subRun(), fEventCounter, tstamp);

  // Put products in the event.
  if (fInputType == "gsimple") {
    art::put_product_in_principal(std::move(mcfluxvec), *outE, fDataProductName);
    art::put_product_in_principal(std::move(mctruthvec), *outE, fDataProductName);
  }
  else if (fInputType == "ntuple") {
    art::put_product_in_principal(std::move(mcfluxvec), *outE, fDataProductName);
    art::put_product_in_principal(std::move(mctruthvec), *outE, fDataProductName);
    art::put_product_in_principal(std::move(gtruthvec), *outE, fDataProductName);
  }
  else {
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Ntuple format " << fInputType << " not supported" << std::endl;
  }

  return true;
}

void TreeReader::endJob() {
  std::cout << "TreeReader job ended" << std::endl;
}

}  // namespace uboone

