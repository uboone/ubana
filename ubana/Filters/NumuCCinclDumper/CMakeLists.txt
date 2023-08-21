cet_make_library(
  SOURCE
  BackTrackerTruthMatch.cxx
  PID.cxx
  LIBRARIES
  PUBLIC
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  ROOT::Core
)

cet_build_plugin(
  NumuCCinklRun3 art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::Filters_NumuCCinclDumper
  ubana::Filters_NumuCCinclDumper_helpers
  ubobj::CRT
  larevt::SpaceChargeService
  lardata::Utilities
  nurandom::RandomUtils_NuRandomService_service
  ROOT::Tree
)

add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
