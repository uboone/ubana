cet_build_plugin(
  NeutrinoSelectionFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)
add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
