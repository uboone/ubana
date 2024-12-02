cet_build_plugin(
  CoarseCCInclusiveFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::Filters_CoarseCCInclusiveFilter_helpers
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
