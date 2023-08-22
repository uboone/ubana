cet_build_plugin(
  SingleTrackFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::Filters_SingleTrackFilter_helpers
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(helpers)
add_subdirectory(Muon)

install_headers()
install_source()
install_fhicl()
