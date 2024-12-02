cet_build_plugin(
  AnodeCathodeTrackAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::PatternFilter_PMAlgs
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  AnodeCathodeTrackFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::PatternFilter_PMAlgs
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
)

install_headers()
install_source()
install_fhicl()
