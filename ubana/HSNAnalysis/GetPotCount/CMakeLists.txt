cet_build_plugin(
  GetNormalizationHistograms art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  GetPotCount art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_source()
install_fhicl()
