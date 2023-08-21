cet_build_plugin(
  HitEfficiencyAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::DetPedestalService
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

add_subdirectory(tools)

#install_headers()
install_fhicl()
install_source()
