cet_build_plugin(
  HsnMcTruthInformation art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  McTruthInformation art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_source()
add_subdirectory(Fcl)
