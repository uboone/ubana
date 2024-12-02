cet_build_plugin(
  CommonOpticalFilter art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Optical
  lardataobj::MCBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  TruthFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::MCBase
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
)

install_headers()
install_fhicl()
install_source()
