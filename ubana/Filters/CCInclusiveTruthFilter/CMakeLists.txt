cet_build_plugin(
  CCInclusiveTruthFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::Utilities
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
)

install_headers()
install_source()
install_fhicl()
