cet_build_plugin(
  InclusiveSinglePhotonGenie art::EDFilter
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::ServiceUtil
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  InclusiveSinglePhoton art::EDFilter
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::ServiceUtil
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
