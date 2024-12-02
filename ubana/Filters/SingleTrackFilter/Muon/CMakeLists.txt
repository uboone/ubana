cet_build_plugin(
  SingleMuonFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::Utilities
  larevt::SpaceChargeService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
