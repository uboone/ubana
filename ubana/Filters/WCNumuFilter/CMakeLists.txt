cet_build_plugin(
  WCNumuFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larpandora::LArPandoraInterface
  larevt::SpaceChargeService
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_fhicl()
install_headers()
install_source()
