cet_build_plugin(
  CTMDataONAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  larevt::SpaceChargeService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_fhicl()
install_source()
