cet_build_plugin(
  AnalysisTree art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
)

add_subdirectory(ObjectTrees)
add_subdirectory(MCTruth)

install_fhicl()
install_source()
