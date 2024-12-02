cet_build_plugin(
  MuCSGainAna art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubcore::Geometry
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  MuCST0Finder art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  MuCSTrackFinder art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  T0TrackCalib art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  TPCObjectCreator art::EDProducer
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  larsim::MCCheater_ParticleInventoryService_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  UBT0Finder art::EDProducer
  LIBRARIES
  PRIVATE
)

# install_headers()
install_fhicl()
install_source()
