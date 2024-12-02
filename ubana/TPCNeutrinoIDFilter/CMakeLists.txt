cet_build_plugin(
  LowLevelNueFilter art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
)

cet_build_plugin(
  OpFlashFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
)

cet_build_plugin(
  TPCNeutrinoIDAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  TPCNeutrinoIDFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  TPCNeutrinoID art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::TPCNeutrinoIDFilter_Algorithms
  larsim::MCCheater_ParticleInventoryService_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

add_subdirectory(Algorithms)

install_headers()
install_source()
install_fhicl()
