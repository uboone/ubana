cet_build_plugin(
  CosmicFlashTaggerAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CosmicFlashTagger art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::CosmicTagging_Algo
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(Algo)
add_subdirectory(job)

install_fhicl()
install_source()
