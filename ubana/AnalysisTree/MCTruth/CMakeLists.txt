cet_build_plugin(
  MCTruthTestAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_make_library(
  LIBRARY_NAME IMCTruthMatching INTERFACE
  SOURCE IMCTruthMatching.h
  LIBRARIES
  INTERFACE
  lardataobj::RecoBase
  nusimdata::SimulationBase
  nug4::ParticleNavigation
  canvas::canvas
)

cet_write_plugin_builder(
  uboone_tool::IMCTruthMatching art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::IMCTruthMatching
)

include(uboone_tool::IMCTruthMatching)

cet_build_plugin(
  AssociationsTruth uboone_tool::IMCTruthMatching
  LIBRARIES
  PRIVATE
  ubana::AnalysisTree_MCTruth_MCTruthBase
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  art::Framework_Principal
)

cet_build_plugin(
  BackTrackerTruth uboone_tool::IMCTruthMatching
  LIBRARIES
  PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_BackTrackerService_service
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
)

add_subdirectory(MCTruthBase)

install_headers()
install_source()
install_fhicl()
