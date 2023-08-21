cet_build_plugin(
  ACPTTrigFlashMatch art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Base
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CosmicRejection art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Base
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  MuCSFlashMatch art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Base
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_make_library(
  LIBRARY_NAME FlashMatchingToolBase INTERFACE
  SOURCE FlashMatchingToolBase_tool.h
  LIBRARIES
  INTERFACE
  ubreco::LLSelectionTool_OpT0Finder_Base
  larpandora::LArPandoraInterface
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::FlashMatchingToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::FlashMatchingToolBase
)

include(uboone_tool::FlashMatchingToolBase)

cet_build_plugin(
  FlashMatchingTool uboone_tool::FlashMatchingToolBase
  LIBRARIES
  PRIVATE
  ubobj::DataOverlay
  ubevt::Database
  larevt::PmtGainService
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
