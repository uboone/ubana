cet_add_compiler_flags(CXX -Wno-unused-but-set-variable)

cet_make_library(
  LIBRARY_NAME AnalysisToolBase INTERFACE
  SOURCE AnalysisToolBase.h
  LIBRARIES
  INTERFACE
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::AnalysisToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::AnalysisToolBase
  lardata::RecoBaseProxy
  lardataobj::RecoBase
)

include(uboone_tool::AnalysisToolBase)

cet_build_plugin(
  BDT uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubana::XGBoost
)

cet_build_plugin(
  CRTApproachAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  lardataobj::AnalysisBase
)

cet_build_plugin(
  CalorimetryAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  larevt::SpaceChargeServices_SpaceChargeServiceStandard_service
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::AnalysisBase
  lardataobj::MCBase
  ROOT::EG
)

cet_build_plugin(
  ContainmentAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  nusimdata::SimulationBase
  lardataobj::MCBase
)

cet_build_plugin(
  CosmicIP uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  DefaultAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::RawData
  lardataobj::MCBase
  nusimdata::SimulationBase
  ROOT::EG
)

cet_build_plugin(
  EventFilter uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  nusimdata::SimulationBase
)

cet_build_plugin(
  EventWeightTree uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
)

cet_build_plugin(
  FlashMatching uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubevt::Utilities
  ubobj::CRT
  larpandora::LArPandoraInterface
  larevt::SpaceChargeService
  lardata::Utilities
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  MCFilter uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
)

cet_build_plugin(
  MCS uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larpandora::LArPandoraInterface
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  nusimdata::SimulationBase
  Eigen3::Eigen
  ROOT::Tree
)

cet_build_plugin(
  NeutrinoEnergy uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  lardata::Utilities
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  NeutrinoTiming uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::PmtGainService
  lardata::Utilities
  lardataobj::RawData
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  ObviousCosmicFlashMatching uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubevt::Database
  ubobj::CRT
  larpandora::LArPandoraInterface
  larevt::SpaceChargeService
  lardata::Utilities
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  PMTNoise uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
)

cet_build_plugin(
  PMTWaveform uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubevt::Utilities
  ubobj::CRT
  larevt::PmtGainService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::AnalysisBase
  lardataobj::MCBase
)

cet_build_plugin(
  Pi0Tagger uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  Pi0TruthAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  SecondShowerTagger uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  lardata::Utilities
  nusimdata::SimulationBase
)

cet_build_plugin(
  ShowerAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  ShowerStartPoint uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  SimPhoton uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceChargeService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  SlicePurCompl uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  TrackAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  ubevt::Database_TPCEnergyCalib
  larreco::RecoAlg
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::MCBase
  nusimdata::SimulationBase
  ROOT::EG
)

cet_build_plugin(
  VertexAnalysis uboone_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  lardata::DetectorPropertiesService
)

install_headers()
install_source()
install_fhicl()
