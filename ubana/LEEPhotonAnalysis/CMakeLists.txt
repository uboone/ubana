cet_make_library(
  SOURCE
  DetectorObjects.cxx
  EnergyHelper.cxx
  FillTreeVariables.cxx
  GeometryHelper.cxx
  ParticleAssociations.cxx
  RecoMCMatching.cxx
  RecoTrueHelper.cxx
  SimchInfo.cxx
  VertexBuilder.cxx
  VertexQuality.cxx
  LIBRARIES
  PRIVATE
  lardataalg::UtilitiesHeaders
  PUBLIC
  ubcore::LLBasicTool_GeoAlgo
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larsim::EventWeight_Base
  larsim::MCCheater_ParticleInventoryService_service
  larevt::DetPedestalService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art_root_io::TFileService_service
  nug4::ParticleNavigation
  ROOT::Tree
)

cet_build_plugin(
  GENIETest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  LEEPhotonAnalysis art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::LEEPhotonAnalysis
  lardata::DetectorPropertiesService
  lardataobj::AnalysisBase
)

cet_build_plugin(
  NCDeltaRadiative art::EDFilter
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  RecoTrueTest art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::LEEPhotonAnalysis
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorClocksService
)

cet_build_plugin(
  SignalFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  flash art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
