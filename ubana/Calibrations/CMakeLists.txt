cet_make_library(
  SOURCE
  CalibrationTPC_Algs.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  lardataobj::RawData
)

cet_build_plugin(
  CalibrationTPC art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::Calibrations
  ubevt::CalData
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CalibrationdEdX art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Database_TPCEnergyCalib
  ubobj::DataOverlay
  larreco::Calorimetry
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  lardataobj::AnalysisBase
)

cet_build_plugin(
  Diffusion art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  ExtractPedestalInformation art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  FilterCrossingMuons art::EDFilter
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  Kplane2 art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  LifetimeQAQC art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larsim::Simulation
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
  ROOT::Minuit
  ROOT::Gpad
)

cet_build_plugin(
  Lifetime art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larsim::Simulation
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  NoiseFilter art::EDProducer
  LIBRARIES
  PRIVATE
  larevt::DetPedestalService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  ValidateLifetime art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larreco::Calorimetry
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  XYZcorrection art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_BackTrackerService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  XYZvalidatioin art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

add_subdirectory(Optical)

FILE(GLOB data_files calibrationList_*.dat )

install_headers()
install_fhicl()
install_source()
install_fw( LIST ${data_files} )
