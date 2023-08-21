cet_build_plugin(
  GammaTruthStudies art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::MCBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  ProtonHitPurity art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::MCBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  ProtonTruthStudies art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::MCBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  SaveSliceHits art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
)

cet_build_plugin(
  SecondShowerPurity art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::MCBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  ShowerMerger art::EDProducer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::RecoBaseProxy
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  lardataobj::MCBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
