cet_build_plugin(
  CellTreeUB art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubevt::Database_TPCEnergyCalib
  ubcore::Geometry_UBOpReadoutMap_service
  ubobj::Trigger
  ubobj::MuCS
  larsim::EventWeight_Base
  larevt::PmtGainService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  WireCellAnaTree art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  larevt::PmtGainService
  larcore::ServiceUtil
  lardataobj::RawData
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
  ROOT::Hist
)

cet_build_plugin(
  WireCellEventWeightTree art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  lardataobj::RawData
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  WireCellPFDump art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataobj::RawData
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)


install_headers()
install_source()

add_subdirectory(job)
add_subdirectory(SinglePhoton)
