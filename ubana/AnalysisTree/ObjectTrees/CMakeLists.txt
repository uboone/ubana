cet_build_plugin(
  HitTreeAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  MCNeutrinoAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  POTSummaryAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larcoreobj::SummaryData
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  UbooneOpticalFilterAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Optical
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  WireTreeAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
