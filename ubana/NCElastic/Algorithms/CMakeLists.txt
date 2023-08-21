cet_make_library(
  SOURCE
  LREvtReconstruction.cxx
  TrackFeatures.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
