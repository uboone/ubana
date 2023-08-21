cet_make_library(
  SOURCE
  SEAviewer.cc
  LIBRARIES
  PUBLIC
  larevt::SpaceChargeServices_SpaceChargeServiceStandard_service
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  ROOT::Gpad
)

install_headers()
install_fhicl()
install_source()
