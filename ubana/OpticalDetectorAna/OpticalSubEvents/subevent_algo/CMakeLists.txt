cet_make_library(
  SOURCE
  CosmicWindowSubEvents.cc
  SubEventModConfig.cc
  SubEventModule.cc
  WaveformData.cc
  pedestal.cc
  scintresponse.cc
  LIBRARIES
  PUBLIC
  ubana::OpticalDetectorAna_OpticalSubEvents_cfdiscriminator_algo
  ubobj::Optical
)

install_headers()
install_source()
