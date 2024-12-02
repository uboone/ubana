cet_make_library(
  SOURCE
  LEEPreCutAlgo.cxx
  LIBRARIES
  PUBLIC
)

cet_build_plugin(
  DLPMTPreCuts art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::OpticalDetectorAna_DLPMTPreCuts
  ubevt::Database
  ubobj::Optical
  larevt::CalibrationDBI_IOVData
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
)

add_subdirectory(fcl)

install_fhicl()
install_headers()
install_source()

