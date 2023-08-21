include_directories( $ENV{SWTRIGGER_INCDIR} )
link_libraries( -L$ENV{SWTRIGGER_LIBDIR} -lFEMBeamTrigger -lSWTriggerBase )

cet_make_library(
  SOURCE
  FEMBasicTriggerAlgo.cxx
  LogicPulseFinder.cxx
  PulseRecoTemplate.cxx
  LIBRARIES
  PUBLIC
)

cet_build_plugin(
  FEMemulator art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::Calibrations_Optical
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  lardataobj::RawData
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  SPEcalibration art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Trigger
  ubana::Calibrations_Optical
  ubana::OpticalDetectorAna_OpticalSubEvents_cfdiscriminator_algo
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  larcore::ServiceUtil
  lardataobj::RawData
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  SWTrigger art::EDFilter
  LIBRARIES
  PRIVATE
  ubobj::Trigger
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  larcore::ServiceUtil
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(fcl)

install_fhicl()
install_headers()
install_source()
