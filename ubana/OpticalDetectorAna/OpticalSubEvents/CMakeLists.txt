cet_build_plugin(
  SubEventBuilder art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::OpticalDetectorAna_OpticalSubEvents_subevent_algo
  ubobj::Optical
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
)

add_subdirectory( cfdiscriminator_algo )
add_subdirectory( subevent_algo )

install_headers()
install_source()

install_fhicl()
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
