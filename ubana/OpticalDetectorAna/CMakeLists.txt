cet_make_library(
  SOURCE
  OpDetWaveformAna.cxx
  LIBRARIES
  PUBLIC
  larana::OpticalDetector_OpHitFinder
)

cet_build_plugin(
  FlashTrigger art::EDFilter
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  OpDigitSaturationCorrection art::EDProducer
  LIBRARIES
  PRIVATE
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  SimpleFlashFinder art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(
  UBBasicOpticalAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::OpticalDetectorAna
  lardata::DetectorClocksService
  larcore::ServiceUtil
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(job)
add_subdirectory(OpticalSubEvents)
add_subdirectory(DLPMTPreCuts)

install_headers()
install_source()
install_fhicl()
