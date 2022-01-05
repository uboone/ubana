include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

cet_enable_asserts()

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
  LIBRARY_NAME     ubana_ParticleID_Algorithms
  LIB_LIBRARIES
  ubana_BasicTool_GeoAlgo
  uboone_UBFlashFinder
  ubreco_SelectionTool_OpT0FinderBase
  ubreco_SelectionTool_OpT0FinderAlgorithms
  larana_ParticleIdentification

  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  ubraw_RawData
  SignalShapingServiceMicroBooNE_service
  larevt_CalibrationDBI_IOVData
  larevt_CalibrationDBI_Providers
  larsim_Simulation
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardataobj_Simulation
  lardata_Utilities
  larpandora_LArPandoraInterface
  nusimdata::SimulationBase
  PANDORASDK
  PANDORAMONITORING
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib::cetlib cetlib_except
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}

  MODULE_LIBRARIES
  ubana_ParticleID_Algorithms
  larsim_MCCheater_BackTracker_service
  larpandora_LArPandoraInterface
  #                   larpandora_LArPandoraContent

  ubana_BasicTool_GeoAlgo
  uboone_UBFlashFinder
  ubreco_SelectionTool_OpT0FinderBase
  ubreco_SelectionTool_OpT0FinderAlgorithms

  ubraw_RawData
  SignalShapingServiceMicroBooNE_service
  larevt_CalibrationDBI_IOVData
  larevt_CalibrationDBI_Providers
  larcore_Geometry
  larcore_Geometry_Geometry_service
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardata_Utilities
  larana_ParticleIdentification
  larsim_Simulation lardataobj_Simulation 
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}
  )

install_fhicl()
install_headers()
install_source()
