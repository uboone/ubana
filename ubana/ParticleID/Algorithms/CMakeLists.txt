include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

cet_enable_asserts()

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
  LIBRARY_NAME     ubana_ubana_ParticleID_Algorithms
  LIB_LIBRARIES
  ubana_BasicTool_GeoAlgo
  uboone_UBFlashFinder
  ubana_SelectionTool_OpT0FinderBase
  ubana_SelectionTool_OpT0FinderAlgorithms
  larana_ParticleIdentification

  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  ubraw_RawData
  SignalShapingServiceMicroBooNE_service
  larevt_CalibrationDBI_IOVData
  larevt_CalibrationDBI_Providers
  larsim_MCCheater_BackTracker
  larsim_Simulation
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardataobj_Simulation
  lardata_Utilities
  larpandora_LArPandoraInterface
  nusimdata_SimulationBase
  ${PANDORASDK}
  ${PANDORAMONITORING}
  ${UBOONECODE_LIB}
  ${LARDATA_LIB}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  #art_Persistency_Common canvas_Persistency_Common
  #art_Persistency_Provenance canvas_Persistency_Provenance
  #art_Utilities canvas_Utilities
  art_Persistency_Common canvas
  art_Utilities canvas
  ${ART_PERSISTENCY_COMMON}
  ${ART_PERSISTENCY_PROVENANCE}
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  ${CETLIB}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_BASIC_LIB_LIST}

  MODULE_LIBRARIES
  ubana_ubana_ParticleID_Algorithms
  larsim_MCCheater_BackTracker
  larpandora_LArPandoraInterface
  #                   larpandora_LArPandoraContent

  ubcore_LLBasicTool_GeoAlgo
  ubreco_UBFlashFinder
  ubana_LLSelectionTool_OpT0Finder_Base
  ubana_LLSelectionTool_OpT0Finder_Algorithms

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
  nusimdata_SimulationBase
  ${UBOONECODE_LIB}
  ${LARDATA_LIB}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  ${ART_PERSISTENCY_PROVENANCE}
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  ${CETLIB}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_BASIC_LIB_LIST}
  )

install_fhicl()
install_headers()
install_source()
