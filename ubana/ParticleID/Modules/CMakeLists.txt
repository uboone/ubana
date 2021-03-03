include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make( BASENAME_ONLY
  LIBRARY_NAME ubana_ParticleID
  LIB_LIBRARIES
  ubana_ParticleID_Algorithms

  MODULE_LIBRARIES 
  ubobj_UBXSec
  uboonecode_uboone_UBXSec_Algorithms
  ubana_ParticleID_Algorithms
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
  larsim_MCCheater_BackTracker
  lardata_Utilities
  larevt_Filters
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  larana_TruncatedMean_Algorithm
  lardata_RecoObjects
  larreco_RecoAlg
  larpandora_LArPandoraInterface
  SignalShapingServiceMicroBooNE_service
  nusimdata_SimulationBase
  ubraw_RawData
  ${LARRECO_LIB}
  ${LARDATA_LIB}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  #art_Persistency_Common canvas_Persistency_Common
  #art_Persistency_Provenance canvas_Persistency_Provenance
  #art_Utilities canvas_Utilities
  art_Persistency_Common canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  ${CETLIB}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_BASIC_LIB_LIST}
  )   

### install_headers()
install_fhicl()
install_source()
