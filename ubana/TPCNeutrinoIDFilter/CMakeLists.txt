# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
  LIBRARY_NAME TPCNeutrinoIDFilter
  LIB_LIBRARIES
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_Utilities
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
  art_Utilities canvas
  ${MF_MESSAGELOGGER}

  ${FHICLCPP}
  cetlib cetlib_except


  MODULE_LIBRARIES
  TPCNeutrinoIDAlgs
  lardataalg_DetectorInfo
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_Utilities
  nusimdata_SimulationBase
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
  art_Utilities canvas
  ${MF_MESSAGELOGGER}

  ${FHICLCPP}
  cetlib cetlib_except


  )

add_subdirectory(Algorithms)

install_headers()
install_source()
install_fhicl()
