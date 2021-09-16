art_make( 

      LIBRARY_NAME     uboonecode_uboone_HyperonProduction_Algorithms

      LIB_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}
      lardataobj_RecoBase
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
      cetlib cetlib_except
      uboonecode_uboone_EventConnectedness
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}

      DICT_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}

      MODULE_LIBRARIES
      larcorealg_Geometry
      ${ROOT_BASIC_LIB_LIST}

)

install_headers()
install_fhicl()
install_source()
