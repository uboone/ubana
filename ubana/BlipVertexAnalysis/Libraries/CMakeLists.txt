FILE( GLOB local_src_libs *.cxx)

art_make_library( LIBRARY_NAME ubana_BlipVertexAnalysis_Libaries
        SOURCE ${local_src_libs}

        LIBRARIES
          art_Persistency_Common canvas
          art_Persistency_Provenance canvas
          art_Utilities canvas
          cetlib cetlib_except
          nusimdata_SimulationBase
          ${LARRECO_LIB}
          ${LARDATA_LIB}
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_FRAMEWORK_SERVICES_BASIC}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
          ${MF_MESSAGELOGGER}
          ${MF_UTILITIES}
          ${FHICLCPP}
          ${ROOT_BASIC_LIB_LIST}

)
