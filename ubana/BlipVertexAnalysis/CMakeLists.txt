link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)

art_make(
          MODULE_LIBRARIES 
          ubana_BlipVertexAnalysis_Libaries
          ubreco_BlipRecoAlg
          lardataobj_RecoBase
          larpandora_LArPandoraInterface

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
install_fhicl()
install_source()

add_subdirectory(jobs)
add_subdirectory(Libraries)
