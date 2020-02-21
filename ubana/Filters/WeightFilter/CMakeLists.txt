art_make( 
          LIBRARY_NAME uboonecode_uboone_Filters
          # Needed for standalone c++ helper files
          LIB_LIBRARIES 

          # Needed for the main module          
          MODULE_LIBRARIES
                    larcorealg_Geometry
                    larreco_RecoAlg
                    larcore_Geometry_Geometry_service
                    larsim_Simulation lardataobj_Simulation
                    lardataobj_RawData
                    lardataobj_RecoBase
                    lardataobj_AnalysisBase
                    lardata_Utilities
                    nusimdata_SimulationBase
                    larsim_EventWeight_Base
                    ${ART_FRAMEWORK_CORE}
                    ${ART_FRAMEWORK_PRINCIPAL}
                    ${ART_FRAMEWORK_SERVICES_REGISTRY}
                    ${ART_FRAMEWORK_SERVICES_BASIC}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                    art_Framework_Services_Registry
                    art_Persistency_Common
                    art_Persistency_Provenance
                    art_Utilities
                    canvas
                    ${MF_MESSAGELOGGER}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    ${Boost_SYSTEM_LIBRARY}
                    ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
install_source()
install_fhicl()
