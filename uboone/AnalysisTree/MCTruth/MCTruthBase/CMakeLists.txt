cet_enable_asserts()

art_make( LIB_LIBRARIES lardataobj_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardata_Utilities
                        nusimdata_SimulationBase
                        canvas_Persistency_Common
                        canvas_Persistency_Provenance
                        canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${Boost_FILESYSTEM_LIBRARY}
                        ${Boost_SYSTEM_LIBRARY}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
        )

install_headers()
install_source()
install_fhicl()
