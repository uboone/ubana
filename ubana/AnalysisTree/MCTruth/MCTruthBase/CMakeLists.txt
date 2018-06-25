cet_enable_asserts()

art_make( LIB_LIBRARIES lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardata_Utilities
                        nusimdata_SimulationBase
                        canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        ${CETLIB}
                        cetlib_except
                        ${Boost_FILESYSTEM_LIBRARY}
                        ${Boost_SYSTEM_LIBRARY}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
        )

install_headers()
install_source()
install_fhicl()
