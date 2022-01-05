cet_enable_asserts()

art_make( LIB_LIBRARIES lardataobj_RecoBase
                        lardataalg_DetectorInfo
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardata_Utilities
                        nusimdata::SimulationBase
                        canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        Boost::filesystem

                        CLHEP::CLHEP
        )

install_headers()
install_source()
install_fhicl()
