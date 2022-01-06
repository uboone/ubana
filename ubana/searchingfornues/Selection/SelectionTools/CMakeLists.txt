cet_enable_asserts()

art_make( TOOL_LIBRARIES lardataobj_RecoBase
                         lardataobj_Simulation
                         lardataobj_MCBase
                         nusimdata::SimulationBase
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         lardata_RecoObjects
                         lardata_Utilities
                         larreco_RecoAlg
                         larpandora_LArPandoraInterface
                         art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         art::Framework_Services_Optional_RandomNumberGenerator_service
                         art::Persistency_Common
                         art::Persistency_Provenance
                         art::Utilities
                         canvas::canvas
                         messagefacility::MF_MessageLogger
                         fhiclcpp::fhiclcpp
                         cetlib::cetlib
                         cetlib_except::cetlib_except
                         Boost::filesystem
                         CLHEP::CLHEP
                         ROOT::EG
        )

install_headers()
install_source()
install_fhicl()
