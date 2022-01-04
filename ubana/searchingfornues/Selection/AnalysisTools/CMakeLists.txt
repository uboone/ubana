cet_enable_asserts()

include_directories( $ENV{EIGEN_INC} )

art_make( TOOL_LIBRARIES lardataobj_RecoBase
                         lardataobj_Simulation
                         lardataobj_MCBase
                         nusimdata::SimulationBase
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         larpandora_LArPandoraInterface
                         lardata_RecoObjects
                         lardata_Utilities
			 larreco_RecoAlg
                         larsim_EventWeight_Base
                         art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         art::Framework_Services_Optional_RandomNumberGenerator_service
                         art::Persistency_Common
                         art::Persistency_Provenance
                         art::Utilities
                         canvas::canvas
			 hep_concurrency
                         messagefacility::MF_MessageLogger
                         fhiclcpp::fhiclcpp
                         cetlib::cetlib
                         cetlib_except::cetlib_except
                         Boost::filesystem
                         CLHEP::CLHEP
                         ${ROOT_BASIC_LIB_LIST}
                         ROOT::EG
                         ROOT::Minuit
                         ubreco_ShowerReco_ProximityClustering_Algorithms
                         uboonecode_uboone_XGBoost
                         art::Persistency_Common
                         art_Framework_Services_System_TriggerNamesService_service
        )

install_headers()
install_source()
install_fhicl()
