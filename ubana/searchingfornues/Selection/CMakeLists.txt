#include_directories( $ENV{PANDORA_INC} )
#include_directories( $ENV{LARPANDORACONTENT_INC} )

add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)
add_subdirectory(job)

art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nug4::ParticleNavigation
                           lardataobj_Simulation
			   lardataobj_RecoBase
			   lardataobj_MCBase	
                           lardata_Utilities
			   nusimdata::SimulationBase
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardataobj_RawData
                           larevt_CalibrationDBI_Providers
                           nurandom::RandomUtils_NuRandomService_service
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
                           CLHEP::CLHEP
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ${ROOT_BASIC_LIB_LIST}
            		   PANDORASDK
			   larpandora_LArPandoraInterface
        )

install_headers()
install_fhicl()
install_source()
