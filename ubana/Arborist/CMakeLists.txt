art_make( BASENAME_ONLY
          LIB_LIBRARIES
			ubraw_RawData
			ubobj_Trigger
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
			ubana_BasicTool_GeoAlgo
            lardata_Utilities
            larevt_Filters
            larpandoraobj
            lardataobj_RawData
			lardataobj_RecoBase
            lardataobj_MCBase
			larcoreobj_SummaryData
			larsim_MCCheater_BackTrackerService_service
			larsim_MCCheater_ParticleInventoryService_service
			larpandora_LArPandoraInterface
            PANDORASDK
            PANDORAMONITORING
            LArPandoraContent 
			nusimdata::SimulationBase
			UB_DATA_TYPES
			Boost::serialization
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP
                        ROOT::Geom 
          MODULE_LIBRARIES
			ubraw_RawData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			ubana_BasicTool_GeoAlgo
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
			lardataobj_RecoBase
            lardataobj_MCBase
            larcoreobj_SummaryData	
			larsim_MCCheater_BackTrackerService_service
			larpandora_LArPandoraInterface
            SignalShapingServiceMicroBooNE_service
			nusimdata::SimulationBase
			ubana_LEEPhotonAnalysis
			UB_DATA_TYPES
			Boost::serialization
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP
                        ROOT::Geom )
                      
install_headers()
install_fhicl()
install_source()
