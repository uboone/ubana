art_make( BASENAME_ONLY

	  EXCLUDE
			plotMV4.C
			plotsupimp.C

          LIB_LIBRARIES
			ubraw_RawData
			ubobj_Trigger
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			ubcore::ubana_BasicTool_GeoAlgo
                        lardata_Utilities
                        larevt_Filters
                        lardataalg_DetectorInfo
                        lardataobj_RawData
			lardataobj_RecoBase
			larcoreobj_SummaryData
			larsim_MCCheater_BackTrackerService_service
			larsim_MCCheater_ParticleInventoryService_service
			larpandora_LArPandoraInterface
			nusimdata::SimulationBase
			UB_DATA_TYPES
			Boost::serialization
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP


          MODULE_LIBRARIES
			ubraw_RawData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			ubcore::ubana_BasicTool_GeoAlgo
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
			lardataobj_RecoBase
			larcoreobj_SummaryData
			larsim_MCCheater_BackTrackerService_service
			larpandora_LArPandoraInterface
			nusimdata::SimulationBase
			ubana_LEEPhotonAnalysis
			UB_DATA_TYPES
			Boost::serialization
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP

        )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
