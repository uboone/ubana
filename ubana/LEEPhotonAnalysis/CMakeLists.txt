link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)

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
			nusimdata_SimulationBase
			UB_DATA_TYPES
			${Boost_SERIALIZATION_LIBRARY}
                        nurandom_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}


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
			nusimdata_SimulationBase
			ubana_LEEPhotonAnalysis
			UB_DATA_TYPES
			${Boost_SERIALIZATION_LIBRARY}
                        nurandom_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}

        )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
