art_make( BASENAME_ONLY
	LIB_LIBRARIES
		larreco_RecoAlg
		larreco_RecoAlg_Cluster3DAlgs
		larsim_Simulation
		lardataobj_RawData
		lardata_RecoObjects
		lardataobj_RecoBase
		lardataobj_MCBase
		lardataobj_AnalysisBase
		larcorealg_Geometry
    larcore_Geometry_Geometry_service
		lardata_Utilities
		nusimdata_SimulationBase
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		${ART_ROOT_IO_TFILESERVICE_SERVICE}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
		${ART_UTILITIES}
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${FHICLCPP}
		cetlib cetlib_except
		${CLHEP}

		${G4_LIB_LIST}

	MODULE_LIBRARIES
		PreSelectAlgorithms
		PreSelectDataObjects
		larreco_RecoAlg
		larreco_RecoAlg_Cluster3DAlgs
		larsim_Simulation
		lardataobj_RawData
		lardata_RecoObjects
		lardataobj_RecoBase
		lardataobj_MCBase
		lardataobj_AnalysisBase
		larcorealg_Geometry
    larcore_Geometry_Geometry_service
		lardata_Utilities
		nusimdata_SimulationBase
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		${ART_ROOT_IO_TFILESERVICE_SERVICE}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
		${ART_UTILITIES}
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${FHICLCPP}
		cetlib cetlib_except
		${CLHEP}

		${G4_LIB_LIST}

	)

add_subdirectory(Algorithms)
add_subdirectory(DataObjects)

install_headers()
install_source()
install_fhicl()
