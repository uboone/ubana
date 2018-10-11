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
		${ART_FRAMEWORK_SERVICES_BASIC}
		${ART_FRAMEWORK_SERVICES_OPTIONAL}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
		${ROOT_BASIC_LIB_LIST}
		${G4_LIB_LIST}
		${Boost_SYSTEM_LIBRARY}
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
		${ART_FRAMEWORK_SERVICES_BASIC}
		${ART_FRAMEWORK_SERVICES_OPTIONAL}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
		${ROOT_BASIC_LIB_LIST}
		${G4_LIB_LIST}
		${Boost_SYSTEM_LIBRARY}
	)

add_subdirectory(Algorithms)
add_subdirectory(DataObjects)
add_subdirectory(Fcl)

install_headers()
install_source()
install_fhicl()
