# basic source code CMakeLists.txt

art_make(BASENAME_ONLY
	MODULE_LIBRARIES
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
		art_Persistency_Common
		art_Persistency_Provenance
		art_Utilities
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${FHICLCPP}
		cetlib cetlib_except
		${CLHEP}
		${ROOT_BASIC_LIB_LIST}
		${G4_LIB_LIST}
	)

install_headers()
install_source()
add_subdirectory(Fcl)