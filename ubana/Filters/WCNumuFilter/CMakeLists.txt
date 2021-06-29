include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make( MODULE_LIBRARIES larcorealg_Geometry
	nusimdata_SimulationBase
	lardataobj_Simulation
	lardataobj_RecoBase
	lardataobj_AnalysisBase
	larcoreobj_SummaryData
	larpandora_LArPandoraInterface
	ubobj_WcpPort
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_PRINCIPAL}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_FRAMEWORK_SERVICES_OPTIONAL}
	${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	art_Persistency_Common 
	art_Persistency_Provenance
	art_Utilities
	canvas
	${MF_MESSAGELOGGER}

	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_GEOM}
	${ROOT_XMLIO}
	${ROOT_GDML}
	${ROOT_EG}
	${ROOT_BASIC_LIB_LIST}

)

install_fhicl()
install_headers()
install_source()
