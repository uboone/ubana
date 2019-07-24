art_make( MODULE_LIBRARIES larcorealg_Geometry
			               larcore_Geometry_Geometry_service
                           larsim_Simulation
                           nug4_ParticleNavigation
                           lardataobj_Simulation
	                       lardata_Utilities
			               larevt_Filters
			               lardataobj_RawData
			               lardataobj_RecoBase
			               lardata_RecoObjects
			               larreco_RecoAlg
			               larpandora_LArPandoraInterface
			               SignalShapingServiceMicroBooNE_service
                           nusimdata_SimulationBase
                           ubraw_RawData
			               ubobj_Trigger
			               ubobj_Optical
                           ${LARRECO_LIB}
                           ${LARDATA_LIB}
                           ${ART_FRAMEWORK_CORE}
			               ${ART_FRAMEWORK_PRINCIPAL}
			               ${ART_FRAMEWORK_SERVICES_REGISTRY}
			               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			               art_Persistency_Common canvas
			               art_Persistency_Provenance canvas
			               art_Utilities canvas
			               ${MF_MESSAGELOGGER}
			               
	                       ${FHICLCPP}
			               cetlib cetlib_except
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
			               ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
