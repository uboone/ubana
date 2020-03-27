art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation
#nug4_ParticleNavigation
                           lardataobj_Simulation
                           lardata_Utilities
                           lardataobj_RecoBase
                           lardata_RecoObjects larreco_Deprecated
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardataobj_RawData
                           larreco_RecoAlg
                           larreco_Calorimetry
                           larevt_CalibrationDBI_Providers
                           larpandora_LArPandoraInterface
                           SignalShapingServiceMicroBooNE_service
                           nusimdata_SimulationBase
                           nutools_RandomUtils_NuRandomService_service
                           #nurandom_RandomUtils_NuRandomService_service
                           ubana_MyClasses
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           #${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           cetlib_except
                           ${CLHEP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
