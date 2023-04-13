
#
simple_plugin( CellTreeUB "module"
               UBOpReadoutMap_service
	       ubobj_Trigger
	       ubobj_MuCS
	       nusimdata_SimulationBase
               lardata_Utilities_DatabaseUtil_service
               lardataobj_RawData
	       lardataobj_RecoBase
	       larevt_CalibrationDBI_Providers
	       larevt_CalibrationDBI_IOVData
	       larcorealg_Geometry
	       larcore_Geometry_Geometry_service
	       larsim_EventWeight_Base
	       ubraw_RawData
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
	       ${MF_MESSAGELOGGER}
	       ${ROOT_BASIC_LIB_LIST}
               BASENAME_ONLY
             )

simple_plugin( WireCellAnaTree "module"
               ubobj_WcpPort
               larcorealg_Geometry
               larcore_Geometry_Geometry_service
               lardataobj_RawData
               lardataobj_RecoBase
               lardataobj_MCBase
               nusimdata_SimulationBase
               lardata_Utilities
               lardataalg_DetectorInfo
               larreco_RecoAlg
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
	       larevt_CalibrationDBI_IOVData
               larevt_CalibrationDBI_Providers
               
               ${FHICLCPP}
               cetlib cetlib_except
               ${ROOT_GEOM}
               ${ROOT_XMLIO}
               ${ROOT_GDML}
               ${ROOT_EG}
               ${ROOT_BASIC_LIB_LIST}
               BASENAME_ONLY
             )

simple_plugin( WireCellEventWeightTree "module"
               ubobj_WcpPort
               larcorealg_Geometry
               larcore_Geometry_Geometry_service
               lardataobj_RawData
               lardataobj_RecoBase
               lardataobj_MCBase
               nusimdata_SimulationBase
               lardata_Utilities
               lardataalg_DetectorInfo
               larreco_RecoAlg
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
               BASENAME_ONLY
             )

simple_plugin( WireCellPFDump "module"
               ubobj_WcpPort
               larcorealg_Geometry
               larcore_Geometry_Geometry_service
               lardataobj_RawData
               lardataobj_RecoBase
               lardataobj_MCBase
               nusimdata_SimulationBase
               lardata_Utilities
               lardataalg_DetectorInfo
               larreco_RecoAlg
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
               BASENAME_ONLY
             )


art_make( EXCLUDE CellTreeUB_module.cc WireCellAnaTree_module.cc WireCellEventWeightTree_module.cc WireCellPFDump_module.cc
          BASENAME_ONLY
          LIBRARY_NAME uboone_MicroBooNEWireCell
          LIB_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
                        lardata_Utilities
                        ubcore_Geometry
			ubraw_RawData
                        UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
			lardataobj_RecoBase
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			larsim_EventWeight_Base
                        ${UBDut}
                        ${PQ}
                        ${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}

          MODULE_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        ubcore_Geometry
			ubraw_RawData
                        UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
			lardataobj_RecoBase
                        larcoreobj_SummaryData
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			larsim_EventWeight_Base
			ubevt_Database
                        uboone_MicroBooNEWireCell
                        ${UBDut}
                        ${PQ}
                        ${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			ubevt_Database
			${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art_Utilities canvas
                              larevt_CalibrationDBI_IOVData
                              lardataobj_RawData
                              larcorealg_Geometry
			      larcore_Geometry_Geometry_service
                              ${ART_FRAMEWORK_SERVICES_REGISTRY}
                              ${ART_FRAMEWORK_CORE}
                              ${ART_FRAMEWORK_PRINCIPAL}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                              art_Persistency_Common canvas
                              art_Persistency_Provenance canvas
                              art_Utilities canvas
                              ${FHICLCPP}
                              ${ROOT_BASIC_LIB_LIST}
)                                                           

install_headers()
install_source()

add_subdirectory(job)
add_subdirectory(SinglePhoton)
