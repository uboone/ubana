
#
simple_plugin( CellTreeUB "module"
               UBOpReadoutMap_service
               lardata_Utilities_DatabaseUtil_service
               lardataobj_RawData
	       lardataobj_RecoBase
	       larevt_CalibrationDBI_Providers
	       larevt_CalibrationDBI_IOVData
	       larcorealg_Geometry
	       larcore_Geometry_Geometry_service
	       uboone_RawData
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
               BASENAME_ONLY
             )

art_make( EXCLUDE CellTreeUB_module.cc
          BASENAME_ONLY
          LIBRARY_NAME uboone_MicroBooNEWireCell
          LIB_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
                        lardata_Utilities
                        uboone_Geometry
			uboone_RawData
                        UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
			lardataobj_RecoBase
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
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
                        uboone_Geometry
			uboone_RawData
                        UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
			lardataobj_RecoBase
                        larcoreobj_SummaryData
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			uboone_Database
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
			uboone_Database
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
