
#
simple_plugin( CellTreeUB "module"
               ubcore::UBOpReadoutMap_service
	       ubobj_Trigger
	       ubobj_MuCS
	       nusimdata::SimulationBase
               lardata_Utilities_DatabaseUtil_service
               lardataobj_RawData
	       lardataobj_RecoBase
	       larevt_CalibrationDBI_Providers
	       larevt_CalibrationDBI_IOVData
	       larcorealg_Geometry
	       larcore_Geometry_Geometry_service
	       larsim_EventWeight_Base
	       ubraw_RawData
               art::Framework_Core
               art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common canvas::canvas
               art::Persistency_Provenance canvas::canvas
               art::Utilities canvas::canvas
	       messagefacility::MF_MessageLogger

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
                        ubcore::UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
			lardataobj_RecoBase
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			larsim_EventWeight_Base
                        ${UBDut}
                        PostgreSQL::PostgreSQL
                        Boost::serialization
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml



          MODULE_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        ubcore_Geometry
			ubraw_RawData
                        ubcore::UBOpReadoutMap_service
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
                        PostgreSQL::PostgreSQL
                        Boost::serialization
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml


	  SERVICE_LIBRARIES
			larevt_CalibrationDBI_Providers
			larevt_CalibrationDBI_IOVData
			ubevt_Database
			art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art::Utilities canvas::canvas
                              larevt_CalibrationDBI_IOVData
                              lardataobj_RawData
                              larcorealg_Geometry
			      larcore_Geometry_Geometry_service
                              art::Framework_Services_Registry
                              art::Framework_Core
                              art::Framework_Principal
                              art_root_io::tfile_support ROOT::Core
                              art_root_io::TFileService_service
                              art::Persistency_Common canvas::canvas
                              art::Persistency_Provenance canvas::canvas
                              art::Utilities canvas::canvas
                              fhiclcpp::fhiclcpp

)

install_headers()
install_source()

add_subdirectory(job)
