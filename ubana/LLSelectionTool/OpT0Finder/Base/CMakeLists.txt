art_make( BASENAME_ONLY
   	  LIBRARY_NAME     ubana_SelectionTool_OpT0FinderBase
	  LIB_LIBRARIES
	    larcorealg_Geometry
	    larcore_Geometry_Geometry_service
            lardataalg_DetectorInfo
        art::Framework_Core
        art::Framework_Services_Registry
        art::Framework_Principal
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        cetlib::cetlib
        CLHEP::CLHEP
        ubcore::ubana_BasicTool_GeoAlgo
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        cetlib_except::cetlib_except
	) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.21.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types
install_headers()
install_source()
