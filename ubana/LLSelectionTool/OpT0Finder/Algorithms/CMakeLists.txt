art_make( BASENAME_ONLY
    	  LIBRARY_NAME     ubana_SelectionTool_OpT0FinderAlgorithms
	  LIB_LIBRARIES
                	  larsim_PhotonPropagation
			  larsim_PhotonPropagation_PhotonVisibilityService_service
			  larcorealg_Geometry
			  ubcore::ubana_BasicTool_GeoAlgo
			  ubreco_SelectionTool_OpT0FinderBase
			  ubobj_Trigger
                          ubcore::UBOpReadoutMap_service
  		          art::Framework_Core
   		          art::Framework_Principal
                          art::Framework_Services_Registry
   		          art_root_io::tfile_support ROOT::Core
     		          art_root_io::TFileService_service
			  art::Persistency_Common canvas::canvas
			  art::Persistency_Provenance canvas::canvas
			  art::Utilities canvas::canvas
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib
                          cetlib_except::cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
