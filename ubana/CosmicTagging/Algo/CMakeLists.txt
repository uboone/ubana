
art_make( BASENAME_ONLY
    	  LIBRARY_NAME     ubana_CosmicTagging_CosmicTaggingAlgorithms
	  LIB_LIBRARIES
                	  larsim_PhotonPropagation
			  larsim_PhotonPropagation_PhotonVisibilityService_service
			  ubcore::ubana_BasicTool_GeoAlgo
			  ubreco_SelectionTool_OpT0FinderBase
                          ubcore::UBOpReadoutMap_service
  		          art::Framework_Core
   		          art::Framework_Principal
                          art::Framework_Services_Registry
   		          art_root_io::tfile_support ROOT::Core
     		          art_root_io::TFileService_service
     		          art::Persistency_Common
                          art::Persistency_Provenance
                          art::Utilities
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
