art_make( BASENAME_ONLY
    	    LIBRARY_NAME     HitCosmicTag_Algorithms
	        LIB_LIBRARIES
                           HitCosmicTag_Base
  		                     art::Framework_Core
   		                     art::Framework_Principal

                           art::Framework_Services_Registry
   		                     art_root_io::tfile_support
     		                   art_root_io::TFileService_service
     		                   art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           ROOT::Physics
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           ROOT::Core
                           ROOT::Eve
	      )

install_headers()
install_source()
