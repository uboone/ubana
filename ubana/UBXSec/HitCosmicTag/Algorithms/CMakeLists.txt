# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
    	    LIBRARY_NAME     HitCosmicTag_Algorithms
	        LIB_LIBRARIES
                           HitCosmicTag_Base
  		                     art::Framework_Core
   		                     art::Framework_Principal

                           art::Framework_Services_Registry
   		                     art_root_io::tfile_support ROOT::Core
     		                   art_root_io::TFileService_service
     		                   art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
	      )

install_headers()
install_source()
