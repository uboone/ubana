set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} ${ROOTLIB} )
art_make( BASENAME_ONLY
   	      LIBRARY_NAME     HitCosmicTag_Base
	      LIB_LIBRARIES
	      fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
	    )
install_headers()
install_source()
