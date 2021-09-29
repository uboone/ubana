# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
    	    LIBRARY_NAME     HitCosmicTag_Algorithms
	        LIB_LIBRARIES
                           HitCosmicTag_Base
  		                     ${ART_FRAMEWORK_CORE}
   		                     ${ART_FRAMEWORK_PRINCIPAL}

                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
   		                     ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
     		                   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
     		                   ${ART_PERSISTENCY_COMMON}
                           ${ART_PERSISTENCY_PROVENANCE}
                           ${ART_UTILITIES}

                           ${FHICLCPP}
                           cetlib cetlib_except
	      )

install_headers()
install_source()
