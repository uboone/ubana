# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

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
  		          ${ART_FRAMEWORK_CORE}
   		          ${ART_FRAMEWORK_PRINCIPAL}

                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
   		          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
     		          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			  art_Persistency_Common canvas
			  art_Persistency_Provenance canvas
			  art_Utilities canvas

                          ${FHICLCPP}
                          cetlib cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
