# Either set these here or set environment variables

#include_directories ( $ENV{DLPMTPRECUTALGO_INCDIR} )
#find_library( DLPMTPRECUTALGO_LIBBASE NAMES LEEPreCutAlgo PATHS $ENV{DLPMTPRECUTALGO_LIBDIR} )

add_subdirectory(fcl)

art_make_library( LIBRARY_NAME uboone_LEEPreCutAlgo
		  SOURCE LEEPreCutAlgo.cxx
		)


simple_plugin( DLPMTPreCuts "module"
	       uboone_LEEPreCutAlgo
	       lardataobj_RecoBase
               larcore_Geometry_Geometry_service
               larcorealg_Geometry
	       ubobj_Optical
	       ubevt_Database
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
	       #${DLPMTPRECUTALGO_LIBBASE}

               BASENAME_ONLY
             )

install_fhicl()
install_headers()
install_source()

