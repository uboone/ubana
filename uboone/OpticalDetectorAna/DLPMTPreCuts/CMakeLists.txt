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
	       ubooneobj
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas_Persistency_Common
               art_Persistency_Provenance canvas_Persistency_Provenance
               art_Utilities canvas_Utilities
	       #${DLPMTPRECUTALGO_LIBBASE}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

install_fhicl()
install_headers()
