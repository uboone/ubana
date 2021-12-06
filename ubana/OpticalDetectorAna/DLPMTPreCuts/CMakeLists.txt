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
	       ubobj_Optical
	       art::Framework_Core
	       art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common canvas::canvas
               art::Persistency_Provenance canvas::canvas
               art::Utilities canvas::canvas
	       #${DLPMTPRECUTALGO_LIBBASE}

               BASENAME_ONLY
             )

install_fhicl()
install_headers()
install_source()

