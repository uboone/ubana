
art_make_library( LIBRARY_NAME uboone_cfdiscriminator
		  SOURCE CFDFire.cc CFDiscConfig.cc
		  LIBRARIES art::Framework_Core

		)

install_headers()
install_source()
