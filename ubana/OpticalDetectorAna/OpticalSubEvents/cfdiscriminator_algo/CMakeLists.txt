
art_make_library( LIBRARY_NAME uboone_cfdiscriminator
		  SOURCE CFDFire.cc CFDiscConfig.cc
		  LIBRARIES ${ART_FRAMEWORK_CORE}

		)

install_headers()
install_source()
