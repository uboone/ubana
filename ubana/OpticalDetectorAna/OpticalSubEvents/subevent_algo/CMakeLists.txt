
art_make_library( LIBRARY_NAME uboone_subeventalgo
		  SOURCE SubEventModule.cc SubEventModConfig.cc
		         scintresponse.cc WaveformData.cc
		  LIBRARIES uboone_cfdiscriminator
		  	    ubobj_Optical
		            ${ART_FRAMEWORK_CORE}

		)

install_headers()
install_source()
