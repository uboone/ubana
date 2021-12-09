
art_make_library( LIBRARY_NAME uboone_subeventalgo
		  SOURCE SubEventModule.cc SubEventModConfig.cc
		         scintresponse.cc WaveformData.cc
		  LIBRARIES uboone_cfdiscriminator
		  	    ubobj_Optical
		            art::Framework_Core

		)

install_headers()
install_source()
