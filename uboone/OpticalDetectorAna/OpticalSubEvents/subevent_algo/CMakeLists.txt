
art_make_library( LIBRARY_NAME uboone_subeventalgo
		  SOURCE SubEventModule.cc SubEventModConfig.cc
		         scintresponse.cc WaveformData.cc
		  LIBRARIES uboone_cfdiscriminator
		  	    ubooneobj_Optical
		            ${ART_FRAMEWORK_CORE}
                  	    ${ROOT_BASIC_LIB_LIST}
		)							   

install_headers()
install_source()
