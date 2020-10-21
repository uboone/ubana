art_make_library( LIBRARY_NAME uboone_hyperontools
                  SOURCE SimParticle.cc
                  LIBRARIES ${ART_FRAMEWORK_CORE}
                            ${ROOT_BASIC_LIB_LIST}
                )

install_headers()
install_source()


