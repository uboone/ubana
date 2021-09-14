art_make( 

      LIBRARY_NAME     uboonecode_uboone_EventConnectedness

      LIB_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}
      lardataobj_RecoBase
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
      cetlib cetlib_except

      )

install_headers()
install_source()
