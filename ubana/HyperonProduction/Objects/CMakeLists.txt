art_make( 

      LIB_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}

      DICT_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}

      MODULE_LIBRARIES
      larcorealg_Geometry
      ${ROOT_BASIC_LIB_LIST}

      )

install_headers()
install_fhicl()
install_source()
