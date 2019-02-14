cet_enable_asserts()
include_directories ( $ENV{GENIE_INC}/GENIE )
file(GLOB INTERFACEFILES *Interface* )
simple_plugin( TreeReader
               "source"
	       TreeReader
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
               BASENAME_ONLY
             )

art_make_library(LIBRARY_NAME TreeInterface 
                 SOURCE	      ${INTERFACEFILES}
                 LIBRARIES
		       nusimdata_SimulationBase
#GENIE STUFF########################################
			xml2
			${LOG4CPP}
			${GALGORITHM}
			${GBASE}
			${GEVGCORE}
			${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
			${GINTERACTION}
			${GMESSENGER}
			${GHEP}
####################################################
                      ${CETLIB}
		      ${FHICLCPP}
                      #${MF_MESSAGELOGGER}
		      ${CLHEP}
                      ${ROOT_BASIC_LIB_LIST}
                      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make( EXCLUDE TreeReader_source.cc ${INTERFACEFILES}
	  BASENAME_ONLY
          LIBRARY_NAME TreeReader
          LIB_LIBRARIES	nusimdata_SimulationBase
	  		TreeInterface
			larcoreobj_SummaryData
	  		${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		      art_Persistency_Common canvas_Persistency_Common canvas_Persistency_Common
                      art_Persistency_Provenance canvas_Persistency_Provenance canvas_Persistency_Provenance
                      art_Utilities canvas_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}

          MODULE_LIBRARIES TreeReader
	  		TreeInterface
                        nusimdata_SimulationBase
			larcoreobj_SummaryData
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#GENIE STUFF########################################
			xml2
			${LOG4CPP}
			${GALGORITHM}
			${GBASE}
			${GEVGCORE}
			${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
			${GINTERACTION}
			${GMESSENGER}
			${GHEP}
####################################################
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )                                                           

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)

