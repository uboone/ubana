if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage -Wno-undefined-var-template )
endif()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )


art_make( 
          EXCLUDE
            ConsolidatedOutputAnalyserAllTracks_module.cc
          # Needed for standalone c++ helper files
          LIB_LIBRARIES 
            PID
            BackTrack
          # Needed for the main module          
          MODULE_LIBRARIES
            PID
            BackTrack
                        larpandora_LArPandoraInterface
                        larcorealg_Geometry
                        larreco_RecoAlg
                        larcore_Geometry_Geometry_service
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardataobj_AnalysisBase
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			            ${ART_FRAMEWORK_PRINCIPAL}
			            ${ART_FRAMEWORK_SERVICES_REGISTRY}
			            ${ART_FRAMEWORK_SERVICES_BASIC}
			            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			            art_Framework_Services_Registry
			            art_Persistency_Common
			            art_Persistency_Provenance
			            art_Utilities
                        canvas
			            ${MF_MESSAGELOGGER}
	                    ${FHICLCPP}
			            cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
                        NuCCHelpers
          )

art_make_library(
   LIBRARY_NAME
       PID
   SOURCE
       PID.cxx
   LIBRARIES
       lardataobj_RawData
       larcorealg_Geometry
       larcore_Geometry_Geometry_service
       art_Persistency_Common canvas
       art_Persistency_Provenance canvas
       art_Utilities canvas
       lardataobj_RecoBase
       lardataobj_AnalysisBase
       lardata_Utilities
       larpandora_LarPandoraEventBuilding
       LArPandoraContent
       larpandora_LArPandoraInterface
       nusimdata_SimulationBase
       LArPandoraContent
       cetlib cetlib_except
       ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}
       ${ROOT_EG}
)

art_make_library(
   LIBRARY_NAME
       BackTrack
   SOURCE
       BackTrackerTruthMatch.cxx
   LIBRARIES
       lardataobj_RawData
       larcorealg_Geometry
       larcore_Geometry_Geometry_service
       ${ART_FRAMEWORK_CORE}
       ${ART_FRAMEWORK_IO_SOURCES}
       ${ART_FRAMEWORK_PRINCIPAL}
       art_Persistency_Common canvas
       art_Persistency_Provenance canvas
       art_Utilities canvas
       lardataobj_RecoBase
       lardataobj_AnalysisBase
       lardata_Utilities
       larpandora_LarPandoraEventBuilding
       LArPandoraContent
       larpandora_LArPandoraInterface
       nusimdata_SimulationBase
       LArPandoraContent
       ${ART_FRAMEWORK_SERVICES_REGISTRY}
       ${MF_MESSAGELOGGER}
       ${FHICLCPP}
       ${ROOT_BASIC_LIB_LIST}
       ${ROOT_GEOM}
       ${Boost_FILESYSTEM_LIBRARY}
       ${Boost_SYSTEM_LIBRARY}
       cetlib cetlib_except
)



add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()