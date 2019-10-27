include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
          LIBRARY_NAME     uboonecode_uboone_UBXSec
          LIB_LIBRARIES
                   uboonecode_uboone_UBXSec_Algorithms
                   BackTrack
                   BrokenTrack
                   RecoTruth
                   Topology
                   PID
                   ubobj_UBXSec
                   ubobj_Trigger
                   HitCosmicTag_Base 
                   HitCosmicTag_Algorithms
                   ubana_BasicTool_GeoAlgo
                   uboone_UBFlashFinder
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
                   ubraw_RawData
                   SignalShapingServiceMicroBooNE_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
                   larsim_MCCheater_BackTrackerService_service
                   larreco_RecoAlg
                   #Simulation
		   #Utilities
		   #TimeService_service
		   #MCBase
                   lardataobj_RawData
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larpandora_LArPandoraInterface
                   nusimdata_SimulationBase
                   ${PANDORASDK}
                   ${PANDORAMONITORING}
                   #larpandora_LArPandoraContent
                   ${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                   art_Persistency_Common canvas
                   art_Persistency_Provenance canvas
                   art_Utilities canvas
                   ${ART_PERSISTENCY_COMMON}
		   ${ART_PERSISTENCY_PROVENANCE}
		   ${ART_UTILITIES}
		   ${MF_MESSAGELOGGER}
		   
		   ${FHICLCPP}
		   cetlib cetlib_except
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}
                   ${ROOT_EG}

          MODULE_LIBRARIES
                   uboonecode_uboone_UBXSec_Algorithms 
                   BackTrack
                   BrokenTrack
                   RecoTruth
                   Topology
                   PID
                   ubobj_UBXSec
		   ubobj_Trigger
                   HitCosmicTag_Base
                   HitCosmicTag_Algorithms
                   #uboonecode_uboone_UBXSec
                   larsim_MCCheater_BackTrackerService_service
		   larpandora_LArPandoraInterface
                   #                   larpandora_LArPandoraContent

		   ubana_BasicTool_GeoAlgo
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
                   uboone_UBFlashFinder

                   ubraw_RawData
                   larreco_RecoAlg
                   SignalShapingServiceMicroBooNE_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larsim_Simulation lardataobj_Simulation 
                   nusimdata_SimulationBase
		   ${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		   ${ART_PERSISTENCY_COMMON}
		   ${ART_PERSISTENCY_PROVENANCE}
		   ${ART_UTILITIES}
		   ${MF_MESSAGELOGGER}
		   
		   ${FHICLCPP}
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

art_make_library(
   LIBRARY_NAME
       BrokenTrack
   SOURCE
       BrokenTrack.cxx
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

art_make_library(
   LIBRARY_NAME
       RecoTruth
   SOURCE
       RecoTruthMCParticle.cxx
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
       ${larevt_SpaceChargeServices}
       ${ART_FRAMEWORK_SERVICES_REGISTRY}
       ${MF_MESSAGELOGGER}
       ${FHICLCPP}
       ${ROOT_BASIC_LIB_LIST}
       ${ROOT_GEOM}
       ${Boost_FILESYSTEM_LIBRARY}
       ${Boost_SYSTEM_LIBRARY}
       cetlib cetlib_except

)

art_make_library(
   LIBRARY_NAME
       Topology
   SOURCE
       Topology.cxx
   LIBRARIES
       cetlib cetlib_except

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

)

install_headers()
install_fhicl()
install_source()

