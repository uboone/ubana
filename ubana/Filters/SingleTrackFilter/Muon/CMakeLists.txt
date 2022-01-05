include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} ${ROOTLIB} )

art_make( BASENAME_ONLY
          LIBRARY_NAME     uboonecode_uboone_UBXSec
          LIB_LIBRARIES
                   uboonecode_uboone_UBXSec_Algorithms
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
                   nusimdata::SimulationBase
                   PANDORASDK
                   PANDORAMONITORING
                   #larpandora_LArPandoraContent
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   
                   art::Persistency_Common canvas
                   art::Persistency_Provenance canvas
                   art::Utilities canvas
                   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
		   ${ROOT_BASIC_LIB_LIST}
                   ROOT::EG

          MODULE_LIBRARIES
                   uboonecode_uboone_UBXSec_Algorithms
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
                   nusimdata::SimulationBase
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
		   ${ROOT_BASIC_LIB_LIST}
                   ROOT::EG
                   )


install_headers()
install_fhicl()
install_source()

