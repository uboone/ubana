
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make(

	  #LIBRARY_NAME  uboonecode_uboone_MyClasses
	  LIB_LIBRARIES

                   ${Boost_SYSTEM_LIBRARY}
		   ${ROOT_BASIC_LIB_LIST}
		   ${CETLIB}

		   art_Framework_Core 
		   art_Framework_Services_Registry 
		   art_Framework_Principal
		   art_Utilities canvas
                   cetlib cetlib_except
		   lardataobj_RawData
		   lardataobj_Simulation
		   larsim_Simulation
                   larsim_MCCheater_BackTrackerService_service
                   larsim_MCCheater_ParticleInventoryService_service
		   larpandora_LArPandoraInterface
		   nusimdata_SimulationBase
		   ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
                   ubraw_RawData


          MODULE_LIBRARIES

                   art_Persistency_Common canvas_Persistency_Common
                   art_Persistency_Provenance canvas_Persistency_Provenance
		   art_Persistency_Common canvas_Persistency_Common
                   art_Utilities canvas_Utilities	
		   larpandora_LArPandoraInterface
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardataobj_RawData
                   lardata_Utilities
                   larsim_Simulation lardataobj_Simulation
                   larsim_MCCheater_BackTrackerService_service
                   larsim_MCCheater_ParticleInventoryService_service
                   nusimdata_SimulationBase
		   ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
                   ubraw_RawData

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
            	   ${MF_UTILITIES}
            	   ${FHICLCPP}
            	   ${CETLIB}
            	   ${ROOT_GEOM}
            	   ${ROOT_XMLIO}
            	   ${ROOT_GDML}
            	   ${ROOT_BASIC_LIB_LIST}

		)

install_headers()
install_fhicl()
install_source()
