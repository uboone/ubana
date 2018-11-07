link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make(BASENAME_ONLY
         LIBRARY_NAME uboonecode_uboone_NCElastic 
         LIB_LIBRARIES 	  		   
           uboonecode_uboone_XGBoost
	         ubcore_Geometry
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation 
           ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
           ubana_AnalysisTree_MCTruth_BackTrackerTruth_tool
           ubana_AnalysisTree_MCTruth_MCTruthBase
           nusimdata_SimulationBase
		       art_Persistency_Common canvas
		       art_Persistency_Provenance canvas
		       art_Utilities canvas
		       ${UBOONECODE_LIB}
		       ${LARDATA_LIB}
		       ${ART_FRAMEWORK_CORE}
		       ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_BASIC}
		       ${ART_FRAMEWORK_SERVICES_REGISTRY}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		       ${MF_MESSAGELOGGER}
		       ${MF_UTILITIES}
		       ${FHICLCPP}
		       cetlib cetlib_except
		       ${ROOT_GEOM}
		       ${ROOT_XMLIO}
		       ${ROOT_GDML}
		       ${ROOT_BASIC_LIB_LIST}
 
         MODULE_LIBRARIES 	  		   
           uboonecode_uboone_XGBoost
           uboonecode_uboone_NCElasticAlgorithms
	         ubcore_Geometry
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation 
           ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
           ubana_AnalysisTree_MCTruth_BackTrackerTruth_tool
           ubana_AnalysisTree_MCTruth_MCTruthBase
           nusimdata_SimulationBase
		       art_Persistency_Common canvas
		       art_Persistency_Provenance canvas
		       art_Utilities canvas
		       ${UBOONECODE_LIB}
		       ${LARDATA_LIB}
		       ${ART_FRAMEWORK_CORE}
		       ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_BASIC}
		       ${ART_FRAMEWORK_SERVICES_REGISTRY}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		       ${MF_MESSAGELOGGER}
		       ${FHICLCPP}
		       cetlib cetlib_except
		       ${ROOT_GEOM}
		       ${ROOT_XMLIO}
		       ${ROOT_GDML}
		       ${ROOT_BASIC_LIB_LIST}
		    )

install_headers()
install_fhicl()
install_source()

add_subdirectory(job)
add_subdirectory(Algorithms)
