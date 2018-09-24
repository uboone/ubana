link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make(BASENAME_ONLY
         LIBRARY_NAME uboonecode_uboone_NCElastic 
         LIB_LIBRARIES 	  		   
           uboonecode_uboone_XGBoost
	  	     larcore_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation 
           nusimdata_SimulationBase
		       art_Persistency_Common canvas_Persistency_Common
		       art_Persistency_Provenance canvas_Persistency_Provenance
		       art_Utilities canvas_Utilities
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
			     ${CETLIB}
		       ${ROOT_GEOM}
		       ${ROOT_XMLIO}
		       ${ROOT_GDML}
		       ${ROOT_BASIC_LIB_LIST}
 
         MODULE_LIBRARIES 	  		   
           uboonecode_uboone_XGBoost
           uboonecode_uboone_NCElasticAlgorithms
	  	     larcore_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation 
           nusimdata_SimulationBase
		       art_Persistency_Common canvas_Persistency_Common
		       art_Persistency_Provenance canvas_Persistency_Provenance
		       art_Utilities canvas_Utilities
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
			     ${CETLIB}
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
