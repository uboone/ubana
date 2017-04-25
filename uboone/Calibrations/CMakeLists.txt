
add_subdirectory(Optical)

art_make( BASENAME_ONLY
	  EXCLUDE Lifetime.C
	  LIBRARY_NAME  uboonecode_uboone_Calibrations
	  LIB_LIBRARIES
                        lardata_Utilities_LArFFT_service
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
      ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
			uboonecode_uboone_Calibrations
			lardataobj_RawData
			CalData
			larcore_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_BackTracker_service
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
                        larreco_Calorimetry
			lardata_Utilities
			larevt_Filters
                        lardata_Utilities_LArFFT_service
			SignalShapingServiceMicroBooNE_service
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_MINUIT}
			${ROOT_MINUIT2}	
        )

FILE(GLOB data_files calibrationList_*.dat )

install_headers()
install_fhicl()
install_source()
install_fw( LIST ${data_files} )
