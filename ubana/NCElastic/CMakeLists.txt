
art_make(BASENAME_ONLY
         LIBRARY_NAME uboonecode_uboone_NCElastic
         LIB_LIBRARIES
           uboonecode_uboone_XGBoost
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation
           nusimdata::SimulationBase
		       art::Persistency_Common canvas::canvas
		       art::Persistency_Provenance canvas::canvas
		       art::Utilities canvas::canvas
		       art::Framework_Core
		       art::Framework_Principal
		       art::Framework_Services_Registry
		       art_root_io::tfile_support ROOT::Core
		       art_root_io::TFileService_service
		       messagefacility::MF_MessageLogger
		       fhiclcpp::fhiclcpp
		       cetlib::cetlib
                       cetlib_except::cetlib_except
		       ROOT::Geom
		       ROOT::XMLIO
		       ROOT::Gdml

         MODULE_LIBRARIES
           uboonecode_uboone_XGBoost
           uboonecode_uboone_NCElasticAlgorithms
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
		       lardataobj_RecoBase
           lardataobj_AnalysisBase
           lardata_Utilities
			     larreco_RecoAlg
			     larpandora_LArPandoraInterface
           larsim_Simulation lardataobj_Simulation
           nusimdata::SimulationBase
		       art::Persistency_Common canvas::canvas
		       art::Persistency_Provenance canvas::canvas
		       art::Utilities canvas::canvas
		       ${UBOONECODE_LIB}
		       art::Framework_Core
		       art::Framework_Principal
		       art::Framework_Services_Registry
		       art_root_io::tfile_support ROOT::Core
		       art_root_io::TFileService_service
		       messagefacility::MF_MessageLogger
		       fhiclcpp::fhiclcpp
		       cetlib::cetlib
                       cetlib_except::cetlib_except
		       ROOT::Geom
		       ROOT::XMLIO
		       ROOT::Gdml

		    )

install_headers()
install_fhicl()
install_source()

add_subdirectory(job)
add_subdirectory(Algorithms)
