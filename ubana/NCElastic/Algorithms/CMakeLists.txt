
art_make(BASENAME_ONLY
         LIBRARY_NAME uboonecode_uboone_NCElasticAlgorithms
         LIB_LIBRARIES
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
	         ubcore_Geometry
           ubcore::UBOpReadoutMap_service
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
	  	     larcorealg_Geometry
		       larcore_Geometry_Geometry_service
	         uboone_Geometry
           ubcore::UBOpReadoutMap_service
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
		       ${LARDATA_LIB}
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
