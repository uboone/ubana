# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
          LIBRARY_NAME     uboonecode_uboone_CosmicTagging
	  LIB_LIBRARIES
		   ubcore::ubana_BasicTool_GeoAlgo
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
                   ubana_CosmicTagging_CosmicTaggingAlgorithms
                   #Simulation
		   #Utilities
		   #TimeService_service
		   #MCBase
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larpandora_LArPandoraInterface
                   ${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   art::Framework_Core
		   art::Framework_Principal

		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger

		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml


          MODULE_LIBRARIES
                   larsim_MCCheater_ParticleInventoryService_service
                   larsim_MCCheater_BackTrackerService_service
		   larpandora_LArPandoraInterface

		   ubcore::ubana_BasicTool_GeoAlgo
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms

                   ubana_CosmicTagging_CosmicTaggingAlgorithms

	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larsim_Simulation lardataobj_Simulation
                   nusimdata::SimulationBase
		   ${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   art::Framework_Core
		   art::Framework_Principal

		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger

		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml

		)

# install_headers()
install_fhicl()
install_source()

add_subdirectory(Algo)
add_subdirectory(job)
