
cet_enable_asserts()

set( hitfinder_tool_lib_list
                        larcorealg_Geometry
                        lardataobj_RecoBase
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        Boost::filesystem
                        CLHEP::CLHEP
    )

art_make( BASENAME_ONLY
          LIBRARY_NAME     uboonecode_uboone_UBXSec_Algorithms
          LIB_LIBRARIES
                   ubobj_UBXSec
                   ubcore::ubana_BasicTool_GeoAlgo
                   ubreco::uboone_UBFlashFinder
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
                   ubraw_RawData
                   ubevt::SignalShapingServiceMicroBooNE_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
                   larsim_MCCheater_BackTrackerService_service
		   larsim_MCCheater_ParticleInventoryService_service
                   larsim_Simulation
                   #Simulation
		   #Utilities
		   #TimeService_service
		   #MCBase
                   lardataobj_RawData
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardataobj_Simulation
                   lardata_Utilities
                   larpandora_LArPandoraInterface
                   nusimdata::SimulationBase
                   PANDORASDK
                   PANDORAMONITORING
                   #larpandora_LArPandoraContent
                   ${UBOONECODE_LIB}
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
                   art::Persistency_Common canvas::canvas
                   art::Persistency_Provenance canvas::canvas
                   art::Utilities canvas::canvas
                   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib
                   cetlib_except::cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml

          TOOL_LIBRARIES
                   ${hitfinder_tool_lib_list}

          MODULE_LIBRARIES
                   ubobj_UBXSec
                   uboonecode_uboone_UBXSec_Algorithms
                   larsim_MCCheater_BackTrackerService_service
		   larpandora_LArPandoraInterface
                   #                   larpandora_LArPandoraContent
		   ubcore::ubana_BasicTool_GeoAlgo
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
                   ubreco::uboone_UBFlashFinder
                   ubraw_RawData
                   ubevt::SignalShapingServiceMicroBooNE_service
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
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib
                   cetlib_except::cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
		)

install_headers()
install_source()

