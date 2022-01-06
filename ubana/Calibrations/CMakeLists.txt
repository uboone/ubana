
cet_remove_compiler_flags(CXX -pedantic)

add_subdirectory(Optical)

art_make( BASENAME_ONLY
	  EXCLUDE Lifetime.C
	  LIBRARY_NAME  ubana_Calibrations
	  LIB_LIBRARIES
                        lardata_Utilities_LArFFT_service
			art::Framework_Core
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger

			cetlib::cetlib
                        cetlib_except::cetlib_except


          MODULE_LIBRARIES
			ubana_Calibrations
                        lardataalg_DetectorInfo
			lardataobj_RawData
			ubevt::CalData
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_ParticleInventoryService_service
      larsim_MCCheater_BackTrackerService_service
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			larreco_RecoAlg
                        larreco_Calorimetry
			lardata_Utilities
			larevt_Filters
                        lardata_Utilities_LArFFT_service
			ubevt::SignalShapingServiceMicroBooNE_service
            larevt_CalibrationDBI_Providers
            nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
        )

FILE(GLOB data_files calibrationList_*.dat )

install_headers()
install_fhicl()
install_source()
install_fw( LIST ${data_files} )
