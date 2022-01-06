
art_make( BASENAME_ONLY
	  LIBRARY_NAME  ubana_OpticalDetectorAna
	  LIB_LIBRARIES
			larana_OpticalDetector
			larana_OpticalDetector_OpHitFinder


	  MODULE_LIBRARIES larcorealg_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                ubcore_Geometry
			ubcore::UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
			larana_OpticalDetector
                        lardataobj_RawData
			lardataobj_RecoBase
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			nusimdata::SimulationBase
			larsim_PhotonPropagation
			ubana_OpticalDetectorAna
                        art::Framework_Core
			art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml

			nusimdata::SimulationBase
        )

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
add_subdirectory(OpticalSubEvents)
add_subdirectory(DLPMTPreCuts)
