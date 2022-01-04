art_make( 
          # Needed for standalone c++ helper files
          LIBRARY_NAME
                        SingleTrackHelpers
          LIB_LIBRARIES 
                        larpandora_LArPandoraInterface
			lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_RecoBase
                        lardataobj_AnalysisBase
                        lardata_Utilities
                        nusimdata::SimulationBase
                        art::Framework_Core
			            art::Framework_Principal
			            art::Framework_Services_Registry
			            art::Framework_Services_Registry
			            art::Persistency_Common
			            art::Persistency_Provenance
			            art::Utilities
                        canvas::canvas
			            messagefacility::MF_MessageLogger
	                    fhiclcpp::fhiclcpp
			            cetlib::cetlib
                                    cetlib_except::cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
          )
          
install_headers()
install_source()
install_fhicl()
