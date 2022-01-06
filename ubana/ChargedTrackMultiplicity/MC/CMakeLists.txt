
art_make(
          MODULE_LIBRARIES larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
                        ubana_AnalysisTree_MCTruth_BackTrackerTruth_tool
			ubana_AnalysisTree_MCTruth_MCTruthBase
	                lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			lardataobj_RecoBase
			lardata_RecoObjects
			larreco_RecoAlg
			larpandora_LArPandoraInterface
			ubevt::SignalShapingServiceMicroBooNE_service
                        nusimdata::SimulationBase
                        ubraw_RawData
			ubobj_Trigger
			ubobj_Optical
                        art::Framework_Core
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

        )

# install_headers()
install_fhicl()
install_source()
