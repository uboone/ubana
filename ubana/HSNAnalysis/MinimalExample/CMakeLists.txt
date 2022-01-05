# basic source code CMakeLists.txt

art_make(BASENAME_ONLY
	LIB_LIBRARIES
		lardataobj_RecoBase
		lardataobj_AnalysisBase
		larcorealg_Geometry
    larcore_Geometry_Geometry_service
		lardata_Utilities
		nusimdata::SimulationBase
		messagefacility::MF_MessageLogger
		
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except
		CLHEP::CLHEP

		${G4_LIB_LIST}
	MODULE_LIBRARIES
		larreco_RecoAlg
		larreco_RecoAlg_Cluster3DAlgs
		larsim_Simulation
		lardataobj_RawData
		lardata_RecoObjects
		lardataobj_RecoBase
		lardataobj_MCBase
		lardataobj_AnalysisBase
		larcorealg_Geometry
    larcore_Geometry_Geometry_service
		lardata_Utilities
		nusimdata::SimulationBase
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except
		CLHEP::CLHEP

		${G4_LIB_LIST}
	)

install_headers()
install_source()
install_fhicl()
