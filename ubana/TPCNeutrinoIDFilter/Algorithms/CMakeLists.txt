# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
	LIBRARY_NAME TPCNeutrinoIDAlgs
	LIB_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
		     art::Persistency_Common canvas::canvas
		     art::Persistency_Provenance canvas::canvas
		     art::Framework_Core
		     art::Framework_Principal
		     art::Framework_Services_Registry
		     art_root_io::tfile_support ROOT::Core
		     art_root_io::TFileService_service
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     art::Utilities canvas::canvas
		     messagefacility::MF_MessageLogger

		     fhiclcpp::fhiclcpp
		     cetlib::cetlib cetlib_except


	MODULE_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
		     art::Persistency_Common canvas::canvas
		     art::Persistency_Provenance canvas::canvas
		     art::Framework_Core
		     art::Framework_Principal
		     art::Framework_Services_Registry
		     art_root_io::tfile_support ROOT::Core
		     art_root_io::TFileService_service
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     art::Utilities canvas::canvas
		     messagefacility::MF_MessageLogger

		     fhiclcpp::fhiclcpp
		     cetlib::cetlib cetlib_except


         )

install_headers()
install_source()
install_fhicl()
