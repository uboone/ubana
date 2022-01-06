# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
	LIBRARY_NAME TPCNeutrinoIDAlgs
	LIB_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
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


	MODULE_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
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


         )

install_headers()
install_source()
install_fhicl()
