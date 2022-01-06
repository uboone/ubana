# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
  LIBRARY_NAME uboonecode_uboone_Filters
	LIB_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
		     art::Persistency_Common
		     art::Persistency_Provenance
		     art::Framework_Core
		     art::Framework_Principal
		     art::Framework_Services_Registry
		     art::Utilities
                     canvas::canvas
		     messagefacility::MF_MessageLogger
		     fhiclcpp::fhiclcpp
		     cetlib::cetlib
                     ROOT::EG
	MODULE_LIBRARIES
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcorealg_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
         uboonecode_uboone_UBXSec_Algorithms
         ubobj_UBXSec
		     art::Persistency_Common
		     art::Persistency_Provenance
		     art::Framework_Core
		     art::Framework_Principal
		     art::Framework_Services_Registry
		     art::Utilities
                     canvas::canvas
		     messagefacility::MF_MessageLogger
		     fhiclcpp::fhiclcpp
		     cetlib::cetlib
                     ROOT::EG
                     CoarseNuCCHelpers
         )

add_subdirectory(helpers)
add_subdirectory(Muon)

install_headers()
install_source()
install_fhicl()
