# basic source code CMakeLists.txt

art_make(
	MODULE_LIBRARIES
			ubana_PatternFilter_PMAlgs
			larcorealg_Geometry
			lardataalg_DetectorInfo
			lardataobj_RecoBase
                        art::Framework_Core
			art_root_io::art_root_io
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
                        ${BOOST_LIB}
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except

)

install_headers()
install_source()
install_fhicl()

# Libraries I didn't need
#			ub_RawDigitAdderAna
#
#			art::Framework_IO_Sources
#                        art::Persistency_Provenance canvas::canvas
#                        messagefacility::MF_MessageLogger
#
#
