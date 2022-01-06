cet_enable_asserts()
find_package(dk2nugenie REQUIRED)
file(GLOB INTERFACEFILES *Interface* )
simple_plugin( TreeReader
               "source"
	       TreeReader
               BASENAME_ONLY
             )

art_make_library(LIBRARY_NAME TreeInterface 
                 SOURCE	      ${INTERFACEFILES}
                 LIBRARIES
		 nusimdata::SimulationBase
#GENIE STUFF########################################
			xml2
			 log4cpp::log4cpp
                 	 GENIE::GFwEG
                 	 GENIE::GFwAlg
                 	 GENIE::GFwMsg
                 	 GENIE::GFwGHEP
                 	 GENIE::GTlFlx
                 	 dk2nu::Tree
                 	 dk2nu::Genie
			log4cpp::log4cpp
####################################################
		      fhiclcpp::fhiclcpp
		      cetlib::cetlib
                      cetlib_except::cetlib_except
		      #messagefacility::MF_MessageLogger
		      CLHEP::CLHEP
                      ROOT::Tree
	        ) 

art_make( EXCLUDE TreeReader_source.cc ${INTERFACEFILES}
	  BASENAME_ONLY
          LIBRARY_NAME TreeReader
          LIB_LIBRARIES	nusimdata::SimulationBase
	  		TreeInterface
			larcoreobj_SummaryData
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
          MODULE_LIBRARIES 
	  		TreeReader
	  		TreeInterface
                        nusimdata::SimulationBase
			larcoreobj_SummaryData
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
#GENIE STUFF########################################
			xml2
			GENIE::GFwEG
                 	GENIE::GFwAlg
                	GENIE::GFwMsg
               		GENIE::GFwGHEP
               		GENIE::GTlFlx
			log4cpp::log4cpp
               		dk2nu::Tree
			dk2nu::Genie
			log4cpp::log4cpp
####################################################
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
		        cetlib::cetlib
                        cetlib_except::cetlib_except
        )                                                           

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
