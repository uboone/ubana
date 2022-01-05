cet_enable_asserts()
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{DK2NUDATA_INC} )
link_directories( $ENV{DK2NUDATA_LIB})
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
                 	 dk2nuTree
                 	 #${DK2NUDATA}
                 	 ${DK2NUGENIE}
			log4cpp::log4cpp
####################################################
		      fhiclcpp::fhiclcpp
		      cetlib::cetlib cetlib_except
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
		      art::Persistency_Common canvas 
                      art::Persistency_Provenance canvas 
                      art::Utilities canvas
			messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        
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
               		dk2nuTree
			${DK2NUGENIE}
			log4cpp::log4cpp
####################################################
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
        )                                                           

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
