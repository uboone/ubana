# find_ups_product is meant to be called from the top CMakeLists.txt file

#add_subdirectory( cfdiscriminator_algo )
#add_subdirectory( subevent_algo )

#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -D__BUILD_ROOT_DICT__" )

art_make( LIBRARY_NAME OpticalSubEvents
	  SUBDIRS cfdiscriminator_algo subevent_algo
	  EXCLUDE SubEventBuilder_module.cc
	  LIB_LIBRARIES
	       ubobj_Optical
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

	)


simple_plugin( SubEventBuilder "module"
	       OpticalSubEvents
	       ubobj_Optical
	       lardataobj_RecoBase
	       larcorealg_Geometry
	       larcore_Geometry_Geometry_service
	       lardata_Utilities
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

               BASENAME_ONLY
             )


install_headers()
install_source()

install_fhicl()
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
