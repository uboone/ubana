# Either set these here or set environment variables

##include_directories ( $ENV{SWTRIGGER_INCDIR} )
cet_find_library( SWTRIG_LIBBASE NAMES SWTriggerBase PATHS $ENV{SWTRIGGER_LIBDIR} )
cet_find_library( SWTRIG_LIBFEMU NAMES FEMBeamTrigger PATHS $ENV{SWTRIGGER_LIBDIR} )

add_subdirectory(fcl)

art_make( LIBRARY_NAME SPEcalibrationCode
          EXCLUDE SPEcalibration_module.cc FEMemulator_module.cc SWTrigger_module.cc CalibrationFit.C
          lardataobj_RawData
          ubcore::UBOpReadoutMap_service
          LIB_LIBRARIES art::Framework_Core
               art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common canvas::canvas
               art::Persistency_Provenance canvas::canvas
               art::Utilities canvas::canvas

               BASENAME_ONLY
        )


simple_plugin( SPEcalibration "module"
               SPEcalibrationCode
               OpticalSubEvents
               ubcore::UBOpReadoutMap_service
               lardataobj_RawData
               art::Framework_Core
               art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common canvas::canvas
               art::Persistency_Provenance canvas::canvas
               art::Utilities canvas::canvas
               ROOT::Tree
               BASENAME_ONLY
             )

simple_plugin( FEMemulator "module"
               SPEcalibrationCode
               OpticalSubEvents
               ubcore::UBOpReadoutMap_service
               lardataobj_RawData
               art::Framework_Core
               art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common canvas::canvas
               art::Persistency_Provenance canvas::canvas
               art::Utilities canvas::canvas
               ROOT::Tree
               BASENAME_ONLY
             )

simple_plugin( SWTrigger "module"
               ubcore::UBOpReadoutMap_service
               ubobj_Trigger
               lardataobj_RawData
               lardataobj_RecoBase
               art::Framework_Core
               art::Framework_IO_Sources
               art::Framework_Principal
               art::Framework_Services_Registry
               art_root_io::tfile_support ROOT::Core
               art_root_io::TFileService_service
               art::Persistency_Common
               art::Persistency_Provenance
               art::Utilities
               SWTRIG_LIBFEMU
               SWTRIG_LIBBASE
               ROOT::Tree
               BASENAME_ONLY
             )

install_fhicl()
install_headers()
install_source()
