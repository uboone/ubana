# Either set these here or set environment variables

##include_directories ( $ENV{SWTRIGGER_INCDIR} )
cet_find_library( SWTRIG_LIBBASE NAMES SWTriggerBase PATHS $ENV{SWTRIGGER_LIBDIR} )
cet_find_library( SWTRIG_LIBFEMU NAMES FEMBeamTrigger PATHS $ENV{SWTRIGGER_LIBDIR} )

add_subdirectory(fcl)

art_make( LIBRARY_NAME SPEcalibrationCode
          EXCLUDE SPEcalibration_module.cc FEMemulator_module.cc SWTrigger_module.cc CalibrationFit.C
          lardataobj_RawData
          ubcore::UBOpReadoutMap_service
          LIB_LIBRARIES ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas

               BASENAME_ONLY
        )


simple_plugin( SPEcalibration "module"
               SPEcalibrationCode
               OpticalSubEvents
               ubcore::UBOpReadoutMap_service
               lardataobj_RawData
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
               ROOT::Tree
               BASENAME_ONLY
             )

simple_plugin( FEMemulator "module"
               SPEcalibrationCode
               OpticalSubEvents
               ubcore::UBOpReadoutMap_service
               lardataobj_RawData
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
               ROOT::Tree
               BASENAME_ONLY
             )

simple_plugin( SWTrigger "module"
               ubcore::UBOpReadoutMap_service
               ubobj_Trigger
               lardataobj_RawData
               lardataobj_RecoBase
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
               ${ART_ROOT_IO_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
               SWTRIG_LIBFEMU
               SWTRIG_LIBBASE
               ROOT::Tree
               BASENAME_ONLY
             )

install_fhicl()
install_headers()
install_source()
