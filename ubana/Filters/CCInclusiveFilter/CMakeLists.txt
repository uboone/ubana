art_make( 
          LIBRARY_NAME uboonecode_uboone_Filters
          # Needed for standalone c++ helper files
          LIB_LIBRARIES 

          # Needed for the main module          
          MODULE_LIBRARIES
                    larpandora_LArPandoraInterface
                    larcorealg_Geometry
                    larreco_RecoAlg
                    larcore_Geometry_Geometry_service
                    larsim_Simulation lardataobj_Simulation
                    lardataobj_RawData
                    lardataobj_RecoBase
                    lardataobj_AnalysisBase
                    lardata_Utilities
                    nusimdata::SimulationBase
                    larsim_EventWeight_Base
                    art::Framework_Core
                    art::Framework_Principal
                    art::Framework_Services_Registry
                    art::Framework_Services_Registry
                    art::Persistency_Common
                    art::Persistency_Provenance
                    art::Utilities
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    cetlib_except::cetlib_except
                    NuCCHelpers
          )

add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
