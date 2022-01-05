if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage -Wno-undefined-var-template )
endif()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make( 
          EXCLUDE
            ConsolidatedOutputAnalyserAllTracks_module.cc
          # Needed for standalone c++ helper files
          LIB_LIBRARIES 
            PID
            BackTrack
          # Needed for the main module          
          MODULE_LIBRARIES
            PID
            BackTrack
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
			            cetlib::cetlib cetlib_except
                        
                        NuCCHelpers
          )

art_make_library(
   LIBRARY_NAME
       PID
   SOURCE
       PID.cxx
   LIBRARIES
       lardataobj_RawData
       larcorealg_Geometry
       larcore_Geometry_Geometry_service
       art::Persistency_Common canvas
       art::Persistency_Provenance canvas
       art::Utilities canvas
       lardataobj_RecoBase
       lardataobj_AnalysisBase
       lardata_Utilities
       larpandora_LArPandoraEventBuilding
       LArPandoraContent
       larpandora_LArPandoraInterface
       nusimdata::SimulationBase
       LArPandoraContent
       cetlib::cetlib cetlib_except
       ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
       ROOT::EG
)

art_make_library(
   LIBRARY_NAME
       BackTrack
   SOURCE
       BackTrackerTruthMatch.cxx
   LIBRARIES
       lardataobj_RawData
       larcorealg_Geometry
       larcore_Geometry_Geometry_service
       art::Framework_Core
       art::Framework_IO_Sources
       art::Framework_Principal
       art::Persistency_Common canvas
       art::Persistency_Provenance canvas
       art::Utilities canvas
       lardataobj_RecoBase
       lardataobj_AnalysisBase
       lardata_Utilities
       larpandora_LArPandoraEventBuilding
       LArPandoraContent
       larpandora_LArPandoraInterface
       nusimdata::SimulationBase
       LArPandoraContent
       art::Framework_Services_Registry
       messagefacility::MF_MessageLogger
       fhiclcpp::fhiclcpp
       ROOT::Geom
       Boost::filesystem
       
       cetlib::cetlib cetlib_except
)



add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
