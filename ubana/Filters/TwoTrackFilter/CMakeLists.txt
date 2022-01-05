include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make( BASENAME_ONLY
          LIBRARY_NAME     NuMuCC0pi0p
          LIB_LIBRARIES
                   #uboonecode_uboone_UBXSec_Algorithms
                   BackTracker
                   RecoTruth
                   Topology
                   PID3pl
                   FiducialVolume
                   ubobj_UBXSec
                   ubobj_Trigger
                   HitCosmicTag_Base 
                   HitCosmicTag_Algorithms
                   ubana_BasicTool_GeoAlgo
                   uboone_UBFlashFinder
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
                   ubraw_RawData
                   SignalShapingServiceMicroBooNE_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
                   larsim_MCCheater_BackTrackerService_service
                   larreco_RecoAlg
                   #Simulation
		   #Utilities
		   #TimeService_service
		   #MCBase
                   lardataobj_RawData
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larpandora_LArPandoraInterface
                   nusimdata::SimulationBase
                   PANDORASDK
                   PANDORAMONITORING
                   #larpandora_LArPandoraContent
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   
                   art::Persistency_Common canvas
                   art::Persistency_Provenance canvas
                   art::Utilities canvas
                   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
                   ROOT::EG

          MODULE_LIBRARIES
                   #uboonecode_uboone_UBXSec_Algorithms 
                   BackTracker
                   RecoTruth
                   Topology
                   PID3pl
                   FiducialVolume
                   ubobj_UBXSec
		   ubobj_Trigger
                   HitCosmicTag_Base
                   HitCosmicTag_Algorithms
                   #uboonecode_uboone_UBXSec
                   larsim_MCCheater_BackTrackerService_service
		   larpandora_LArPandoraInterface
                   #                   larpandora_LArPandoraContent

		   ubana_BasicTool_GeoAlgo
		   ubreco_SelectionTool_OpT0FinderBase
		   ubreco_SelectionTool_OpT0FinderAlgorithms
                   uboone_UBFlashFinder

                   ubraw_RawData
                   larreco_RecoAlg
                   SignalShapingServiceMicroBooNE_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
	  	   larcorealg_Geometry
		   larcore_Geometry_Geometry_service
		   lardataobj_RecoBase
                   lardataobj_AnalysisBase
                   lardata_Utilities
                   larsim_Simulation lardataobj_Simulation 
                   nusimdata::SimulationBase
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   messagefacility::MF_MessageLogger
		   
		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
                   ROOT::EG
                   )

art_make_library(
   LIBRARY_NAME
       BackTracker
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


art_make_library(
   LIBRARY_NAME
       RecoTruth
   SOURCE
       RecoTruthMCParticle.cxx
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

art_make_library(
   LIBRARY_NAME
       Topology
   SOURCE
       Topology.cxx
   LIBRARIES
       cetlib::cetlib cetlib_except

)

art_make_library(
   LIBRARY_NAME
       PID3pl
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

)

art_make_library(
   LIBRARY_NAME
       FiducialVolume
   SOURCE
       FiducialVolume.cxx
   LIBRARIES
       ubana_BasicTool_GeoAlgo
       lardataobj_RawData
       larcorealg_Geometry
       larcore_Geometry_Geometry_service
       art::Framework_Core
       art::Framework_IO_Sources
       art::Framework_Principal
       art::Persistency_Common canvas
       art::Persistency_Provenance canvas
       art::Utilities canvas
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
       art::Framework_Services_Registry
       messagefacility::MF_MessageLogger
       fhiclcpp::fhiclcpp
       ROOT::Geom
       Boost::filesystem
       
       cetlib::cetlib cetlib_except

)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

