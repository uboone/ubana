cet_build_plugin(
  ACPTTagger art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::UBFlashFinder
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CandidateConsistency art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubana::UBXSec_HitCosmicTag_Algorithms
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubcore::LLBasicTool_GeoAlgo
  ubobj::UBXSec
  ubevt::Utilities_SignalShapingServiceMicroBooNE_service
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CosmicFlashMatch art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubcore::Geometry
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  CosmicTaggerAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  CosmicTrackHitTagger art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubobj::Trigger
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardataobj::RecoBase
)

cet_build_plugin(
  DeDxAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  larevt::DetPedestalService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  FlashMatchCalib art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubcore::Geometry
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  GeoCosmicTagger art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::UBXSec
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  InTimeStudy art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  NeutrinoFlashMatch art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubcore::Geometry
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  NeutrinoMCFlash art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  PhotonActivity art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  RecoTrueMatcher art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  SimpleAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::LLBasicTool_GeoAlgo
  ubobj::UBXSec
  larsim::EventWeight_Base
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  StoppingMuonTagger art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_HitCosmicTag_Algorithms
  ubana::Utilities
  ubreco::UBFlashFinder
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  TPCObjectMaker art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubcore::LLBasicTool_GeoAlgo
  ubobj::UBXSec
  larpandora::LArPandoraInterface
  larsim::MCCheater_BackTrackerService_service
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  UBXSec art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::UBXSec_Algorithms
  ubreco::UBFlashFinder
  ubevt::Utilities
  ubobj::UBXSec
  ubobj::Trigger
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nug4::ParticleNavigation
  art_root_io::TFileService_service
  ROOT::Tree
  ROOT::EG
)

install_headers()
install_fhicl()
install_source()
