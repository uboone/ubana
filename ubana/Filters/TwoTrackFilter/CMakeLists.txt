cet_make_library(
  SOURCE
  BackTrackerTruthMatch.cxx
  FiducialVolume.cxx
  PID.cxx
  RecoTruthMCParticle.cxx
  Topology.cxx
  TrackDirection.cxx
  LIBRARIES
  PUBLIC
  ubcore::LLBasicTool_GeoAlgo
  larpandora::LArPandoraInterface
  larsim::MCCheater_BackTrackerService_service
  larevt::DetPedestalService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art_root_io::TFileService_service
  canvas::canvas
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  ROOT::Physics
)

cet_build_plugin(
  OneMuOnePFlat art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::Filters_TwoTrackFilter
  ubobj::CRT
  larreco::RecoAlg
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  nurandom::RandomUtils_NuRandomService_service
  ROOT::Tree
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

