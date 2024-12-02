cet_make_library(
  LIBRARY_NAME IHitEfficiencyHistogramTool INTERFACE
  SOURCE IHitEfficiencyHistogramTool.h
  LIBRARIES
  INTERFACE
  art_root_io::TFileService_service
)

cet_write_plugin_builder(
  uboone_tool::IHitEfficiencyHistogramTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::IHitEfficiencyHistogramTool
)

include(uboone_tool::IHitEfficiencyHistogramTool)

cet_build_plugin(
  ShowerHitEfficiencyAnalysis uboone_tool::IHitEfficiencyHistogramTool
  LIBRARIES
  PRIVATE
  larevt::ChannelStatusService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  TrackHitEfficiencyAnalysis uboone_tool::IHitEfficiencyHistogramTool
  LIBRARIES
  PRIVATE
  larevt::ChannelStatusService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  nusimdata::SimulationBase
  Eigen3::Eigen
  ROOT::Tree
)

install_headers()
install_source()
install_fhicl()
