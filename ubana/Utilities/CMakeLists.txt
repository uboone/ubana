cet_make_library(
  SOURCE
  FiducialVolume.cxx
  LIBRARIES
  PUBLIC
  ubcore::LLBasicTool_GeoAlgo
  larcorealg::GeoAlgo
  larcoreobj::SimpleTypesAndConstants
  fhiclcpp::fhiclcpp
  ROOT::MathCore
)

install_headers()
install_source()
