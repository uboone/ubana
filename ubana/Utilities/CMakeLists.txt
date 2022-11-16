art_make_library(
   SOURCE
       FiducialVolume.cxx
   LIBRARIES
       ubcore::ubana_BasicTool_GeoAlgo
       larcoreobj::SimpleTypesAndConstants
       fhiclcpp::fhiclcpp
)

install_headers()
install_source()
