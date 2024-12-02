cet_make_library(
  SOURCE
  IncompatibilityChecker.cxx
  LIBRARIES
  PUBLIC
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  fhiclcpp::fhiclcpp
)

install_headers()
#install_fhicl()
install_source()
