cet_build_plugin(
  ReBooNE art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::ReBooNE_Utils
  lardata::Utilities
  nusimdata::SimulationBase
  dk2nu::Tree
  ROOT::Tree
)


add_subdirectory(Utils)


install_headers()
install_fhicl()
install_source()
