cet_build_plugin(
  ReDk2Nu art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  nusimdata::SimulationBase
  dk2nu::Tree
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
