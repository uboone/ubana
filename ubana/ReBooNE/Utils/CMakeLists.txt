cet_make_library(
  SOURCE
  BooNENtuple.cc
  BeamNtuple.cc
  BooNEInfo.cc
  LIBRARIES
  PUBLIC
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  canvas::canvas
)

art_dictionary(
  DICTIONARY_LIBRARIES
  lardataobj::RecoBase
  nusimdata::SimulationBase
)


install_headers()
install_fhicl()
install_source()
