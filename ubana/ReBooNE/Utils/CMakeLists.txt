cet_make_library(
  SOURCE
  BooNENtuple.cc
  BeamNtuple.cc
  LIBRARIES
  PUBLIC
  ubobj::ReBooNE
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  canvas::canvas
)

install_headers()
install_source()
