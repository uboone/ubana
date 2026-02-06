cet_make_library(
  SOURCE
  GSimpleInterface.cxx
  NTupleInterface.cxx
  TreeReader.cxx
  LIBRARIES
  PUBLIC
  larcoreobj::SummaryData
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_IO_Sources
  art_root_io::TFileService_service
  GENIE::GTlFlx
  ROOT::Tree
)

cet_build_plugin(
  TreeReader art::source
  LIBRARIES
  PRIVATE
  ubana::TreeReader
)

add_subdirectory(job)

install_headers()
install_source()
install_fhicl()
install_scripts()
