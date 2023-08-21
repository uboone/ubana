cet_build_plugin(
  Arborist art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
