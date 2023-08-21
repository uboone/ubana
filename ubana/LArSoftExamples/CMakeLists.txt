cet_build_plugin(
  UsingAssnsAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  UsingHandlesAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

# install_headers()
install_fhicl()
install_source()
