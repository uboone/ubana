cet_build_plugin(
  HsnFinder art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::HSNAnalysis_HsnFinder_Algorithms
)

add_subdirectory(Algorithms)
add_subdirectory(DataObjects)

install_headers()
install_source()
install_fhicl()
