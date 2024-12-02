cet_build_plugin(
  MCC9Pi0Filter art::EDFilter
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
  lardata::RecoBaseProxy
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  MCC9SliceIDFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
  lardata::RecoBaseProxy
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  NuCCanalyzer art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::Filters_CCInclusiveFilter_helpers
)

cet_build_plugin(
  NuCCfilter art::EDFilter
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  lardataobj::RecoBase
)

cet_build_plugin(
  NuCCproducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::Filters_CCInclusiveFilter_helpers
)

add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
