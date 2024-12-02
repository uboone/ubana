cet_build_plugin(
  DecisionTreeID art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::NCElastic_Algorithms
  ubana::XGBoost
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  QEEventSelectionFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::NCElastic_Algorithms
  larsim::EventWeight_Base
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  QEEventSelection art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::NCElastic_Algorithms
  larsim::EventWeight_Base
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)


install_headers()
install_fhicl()
install_source()

add_subdirectory(job)
add_subdirectory(Algorithms)
