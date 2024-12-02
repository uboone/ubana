cet_build_plugin(
  CalibratedEdx art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
)

cet_build_plugin(
  ParticleIdValidationPlots art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubana::ParticleID_Algorithms
  ubobj::UBXSec
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  ParticleId art::EDProducer
  LIBRARIES
  PRIVATE
  ubana::ParticleID_Algorithms
  larana::ParticleIdentification
  larana::TruncatedMean_Algorithm
  larreco::RecoAlg
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
)

### install_headers()
install_fhicl()
install_source()
