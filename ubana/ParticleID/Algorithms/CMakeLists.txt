cet_make_library(
  SOURCE
  Bragg_Likelihood_Estimator.cxx
  FiducialVolume.cxx
  GetDaughterTracksShowers.cxx
  KernelDensityEstimator.cxx
  PIDA.cxx
  Theory_dEdx_resrange.cxx
  LIBRARIES
  PUBLIC
  lardataobj::RecoBase
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::MathCore
)

install_fhicl()
install_headers()
install_source()
