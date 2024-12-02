cet_make_library(
  SOURCE
  EnergyHelper.cc
  PandoraInterfaceHelper.cxx
  TrackHelper.cc
  LIBRARIES
  PUBLIC
  larpandora::LArPandoraInterface
  larreco::RecoAlg
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
)

install_headers()
install_source()
install_fhicl()
