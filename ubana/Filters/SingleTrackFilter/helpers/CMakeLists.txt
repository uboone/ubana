cet_make_library(
  SOURCE
  PandoraInterfaceHelper.cxx
  LIBRARIES
  PUBLIC
  larpandora::LArPandoraInterface
  larevt::SpaceChargeService
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardataobj::RecoBase
  art::Framework_Principal
)

install_headers()
install_source()
install_fhicl()
