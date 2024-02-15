cet_build_plugin(
  NCRadiativeResonant art::EDFilter
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  SinglePhoton art::EDFilter
  LIBRARIES
  PRIVATE
  ubana::SinglePhotonAnalysis_SEAview
  ubobj::CRT
  larpandora::LArPandoraInterface
  larsim::EventWeight_Base
  larevt::SpaceChargeServices_SpaceChargeServiceStandard_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

                      
add_subdirectory(SEAview)
add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
