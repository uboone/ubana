art_make( BASENAME_ONLY
  LIBRARY_NAME ubana_ParticleID
  LIB_LIBRARIES
  ubana_ParticleID_Algorithms

  MODULE_LIBRARIES 
  ubobj_UBXSec
  uboonecode_uboone_UBXSec_Algorithms
  ubana_ParticleID_Algorithms
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
  larsim_MCCheater_BackTracker
  lardata_Utilities
  larevt_Filters
  larana_TruncatedMean_Algorithm
  larana_ParticleIdentification
  lardata_RecoObjects
  larreco_RecoAlg
  larpandora_LArPandoraInterface
  SignalShapingServiceMicroBooNE_service
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  nusimdata::SimulationBase
  ubraw_RawData
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  art::Utilities
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  )   

### install_headers()
install_fhicl()
install_source()
