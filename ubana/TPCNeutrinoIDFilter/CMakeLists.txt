# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
  LIBRARY_NAME TPCNeutrinoIDFilter
  LIB_LIBRARIES
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_Utilities
  art::Persistency_Common
  art::Persistency_Provenance
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support ROOT::Core
  art_root_io::TFileService_service
  art::Utilities
  canvas::canvas
  messagefacility::MF_MessageLogger

  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except


  MODULE_LIBRARIES
  TPCNeutrinoIDAlgs
  lardataalg_DetectorInfo
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_Utilities
  nusimdata::SimulationBase
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support ROOT::Core
  art_root_io::TFileService_service
  art::Utilities
  canvas::canvas
  messagefacility::MF_MessageLogger

  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except


  )

add_subdirectory(Algorithms)

install_headers()
install_source()
install_fhicl()
