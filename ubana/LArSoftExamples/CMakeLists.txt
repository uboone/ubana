# CMakeLists are the devil and anybody who tells you they understand them is a liar
# This one is stolen from the AnalysisTree directory, which is a good bet for
# a first try.

art_make(
  MODULE_LIBRARIES
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation
    nug4::ParticleNavigation
    lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    lardataobj_RawData
    lardataobj_RecoBase
    lardata_RecoObjects
    larreco_RecoAlg
    larpandora_LArPandoraInterface
    ubevt::SignalShapingServiceMicroBooNE_service
    nusimdata::SimulationBase
    ubraw_RawData
    ubobj_Trigger
    ubobj_Optical
    ${LARRECO_LIB}
    ${LARDATA_LIB}
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support ROOT::Core
    art_root_io::TFileService_service
    art::Persistency_Common canvas::canvas
    art::Persistency_Provenance canvas::canvas
    art::Utilities canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::Geom
    ROOT::XMLIO
    ROOT::Gdml

  )

# install_headers()
install_fhicl()
install_source()
