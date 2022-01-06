
art_make( MODULE_LIBRARIES larcorealg_Geometry
	nusimdata::SimulationBase
	lardataobj_Simulation
	lardataobj_RecoBase
	lardataobj_AnalysisBase
	larcoreobj_SummaryData
	larpandora_LArPandoraInterface
	ubobj_WcpPort
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common 
	art::Persistency_Provenance
	art::Utilities
	canvas::canvas
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib::cetlib
        cetlib_except::cetlib_except
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	ROOT::EG
)

install_fhicl()
install_headers()
install_source()
