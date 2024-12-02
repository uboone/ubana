cet_make_library(
  LIBRARY_NAME SelectionToolBase INTERFACE
  SOURCE SelectionToolBase.h
  LIBRARIES
  INTERFACE
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::SelectionToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::SelectionToolBase
  lardata::RecoBaseProxy
  lardataobj::RecoBase
)

include(uboone_tool::SelectionToolBase)

cet_build_plugin(
  CC0piNpSelection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
  ubevt::Database_TPCEnergyCalib
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  ROOT::EG
)

cet_build_plugin(
  CRTApproachSelection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
  lardataobj::AnalysisBase
)

cet_build_plugin(
  EmptySelection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  NuMuSelection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
  ubevt::Database_TPCEnergyCalib
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  ROOT::EG
)

cet_build_plugin(
  Pi0Selection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  SelectionExample uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  ShowerSelection uboone_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
)

install_headers()
install_source()
install_fhicl()
