cet_make_library(
  SOURCE
  AltNuMuCCInclusiveAlg.cxx
  ChargedTrackMultiplicityAlg.cxx
  Cluster2DNuAlg.cxx
  ModNuMuCCInclusiveAlg.cxx
  NeutrinoIDAlgFactory.cxx
  NuMuCCInclusiveAlg.cxx
  NuMuCCSelectionIIAlg.cxx
  NuMuCCSelectionIIAlgMCC7.cxx
  TrackPairPlusVertexAlg.cxx
  LIBRARIES
  PUBLIC
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  larcorealg::Geometry
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
)

install_headers()
install_source()
install_fhicl()
