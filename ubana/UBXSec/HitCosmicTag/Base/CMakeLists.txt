cet_make_library(
  SOURCE
  BaseAlgorithm.cxx
  BaseDqDsCalculatorAlgo.cxx
  BaseDqDsSmootherAlgo.cxx
  BaseHitOrdererAlgo.cxx
  BaseHitSmootherAlgo.cxx
  BaseLocalLinearityCalculatorAlgo.cxx
  BaseStartHitFinderAlgo.cxx
  CosmicTagManager.cxx
  CustomAlgoFactory.cxx
  DqDsCalculatorFactory.cxx
  DqDsSmootherFactory.cxx
  HitCosmicTagException.cxx
  HitCosmicTagLogger.cxx
  HitOrdererFactory.cxx
  HitSmootherFactory.cxx
  LocalLinearityCalculatorFactory.cxx
  StartHitFinderFactory.cxx
  Tools.cxx
  LIBRARIES
  PUBLIC
  fhiclcpp::fhiclcpp
)

install_headers()
install_source()
