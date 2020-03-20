# For saving GENIE event records to the NuCCanalyzer module output
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Determine the major version number of GENIE based on the GENIE_VERSION
# variable set by the genie ups product. Define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think
# that we can access the result of that version check in the ubsim build.
parse_ups_version( ${GENIE_VERSION} )
if( ${VMAJ} LESS 3 )

  add_definitions(-DGENIE_PRE_R3)

  # GENIE v2 library list
  set(GENIE_LIB_LIST  ${GALGORITHM}
                    ${GBARYONRESONANCE}
                    ${GBASE}
                    ${GBODEKYANG}
                    ${GCHARM}
                    ${GCOH}
                    ${GDFRC}
                    ${GDIS}
                    ${GCROSSSECTIONS}
                    ${GDECAY}
                    ${GELAS}
                    ${GELFF}
                    ${GHEP}
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GEVGDRIVERS}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GINTERACTION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GMESSENGER}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNTUPLE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GUTILS}
                    ${GGEO}
                    ${GFLUXDRIVERS}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )

  #MESSAGE("--UBSIM-- GENIE_VERSION=${VMAJ} ${VMIN} -DGENIE_PRE_R3")
else()

  # GENIE v3 library list
  set (GENIE_LIB_LIST ${GFWMSG}
                      ${GFWREG}
                      ${GFWALG}
                      ${GFWINT}
                      ${GFWGHEP}
                      ${GFWNUM}
                      ${GSL}  # FWNUM relies on GSL
                      ${GFWUTL}
                      ${GFWPARDAT}
                      ${GFWEG}
                      ${GFWNTP}
                      ${GPHXSIG}
                      ${GPHPDF}
                      ${GPHNUCLST}
                      ${GPHCMN}
                      ${GPHDCY}
                      ${GPHHADTRANSP}
                      ${GPHHADNZ}
                      ${GPHDEEX}
                      ${GPHAMNGXS}
                      ${GPHAMNGEG}
                      ${GPHCHMXS}
                      ${GPHCOHXS}
                      ${GPHCOHEG}
                      ${GPHDISXS}
                      ${GPHDISEG}
                      ${GPHDFRCXS}
                      ${GPHDFRCEG}
                      ${GPHGLWRESXS}
                      ${GPHGLWRESEG}
                      ${GPHIBDXS}
                      ${GPHIBDEG}
                      ${GPHMNUCXS}
                      ${GPHMNUCEG}
                      ${GPHMEL}
                      ${GPHNUELXS}
                      ${GPHNUELEG}
                      ${GPHQELXS}
                      ${GPHQELEG}
                      ${GPHRESXS}
                      ${GPHRESEG}
                      ${GPHSTRXS}
                      ${GPHSTREG}
                      ${GPHNDCY}
                      ${GTLGEO}
                      ${GTLFLX}
                      ${GRWFWK}
                      ${GRWIO}
                      ${GRWCLC})

  #MESSAGE("--UBSIM-- GENIE_VERSION=${VMAJ} ${VMIN} not defining -DGENIE_PRE_R3")
endif()



art_make( 
          LIBRARY_NAME uboonecode_uboone_Filters
          # Needed for standalone c++ helper files
          LIB_LIBRARIES 

          # Needed for the main module          
          MODULE_LIBRARIES
                    larpandora_LArPandoraInterface
                    larcorealg_Geometry
                    larreco_RecoAlg
                    larcore_Geometry_Geometry_service
                    larsim_Simulation lardataobj_Simulation
                    lardataobj_RawData
                    lardataobj_RecoBase
                    lardataobj_AnalysisBase
                    lardata_Utilities
                    nusimdata_SimulationBase
                    larsim_EventWeight_Base
                    ${ART_FRAMEWORK_CORE}
                    ${ART_FRAMEWORK_PRINCIPAL}
                    ${ART_FRAMEWORK_SERVICES_REGISTRY}
                    ${ART_FRAMEWORK_SERVICES_BASIC}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                    art_Framework_Services_Registry
                    art_Persistency_Common
                    art_Persistency_Provenance
                    art_Utilities
                    canvas
                    ${MF_MESSAGELOGGER}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    ${Boost_SYSTEM_LIBRARY}
                    ${ROOT_BASIC_LIB_LIST}
                    NuCCHelpers
                    ${GENIE_LIB_LIST}
                    nutools_EventGeneratorBase_GENIE
          )

add_subdirectory(job)
add_subdirectory(helpers)

install_headers()
install_source()
install_fhicl()
