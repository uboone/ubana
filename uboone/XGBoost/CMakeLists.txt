link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

include_directories (
    xgboost 
    dmlc/include
    rabit/include 
    )

file(GLOB SOURCES 
    xgboost/src/c_api/*.cc
    xgboost/src/common/*.cc
    xgboost/src/data/*.cc
    xgboost/src/gbm/*.cc
    xgboost/src/metric/*.cc
    xgboost/src/objective/*.cc
    xgboost/src/tree/*.cc
    xgboost/src/*.cc
    rabit/*.cc
    rabit/src/*.cc
    dmlc/*.cc
    dmlc/src/*.cc
    dmlc/src/data/*.cc
    dmlc/src/io/*.cc
)


art_make_library( LIBRARY_NAME  uboonecode_uboone_XGBoost
                  SOURCE ${SOURCES}
                  LIBRARIES
		                ${CETLIB}
                    
                )

install_headers()
install_source()
