file(GLOB Hadr01src ${PROJECT_SOURCE_DIR}/geant4/*.cc)

ADD_LIBRARY(pctGeant4 pctGeant4.cxx ${Hadr01src})

target_link_libraries(pctGeant4 ${Geant4_LIBRARIES} )

