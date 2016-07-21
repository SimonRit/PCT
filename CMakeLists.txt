# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

# This project is designed to be built outside the Insight source tree.
PROJECT(ProtonCT)

#=========================================================
# Help function to debug CMake
MACRO (DD in)
  MESSAGE(${in}=${${in}})
ENDMACRO(DD)

#=========================================================
# CMake packages
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake"
    CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
    ${CMAKE_MODULE_PATH})

FIND_PACKAGE(Gengetopt REQUIRED)

FIND_PACKAGE(RTK REQUIRED)
INCLUDE(${RTK_USE_FILE})

IF(NOT WIN32)
  FIND_PACKAGE(ROOT REQUIRED)
  INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
ENDIF(NOT WIN32)

IF(NOT ITK_USE_FFTWD OR NOT ITK_USE_FFTWF)
    MESSAGE(FATAL_ERROR "You must recompile ITK with ITK_USE_FFTWF and ITK_USE_FFTWD. Then, recompile RTK first and pct second.")
ENDIF()

#=========================================================
# Geant4 stuff
OPTION(PCT_WITH_GEANT4 "Use Geant4" OFF)
IF(PCT_WITH_GEANT4)
  FIND_PACKAGE(Geant4 REQUIRED)
  INCLUDE(${Geant4_USE_FILE})
  ADD_SUBDIRECTORY(geant4)
  ADD_DEFINITIONS(-DPCT_GEANT4)
  INCLUDE_DIRECTORIES(geant4)
  STRING(REGEX REPLACE "-ansi"     "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  STRING(REGEX REPLACE "-pedantic" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
ELSE()
  FIND_PACKAGE(CLHEP 2.0 REQUIRED)
  INCLUDE_DIRECTORIES(${CLHEP_INCLUDE_DIRS})
ENDIF()

#=========================================================
# We remove the warning option "shadow", because there are tons of
# such warning related to clhep/g4 system of units.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")

#=========================================================
# Executables
MACRO (ADD_PCT_EXECUTABLE EXE_NAME)
  WRAP_GGO(${EXE_NAME}_GGO_C ${EXE_NAME}.ggo)
  ADD_EXECUTABLE(${EXE_NAME} ${EXE_NAME}.cxx ${${EXE_NAME}_GGO_C} ${ARGN})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${RTK_LIBRARIES})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${ITK_LIBRARIES})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${CLHEP_LIBRARIES})
  IF(PCT_WITH_GEANT4)
    TARGET_LINK_LIBRARIES(${EXE_NAME} pctGeant4)
  ENDIF(PCT_WITH_GEANT4)
ENDMACRO (ADD_PCT_EXECUTABLE)

# Executables
IF(NOT WIN32)
  ADD_PCT_EXECUTABLE(pctpairprotons pctpairprotons.cxx)
  TARGET_LINK_LIBRARIES(pctpairprotons ${ROOT_LIBRARIES})
  ADD_PCT_EXECUTABLE(pctpairprotonsLomaLinda pctpairprotonsLomaLinda.cxx)
  TARGET_LINK_LIBRARIES(pctpairprotonsLomaLinda ${ROOT_LIBRARIES})
  ADD_PCT_EXECUTABLE(pctpairprotonsLMU_IMPCT pctpairprotonsLMU_IMPCT.cxx)
  TARGET_LINK_LIBRARIES(pctpairprotonsLMU_IMPCT ${ROOT_LIBRARIES})
ENDIF(NOT WIN32)
ADD_PCT_EXECUTABLE(pctpaircuts pctpaircuts.cxx)
TARGET_LINK_LIBRARIES(pctpaircuts ${ROOT_LIBRARIES} RooFit RooFitCore Minuit)
ADD_PCT_EXECUTABLE(pctpairgeometry pctpairgeometry.cxx)
ADD_PCT_EXECUTABLE(pctmostlikelypath pctmostlikelypath.cxx)
ADD_PCT_EXECUTABLE(pctbinning pctbinning.cxx)
ADD_PCT_EXECUTABLE(pctfillholl pctfillholl.cxx)
ADD_PCT_EXECUTABLE(pctbackprojectionbinning pctbackprojectionbinning.cxx)
ADD_PCT_EXECUTABLE(pctzengbackprojections pctzengbackprojections.cxx)
ADD_PCT_EXECUTABLE(pctparkershortscanweighting pctparkershortscanweighting.cxx)
ADD_PCT_EXECUTABLE(pctprojections pctprojections.cxx)
ADD_PCT_EXECUTABLE(pctbackprojections pctbackprojections.cxx)
ADD_PCT_EXECUTABLE(pctfdktwodweights pctfdktwodweights.cxx)
ADD_PCT_EXECUTABLE(pctfdk pctfdk.cxx)
ADD_PCT_EXECUTABLE(pctschulte pctschulte.cxx)
ADD_PCT_EXECUTABLE(pctswapcoordinates pctswapcoordinates.cxx)
ADD_PCT_EXECUTABLE(pctpairarithm pctpairarithm.cxx)
IF(PCT_WITH_GEANT4)
  ADD_PCT_EXECUTABLE(pctdiffxsection pctdiffxsection.cxx)
ENDIF(PCT_WITH_GEANT4)

