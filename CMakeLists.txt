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

FIND_PACKAGE(ROOT REQUIRED)
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})

FIND_PACKAGE(CLHEP 2.0 REQUIRED)
INCLUDE_DIRECTORIES(${CLHEP_INCLUDE_DIRS})
#INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/CLHEP")

#=========================================================
# Executables
MACRO (ADD_PCT_EXECUTABLE EXE_NAME)
  WRAP_GGO(${EXE_NAME}_GGO_C ${EXE_NAME}.ggo)
  ADD_EXECUTABLE(${EXE_NAME} ${EXE_NAME}.cxx ${${EXE_NAME}_GGO_C} ${ARGN})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${RTK_LIBRARIES})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${ITK_LIBRARIES})
  TARGET_LINK_LIBRARIES(${EXE_NAME} ${ROOT_LIBRARIES})
ENDMACRO (ADD_PCT_EXECUTABLE)

# Executable
ADD_PCT_EXECUTABLE(pctpairprotons pctpairprotons.cxx)
ADD_PCT_EXECUTABLE(pctpaircuts pctpaircuts.cxx)
ADD_PCT_EXECUTABLE(pctprojection pctprojection.cxx)
ADD_PCT_EXECUTABLE(pctbackprojections pctbackprojections.cxx)
ADD_PCT_EXECUTABLE(pctfdktwodweights pctfdktwodweights.cxx)
ADD_PCT_EXECUTABLE(pctfdk pctfdk.cxx)
