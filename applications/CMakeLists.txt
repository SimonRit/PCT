find_package(Gengetopt)


#-----------------------------------------------------------------------------
# Setup RPATH :
#
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# add path to PCT libraries to the install path
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${PCT_INSTALL_LIB_DIR}")

#-----------------------------------------------------------------------------
# Find ITK.
# Required to include ITK_USE_FILE in order to Register IO factories
# Force requested modules to be PCT dependencies only, otherwise all
# available factories will try to register themselves.
if (NOT ITK_DIR)
  set(ITK_DIR ${ITK_BINARY_DIR}/CMakeTmp)
endif()
find_package(ITK REQUIRED COMPONENTS ${ITK_MODULE_PCT_DEPENDS})
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Executables
add_subdirectory(pctfdktwodweights)
add_subdirectory(pctpaircuts)
add_subdirectory(pctparkershortscanweighting)
add_subdirectory(pctbackprojectionbinning)
add_subdirectory(pctfillholl)
add_subdirectory(pctpairgeometry)
add_subdirectory(pctprojections)
add_subdirectory(pctbinning)
add_subdirectory(pctlluconverter)
add_subdirectory(pctpairprotons)
add_subdirectory(pctzengbackprojections)
add_subdirectory(pctdiffxsection)
add_subdirectory(pctmostlikelypath)
add_subdirectory(pctpairprotonsLMU_IMPCT)
add_subdirectory(pctfdk)
add_subdirectory(pctpairarithm)
add_subdirectory(pctpairprotonsLomaLinda)

#=========================================================

