#=========================================================
# pCT : pronton CT
#=========================================================
# Respect the CMAKE_CXX_STANDARD flags when building for
# ITKv5 or C++11.
cmake_minimum_required(VERSION 3.9.5 FATAL_ERROR)

## Only policies introduced after the cmake_minimum_required
## version need to explicitly be set to NEW.
## Refer to https://cmake.org/cmake/help/v3.11/manual/cmake-policies.7.html
set(CMAKE_POLICIES
  CMP0015
  CMP0070
  CMP0071
  CMP0072)

foreach(p ${CMAKE_POLICIES})
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

#=========================================================
# Help function to debug CMake
macro (DD in)
  message(${in}=[${${in}}])
endmacro()
#=========================================================

project(ProtonCT)

## PCT Version
set(PCT_VERSION_MAJOR "1")
set(PCT_VERSION_MINOR "0")
set(PCT_VERSION_PATCH "0")
set(PCT_VERSION_STRING "${PCT_VERSION_MAJOR}.${PCT_VERSION_MINOR}")

set(PCT_LIBRARIES PCT)

## Default to release
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

# CMake packages
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake"
    CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
    ${CMAKE_MODULE_PATH})

IF(NOT WIN32)
  FIND_PACKAGE(ROOT REQUIRED)
  INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
ENDIF(NOT WIN32)

#=========================================================
# Installation variables
#=========================================================
if(NOT PCT_INSTALL_RUNTIME_DIR)
  set(PCT_INSTALL_RUNTIME_DIR bin)
endif()
if(NOT PCT_INSTALL_LIB_DIR)
  set(PCT_INSTALL_LIB_DIR lib)
endif()
if(NOT PCT_INSTALL_ARCHIVE_DIR)
  set(PCT_INSTALL_ARCHIVE_DIR lib)
endif()
if(NOT PCT_INSTALL_INCLUDE_DIR)
  set(PCT_INSTALL_INCLUDE_DIR include/PCT)
endif()
if(NOT PCT_INSTALL_PACKAGE_DIR)
  set(PCT_INSTALL_PACKAGE_DIR "lib/cmake/PCT")
endif()

#Set position independent code for Unix (-fPIC)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#=========================================================
# Remove some MS Visual c++ flags
if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
endif()

#=========================================================
# Remove some Intel compiler warnings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  string(REPLACE "-Wno-unused-parameter" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  if(WIN32)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qwd1268")
  else()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1268")
  endif()
endif()

#=========================================================
# Allow for multiple CL.EXE to write to the same .PDB file

if(WIN32)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FS")
endif()


# --------------------------------------------------------
# Find ITK (required)
if(NOT ITK_SOURCE_DIR)
  find_package(ITK 5.0.0 REQUIRED COMPONENTS RTK)
endif()

# --------------------------------------------------------
# Shared libraries option
if(NOT ITK_SOURCE_DIR)
  set(PCT_BUILD_SHARED_LIBS ${ITK_BUILD_SHARED})
else()
  set(PCT_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
endif()

# ----------------------------------------------------------------------------
# Set PCT_DATA_ROOT
# Internally used by XRadRawToAttenuationImageFilter to set the path for its
# flat and dark headers.
if(NOT ITK_SOURCE_DIR)
  set(PCT_DATA_ROOT ${CMAKE_BINARY_DIR}/ExternalData/test CACHE PATH "Path of the data root" FORCE)
else()
  file(RELATIVE_PATH PCT_RELATIVE_PATH ${ITK_SOURCE_DIR} ${PCT_SOURCE_DIR})
  set(PCT_DATA_ROOT ${ExternalData_BINARY_ROOT}/${PCT_RELATIVE_PATH}/test CACHE PATH "Path of the data root" FORCE)
endif()
mark_as_advanced(PCT_DATA_ROOT)

# --------------------------------------------------------
# Initialization
if(NOT ITK_SOURCE_DIR)
  include(itk-module-init.cmake)
endif()

# Propagate cmake options in a header file
configure_file(${PCT_SOURCE_DIR}/pctConfiguration.h.in
  ${PCT_BINARY_DIR}/pctConfiguration.h)

#=========================================================
# If choose to build documentation, then search for Doxygen executables.
option(PCT_BUILD_DOXYGEN "Build Doxygen Documentation" OFF)
if(PCT_BUILD_DOXYGEN)
  add_subdirectory(documentation/Doxygen)
endif()

# Setup build locations.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PCT_BINARY_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PCT_BINARY_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PCT_BINARY_DIR}/lib)
endif()

# --------------------------------------------------------
# Setup KWStyle from ITK
if(ITK_USE_KWSTYLE)
  set(WORKING_DIR "${CMAKE_SOURCE_DIR}")
  if(ITK_SOURCE_DIR)
    set(WORKING_DIR "${ITK_SOURCE_DIR}")
  endif()

  find_package( Git )
  if(GIT_FOUND AND EXISTS "${WORKING_DIR}/.git/config")
    execute_process( COMMAND ${GIT_EXECUTABLE} config hooks.KWStyle.path
      "${KWSTYLE_EXECUTABLE}"
      WORKING_DIRECTORY ${WORKING_DIR} )
  endif()
endif()

#=========================================================
# ITKCudaCommon
#=========================================================
if(PCT_USE_CUDA)
  if(NOT TARGET ITKCudaCommon)
    add_subdirectory(utilities/ITKCudaCommon)
  endif()
endif()
list(APPEND PCT_INCLUDE_DIRS
  ${ITKCudaCommon_INCLUDE_DIRS})
#=========================================================
# Include directories
#=========================================================
list(APPEND PCT_INCLUDE_DIRS
  ${PCT_BINARY_DIR})

# Export library linking directories (used by LIBPATH)
set(PCT_SYSTEM_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PCT_SYSTEM_INCLUDE_DIRS ${PCT_INCLUDE_DIRS})

#=========================================================
# Generate PCTConfig.cmake for the build tree.
set(PCT_MODULE_PATH_CONFIG ${CMAKE_MODULE_PATH})
set(PCT_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(PCT_EXPORT_CODE_BUILD "
# The PCT version number
set(PCT_VERSION_MAJOR ${PCT_VERSION_MAJOR})
set(PCT_VERSION_MINOR ${PCT_VERSION_MINOR})
set(PCT_VERSION_PATCH ${PCT_VERSION_PATCH})

# Whether the compiled version of PCT uses CUDA
set(PCT_USE_CUDA ${PCT_USE_CUDA})

if(${PCT_USE_CUDA})
  set(CMAKE_MODULE_PATH ${PCT_MODULE_PATH_CONFIG} ${CMAKE_MODULE_PATH})
  # if we are using CUDA, make sure CUDA libraries are available
  find_package(CUDA_wrap REQUIRED)
  include_directories(${CUDA_INCLUDE_DIRS})
  set(PCT_CUDA_PROJECTIONS_SLAB_SIZE \"16\" CACHE STRING \"Number of projections processed simultaneously in CUDA forward and back projections\")
endif()
")

set(PCT_EXPORT_CODE_INSTALL "
# The PCT version number
set(PCT_VERSION_MAJOR ${PCT_VERSION_MAJOR})
set(PCT_VERSION_MINOR ${PCT_VERSION_MINOR})
set(PCT_VERSION_PATCH ${PCT_VERSION_PATCH})

# Whether the compiled version of PCT uses CUDA
set(PCT_USE_CUDA ${PCT_USE_CUDA})

if(${PCT_USE_CUDA})
  # if we are using CUDA, make sure CUDA libraries are available
  find_package(CUDA_wrap REQUIRED)
  set(PCT_CUDA_PROJECTIONS_SLAB_SIZE \"16\" CACHE STRING \"Number of projections processed simultaneously in CUDA forward and back projections\")
endif()
")

#=========================================================
# Configure and build ITK external module
#=========================================================
if(NOT ITK_SOURCE_DIR)
  if(NOT EXISTS ${ITK_CMAKE_DIR}/ITKModuleMacros.cmake)
    message(FATAL_ERROR "Modules can only be built against an ITK build tree; they cannot be built against an ITK install tree.")
  endif()

  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  itk_module_impl()
endif()

# Install configuration file
install(FILES ${PCT_BINARY_DIR}/PCTConfiguration.h DESTINATION ${PCT_INSTALL_INCLUDE_DIR})
install(FILES ${PCT_SOURCE_DIR}/cmake/FindCUDA_wrap.cmake
              ${PCT_SOURCE_DIR}/cmake/FindGengetopt.cmake
              ${PCT_SOURCE_DIR}/cmake/GetGitRevisionDescription.cmake
              ${PCT_SOURCE_DIR}/cmake/GetGitRevisionDescription.cmake.in
              ${PCT_SOURCE_DIR}/cmake/has_cuda_gpu.cxx
              ${PCT_SOURCE_DIR}/cmake/nvcc-check.cmake
        DESTINATION ${ITK_INSTALL_PACKAGE_DIR})
if(NOT ITK_SOURCE_DIR)
  install(CODE "MESSAGE(FATAL_ERROR \"Cannot install, PCT is compiled separately from ITK. Installation is only functional if PCT is compiled within ITK.\")")
endif()

IF(NOT ITK_USE_FFTWD OR NOT ITK_USE_FFTWF)
    MESSAGE(FATAL_ERROR "You must recompile ITK with ITK_USE_FFTWF and ITK_USE_FFTWD. Then, recompile PCT.")
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
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

#=========================================================
# We remove the warning option "shadow", because there are tons of
# such warning related to clhep/g4 system of units.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")

#=========================================================
# Executables
option(PCT_BUILD_APPLICATIONS "Build pCT applications" ON)
if(PCT_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif(PCT_BUILD_APPLICATIONS)

