#
# Find the packages required by this module
#
list(APPEND CMAKE_MODULE_PATH ${PCT_SOURCE_DIR}/cmake)

find_package(CUDA_wrap QUIET)
if(CUDA_FOUND)
  if(${CUDA_VERSION} VERSION_LESS 8.0)
    message(WARNING "CUDA version ${CUDA_VERSION} is not supported by PCT.")
    set(PCT_USE_CUDA_DEFAULT OFF)
  else()
    set(PCT_USE_CUDA_DEFAULT ON)
  endif()
else()
  set(PCT_USE_CUDA_DEFAULT OFF)
endif()
option(PCT_USE_CUDA "Use CUDA for PCT" ${PCT_USE_CUDA_DEFAULT})

if(PCT_USE_CUDA)
  if(NOT CUDA_FOUND)
    find_package(CUDA_wrap REQUIRED)
  endif()
  set(PCT_CUDA_PROJECTIONS_SLAB_SIZE "16" CACHE STRING "Number of projections processed simultaneously in CUDA forward and back projections")
endif()

