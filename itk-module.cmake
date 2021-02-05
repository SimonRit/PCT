set(DOCUMENTATION "")

# -----------------------------------------
#  Required Modules to build PCT library :
set(PCT_IO_DEPENDS
  ITKIOCSV
  ITKIOGDCM
  ITKIOMeta
  ITKIORAW
  ITKIOTIFF
  ITKIOXML
  )

set(PCT_DEPENDS
  ITKCommon
  ITKConvolution
  ITKFFT
  ITKOptimizers
  ITKRegistrationCommon
  ITKSmoothing
  ITKImageNoise
  RTK
  ${PCT_IO_DEPENDS}
  )

# -----------------------------------------
#  Required Modules to build PCT tests :
set(PCT_TEST_DEPENDS
  ITKTestKernel)

# # -----------------------------------------
# # CUDA optional dependencies
if(ITK_SOURCE_DIR)
  if(${PCT_USE_CUDA})
    list(APPEND PCT_DEPENDS ITKCudaCommon)
  endif()
endif()

#=========================================================
# Module PCT
#=========================================================
itk_module(PCT
  ENABLE_SHARED
  EXCLUDE_FROM_DEFAULT
  DEPENDS
    ${PCT_DEPENDS}
  TEST_DEPENDS
    ${PCT_TEST_DEPENDS}
  DESCRIPTION
    "${DOCUMENTATION}"
  )
