itk_wrap_module(PCT)
itk_auto_load_submodules()
itk_end_wrap_module()

list(
  APPEND
  PCT_PYTHON_APP
  ${PCT_SOURCE_DIR}/applications/pctfdktwodweights/pctfdktwodweights.py
  ${PCT_SOURCE_DIR}/applications/pctpairprotons/pctpairprotons.py
  ${PCT_SOURCE_DIR}/applications/pctweplfit/pctweplfit.py
  ${PCT_SOURCE_DIR}/applications/pctlomalinda/pctlomalinda.py
  ${PCT_SOURCE_DIR}/applications/pctfilterprotons/pctfilterprotons.py
  ${PCT_SOURCE_DIR}/gate/protonct.py
)

# Install python applications
wrap_itk_python_bindings_install(/itk "RTK"  ${PCT_PYTHON_APP})

# Copy python applications to the ITK wrapping directory to ensure they can be imported in tests.
# This directory is added to the PYTHONPATH by itk_python_add_test.
if(ITK_DIR)
  set(itk_wrap_python_binary_dir "${ITK_DIR}/Wrapping/Generators/Python")
else()
  set(itk_wrap_python_binary_dir "${ITK_BINARY_DIR}/Wrapping/Generators/Python")
endif()
file(COPY ${PCT_PYTHON_APP} DESTINATION "${itk_wrap_python_binary_dir}")
