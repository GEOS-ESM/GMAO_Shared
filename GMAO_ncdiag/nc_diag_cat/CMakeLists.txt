esma_set_this ()

set (SRCS
  gen_diff_size_netcdf.f90 kinds.F90 ncdc_climsg.F90 ncdc_cli_process.F90
  ncdc_data.F90
  ncdc_dims.F90 ncdc_metadata.F90 ncdc_realloc.F90 ncdc_state.F90
  ncdc_types.f90 ncdc_vars.F90 netcdf_unlimdims.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES nc_diag_write NetCDF::NetCDF_Fortran
  )
target_compile_definitions (${this} PRIVATE _REAL8_ ERROR_TRACEBACK QUIET)

ecbuild_add_executable (TARGET test_nc_unlimdims.x SOURCES test_nc_unlimdims.F90)
target_link_libraries (test_nc_unlimdims.x PUBLIC ${this} esmf)

ecbuild_add_executable (TARGET nc_diag_cat.x SOURCES nc_diag_cat.F90)
target_link_libraries (nc_diag_cat.x PUBLIC ${this} esmf)
