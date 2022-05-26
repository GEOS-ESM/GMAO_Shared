esma_set_this ()

set (SRCS
  kinds.F90 nc_diag_read_mod.F90 ncdr_alloc_assert.f90 ncdr_attrs.f90 ncdr_attrs_fetch.f90
  ncdr_check.f90 ncdr_climsg.F90 ncdr_dims.f90 ncdr_global_attrs.f90 ncdr_global_attrs_fetch.f90
  ncdr_realloc_mod.F90 ncdr_state.f90 ncdr_types.f90 ncdr_vars.f90 ncdr_vars_fetch.f90 netcdf_unlimdims.F90
  )

esma_add_library (${this}
  SRCS ${SRCS} 
  DEPENDENCIES NetCDF::NetCDF_Fortran
  )

target_compile_definitions (${this} PRIVATE _REAL8_ ERROR_TRACEBACK)

