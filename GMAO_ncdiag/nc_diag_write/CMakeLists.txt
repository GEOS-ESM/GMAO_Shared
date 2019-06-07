esma_set_this ()

set (SRCS
  kinds.F90 nc_diag_write_mod.F90 ncdw_chaninfo.F90 ncdw_ciresize.F90
  ncdw_climsg.F90 ncdw_data2d.F90 ncdw_dresize.F90 ncdw_lheader.F90
  ncdw_metadata.F90 ncdw_mresize.F90 ncdw_realloc.F90 ncdw_state.f90
  ncdw_strarrutils.F90 ncdw_types.F90 ncdw_varattr.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  INCLUDES ${INC_NETCDF}
  )

target_compile_definitions (${this} PRIVATE _REAL8_ ERROR_TRACEBACK ENABLE_ACTION_MSGS)

