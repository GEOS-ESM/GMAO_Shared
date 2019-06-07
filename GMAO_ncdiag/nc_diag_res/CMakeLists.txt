esma_set_this ()

set (SRCS
  nc_diag_res.f90 ncdres_climsg.F90
  nc_diag_fson/src/nc_diag_fson.f90 nc_diag_fson/src/ncdf_value_m.f90
  nc_diag_fson/src/ncdf_string_m.f90 nc_diag_fson/src/ncdf_path_m.f90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  )
