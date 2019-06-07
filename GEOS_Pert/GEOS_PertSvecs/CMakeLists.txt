esma_set_this ()

set (SRCS
  svecs_cf.F90
  svecs_norms.F90
  svecs_proj.F90
  svecs_tlmadm.F90
  svecs_utils.F90
  svecs_weights.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_cfio_r4 MAPL_Base GMAO_mpeu GEOS_PertShared
  INCLUDES ${INC_ESMF})
