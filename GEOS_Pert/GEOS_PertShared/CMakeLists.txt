esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_PertSharedMod.F90
  DEPENDENCIES MAPL_cfio_r4 MAPL_Base GMAO_mpeu GMAO_hermes
  INCLUDES ${INC_ESMF})

