esma_set_this ()

set (SRCS
  VegType_Mod.F90
  VegType_io.c
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GMAO_gfio_r4
  INCLUDES ${INC_ESMF})

