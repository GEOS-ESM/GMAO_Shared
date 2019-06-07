esma_set_this ()

esma_add_library (${this}
  SRCS adBuffer.f adStack.c tapenade_iter.F90
  DEPENDENCIES MAPL_cfio_r4
  INCLUDES ${INC_ESMF})

