esma_set_this (OVERRIDE putil)

set (SRCS
  pcmout.f  pcvout.f  pdmout.f  pdvout.f  pivout.f  psmout.f  psvout.f  pzmout.f  pzvout.f
  )

esma_add_library (${this}
  SRCS ${SRCS}
  )
