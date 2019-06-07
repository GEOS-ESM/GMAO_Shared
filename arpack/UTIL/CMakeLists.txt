esma_set_this (OVERRIDE util)

set (SRCS
  cmout.f  dmout.f  icnteq.f  iset.f   ivout.f   smout.f  zmout.f 
  cvout.f  dvout.f  icopy.f   iswap.f  second.f  svout.f  zvout.f
  )

esma_add_library (${this}
  SRCS ${SRCS}
)
