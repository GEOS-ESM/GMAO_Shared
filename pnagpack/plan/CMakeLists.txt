esma_set_this (OVERRIDE planso)

set (SRCS clanc.f  ddotmpi.f plandr2.f  plandr.f  planso.f  ppurge.f  pstartv.f  pstpone.f)
esma_add_library(${this} SRCS ${SRCS})
set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
