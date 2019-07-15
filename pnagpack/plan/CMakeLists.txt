esma_set_this (OVERRIDE planso)

set (SRCS clanc.f  ddotmpi.f plandr2.f  plandr.f  planso.f  ppurge.f  pstartv.f  pstpone.f)

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${FREAL8} ${BIG_ENDIAN} ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

esma_add_library(${this} SRCS ${SRCS})
set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
