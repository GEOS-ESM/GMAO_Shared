esma_set_this (OVERRIDE lanso)

set (SRCS dsort2.f  lanso.f   machar.f  ortbnd.f  pythag.f  ritvec.f  stpone.f  tqlb.f 
    datx.f  landr.f   loadev.f   purge.f   random.f  startv.f  tql2.f store.f deteps.f)

esma_add_library(${this} SRCS ${SRCS})
