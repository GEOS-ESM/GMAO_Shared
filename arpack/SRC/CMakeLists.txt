esma_set_this (OVERRIDE arpack)

set (SRCS
  cgetv0.f  cstatn.f            dsconv.f  dstqrb.f  snconv.f  sseigt.f sngets.f 
  cnaitr.f  dnconv.f  dseigt.f  sneigh.f  ssesrt.f  zsortc.f 
  cnapps.f  dneigh.f  dsesrt.f  sgetv0.f  sneupd.f  sseupd.f  zgetv0.f   zstatn.f 
  cnaup2.f  dgetv0.f  dneupd.f  dseupd.f  slaqrb.f  sngets.f  ssgets.f  znaitr.f 
  cnaupd.f  dlaqrb.f  dngets.f  dsgets.f  snaitr.f  ssaitr.f  ssortc.f  znapps.f 
  cneigh.f  dnaitr.f  dsaitr.f  dsortc.f  snapps.f  ssapps.f  ssortr.f  znaup2.f 
  cneupd.f  dnapps.f  dsapps.f  dsortr.f  snaup2.f  ssaup2.f  sstatn.f  znaupd.f 
  cngets.f  dnaup2.f  dsaup2.f  dstatn.f  snaupd.f  ssaupd.f  sstats.f  zneigh.f 
  csortc.f  dnaupd.f  dsaupd.f  dstats.f            ssconv.f  sstqrb.f  zneupd.f 
  sspg.f    dspg.f
  )

esma_add_library (${this}
  SRCS ${SRCS}
  )


target_include_directories (${this} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
