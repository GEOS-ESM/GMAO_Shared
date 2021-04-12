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

# This is equivalent to FOPT=$(FOPT2) in GNU Make
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   string (REPLACE "${FOPT3}" "${FOPT2} ${BIG_ENDIAN}" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${OPTREPORT0}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${FTZ}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${ALIGN_ALL}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${NO_ALIAS}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
endif ()
esma_add_library (${this}
  SRCS ${SRCS}
  )

target_include_directories (${this} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
