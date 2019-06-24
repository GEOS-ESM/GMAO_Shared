esma_set_this()

add_definitions (-DHERMES -DDEBUG_SHAVE -Dunix -D__unix__)
include_directories (${include_MAPL_Base} ${include_GEOS_Shared})
include_directories (${include_GMAO_gfio_r8})

string (REPLACE " " ";" flags ${FREAL8})
add_compile_options(${flags})

set (srcs
  m_ana2dyn.F m_const.f90 m_dyn.f90 m_dyn2dyn.f90 m_dynp.F90
  m_insitu.F m_interp.F90 m_lp.F m_maph.F90 m_maph_pert.F90 m_mapz.F m_mapz_pert.F90
  m_qsat.F m_remap.f90 m_set_eta.F90 m_vdc2vdc.F90 progsig.f90 m_topo_remap.F90
  m_interpack.F m_interpack_ad.F m_ec_set_eta.F90 m_ecdyn.f90 m_tick.f90
  m_dyn_util.F90 m_spline.f90
  )

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu GMAO_gfio_r8)
target_compile_options(${this} PRIVATE ${flags})

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM}")

if (EXTENDED_SOURCE)
  set_source_files_properties (m_interpack.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_interpack_ad.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_ana2dyn.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_insitu.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_set_eta.F90 PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()
set_source_files_properties (GFIO_Shave.f90 PROPERTIES COMPILE_FLAGS ${PP})
set_source_files_properties (GFIO_remap.f90 PROPERTIES COMPILE_FLAGS ${PP})

ecbuild_add_executable (TARGET dyndims.x SOURCES dyndims.f LIBS ${this})
ecbuild_add_executable (TARGET GFIO_remap.x SOURCES GFIO_remap.f90 LIBS ${this} MAPL_cfio_r8)
ecbuild_add_executable (TARGET GFIO_Shave.x SOURCES GFIO_Shave.f90 LIBS ${this} MAPL_cfio_r8)
ecbuild_add_executable (TARGET echorc.x SOURCES echorc.f90 LIBS ${this})
ecbuild_add_executable (TARGET fvsubset.x SOURCES fvsubset.F90 LIBS ${this} MAPL_cfio_r8)




