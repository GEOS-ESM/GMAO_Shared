esma_set_this()

add_definitions (-Dopenmp)

# NOTE: The concommitant GNUmakefile for this directory says we should
# use the r8i4 versions of the NCEP libraries. Unfortunately, if you do
# this you break the use of NCEP_sp in GEOSmkiau_GridComp which requires
# the r4i4 interface. Until that is converted to r8, this routine must
# use r4i4 to get zero-diff with the GNU Make build
# The "correct" lines should be:
#
#include_directories (${include_NCEP_sp_r4i4})
#include_directories (${include_NCEP_w3_r8i4})
#esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu NCEP_sp_r4i4 NCEP_w3_r8i4 NCEP_sfcio GMAO_hermes NCEP_sigio)

include_directories (${include_GMAO_mpeu})
include_directories (${include_NCEP_sp_r4i4})
include_directories (${include_NCEP_w3_r8i4})
include_directories (${include_NCEP_sfcio})
include_directories (${include_GMAO_hermes})
include_directories (${include_NCEP_sigio})
include_directories (${include_${this}})

#target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>)

set (srcs
   atod.f
   blend_uvsh.F90
   dtoa.F
   dynamics_lattice_module.F
   m_elliptic.f
   fft.F
   fv2fv.F90
   fv2ll.F
   gauss_lat_nmc.f
   get_ncep_phis.f
   get_ncep_tsx.f
   m_ec2fv.F
   m_hinterp.F
   interp_h.F
   ll2fv.F
   m_InterleaveScatterer.F90
   m_InterleaveScattererComm.F90
   m_Interleaved.F90
   m_InterleavedComm.F90
   m_SubdomainDistributor.F90
   m_SubdomainDistributorComm.F90
   m_checksums.F90
   m_daInterp.F90
   m_fgInterp.F90
   m_fvGrid.F90
   m_fvGridHeader.F90
   m_fvGridThickness.F90
   m_fvGriddedState.F90
   m_fvSurface.F90
   m_geosapi.F90
   m_gg.F90
   m_ggGradient.F90
   m_ggGradientSP.F90
   m_gsiGriddedState.F90
   m_interleavedObject.F90
   m_llInterp.F90
   m_massadj.F90
   m_ncSurface.F90
   m_ncepphis.F90
   m_parGrADS.F90
   m_ppInterp.F90
   m_simpleGrid.F90
   m_ss.f90
   m_swapij.F90
   m_undef.F90
   m_xform.f90
   make_surf.f90
   mappm.F
   mpi_util.F
   mymalloc.F
   qsat.f
   set_regrid.f90
   tick.f
   timer.F
   util.f90
   gmap.F
   )

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM} ${MISMATCH}")
set (CMAKE_Fortran_FLAGS_DEBUG   "${GEOS_Fortran_FLAGS_DEBUG}")

ecbuild_add_executable(TARGET dyn_divo.x SOURCES dyn_divo.F90 LIBS ${this} FFLAGS "${FREAL8}")
ecbuild_add_executable(TARGET ec2fv.x SOURCES ec2fv.F90 LIBS ${this} FFLAGS "${FREAL8}")
ecbuild_add_executable(TARGET fv2ss.x SOURCES fv2ss.f90 LIBS ${this} FFLAGS "${FREAL8}")
ecbuild_add_executable(TARGET g4tog5.x SOURCES g4tog5.f90 LIBS ${this} FFLAGS "${FREAL8}")
ecbuild_add_executable(TARGET ss2fv.x SOURCES ss2fv.f90 LIBS ${this} FFLAGS "${FREAL8}")

install (PROGRAMS ec2fv.csh DESTINATION bin)

#This needs to be worked on
ecbuild_add_executable(TARGET ss2gg.x SOURCES ss2gg.f90 LIBS NCEP_w3_r4i4 NCEP_sp_r4i4 NCEP_bufr_r4i4 NCEP_bacio_r4i4 NCEP_sigio)

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu NCEP_sp_r4i4 NCEP_bufr_r4i4  NCEP_w3_r8i4 NCEP_sfcio GMAO_hermes NCEP_sigio ${MKL_LIBRARIES})
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})
if (EXTENDED_SOURCE)
   set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()

string (REPLACE " " ";" flags ${FREAL8})
foreach (flag ${flags})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()
