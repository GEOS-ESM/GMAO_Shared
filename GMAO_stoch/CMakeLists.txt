esma_set_this()

# NOTE: The concommitant GNUmakefile for this directory says we should
# use the r8i4 versions of the NCEP libraries. Unfortunately, if you do
# this you break the use of NCEP_sp in GEOSmkiau_GridComp which requires
# the r4i4 interface. Until that is converted to r8, this routine must
# use r4i4 to get zero-diff with the GNU Make build.
# The "correct" lines should be:
#
#include_directories (${include_NCEP_sp_r8i4})
#target_link_libraries (${this} GMAO_mpeu GMAO_gfio_r8 GMAO_transf MAPL NCEP_sfcio NCEP_sp_r8i4 MAPL_cfio_r4)

include_directories (${include_NCEP_sigio})
include_directories (${include_NCEP_sfcio})
include_directories (${include_GMAO_hermes})
include_directories (${include_GMAO_transf})
include_directories (${include_GMAO_mpeu})
include_directories (${INC_NETCDF})
include_directories (${INC_ESMF})
include_directories (${include_NCEP_sp_r4i4})
include_directories (${include_${this}})


set (srcs
   epslon.f90
   gauss2grd.f90
   get_ls_node.f
   get_pattern.f90
   glats.f
   grd_xform.F90
   gozrineo.f
   initialize.F90
   mersenne_twister.f
   mod_four.F90
   mod_uvdz.F90
   mod_param.F90
   mpi_more.f
   num_parthd.f90
   patterngenerator.f90
   pln2eo.f
   setlats.f
   stochmain.f90
   stoch_data.f90
   stoch_module.F90
   wrap_fftpack.f
   )

 esma_add_library(${this}
   SRCS ${srcs}
   DEPENDENCIES GMAO_mpeu GMAO_gfio_r8 GMAO_transf MAPL NCEP_sfcio NCEP_sp_r4i4
   )
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

string (REPLACE " " ";" flags ${FREAL8})
foreach (flag ${flags})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()

set (CMAKE_Fortran_FLAGS_RELEASE "-O1 ${BIG_ENDIAN} ${common_Fortran_flags} ${EXTENDED_SOURCE} ${TRACEBACK} ${ALIGNCOM}")

