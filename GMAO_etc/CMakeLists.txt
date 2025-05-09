esma_set_this()

# Add a dummy custom target because GMAO_etc depends
# on Perl modules in GMAO_perllib

add_custom_target(${this} ALL)
add_dependencies(${this} GMAO_perllib)

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

ecbuild_add_executable (
  TARGET argopt
  SOURCES argopt.c
  )

ecbuild_add_executable (
  TARGET fresolv.x
  SOURCES fresolv.f90
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET gethdftinc.x
  SOURCES gethdftinc.f
  LIBS GMAO_mpeu NetCDF::NetCDF_Fortran Threads::Threads
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_ci.x
  SOURCES zeit_ci.f zeit_usage.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_co.x
  SOURCES zeit_co.f zeit_usage.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_pr.x
  SOURCES zeit_pr.f zeit_usage.f zeit_read.f zeit_jday.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET mkdrstdate.x
  SOURCES mkdrstdate.f
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET rst_date
  SOURCES rst_date.f
  FFLAGS ${FREAL8})

##### Installation Section #####

## bin ##

file(GLOB perl_files CONFIGURE_DEPENDS *.pl)
file(GLOB perl_packages CONFIGURE_DEPENDS *.pm)
file(GLOB python_files CONFIGURE_DEPENDS *.py)

set(perlscripts
   esma_mpirun
   g5_modules_perl_wrapper
   jobIDfilter
   lsfpL
   mkdir_remote
   n4zip
   now
   num_days_in_month
   pasta
   pesto
   prund.pl
   r_dist
   require_obsys
   rndasfiles
   tick
   token_resolve
   vED
   stripname
   )

set(otherbin
   Functions.csh
   pyrob
   pyrob_CF
   RmShmKeys_sshmpi.csh
   rmshmkeyhere.sh
   protex
   )

install (
   PROGRAMS ${perl_files} ${perl_packages} ${python_files} ${perlscripts} ${otherbin}
   DESTINATION bin
   )

set(perlscripts_with_dasperl
   bacon
   fvgetf
   grads2grib
   rget
   rst4ncepf
   rst2rcyc
   rput
   )

set (DASPERL /usr/bin/perl)
foreach (script ${perlscripts_with_dasperl})
   configure_file(${script} ${script} @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin)
endforeach ()

## etc ##

file(GLOB obsys_rc CONFIGURE_DEPENDS obsys*.rc)
file(GLOB templates CONFIGURE_DEPENDS *.tmpl)

install (
   FILES ${obsys_rc} ${templates} pesto.arc PUBLICTAG
   DESTINATION etc
   )

## grads ##

set (grads_files
   lats4d.gs
   lats4d
   info.gs
   )

install (PROGRAMS ${grads_files}
   DESTINATION lib/grads
   )
