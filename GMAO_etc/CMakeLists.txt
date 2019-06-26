set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

ecbuild_add_executable (
  TARGET argopt
  SOURCES argopt.c)

ecbuild_add_executable (
  TARGET fresolv.x
  SOURCES fresolv.f90
  LIBS GMAO_mpeu)

ecbuild_add_executable (
  TARGET gethdftinc.x
  SOURCES gethdftinc.f
  LIBS GMAO_mpeu ${NETCDF_LIBRARIES} Threads::Threads
  INCLUDES ${INC_NETCDF})

ecbuild_add_executable (
  TARGET zeit_ci.x
  SOURCES zeit_ci.f zeit_usage.f
  LIBS GMAO_mpeu)

ecbuild_add_executable (
  TARGET zeit_co.x
  SOURCES zeit_co.f zeit_usage.f
  LIBS GMAO_mpeu)

ecbuild_add_executable (
  TARGET zeit_pr.x
  SOURCES zeit_pr.f zeit_usage.f zeit_read.f zeit_jday.f
  LIBS GMAO_mpeu)

ecbuild_add_executable (
  TARGET mkdrstdate.x
  SOURCES mkdrstdate.f)

ecbuild_add_executable (
  TARGET rst_date
  SOURCES rst_date.f)

##### Installation Section #####

## bin ##

file(GLOB perl_files *.pl)
file(GLOB perl_packages *.pm)
file(GLOB python_files *.py)

set(perlscripts
   acquire
   acquire_obsys
   bacon
   esma_mpirun
   fvgetf
   g5_modules_perl_wrapper
   grads2grib
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
   rget
   rndasfiles
   rput
   rst2rcyc
   rst4ncepf
   tick
   token_resolve
   vED
   stripname
   )

set(otherbin
   Functions.csh
   pyrob
   RmShmKeys_sshmpi.csh
   rmshmkeyhere.sh
   )

install (
   FILES ${perl_files} ${perl_packages} ${python_files} ${perlscripts} ${otherbin}
   DESTINATION bin
   )

## etc ##

file(GLOB obsys_rc obsys*.rc)
file(GLOB templates *.tmpl)

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

install (FILES ${grads_files}
   DESTINATION lib/grads
   )
