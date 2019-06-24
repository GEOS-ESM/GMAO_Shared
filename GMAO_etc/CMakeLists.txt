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


set(exec_scripts
    Manipulate_time.pm
    tick
    stripname
    RmShmKeys_sshmpi.csh
    rmshmkeyhere.sh
    esma_mpirun
    GMAO_utils.pm
   )

foreach(exec ${exec_scripts})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${exec} ${CMAKE_BINARY_DIR}/bin/${exec} COPYONLY)
endforeach()
