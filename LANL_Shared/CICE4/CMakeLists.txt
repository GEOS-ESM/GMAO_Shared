esma_set_this()
get_filename_component(this "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
set (srcs
  csm_share/shr_orb_mod.F90
  drivers/cice4/ice_constants.F90
  drivers/cice4/ice_prescribed_mod.F90
  mpi/ice_boundary.F90
  mpi/ice_broadcast.F90
  mpi/ice_communicate.F90
  mpi/ice_exit.F90
  mpi/ice_gather_scatter.F90
  mpi/ice_global_reductions.F90
  mpi/ice_timers.F90
  source/ice_age.F90
  source/ice_atmo.F90
  source/ice_blocks.F90
  source/ice_calendar.F90
  source/ice_diagnostics.F90
  source/ice_distribution.F90
  source/ice_domain.F90
  source/ice_domain_size.F90
  source/ice_dyn_evp.F90
  source/ice_fileunits.F90
  source/ice_flux.F90
  source/ice_forcing.F90
  source/ice_grid.F90
  source/ice_history.F90
  source/ice_init.F90
  source/ice_itd.F90
  source/ice_kinds_mod.F90
  source/ice_mechred.F90
  source/ice_meltpond.F90
  source/ice_ocean.F90
  source/ice_orbital.F90
  source/ice_read_write.F90
  source/ice_restart.F90
  source/ice_restoring.F90
  source/ice_shortwave.F90
  source/ice_spacecurve.F90
  source/ice_state.F90
  source/ice_step_mod.F90
  source/ice_therm_itd.F90
  source/ice_therm_vertical.F90
  source/ice_transport_driver.F90
  source/ice_transport_remap.F90
  source/ice_work.F90
  )

esma_add_library(${this}
                 SRCS ${srcs}
                 DEPENDENCIES MAPL
                 TYPE SHARED
                 )

target_compile_definitions (${this} PRIVATE coupled GEOS MODIFY_TOPOGRAPHY USE_R8)

if (CMAKE_SYSTEM_NAME MATCHES Linux)
   target_compile_definitions(${this} PRIVATE LINUX)
endif ()

# It turns out the new vectorized Release flags for Intel cause a segfault
# in CICE4. For now, we workaround this by restoring the old non-vectorized
# flags here when building as Release
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Release)
  set (CMAKE_Fortran_FLAGS_RELEASE  "${GEOS_Fortran_FLAGS_NOVECT}")
endif ()

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
if ( BUILD_MIT_OCEAN )
  add_compile_definitions(BUILD_MIT_OCEAN)
endif ()
  
ecbuild_add_executable(TARGET makdep SOURCES bld/makdep.c)
