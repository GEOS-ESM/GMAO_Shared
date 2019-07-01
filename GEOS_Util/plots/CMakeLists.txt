set (zonal_prog_name
   zonal_${CMAKE_SYSTEM_NAME}.x
   )

ecbuild_add_executable (TARGET ${zonal_prog_name} SOURCES zonal.f LIBS MAPL_Base)
if (EXTENDED_SOURCE)
  set_source_files_properties(zonal.f PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()

set (plots_progs
   quickplot
   quickstat
   configure
   landscape.script
   portrait.script
   moveplot
   movestat
   get_exports
   )

install (
   PROGRAMS ${plots_progs}
   DESTINATION plots
   )
