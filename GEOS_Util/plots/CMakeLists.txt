set (zonal_prog_name
   zonal_${CMAKE_SYSTEM_NAME}.x
   )

ecbuild_add_executable (TARGET ${zonal_prog_name} SOURCES zonal.f LIBS MAPL_Base)
if (EXTENDED_SOURCE)
  set_source_files_properties(zonal.f PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()

set (regrid2_prog_name
   regrid2_${CMAKE_SYSTEM_NAME}.x
   )

ecbuild_add_executable (TARGET ${regrid2_prog_name} SOURCES grads_util/regrid2.f)
if (EXTENDED_SOURCE)
  set_source_files_properties(regrid2.f PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
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
   chckhist.new
   chckhist.old
	chckrc
   zonal
   )

install (
   PROGRAMS ${plots_progs}
   DESTINATION plots
   )

# From https://stackoverflow.com/questions/7787823/cmake-how-to-get-the-name-of-all-subdirectories-of-a-directory

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

# Using MESSAGE_NEVER to not overwhelm the installation
# messages with all these installs
install (
   DIRECTORY ${SUBDIRS}
   DESTINATION plots
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
