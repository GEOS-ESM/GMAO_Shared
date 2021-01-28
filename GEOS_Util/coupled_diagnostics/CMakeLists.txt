# From https://stackoverflow.com/questions/7787823/cmake-how-to-get-the-name-of-all-subdirectories-of-a-directory

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children CONFIGURE_DEPENDS RELATIVE ${curdir} ${curdir}/*)
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
   DESTINATION coupled_diagnostics
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
