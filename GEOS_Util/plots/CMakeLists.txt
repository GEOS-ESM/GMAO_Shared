ecbuild_add_executable (TARGET zonal.x SOURCES zonal.f LIBS MAPL_Base)
if (EXTENDED_SOURCE)
  set_source_files_properties(zonal.f PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quickplot ${CMAKE_CURRENT_BINARY_DIR}/quickplot COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure ${CMAKE_CURRENT_BINARY_DIR}/configure COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/landscape.script ${CMAKE_CURRENT_BINARY_DIR}/landscape.script COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/portrait.script ${CMAKE_CURRENT_BINARY_DIR}/portrait.script COPYONLY)
