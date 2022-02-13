esma_set_this()

#set ( srcs
#  )

#esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL)
#target_link_libraries (${this} PUBLIC ${NETCDF_LIBRARIES})
#target_include_directories (${this} PUBLIC ${INC_NETCDF})


file(GLOB cshscripts CONFIGURE_DEPENDS *.csh)
file(GLOB perlscripts CONFIGURE_DEPENDS *.pl)
install(PROGRAMS ${cshscripts} ${perlscripts} DESTINATION bin)
