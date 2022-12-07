esma_set_this()

# add a dummy target so that GMAO_etc can depend
# on GMAO_perllib

add_custom_target(${this} ALL)

set(perldirs
   ESMA
   GrADS
   Net
   Time
   XML
   )

foreach(perldir ${perldirs})
   install (
      DIRECTORY ${perldirs}
      DESTINATION bin
      USE_SOURCE_PERMISSIONS
      )
endforeach()
