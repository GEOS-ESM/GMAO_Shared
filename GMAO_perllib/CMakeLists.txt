esma_set_this()

# add a dummy target so that GMAO_etc can depend
# on GMAO_perllib

add_custom_target(${this} ALL)

file(GLOB perl_files CONFIGURE_DEPENDS *.pl)
# We must remove two files because of DASPERL
list(REMOVE_ITEM perl_files
   ${CMAKE_CURRENT_SOURCE_DIR}/Err_Log.pl
   ${CMAKE_CURRENT_SOURCE_DIR}/fndate.pl
   )

file(GLOB perl_packages CONFIGURE_DEPENDS *.pm)

set(perlscripts
   acquire
   acquire_obsys
   )

install (
   PROGRAMS ${perl_files} ${perl_packages} ${perlscripts}
   DESTINATION bin
   )

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

set(perlscripts_with_dasperl
   Err_Log.pl
   fndate.pl
   )
   
set (DASPERL /usr/bin/perl)
foreach (script ${perlscripts_with_dasperl})
   configure_file(${script} ${script} @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin)
endforeach ()
