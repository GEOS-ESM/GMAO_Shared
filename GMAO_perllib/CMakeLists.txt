set(perldirs
   ESMA
   GrADS
   Net
   Time
   XML
   )

foreach(perldir ${perldirs})
   install (DIRECTORY ${perldirs}
      DESTINATION bin
      )
endforeach()
