
file(GLOB test_scripts CONFIGURE_DEPENDS ./*.*)
install(
   PROGRAMS ${test_scripts} 
   DESTINATION bin)
