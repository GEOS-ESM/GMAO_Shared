
file(GLOB test_scripts CONFIGURE_DEPENDS ./*.py ./*.yaml)
install(
   PROGRAMS ${test_scripts} 
   DESTINATION bin)
