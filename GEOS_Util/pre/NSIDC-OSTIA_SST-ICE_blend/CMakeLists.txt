esma_set_this()

set (srcs
  bin2bin.F90
  check.F90
  fill_Land.F90
  hflip.F90
  interp_to_eight_deg.F90
  read_input.F90
  read_input_quart.F90
  read_Ostia.F90
  read_Ostia_quart.F90
  read_Reynolds.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)
target_link_libraries (${this} PUBLIC ${NETCDF_LIBRARIES})
target_include_directories (${this} PUBLIC ${INC_NETCDF})

ecbuild_add_executable (
  TARGET regrid_forcing.x
  SOURCES regrid_forcing.F90
  LIBS ${this})
set_target_properties(regrid_forcing.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})

ecbuild_add_executable (
  TARGET sst_sic_EIGTHdeg.x
  SOURCES proc_SST_FRACI.F90
  LIBS ${this})
set_target_properties(sst_sic_EIGTHdeg.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})

ecbuild_add_executable (
  TARGET sst_sic_QUARTdeg.x
  SOURCES proc_SST_FRACI_ostia_quart.F90
  LIBS ${this})
set_target_properties(sst_sic_QUARTdeg.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})

file(GLOB cshscripts *.csh)
file(GLOB perlscripts *.pl)
install(PROGRAMS ${cshscripts} ${perlscripts} DESTINATION bin)
