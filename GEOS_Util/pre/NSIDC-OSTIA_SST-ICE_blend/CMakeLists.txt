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

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL)
target_link_libraries (${this} PUBLIC ${NETCDF_LIBRARIES})
target_include_directories (${this} PUBLIC ${INC_NETCDF})

ecbuild_add_executable (
  TARGET regrid_forcing.x
  SOURCES regrid_forcing.F90
  LIBS ${this})

ecbuild_add_executable (
  TARGET regrid_forcing_esmf.x
  SOURCES regrid_forcing_esmf.F90
  LIBS ${this})

ecbuild_add_executable (
  TARGET sst_sic_EIGTHdeg.x
  SOURCES proc_SST_FRACI.F90
  LIBS ${this})

ecbuild_add_executable (
  TARGET sst_sic_QUARTdeg.x
  SOURCES proc_SST_FRACI_ostia_quart.F90
  LIBS ${this})

ecbuild_add_executable (
  TARGET lake_sst_sic_EIGTHdeg.x
  SOURCES lake_data_EIGTHdeg.F90
  LIBS ${this})

if (F2PY_FOUND)
   add_f2py_module(read_ops_bcs
      SOURCES read_bin.f90
      DESTINATION bin
      INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
      )
   add_dependencies(read_ops_bcs ${this})
endif (F2PY_FOUND)

file(GLOB cshscripts CONFIGURE_DEPENDS *.csh)
file(GLOB perlscripts CONFIGURE_DEPENDS *.pl)
install(PROGRAMS ${cshscripts} ${perlscripts} DESTINATION bin)
