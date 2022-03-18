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

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL NetCDF::NetCDF_Fortran)

ecbuild_add_executable (
  TARGET regrid_forcing.x
  SOURCES regrid_forcing.F90
  LIBS ${this})
set_target_properties(regrid_forcing.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

ecbuild_add_executable (
  TARGET regrid_forcing_esmf.x
  SOURCES regrid_forcing_esmf.F90
  LIBS ${this})
set_target_properties(regrid_forcing_esmf.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

ecbuild_add_executable (
  TARGET sst_sic_EIGTHdeg.x
  SOURCES proc_SST_FRACI.F90
  LIBS ${this})
set_target_properties(sst_sic_EIGTHdeg.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

ecbuild_add_executable (
  TARGET sst_sic_QUARTdeg.x
  SOURCES proc_SST_FRACI_ostia_quart.F90
  LIBS ${this})
set_target_properties(sst_sic_QUARTdeg.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

ecbuild_add_executable (
  TARGET lake_sst_sic_EIGTHdeg.x
  SOURCES lake_data_EIGTHdeg.F90
  LIBS ${this})
set_target_properties(lake_sst_sic_EIGTHdeg.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

ecbuild_add_executable (
  TARGET reynolds_QUART.x
  SOURCES proc_SST_FRACI_reynolds_quart.F90
  LIBS ${this})
set_target_properties(reynolds_QUART.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

if (USE_F2PY)
   find_package(F2PY2)
   if (F2PY2_FOUND)
      esma_add_f2py2_module(read_ops_bcs
         SOURCES read_bin.f90
         DESTINATION bin
         INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
         )
      add_dependencies(read_ops_bcs ${this})
   endif ()
endif ()

file(GLOB cshscripts CONFIGURE_DEPENDS *.csh)
file(GLOB perlscripts CONFIGURE_DEPENDS *.pl)
install(PROGRAMS ${cshscripts} ${perlscripts} DESTINATION bin)
