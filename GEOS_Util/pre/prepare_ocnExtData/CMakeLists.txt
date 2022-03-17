esma_set_this()

set ( srcs
  sst_ice_helpers.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL NetCDF::NetCDF_Fortran)

set (pyscripts
  plot_binary_sst_ice.py
  plot_sst_ice.py
  ocnExtData_utils.py
  write_kPAR.py
  )

ecbuild_add_executable (
  TARGET sst_sic_eight.x
  SOURCES proc_SST_FRACI_eight.F90
  LIBS ${this})
set_target_properties(sst_sic_eight.x PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

if (USE_F2PY)
   find_package(F2PY2)
   if (F2PY2_FOUND)
      esma_add_f2py2_module(read_sst_ice_bcs
         SOURCES read_sst_ice_bin.f90
         DESTINATION bin
         INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
         )
      add_dependencies(read_sst_ice_bcs ${this})
   endif ()
endif ()

file(GLOB cshscripts CONFIGURE_DEPENDS *.csh)
file(GLOB perlscripts CONFIGURE_DEPENDS *.pl)
install(PROGRAMS ${cshscripts} ${perlscripts} ${pyscripts} DESTINATION bin)
