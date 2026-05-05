esma_set_this()

set (srcs
  GEOS_Utilities.F90 DragCoefficients.F90 GEOS.F90
  tropovars.F90 surfacelayer.F90 GEOS_TopoGet.F90
  GEOS_GmapMod.F90 GEOS_RemapMod.F90 G3_MPI_Util_Mod.F90 windfix.F90
  PWSSSP.F sphere.F calcdbz.F OVP.F90
  getco2.F90 atmOceanIntLayer.F90
  )

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MAPL ESMF::ESMF NetCDF::NetCDF_Fortran)

# special cases
set_source_files_properties(sphere.F PROPERTIES COMPILE_FLAGS "${FREAL8}")

if (EXTENDED_SOURCE)
   set_source_files_properties (PWSSSP.F PROPERTIES COMPILE_FLAGS "${FREAL8} ${EXTENDED_SOURCE}")
else ()
   set_source_files_properties (PWSSSP.F PROPERTIES COMPILE_FLAGS "${FREAL8}")
endif ()



