esma_check_if_debug()

option(HERMES_LIGHT "Build light-weighted GMAO_shared library" OFF)

esma_add_subdirectories(GEOS_Shared)
esma_add_subdirectories(GMAO_hermes)

if ( NOT HERMES_LIGHT )

   esma_add_subdirectories(
     GMAO_mpeu
     GMAO_pilgrim
     GMAO_etc
     GEOS_Util 
     LANL_Shared
     GMAO_perllib
     GMAO_transf
     GMAO_stoch
   
     GEOS_Pert
     GMAO_mfhdf3
     GMAO_gfioeos
     GMAO_ods
     GMAO_psas
     GMAO_gems
     GMAO_iret
     GMAO_pyobs
     GMAO_radmon
   
     GMAO_ncdiag
     arpack
     pnagpack
     )
   
   # Special case - GMAO_gfio is built twice with two different precisions.
   if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GMAO_gfio)
     add_subdirectory (GMAO_gfio GMAO_gfio_r4)
     add_subdirectory (GMAO_gfio GMAO_gfio_r8)
     add_dependencies (GMAO_gfio_r4 GMAO_gfio_r8)
   endif ()
   
endif()
