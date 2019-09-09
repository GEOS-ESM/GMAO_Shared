esma_set_this()


set(srcs
  dynamics_lattice_module.F mpi_util.F mymalloc.F PWSSSP.F windfix.F
  timer.F gmap.F alias.F fft.F
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL_Base GMAO_pFIO GMAO_hermes)
target_compile_definitions(${this} PRIVATE mpi)


esma_add_library (post_nompi
  SRCS dynamics_lattice_module.F mpi_util.F mymalloc.F timer.F
  DEPENDENCIES MAPL_Base GMAO_pFIO GMAO_hermes)
target_compile_definitions(post_nompi PRIVATE mpi)

foreach (basename
    makeiau hdf2rs gg2fv eta2prs eta2rst rs_hinterp
    time_ave convert_eta rs2hdf
    ec_prs2fv flat2hdf fvrst gg2eta ec_prs2eta
    merra2scm rsg3_vinterp era5_prs2eta)
 ecbuild_add_executable (TARGET ${basename}.x SOURCES ${basename}.F DEFINITIONS mpi)
endforeach ()
ecbuild_add_executable (TARGET binarytile.x SOURCES binarytile.F90)
ecbuild_add_executable (TARGET rs_numtiles.x SOURCES rs_numtiles.F90)
ecbuild_add_executable (TARGET rs_vinterp.x SOURCES rs_vinterp.F90)
ecbuild_add_executable (TARGET rs_vinterp_scm.x SOURCES rs_vinterp_scm.F90)
ecbuild_add_executable (TARGET stats.x SOURCES stats.F90)
if (DISABLE_FIELD_WIDTH_WARNING)
  set_target_properties (stats.x PROPERTIES COMPILE_FLAGS ${DISABLE_FIELD_WIDTH_WARNING})
endif ()
ecbuild_add_executable (TARGET convert_aerosols.x SOURCES convert_aerosols.F)

file(GLOB perlscripts *.pl)
install(
   PROGRAMS ${perlscripts} ec2grd.csh
   DESTINATION bin)

set_property(SOURCE rsg3_vinterp.F APPEND_STRING PROPERTY COMPILE_FLAGS "${FREAL8} ${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE rs_vinterp.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FREAL8} ${BYTERECLEN}")
set_property(SOURCE rs_vinterp_scm.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FREAL8} ${BYTERECLEN}")

set_property(SOURCE dynamics_lattice_module.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE mpi_util.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE mymalloc.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE windfix.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE timer.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE gmap.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE fft.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE convert_eta.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE ec_prs2eta.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE era5_prs2eta.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE ec_prs2fv.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE eta2prs.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE eta2rst.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE flat2hdf.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE fvrst.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE gg2eta.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE gg2fv.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE hdf2rs.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE makeiau.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE merra2scm.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE rs2hdf.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE rs_hinterp.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE time_ave.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE convert_aerosols.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN} ${EXTENDED_SOURCE}")
set_property(SOURCE PWSSSP.F APPEND_STRING PROPERTY COMPILE_FLAGS "${FREAL8} ${EXTENDED_SOURCE} ${BYTERECLEN}")
set_property(SOURCE alias.F APPEND_STRING PROPERTY COMPILE_FLAGS "${BYTERECLEN}")

target_link_libraries (convert_eta.x ${this} GMAO_hermes MAPL_Base)
target_link_libraries (eta2prs.x ${this} GMAO_gfio_r4 MAPL_Base GEOS_Shared)
target_link_libraries (ec_prs2eta.x ${this} GMAO_gfio_r4 MAPL_Base)
target_link_libraries (era5_prs2eta.x ${this} GMAO_gfio_r4 MAPL_Base)
target_link_libraries (ec_prs2fv.x ${this} GMAO_gfio_r4 MAPL_Base)
target_link_libraries (stats.x post_nompi GMAO_gfio_r4 MAPL_Base ${ESMF_LIBRARIES})
target_link_libraries (flat2hdf.x post_nompi ${this} GMAO_gfio_r4)
target_link_libraries (gg2fv.x post)
target_link_libraries (rs2hdf.x post)
target_link_libraries (time_ave.x ${this} GMAO_gfio_r4 MAPL_Base )

target_link_libraries (rsg3_vinterp.x ${this} GMAO_hermes)
target_link_libraries (rs_vinterp.x ${this} GMAO_hermes)
target_link_libraries (eta2rst.x ${this} GMAO_hermes ${MPI_Fortran_LIBRARIES})
target_link_libraries (fvrst.x GMAO_gfio_r4 MAPL_Base GMAO_pFIO)
set_target_properties (fvrst.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
target_link_libraries (gg2eta.x ${this} GMAO_gfio_r4)
target_link_libraries (gg2fv.x ${this} GMAO_gfio_r4)
target_link_libraries (hdf2rs.x ${this} GMAO_gfio_r4 ${MPI_Fortran_LIBRARIES})
target_link_libraries (makeiau.x ${this} GMAO_gfio_r4 ${MPI_Fortran_LIBRARIES})
target_link_libraries (merra2scm.x ${this} GMAO_gfio_r4 ${MPI_Fortran_LIBRARIES})
target_link_libraries (rs2hdf.x ${this} GMAO_gfio_r4)
target_link_libraries (rs_hinterp.x ${this} MAPL_Base GMAO_pFIO)
set_target_properties (rs_hinterp.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

target_link_libraries (rs_numtiles.x MAPL_Base GMAO_pFIO)
set_target_properties (rs_numtiles.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

# Unit tests?
#add_executable (checkdate.x checkdate.F)
#add_executable (swapendian_FV.x swapendian_FV.F90)
#add_executable (swapendian_RST.x swapendian_RST.F90)
#target_link_libraries (swapendian_FV.x)
#target_link_libraries (swapendian_RST.x)

target_link_libraries (rs_vinterp_scm.x ${this} GMAO_hermes)
target_link_libraries (convert_aerosols.x ${this} GMAO_gfio_r8 MAPL_cfio_r8 GMAO_hermes)


include_directories (${include_MAPL_Base})
include_directories (${include_GMAO_hermes})
include_directories (${include_GEOS_Shared})
include_directories (${include_GMAO_pFIO})

set (post_progs
   gcmpost.script
   tick
   get_streams
   get_archive
   getdates
   gcmclim.script
   check_clim
   )

install (
   PROGRAMS ${post_progs}
   DESTINATION post
   )

set (post_rc
   plot.rc
   post.rc
   stats.rc
   time_ave.rc
   )

install (
   FILES ${post_rc}
   DESTINATION post
   )

install (
   FILES stats.rc time_ave.rc
   DESTINATION etc
   )
