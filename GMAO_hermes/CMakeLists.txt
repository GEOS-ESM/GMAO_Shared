esma_set_this()

add_definitions (-DHERMES -DDEBUG_SHAVE -Dunix -D__unix__)
string (REPLACE " " ";" flags ${FREAL8})
add_compile_options(${flags})
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM} ${MISMATCH}")

if ( NOT HERMES_LIGHT )

   include_directories (${include_GEOS_Shared})
   include_directories (${include_GMAO_gfio_r8})
   
   set (srcs
     m_ana2dyn.F m_const.f90 m_dyn.f90 m_dyn2dyn.f90 m_dynp.F90
     m_insitu.F m_interp.F90 m_lp.F m_maph.F90 m_maph_pert.F90 m_mapz.F m_mapz_pert.F90
     m_qsat.F m_remap.f90 m_set_eta.F90 m_vdc2vdc.F90 progsig.f90 m_topo_remap.F90
     m_interpack.F m_interpack_ad.F m_ec_set_eta.F90 m_ecdyn.f90 m_tick.f90
     m_dyn_util.F90 m_spline.f90
     )
   
   esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu GMAO_gfio_r8)
   target_compile_options(${this} PRIVATE ${flags})
   
   if (EXTENDED_SOURCE)
     set_source_files_properties (m_interpack.F    PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (m_interpack_ad.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (m_ana2dyn.F      PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (m_insitu.F       PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (m_set_eta.F90    PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (ibc_upd.F        PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (gcmbkg2ana.F     PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
     set_source_files_properties (diag2dyn.F       PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
   endif ()
   set_source_files_properties (GFIO_Shave.f90 PROPERTIES COMPILE_FLAGS ${PP})
   set_source_files_properties (GFIO_remap.f90 PROPERTIES COMPILE_FLAGS ${PP})
   
   ecbuild_add_executable (TARGET GFIO_remap.x     SOURCES GFIO_remap.f90     LIBS ${this} MAPL_cfio_r8)
   ecbuild_add_executable (TARGET GFIO_Shave.x     SOURCES GFIO_Shave.f90     LIBS ${this} MAPL_cfio_r8)
   ecbuild_add_executable (TARGET lcv2prs.x        SOURCES lcv2prs.F90        LIBS ${this} MAPL_cfio_r8 MAPL GEOS_Shared)
   
   ecbuild_add_executable (TARGET ana5sfc.x        SOURCES ana5sfc.F90        LIBS ${this})
   ecbuild_add_executable (TARGET diag2dyn.x       SOURCES diag2dyn.F         LIBS ${this})
   ecbuild_add_executable (TARGET drs2dyn.x        SOURCES drs2dyn.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn2drs.x        SOURCES dyn2drs.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn2dyn.x        SOURCES dyn2dyn.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn_boot.x       SOURCES dyn_boot.f90       LIBS ${this})
   ecbuild_add_executable (TARGET dyn2prs.x        SOURCES dyn2prs.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn2rs5.x        SOURCES dyn2rs5.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn52dyn.x       SOURCES dyn52dyn.f90       LIBS ${this})
   ecbuild_add_executable (TARGET dyndiff.x        SOURCES dyndiff.f90        LIBS ${this} GEOS_Shared)
   ecbuild_add_executable (TARGET dyndims.x        SOURCES dyndims.f          LIBS ${this})
   ecbuild_add_executable (TARGET dynp.x           SOURCES dynp.f90           LIBS ${this})
   ecbuild_add_executable (TARGET echorc.x         SOURCES echorc.f90         LIBS ${this})
   ecbuild_add_executable (TARGET fvsubset.x       SOURCES fvsubset.F90       LIBS ${this})
   ecbuild_add_executable (TARGET getgfiodim.x     SOURCES getgfiodim.f90     LIBS ${this})
   ecbuild_add_executable (TARGET ibc_upd.x        SOURCES ibc_upd.F          LIBS ${this})
   ecbuild_add_executable (TARGET maph.x           SOURCES maph.f90           LIBS ${this})
   ecbuild_add_executable (TARGET maph_pert.x      SOURCES maph_pert.f90      LIBS ${this})
   ecbuild_add_executable (TARGET mapz.x           SOURCES mapz.f             LIBS ${this})
   ecbuild_add_executable (TARGET recalcsfc.x      SOURCES recalcsfc.F90      LIBS ${this})
   ecbuild_add_executable (TARGET reset_time.x     SOURCES reset_time.f       LIBS ${this})
   ecbuild_add_executable (TARGET rs52dyn.x        SOURCES rs52dyn.f90        LIBS ${this})
   ecbuild_add_executable (TARGET gcmbkg2ana.x     SOURCES gcmbkg2ana.F       LIBS ${this})
   ecbuild_add_executable (TARGET dyndot.x         SOURCES dyndot.f90         LIBS ${this})
   ecbuild_add_executable (TARGET dyn_pert_remap.x SOURCES dyn_pert_remap.F90 LIBS ${this})
   ecbuild_add_executable (TARGET dyn_recenter.x   SOURCES dyn_recenter.f90   LIBS ${this})
   ecbuild_add_executable (TARGET dyn_rankhist.x   SOURCES dyn_rankhist.F90   LIBS ${this})
   ecbuild_add_executable (TARGET dyn2real_eta.x   SOURCES dyn2real_eta.f90   LIBS ${this})
   ecbuild_add_executable (TARGET dyn_ncf2dyn.x    SOURCES dyn_ncf2dyn.f90    LIBS ${this})
   ecbuild_add_executable (TARGET dyn_iupd.x       SOURCES dyn_iupd.F90       LIBS ${this})
   ecbuild_add_executable (TARGET dyn_cov.x        SOURCES dyn_cov.f90        LIBS ${this})
   ecbuild_add_executable (TARGET dyn_efsens.x     SOURCES dyn_efsens.f90     LIBS ${this})
   ecbuild_add_executable (TARGET eta_echo.x       SOURCES eta_echo.f90       LIBS ${this})
   ecbuild_add_executable (TARGET dyn_inflate.x    SOURCES dyn_inflate.f90    LIBS ${this})
   
   install(PROGRAMS cnv2prs.pl ana5sfc DESTINATION bin)

else ()
   set (srcs
     m_set_eta.F90 shared_topo_remap.F90
     )
   esma_add_library(${this} SRCS ${srcs} )
   target_compile_options(${this} PRIVATE ${flags})
endif ()

set(plain_resource_files
   ana.prs.rc
   dyn2dyn.rc
   dyn2prs.rc
   fcst_fv2prs.lcv.rc
   fcst_fv2prs.rc
   fv2prs.lcv.rc
   fv2prs.rc
   )

foreach (rcfile ${plain_resource_files})
   install(FILES ${rcfile} DESTINATION etc)
endforeach ()

set(resource_files_needing_processing
   const_2d_asm_Nx.rc
   const_2d_chm_Fx.rc
   inst1_2d_int_Nx.rc
   inst1_2d_int_NxDM.rc
   inst1_2d_int_NxMM.rc
   inst2d_met_x.rc
   inst3_3d_asm_Cp.rc
   inst3_3d_asm_CpDM.rc
   inst3_3d_asm_CpMM.rc
   inst3_3d_chm_Ne.rc
   inst3_3d_chm_Nv.rc
   inst3d_met_p.rc
   inst6_3d_ana_Np.rc
   inst6_3d_ana_NpDM.rc
   inst6_3d_ana_NpMM.rc
   inst6_3d_ana_Nv.rc
   tavg1_2d_flx_Nx.rc
   tavg1_2d_flx_NxDM.rc
   tavg1_2d_flx_NxMM.rc
   tavg1_2d_int_Nx.rc
   tavg1_2d_int_NxDM.rc
   tavg1_2d_int_NxMM.rc
   tavg1_2d_lnd_Nx.rc
   tavg1_2d_lnd_NxDM.rc
   tavg1_2d_lnd_NxMM.rc
   tavg1_2d_rad_Nx.rc
   tavg1_2d_rad_NxDM.rc
   tavg1_2d_rad_NxMM.rc
   tavg1_2d_slv_Nx.rc
   tavg1_2d_slv_NxDM.rc
   tavg1_2d_slv_NxMM.rc
   tavg2d_met_x.rc
   tavg3_2d_chm_Fx.rc
   tavg3_2d_chm_FxDM.rc
   tavg3_2d_chm_FxMM.rc
   tavg3_2d_chm_Nx.rc
   tavg3_3d_chm_Fe.rc
   tavg3_3d_chm_Fv.rc
   tavg3_3d_chm_Ne.rc
   tavg3_3d_chm_Nv.rc
   tavg3_3d_cld_Cp.rc
   tavg3_3d_cld_CpDM.rc
   tavg3_3d_cld_CpMM.rc
   tavg3_3d_mst_Cp.rc
   tavg3_3d_mst_CpDM.rc
   tavg3_3d_mst_CpMM.rc
   tavg3_3d_odt_Cp.rc
   tavg3_3d_odt_CpDM.rc
   tavg3_3d_odt_CpMM.rc
   tavg3_3d_qdt_Cp.rc
   tavg3_3d_qdt_CpDM.rc
   tavg3_3d_qdt_CpMM.rc
   tavg3_3d_rad_Cp.rc
   tavg3_3d_rad_CpDM.rc
   tavg3_3d_rad_CpMM.rc
   tavg3_3d_tdt_Cp.rc
   tavg3_3d_tdt_CpDM.rc
   tavg3_3d_tdt_CpMM.rc
   tavg3_3d_trb_Cp.rc
   tavg3_3d_trb_CpDM.rc
   tavg3_3d_trb_CpMM.rc
   tavg3_3d_udt_Cp.rc
   tavg3_3d_udt_CpDM.rc
   tavg3_3d_udt_CpMM.rc
   tavg3d_cld_p.rc
   tavg3d_cld_v.rc
   tavg3d_dyn_p.rc
   tavg3d_dyn_v.rc
   tavg3d_met_e.rc
   tavg3d_mis_e.rc
   tavg3d_mst_p.rc
   tavg3d_mst_v.rc
   tavg3d_prs_v.rc
   tavg3d_tmp_p.rc
   tavg3d_tmp_v.rc
   tavg3d_trp_e.rc
   tavg3d_wnd_p.rc
   tavg3d_wnd_v.rc
   )

set(PUBLICTAG "GEOS-${CMAKE_PROJECT_VERSION}")
set(CVSTAG ${PUBLICTAG})
foreach (rcfile ${resource_files_needing_processing})
   configure_file(${rcfile} ${rcfile} @ONLY)
   install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${rcfile} DESTINATION etc)
endforeach ()
