!=======================================================================
!BOP
!
! !MODULE: ice_flux - flux variable declarations: coupler, diagnostic and
!          internal
!
! !DESCRIPTION:
!
! Flux variable declarations; these include fields sent from the coupler
! ("in"), sent to the coupler ("out"), written to diagnostic history files
! ("diagnostic"), and used internally ("internal").
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author Elizabeth C. Hunke, LANL
!
! 2004: Block structure added by William Lipscomb
!       Swappped, revised, and added some subroutines
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_flux
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &

       ! in from atmos (if .not.calc_strair)  
         strax   , & ! wind stress components (N/m^2)
         stray   , & ! 

       ! in from ocean
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty , & ! sea surface slope, y-direction

       ! out to atmosphere (if calc_strair)
         strairxT, & ! stress on ice by air, x-direction
         strairyT, & ! stress on ice by air, y-direction

       ! out to ocean          T-cell (kg/m s^2)
       ! Note, CICE_IN_NEMO uses strocnx and strocny for coupling
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

       ! diagnostic

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         sig1    , & ! principal stress component
         sig2    , & ! principal stress component
         aiceu   , & ! ice concentration on U grid 
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         strairxTsm, & ! smoothed stress on ice by air, x-direction
         strairyTsm, & ! smoothed stress on ice by air, y-direction
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strtltx , & ! stress due to sea surface slope, x-direction
         strtlty , & ! stress due to sea surface slope, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         daidtd  , & ! ice area tendency due to transport   (1/s)
         dvidtd  , & ! ice volume tendency due to transport (m/s)
         vice0   , & ! 
         vice1   , & ! 
         dardg1dt, & ! rate of area loss by ridging ice (1/s)
         dardg2dt, & ! rate of area gain by new ridges (1/s)
         dvirdgdt, & ! rate of ice volume ridged (m/s)
         divuocn,  & ! divergence of ocean surface currents (1/s)
         fakediv,  & ! divergence of a specified velocity fields (1/s)
         tauodiv,  & ! divergence of ocean-ice stress 
         tauadiv,  & ! divergence of air-ice stress
         tau1div,  & ! divergence of air-ice stress (no aice scaling)
         opening     ! rate of opening due to divergence/shear (1/s)

      real (kind=dbl_kind), allocatable, dimension (:,:,:,:) :: &
         hiceflxe, & ! ice volume transports across E cell edges 
         hiceflxn    ! ice volume transports across N cell edges 

       ! restart

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
       ! ice stress tensor in each corner of T cell (kg/s^2)
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      logical (kind=log_kind), &
         allocatable, dimension (:,:,:) :: &
         iceumask   ! ice extent mask (U-cell)

       ! internal

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         prs_sig  , & ! replacement pressure, for stress calc
         fm           ! Coriolis param. * mass in U-cell (kg/s)

      !-----------------------------------------------------------------
      ! Thermodynamic component
      !-----------------------------------------------------------------

       ! in from atmosphere (if calc_Tsfc)

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         zlvl    , & ! atm level height (m)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         wind    , & ! wind speed (m/s)
         potT    , & ! air potential temperature  (K)
         Tair    , & ! air temperature  (K)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         flw         ! incoming longwave radiation (W/m^2)

       ! in from atmosphere (if .not. Tsfc_calc)
       ! required for coupling to HadGEM3
       ! NOTE: when in CICE_IN_NEMO mode, these are gridbox mean fields,
       ! not per ice area. When in standalone mode, these are per ice area.

      real (kind=dbl_kind), & 
         allocatable, dimension (:,:,:,:) :: &
         fsurfn_f   , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f, & ! downward cond flux at top surface (W m-2)
         flatn_f        ! latent heat flux (W m-2)

       ! in from atmosphere

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow       ! snowfall rate (kg/m^2 s)

       ! in from ocean

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         frzmlt  , & ! freezing/melting potential (W/m^2)
         Tf      , & ! freezing temperature (C)
         qdp     , & ! deep ocean heat flux (W/m^2), negative upward
         hmix        ! mixed layer depth (m)

      character (char_len) :: &
         Tfrzpt      ! ocean freezing temperature formulation
                     ! 'constant' (-1.8C), 'linear_S'

       ! out to atmosphere (if calc_Tsfc)
       ! note Tsfc is in ice_state.F

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         fsens   , & ! sensible heat flux (W/m^2)
         flat    , & ! latent heat flux   (W/m^2)
         fswabs  , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         flwout  , & ! outgoing longwave radiation (W/m^2)
         Tref    , & ! 2m atm reference temperature (K)
         Qref    , & ! 2m atm reference spec humidity (kg/kg)
         evap        ! evaporative water flux (kg/m^2/s)

       ! albedos aggregated over categories (if calc_Tsfc)
      real (kind=dbl_kind), allocatable, dimension(:,:,:) :: &
         alvdr   , & ! visible, direct   (fraction)
         alidr   , & ! near-ir, direct   (fraction)
         alvdf   , & ! visible, diffuse  (fraction)
         alidf   , & ! near-ir, diffuse  (fraction)
         ! grid-box-mean versions
         alvdr_gbm, & ! visible, direct   (fraction)
         alidr_gbm, & ! near-ir, direct   (fraction)
         alvdf_gbm, & ! visible, diffuse  (fraction)
         alidf_gbm, & ! near-ir, diffuse  (fraction)
         ! components for history
         albice   , & ! bare ice albedo
         albsno   , & ! snow albedo
         albpnd   , & ! melt pond albedo
         albcnt       ! counter for zenith angle

       ! out to ocean 
       ! (Note CICE_IN_NEMO does not use these for coupling.  
       !  It uses fresh_gbm,fsalt_gbm,fhocn_gbm and fswthru_gbm)
      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         fresh   , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt   , & ! salt flux to ocean (kg/m^2/s)
         fhocn   , & ! net heat flux to ocean (W/m^2)
         fswthru     ! shortwave penetrating to ocean (W/m^2)

       ! internal

      real (kind=dbl_kind), &
         allocatable, dimension (:,:,:) :: &
         fswfac  , & ! for history
         scale_factor! scaling factor for shortwave components

      logical (kind=log_kind) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

#ifdef GEOS
      real (kind=dbl_kind)   :: &
          ice_ref_salinity
#endif

      !-----------------------------------------------------------------
      ! quantities passed from ocean mixed layer to atmosphere
      ! (for running with CAM)
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         strairx_ocn , & ! stress on ocean by air, x-direction
         strairy_ocn , & ! stress on ocean by air, y-direction
         fsens_ocn   , & ! sensible heat flux (W/m^2)
         flat_ocn    , & ! latent heat flux   (W/m^2)
         flwout_ocn  , & ! outgoing longwave radiation (W/m^2)
         evap_ocn    , & ! evaporative water flux (kg/m^2/s)
         alvdr_ocn   , & ! visible, direct   (fraction)
         alidr_ocn   , & ! near-ir, direct   (fraction)
         alvdf_ocn   , & ! visible, diffuse  (fraction)
         alidf_ocn   , & ! near-ir, diffuse  (fraction)
         Tref_ocn    , & ! 2m atm reference temperature (K)
         Qref_ocn        ! 2m atm reference spec humidity (kg/kg)

      !-----------------------------------------------------------------
      ! diagnostic
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         fsurf , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop,&! top surface conductive flux        (W/m^2)
         congel, & ! basal ice growth         (m/step-->cm/day)
         frazil, & ! frazil ice growth        (m/step-->cm/day)
         snoice, & ! snow-ice formation       (m/step-->cm/day)
         meltt , & ! top ice melt             (m/step-->cm/day)
         melts , & ! snow melt                (m/step-->cm/day)
         meltb , & ! basal ice melt           (m/step-->cm/day)
         meltl , & ! lateral ice melt         (m/step-->cm/day)
         daidtt, & ! ice area tendency thermo.   (s^-1)
         dvidtt, & ! ice volume tendency thermo. (m/s)
         mlt_onset, &! day of year that sfc melting begins
         frz_onset   ! day of year that freezing begins (congel or frazil)
         
      real (kind=dbl_kind), & 
         allocatable, dimension (:,:,:,:) :: &
         fsurfn,   & ! category fsurf
         fcondtopn,& ! category fcondtop
         flatn       ! cagegory latent heat flux

      ! As above but these remain grid box mean values i.e. they are not
      ! divided by aice at end of ice_dynamics.  These are used in
      ! CICE_IN_NEMO for coupling and also for generating
      ! ice diagnostics and history files as these are more accurate. 
      ! (The others suffer from problem of incorrect values at grid boxes
      !  that change from an ice free state to an icy state.)
    
      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         fresh_gbm, & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_gbm, & ! salt flux to ocean (kg/m^2/s)
         fhocn_gbm, & ! net heat flux to ocean (W/m^2)
         fswthru_gbm  ! shortwave penetrating to ocean (W/m^2)

      !-----------------------------------------------------------------
      ! internal
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension (:,:,:) :: &
         rside   , & ! fraction of ice that melts laterally
         fsw     , & ! incoming shortwave radiation (W/m^2)
         coszen  , & ! cosine solar zenith angle, < 0 for sun below horizon 
         rdg_conv, & ! convergence term for ridging (1/s)
         rdg_shear   ! shear term for ridging (1/s)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: alloc_flux - allocate ice_flux 
!
! !INTERFACE:
!
      subroutine alloc_flux
!
! !DESCRIPTION:
!
! Allocate ice_flux state
!
! !REVISION HISTORY:
!
! author Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         allocate(strax(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stray(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(uocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(ss_tltx(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(ss_tlty(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(strairxT(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairyT(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairxTsm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairyTsm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(strocnxT(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strocnyT(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(sig1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(sig2(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(aiceu(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairx(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairy(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strocnx(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strocny(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strtltx(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strtlty(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strintx(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strinty(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(daidtd(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(dvidtd(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vice0(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vice1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(dardg1dt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(dardg2dt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(dvirdgdt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(opening(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(divuocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fakediv(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(tauodiv(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(tauadiv(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(tau1div(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(hiceflxe(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(hiceflxn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

         allocate(stressp_1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressp_2(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressp_3(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressp_4(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressm_1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressm_2(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressm_3(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stressm_4 (nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stress12_1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stress12_2(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stress12_3(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(stress12_4(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(iceumask(nx_block,ny_block,max_blocks), source=.false.)

         allocate(prs_sig(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(zlvl(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(uatm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vatm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(wind(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(potT(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Tair(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Qa(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(rhoa(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(swvdr(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(swvdf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(swidr(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(swidf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(flw(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(fsurfn_f(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(fcondtopn_f(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(flatn_f(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

         allocate(frain(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fsnow(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(sss(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(sst(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(frzmlt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Tf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(qdp(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(hmix(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(fsens(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(flat(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fswabs(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(flwout(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Tref(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Qref(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(evap(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(alvdr(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidr(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alvdf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alvdr_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidr_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alvdf_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidf_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(albice(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(albsno(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(albpnd(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(albcnt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(fresh(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fsalt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fhocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fswthru(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(fswfac(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(scale_factor(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(strairx_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strairy_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fsens_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(flat_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(flwout_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(evap_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alvdr_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidr_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alvdf_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(alidf_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Tref_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(Qref_ocn(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(fsurf(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fcondtop(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(congel(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(frazil(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(snoice(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(meltt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(melts(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(meltb(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(meltl(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(daidtt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(dvidtt(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(mlt_onset(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(frz_onset(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         
         allocate(fsurfn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(fcondtopn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(flatn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

         allocate(fresh_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fsalt_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fhocn_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fswthru_gbm(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(rside(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(fsw(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(coszen(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(rdg_conv(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(rdg_shear(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

      end subroutine alloc_flux

!=======================================================================
!BOP
!
! !IROUTINE: dealloc_flux - deallocate ice_flux 
!
! !INTERFACE:
!
      subroutine dealloc_flux
!
! !DESCRIPTION:
!
! Deallocate ice_flux state
!
! !REVISION HISTORY:
!
! author Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         deallocate(strax)
         deallocate(stray)

         deallocate(uocn)
         deallocate(vocn)
         deallocate(ss_tltx)
         deallocate(ss_tlty)

         deallocate(strairxT)
         deallocate(strairyT)

         deallocate(strairxTsm)
         deallocate(strairyTsm)

         deallocate(strocnxT)
         deallocate(strocnyT)

         deallocate(sig1)
         deallocate(sig2)
         deallocate(aiceu)
         deallocate(strairx)
         deallocate(strairy)
         deallocate(strocnx)
         deallocate(strocny)
         deallocate(strtltx)
         deallocate(strtlty)
         deallocate(strintx)
         deallocate(strinty)
         deallocate(daidtd)
         deallocate(dvidtd)
         deallocate(vice0)
         deallocate(vice1)
         deallocate(hiceflxe)
         deallocate(hiceflxn)
         deallocate(dardg1dt)
         deallocate(dardg2dt)
         deallocate(dvirdgdt)
         deallocate(opening)
         deallocate(divuocn)
         deallocate(fakediv)
         deallocate(tauadiv)
         deallocate(tauodiv)
         deallocate(tau1div)

         deallocate(stressp_1)
         deallocate(stressp_2)
         deallocate(stressp_3)
         deallocate(stressp_4)
         deallocate(stressm_1)
         deallocate(stressm_2)
         deallocate(stressm_3)
         deallocate(stressm_4 )
         deallocate(stress12_1)
         deallocate(stress12_2)
         deallocate(stress12_3)
         deallocate(stress12_4)

         deallocate(iceumask)

         deallocate(prs_sig)
         deallocate(fm)

         deallocate(zlvl)
         deallocate(uatm)
         deallocate(vatm)
         deallocate(wind)
         deallocate(potT)
         deallocate(Tair)
         deallocate(Qa)
         deallocate(rhoa)
         deallocate(swvdr)
         deallocate(swvdf)
         deallocate(swidr)
         deallocate(swidf)
         deallocate(flw)

         deallocate(fsurfn_f)
         deallocate(fcondtopn_f)
         deallocate(flatn_f)

         deallocate(frain)
         deallocate(fsnow)

         deallocate(sss)
         deallocate(sst)
         deallocate(frzmlt)
         deallocate(Tf)
         deallocate(qdp)
         deallocate(hmix)

         deallocate(fsens)
         deallocate(flat)
         deallocate(fswabs)
         deallocate(flwout)
         deallocate(Tref)
         deallocate(Qref)
         deallocate(evap)

         deallocate(alvdr)
         deallocate(alidr)
         deallocate(alvdf)
         deallocate(alidf)
         deallocate(alvdr_gbm)
         deallocate(alidr_gbm)
         deallocate(alvdf_gbm)
         deallocate(alidf_gbm)
         deallocate(albice)
         deallocate(albsno)
         deallocate(albpnd)
         deallocate(albcnt)

         deallocate(fresh)
         deallocate(fsalt)
         deallocate(fhocn)
         deallocate(fswthru)

         deallocate(fswfac)
         deallocate(scale_factor)

         deallocate(strairx_ocn)
         deallocate(strairy_ocn)
         deallocate(fsens_ocn)
         deallocate(flat_ocn)
         deallocate(flwout_ocn)
         deallocate(evap_ocn)
         deallocate(alvdr_ocn)
         deallocate(alidr_ocn)
         deallocate(alvdf_ocn)
         deallocate(alidf_ocn)
         deallocate(Tref_ocn)
         deallocate(Qref_ocn)

         deallocate(fsurf)
         deallocate(fcondtop)
         deallocate(congel)
         deallocate(frazil)
         deallocate(snoice)
         deallocate(meltt)
         deallocate(melts)
         deallocate(meltb)
         deallocate(meltl)
         deallocate(daidtt)
         deallocate(dvidtt)
         deallocate(mlt_onset)
         deallocate(frz_onset)

         deallocate(fsurfn)
         deallocate(fcondtopn)
         deallocate(flatn)

         deallocate(fresh_gbm)
         deallocate(fsalt_gbm)
         deallocate(fhocn_gbm)
         deallocate(fswthru_gbm)
         deallocate(rside)

         deallocate(fsw)
         deallocate(coszen)
         deallocate(rdg_conv)
         deallocate(rdg_shear)

      end subroutine dealloc_flux

!=======================================================================
!BOP
!
! !IROUTINE: init_coupler_flux - initialize fluxes exchanged with coupler
!
! !INTERFACE:
!
      subroutine init_coupler_flux
!
! !DESCRIPTION:
!
! Initialize all fluxes exchanged with flux coupler
! and some data-derived fields
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, n, iblk

      logical (kind=log_kind), parameter ::     & 
!         l_winter = .true.   ! winter/summer default switch
         l_winter = .false.   ! winter/summer default switch

      real (kind=dbl_kind) :: fcondtopn_d(6), fsurfn_d(6)

      data fcondtopn_d / -50.0_dbl_kind,-17.0_dbl_kind,-12.0_dbl_kind, &
                          -9.0_dbl_kind, -7.0_dbl_kind, -3.0_dbl_kind /
      data fsurfn_d    /  0.20_dbl_kind, 0.15_dbl_kind, 0.10_dbl_kind, &
                          0.05_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind /


#ifdef GEOS
      Tfrzpt = 'linear_S'
#endif
      !-----------------------------------------------------------------
      ! fluxes received from atmosphere
      !-----------------------------------------------------------------
      zlvl  (:,:,:) = c10             ! atm level height (m)
      rhoa  (:,:,:) = 1.3_dbl_kind    ! air density (kg/m^3)
      uatm  (:,:,:) = c5              ! wind velocity    (m/s)
      vatm  (:,:,:) = c5
      strax (:,:,:) = 0.05_dbl_kind
      stray (:,:,:) = 0.05_dbl_kind
      if (l_winter) then
         !typical winter values
         potT  (:,:,:) = 253.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 253.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0006_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         flw   (:,:,:) = c180            ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         fsnow (:,:,:) = 4.0e-6_dbl_kind ! snowfall rate (kg/m2/s)
         do n = 1, ncat                  ! conductive heat flux (W/m^2)
            fcondtopn_f(:,:,n,:) = fcondtopn_d(n)
         enddo
         fsurfn_f = fcondtopn_f          ! surface heat flux (W/m^2)
         flatn_f(:,:,:,:) = c0           ! latent heat flux (kg/m2/s)
      else
         !typical summer values
         potT  (:,:,:) = 273.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 273.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0035_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         flw   (:,:,:) = 280.0_dbl_kind  ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         fsnow (:,:,:) = c0              ! snowfall rate (kg/m2/s)
         do n = 1, ncat                  ! surface heat flux (W/m^2)
            fsurfn_f(:,:,n,:) = fsurfn_d(n)
         enddo
         fcondtopn_f(:,:,:,:) = 0.0_dbl_kind ! conductive heat flux (W/m^2)
         flatn_f(:,:,:,:) = -2.0_dbl_kind    ! latent heat flux (W/m^2)
      endif !     l_winter

      !-----------------------------------------------------------------
      ! fluxes received from ocean
      !-----------------------------------------------------------------

      ss_tltx(:,:,:)= c0              ! sea surface tilt (m/m)
      ss_tlty(:,:,:)= c0
      uocn  (:,:,:) = c0              ! surface ocean currents (m/s)
      vocn  (:,:,:) = c0
      frzmlt(:,:,:) = c0              ! freezing/melting potential (W/m^2)
      sss   (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (o/oo)
      if (trim(Tfrzpt) == 'constant') then
         Tf    (:,:,:) = -1.8_dbl_kind   ! freezing temp (C)     
      else ! default:  Tfrzpt = 'linear_S'
         Tf    (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
      endif
#ifndef CICE_IN_NEMO
      sst   (:,:,:) = Tf(:,:,:)       ! sea surface temp (C)
#endif
      qdp   (:,:,:) = c0              ! deep ocean heat flux (W/m^2)
      hmix  (:,:,:) = c20             ! ocean mixed layer depth

      !-----------------------------------------------------------------
      ! fluxes sent to atmosphere
      !-----------------------------------------------------------------

!echmod - for rectangular grid tests without thermo
!      strairxT(:,:,:) = 0.15_dbl_kind
!      strairyT(:,:,:) = 0.15_dbl_kind

      strairxT(:,:,:) = c0            ! wind stress, T grid
      strairyT(:,:,:) = c0
!echmod
      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
      flwout  (:,:,:) = -stefan_boltzmann*Tffresh**4   
                        ! in case atm model diagnoses Tsfc from flwout
      evap    (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0
      alvdr   (:,:,:) = c0
      alidr   (:,:,:) = c0
      alvdf   (:,:,:) = c0
      alidf   (:,:,:) = c0

      !-----------------------------------------------------------------
      ! fluxes sent to ocean
      !-----------------------------------------------------------------

      strocnxT(:,:,:) = c0    ! ice-ocean stress, x-direction (T-cell)
      strocnyT(:,:,:) = c0    ! ice-ocean stress, y-direction (T-cell)
      fresh   (:,:,:) = c0
      fsalt   (:,:,:) = c0
      fhocn   (:,:,:) = c0
      fswthru (:,:,:) = c0

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

      fsw     (:,:,:) = c0            ! shortwave radiation (W/m^2)
      scale_factor(:,:,:) = c1        ! shortwave scaling factor 
      wind    (:,:,:) = sqrt(uatm(:,:,:)**2 &
                           + vatm(:,:,:)**2)  ! wind speed, (m/s)

      end subroutine init_coupler_flux

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_atm - initialize all atmospheric fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_atm
!
! !DESCRIPTION:
!
! Initialize some fluxes sent to coupler for use by the atm model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      !-----------------------------------------------------------------
      ! initialize albedo and fluxes
      !-----------------------------------------------------------------

      strairxT(:,:,:) = c0      ! wind stress, T grid
      strairyT(:,:,:) = c0
      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
      flwout  (:,:,:) = c0
      evap    (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0

      end subroutine init_flux_atm

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_ocn - initialize ocean fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_ocn
!
! !DESCRIPTION:
!
! Initialize some fluxes sent to coupler for use by the ocean model
!
! NOTE: These fluxes should be initialized immediately after the
!       call to the coupler.  The atmospheric fluxes can be initialized
!       at the beginning of the following time step because they are
!       not modified by any subroutines between the call to_coupler
!       and the end of the time step.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------

      fresh  (:,:,:)  = c0
      fsalt  (:,:,:)  = c0
      fhocn  (:,:,:)  = c0
      fswthru(:,:,:)  = c0

      end subroutine init_flux_ocn

!=======================================================================
!BOP
!
! !IROUTINE: init_history_therm - initialize thermo history fields
!
! !INTERFACE:
!
      subroutine init_history_therm
!
! !DESCRIPTION:
!
! Initialize thermodynamic fields written to history files.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_state, only: aice, vice

      fsurf  (:,:,:) = c0
      fcondtop(:,:,:)= c0
      congel (:,:,:) = c0
      frazil (:,:,:) = c0
      snoice (:,:,:) = c0
      meltt  (:,:,:) = c0
      melts  (:,:,:) = c0
      meltb  (:,:,:) = c0
      meltl  (:,:,:) = c0
      daidtt (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtt (:,:,:) = vice(:,:,:) ! temporary initial volume
      fsurfn    (:,:,:,:) = c0
      fcondtopn (:,:,:,:) = c0
      flatn     (:,:,:,:) = c0
      fresh_gbm  (:,:,:) = c0
      fsalt_gbm  (:,:,:) = c0
      fhocn_gbm  (:,:,:) = c0
      fswthru_gbm(:,:,:) = c0
      albice (:,:,:) = c0
      albsno (:,:,:) = c0
      albpnd (:,:,:) = c0
      
      end subroutine init_history_therm

!=======================================================================
!BOP
!
! !IROUTINE: init_history_dyn - initialize dynamic history fields
!
! !INTERFACE:
!
      subroutine init_history_dyn
!
! !DESCRIPTION:
!
! Initialize dynamic fields written to history files.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_state, only: aice, vice
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      sig1    (:,:,:) = c0
      sig2    (:,:,:) = c0
      strocnx (:,:,:) = c0
      strocny (:,:,:) = c0
      strairx (:,:,:) = c0
      strairy (:,:,:) = c0
      strtltx (:,:,:) = c0
      strtlty (:,:,:) = c0
      strintx (:,:,:) = c0
      strinty (:,:,:) = c0
      dardg1dt(:,:,:) = c0
      dardg2dt(:,:,:) = c0
      dvirdgdt(:,:,:) = c0
      opening (:,:,:) = c0
      divuocn (:,:,:) = c0
      fakediv (:,:,:) = c0
      tauodiv (:,:,:) = c0
      tauadiv (:,:,:) = c0
      tau1div (:,:,:) = c0
      daidtd  (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtd  (:,:,:) = vice(:,:,:) ! temporary initial volume
      fm      (:,:,:) = c0
      prs_sig (:,:,:) = c0

      end subroutine init_history_dyn

!=======================================================================
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_fluxes (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               aicen,                &    
                               flw,      coszn,      &
                               strairxn, strairyn,   &
                               fsurfn,   fcondtopn,  &  
                               fsensn,   flatn,      & 
                               fswabsn,  flwoutn,    &
                               evapn,                &
                               Trefn,    Qrefn,      &
                               freshn,   fsaltn,     &
                               fhocnn,   fswthrun,   &
                               strairxT, strairyT,   &  
                               fsurf,    fcondtop,   &
                               fsens,    flat,       & 
                               fswabs,   flwout,     &
                               evap,                 & 
                               Tref,     Qref,       &
                               fresh,    fsalt,      & 
                               fhocn,    fswthru,    &
                               melttn, meltsn, meltbn, congeln, snoicen, &
                               meltt,  melts,  &
                               meltb,                       &
                               congel,  snoice)
                               

!
! !DESCRIPTION:
!
! Aggregate flux information from all ice thickness categories
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &
          aicen   , & ! concentration of ice
          flw     , & ! downward longwave flux          (W/m**2)
          coszn   , & ! cosine of solar zenith angle 
          strairxn, & ! air/ice zonal  strss,           (N/m**2)
          strairyn, & ! air/ice merdnl strss,           (N/m**2)
          fsurfn  , & ! net heat flux to top surface    (W/m**2)
          fcondtopn,& ! downward cond flux at top sfc   (W/m**2)
          fsensn  , & ! sensible heat flx               (W/m**2)
          flatn   , & ! latent   heat flx               (W/m**2)
          fswabsn , & ! shortwave absorbed heat flx     (W/m**2)
          flwoutn , & ! upwd lw emitted heat flx        (W/m**2)
          evapn   , & ! evaporation                     (kg/m2/s)
          Trefn   , & ! air tmp reference level         (K)
          Qrefn   , & ! air sp hum reference level      (kg/kg)
          freshn  , & ! fresh water flux to ocean       (kg/m2/s)
          fsaltn  , & ! salt flux to ocean              (kg/m2/s)
          fhocnn  , & ! actual ocn/ice heat flx         (W/m**2)
          fswthrun, & ! sw radiation through ice bot    (W/m**2)
          melttn  , & ! top ice melt                    (m)
          meltbn  , & ! bottom ice melt                 (m)
          meltsn  , & ! snow melt                       (m)
          congeln , & ! congelation ice growth          (m)
          snoicen     ! snow-ice growth                 (m)
           
      ! cumulative fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          strairxT, & ! air/ice zonal  strss,           (N/m**2)
          strairyT, & ! air/ice merdnl strss,           (N/m**2)
          fsurf   , & ! net heat flux to top surface    (W/m**2)
          fcondtop, & ! downward cond flux at top sfc   (W/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          meltt   , & ! top ice melt                    (m)
          meltb   , & ! bottom ice melt                 (m)
          melts   , & ! snow melt                       (m)
          congel  , & ! congelation ice growth          (m)
          snoice      ! snow-ice growth                 (m)
!
!EOP
!
      integer (kind=int_kind) :: &
          ij, i, j    ! horizontal indices

      !-----------------------------------------------------------------
      ! Merge fluxes
      ! NOTE: The albedo is aggregated only in cells where ice exists
      !       and (for the delta-Eddington scheme) where the sun is above
      !       the horizon. 
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! atmo fluxes

         strairxT (i,j)  = strairxT(i,j) + strairxn(i,j)*aicen(i,j)
         strairyT (i,j)  = strairyT(i,j) + strairyn(i,j)*aicen(i,j)
         fsurf    (i,j)  = fsurf   (i,j) + fsurfn  (i,j)*aicen(i,j)
         fcondtop (i,j)  = fcondtop(i,j) + fcondtopn(i,j)*aicen(i,j) 
         fsens    (i,j)  = fsens   (i,j) + fsensn  (i,j)*aicen(i,j)
         flat     (i,j)  = flat    (i,j) + flatn   (i,j)*aicen(i,j)
         fswabs   (i,j)  = fswabs  (i,j) + fswabsn (i,j)*aicen(i,j)
         flwout   (i,j)  = flwout  (i,j) &
             + (flwoutn(i,j) - (c1-emissivity)*flw(i,j))*aicen(i,j)
         evap     (i,j)  = evap    (i,j) + evapn   (i,j)*aicen(i,j)
         Tref     (i,j)  = Tref    (i,j) + Trefn   (i,j)*aicen(i,j)
         Qref     (i,j)  = Qref    (i,j) + Qrefn   (i,j)*aicen(i,j)

         ! ocean fluxes

         fresh    (i,j) = fresh    (i,j) + freshn  (i,j)*aicen(i,j)
         fsalt    (i,j) = fsalt    (i,j) + fsaltn  (i,j)*aicen(i,j)
         fhocn    (i,j) = fhocn    (i,j) + fhocnn  (i,j)*aicen(i,j)
         fswthru  (i,j) = fswthru  (i,j) + fswthrun(i,j)*aicen(i,j)

         ! ice/snow thickness

         meltt    (i,j) = meltt    (i,j) + melttn  (i,j)*aicen(i,j)
         meltb    (i,j) = meltb    (i,j) + meltbn  (i,j)*aicen(i,j)
         melts    (i,j) = melts    (i,j) + meltsn  (i,j)*aicen(i,j)
         congel   (i,j) = congel   (i,j) + congeln (i,j)*aicen(i,j)
         snoice   (i,j) = snoice   (i,j) + snoicen (i,j)*aicen(i,j)

      enddo                     ! ij
      
      end subroutine merge_fluxes

!=======================================================================
!BOP
!
! !IROUTINE: scale_fluxes
!
! !DESCRIPTION:
!
!  Divide ice fluxes by ice area before sending them to the
!  coupler, since the coupler multiplies by ice area.
!
! !INTERFACE:
!
      subroutine scale_fluxes (nx_block, ny_block, &
                               tmask,              &
                               aice,     Tf,       &
                               Tair,     Qa,       &
                               strairxT, strairyT, &
                               fsens,    flat,     &
                               fswabs,   flwout,   &
                               evap,               &
                               Tref,     Qref,     &
                               fresh,    fsalt,    &
                               fhocn,    fswthru,  &
                               alvdr,    alidr,    &
                               alvdf,    alidf)
!
! !REVISION HISTORY:
!
! authors: C.M.Bitz, William H. Lipscomb
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block    ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask     ! land/boundary mask, thickness (T-cell)


      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice    , & ! fractional ice area
          Tf      , & ! freezing temperature            (C)
          Tair    , & ! surface air temperature         (K)
          Qa          ! sfc air specific humidity       (kg/kg)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          strairxT, & ! air/ice zonal  stress           (N/m**2)
          strairyT, & ! air/ice merdnl stress           (N/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          alvdr   , & ! visible, direct   (fraction)
          alidr   , & ! near-ir, direct   (fraction)
          alvdf   , & ! visible, diffuse  (fraction)
          alidf       ! near-ir, diffuse  (fraction)
!
!EOP
!
      real (kind=dbl_kind) :: ar   ! 1/aice

      integer (kind=int_kind) :: &
          i, j    ! horizontal indices


!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j) .and. aice(i,j) > c0) then
            ar = c1 / aice(i,j)
            strairxT(i,j) = strairxT(i,j) * ar
            strairyT(i,j) = strairyT(i,j) * ar
            fsens   (i,j) = fsens   (i,j) * ar
            flat    (i,j) = flat    (i,j) * ar
            fswabs  (i,j) = fswabs  (i,j) * ar
            flwout  (i,j) = flwout  (i,j) * ar
            evap    (i,j) = evap    (i,j) * ar
            Tref    (i,j) = Tref    (i,j) * ar
            Qref    (i,j) = Qref    (i,j) * ar
            fresh   (i,j) = fresh   (i,j) * ar
            fsalt   (i,j) = fsalt   (i,j) * ar
            fhocn   (i,j) = fhocn   (i,j) * ar
            fswthru (i,j) = fswthru (i,j) * ar
            alvdr   (i,j) = alvdr   (i,j) * ar
            alidr   (i,j) = alidr   (i,j) * ar
            alvdf   (i,j) = alvdf   (i,j) * ar
            alidf   (i,j) = alidf   (i,j) * ar
         else                   ! zero out fluxes
            strairxT(i,j) = c0
            strairyT(i,j) = c0
            fsens   (i,j) = c0
            flat    (i,j) = c0
            fswabs  (i,j) = c0
            flwout  (i,j) = -stefan_boltzmann *(Tf(i,j) + Tffresh)**4
               ! to make upward longwave over ocean reasonable for history file
            evap    (i,j) = c0
            Tref    (i,j) = Tair(i,j)
            Qref    (i,j) = Qa  (i,j)
            fresh   (i,j) = c0
            fsalt   (i,j) = c0
            fhocn   (i,j) = c0
            fswthru (i,j) = c0
            alvdr   (i,j) = c0  ! zero out albedo where ice is absent
            alidr   (i,j) = c0
            alvdf   (i,j) = c0 
            alidf   (i,j) = c0
         endif                  ! tmask and aice > 0
      enddo                     ! i
      enddo                     ! j
      
      end subroutine scale_fluxes

#ifdef GEOS

      subroutine smooth_windstress(strxin, stryin, strxout, stryout) 

      use ice_boundary
      use ice_domain
      use ice_state,   only: aice
      use ice_grid,    only: t2ugrid_vector, ANGLET, dxu, dyu, tmask, &
                             umask, hm
      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
         strxin, stryin
      real (kind=dbl_kind), &
         dimension(:,:,:), intent(out) :: &
              strxout, stryout

         
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
            workx, worky

      integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: &
         niters = 10

      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              if(tmask(i,j,iblk)) then
                strxout(i,j,iblk) = strxin(i1,j1)
                stryout(i,j,iblk) = stryin(i1,j1)
              else
                strxout(i,j,iblk) = c0 
                stryout(i,j,iblk) = c0 
              endif  
         enddo !i
         enddo !j

      enddo                     ! iblk
      !call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strxout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (stryout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      !call ice_timer_stop(timer_bound)

      do n=1, niters  

      workx(:,:,:) = strxout(:,:,:) 
      worky(:,:,:) = stryout(:,:,:) 

      strxout(:,:,:) = c0
      stryout(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
             if(tmask(i,j,iblk) .and. aice(i,j,iblk) > puny) then 
             k = 0 
             do j1=j-1,j+1
             do i1=i-1,i+1
                if(tmask(i1,j1,iblk) .and. aice(i,j,iblk) > puny) then
                  strxout(i,j,iblk) = strxout(i,j,iblk) + &
                                      workx(i1,j1,iblk)
                  stryout(i,j,iblk) = stryout(i,j,iblk) + &
                                      worky(i1,j1,iblk)
                  k = k+1 
                endif
             enddo !i1
             enddo !j1
             strxout(i,j,iblk) = strxout(i,j,iblk) / &
                                  real(k, kind=dbl_kind)
             stryout(i,j,iblk) = stryout(i,j,iblk) / & 
                                  real(k, kind=dbl_kind)
             endif
         enddo !i
         enddo !j
      enddo                     ! iblk

      !call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strxout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (stryout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      !call ice_timer_stop(timer_bound)

      enddo                     ! n
          

      end subroutine smooth_windstress 

      subroutine transformA2B(strxin, stryin, strxout, stryout, undef) 

      use ice_boundary
      use ice_domain
      use ice_state,   only: aice
      use ice_grid,    only: t2ugrid_vector, ANGLET, dxu, dyu, tmask, &
                             umask, hm
      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
         strxin, stryin
      real (kind=dbl_kind), &
         dimension(:,:,:), intent(out) :: &
              strxout, stryout
      real (kind=dbl_kind), intent(in) :: &
             undef

         
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
            workx, worky

      integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n

      type (block) :: &
         this_block           ! block information for current block


      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              if(tmask(i,j,iblk)) then
                strxout(i,j,iblk) = strxin(i1,j1)
                stryout(i,j,iblk) = stryin(i1,j1)
              else
                strxout(i,j,iblk) = undef 
                stryout(i,j,iblk) = undef
              endif  
         enddo !i
         enddo !j

      enddo                     ! iblk
      !call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strxout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (stryout,           halo_info, &
                           field_loc_center,  field_type_scalar)
      !call ice_timer_stop(timer_bound)


      workx(:,:,:) = strxout(:,:,:) 
      worky(:,:,:) = stryout(:,:,:) 

      strxout(:,:,:) = c0
      stryout(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
             if(umask(i,j,iblk)) then 
             k = 0 
             do j1=j,j+1
             do i1=i,i+1
                if(workx(i,j,iblk) /= undef) then
                  strxout(i,j,iblk) = strxout(i,j,iblk) + &
                                      workx(i1,j1,iblk)
                  stryout(i,j,iblk) = stryout(i,j,iblk) + &
                                      worky(i1,j1,iblk)
                  k = k+1 
                endif
             enddo !i1
             enddo !j1
             strxout(i,j,iblk) = strxout(i,j,iblk) / &
                                  real(k, kind=dbl_kind)
             stryout(i,j,iblk) = stryout(i,j,iblk) / & 
                                  real(k, kind=dbl_kind)
             endif !umask
         enddo !i
         enddo !j
      enddo                     ! iblk
          
      end subroutine transformA2B 

      subroutine set_atm_fluxes(strx, stry, uw, vw, ssh, undef)

      use ice_boundary
      use ice_domain
      use ice_state,   only: aice
      use ice_work,    only: work1 
      use ice_grid,    only: t2ugrid_vector, ANGLET, dxu, dyu, tmask, &
                             umask, hm

      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
         strx, stry
      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
              uw, vw, ssh
      !logical (kind=log_kind), intent(in) :: smooth 
      real (kind=real_kind), intent(in) :: &
            undef
 

      integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
      type (block) :: &
         this_block           ! block information for current block
      real (kind=real_kind) :: &
          workx, worky            

    
        !if(smooth) then  
        !   call smooth_windstress(strx, stry, tauxsm, tauysm) 
        !endif 
        !call transformA2B(strx, stry, tauxsm, tauysm, &
        !                  real(undef,kind=dbl_kind)) 

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              if(tmask(i,j,iblk)) then
                work1(i,j,iblk) = ssh(i1,j1)  
                workx        = strx(i1,j1)
                worky        = stry(i1,j1)
                strax(i,j,iblk) = workx
                stray(i,j,iblk) = worky
                strairxT(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                                   + worky*sin(ANGLET(i,j,iblk))   ! note strax, stray, wind
                strairyT(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & !  are on the T-grid here
                                   - workx*sin(ANGLET(i,j,iblk))
                !**** comment out for now 
                !**** will need to uncomment out later 
                !strairxT(i,j,iblk) = strairxT(i,j,iblk) * aice(i,j,iblk)
                !strairyT(i,j,iblk) = strairyT(i,j,iblk) * aice(i,j,iblk)
              else
                strax(i,j,iblk)    =  c0 
                stray(i,j,iblk)    =  c0 
                strairxT(i,j,iblk) =  c0 
                strairyT(i,j,iblk) =  c0 
                work1(i,j,iblk)    =  c0
              endif  
              if(umask(i,j,iblk)) then
                !*** uw and vw from MOM are along native grid directions 
                !*** and on B-grid(velocity point) already 
                uocn(i,j,iblk)  =  uw(i1,j1)
                vocn(i,j,iblk)  =  vw(i1,j1)
#if 0
                workx           =  tauxsm(i,j,iblk)
                worky           =  tauysm(i,j,iblk)
                strax(i,j,iblk) =  workx
                stray(i,j,iblk) =  worky
                strairxT(i,j,iblk) = workx*cos(ANGLE(i,j,iblk)) & ! convert to POP grid
                                   + worky*sin(ANGLE(i,j,iblk))   ! note strax, stray, wind
                strairyT(i,j,iblk) = worky*cos(ANGLE(i,j,iblk)) & !  are on the T-grid here
                                   - workx*sin(ANGLE(i,j,iblk))
                !**** will multiply by aiu later 
                !strairxT(i,j,iblk) = strairxT(i,j,iblk) * aice(i,j,iblk)
                !strairyT(i,j,iblk) = strairyT(i,j,iblk) * aice(i,j,iblk)
#endif
              else
                uocn(i,j,iblk)  =  c0 
                vocn(i,j,iblk)  =  c0 
#if 0
                strax(i,j,iblk) =  c0 
                stray(i,j,iblk) =  c0 
                strairxT(i,j,iblk) =  c0 
                strairyT(i,j,iblk) =  c0 
#endif
              endif 
           enddo !i
           enddo !j
         enddo ! nblocks 

         call ice_HaloUpdate (work1,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (uocn,             halo_info, &
                              field_loc_NEcorner, field_type_vector)
         call ice_HaloUpdate (vocn,             halo_info, &
                              field_loc_NEcorner, field_type_vector)

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              if(umask(i,j,iblk)) then
                ss_tltx(i,j,iblk) = p5*(work1(i+1,j+1,iblk)-work1(i,j+1,iblk)  &
                                       +work1(i+1,j  ,iblk)-work1(i,j  ,iblk)) &
                                       /dxu(i,j,iblk) 
                ss_tlty(i,j,iblk) = p5*(work1(i+1,j+1,iblk)+work1(i,j+1,iblk)  &
                                       -work1(i+1,j  ,iblk)-work1(i,j  ,iblk)) &
                                       /dyu(i,j,iblk) 
              else
                ss_tltx(i,j,iblk) = c0 
                ss_tlty(i,j,iblk) = c0 
              endif  
           enddo !i
           enddo !j
         enddo ! nblocks 

         strairxTsm(:,:,:) = strax(:,:,:)
         strairyTsm(:,:,:) = stray(:,:,:)
 
      end subroutine set_atm_fluxes

      subroutine init_stress_tensor(iceu, strcomp)

      use ice_boundary
      use ice_domain
      use ice_blocks
      use ice_work,    only: work1

      !logical (kind=log_kind), dimension (:,:), & 
      real (kind=real_kind), dimension (:,:), & 
         intent(in) :: &
         iceu    ! ice extent mask (U-cell)

      real (kind=dbl_kind), &
         dimension(:,:,:), intent(in) :: &
         strcomp

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              work1(i,j,iblk)    = real(iceu(i1,j1), kind=dbl_kind)
           enddo 
           enddo 
         enddo 
         call ice_HaloUpdate (work1,           halo_info, &
                           field_loc_center,   field_type_scalar)

         iceumask(:,:,:) = .false.
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
            enddo
            enddo
         enddo

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              stressp_1(i,j,iblk)   = strcomp(i1,j1, 1)
              stressp_2(i,j,iblk)   = strcomp(i1,j1, 2)
              stressp_3(i,j,iblk)   = strcomp(i1,j1, 3)
              stressp_4(i,j,iblk)   = strcomp(i1,j1, 4)
              stressm_1(i,j,iblk)   = strcomp(i1,j1, 5)
              stressm_2(i,j,iblk)   = strcomp(i1,j1, 6)
              stressm_3(i,j,iblk)   = strcomp(i1,j1, 7)
              stressm_4(i,j,iblk)   = strcomp(i1,j1, 8)
              stress12_1(i,j,iblk)  = strcomp(i1,j1, 9)
              stress12_2(i,j,iblk)  = strcomp(i1,j1,10)
              stress12_3(i,j,iblk)  = strcomp(i1,j1,11)
              stress12_4(i,j,iblk)  = strcomp(i1,j1,12)
           enddo 
           enddo 
         enddo 
 
      end subroutine init_stress_tensor

      subroutine   gather_scatter_stress

      use ice_communicate,    only: master_task
      use ice_domain
      use ice_work,           only: work_g1, work_g2
      use ice_gather_scatter, only: gather_global, scatter_global_stress

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n, im,jm
       integer (kind=int_kind) :: &
         ig, jg
       type (block) :: &
         this_block           ! block information for current block

!#ifdef GEOS
#if 0
      do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
      do j = jlo, jhi
      do i = ilo, ihi
         ig = this_block%i_glob(i)
         jg = this_block%j_glob(j)
         !if ( ig == 180 .and. jg == 410) then
         if ( ig == 540 .and. jg == 410) then
            print*, 'before global scatter'
            print*, 'Global i, j,', ig, jg
            print*, 'stressp_1,2',  stressp_1(i,j,iblk), stressp_2(i,j,iblk)
            print*, 'stressp_3,4',  stressp_3(i,j,iblk), stressp_4(i,j,iblk)
            print*, 'stressm_1,2',  stressm_1(i,j,iblk), stressm_2(i,j,iblk)
            print*, 'stressm_3,4',  stressm_3(i,j,iblk), stressm_4(i,j,iblk)
            print*, 'stress12_1,2', stress12_1(i,j,iblk), stress12_2(i,j,iblk)
            print*, 'stress12_3,4', stress12_3(i,j,iblk), stress12_4(i,j,iblk)
            print*, 'NE neighbor i, j,',  this_block%i_glob(i+1), &
                                          this_block%j_glob(j+1)
            print*, 'stressp_1,2',  stressp_1(i+1,j+1,iblk), stressp_2(i+1,j+1,iblk)
            print*, 'stressp_3,4',  stressp_3(i+1,j+1,iblk), stressp_4(i+1,j+1,iblk)
            print*, 'stressm_1,2',  stressm_1(i+1,j+1,iblk), stressm_2(i+1,j+1,iblk)
            print*, 'stressm_3,4',  stressm_3(i+1,j+1,iblk), stressm_4(i+1,j+1,iblk)
            print*, 'stress12_1,2', stress12_1(i+1,j+1,iblk), stress12_2(i+1,j+1,iblk)
            print*, 'stress12_3,4', stress12_3(i+1,j+1,iblk), stress12_4(i+1,j+1,iblk)
            print*, 'N neighbor i, j,',  this_block%i_glob(i), &
                                         this_block%j_glob(j+1)
            print*, 'stressp_1,2',  stressp_1(i,j+1,iblk), stressp_2(i,j+1,iblk)
            print*, 'stressp_3,4',  stressp_3(i,j+1,iblk), stressp_4(i,j+1,iblk)
            print*, 'stressm_1,2',  stressm_1(i,j+1,iblk), stressm_2(i,j+1,iblk)
            print*, 'stressm_3,4',  stressm_3(i,j+1,iblk), stressm_4(i,j+1,iblk)
            print*, 'stress12_1,2', stress12_1(i,j+1,iblk), stress12_2(i,j+1,iblk)
            print*, 'stress12_3,4', stress12_3(i,j+1,iblk), stress12_4(i,j+1,iblk)
            print*, 'E neighbor i, j,',  this_block%i_glob(i+1), &
                                         this_block%j_glob(j)
            print*, 'stressp_1,2',  stressp_1(i+1,j,iblk), stressp_2(i+1,j,iblk)
            print*, 'stressp_3,4',  stressp_3(i+1,j,iblk), stressp_4(i+1,j,iblk)
            print*, 'stressm_1,2',  stressm_1(i+1,j,iblk), stressm_2(i+1,j,iblk)
            print*, 'stressm_3,4',  stressm_3(i+1,j,iblk), stressm_4(i+1,j,iblk)
            print*, 'stress12_1,2', stress12_1(i+1,j,iblk), stress12_2(i+1,j,iblk)
            print*, 'stress12_3,4', stress12_3(i+1,j,iblk), stress12_4(i+1,j,iblk)
         endif
      enddo
      enddo
      enddo
#endif
      allocate (work_g1(nx_global,ny_global), &
                work_g2(nx_global,ny_global))

      call gather_global(work_g1, stressp_1, master_task, distrb_info)
      call gather_global(work_g2, stressp_3, master_task, distrb_info)

      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call gather_global(work_g1, stressp_2, master_task, distrb_info)
      call gather_global(work_g2, stressp_4, master_task, distrb_info)

      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call gather_global(work_g1, stressm_1, master_task, distrb_info)
      call gather_global(work_g2, stressm_3, master_task, distrb_info)

      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call gather_global(work_g1, stressm_2, master_task, distrb_info)
      call gather_global(work_g2, stressm_4, master_task, distrb_info)

      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call gather_global(work_g1, stress12_1, master_task, distrb_info)
      call gather_global(work_g2, stress12_3, master_task, distrb_info)

      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call gather_global(work_g1, stress12_2, master_task, distrb_info)
      call gather_global(work_g2, stress12_4, master_task, distrb_info)

      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
           !-----------------------------------------------------------------
           ! Ensure unused stress values in west and south ghost cells are 0
           !-----------------------------------------------------------------

           do j = 1, nghost
           do i = 1, nx_block
             stressp_1 (i,j,iblk) = c0
             stressp_2 (i,j,iblk) = c0
             stressp_3 (i,j,iblk) = c0
             stressp_4 (i,j,iblk) = c0
             stressm_1 (i,j,iblk) = c0
             stressm_2 (i,j,iblk) = c0
             stressm_3 (i,j,iblk) = c0
             stressm_4 (i,j,iblk) = c0
             stress12_1(i,j,iblk) = c0
             stress12_2(i,j,iblk) = c0
             stress12_3(i,j,iblk) = c0
             stress12_4(i,j,iblk) = c0
           enddo
           enddo
           do j = 1, ny_block
           do i = 1, nghost
             stressp_1 (i,j,iblk) = c0
             stressp_2 (i,j,iblk) = c0
             stressp_3 (i,j,iblk) = c0
             stressp_4 (i,j,iblk) = c0
             stressm_1 (i,j,iblk) = c0
             stressm_2 (i,j,iblk) = c0
             stressm_3 (i,j,iblk) = c0
             stressm_4 (i,j,iblk) = c0
             stress12_1(i,j,iblk) = c0
             stress12_2(i,j,iblk) = c0
             stress12_3(i,j,iblk) = c0
             stress12_4(i,j,iblk) = c0
           enddo
           enddo

      enddo

!#ifdef GEOS
#if 0
      do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
      do j = jlo, jhi
      do i = ilo, ihi
         ig = this_block%i_glob(i)
         jg = this_block%j_glob(j)
         !if ( ig == 180 .and. jg == 410) then
         if ( ig == 540 .and. jg == 410) then
            print*, 'after global scatter'
            print*, 'Global i, j,', ig, jg
            print*, 'stressp_1,2',  stressp_1(i,j,iblk), stressp_2(i,j,iblk)
            print*, 'stressp_3,4',  stressp_3(i,j,iblk), stressp_4(i,j,iblk)
            print*, 'stressm_1,2',  stressm_1(i,j,iblk), stressm_2(i,j,iblk)
            print*, 'stressm_3,4',  stressm_3(i,j,iblk), stressm_4(i,j,iblk)
            print*, 'stress12_1,2', stress12_1(i,j,iblk), stress12_2(i,j,iblk)
            print*, 'stress12_3,4', stress12_3(i,j,iblk), stress12_4(i,j,iblk)
            print*, 'NE neighbor i, j,',  this_block%i_glob(i+1), &
                                          this_block%j_glob(j+1)
            print*, 'stressp_1,2',  stressp_1(i+1,j+1,iblk), stressp_2(i+1,j+1,iblk)
            print*, 'stressp_3,4',  stressp_3(i+1,j+1,iblk), stressp_4(i+1,j+1,iblk)
            print*, 'stressm_1,2',  stressm_1(i+1,j+1,iblk), stressm_2(i+1,j+1,iblk)
            print*, 'stressm_3,4',  stressm_3(i+1,j+1,iblk), stressm_4(i+1,j+1,iblk)
            print*, 'stress12_1,2', stress12_1(i+1,j+1,iblk), stress12_2(i+1,j+1,iblk)
            print*, 'stress12_3,4', stress12_3(i+1,j+1,iblk), stress12_4(i+1,j+1,iblk)
            print*, 'N neighbor i, j,',  this_block%i_glob(i), &
                                         this_block%j_glob(j+1)
            print*, 'stressp_1,2',  stressp_1(i,j+1,iblk), stressp_2(i,j+1,iblk)
            print*, 'stressp_3,4',  stressp_3(i,j+1,iblk), stressp_4(i,j+1,iblk)
            print*, 'stressm_1,2',  stressm_1(i,j+1,iblk), stressm_2(i,j+1,iblk)
            print*, 'stressm_3,4',  stressm_3(i,j+1,iblk), stressm_4(i,j+1,iblk)
            print*, 'stress12_1,2', stress12_1(i,j+1,iblk), stress12_2(i,j+1,iblk)
            print*, 'stress12_3,4', stress12_3(i,j+1,iblk), stress12_4(i,j+1,iblk)
            print*, 'E neighbor i, j,',  this_block%i_glob(i+1), &
                                         this_block%j_glob(j)
            print*, 'stressp_1,2',  stressp_1(i+1,j,iblk), stressp_2(i+1,j,iblk)
            print*, 'stressp_3,4',  stressp_3(i+1,j,iblk), stressp_4(i+1,j,iblk)
            print*, 'stressm_1,2',  stressm_1(i+1,j,iblk), stressm_2(i+1,j,iblk)
            print*, 'stressm_3,4',  stressm_3(i+1,j,iblk), stressm_4(i+1,j,iblk)
            print*, 'stress12_1,2', stress12_1(i+1,j,iblk), stress12_2(i+1,j,iblk)
            print*, 'stress12_3,4', stress12_3(i+1,j,iblk), stress12_4(i+1,j,iblk)
         endif
      enddo
      enddo
      enddo
#endif

      end subroutine  gather_scatter_stress


      subroutine get_bottom_stresses(tauxb, tauyb)

      use ice_domain
      use ice_grid,        only:  tmask, ANGLET

      real (kind=real_kind), &
      !real (kind=dbl_kind), &
         dimension(:,:), intent(out) :: &
         tauxb, tauyb

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block
       real (kind=dbl_kind) :: &
          workx, worky 


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
             i1 = i-nghost
             j1 = j-nghost  
             if(tmask(i,j,iblk)) then
               workx = strocnxT(i,j,iblk)
               worky = strocnyT(i,j,iblk)
               !*** need to rotate them back onto lon-lat grid
               !*** MOM will rotate them again onto its native (tripolar) grid   
               tauxb(i1,j1)  =  workx*cos(ANGLET(i,j,iblk)) &
                              - worky*sin(ANGLET(i,j,iblk))  
               tauyb(i1,j1)  =  worky*cos(ANGLET(i,j,iblk)) &
                              + workx*sin(ANGLET(i,j,iblk))   
             else
               tauxb(i1,j1) = 0.0 
               tauyb(i1,j1) = 0.0
             endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_bottom_stresses

      subroutine get_smooth_windstress(taux, tauy)

      use ice_domain
      use ice_grid,        only:  tmask, ANGLET

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         taux, tauy

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
             i1 = i-nghost
             j1 = j-nghost  
             if(tmask(i,j,iblk)) then
               taux(i1,j1) = real(strairxTsm(i,j,iblk), kind=real_kind)  
               tauy(i1,j1) = real(strairyTsm(i,j,iblk), kind=real_kind)  
             endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_smooth_windstress

      subroutine get_ocn_fluxes(fre, sal, hea)

      use ice_domain

      real (kind=real_kind), &
         dimension(:,:), optional, intent(inout) :: &
         fre, sal, hea

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
             i1 = i-nghost
             j1 = j-nghost  
             if(present(fre)) fre(i1,j1) = fresh(i,j,iblk) 
             if(present(sal)) sal(i1,j1) = fsalt(i,j,iblk) 
             if(present(hea)) hea(i1,j1) = fhocn(i,j,iblk) 
           enddo 
           enddo 
         enddo 
 
      end subroutine get_ocn_fluxes

      subroutine get_strair_term(strx, stry)

      use ice_domain
      use ice_grid,        only:  umask

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         strx, stry


       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block
       real (kind=dbl_kind) :: &
          workx, worky 


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              if(umask(i,j,iblk)) then
                strx(i1,j1) = strairx(i,j,iblk)
                stry(i1,j1) = strairy(i,j,iblk)
              endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_strair_term

      subroutine get_strtilt_term(strx, stry, undef)

      use ice_domain
      use ice_grid,        only:  umask

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         strx, stry
      real (kind=real_kind), intent(in) :: &
         undef


       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block
       real (kind=dbl_kind) :: &
          workx, worky 


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              if(umask(i,j,iblk)) then
                strx(i1,j1) = strtltx(i,j,iblk)
                stry(i1,j1) = strtlty(i,j,iblk)
              else
                strx(i1,j1) = undef
                stry(i1,j1) = undef
              endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_strtilt_term

      subroutine get_strint_term(strx, stry, undef)

      use ice_domain
      use ice_grid,        only:  umask

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         strx, stry
      real (kind=real_kind), intent(in) :: &
         undef


       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block
       real (kind=dbl_kind) :: &
          workx, worky 


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              if(umask(i,j,iblk)) then
                strx(i1,j1) = strintx(i,j,iblk)
                stry(i1,j1) = strinty(i,j,iblk)
              else
                strx(i1,j1) = undef
                stry(i1,j1) = undef
              endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_strint_term

      subroutine get_ice_flux_vars(v0, v1)
                                 

      use ice_domain
      use ice_grid,        only:  tmask, umask, uoceanfr

      real (kind=real_kind), &
         dimension(:,:), optional, intent(out) :: &
          v0, v1

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              if(tmask(i,j,iblk)) then
                if(present(v0)) then
                    v0(i1,j1) = vice0(i,j,iblk)
                endif
                if(present(v1)) then
                    v1(i1,j1) = vice1(i,j,iblk)
                endif
              endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_ice_flux_vars

      subroutine get_momentum_terms(divuo, fakedivo, tauodivo, &
                                  tauadivo, tau1divo, &
                                   aiuo,              & 
                                   strcorxo, strcoryo, &
                                   strocnxo, strocnyo)

      use ice_domain
      use ice_state,       only:  uvel, vvel  
      use ice_grid,        only:  tmask, umask, uoceanfr

      real (kind=real_kind), &
         dimension(:,:), optional, intent(out) :: &
          divuo, fakedivo, tauodivo, aiuo, tauadivo, tau1divo
      real (kind=real_kind), &
         dimension(:,:), optional, intent(out) :: &
           strcorxo, strcoryo, strocnxo, strocnyo


       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block
       real (kind=dbl_kind) :: &
          workx, worky 

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              if(tmask(i,j,iblk)) then
                if(present(divuo)) divuo(i1,j1) = divuocn(i,j,iblk)
                if(present(fakedivo)) then
                    fakedivo(i1,j1) = fakediv(i,j,iblk)
                endif
                if(present(tauadivo)) then
                    tauadivo(i1,j1) = tauadiv(i,j,iblk)
                endif
                if(present(tau1divo)) then
                    tau1divo(i1,j1) = tau1div(i,j,iblk)
                endif
                if(present(tauodivo)) then
                    tauodivo(i1,j1) = tauodiv(i,j,iblk)
                endif
              endif
              if(umask(i,j,iblk)) then
                 if(present(strcorxo)) then
                   strcorxo(i1,j1) =  fm(i,j,iblk) * vvel(i,j,iblk)
                 endif  
                 if(present(strcoryo)) then
                   strcoryo(i1,j1) = -fm(i,j,iblk) * uvel(i,j,iblk)
                 endif  
                 if(present(strocnxo)) then
                   strocnxo(i1,j1) = uoceanfr(i,j,iblk) * strocnx(i,j,iblk) 
                 endif  
                 if(present(strocnyo)) then
                   strocnyo(i1,j1) = uoceanfr(i,j,iblk) * strocny(i,j,iblk) 
                 endif  
                 if(present(aiuo)) then
                   aiuo(i1,j1) = aiceu(i,j,iblk) 
                 endif  
              endif
           enddo 
           enddo 
         enddo 
 
      end subroutine get_momentum_terms

      subroutine get_ice_tendencies(at, vt)

      use ice_domain

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         at, vt

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              at(i1,j1)  = daidtd(i,j,iblk) 
              vt(i1,j1)  = dvidtd(i,j,iblk)  
           enddo 
           enddo 
         enddo 
 
      end subroutine get_ice_tendencies

      subroutine get_ice_ridge_tendencies(vt, rindex)

      use ice_domain

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         vt

       integer (kind=int_kind), intent(in) :: &
           rindex

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              select case (rindex)
              case (1) 
                 vt(i1,j1)  = dardg1dt(i,j,iblk)  
              case (2) 
                 vt(i1,j1)  = dardg2dt(i,j,iblk)  
              case (3) 
                 vt(i1,j1)  = dvirdgdt(i,j,iblk)  
              case (4) 
                 vt(i1,j1)  = opening(i,j,iblk)  
              end select
           enddo 
           enddo 
         enddo 
 
      end subroutine get_ice_ridge_tendencies

      subroutine get_ice_transport_flux(ve, vn)

      use ice_domain

      real (kind=real_kind), &
         dimension(:,:,:), intent(out) :: &
         ve, vn

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block


        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost  
              ve(i1,j1,:)  = hiceflxe(i,j,:,iblk)  
              vn(i1,j1,:)  = hiceflxn(i,j,:,iblk)  
           enddo 
           enddo 
         enddo 
 
      end subroutine get_ice_transport_flux


      subroutine finalize_stress_tensor(iceu, strcomp)

      use ice_communicate,    only: master_task
      use ice_domain

      !logical (kind=log_kind), dimension (:,:), & 
      real(kind=real_kind), dimension (:,:), & 
         intent(out) :: &
         iceu    ! ice extent mask (U-cell)

      real (kind=dbl_kind), &
         dimension(:,:,:), intent(out) :: &
         strcomp

       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block



        !work1(:,:,:) = c0

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i-nghost
              j1 = j-nghost
              !work1(i,j,iblk) = c0
              !if (iceumask(i,j,iblk)) work1(i,j,iblk) = c1              
              !iceu(i1,j1)    = real(work1(i,j,iblk), kind=real_kind)
              iceu(i1,j1) = c0   
              if (iceumask(i,j,iblk)) &
                  iceu(i1,j1)    = real(c1,  kind=real_kind)        
              strcomp(i1,j1, 1)  = stressp_1(i,j,iblk)   
              strcomp(i1,j1, 2)  = stressp_2(i,j,iblk)   
              strcomp(i1,j1, 3)  = stressp_3(i,j,iblk)    
              strcomp(i1,j1, 4)  = stressp_4(i,j,iblk)   
              strcomp(i1,j1, 5)  = stressm_1(i,j,iblk)   
              strcomp(i1,j1, 6)  = stressm_2(i,j,iblk)    
              strcomp(i1,j1, 7)  = stressm_3(i,j,iblk)    
              strcomp(i1,j1, 8)  = stressm_4(i,j,iblk)    
              strcomp(i1,j1, 9)  = stress12_1(i,j,iblk)   
              strcomp(i1,j1,10)  = stress12_2(i,j,iblk)   
              strcomp(i1,j1,11)  = stress12_3(i,j,iblk)   
              strcomp(i1,j1,12)  = stress12_4(i,j,iblk)   
           enddo 
           enddo 
         enddo 

 
      end subroutine finalize_stress_tensor

#endif

!=======================================================================

      end module ice_flux

!=======================================================================
