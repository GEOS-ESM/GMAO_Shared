!=========================================================================
!BOP
!
! !MODULE: ice_therm_vertical - thermo calculations before call to coupler
!
! !DESCRIPTION:
!
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
!
! !INTERFACE:
!
      module ice_therm_vertical
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size, only: ncat, nilyr, nslyr, ntilyr, ntslyr, ntrcr
      use ice_constants
      use ice_fileunits, only: nu_diag
      use ice_age, only: tr_iage
!
!EOP
!
      implicit none
      save

      real (kind=dbl_kind), parameter :: &
         saltmax = 3.2_dbl_kind,  & ! max salinity at ice base (ppt)
         !*** if dEdd scheme is used, we must choose a smaller hs_min
         !*** because dEdd always computes a non-zero SSWABS which is included
         !*** in FSWINT. if hs_min is too big (>1.e-3), l_snow will be 
         !*** false and SSWABS will not be included in the rhs of the tridiagonal
         !*** system. hence there is inconsitency when computing the energy change
         !*** rate and energy input which causes non-convergence in 
         !*** temperature_change routine.
         !*** the non-convergence is due to condition (5) NOT being satified;
         !*** the flux error is extactly the value of SSWABS, hence if SSWABS 
         !*** is > 0.9*ferrmax, then temperature iteration will never meet
         !*** condition (5)
         !*** the vanilla version of CICE uses 1.e-4 
         !***  
         hs_min = 1.e-4_dbl_kind, & ! min snow thickness for computing Tsno (m)
         betak   = 0.13_dbl_kind, & ! constant in formula for k (W m-1 ppt-1)
         kimin   = 0.10_dbl_kind    ! min conductivity of saline ice (W m-1 deg-1)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         salin       , & ! salinity (ppt)   
         Tmlt            ! melting temp, -depressT * salinity
                         ! nilyr + 1 index is for bottom surface

      real (kind=dbl_kind) :: &
         ustar_scale     ! scaling for ice-ocean heat flux

      !*** move into namelist
      real (kind=dbl_kind) :: &
         ustar_min 

      real (kind=dbl_kind), parameter, private :: &
#ifdef GEOS
#ifdef USE_R8
         ferrmax = 1.0e-3_dbl_kind    ! max allowed energy flux error (W m-2)
#else
         ! have to increase to 1.0 for the single precision  to work
         ferrmax = 1.0e-2_dbl_kind       ! max allowed energy flux error (W m-2)
#endif 
#else
         ferrmax = 1.0e-3_dbl_kind    ! max allowed energy flux error (W m-2)
#endif
                                      ! recommend ferrmax < 0.01 W m-2

      character (char_len) :: stoplabel
 
      character (char_len) :: &
         conduct         ! 'MU71' or 'bubbly'

      logical (kind=log_kind) :: &
         l_brine         ! if true, treat brine pocket effects

      logical (kind=log_kind) :: &
         heat_capacity, &! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics
         calc_Tsfc       ! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE

#ifdef GEOS
      real (kind=dbl_kind)  :: &
         ksno
#endif

      integer (kind=int_kind), parameter :: &
         DIRICHLET = 1,                     &
         NEUMANN   = 2

      character (char_len)  ::              &
         top_bc  = 'flux'  ! default: 'flux', i.e., Neumann BC, this
                           !          turns on the old behavior when
                           !          calc_Tsfc = .false. 
                           !          'mixed', mixed Dirichlet/Neumann BC
                           !          when not converging or Tsfc = 0C
                           !          switch to Dirichlet BC      


!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: alloc_therm_vertical - allocate ice_therm_vertical arrays
!
! !INTERFACE:
!
      subroutine alloc_therm_vertical
!
! !DESCRIPTION:
!
! Allocate therm_vertical arrays
!
! !REVISION HISTORY:
!
! author  Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         allocate(salin(nilyr+1), source=0.0_dbl_kind)
         allocate(Tmlt (nilyr+1), source=0.0_dbl_kind)

      end subroutine alloc_therm_vertical

!=======================================================================
!
!BOP
!
! !IROUTINE: dealloc_therm_vertical - deallocate ice_therm_vertical arrays
!
! !INTERFACE:
!
      subroutine dealloc_therm_vertical
!
! !DESCRIPTION:
!
! Deallocate therm_vertical arrays
!
! !REVISION HISTORY:
!
! author  Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         deallocate(salin)
         deallocate(Tmlt )

      end subroutine dealloc_therm_vertical

!=======================================================================
!BOP
!
! !ROUTINE: thermo_vertical - driver for pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and atmospheric fluxes.
!
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!
! !INTERFACE:


      subroutine thermo_vertical (nx_block,    ny_block,  &
                                  dt,          icells,    &
                                  indxi,       indxj,     &
                                  aicen,       trcrn,     &
                                  vicen,       vsnon,     &
                                  eicen,       esnon,     &
                                  flw,         potT,      &
                                  Qa,          rhoa,      &
                                  fsnow,                  &
                                  fbot,        Tbot,      &
                                  lhcoef,      shcoef,    &
                                  fswsfc,      fswint,    &
                                  fswthrun,               &
                                  Sswabs,      Iswabs,    &
                                  fsurfn,      fcondtopn, &
                                  fsensn,      flatn,     &
                                  fswabsn,     flwoutn,   &
                                  evapn,       freshn,    &
                                  fsaltn,      fhocnn,    &
                                  meltt,       melts,     &
                                  meltb,                  &
                                  congel,      snoice,    &
#ifdef GEOS
                                  DFSDT,DSHDT,DLHDT,DLWDT,&
                                  tlat, tlon, observe,    &
                                  fcondbotl,  sblx,       &            
                                  fcondtopn_repar,        &
#endif
                                  mlt_onset,   frz_onset, &
                                  yday,        l_stop,    &
#ifdef GEOS
                                  istop,       jstop,     &
                                  datm)
#else                                  
                                  istop,       jstop)
#endif

! 
! !USES:
!
      use ice_itd, only: ilyr1, slyr1, ilyrn, slyrn
      use ice_state, only: nt_Tsfc, nt_iage
#ifdef GEOS
      use ice_flux, only: ice_ref_salinity
#endif
#ifdef GEOS
      integer, parameter :: my_task=0, master_task=0, istep1=1
#else
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only: istep1
      use ice_exit
      use ice_ocean
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! ice state variables
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      ! input from atmosphere
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         flw     , & ! incoming longwave radiation (W/m^2)
         potT    , & ! air potential temperature  (K) 
         Qa      , & ! specific humidity (kg/kg) 
         rhoa    , & ! air density (kg/m^3) 
         fsnow   , & ! snowfall rate (kg m-2 s-1)
         shcoef  , & ! transfer coefficient for sensible heat
         lhcoef      ! transfer coefficient for latent heat

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fswsfc  , & ! SW absorbed at ice/snow surface (W m-2)
         fswint  , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun    ! SW through ice to ocean         (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         Sswabs      ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         Iswabs      ! SW radiation absorbed in ice layers (W m-2)

      ! input from ocean
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fbot    , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot        ! ice bottom surface temperature (deg C)

      ! coupler fluxes to atmosphere
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         fsensn  , & ! sensible heat flux (W/m^2) 
         fswabsn , & ! shortwave flux absorbed in ice and ocean (W/m^2) 
         flwoutn , & ! outgoing longwave radiation (W/m^2) 
         evapn       ! evaporative water flux (kg/m^2/s) 

      ! Note: these are intent out if calc_Tsfc = T, otherwise intent in
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         flatn    , & ! latent heat flux   (W/m^2) 
         fsurfn   , & ! net flux to top surface, excluding fcondtopn
         fcondtopn    ! downward cond flux at top surface (W m-2)

      ! coupler fluxes to ocean
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         freshn  , & ! fresh water flux to ocean (kg/m^2/s)
         fsaltn  , & ! salt flux to ocean (kg/m^2/s)
         fhocnn      ! net heat flux to ocean (W/m^2) 

#ifdef GEOS
      ! derivative for GEOS-5 implementation
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in):: &
         DFSDT   , & ! derivative of total heat flux (W/m^2) 
         DSHDT   , & ! derivative of sensible heat flux (W/m^2) 
         DLHDT   , & ! derivative of latent heat flux (W/m^2) 
         DLWDT       ! derivative of upward longwave flux (W/m^2) 
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in):: &
         tlat, tlon
      logical, dimension (nx_block,ny_block), intent(in):: &
         observe
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         fcondbotl,             &
         fcondtopn_repar,       &  
         sblx
      logical (kind=log_kind), optional, intent(in) :: &
         datm
#endif

      ! diagnostic fields
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout):: &
         meltt    , & ! top ice melt             (m/step-->cm/day) 
         melts    , & ! snow melt                (m/step-->cm/day) 
         meltb    , & ! basal ice melt           (m/step-->cm/day) 
         congel   , & ! basal ice growth         (m/step-->cm/day) 
         snoice   , & ! snow-ice formation       (m/step-->cm/day) 
         mlt_onset, & ! day of year that sfc melting begins 
         frz_onset    ! day of year that freezing begins (congel or frazil) 

      real (kind=dbl_kind), intent(in) :: &
         yday      ! day of year

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, print diagnostics and abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where code aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, m     , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         k           , & ! ice layer index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      real (kind=dbl_kind) :: &
         dhi         , & ! change in ice thickness
         dhs             ! change in snow thickness

! 2D state variables (thickness, temperature, enthalpy)

      real (kind=dbl_kind), dimension (icells) :: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         Tsf         , & ! ice/snow top surface temp, same as Tsfcn (deg C)
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hsn_new     , & ! thickness of new snow (m)
         worki       , & ! local work array
         works           ! local work array

      real (kind=dbl_kind), dimension (icells) :: &
         fcondtopn_save

      logical (kind=log_kind), dimension (icells) :: &
         flag_mixed

      real (kind=dbl_kind) :: &
         tmp

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         qin         , & ! ice layer enthalpy, qin < 0 (J m-3)
         Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         qsn         , & ! snow layer enthalpy, qsn < 0 (J m-3)
         Tsn             ! internal snow layer temperatures

! other 2D flux and energy variables

      real (kind=dbl_kind), dimension (icells) :: &
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         einit       , & ! initial energy of melting (J m-2)
         efinal          ! final energy of melting (J m-2)

! ech: the size of these arrays should be reduced to icells
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Tsfcn, & ! temperature of ice/snow top surface  (C)
         iage     ! ice age (s)

#ifdef GEOS
      logical (kind=log_kind) :: &
         atmos_forcing_specified
#endif
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

#ifdef GEOS
      if(present(datm)) then
          atmos_forcing_specified = datm
      else
          atmos_forcing_specified = .false.
      endif
#endif

      do j=1, ny_block
      do i=1, nx_block
         fswabsn(i,j) = c0
        ! the 3 fluxes below are passed by atm in GEOS-5 mode
        ! they should NOT be set to zero, although declared as OUT  
#ifdef GEOS 
!#ifdef DATAATM
         if(atmos_forcing_specified) then  
            fsensn (i,j) = c0
            flwoutn(i,j) = c0
            evapn  (i,j) = c0
         endif
!#endif
#endif
         freshn (i,j) = c0
         fsaltn (i,j) = c0
         fhocnn (i,j) = c0

         fcondtopn_repar (i,j) = c0

         meltt  (i,j) = c0
         meltb  (i,j) = c0
         melts  (i,j) = c0
         congel (i,j) = c0
         snoice (i,j) = c0

         Tsfcn(i,j) = trcrn(i,j,nt_Tsfc)
         if (tr_iage) iage(i,j) = trcrn(i,j,nt_iage)
      enddo
      enddo

      if (calc_Tsfc) then
         do j=1, ny_block
         do i=1, nx_block
           ! see note above 
#ifdef GEOS 
!#ifdef DATAATM
         if(atmos_forcing_specified) then  
            flatn    (i,j) = c0
            fsurfn   (i,j) = c0
         endif
!#endif
#endif
            fcondtopn(i,j) = c0
         enddo
         enddo
      endif

      !-----------------------------------------------------------------
      ! Compute variables needed for vertical thermo calculation
      !-----------------------------------------------------------------

      call init_vertical_profile (nx_block,     ny_block,     &
                                  my_task,      istep1,       &
                                  icells,                     &
                                  indxi,        indxj,        &
                                  aicen(:,:),                 &
                                  vicen(:,:),   vsnon(:,:),   &
                                  Tsfcn(:,:),                 &
                                  eicen(:,:,:), esnon(:,:,:), &
                                  hin,          hilyr,        &
                                  hsn,          hslyr,        &
                                  qin,          Tin,          &
                                  qsn,          Tsn,          &
                                  Tsf,          einit,        &
#ifdef GEOS
#ifndef USE_R8
                                  esumn,                      &  
#endif
                                  tlat, tlon,                 &
#endif
                                  l_stop,                     &
                                  istop,        jstop)

      if (l_stop) return

      do ij = 1, icells
         ! Save initial ice and snow thickness (for fresh and fsalt)
         worki(ij) = hin(ij)
         works(ij) = hsn(ij)
      enddo
   

      !-----------------------------------------------------------------
      ! Compute new surface temperature and internal ice and snow
      !  temperatures.
      !-----------------------------------------------------------------

      if (heat_capacity) then   ! usual case

         call temperature_changes(nx_block,      ny_block, &
                                  my_task,       istep1,   &
                                  dt,            icells,   & 
                                  indxi,         indxj,    &
                                  rhoa,          flw,      &
                                  potT,          Qa,       &
                                  shcoef,        lhcoef,   &
                                  fswsfc,        fswint,   &
                                  fswthrun,      Sswabs,   &
                                  Iswabs,                  &
                                  fcondtopn_save,          &
                                  flag_mixed,              &
                                  hilyr,         hslyr,    &
                                  qin,           Tin,      &
                                  qsn,           Tsn,      &
                                  Tsf,           Tbot,     &
                                  fsensn,        flatn,    &
                                  fswabsn,       flwoutn,  &
                                  fsurfn,                  &
                                  fcondtopn,     fcondbot, &
#ifdef GEOS
                                  dfsdt, dshdt, dlhdt, dlwdt, &
                                  observe,                 &
                                  atmos_forcing_specified, &
#endif
                                  einit,         l_stop,   &
                                  istop,         jstop)

      else

         if (calc_Tsfc) then       

            call zerolayer_temperature(nx_block,      ny_block, &
                                       my_task,       istep1,   &
                                       dt,            icells,   & 
                                       indxi,         indxj,    &
                                       rhoa,          flw,      &
                                       potT,          Qa,       &
                                       shcoef,        lhcoef,   &
                                       fswsfc,        fswthrun, &
                                       hilyr,         hslyr,    &
                                       Tsf,           Tbot,     &
                                       fsensn,        flatn,    &
                                       fswabsn,       flwoutn,  &
                                       fsurfn,                  &
                                       fcondtopn,     fcondbot, &
                                       l_stop,                  &
                                       istop,         jstop)

         else

            !------------------------------------------------------------
            ! Set fcondbot = fcondtop for zero layer thermodynamics
            ! fcondtop is set in call to set_sfcflux in step_therm1
            !------------------------------------------------------------

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               fcondbot(ij)  = fcondtopn(i,j)   ! zero layer         
            enddo
      
         endif      ! calc_Tsfc

      endif         ! heat_capacity

      if (l_stop) return

#ifdef GEOS
      fcondbotl(1,1) = fcondbot(1)

      if (.not. calc_Tsfc .and. top_bc == 'mixed') then
          do ij = 1, icells
              i = indxi(ij)
              j = indxj(ij)
              if (Tsf(ij) < c0 .and. flag_mixed(ij)) then
                  tmp =  fcondtopn(i,j)
                  fcondtopn(i,j) = fcondtopn_save(ij)
                  fcondtopn_save(ij) = tmp
              endif
          enddo
      endif
#endif

      !-----------------------------------------------------------------
      ! Compute growth and/or melting at the top and bottom surfaces.
      ! Add new snowfall.
      ! Repartition ice into equal-thickness layers, conserving energy.
      !-----------------------------------------------------------------

      call thickness_changes(nx_block,     ny_block, &
                             dt,                     &
                             yday,         icells,   &
                             indxi,        indxj,    &
                             efinal,                 &
                             hin,          hilyr,    &
                             hsn,          hslyr,    &
                             qin,          qsn,      &
                             fbot,         Tbot,     &
                             flatn,        fsurfn,   &
                             fcondtopn,    fcondbot, &
                             fsnow,        hsn_new,  &
                             fhocnn,       evapn,    &
                             meltt,        melts,    &
                             meltb,        iage,     &
                             congel,       snoice,   &
#ifdef GEOS
                             observe,      sblx,     &
#endif
                             mlt_onset,    frz_onset)

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy
      ! to the net energy input
      !-----------------------------------------------------------------

      if (.not. calc_Tsfc .and. top_bc == 'mixed') then
          do ij = 1, icells
              i = indxi(ij)
              j = indxj(ij)
              if (Tsf(ij) < c0 .and. flag_mixed(ij)) then
                  tmp =  fcondtopn(i,j)
                  fcondtopn(i,j) = fcondtopn_save(ij)
                  fcondtopn_save(ij) = tmp
              endif
              if (Tsf(ij) == c0 .and. flag_mixed(ij)) then
                   if (fsurfn(i,j) < fcondtopn(i,j)) then
                      fhocnn(i,j) = fhocnn(i,j) -   &
                                (fcondtopn(i,j)-fsurfn(i,j)) 
                   endif
              endif
          enddo
      endif

      call conservation_check_vthermo(nx_block, ny_block, &
                                      my_task,  istep1,   &
                                      dt,       icells,   &
                                      indxi,    indxj,    &
                                      fsurfn,   flatn,    &
                                      fhocnn,   fswint,   &
                                      fsnow,              &
                                      einit,    efinal,   &
#ifdef GEOS
                                      observe,            & 
                                      fcondtopn,          & 
                                      Tsf,                & 
                                      flag_mixed,         &         
#endif
                                      l_stop,             &
                                      istop,    jstop)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Compute fluxes of water and salt from ice to ocean.
      ! evapn < 0 => sublimation, evapn > 0 => condensation
      !-----------------------------------------------------------------
      ! the reason evapn gets in here is because freshwater
      ! flux computation is based on total ice/snow thickness
      ! changes and these changes include water from atmosphere
      ! due to sublimation/condensation so it needs to be
      ! substracted.    
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
            
         dhi = hin(ij) - worki(ij)
         dhs = hsn(ij) - works(ij)
               
         freshn(i,j) = evapn(i,j) - &
                       (rhoi*dhi + rhos*(dhs-hsn_new(ij))) / dt
         fsaltn(i,j) = -rhoi*dhi*ice_ref_salinity*p001/dt

      enddo                     ! ij

      !-----------------------------------------------------------------
      !  Given the vertical thermo state variables (hin, hsn, Tsf,
      !   qin, qsn,), compute the new ice state variables (vicen, vsnon,
      !   Tsfcn, eicen, esnon).
      !-----------------------------------------------------------------

      call update_state_vthermo(nx_block,     ny_block,   &
                                icells,                   &
                                indxi,        indxj,      &
                                Tbot,         Tsf,        &     
                                hin,          hsn,        &
                                qin,          qsn,        &
                                aicen(:,:),               &
                                vicen(:,:),   vsnon(:,:), &
                                Tsfcn(:,:),               &
                                eicen(:,:,:), esnon(:,:,:))

      !-----------------------------------------------------------------
      !  Repartition surplus top conduvtive flux (if any) to the base
      !  of the ice and added to ice-ocean heat flux  
      !  This is done also in HadGEM3-GC3.1
      !-----------------------------------------------------------------

      if (.not. calc_Tsfc .and. top_bc == 'mixed') then
          do ij = 1, icells
              i = indxi(ij)
              j = indxj(ij)
              if (Tsf(ij) < c0 .and. flag_mixed(ij)) then
                  fhocnn(i,j) =  fhocnn(i,j) + fcondtopn_save(ij) &
                                 - fcondtopn(i,j)
                  fcondtopn_repar(i,j) = fcondtopn_save(ij) &
                                         - fcondtopn(i,j) 
              endif
          enddo
      endif

      !-----------------------------------------------------------------
      ! Reload tracer array
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         trcrn(i,j,nt_Tsfc) = Tsfcn(i,j)
         if (tr_iage) trcrn(i,j,nt_iage) = iage(i,j)
      enddo
      enddo

      end subroutine thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: init_thermo_vertical - initialize salinity and melting temp
!
! !DESCRIPTION:
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_thermo_vertical
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind, &
         min_salin = 0.1_dbl_kind  ! threshold for brine pocket treatment 

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      ! Set l_brine to false for zero layer thermodynamics
      !-----------------------------------------------------------------

      if (saltmax > min_salin .and. heat_capacity) then
         l_brine = .true.
      else
         l_brine = .false.
      endif
      
      write(nu_diag,*) 'top surface BC type: ', trim(top_bc)

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            salin(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
!            salin(k)=saltmax ! for isosaline ice
            Tmlt(k) = -salin(k)*depressT
         enddo
         salin(nilyr+1) = saltmax
         Tmlt(nilyr+1) = -salin(nilyr+1)*depressT
      else
         do k = 1, nilyr+1
            salin(k) = c0
            Tmlt(k) = c0
         enddo
      endif

      ustar_scale = c1           ! for nonzero currents
      !ustar_scale = c10           ! for zero currents (uncoupled to ocean)

      end subroutine init_thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: calculate_Tin_from_qin  - calculate internal ice temperatures
!
! !DESCRIPTION:
!
!  Compute the internal ice temperatures from enthalpy using
!  quadratic formula
!
! !REVISION HISTORY:
!
! !INTERFACE:
!
      function calculate_Tin_from_qin (qin, Tmltk) &
               result(Tin)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         qin   , &              ! enthalpy
         Tmltk                  ! melting temperature at one level
!
! !OUTPUT PARAMETERS
!
     real (kind=dbl_kind) :: &
         Tin                 ! internal temperature
!
!EOP
!
      real (kind=dbl_kind) :: &
         aa1,bb1,cc1         ! quadratic solvers


      if (l_brine) then
         aa1 = cp_ice
         bb1 = (cp_ocn-cp_ice)*Tmltk - qin/rhoi - Lfresh
         cc1 = Lfresh * Tmltk
         Tin =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                         (c2*aa1)

      else                ! fresh ice
         Tin = (Lfresh + qin/rhoi) / cp_ice
      endif

      end function calculate_Tin_from_qin

!=======================================================================
!BOP
!
! !ROUTINE: frzmlt_bottom_lateral - bottom and lateral heat fluxes
!
! !DESCRIPTION:
!
! Adjust frzmlt to account for changes in fhocn since from_coupler.
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! !REVISION HISTORY:
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine frzmlt_bottom_lateral (nx_block, ny_block, &
                                        ilo, ihi, jlo, jhi, &
                                        dt,                 &
#ifdef GEOS
                                        observe,            &
#endif
                                        aice,     frzmlt,   &
                                        eicen,    esnon,    &
                                        sst,      Tf,       &
                                        strocnxT, strocnyT, &
                                        Tbot,     fbot,     &
                                        rside)
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1

! !INPUT/OUTPUT PARAMETERS:

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step

#ifdef GEOS
      logical, dimension (nx_block,ny_block), intent(in):: &
         observe
#endif

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         aice    , & ! ice concentration
         frzmlt  , & ! freezing/melting potential (W/m^2)
         sst     , & ! sea surface temperature (C)
         Tf      , & ! freezing temperature (C)
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntilyr), &
         intent(in) :: &
         eicen       ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntslyr), &
         intent(in) :: &
         esnon       ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: &
         Tbot    , & ! ice bottom surface temperature (deg C)
         fbot    , & ! heat flux to ice bottom  (W/m^2)
         rside       ! fraction of ice that melts laterally
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k              , & ! layer index
         ij             , & ! horizontal index, combines i and j loops
         imelt              ! number of cells with ice melting

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      real (kind=dbl_kind), dimension (:), allocatable :: &
         etot    , & ! total energy in column
         fside       ! lateral heat flux (W/m^2)

      real (kind=dbl_kind) :: &
         deltaT    , & ! SST - Tbot >= 0
         ustar     , & ! skin friction velocity for fbot (m/s)
         wlat      , & ! lateral melt rate (m/s)
         xtmp          ! temporary variable

      ! Parameters for bottom melting

      ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut

      real (kind=dbl_kind), parameter :: &
         cpchr = -cp_ocn*rhow*0.006_dbl_kind
!#ifdef GEOS
!         !*** CCSM values
!         ustar_min = 1.e-3_dbl_kind
!#else
!         ustar_min = 5.e-3_dbl_kind
!#endif


      ! Parameters for lateral melting

      real (kind=dbl_kind), parameter :: &
         floediam = 300.0_dbl_kind, & ! effective floe diameter (m)
         alpha    = 0.66_dbl_kind , & ! constant from Steele (unitless)
         m1 = 1.6e-6_dbl_kind     , & ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
         m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      do j = 1, ny_block
      do i = 1, nx_block
         rside(i,j) = c0
         Tbot (i,j) = Tf(i,j)
         fbot (i,j) = c0
      enddo
         !*** CCSM values
      enddo

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      imelt = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j) > puny .and. frzmlt(i,j) < c0) then ! ice can melt
            imelt = imelt + 1
            indxi(imelt) = i
            indxj(imelt) = j
         endif
      enddo                     ! i
      enddo                     ! j

      allocate(etot (imelt))
      allocate(fside(imelt))

      do ij = 1, imelt  ! cells where ice can melt
         i = indxi(ij)
         j = indxj(ij)

         fside(ij) = c0

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst(i,j)-Tbot(i,j)),c0)

         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2
         ustar = sqrt (sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow)
         ustar = max (ustar,ustar_min*ustar_scale)
#ifdef GEOS
         !if(observe(i,j)) then
         !   print*, ' ustar = ', ustar, 'stress star = ', &
         !         sqrt(sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow) 
         !endif
#endif

         fbot(i,j) = cpchr * deltaT * ustar ! < 0
         fbot(i,j) = max (fbot(i,j), frzmlt(i,j)) ! frzmlt < fbot < 0

!!! uncomment to use all frzmlt for standalone runs
!!!         fbot(i,j) = min (c0, frzmlt(i,j))

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Steele (1992): JGR, 97, 17,729-17,738
      !-----------------------------------------------------------------

         wlat = m1 * deltaT**m2 ! Maykut & Perovich
         rside(i,j) = wlat*dt*pi/(alpha*floediam) ! Steele
         rside(i,j) = max(c0,min(rside(i,j),c1))

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

      do n = 1, ncat

         do ij = 1, imelt
            etot(ij) = c0
         enddo

         ! melting energy/unit area in each column, etot < 0

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, imelt
               i = indxi(ij)
               j = indxj(ij)
               etot(ij) = etot(ij) + esnon(i,j,slyr1(n)+k-1)
            enddo               ! ij
         enddo

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, imelt
               i = indxi(ij)
               j = indxj(ij)
               etot(ij) = etot(ij) + eicen(i,j,ilyr1(n)+k-1)
            enddo               ! ij
         enddo                  ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, imelt
            i = indxi(ij)
            j = indxj(ij)
            ! lateral heat flux
            fside(ij) = fside(ij) + rside(i,j)*etot(ij)/dt ! fside < 0
         enddo                  ! ij

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu


      do ij = 1, imelt
         i = indxi(ij)
         j = indxj(ij)

         xtmp = frzmlt(i,j)/(fbot(i,j) + fside(ij) + puny) 
         xtmp = min(xtmp, c1)
         fbot(i,j)  = fbot(i,j)  * xtmp
         rside(i,j) = rside(i,j) * xtmp
      enddo                     ! ij

      deallocate(etot)
      deallocate(fside)

      end subroutine frzmlt_bottom_lateral

!=======================================================================
!BOP
!
! !ROUTINE: init_vertical_profile - initial thickness, enthalpy, temperature
!
! !DESCRIPTION:
!
! Given the state variables (vicen, vsnon, eicen, esnon, Tsfcn),
! compute variables needed for the vertical thermodynamics
! (hin, hsn, qin, qsn, Tin, Tsn, Tsf).
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine init_vertical_profile(nx_block, ny_block, &
                                       my_task,  istep1,   &
                                       icells,             &
                                       indxi,    indxj,    &
                                       aicen,    vicen,    &
                                       vsnon,    Tsfcn,    &
                                       eicen,    esnon,    &
                                       hin,      hilyr,    &
                                       hsn,      hslyr,    &
                                       qin,      Tin,      &
                                       qsn,      Tsn,      &
                                       Tsf,      einit,    &
#ifdef GEOS
#ifndef USE_R8 
                                       esumn,              & 
#endif
                                       tlat, tlon,         &
#endif
                                       l_stop,             &
                                       istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task           , & ! task number (diagnostic only)
         istep1            , & ! time step index (diagnostic only)
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon , & ! volume per unit area of snow         (m)
         Tsfcn     ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr), &
         intent(in) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr), &
         intent(in) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(icells), intent(out):: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         Tsf         , & ! ice/snow surface temperature, Tsfcn
         einit           ! initial energy of melting (J m-2)

      real (kind=dbl_kind), dimension(icells), intent(out):: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(out) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(out) :: &
         qsn         , & ! snow enthalpy
         Tsn             ! snow temperature
#ifdef GEOS
#ifndef USE_R8
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         esumn 
#endif
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in):: &
         tlat, tlon
#endif
      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         Tmin = -100._dbl_kind ! min allowed internal temperature (deg C)

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         rnslyr,        & ! real(nslyr)
         aa1, bb1, cc1, & ! terms in quadratic formula
         Tmax             ! maximum allowed snow/ice temperature (deg C)

      logical (kind=log_kind) :: &   ! for vector-friendly error checks
         tsno_high   , & ! flag for Tsn > Tmax
         tice_high   , & ! flag for Tin > Tmlt
         tsno_low    , & ! flag for Tsn < Tmin
         tice_low        ! flag for Tin < Tmin

#ifdef GEOS
      real (kind=dbl_kind) :: &
         mypuny
#endif
      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         esn          ! snow enthalpy
      real (kind=dbl_kind), dimension (icells) :: &
         vsn          ! snow temperature
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rnslyr = real(nslyr,kind=dbl_kind)

      do ij = 1, icells
         einit(ij) = c0
      enddo

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

#ifdef GEOS
      mypuny = 1.e-4_dbl_kind
      !mypuny = puny
#endif
      !-----------------------------------------------------------------
      ! Load arrays for vertical thermo calculation.
      !-----------------------------------------------------------------
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness
      ! Initialize internal energy
      !-----------------------------------------------------------------

         Tsf(ij)    = Tsfcn(i,j)
         hin(ij)    = vicen(i,j) / aicen(i,j)
         hsn(ij)    = vsnon(i,j) / aicen(i,j)
         hilyr(ij)    = hin(ij) / real(nilyr,kind=dbl_kind)
         hslyr(ij)    = hsn(ij) / rnslyr
         vsn(ij)    = vsnon(i,j)
         do k = 1, nslyr
           esn(ij,k) = esnon(i,j,k)
         enddo
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Snow enthalpy and maximum allowed snow temperature
      ! If heat_capacity = F, qsn and Tsn are never used.
      !-----------------------------------------------------------------

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      !
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

            if (hslyr(ij) > hs_min/rnslyr .and. heat_capacity) then
               ! qsn, esnon < 0              
               qsn  (ij,k) = esnon(i,j,k)*rnslyr/vsnon(i,j) 
#ifdef GEOS
               Tmax = -qsn(ij,k)*mypuny*rnslyr / &
#else
               Tmax = -qsn(ij,k)*puny*rnslyr / &
#endif
                       (rhos*cp_ice*vsnon(i,j))
            else
               qsn  (ij,k) = -rhos * Lfresh
               Tmax = puny
            endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.
      ! Note: qsn <= -rhos*Lfresh, so Tsn <= 0.
      !-----------------------------------------------------------------
            Tsn(ij,k) = (Lfresh + qsn(ij,k)/rhos)/cp_ice
  
      !-----------------------------------------------------------------
      ! Check for Tsn > Tmax (allowing for roundoff error) and Tsn < Tmin.
      !-----------------------------------------------------------------
            if (Tsn(ij,k) > Tmax) then
               tsno_high = .true.
            elseif (Tsn(ij,k) < Tmin) then
               tsno_low  = .true.
            endif

         enddo                  ! ij
      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! If Tsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

      if (tsno_high .and. heat_capacity) then
         do k = 1, nslyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (hslyr(ij) > hs_min/rnslyr) then
#ifdef GEOS
                  Tmax = -qsn(ij,k)*mypuny*rnslyr / &
#else
                  Tmax = -qsn(ij,k)*puny*rnslyr / &
#endif
                           (rhos*cp_ice*vsnon(i,j))
               else
                  Tmax = puny
               endif

               if (Tsn(ij,k) > Tmax) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, Tsn > Tmax'
                  write(nu_diag,*) 'Tsn=',Tsn(ij,k)
                  write(nu_diag,*) 'Tmax=',Tmax
                  write(nu_diag,*) 'istep1, my_task, i, j:', &
                                    istep1, my_task, i, j
                  write(nu_diag,*) 'qsn',qsn(ij,k)
                  write(nu_diag,*) 'hsn', hsn(ij)
                  write(nu_diag,*) 'esn', esn(ij,k)
                  write(nu_diag,*) 'vsn', vsn(ij)
                  write(nu_diag,*) 'Tsf=', Tsf(ij)
                  write(nu_diag,*) 'Lfresh=', Lfresh
                  write(nu_diag,*) 'rhos=', rhos
                  write(nu_diag,*) 'cp_ice=', cp_ice
#ifdef GEOS
                  write(nu_diag,*) 'lat= ', tlat(i,j), &
                                   'lon= ', tlon(i,j)
#endif
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif

            enddo               ! ij
         enddo                  ! nslyr
      endif                     ! tsno_high

      if (tsno_low .and. heat_capacity) then
         do k = 1, nslyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (Tsn(ij,k) < Tmin) then ! allowing for roundoff error
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, Tsn < Tmin'
                  write(nu_diag,*) 'Tsn=', Tsn(ij,k)
                  write(nu_diag,*) 'Tmin=', Tmin
                  write(nu_diag,*) 'istep1, my_task, i, j:', &
                                    istep1, my_task, i, j
                  write(nu_diag,*) 'qsn', qsn(ij,k)
                  write(nu_diag,*) 'hsn', hsn(ij)
                  write(nu_diag,*) 'esn', esn(ij,k)
                  write(nu_diag,*) 'vsn', vsn(ij)
                  write(nu_diag,*) 'Tsf=', Tsf(ij)
#ifdef GEOS
                  write(nu_diag,*) 'lat= ', tlat(i,j), &
                                   'lon= ', tlon(i,j)
#endif
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif

            enddo               ! ij
         enddo                  ! nslyr
      endif                     ! tsno_low

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells

            if (Tsn(ij,k) > c0) then   ! correct roundoff error
               Tsn(ij,k) = c0
               qsn(ij,k) = -rhos*Lfresh
            endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------
!#ifndef GEOS
!            einit(ij) = einit(ij) + hslyr(ij)*qsn(ij,k)
!#elif defined(USE_R8)    
            einit(ij) = einit(ij) + hslyr(ij)*qsn(ij,k)
!#endif

         enddo                  ! ij
      enddo                     ! nslyr

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute ice enthalpy
      ! If heat_capacity = F, qin and Tin are never used.
      !-----------------------------------------------------------------
            ! qin, eicen < 0
            qin(ij,k) = eicen(i,j,k)*real(nilyr,kind=dbl_kind) &
                        /vicen(i,j)  

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

            if (l_brine) then
               aa1 = cp_ice
               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin(ij,k)/rhoi - Lfresh 
               cc1 = Lfresh * Tmlt(k)
               Tin(ij,k) =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                             (c2*aa1)
               Tmax = Tmlt(k)

            else                ! fresh ice
               Tin(ij,k) = (Lfresh + qin(ij,k)/rhoi) / cp_ice
#ifdef GEOS
               Tmax = -qin(ij,k)*mypuny/(rhos*cp_ice*vicen(i,j))
#else
               Tmax = -qin(ij,k)*puny/(rhos*cp_ice*vicen(i,j))
#endif
                         ! as above for snow
            endif

      !-----------------------------------------------------------------
      ! Check for Tin > Tmax and Tin < Tmin
      !-----------------------------------------------------------------
            if (Tin(ij,k) > Tmax) then
               tice_high = .true.
            elseif (Tin(ij,k) < Tmin) then
               tice_low  = .true.
            endif

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! If Tin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

         if (tice_high .and. heat_capacity) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (l_brine) then
                  Tmax = Tmlt(k)
               else             ! fresh ice
#ifdef GEOS
                  Tmax = -qin(ij,k)*mypuny/(rhos*cp_ice*vicen(i,j))
#else
                  Tmax = -qin(ij,k)*puny/(rhos*cp_ice*vicen(i,j))
#endif
               endif

               if (Tin(ij,k) > Tmax) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, T > Tmax, layer', k
                  write(nu_diag,*) 'Tin=',Tin(ij,k),', Tmax=',Tmax
                  write(nu_diag,*) 'istep1, my_task, i, j:', &
                                    istep1, my_task, i, j
                  write(nu_diag,*) 'qin',qin(ij,k)
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif
            enddo               ! ij
         endif                  ! tice_high

         if (tice_low .and. heat_capacity) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (Tin(ij,k) < Tmin) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo T < Tmin, layer', k
                  write(nu_diag,*) 'Tin =', Tin(ij,k)
                  write(nu_diag,*) 'Tmin =', Tmin
                  write(nu_diag,*) 'istep1, my_task, i, j:', &
                                    istep1, my_task, i, j
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif
            enddo               ! ij
         endif                  ! tice_low

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells

            if (Tin(ij,k) > c0) then ! correct roundoff error
               Tin(ij,k) = c0
               qin(ij,k) = -rhoi*Lfresh
            endif
!#ifndef GEOS            
!            einit(ij) = einit(ij) + hilyr(ij)*qin(ij,k) 
!#elif defined(USE_R8)    
            einit(ij) = einit(ij) + hilyr(ij)*qin(ij,k) 
!#endif
         enddo                  ! ij

      enddo                     ! nilyr

!#ifdef GEOS            
!#ifndef USE_R8
!      do ij = 1, icells
!         i = indxi(ij)
!         j = indxj(ij)
!         einit(ij) = esumn(i,j)
!      enddo      
!#endif
!#endif

      end subroutine init_vertical_profile

!=======================================================================
!BOP
!
! !ROUTINE: temperature_changes  - new vertical temperature profile
!
! !DESCRIPTION:
!
! Compute new surface temperature and internal ice and snow
! temperatures.  Include effects of salinity on sea ice heat
! capacity in a way that conserves energy (Bitz and Lipscomb, 1999).
!
! New temperatures are computed iteratively by solving a tridiagonal
! system of equations; heat capacity is updated with each iteration.
! Finite differencing is backward implicit.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669-15,677.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine temperature_changes (nx_block, ny_block, &
                                      my_task,  istep1,   &
                                      dt,       icells,   & 
                                      indxi,    indxj,    &
                                      rhoa,     flw,      &
                                      potT,     Qa,       &
                                      shcoef,   lhcoef,   &
                                      fswsfc,   fswint,   &
                                      fswthrun, Sswabs,   &
                                      Iswabs,             &
                                      fcondtopn_save,     &
                                      flag_mixed,         &
                                      hilyr,    hslyr,    &
                                      qin,      Tin,      &
                                      qsn,      Tsn,      &
                                      Tsf,      Tbot,     &
                                      fsensn,   flatn,    &
                                      fswabsn,  flwoutn,  &
                                      fsurfn,             &
                                      fcondtopn,fcondbot, &
#ifdef GEOS
                                      dfsdt,dshdt,dlhdt,dlwdt, &
                                      observe,            &   
                                      atmos_forcing_specified, &
#endif
                                      einit,    l_stop,   &
                                      istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task     , & ! task number (diagnostic only)
         istep1      , & ! time step index (diagnostic only)
         icells          ! number of cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun        ! SW through ice to ocean         (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         einit           ! initial energy of melting (J m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(out) :: &
         fcondtopn_save   ! save initial top conductive flux 

      logical (kind=log_kind), dimension (icells), intent(out) :: &
         flag_mixed   ! 

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fswabsn     , & ! shortwave absorbed by ice (W m-2)
         flwoutn         ! upward LW at surface (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout):: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! internal snow layer temperatures

#ifdef GEOS
      ! derivative for GEOS-5 implementation
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in):: &
         DFSDT   , & ! derivative of total heat flux (W/m^2) 
         DSHDT   , & ! derivative of sensible heat flux (W/m^2) 
         DLHDT   , & ! derivative of latent heat flux (W/m^2) 
         DLWDT       ! derivative of upward longwave flux (W/m^2) 
      logical, dimension (nx_block,ny_block), intent(in):: &
         observe
      logical (kind=log_kind),   intent(in) :: &
         atmos_forcing_specified
#endif

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      integer (kind=int_kind), parameter :: &
         nitermax = 50   ! max number of iterations in temperature solver

      integer (kind=int_kind) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), parameter :: &
         Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k           , & ! ice layer index
         niter           ! iteration counter in temperature solver

      integer (kind=int_kind) :: &
         isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed

      integer (kind=int_kind), dimension (icells) :: &
         tbc_type        ! Neumann or Dirichlet

      real (kind=dbl_kind), dimension (:), allocatable :: &
         Tsf_start   , & ! Tsf at start of iteration
         dTsf        , & ! Tsf - Tsf_start
         dTi1        , & ! Ti1(1) - Tin_start(1)
         dfsurf_dT   , & ! derivative of fsurf wrt Tsf
         avg_Tsi     , & ! = 1. if new snow/ice temps avg'd w/starting temps
         enew            ! new energy of melting after temp change (J m-2)

      real (kind=dbl_kind), dimension (icells) :: &
         dTsf_prev   , & ! dTsf from previous iteration
         dTi1_prev   , & ! dTi1 from previous iteration
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT  , & ! deriv of flwout wrt Tsf (W m-2 deg-1)
         dt_rhoi_hlyr    ! dt/(rhoi*hilyr)

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         Tin_init    , & ! Tin at beginning of time step
         Tin_start       ! Tin at start of iteration

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         Tsn_init    , & ! Tsn at beginning of time step
         Tsn_start   , & ! Tsn at start of iteration
         etas            ! dt / (rho * cp * h) for snow layers

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         etai        , & ! dt / (rho * cp * h) for ice layers
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs         , & ! rhs of tri-diagonal matrix equation
         Tmat            ! matrix output temperatures

      real (kind=dbl_kind), dimension(icells,nilyr+nslyr+1):: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)

      real (kind=dbl_kind) :: &
         ci          , & ! specific heat of sea ice (J kg-1 deg-1)
         avg_Tsf     , & ! = 1. if Tsf averaged w/Tsf_start, else = 0.
         ferr        , & ! energy conservation error (W m-2)
         Iswabs_tmp  , & ! energy to melt through fraction frac of layer
         Sswabs_tmp  , & ! same for snow
         frac        , & ! fraction of layer that can be melted through
         dTemp           ! minimum temperature difference for absorption
#ifdef GEOS
      real (kind=dbl_kind) :: &
          enewsave
#endif
      logical (kind=log_kind), dimension (icells) :: &
         converged      ! = true when local solution has converged

      logical (kind=log_kind) :: &
         all_converged  ! = true when all cells have converged

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      nmat = nslyr + nilyr + 1  ! matrix dimension

      all_converged   = .false.

      do ij = 1, icells

         converged (ij) = .false.
         l_snow    (ij) = .false.
         l_cold    (ij) = .true.
         tbc_type  (ij) = NEUMANN ! default Neumann
         fcondbot  (ij) = c0
         dTsf_prev (ij) = c0
         dTi1_prev (ij) = c0
         dfsens_dT (ij) = c0
         dflat_dT  (ij) = c0
         dflwout_dT(ij) = c0  
         dt_rhoi_hlyr(ij) = dt / (rhoi*hilyr(ij))  ! hilyr > 0
         if (hslyr(ij) > hs_min/real(nslyr,kind=dbl_kind)) &
            l_snow(ij) = .true.
      enddo                     ! ij

      do k = 1, nslyr
         do ij = 1, icells
            Tsn_init (ij,k) = Tsn(ij,k) ! beginning of time step
            Tsn_start(ij,k) = Tsn(ij,k) ! beginning of iteration
            if (l_snow(ij)) then
               etas(ij,k) = dt/(rhos*cp_ice*hslyr(ij))
            else
               etas(ij,k) = c0
            endif
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nilyr
         do ij = 1, icells
            Tin_init (ij,k) = Tin(ij,k)   ! beginning of time step
            Tin_start(ij,k) = Tin(ij,k)   ! beginning of iteration
         enddo
      enddo

      if(observe(1,1)) then
          write(nu_diag,*) 'in t_c a:', fsurfn(1,1)
          write(nu_diag,*) 'in t_c b:', fcondtopn(1,1)
          write(nu_diag,*) 'in t_c c:', flatn(1,1)
      endif

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during
      !  subsequent iterations).
      ! Ice and snow interfaces are combined into one array (kh) to
      !  simplify the logic.
      !-----------------------------------------------------------------

      call conductivity (nx_block, ny_block,         &
                         l_snow,   icells,           &
                         indxi,    indxj,    indxij, &
                         hilyr,    hslyr,            &
                         Tin,      kh )

!#ifndef GEOS
      !-----------------------------------------------------------------
      ! Check for excessive absorbed solar radiation that may result in
      ! temperature overshoots. Convergence is particularly difficult
      ! if the starting temperature is already very close to the melting 
      ! temperature and extra energy is added.   In that case, or if the
      ! amount of energy absorbed is greater than the amount needed to
      ! melt through a given fraction of a layer, we put the extra 
      ! energy into the surface.
      ! NOTE: This option is not available if the atmosphere model
      !       has already computed fsurf.  (Unless we adjust fsurf here)
      !-----------------------------------------------------------------
!mclaren: Should there be an if calc_Tsfc statement here then?? 
      if ( calc_Tsfc ) then
      !frac = c1 - puny
      !dTemp = p01
      frac = 0.9_dbl_kind    !from CICE 5.0
      dTemp = 0.02_dbl_kind  !from CICE 5.0
      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            Iswabs_tmp = c0
            if (Tin_init(ij,k) <= Tmlt(k) - dTemp) then
               if (l_brine) then
                  !ci = cp_ice - Lfresh / Tin_init(ij,k)
                  ci = cp_ice - Lfresh * Tmlt(k) / (Tin_init(ij,k)**2) ! from CICE 5.0
                  Iswabs_tmp = min(Iswabs(i,j,k), &
                     frac*(Tmlt(k)-Tin_init(ij,k))*ci/dt_rhoi_hlyr(ij))
               else
                  ci = cp_ice
                  Iswabs_tmp = min(Iswabs(i,j,k), &
                     frac*(       -Tin_init(ij,k))*ci/dt_rhoi_hlyr(ij))
               endif
            endif

            fswsfc(i,j)   = fswsfc(i,j) + (Iswabs(i,j,k) - Iswabs_tmp)
#ifdef GEOS
            fsurfn(i,j)   = fsurfn(i,j) + (Iswabs(i,j,k) - Iswabs_tmp)
#endif
            fswint(i,j)   = fswint(i,j) - (Iswabs(i,j,k) - Iswabs_tmp)
            Iswabs(i,j,k) = Iswabs_tmp

         enddo
      enddo

      do k = 1, nslyr
         do ij = 1, icells
            if (l_snow(ij)) then
               i = indxi(ij)
               j = indxj(ij)

               Sswabs_tmp = c0
               if (Tsn_init(ij,k) <= -dTemp) then
                  Sswabs_tmp = min(Sswabs(i,j,k), &
                          -frac*Tsn_init(ij,k)/etas(ij,k))
               endif

               fswsfc(i,j)   = fswsfc(i,j) + (Sswabs(i,j,k) - Sswabs_tmp)
#ifdef GEOS
               fsurfn(i,j)   = fsurfn(i,j) + (Sswabs(i,j,k) - Sswabs_tmp)
#endif
               fswint(i,j)   = fswint(i,j) - (Sswabs(i,j,k) - Sswabs_tmp)
               Sswabs(i,j,k) = Sswabs_tmp

            endif
         enddo
      enddo
      endif

!lipscomb - This could be done in the shortwave module instead.
!           (Change zerolayer routine also)
      ! absorbed shortwave flux for coupler

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         fswabsn(i,j) = fswsfc(i,j) + fswint(i,j) + fswthrun(i,j)
      enddo



#ifdef GEOS
      if(.not. atmos_forcing_specified) then
         dflat_dT(1)   = dlhdt(1,1)
         dfsens_dT(1)  = dshdt(1,1)
         dflwout_dT(1) = dlwdt(1,1)
      endif
#endif

#ifdef GEOS
      isolve = 0
      do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
          fcondtopn_save(ij) = fcondtopn(i,j) 
          flag_mixed(ij) = .false. 
          isolve = isolve + 1
          indxii(isolve) = i
          indxjj(isolve) = j
          indxij(isolve) = ij
      enddo               ! ij
      allocate(   sbdiag(isolve,nilyr+nslyr+1))
      allocate(     diag(isolve,nilyr+nslyr+1))
      allocate(   spdiag(isolve,nilyr+nslyr+1))
      allocate(      rhs(isolve,nilyr+nslyr+1))
      allocate(     Tmat(isolve,nilyr+nslyr+1))
      allocate(     etai(isolve,nilyr))
      allocate(Tsf_start(isolve))
      allocate(     dTsf(isolve))
      allocate(dfsurf_dT(isolve))
      allocate(  avg_Tsi(isolve))
      allocate(     enew(isolve))
      allocate(     dTi1(isolve))
#endif


      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal energy error.
      !-----------------------------------------------------------------

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Identify cells, if any, where calculation has not converged.
      !-----------------------------------------------------------------

         if (all_converged) then  ! thermo calculation is done
            exit
#ifndef GEOS
         else                     ! identify cells not yet converged
            isolve = 0
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (.not.converged(ij)) then
                  isolve = isolve + 1
                  indxii(isolve) = i
                  indxjj(isolve) = j
                  indxij(isolve) = ij
               endif
            enddo               ! ij
#endif
         endif

#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'iter = ', niter, ' starts'
             endif
#endif
#endif
      !-----------------------------------------------------------------
      ! Allocate and initialize
      !-----------------------------------------------------------------

#ifndef GEOS
         allocate(   sbdiag(isolve,nilyr+nslyr+1))
         allocate(     diag(isolve,nilyr+nslyr+1))
         allocate(   spdiag(isolve,nilyr+nslyr+1))
         allocate(      rhs(isolve,nilyr+nslyr+1))
         allocate(     Tmat(isolve,nilyr+nslyr+1))
         allocate(     etai(isolve,nilyr))
         allocate(Tsf_start(isolve))
         allocate(     dTsf(isolve))
         allocate(dfsurf_dT(isolve))
         allocate(  avg_Tsi(isolve))
         allocate(     enew(isolve))
         allocate(     dTi1(isolve))
#endif

         all_converged = .true.

         do ij = 1, isolve
            m = indxij(ij)
            converged(m)  = .true.
            dfsurf_dT(ij) = c0
            avg_Tsi  (ij) = c0
            enew     (ij) = c0
         enddo

      !-----------------------------------------------------------------
      ! Update specific heat of ice layers.
      ! To ensure energy conservation, the specific heat is a function of
      ! both the starting temperature and the (latest guess for) the
      ! final temperature.
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, isolve
               m = indxij(ij)
               i = indxii(ij)
               j = indxjj(ij)

               if (l_brine) then
                  ci = cp_ice - Lfresh*Tmlt(k) /  &
                                (Tin(m,k)*Tin_init(m,k))
               else
                  ci = cp_ice
               endif
               etai(ij,k) = dt_rhoi_hlyr(m) / ci

            enddo
         enddo

         if (calc_Tsfc) then

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

            !if(isolve>0)then
#ifdef GEOS
!#ifndef DATAATM
      if(.not. atmos_forcing_specified) then
               !*** note dfsurf_dT will be a constant for each iteration
               dfsurf_dT(1)  = dflwout_dT(1) + dfsens_dT(1) + dflat_dT(1)
      endif
              ! if(observe(1,1)) then
              !   write(nu_diag,*) 'set iter = ', niter, &
              !                    'dfurf_dt = ', dfsurf_dT(1)
              !   write(nu_diag,*) 'set iter = ', niter, &
              !                    'dflwout_dT = ', dflwout_dT(1), &
              !                    'dfsens_dT = ', dfsens_dT(1), &
              !                    'dflat_dT = ', dflat_dT(1)    
              ! endif
!#endif
#endif

#ifdef GEOS
!#ifdef DATAATM
      if(atmos_forcing_specified) then
               call surface_fluxes( &
                              nx_block,    ny_block,          &
                              isolve,      icells,            &
                              indxii,      indxjj,    indxij, &
                              Tsf,         fswsfc,            &
                              rhoa,        flw,               &
                              potT,        Qa,                &
                              shcoef,      lhcoef,            &
                              flwoutn,     fsensn,            &
                              flatn,       fsurfn,            &
                              dflwout_dT,  dfsens_dT,         &
                              dflat_dT,    dfsurf_dT)
      endif
!#endif
#endif
            !endif

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtopn.
      ! If fsurfn < fcondtopn and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            if (l_snow(m)) then
               fcondtopn(i,j) = kh(m,1) * (Tsf(m) - Tsn(m,1))
            else
               fcondtopn(i,j) = kh(m,1+nslyr) * (Tsf(m) - Tin(m,1))
            endif

            !if (fsurfn(i,j) < fcondtopn(i,j)) &
            !     Tsf(m) = min (Tsf(m), -puny)
            ! new code from v5.1
            if (Tsf(m) >= c0 .and. fsurfn(i,j) < fcondtopn(i,j)) &
                 Tsf(m) = -puny

      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------

            Tsf_start(ij) = Tsf(m)

            !if (Tsf(m) <= -puny) then
            if (Tsf(m) < c0) then
               l_cold(m) = .true.
            else
               l_cold(m) = .false.
            endif
          enddo                  ! ij
#ifdef GEOS
#ifdef DEBUG
       if(observe(1,1)) then
          write(nu_diag,*) 'iter = ', niter, ' before get_matrix'
       endif
#endif
#endif
      !-----------------------------------------------------------------
      ! Compute elements of tridiagonal matrix.
      !-----------------------------------------------------------------

            call get_matrix_elements_calc_Tsfc &
                                  (nx_block, ny_block,         &
                                   isolve,   icells,           &
                                   indxii,   indxjj,   indxij, &
#ifdef GEOS
                                   observe,                    &
#endif
                                   l_snow,   l_cold,           &
                                   Tsf,      Tbot,             &
                                   fsurfn,   dfsurf_dT,        &
                                   Tin_init, Tsn_init,         &
                                   kh,       Sswabs,           &
                                   Iswabs,                     &
                                   etai,     etas,             &
                                   sbdiag,   diag,             &
                                   spdiag,   rhs)

#ifdef GEOS
#ifdef DEBUG
       if(observe(1,1)) then
          write(nu_diag,*) 'iter = ', niter, ' after get_matrix'
       endif
#endif
#endif
         else

             if (top_bc == 'mixed') then
               do ij = 1, isolve
                  i = indxii(ij)
                  j = indxjj(ij)
                  m = indxij(ij)
                  if (Tsf(m) < c0) then
                     l_cold(m) = .true.
                  else
                     l_cold(m) = .false.
                     tbc_type(m) = DIRICHLET
                     flag_mixed(ij) = .true.
                  endif
               enddo
            endif
 
            call get_matrix_elements_know_Tsfc &
                                  (nx_block, ny_block,         &
                                   isolve,   icells,           &
                                   indxii,   indxjj,   indxij, &
                                   l_snow,   Tbot,             &
                                   Tin_init, Tsn_init,         &
                                   kh,       Sswabs,           &
                                   Iswabs,                     &
                                   etai,     etas,             &
                                   sbdiag,   diag,             &
                                   spdiag,   rhs,              &
                                   tbc_type, Tsf,              &
                                   fcondtopn)
         endif  ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------
#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'iter = ', niter
               write(nu_diag,*) 'diag =  ', &
                                (diag(1,k), k=1,nilyr+nslyr+1)
               write(nu_diag,*) 'sbdiag =  ', &
                                (sbdiag(1,k), k=1,nilyr+nslyr+1)
               write(nu_diag,*) 'spdiag =  ', &
                                (spdiag(1,k), k=1,nilyr+nslyr+1)
               write(nu_diag,*) 'rhs =  ', &
                                (rhs(1,k), k=1,nilyr+nslyr+1)
             endif
#endif
#endif

         call tridiag_solver (nx_block, ny_block, &
                              isolve,   icells,   &
                              indxii,   indxjj,   &
                              nmat,     sbdiag,   &
                              diag,     spdiag,   &
                              rhs,      Tmat)

#ifdef GEOS
            !if(observe(1,1)) then
            !     write(nu_diag,*) 'after tridiag iter = ', niter, &
            !                      'dfurf_dt = ', dfsurf_dT(1)
            ! endif
#endif
      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Five conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurfn, must be greater than or equal to the downward
      !        conductive flux, fcondtopn.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        within the prescribed error ferrmax.
      !
      ! For briny ice (the standard case), Tsn and Tin are limited
      !  to prevent them from exceeding their melting temperatures.
      !  (Note that the specific heat formula for briny ice assumes
      !  that T < Tmlt.)  
      ! For fresh ice there is no limiting, since there are cases
      !  when the only convergent solution has Tsn > 0 and/or Tin > 0.
      !  Above-zero temperatures are then reset to zero (with melting 
      !  to conserve energy) in the thickness_changes subroutine.
      !-----------------------------------------------------------------

         if (calc_Tsfc) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tsf from matrix solution
      !-----------------------------------------------------------------

            if (l_cold(m)) then
               if (l_snow(m)) then
                  Tsf(m) = Tmat(ij,1)
               else
                  Tsf(m) = Tmat(ij,1+nslyr)
               endif
            else                ! melting surface
               Tsf(m) = c0
            endif


      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------

            dTsf(ij) = Tsf(m) - Tsf_start(ij)
            avg_Tsf  = c0

#ifdef GEOS
            ! if(observe(1,1)) then
            !   write(nu_diag,*) 'iter = ', niter, 'Tsf = ', Tsf(m), &
            !                    'fsurfn = ', fsurfn(1,1), & 
            !                    'df_dt = ', dfsurf_dT(ij), &
            !                    'dTsf = ', dTsf(ij)
            !                         
            ! endif
#endif
      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0, then average Tsn and Tin to force
      ! internal temps below their melting temps.
      !-----------------------------------------------------------------

            if (Tsf(m) > puny) then
               Tsf(m) = c0
               dTsf(ij) = -Tsf_start(ij)
               if (l_brine) avg_Tsi(ij) = c1   ! avg with starting temp
               converged(m) = .false.
               all_converged = .false.

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      !-----------------------------------------------------------------

            elseif (niter > 1 &                ! condition (2)
              .and. Tsf_start(ij) <= -puny &
              .and. abs(dTsf(ij)) > puny &
              .and. abs(dTsf_prev(m)) > puny &
              .and. -dTsf(ij)/(dTsf_prev(m)+puny*puny) > p5) then

               if (l_brine) then ! average with starting temp
                  avg_Tsf  = c1    
                  avg_Tsi(ij) = c1
               endif
               dTsf(ij) = p5 * dTsf(ij)
               converged(m) = .false.
               all_converged = .false.
            endif

!!!            dTsf_prev(m) = dTsf(ij)

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
            Tsf(m)  = Tsf(m) &
                      + avg_Tsf * p5 * (Tsf_start(ij) - Tsf(m))

          enddo  ! ij

         endif   ! calc_Tsfc

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tsn from matrix solution
      !-----------------------------------------------------------------

               if (l_snow(m)) then
                  Tsn(m,k) = Tmat(ij,k+1)
               else
                  Tsn(m,k) = c0
               endif
               if (l_brine) Tsn(m,k) = min(Tsn(m,k), c0)

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new snow layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tsn(m,k) = Tsn(m,k) &
                         + avg_Tsi(ij)*p5*(Tsn_start(m,k)-Tsn(m,k))

      !-----------------------------------------------------------------
      ! Compute qsn and increment new energy.
      !-----------------------------------------------------------------
               qsn(m,k) = -rhos * (Lfresh - cp_ice*Tsn(m,k))
               enew(ij) = enew(ij) + hslyr(m) * qsn(m,k)

#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'iter = ', niter, 'k = ', k
               write(nu_diag,*) 'qsn = ', qsn(m,k), &
                                'hslyr = ', hslyr(m), & 
                                'enew = ', enew(ij)  
             endif
#endif
#endif
               Tsn_start(m,k) = Tsn(m,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nslyr

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tin from matrix solution
      !-----------------------------------------------------------------
               Tin(m,k) = Tmat(ij,k+1+nslyr)
               if (l_brine) Tin(m,k) = min(Tin(m,k), Tmlt(k))

      !-----------------------------------------------------------------
      ! Condition 2b: check for oscillating Tin(1)
      ! If oscillating, average all ice temps to increase rate of convergence.
      !-----------------------------------------------------------------

               if (k==1 .and. .not.calc_Tsfc) then
                  dTi1(ij) = Tin(m,k) - Tin_start(m,k)

                  if (niter > 1 &                    ! condition 2b    
                      .and. abs(dTi1(ij)) > puny &
                      .and. abs(dTi1_prev(m)) > puny &
                      .and. -dTi1(ij)/(dTi1_prev(m)+puny*puny) > p5) then

                     if (l_brine) avg_Tsi(ij) = c1
                     dTi1(ij) = p5 * dTi1(ij)
                     converged(m) = .false.
                     all_converged = .false.
                  endif
                  dTi1_prev(m) = dTi1(ij)
               endif   ! k = 1 .and. calc_Tsfc = F

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new ice layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tin(m,k) = Tin(m,k) &
                         + avg_Tsi(ij)*p5*(Tin_start(m,k)-Tin(m,k))

      !-----------------------------------------------------------------
      ! Compute qin and increment new energy.
      !-----------------------------------------------------------------
               qin(m,k) = -rhoi * (cp_ice*(Tmlt(k)-Tin(m,k)) &
                                   + Lfresh*(c1-Tmlt(k)/Tin(m,k)) &
                                   - cp_ocn*Tmlt(k))
               enew(ij) = enew(ij) + hilyr(m) * qin(m,k)
#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'iter = ', niter, 'k = ', k
               write(nu_diag,*) 'Tin = ', Tin(m,k), &
                                'qin = ', qin(m,k)
               write(nu_diag,*) 'Tmlt = ', Tmlt(k), &
                                'Tmlt-Tin = ', Tmlt(k)-Tin(m,k)
               write(nu_diag,*) 'hilyr = ', hilyr(m), & 
                                'enew = ', enew(ij)  
             endif
#endif
#endif

               Tin_start(m,k) = Tin(m,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nilyr

         if (calc_Tsfc) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf(ij)) > Tsf_errmax) then
               converged(m) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurfn < fcondtopn with Tsf >= 0
      !-----------------------------------------------------------------

            fsurfn(i,j) = fsurfn(i,j) + dTsf(ij)*dfsurf_dT(ij)
            if (l_snow(m)) then
               fcondtopn(i,j) = kh(m,1) * (Tsf(m)-Tsn(m,1))
            else
               fcondtopn(i,j) = kh(m,1+nslyr) * (Tsf(m)-Tin(m,1))
            endif

            !if (Tsf(m) > -puny .and. fsurfn(i,j) < fcondtopn(i,j)) then
            if (Tsf(m) >= c0 .and. fsurfn(i,j) < fcondtopn(i,j)) then
               converged(m) = .false.
               all_converged = .false.
            endif

            dTsf_prev(m) = dTsf(ij)

          enddo                  ! ij
         endif                   ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Condition 5: check for energy conservation error
      ! Change in internal ice energy should equal net energy input.
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)


            if (.not. calc_Tsfc .and. top_bc == 'mixed') then
               if (tbc_type(m) == DIRICHLET) then
                  if (l_snow(m)) then
                     fcondtopn(i,j) = kh(m,1) * (Tsf(m)-Tsn(m,1))
                 !    print*, 'iter, snow temp at k = 1 ', niter, Tsn(m,1), kh(m,1)
                  else
                     fcondtopn(i,j) = kh(m,1+nslyr) * (Tsf(m)-Tin(m,1))
                 !    print*, 'iter, ice temp at k = 1 ', niter, Tin(m,1), kh(m,1+nslyr)
                  endif
                 ! print*, 'iter, Fcond, Fsurf: ', niter, fcondtopn(1,1), fsurfn(1,1)
               endif
            endif

            fcondbot(m) = kh(m,1+nslyr+nilyr) * &
                           (Tin(m,nilyr)   - Tbot(i,j))

            ferr = abs((enew(ij)-einit(m))/dt &
                 - (fcondtopn(i,j) - fcondbot(m) + fswint(i,j)) )
#ifdef GEOS
            enewsave = enew(ij)
#endif

            ! factor of 0.9 allows for roundoff errors later
            if (ferr > 0.9_dbl_kind*ferrmax) then         ! condition (5)
               converged(m) = .false.
               all_converged = .false.
            endif

            if (.not. calc_Tsfc .and. top_bc == 'mixed') then
               if (tbc_type(m) == NEUMANN .and. niter > 10) then
                  tbc_type(m) = DIRICHLET
                  flag_mixed(ij) = .true.
               endif
            endif

         enddo                  ! ij

#ifdef GEOS
      !if(.not. atmos_forcing_specified .or. calc_Tsfc) then
      !   flwoutn(1,1) = flwoutn(1,1) + dflwout_dT(1)*dTsf(1)
      !   fsensn(1,1)  = fsensn(1,1)  + dfsens_dT(1)*dTsf(1)
      !   flatn(1,1)   = flatn(1,1)   + dflat_dT(1)*dTsf(1)
      !endif
         !*** TODO, the follow should be uncommented?
         !*** should fsurfn be updated in each iteration?  
         !*** it turns out that the following should never be done here
         !*** otherwise, cice will fail with energy conservation error
         !fsurfn(1,1)  = fsurfn(1,1)  + dfsurf_dT(1)*dTsf(1)
#ifdef DEBUG
         if(observe(1,1)) then
            write(nu_diag,*) 'iter = ', niter, 'dTsf = ', dTsf(1)
            write(nu_diag,*) 'dflwout_dT = ', dflwout_dT(1), &
                             'dfsens_dT = ', dfsens_dT(1), & 
                             'dflat_dT = ', dflat_dT(1) 
          endif
#endif
#endif

#ifndef GEOS
         deallocate(sbdiag)
         deallocate(diag)
         deallocate(spdiag)
         deallocate(rhs)
         deallocate(Tmat)
         deallocate(etai)
         deallocate(Tsf_start)
         deallocate(dTsf)
         deallocate(dfsurf_dT)
         deallocate(avg_Tsi)
         deallocate(enew)
         deallocate(dTi1)
#endif

      enddo                     ! temperature iteration niter

#ifdef GEOS
      deallocate(sbdiag)
      deallocate(diag)
      deallocate(spdiag)
      deallocate(rhs)
      deallocate(Tmat)
      deallocate(etai)
      deallocate(Tsf_start)
      deallocate(dTsf)
      deallocate(dfsurf_dT)
      deallocate(avg_Tsi)
      deallocate(enew)
      deallocate(dTi1)
#endif

!#ifdef GEOS
      if(observe(1,1)) then
            write(nu_diag,*) 't_c: flatn: ', flatn(i,j)
      endif
!#endif
      if (.not.all_converged) then

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
            if (.not.converged(ij)) then
               write(nu_diag,*) 'Thermo iteration does not converge,', &
                                'istep1, my_task, i, j:', &
                                 istep1, my_task, i, j
               write(nu_diag,*) 'Ice thickness:',  hilyr(ij)*nilyr
               write(nu_diag,*) 'Snow thickness:', hslyr(ij)*nslyr
               write(nu_diag,*) 'dTsf, Tsf_errmax:',dTsf_prev(ij), &
                                 Tsf_errmax
               write(nu_diag,*) 'Tsf:', Tsf(ij)
               write(nu_diag,*) 'Tbot:', Tbot(i,j)
               write(nu_diag,*) 'fsurf:', fsurfn(i,j)
               write(nu_diag,*) 'fcondtop, fcondbot, fswint', &
                                 fcondtopn(i,j), fcondbot(ij), fswint(i,j)
               write(nu_diag,*) 'fswsfc, fswthrun', &
                                 fswsfc(i,j), fswthrun(i,j)
               write(nu_diag,*) 'Flux conservation error =', ferr
#ifdef GEOS
               write(nu_diag,*) 'enew =', enewsave
               write(nu_diag,*) 'einit =', einit(1)
               write(nu_diag,*) 'enew-einit =', enewsave-einit(1)
               write(nu_diag,*) '(enew-einit)/dt =', &
                                 (enewsave-einit(1))/dt
#endif
               write(nu_diag,*) 'Initial snow temperatures:'
               write(nu_diag,*) (Tsn_init(ij,k),k=1,nslyr)
               write(nu_diag,*) 'Initial ice temperatures:'
               write(nu_diag,*) (Tin_init(ij,k),k=1,nilyr)
               write(nu_diag,*) 'Final snow temperatures:'
               write(nu_diag,*) (Tsn(ij,k),k=1,nslyr)
               write(nu_diag,*) 'Final ice temperatures:'
               write(nu_diag,*) (Tin(ij,k),k=1,nilyr)
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo                  ! ij
      endif                     ! all_converged

#ifdef GEOS
!#ifdef DATAATM
      if(atmos_forcing_specified) then
      if (calc_Tsfc) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            ! update fluxes that depend on Tsf
            flwoutn(i,j) = flwoutn(i,j) + dTsf_prev(ij) * dflwout_dT(ij)
            fsensn(i,j)  = fsensn(i,j)  + dTsf_prev(ij) * dfsens_dT(ij)
            flatn(i,j)   = flatn(i,j)   + dTsf_prev(ij) * dflat_dT(ij)

         enddo                     ! ij
      endif                        ! calc_Tsfc
      endif                        ! calc_Tsfc
!#endif
#endif

      end subroutine temperature_changes

!=======================================================================
!BOP
!
! !ROUTINE: conductivity - compute ice thermal conductivity
!
! !DESCRIPTION:
!
! Compute thermal conductivity at interfaces (held fixed during
!  the subsequent iteration).
!
! NOTE: Ice conductivity must be >= kimin
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conductivity (nx_block, ny_block,         &
                               l_snow,   icells,           &
                               indxi,    indxj,    indxij, &
                               hilyr,    hslyr,            &
                               Tin,      kh)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells          ! number of cells with aicen > puny

      logical (kind=log_kind), dimension(icells), &
         intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hilyr       , & ! ice layer thickness (same for all ice layers)
         hslyr           ! snow layer thickness (same for all snow layers)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (icells,nilyr+nslyr+1), &
         intent(out) :: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k               ! vertical index

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         kilyr           ! thermal cond at ice layer midpoints (W m-1 deg-1)

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         kslyr           ! thermal cond at snow layer midpoints (W m-1 deg-1)

      ! interior snow layers (simple for now, but may be fancier later)
      do k = 1, nslyr
         do ij = 1, icells
            kslyr(ij,k) = ksno
         enddo
      enddo                     ! nslyr

      ! interior ice layers
      if (conduct == 'MU71') then 
        ! Maykut and Untersteiner 1971 form (with Wettlaufer 1991 constants)
        do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            kilyr(ij,k) = kice + betak*salin(k)/min(-puny,Tin(ij,k))
            kilyr(ij,k) = max (kilyr(ij,k), kimin)
         enddo
        enddo                     ! nilyr
      else
      ! Pringle et al JGR 2007 'bubbly brine'
         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               kilyr(ij,k) = (2.11_dbl_kind - 0.011_dbl_kind*Tin(ij,k) &
                            + 0.09_dbl_kind*salin(k)/min(-puny,Tin(ij,k))) &
                            * rhoi / 917._dbl_kind
               kilyr(ij,k) = max (kilyr(ij,k), kimin)
            enddo
         enddo                     ! nilyr
      endif ! conductivity

      ! top snow interface, top and bottom ice interfaces
      do ij = 1, icells
         ! top of snow layer; top surface of top ice layer
         if (l_snow(ij)) then
            kh(ij,1)       = c2 * kslyr(ij,1) / hslyr(ij)
            kh(ij,1+nslyr) = c2 * kslyr(ij,nslyr) * kilyr(ij,1) / &
                             ( kslyr(ij,nslyr)*hilyr(ij) +  &
                               kilyr(ij,1    )*hslyr(ij) )
         else
            kh(ij,1)       = c0
            kh(ij,1+nslyr) = c2 * kilyr(ij,1) / hilyr(ij)
         endif

         ! bottom surface of bottom ice layer
         kh(ij,1+nslyr+nilyr) = c2 * kilyr(ij,nilyr) / hilyr(ij)

      enddo                     ! ij

      ! interior snow interfaces

      if (nslyr > 1) then
         do k = 2, nslyr
            do ij = 1, icells
               if (l_snow(ij)) then
                  kh(ij,k) = c2 * kslyr(ij,k-1) * kslyr(ij,k) / &
                            ((kslyr(ij,k-1) + kslyr(ij,k))*hslyr(ij))
               else
                  kh(ij,k) = c0
               endif
            enddo                  ! ij
         enddo                     ! nilyr
      endif ! nslyr > 1

      ! interior ice interfaces
      do k = 2, nilyr
         do ij = 1, icells
            kh(ij,k+nslyr) = c2 * kilyr(ij,k-1) * kilyr(ij,k) / &
                            ((kilyr(ij,k-1) + kilyr(ij,k))*hilyr(ij))
         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine conductivity


!=======================================================================
!BOP
!
! !ROUTINE: calculate_ki_from_Tin  - calculate ice thermal conductivity
!
! !DESCRIPTION:
!
!  Compute the ice thermal conductivity
!
! !REVISION HISTORY:
!
! !INTERFACE:
!
      function calculate_ki_from_Tin (Tink, salink) &
               result(ki)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         Tink   , &             ! ice layer temperature
         salink                 ! salinity at one level
!
! !OUTPUT PARAMETERS
!
     real (kind=dbl_kind) :: &
         ki                     ! ice conductivity
!
!EOP
!
      if (conduct == 'MU71') then
         ! Maykut and Untersteiner 1971 form (with Wettlaufer 1991 constants)
         ki = kice + betak*salink/min(-puny,Tink)
      else
         ! Pringle et al JGR 2007 'bubbly brine'
         ki = (2.11_dbl_kind - 0.011_dbl_kind*Tink &
             + 0.09_dbl_kind*salink/min(-puny,Tink)) &
             * rhoi / 917._dbl_kind
      endif

      ki = max (ki, kimin)

      end function calculate_ki_from_Tin


!=======================================================================
!BOP
!
! !ROUTINE: surface_fluxes - surface radiative and turbulent fluxes
!
! !DESCRIPTION:
!
! Compute radiative and turbulent fluxes and their derivatives
! with respect to Tsf.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine surface_fluxes (nx_block,   ny_block,          &
                                 isolve,     icells,            &
                                 indxii,     indxjj,    indxij, &
                                 Tsf,        fswsfc,            &
                                 rhoa,       flw,               &
                                 potT,       Qa,                &
                                 shcoef,     lhcoef,            &
                                 flwoutn,    fsensn,            &
                                 flatn,      fsurfn,            &
                                 dflwout_dT, dfsens_dT,         &
                                 dflat_dT,   dfsurf_dT)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxij          ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)

      real (kind=dbl_kind), dimension (isolve), &
         intent(inout) :: &
         dfsurf_dT       ! derivative of fsurfn wrt Tsf
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         TsfK        , & ! ice/snow surface temperature (K)
         Qsfc        , & ! saturated surface specific humidity (kg/kg)
         dQsfcdT     , & ! derivative of Qsfc wrt surface temperature
         qsat        , & ! the saturation humidity of air (kg/m^3)
         flwdabs     , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar          ! 1/TsfK

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         i = indxii(ij)         ! NOTE: not indxi and indxj
         j = indxjj(ij)
         m = indxij(ij)

         ! ice surface temperature in Kelvin
         TsfK = Tsf(m) + Tffresh
         tmpvar = c1/TsfK

         ! saturation humidity
         qsat    = qqqice * exp(-TTTice*tmpvar)
         Qsfc    = qsat / rhoa(i,j)
         dQsfcdT = TTTice * tmpvar*tmpvar * Qsfc

         ! longwave radiative flux
         flwdabs =  emissivity * flw(i,j)
         flwoutn(i,j) = -emissivity * stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsensn(i,j) = shcoef(i,j) * (potT(i,j) - TsfK)
         flatn(i,j)  = lhcoef(i,j) * (Qa(i,j) - Qsfc)

         ! derivatives wrt surface temp
         dflwout_dT(m) = - emissivity*stefan_boltzmann * c4*TsfK**3
         dfsens_dT(m)  = - shcoef(i,j)
         dflat_dT(m)   = - lhcoef(i,j) * dQsfcdT

         fsurfn(i,j) = fswsfc(i,j) + flwdabs + flwoutn(i,j) &
                     + fsensn(i,j) + flatn(i,j)
         dfsurf_dT(ij) = dflwout_dT(m) &
                         + dfsens_dT(m) + dflat_dT(m)

      enddo                     ! ij

      end subroutine surface_fluxes

!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is being computed.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 
!
! !INTERFACE:
!
      subroutine get_matrix_elements_calc_Tsfc &
                                     (nx_block, ny_block,         &
                                      isolve,   icells,           &
                                      indxii,   indxjj,   indxij, &
#ifdef GEOS
                                      observe,                    &   
#endif
                                      l_snow,   l_cold,           &
                                      Tsf,      Tbot,             &
                                      fsurfn,   dfsurf_dT,        &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells), &
         intent(in) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed
#ifdef GEOS
      logical, dimension (nx_block,ny_block), intent(in):: &
         observe
#endif
      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         Tsf             ! ice/snow top surface temp (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn (W/m^2)
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (isolve), intent(in) :: &
         dfsurf_dT       ! derivative of fsurf wrt Tsf

      real (kind=dbl_kind), dimension (isolve,nilyr), &
         intent(in) :: &
         etai            ! dt / (rho*cp*h) for ice layers

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin_init        ! ice temp at beginning of time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         Iswabs          ! absorbed SW flux in ice layers

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (icells,nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (isolve,nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         do ij = 1, isolve
            sbdiag(ij,k) = c0
            diag  (ij,k) = c1
            spdiag(ij,k) = c0
            rhs   (ij,k) = c0
         enddo
      enddo
            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Tsf equation for case of cold surface (with or without snow)
      !-----------------------------------------------------------------
            if (l_cold(m)) then
               if (l_snow(m)) then
                  k = 1
               else                ! no snow
                  k = 1 + nslyr
               endif
               kr = k
            
               sbdiag(ij,kr) = c0
               diag  (ij,kr) = dfsurf_dT(ij) - kh(m,k)
               spdiag(ij,kr) = kh(m,k)
               rhs   (ij,kr) = dfsurf_dT(ij)*Tsf(m) - fsurfn(i,j)
#ifdef GEOS
#ifdef DEBUG
            if(observe(1,1)) then
               write(nu_diag, *) 'A. kr = ', kr, ' rhs = ', rhs(ij,kr), &
                                 'Tsf = ', Tsf(m) 
               write(nu_diag, *) 'dfsurf_dT = ', dfsurf_dT(ij), &
                                 'fsurfn = ',fsurfn(i,j)  
            endif
#endif
#endif
            endif                  ! l_cold

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!           k = 1
!           kr = 2

            if (l_snow(m)) then
               if (l_cold(m)) then
                  sbdiag(ij,2) = -etas(m,1) * kh(m,1)
                  spdiag(ij,2) = -etas(m,1) * kh(m,2)
                  diag  (ij,2) = c1 &
                                + etas(m,1) * (kh(m,1) + kh(m,2))
                  rhs   (ij,2) = Tsn_init(m,1) &
                                + etas(m,1) * Sswabs(i,j,1)
               else                ! melting surface
                  sbdiag(ij,2) = c0
                  spdiag(ij,2) = -etas(m,1) * kh(m,2)
                  diag  (ij,2) = c1 &
                                + etas(m,1) * (kh(m,1) + kh(m,2))
                  rhs   (ij,2) = Tsn_init(m,1) &
                                + etas(m,1)*kh(m,1)*Tsf(m) &
                                + etas(m,1) * Sswabs(i,j,1)
               endif               ! l_cold
            endif                  ! l_snow
#ifdef GEOS
#ifdef DEBUG
            if(observe(1,1)) then
               write(nu_diag, *) 'B. rhs = ', rhs(ij,2), &
                                 'Tsn_init = ', Tsn_init(m,1) 
               write(nu_diag, *) 'etas = ', etas(m,1), &
                                 'Sswabs = ', Sswabs(i,j,1) 
               write(nu_diag, *) 'kh = ', kh(m,1), &
                                 'Tsf = ', Tsf(m)
            endif
#endif
#endif

         enddo                    ! ij

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etas(m,k) * kh(m,k)
                  spdiag(ij,kr) = -etas(m,k) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                               + etas(m,k) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tsn_init(m,k) &
                               + etas(m,k) * Sswabs(i,j,k)
               endif
            enddo               ! ij
         enddo                  ! nslyr

      endif                     ! nslyr > 1


      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m) .or. l_cold(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki) &
                                 + etai(ij,ki)*Iswabs(i,j,ki)
               else    ! no snow, warm surface
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki) &
                                 + etai(ij,ki)*Iswabs(i,j,ki) &
                                 + etai(ij,ki)*kh(m,k)*Tsf(m)
               endif
            
            enddo    ! ij

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
      
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)
            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = c0
            diag  (ij,kr) = c1  &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki) &
                           + etai(ij,ki)*kh(m,k+1)*Tbot(i,j)
#ifdef GEOS
#ifdef DEBUG
            if(observe(1,1)) then
               write(nu_diag, *) 'rhs = ', rhs(ij,kr), &
                                 'Tin_init = ', Tin_init(m,ki) 
               write(nu_diag, *) 'etai = ', etai(ij,ki), &
                                 'Iswabs = ', Iswabs(i,j,ki) 
               write(nu_diag, *) 'kh = ', kh(m,k+1), &
                                 'Tbot = ', Tbot(i,j)
            endif
#endif
#endif

         enddo                   ! ij
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m) .or. l_cold(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               else   ! no snow, warm surface
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k)*Tsf(m)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               endif
            enddo                  ! ij
 
      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
            diag  (ij,kr) = c1 &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki)

         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine get_matrix_elements_calc_Tsfc

!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is already known.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 
!
! !INTERFACE:
!
      subroutine get_matrix_elements_know_Tsfc &
                                     (nx_block, ny_block,         &
                                      isolve,   icells,           &
                                      indxii,   indxjj,   indxij, &
                                      l_snow,   Tbot,             &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      bnd_type, Tsf,              &
                                      fcondtopn)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells), &
         intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (isolve,nilyr), &
         intent(in) :: &
         etai            ! dt / (rho*cp*h) for ice layers

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin_init        ! ice temp at beginning of time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         Iswabs          ! absorbed SW flux in ice layers

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (icells,nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (isolve,nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      integer (kind=int_kind), dimension (icells), &
         intent(in), optional ::                   &
         bnd_type          ! top boundary condition type
                           ! 1: Dirichlet
                           ! 2: Neumann (default)

      real (kind=dbl_kind), dimension (icells), intent(in),  &
         optional :: &
         Tsf             ! ice/snow top surface temp (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in),  &
         optional :: &
         fcondtopn       ! conductive flux at top sfc, positive down (W/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      integer (kind=int_kind), dimension(icells) :: &
         bnd_typ

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         do ij = 1, isolve
            sbdiag(ij,k) = c0
            diag  (ij,k) = c1
            spdiag(ij,k) = c0
            rhs   (ij,k) = c0
         enddo
      enddo

      if (present(bnd_type)) then
         bnd_typ(:) = bnd_type(:)
      else
         bnd_typ(:) = NEUMANN
      endif

            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!        k = 1
!        kr = 2

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            if (l_snow(m)) then
               sbdiag(ij,2) = c0
               spdiag(ij,2) = -etas(m,1) * kh(m,2)
               if(bnd_typ(m) == NEUMANN) then
                  diag  (ij,2) = c1                                 &
                             + etas(m,1) * kh(m,2)
                  rhs   (ij,2) = Tsn_init(m,1)                      &
                             + etas(m,1) * Sswabs(i,j,1)         &
                             + etas(m,1) * fcondtopn(i,j)
               elseif(bnd_typ(m) == DIRICHLET) then    ! melting surface or not converging
                  diag  (ij,2) = c1 &
                                 + etas(m,1) * (kh(m,1) + kh(m,2))
                  rhs   (ij,2) = Tsn_init(m,1) &
                                 + etas(m,1)*kh(m,1)*Tsf(m) &
                                 + etas(m,1) * Sswabs(i,j,1)
               endif   
            endif   ! l_snow
         enddo   ! ij

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etas(m,k) * kh(m,k)
                  spdiag(ij,kr) = -etas(m,k) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                               + etas(m,k) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tsn_init(m,k) &
                               + etas(m,k) * Sswabs(i,j,k)
               endif
            enddo               ! ij
         enddo                  ! nslyr

      endif                     ! nslyr > 1


      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then

                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1                                &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                    &
                                 + etai(ij,ki) * Iswabs(i,j,ki)
               elseif(bnd_typ(m) == NEUMANN) then
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1                                &
                                 + etai(ij,ki) * kh(m,k+1)
                  rhs   (ij,kr) = Tin_init(m,ki)                    &
                                 + etai(ij,ki) * Iswabs(i,j,ki)       &
                                 + etai(ij,ki) * fcondtopn(i,j)
               elseif(bnd_typ(m) == DIRICHLET) then     ! no snow, warm surface
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki) &
                                 + etai(ij,ki)*Iswabs(i,j,ki) &
                                 + etai(ij,ki)*kh(m,k)*Tsf(m)

               endif  ! l_snow
            enddo   ! ij

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
      
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)
            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = c0
            diag  (ij,kr) = c1  &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki) &
                           + etai(ij,ki)*kh(m,k+1)*Tbot(i,j)

         enddo                   ! ij
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               else
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * kh(m,k+1)
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * fcondtopn(i,j)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               endif
            enddo                     ! ij

      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
            diag  (ij,kr) = c1 &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki)

         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine get_matrix_elements_know_Tsfc

!=======================================================================
!BOP
!
! !ROUTINE: tridiag_solver - tridiagonal matrix solver
!
! !DESCRIPTION:
!
! Tridiagonal matrix solver--used to solve the implicit vertical heat
! equation in ice and snow
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine tridiag_solver (nx_block, ny_block, &
                                 isolve,   icells,   &
                                 indxii,   indxjj,   &
                                 nmat,     sbdiag,   &
                                 diag,     spdiag,   &
                                 rhs,      xout)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (isolve,nmat), &
           intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (isolve,nmat), &
           intent(inout) :: &
         xout            ! solution vector
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! row counter

      real (kind=dbl_kind), dimension (isolve) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(isolve,nilyr+nslyr+1):: &
         wgamma          ! temporary matrix variable

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         wbeta(ij) = diag(ij,1)
         xout(ij,1) = rhs(ij,1) / wbeta(ij)
      enddo                     ! ij

      do k = 2, nmat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            wgamma(ij,k) = spdiag(ij,k-1) / wbeta(ij)
            wbeta(ij) = diag(ij,k) - sbdiag(ij,k)*wgamma(ij,k)
            xout(ij,k) = (rhs(ij,k) - sbdiag(ij,k)*xout(ij,k-1)) &
                         / wbeta(ij)
         enddo                  ! ij
      enddo                     ! k

      do k = nmat-1, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            xout(ij,k) = xout(ij,k) - wgamma(ij,k+1)*xout(ij,k+1)
         enddo                  ! ij
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================
!BOP
!
! !ROUTINE: zerolayer_temperature  - new surface temperature calculation
!
! !DESCRIPTION:
!
! Compute new surface temperature using zero layer model of Semtner
! (1976).
!
! New temperatures are computed iteratively by solving a
! surface flux balance equation (i.e. net surface flux from atmos
! equals conductive flux from the top to the bottom surface).
!
! !REVISION HISTORY:
!
! author:  Alison McLaren, Met Office
!         (but largely taken from temperature_changes)
!
! !INTERFACE:
!
      subroutine zerolayer_temperature(nx_block, ny_block, &
                                       my_task,  istep1,   &
                                       dt,       icells,   & 
                                       indxi,    indxj,    &
                                       rhoa,     flw,      &
                                       potT,     Qa,       &
                                       shcoef,   lhcoef,   &
                                       fswsfc,   fswthrun, &
                                       hilyr,    hslyr,    &
                                       Tsf,      Tbot,     &
                                       fsensn,   flatn,    &
                                       fswabsn,  flwoutn,  &
                                       fsurfn,             &
                                       fcondtopn,fcondbot, &
                                       l_stop,             &
                                       istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task     , & ! task number (diagnostic only)
         istep1      , & ! time step index (diagnostic only)
         icells          ! number of cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surface temperature (deg C)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun        ! SW through ice to ocean         (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         fswabsn     , & ! shortwave flux absorbed by ice (W/m-2) 
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout):: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      logical (kind=log_kind), parameter :: &
         l_zerolayerchecks = .true.

      integer (kind=int_kind), parameter :: &
         nitermax = 50   ! max number of iterations in temperature solver

      real (kind=dbl_kind), parameter :: &
         Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         niter           ! iteration counter in temperature solver

      integer (kind=int_kind) :: &
         isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxij          ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (:), allocatable :: &
         Tsf_start   , & ! Tsf at start of iteration
         dTsf        , & ! Tsf - Tsf_start
         dfsurf_dT       ! derivative of fsurfn wrt Tsf

      real (kind=dbl_kind), dimension (icells) :: &
         dTsf_prev   , & ! dTsf from previous iteration
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)

      real (kind=dbl_kind), dimension (:), allocatable :: &
         heff        , & ! effective ice thickness (m)
                         ! ( hice + hsno*kseaice/ksnow)
         kh          , & ! effective conductivity
         diag        , & ! diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix equation

      real (kind=dbl_kind) :: &
         kratio      , & ! ratio of ice and snow conductivies
         avg_Tsf         ! = 1. if Tsf averaged w/Tsf_start, else = 0.

      logical (kind=log_kind), dimension (icells) :: &
         converged      ! = true when local solution has converged

      logical (kind=log_kind) :: &
         all_converged  ! = true when all cells have converged

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      all_converged   = .false.

      do ij = 1, icells
         fcondbot(ij) = c0

         converged (ij) = .false.

         dTsf_prev (ij) = c0

      enddo                     ! ij
      
      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal temperature
      ! change.
      !-----------------------------------------------------------------

      do niter = 1, nitermax

         if (all_converged) then  ! thermo calculation is done
            exit
         else                     ! identify cells not yet converged
            isolve = 0
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (.not.converged(ij)) then
                  isolve = isolve + 1
                  indxii(isolve) = i
                  indxjj(isolve) = j
                  indxij(isolve) = ij
               endif
            enddo               ! ij
         endif

         allocate(     diag(isolve))
         allocate(      rhs(isolve))
         allocate(       kh(isolve))
         allocate(     heff(isolve))
         allocate(Tsf_start(isolve))
         allocate(     dTsf(isolve))
         allocate(dfsurf_dT(isolve))

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

         call surface_fluxes (nx_block,    ny_block,          &
                              isolve,      icells,            &
                              indxii,      indxjj,    indxij, &
                              Tsf,         fswsfc,            &
                              rhoa,        flw,               &
                              potT,        Qa,                &
                              shcoef,      lhcoef,            &
                              flwoutn,     fsensn,            &
                              flatn,       fsurfn,            &
                              dflwout_dT,  dfsens_dT,         &
                              dflat_dT,    dfsurf_dT)

      !-----------------------------------------------------------------
      ! Compute effective ice thickness (includes snow) and thermal 
      ! conductivity 
      !-----------------------------------------------------------------

         kratio = kseaice/ksno
 
         do ij = 1, isolve
             m = indxij(ij)
   
             heff(ij) = hilyr(m) + kratio * hslyr(m)
             kh(ij) = kseaice / heff(ij)        
         enddo                     ! ij


!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtopn.
      ! If fsurfn < fcondtopn and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            fcondtopn(i,j) = kh(ij) * (Tsf(m) - Tbot(i,j))

            if (fsurfn(i,j) < fcondtopn(i,j)) &
                 Tsf(m) = min (Tsf(m), -puny)

      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------

            Tsf_start(ij) = Tsf(m)

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Solve surface balance equation to obtain the new temperatures.
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            diag(ij)  = dfsurf_dT(ij) - kh(ij)
            rhs(ij)   = dfsurf_dT(ij)*Tsf(m) - fsurfn(i,j)   &
                        - kh(ij)*Tbot(i,j)
            Tsf(m)  = rhs(ij) / diag(ij)

         enddo

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Four conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurfn, must be greater than or equal to the downward
      !        conductive flux, fcondtopn.
      !-----------------------------------------------------------------

         ! initialize global convergence flag
         all_converged = .true.

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------

            converged(m) = .true.
            dTsf(ij) = Tsf(m) - Tsf_start(ij)
            avg_Tsf      = c0

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0 and leave converged=.true.
      !-----------------------------------------------------------------

            if (Tsf(m) > puny) then
               Tsf(m) = c0
               dTsf(ij) = -Tsf_start(ij)

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      ! It is possible that this may never occur.
      !-----------------------------------------------------------------

            elseif (niter > 1 &                ! condition (2)
              .and. Tsf_start(ij) <= -puny &
              .and. abs(dTsf(ij)) > puny &
              .and. abs(dTsf_prev(m)) > puny &
              .and. -dTsf(ij)/(dTsf_prev(m)+puny*puny) > p5) then

               avg_Tsf  = c1  ! average with starting temp  
               dTsf(ij) = p5 * dTsf(ij)
               converged(m) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
            Tsf(m)  = Tsf(m) &
                      + avg_Tsf * p5 * (Tsf_start(ij) - Tsf(m))

         enddo  ! ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf(ij)) > Tsf_errmax) then
               converged(m) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurfn < fcondtopn with Tsf > 0
      !-----------------------------------------------------------------

            fsurfn(i,j) = fsurfn(i,j) + dTsf(ij)*dfsurf_dT(ij)
            fcondtopn(i,j) = kh(ij) * (Tsf(m)-Tbot(i,j))

            if (Tsf(m) > -puny .and. fsurfn(i,j) < fcondtopn(i,j)) then
               converged(m) = .false.
               all_converged = .false.
            endif

            fcondbot(m) = fcondtopn(i,j)

            dTsf_prev(m) = dTsf(ij)

         enddo                  ! ij

         deallocate(diag)
         deallocate(rhs)
         deallocate(kh)
         deallocate(heff)
         deallocate(Tsf_start)
         deallocate(dTsf)
         deallocate(dfsurf_dT)

      enddo                     ! temperature iteration niter

      if (.not.all_converged) then

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
            if (.not.converged(ij)) then
               write(nu_diag,*) 'Thermo iteration does not converge,', &
                                'istep1, my_task, i, j:', &
                                 istep1, my_task, i, j
               write(nu_diag,*) 'Ice thickness:',  hilyr(ij)*nilyr
               write(nu_diag,*) 'Snow thickness:', hslyr(ij)*nslyr
               write(nu_diag,*) 'dTsf, Tsf_errmax:',dTsf_prev(ij), &
                                 Tsf_errmax
               write(nu_diag,*) 'Tsf:', Tsf(ij)
               write(nu_diag,*) 'fsurfn:', fsurfn(i,j)
               write(nu_diag,*) 'fcondtopn, fcondbot', &
                                 fcondtopn(i,j), fcondbot(ij)
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo                  ! ij
      endif                     ! all_converged

      !-----------------------------------------------------------------
      ! Check that if Tsfc < 0, then fcondtopn = fsurfn
      !-----------------------------------------------------------------

      if (l_zerolayerchecks) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            
            if (Tsf(ij) < c0 .and. & 
                  abs(fcondtopn(i,j)-fsurfn(i,j)) > puny) then

               write(nu_diag,*) 'fcondtopn does not equal fsurfn,', &
                                'istep1, my_task, i, j:', &
                                 istep1, my_task, i, j
               write(nu_diag,*) 'Tsf=',Tsf(ij)
               write(nu_diag,*) 'fcondtopn=',fcondtopn(i,j)
               write(nu_diag,*) 'fsurfn=',fsurfn(i,j)
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo                  ! ij
      endif                     ! l_zerolayerchecks


!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! update fluxes that depend on Tsf
         flwoutn(i,j) = flwoutn(i,j) + dTsf_prev(ij) * dflwout_dT(ij)
         fsensn(i,j)  = fsensn(i,j)  + dTsf_prev(ij) * dfsens_dT(ij)
         flatn(i,j)   = flatn(i,j)   + dTsf_prev(ij) * dflat_dT(ij)

         ! absorbed shortwave flux for coupler
         fswabsn(i,j) = fswsfc(i,j) + fswthrun(i,j)

      enddo                     ! ij

      end subroutine zerolayer_temperature

!=======================================================================
!BOP
!
! !ROUTINE: thickness changes - top and bottom growth/melting
!
! !DESCRIPTION:
!
! Compute growth and/or melting at the top and bottom surfaces.
! Convert snow to ice if necessary.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine thickness_changes (nx_block,  ny_block, &
                                    dt,                  &
                                    yday,      icells,   &
                                    indxi,     indxj,    &
                                    efinal,              & 
                                    hin,       hilyr,    &
                                    hsn,       hslyr,    &
                                    qin,       qsn,      &
                                    fbot,      Tbot,     &
                                    flatn,     fsurfn,   &
                                    fcondtopn, fcondbot, &
                                    fsnow,     hsn_new,  &
                                    fhocnn,    evapn,    &
                                    meltt,     melts,    &
                                    meltb,     iage,     &
                                    congel,    snoice,   &  
#ifdef GEOS
                                    observe,   sblx,     &
#endif
                                    mlt_onset, frz_onset)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         yday            ! day of the year

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot        , & ! ice bottom surface temperature (deg C)
         fsnow       , & ! snowfall rate (kg m-2 s-1)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(inout) :: &
         qsn             ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         iage        , & ! ice age (s)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         hin         , & ! total ice thickness (m)
         hsn             ! total snow thickness (m)

#ifdef GEOS
      logical, dimension (nx_block,ny_block), intent(in):: &
         observe
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         sblx  
#endif

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         efinal          ! final energy of melting (J m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         fhocnn      , & ! fbot, corrected for any surplus energy (W m-2)
         evapn           ! ice/snow mass sublimated/condensed (kg m-2 s-1)

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         hsn_new         ! thickness of new snow (m)
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         qbotmax = -p5*rhoi*Lfresh  ! max enthalpy of ice growing at bottom

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! vertical index

      real (kind=dbl_kind), dimension (icells) :: &
         esub        , & ! energy for sublimation, > 0    (J m-2)
         econ        , & ! energy for condensation, < 0   (J m-2)
         etop_mlt    , & ! energy for top melting, > 0    (J m-2)
         ebot_mlt    , & ! energy for bottom melting, > 0 (J m-2)
         ebot_gro        ! energy for bottom growth, < 0  (J m-2)

      real (kind=dbl_kind) :: &
         dhi         , & ! change in ice thickness
         dhs         , & ! change in snow thickness
         Ti          , & ! ice temperature
         Ts          , & ! snow temperature
         qbot        , & ! enthalpy of ice growing at bottom surface (J m-3)
         qsub        , & ! energy/unit volume to sublimate ice/snow (J m-3)
         hqtot       , & ! sum of h*q for two layers
         wk1         , & ! temporary variable
         qsnew       , & ! enthalpy of new snow (J m-3)
         hstot           ! snow thickness including new snow (m)

      real (kind=dbl_kind), dimension (icells,nilyr+1) :: &
         zi1         , & ! depth of ice layer boundaries (m)
         zi2             ! adjusted depths, with equal hilyr (m)

      real (kind=dbl_kind), dimension (icells,nslyr+1) :: &
         zs1         , & ! depth of snow layer boundaries (m)
         zs2             ! adjusted depths, with equal hslyr (m)

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         dzi             ! ice layer thickness after growth/melting

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         dzs             ! snow layer thickness after growth/melting

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
#ifdef GEOS
#ifdef DEBUG
         if(observe(1,1)) then
            write(nu_diag,*) 'entering thickness change'
         endif
#endif
#endif

      hsn_new (:) = c0

      do k = 1, nilyr
         do ij = 1, icells
            dzi(ij,k) = hilyr(ij)
         enddo
      enddo

      do k = 1, nslyr
         do ij = 1, icells
            dzs(ij,k) = hslyr(ij)
         enddo
      enddo

      !-----------------------------------------------------------------
      ! For l_brine = false (fresh ice), check for temperatures > 0.
      !  Melt ice or snow as needed to bring temperatures back to 0.
      ! For l_brine = true, this should not be necessary.
      !-----------------------------------------------------------------

      if (.not. l_brine) then 

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells

               Ts = (Lfresh + qsn(ij,k)/rhos) / cp_ice
               if (Ts > c0) then
                  dhs = cp_ice*Ts*dzs(ij,k) / Lfresh
                  dzs(ij,k) = dzs(ij,k) - dhs
                  qsn(ij,k) = -rhos*Lfresh
               endif
            enddo
         enddo

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells

               Ti = (Lfresh + qin(ij,k)/rhoi) / cp_ice
               if (Ti > c0) then
                  dhi = cp_ice*Ti*dzi(ij,k) / Lfresh
                  dzi(ij,k) = dzi(ij,k) - dhi
                  qin(ij,k) = -rhoi*Lfresh
               endif
            enddo               ! ij
         enddo                  ! k

      endif                     ! .not. l_brine

      !-----------------------------------------------------------------
      ! Compute energy available for sublimation/condensation, top melt,
      ! and bottom growth/melt.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         wk1 = -flatn(i,j) * dt
         esub(ij) = max(wk1, c0)     ! energy for sublimation, > 0
         econ(ij) = min(wk1, c0)     ! energy for condensation, < 0

         wk1 = (fsurfn(i,j) - fcondtopn(i,j)) * dt
         etop_mlt(ij) = max(wk1, c0)           ! etop_mlt > 0

         wk1 = (fcondbot(ij) - fbot(i,j)) * dt
         ebot_mlt(ij) = max(wk1, c0)           ! ebot_mlt > 0
         ebot_gro(ij) = min(wk1, c0)           ! ebot_gro < 0

         !--------------------------------------------------------------
         ! Condensation (evapn > 0)
         ! Note: evapn here has unit of kg/m^2.  Divide by dt later.
         !--------------------------------------------------------------

         evapn   (i,j) = c0          ! initialize
#ifdef GEOS
         sblx    (i,j) = c0          ! initialize
#endif

         if (hsn(ij) > puny) then   
#ifdef GEOS
            ! add snow with enthalpy -rLfs (0C)
            dhs = econ(ij) / (-rLfs - rhos*Lvap) ! econ < 0, dhs > 0
#else
            ! add snow with enthalpy qsn(ij,1)
            dhs = econ(ij) / (qsn(ij,1) - rhos*Lvap) ! econ < 0, dhs > 0
#endif
#ifdef GEOS
            hstot = dzs(ij,1) + dhs
            ! adjust top layer snow enthalpy b.c. we added them at 0C
            qsnew = -rLfs
            if (hstot > puny) then
               qsn(ij,1) =  (dzs(ij,1) * qsn(ij,1) &
                          + dhs * qsnew) / hstot
               ! avoid roundoff errors
               qsn(ij,1) = min(qsn(ij,1), -rLfs)
            endif
#endif
            dzs(ij,1) = dzs(ij,1) + dhs
            evapn(i,j) = evapn(i,j) + dhs*rhos
         else         
#ifdef GEOS
            ! add ice with enthalpy -rLfi (fresh at 0C)
            dhi = econ(ij) / (-rLfi - rhoi*Lvap) ! econ < 0, dhi > 0
#else
           ! add ice with enthalpy qin(ij,1)
            dhi = econ(ij) / (qin(ij,1) - rhoi*Lvap) ! econ < 0, dhi > 0
#endif
#ifdef GEOS
            ! adjust top layer ice enthalpy b.c. we added them at 0C
            qsnew = -rLfi
            hqtot = dzi(ij,1)*qin(ij,1) + dhi*qsnew
#endif
            dzi(ij,1) = dzi(ij,1) + dhi
#ifdef GEOS
            if (dzi(ij,1) > puny) &
                qin(ij,1) = hqtot / dzi(ij,1)
#endif
            evapn(i,j) = evapn(i,j) + dhi*rhoi
         endif

         !--------------------------------------------------------------
         ! Grow ice (bottom)
         !--------------------------------------------------------------

         ! enthalpy of new ice growing at bottom surface
         if (heat_capacity) then
           qbot = -rhoi * (cp_ice * (Tmlt(nilyr+1)-Tbot(i,j)) &
                         + Lfresh * (c1-Tmlt(nilyr+1)/Tbot(i,j)) &
                        ! + Lfresh * (c1-Tmlt(nilyr+1)/Tbot(i,j)))
          !*** CCSM has the line below commented
                        - cp_ocn * Tmlt(nilyr+1))
           qbot = min (qbot, qbotmax)      ! in case Tbot is close to Tmlt
         else   ! zero layer
           qbot = -rhoi * Lfresh
         endif

         dhi  = ebot_gro(ij) / qbot     ! dhi > 0
#ifdef GEOS
         !if(observe(1,1)) then
         !   write(nu_diag,*) 'dhi = ', dhi, ' qbot = ', qbot
         !   write(nu_diag,*) 'ebot_gro = ', ebot_gro(ij)
         !   write(nu_diag,*) 'fcondbot = ', fcondbot(ij), ' fbot = ', fbot(i,j)
         !endif
#endif

         hqtot = dzi(ij,nilyr)*qin(ij,nilyr) + dhi*qbot
         dzi(ij,nilyr) = dzi(ij,nilyr) + dhi

         if (dzi(ij,nilyr) > puny) &
              qin(ij,nilyr) = hqtot / dzi(ij,nilyr)

         ! update ice age due to freezing (new ice age = dt)
!         if (tr_iage) &
!            iage(i,j) = (iage(i,j)*hin(ij) + dt*dhi) / (hin(ij) + dhi)

         ! history diagnostics
         congel(i,j) = congel(i,j) + dhi
         if (dhi > puny .and. frz_onset(i,j) < puny) &
                 frz_onset(i,j) = yday

      enddo                     ! ij

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

         !--------------------------------------------------------------
         ! Sublimation of snow (evapn < 0)
         !--------------------------------------------------------------
#ifdef GEOS
            qsub = -rLfs - rhos*Lvap ! qsub < 0
#else
            qsub = qsn(ij,k) - rhos*Lvap ! qsub < 0
#endif
            dhs  = max (-dzs(ij,k), esub(ij)/qsub)  ! esub > 0, dhs < 0
#ifdef GEOS
            sblx(i,j) = sblx(i,j) + (-dhs)*min(qsn(ij,k) - (-rLfs), c0) ! sblx < 0 (J m-2) 
#endif
            dzs(ij,k) = dzs(ij,k) + dhs
            esub(ij) = esub(ij) - dhs*qsub
            esub(ij) = max(esub(ij), c0)   ! in case of roundoff error
            evapn(i,j) = evapn(i,j) + dhs*rhos

         !--------------------------------------------------------------
         ! Melt snow (top)
         !--------------------------------------------------------------

            dhs = max(-dzs(ij,k), etop_mlt(ij)/qsn(ij,k))
            dzs(ij,k) = dzs(ij,k) + dhs         ! qsn < 0, dhs < 0
            etop_mlt(ij) = etop_mlt(ij) - dhs*qsn(ij,k)
            etop_mlt(ij) = max(etop_mlt(ij), c0) ! in case of roundoff error

            ! history diagnostics
            if (dhs < -puny .and. mlt_onset(i,j) < puny) &
               mlt_onset(i,j) = yday
            melts(i,j) = melts(i,j) - dhs

         enddo                  ! ij
      enddo                     ! nslyr

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

         !--------------------------------------------------------------
         ! Sublimation of ice (evapn < 0)
         !--------------------------------------------------------------

#ifdef GEOS
            qsub = -rLfi - rhoi*Lvap                  ! qsub < 0
#else
            qsub = qin(ij,k) - rhoi*Lvap              ! qsub < 0
#endif
            dhi  = max (-dzi(ij,k), esub(ij)/qsub) ! esub < 0, dhi < 0
#ifdef GEOS
            sblx(i,j) = sblx(i,j) + (-dhi)*(qin(ij,k) - (-rLfi)) ! sblx can be v+- (J m-2) 
#endif
            dzi(ij,k) = dzi(ij,k) + dhi
            esub(ij) = esub(ij) - dhi*qsub
            esub(ij) = max(esub(ij), c0)
            evapn(i,j) = evapn(i,j) + dhi*rhoi

         !--------------------------------------------------------------
         ! Melt ice (top)
         !--------------------------------------------------------------

            dhi = max(-dzi(ij,k), etop_mlt(ij)/qin(ij,k))
            dzi(ij,k) = dzi(ij,k) + dhi         ! qin < 0, dhi < 0
            etop_mlt(ij) = etop_mlt(ij) - dhi*qin(ij,k)
            etop_mlt(ij) = max(etop_mlt(ij), c0)

            ! history diagnostics
            if (dhi < -puny .and. mlt_onset(i,j) < puny) &
                 mlt_onset(i,j) = yday
            meltt(i,j) = meltt(i,j) - dhi

         enddo                  ! ij
      enddo                     ! nilyr

      do k = nilyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

         !--------------------------------------------------------------
         ! Melt ice (bottom)
         !--------------------------------------------------------------

            dhi = max(-dzi(ij,k), ebot_mlt(ij)/qin(ij,k))
            dzi(ij,k) = dzi(ij,k) + dhi         ! qin < 0, dhi < 0
            ebot_mlt(ij) = ebot_mlt(ij) - dhi*qin(ij,k)
            ebot_mlt(ij) = max(ebot_mlt(ij), c0)

            ! history diagnostics
            meltb(i,j) = meltb(i,j) - dhi

         enddo                  ! ij
      enddo                     ! nilyr

      do k = nslyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells

         !--------------------------------------------------------------
         ! Melt snow (only if all the ice has melted)
         !--------------------------------------------------------------

            dhs = max(-dzs(ij,k), ebot_mlt(ij)/qsn(ij,k))
            dzs(ij,k) = dzs(ij,k) + dhs         ! qsn < 0, dhs < 0
            ebot_mlt(ij) = ebot_mlt(ij) - dhs*qsn(ij,k)
            ebot_mlt(ij) = max(ebot_mlt(ij), c0)

         enddo                  ! ij
      enddo                     ! nslyr


      !-----------------------------------------------------------------
      ! Compute heat flux used by the ice (<=0).
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         fhocnn(i,j) = fbot(i,j) &
                     + (esub(ij) + etop_mlt(ij) + ebot_mlt(ij))/dt
      enddo

!---!-----------------------------------------------------------------
!---! Add new snowfall at top surface.
!---!-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !----------------------------------------------------------------
      ! NOTE: If heat flux diagnostics are to work, new snow should
      !       have T = 0 (i.e. q = -rhos*Lfresh) and should not be
      !       converted to rain.
      !----------------------------------------------------------------

         if (fsnow(i,j) > c0) then

            hsn_new(ij) = fsnow(i,j)/rhos * dt
            qsnew = -rhos*Lfresh
            hstot = dzs(ij,1) + hsn_new(ij)

            if (hstot > c0) then
               qsn(ij,1) =  (dzs(ij,1) * qsn(ij,1) &
                          + hsn_new(ij) * qsnew) / hstot
               ! avoid roundoff errors
               qsn(ij,1) = min(qsn(ij,1), -rhos*Lfresh)

               dzs(ij,1) = hstot
            endif
         endif

    !-----------------------------------------------------------------
    ! Find the new ice and snow thicknesses.
    !-----------------------------------------------------------------

         hin(ij) = c0
         hsn(ij) = c0
      enddo

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            hin(ij) = hin(ij) + dzi(ij,k)
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            hsn(ij) = hsn(ij) + dzs(ij,k)
         enddo                  ! ij
      enddo                     ! k

     !*** freeboard call is moved after linear_itd call in CCSM
     !*** we do the same here 
#ifndef GEOS
    !-------------------------------------------------------------------
    ! Convert snow to ice if snow lies below freeboard.
    !-------------------------------------------------------------------

      call freeboard (nx_block, ny_block, &
                      icells,             &
                      indxi,    indxj,    &
                      dt,                 &
                      snoice,             &
                      iage,               &
                      hin,      hsn,      &
                      qin,      qsn,      &
                      dzi,      dzs)
#endif

!---!-------------------------------------------------------------------
!---! Repartition the ice and snow into equal-thickness layers,
!---! conserving energy.
!---!-------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute desired layer thicknesses.
      !-----------------------------------------------------------------

      do ij = 1, icells
 
         if (hin(ij) > c0) then
            hilyr(ij) = hin(ij) / real(nilyr,kind=dbl_kind)
         else
            hin(ij) = c0
            hilyr(ij) = c0
         endif
         if (hsn(ij) > c0) then
            hslyr(ij) = hsn(ij) / real(nslyr,kind=dbl_kind)
         else
            hsn(ij) = c0
            hslyr(ij) = c0
         endif

      !-----------------------------------------------------------------
      ! Compute depths zi1 of old layers (unequal thickness).
      ! Compute depths zi2 of new layers (equal thickness).
      !-----------------------------------------------------------------

         zi1(ij,1) = c0
         zi1(ij,1+nilyr) = hin(ij)
 
         zi2(ij,1) = c0
         zi2(ij,1+nilyr) = hin(ij)

      enddo   ! ij

#ifdef GEOS
#ifdef DEBUG
      do k = 1, nilyr
         do ij = 1, icells
             if(observe(1,1)) then
               write(nu_diag,*) 'before adjust ice ','k = ', k
               write(nu_diag,*) 'qin = ', qin(ij,k)
             endif
         enddo                  ! ij
      enddo                     ! k
#endif
#endif

      if (heat_capacity) then

         do k = 1, nilyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               zi1(ij,k+1) = zi1(ij,k) + dzi(ij,k)
               zi2(ij,k+1) = zi2(ij,k) + hilyr(ij)
            end do
         enddo

        !-----------------------------------------------------------------
        ! Conserving energy, compute the enthalpy of the new equal layers.
        !-----------------------------------------------------------------

        call adjust_enthalpy (nx_block, ny_block, &
                              nilyr,    icells,   &
                              indxi,    indxj,    &
                              zi1,      zi2,      &
                              hilyr,    hin,      &
                              qin)

      else ! zero layer (nilyr=1)

         do ij = 1, icells
            qin(ij,1) = -rhoi * Lfresh
            qsn(ij,1) = -rhos * Lfresh
         end do
       
      endif

#ifdef GEOS
#ifdef DEBUG
      do k = 1, nilyr
         do ij = 1, icells
             if(observe(1,1)) then
               write(nu_diag,*) 'after adjust ice ','k = ', k
               write(nu_diag,*) 'qin = ', qin(ij,k)
             endif
         enddo                  ! ij
      enddo                     ! k
#endif
#endif

      if (nslyr > 1) then

      !-----------------------------------------------------------------
      ! Compute depths zs1 of old layers (unequal thickness).
      ! Compute depths zs2 of new layers (equal thickness).
      !-----------------------------------------------------------------

         do ij = 1, icells
            zs1(ij,1) = c0
            zs1(ij,1+nslyr) = hsn(ij)

            zs2(ij,1) = c0
            zs2(ij,1+nslyr) = hsn(ij)
         enddo

         do k = 1, nslyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               zs1(ij,k+1) = zs1(ij,k) + dzs(ij,k)
               zs2(ij,k+1) = zs2(ij,k) + hslyr(ij)
            end do
         enddo

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers.
      !-----------------------------------------------------------------

         call adjust_enthalpy (nx_block, ny_block, &
                               nslyr,    icells,   &
                               indxi,    indxj,    &
                               zs1,      zs2,      &
                               hslyr,    hsn,      &
                               qsn)

      endif   ! nslyr > 1

      !-----------------------------------------------------------------
      ! Compute final ice-snow energy, including the energy of
      !  sublimated/condensed ice.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
#ifdef GEOS
         efinal(ij) = -evapn(i,j)*Lvap + sblx(i,j)
#else
         efinal(ij) = -evapn(i,j)*Lvap
#endif
#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'evap  = ', evapn(i,j)
               write(nu_diag,*) 'Lvap  = ', Lvap
               write(nu_diag,*) 'efinal = ', efinal(ij)  
             endif
#endif
#endif
         evapn(i,j) =  evapn(i,j)/dt
      enddo

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            efinal(ij) = efinal(ij) + hslyr(ij)*qsn(ij,k)
#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'sno ','k = ', k
               write(nu_diag,*) 'qsn = ', qsn(ij,k), &
                                'hslyr = ', hslyr(ij), & 
                                'efinal = ', efinal(ij)  
             endif
#endif
#endif
         enddo                  ! ij
      enddo

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            efinal(ij) = efinal(ij) + hilyr(ij)*qin(ij,k)
#ifdef GEOS
#ifdef DEBUG
             if(observe(1,1)) then
               write(nu_diag,*) 'ice ','k = ', k
               write(nu_diag,*) 'qin = ', qin(ij,k), &
                                'hilyr = ', hilyr(ij), & 
                                'efinal = ', efinal(ij)  
             endif
#endif
#endif
         enddo                  ! ij
      enddo                     ! k

      end subroutine thickness_changes

!=======================================================================
!BOP
!
! !ROUTINE: freeboard - snow-ice conversion
!
! !DESCRIPTION:
!
! If there is enough snow to lower the ice/snow interface below
! sea level, convert enough snow to ice to bring the interface back
! to sea level.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine freeboard (nx_block, ny_block, &
                            icells,             &
                            indxi,    indxj,    &
                            dt,                 &
                            snoice,             &
                            iage,               &
                            hin,      hsn,      &
                            qin,      qsn,      &
                            dzi,      dzs)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         snoice  , & ! snow-ice formation       (m/step-->cm/day)
         iage        ! snow thickness (m)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         hin     , & ! ice thickness (m)
         hsn         ! snow thickness (m)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         qin         ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         dzi         ! ice layer thicknesses (m)

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         qsn         ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(inout) :: &
         dzs         ! snow layer thicknesses (m)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! vertical index

      real (kind=dbl_kind), dimension (icells) :: &
         dhin        , & ! change in ice thickness (m)
         dhsn        , & ! change in snow thickness (m)
         hqs             ! sum of h*q for snow (J m-2)

      real (kind=dbl_kind) :: &
         wk1         , & ! temporary variable
         dhs             ! snow to remove from layer (m)

      !-----------------------------------------------------------------
      ! Determine whether snow lies below freeboard.
      !-----------------------------------------------------------------

      do ij = 1, icells

         dhin(ij) = c0
         dhsn(ij) = c0
         hqs (ij) = c0

         wk1 = hsn(ij) - hin(ij)*(rhow-rhoi)/rhos

         if (wk1 > puny .and. hsn(ij) > puny) then  ! snow below freeboard
            dhsn(ij) = min(wk1*rhoi/rhow, hsn(ij)) ! snow to remove
            dhin(ij) = dhsn(ij) * rhos/rhoi        ! ice to add
         endif
      enddo

      !-----------------------------------------------------------------
      ! Adjust snow layer thickness.
      ! Compute energy to transfer from snow to ice.
      !-----------------------------------------------------------------

      do k = nslyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            if (dhin(ij) > puny) then
               dhs = min(dhsn(ij), dzs(ij,k)) ! snow to remove from layer
               hsn(ij) = hsn(ij) - dhs
               dzs(ij,k) = dzs(ij,k) - dhs
               dhsn(ij) = dhsn(ij) - dhs
               dhsn(ij) = max(dhsn(ij),c0)
               hqs(ij) = hqs(ij) + dhs * qsn(ij,k)
            endif               ! dhin > puny
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Transfer volume and energy from snow to top ice layer.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (dhin(ij) > puny) then
            ! update ice age due to freezing (new ice age = dt)
!            if (tr_iage) &
!               iage(i,j) = (iage(i,j)*hin(ij)+dt*dhin(ij))/(hin(ij)+dhin(ij))

            wk1 = dzi(ij,1) + dhin(ij)
            hin(ij) = hin(ij) + dhin(ij)
            qin(ij,1) = (dzi(ij,1)*qin(ij,1) + hqs(ij)) / wk1
            dzi(ij,1) = wk1

            ! history diagnostic
            snoice(i,j) = snoice(i,j) + dhin(ij)
         endif               ! dhin > puny

      enddo                  ! ij

      end subroutine freeboard

!=======================================================================
!BOP
!
! !ROUTINE: adjust_enthalpy -- enthalpy of new layers
!
! !DESCRIPTION:
!
! Conserving energy, compute the new enthalpy of equal-thickness ice
! or snow layers.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine adjust_enthalpy (nx_block, ny_block, &
                                  nlyr,     icells,   &
                                  indxi,    indxj,    &
                                  z1,       z2,       &
                                  hlyr,     hn,       &
                                  qn)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nlyr              , & ! number of layers (nilyr or nslyr)
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (icells,nlyr+1), &
         intent(in) :: &
         z1          , & ! interface depth for old, unequal layers (m)
         z2              ! interface depth for new, equal layers (m)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hlyr            ! new layer thickness (m)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hn              ! total thickness (m)

      real (kind=dbl_kind), dimension (icells,nlyr), &
         intent(inout) :: &
         qn              ! layer enthalpy (J m-3)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k, k1, k2       ! vertical indices

      real (kind=dbl_kind) :: &
         hovlp           ! overlap between old and new layers (m)

      real (kind=dbl_kind), dimension (icells) :: &
         rhlyr           ! 1./hlyr

      real (kind=dbl_kind), dimension (icells,nlyr) :: &
         hq              ! h * q for a layer

      !-----------------------------------------------------------------
      ! Compute reciprocal layer thickness.
      !-----------------------------------------------------------------

      do ij = 1, icells
         rhlyr(ij) = c0
         if (hn(ij) > puny) rhlyr(ij) = c1 / hlyr(ij)
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute h*q for new layers (k2) given overlap with old layers (k1)
      !-----------------------------------------------------------------

      do k2 = 1, nlyr

         do ij = 1, icells
            hq(ij,k2) = c0
         enddo

         do k1 = 1, nlyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               hovlp = min (z1(ij,k1+1), z2(ij,k2+1)) &
                     - max (z1(ij,k1),   z2(ij,k2))
               hovlp = max (hovlp, c0)

               hq(ij,k2) = hq(ij,k2) + hovlp*qn(ij,k1)
            enddo               ! ij
         enddo                  ! kold
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Compute new enthalpies.
      !-----------------------------------------------------------------

      do k = 1, nlyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            qn(ij,k) = hq(ij,k) * rhlyr(ij)
         enddo                  ! ij
      enddo                     ! k

      end subroutine adjust_enthalpy

!=======================================================================
!BOP
!
! !ROUTINE: conservation_check_vthermo - energy conservation check
!
! !DESCRIPTION:
!
! Check for energy conservation by comparing the change in energy
! to the net energy input.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conservation_check_vthermo(nx_block, ny_block, &
                                            my_task,  istep1,   &
                                            dt,       icells,   &
                                            indxi,    indxj,    &
                                            fsurfn,   flatn,    &
                                            fhocnn,   fswint,   &
                                            fsnow,              &
                                            einit,    efinal,   &
#ifdef GEOS
                                            observe,            &
                                            fcondtopn,          &
                                            Tsf,                &  
                                            flag_mixed,         &   
#endif
                                            l_stop,             &
                                            istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task         , & ! task number (diagnostic only)
         istep1          , & ! time step index (diagnostic only)
         icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         flatn       , & ! surface downward latent heat (W m-2)
         fhocnn      , & ! fbot, corrected for any surplus energy
         fswint      , & ! SW absorbed in ice interior, below surface (W m-2)
         fsnow           ! snowfall rate (kg m-2 s-1)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fcondtopn


      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         einit       , & ! initial energy of melting (J m-2)
         efinal          ! final energy of melting (J m-2)

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

#ifdef GEOS
      logical (kind=log_kind), dimension (nx_block,ny_block), & 
          intent(in) :: &
         observe

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         Tsf

      logical (kind=log_kind), dimension (icells), intent(in) :: &
         flag_mixed
#endif

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij              ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         einp        , & ! energy input during timestep (J m-2)
         ferr            ! energy conservation error (W m-2)

      !----------------------------------------------------------------
      ! If energy is not conserved, print diagnostics and exit.
      !----------------------------------------------------------------
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (flag_mixed(ij) .and. Tsf(ij) < c0) then
            einp = (fcondtopn(i,j) - flatn(i,j) + fswint(i,j) - fhocnn(i,j) &
                    - fsnow(i,j)*Lfresh) * dt
         else  
            einp = (fsurfn(i,j) - flatn(i,j) + fswint(i,j) - fhocnn(i,j) &
                    - fsnow(i,j)*Lfresh) * dt
         endif  
         ferr = abs(efinal(ij)-einit(ij)-einp) / dt
#ifdef GEOS
#ifdef DEBUG
         if(observe(1,1)) then
         write(nu_diag,*) 'Flux error (W/m^2) =', ferr
         write(nu_diag,*) 'Energy error (J) =', ferr*dt
         write(nu_diag,*) 'Initial energy =', einit(ij)
         write(nu_diag,*) 'Final energy =', efinal(ij)
         write(nu_diag,*) 'efinal - einit =', &
                           efinal(ij)-einit(ij)
         write(nu_diag,*) 'Input energy =', einp
         endif
#endif
#endif
         if (ferr > ferrmax) then
            l_stop = .true.
            istop = i
            jstop = j

      !-----------------------------------------------------------------
      ! Note that fsurf - flat = fsw + flw + fsens; i.e., the latent
      ! heat is not included in the energy input, since (efinal - einit)
      ! is the energy change in the system ice + vapor, and the latent
      ! heat lost by the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------

!         einp = (fsurfn(i,j) - flatn(i,j) + fswint(i,j) - fhocnn(i,j) &
!                -fsnow(i,j)*Lfresh) * dt
!         ferr = abs(efinal(ij)-einit(ij)-einp) / dt

         write(nu_diag,*) 'Thermo energy conservation error'
         write(nu_diag,*) 'istep1, my_task, i, j:', &
                           istep1, my_task, i, j
         write(nu_diag,*) 'Flux error (W/m^2) =', ferr
         write(nu_diag,*) 'Energy error (J) =', ferr*dt
         write(nu_diag,*) 'Initial energy =', einit(ij)
         write(nu_diag,*) 'Final energy =', efinal(ij)
         write(nu_diag,*) 'efinal - einit =', &
                           efinal(ij)-einit(ij)
         write(nu_diag,*) 'Input energy =', einp
         write(nu_diag,*) 'fsnow =', Lfresh*fsnow(i,j), &
                          'fhocnn = ',fhocnn(i,j) 
         write(nu_diag,*) 'flatn =', flatn(i,j), 'fswint = ',fswint(i,j)
         write(nu_diag,*) 'fsurfn =', fsurfn(i,j)
         write(nu_diag,*) 'fcondtopn =', fcondtopn(i,j)
         write(nu_diag,*) 'Tsf =', Tsf(ij)
         if(flag_mixed(ij)) then
            write(nu_diag,*) 'Dirichelet BC used'
         endif 

         return
         endif
      enddo

      end subroutine conservation_check_vthermo

!=======================================================================
!BOP
!
! !ROUTINE: update_state_vthermo - new state variables
!
! !DESCRIPTION:
!
! Given the vertical thermo state variables (hin, hsn, qin,
!  qsn, Tsf), compute the new ice state variables (vicen, vsnon,
!  eicen, esnon, Tsfcn).
! Zero out state variables if ice has melted entirely.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine update_state_vthermo (nx_block, ny_block, &
                                       icells,             &
                                       indxi,    indxj,    &
                                       Tf,       Tsf,      &
                                       hin,      hsn,      &
                                       qin,      qsn,      &
                                       aicen,    vicen,    &
                                       vsnon,    Tsfcn,    &
                                       eicen,    esnon)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         Tf              ! freezing temperature (C)

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         qsn             ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         Tsfcn           ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         eicen           ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         esnon           ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! ice layer index

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hin(ij) > c0) then
            ! aicen is already up to date
            vicen(i,j) = aicen(i,j) * hin(ij)
            vsnon(i,j) = aicen(i,j) * hsn(ij)
            Tsfcn(i,j) = Tsf(ij)
         else  ! (hin(ij) == c0)
            aicen(i,j) = c0
            vicen(i,j) = c0
            vsnon(i,j) = c0
            Tsfcn(i,j) = Tf(i,j)
         endif

      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hin(ij) > c0) then
               eicen(i,j,k) = qin(ij,k) * vicen(i,j) &
                                          /real(nilyr,kind=dbl_kind)
            else
               eicen(i,j,k) = c0
            endif

         enddo                  ! ij
      enddo                     ! nilyr

      do k = 1, nslyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hin(ij) > c0) then
               esnon(i,j,k) = qsn(ij,k) * vsnon(i,j) &
                                          /real(nslyr,kind=dbl_kind)
            else
               esnon(i,j,k) = c0
            endif

         enddo                  ! ij
      enddo                     ! nslyr

      end subroutine update_state_vthermo

#ifdef GEOS
      subroutine get_is_interface_temp(tsnint, undef)

      use ice_boundary
      use ice_domain
      use ice_domain_size
      use ice_blocks
      use ice_constants
      use ice_grid,        only: tmask
      use ice_state,       only: aicen, vicen, vsnon, eicen, esnon 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!


      real (kind=real_kind), dimension(:,:),  intent(out) :: & 
         tsnint 

      real (kind=real_kind),  intent(in) :: &
         undef
!
!EOP
!

      real (kind=dbl_kind) :: &
         rnslyr,        & ! real(nslyr)
         aa1, bb1, cc1    ! terms in quadratic formula
         
      
      real (kind=dbl_kind) :: &
         hsn, hin, hslyr, hilyr, tint, qsn, Tsn, qin, Tin 


       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n

       type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:) :: &
         wsn, aiceall !

      logical (kind=log_kind), allocatable, &
         dimension (:,:,:) :: &
         havesno !

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(wsn(nx_block,ny_block,max_blocks))
      allocate(aiceall(nx_block,ny_block,max_blocks))
      allocate(havesno(nx_block,ny_block,max_blocks))

      rnslyr = real(nslyr,kind=dbl_kind)
      wsn(:,:,:) = c0
      havesno(:,:,:) = .true.

      aiceall = sum(aicen, dim=3)
      
      do n=1,ncat

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
           do j = jlo, jhi
           do i = ilo, ihi
              i1 = i - nghost
              j1 = j - nghost
              if(tmask(i,j,iblk) .and. aicen(i,j,n,iblk) > puny) then
                 hsn   = vsnon(i,j,n,iblk) / aicen(i,j,n,iblk)
                 hin   = vicen(i,j,n,iblk) / aicen(i,j,n,iblk)
                 hilyr = p5 * hin / real(nilyr,kind=dbl_kind)
                 hslyr = p5 * hsn / rnslyr
                 if(hsn > hs_min) then
                    ! only bottom snow layer temp. needed 
                    k   = nslyr
                    qsn = esnon(i,j,(n-1)*nslyr+k,iblk)*rnslyr/vsnon(i,j,n,iblk) 
                    Tsn = (Lfresh + qsn/rhos)/cp_ice
                    ! only top ice layer temp. needed 
                    k = 1
                    qin = eicen(i,j,(n-1)*nilyr+k,iblk)*real(nilyr,kind=dbl_kind) &
                                /vicen(i,j,n,iblk)  
                    if (l_brine) then
                        aa1 = cp_ice
                        bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin/rhoi - Lfresh 
                        cc1 = Lfresh * Tmlt(k)
                        Tin =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                             (c2*aa1)

                    else                ! fresh ice
                         Tin = (Lfresh + qin/rhoi) / cp_ice
                    endif
                    ! do linear interpolation b.w. bottom snow layer temp. and 
                    ! top ice layer temp.
                    !tsnint(i1,j1,n) = (Tin * hslyr + Tsn * hilyr) &
                    tint = (Tin * hslyr + Tsn * hilyr) &
                                   / (hilyr + hslyr) 
                    wsn(i,j,iblk) = wsn(i,j,iblk) + aicen(i,j,n,iblk) * tint
                    if(tint > 1.e5) then
                      print*, n, i, j 
                      print*, hilyr, hslyr   
                      print*, Tin, Tsn
                      print*, hin, hsn
                      print*, esnon(i,j,(n-1)*nslyr+k,iblk), eicen(i,j,(n-1)*nilyr+k,iblk) 
                      print*, vicen(i,j,n,iblk), vsnon(i,j,n,iblk), aicen(i,j,n,iblk)
                      stop
                    endif
                 else ! snow depth too small, report missing value
                    !tsnint(i1,j1,n) = undef 
                    havesno(i,j,iblk) = .false.
                 endif                
              else
                 !if(aicen(i,j,n,iblk) <= puny) tsnint(i1,j1,n) = undef 
                 havesno(i,j,iblk) = .false. 
              endif
           enddo 
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
              i1 = i - nghost
              j1 = j - nghost
              if(havesno(i,j,iblk)) then
                 if(aiceall(i,j,iblk) > puny) then
                    tsnint(i1,j1) = wsn(i,j,iblk) / aiceall(i,j,iblk)
                    tsnint(i1,j1) = tsnint(i1,j1) + Tffresh
                 else
                    tsnint(i1,j1) = undef     
                 endif
              else
                 tsnint(i1,j1) = undef     
              endif                     
           enddo 
           enddo 
      enddo 

      deallocate(wsn)
      deallocate(aiceall)
      deallocate(havesno)

      end subroutine get_is_interface_temp 


      subroutine diagnose_internal_ice_temp(vicen, eicen, Tin)


      use ice_itd, only: ilyr1

      real (kind=dbl_kind), intent(in) :: &
            vicen(:), eicen(:,:)  

      real (kind=real_kind), intent(out) :: &
            Tin(:)  

      integer  ::  &
          n, k
      real (kind=dbl_kind) :: &
          qn, tn

      do n = 1, ncat
         if (vicen(n) > puny) then
             do k = 1, nilyr
                 !qn = eicen(ilyr1(n)+k-1) &
                 qn = eicen(k,n) &
                     * real(nilyr,kind=dbl_kind)/vicen(n)  
                 tn = calculate_Tin_from_qin(qn,Tmlt(k))
                 Tin(ilyr1(n)+k-1) = real(tn, kind=real_kind)
             enddo
         endif
      enddo
 
      end subroutine diagnose_internal_ice_temp
#endif
!=======================================================================

      end module ice_therm_vertical

!=======================================================================
