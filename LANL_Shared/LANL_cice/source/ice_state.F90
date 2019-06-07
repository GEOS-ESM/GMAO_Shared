!=======================================================================
!BOP
!
! !MODULE: ice_state - primary state variables
!
! !DESCRIPTION:
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! eicen(i,j,k)         eice(i,j)           J/m^2
! esnon(i,j,k)         esno(i,j)           J/m^2
! trcrn(i,j,it,n)      trcr(i,j,it)        
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of
! m^3, and eice/esno have units of J/m^2 instead of J.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., ice_therm_vertical).
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_state
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_blocks
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension(:,:,:) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         eice  , & ! energy of melt. of ice           (J/m^2)
         esno      ! energy of melt. of snow layer    (J/m^2)

      real (kind=dbl_kind), allocatable, &
         dimension(:,:,:,:) :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow (C)
                   ! 2: meltpond volume                 (m)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension (:,:,:):: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:,:) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:,:) :: &
         apondn , & ! concentration of ponds
         hpondn     ! pond depth         (m)

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:,:,:) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      integer (kind=int_kind), allocatable, dimension (:) :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:,:) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), allocatable, &
         dimension (:,:,:,:) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      !-----------------------------------------------------------------
      ! indices for tracers
      ! The maximum index should be no greater than ntrcr 
      ! (ice_domain_size) to prevent array out-of-bounds errors.
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
         nt_Tsfc  =  1, & ! ice/snow surface temperature
         nt_iage  =  2, & ! volume-weighted ice age
         nt_volpn =  3    ! melt pond volume

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension(:,:,:) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step 
      !-----------------------------------------------------------------

      real (kind=dbl_kind), allocatable, dimension(:,:,:) :: &
         aice_init       ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), allocatable, &
         dimension(:,:,:,:) :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init      ! initial ice volume (m), for linear ITD

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: alloc_state - allocate ice_state variables
!
! !INTERFACE:
!
      subroutine alloc_state
!
! !DESCRIPTION:
!
! Allocate ice_state module variables
!
! !REVISION HISTORY:
!
! author: Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         allocate(aice(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vice(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vsno(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(eice(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(esno(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(trcr(nx_block,ny_block,ntrcr,max_blocks), source=0.0_dbl_kind)


         allocate(aice0(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(aicen(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(vicen(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(vsnon(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

         allocate(apondn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(hpondn(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

         allocate(trcrn(nx_block,ny_block,ntrcr,ncat,max_blocks), source=0.0_dbl_kind)


         allocate(eicen(nx_block,ny_block,ntilyr,max_blocks), source=0.0_dbl_kind)

         allocate(esnon(nx_block,ny_block,ntslyr,max_blocks), source=0.0_dbl_kind)

         allocate(uvel(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(vvel(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(divu(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(shear(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
         allocate(strength(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(aice_init(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

         allocate(aicen_init(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)
         allocate(vicen_init(nx_block,ny_block,ncat,max_blocks), source=0.0_dbl_kind)

     end subroutine alloc_state

!=======================================================================
!BOP
!
! !IROUTINE: dealloc_state - deallocate ice_state variables
!
! !INTERFACE:
!
      subroutine dealloc_state
!
! !DESCRIPTION:
!
! Deallocate ice_state module variables
!
! !REVISION HISTORY:
!
! author: Matthew A. Thompson, NASA/GMAO
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

         deallocate(aice)
         deallocate(vice)
         deallocate(vsno)
         deallocate(eice)
         deallocate(esno)

         deallocate(trcr)

         deallocate(aice0)

         deallocate(aicen)
         deallocate(vicen)
         deallocate(vsnon)

         deallocate(apondn)
         deallocate(hpondn)

         deallocate(trcrn)

         deallocate(eicen)

         deallocate(esnon)

         deallocate(uvel)
         deallocate(vvel)
         deallocate(divu)
         deallocate(shear)
         deallocate(strength)

         deallocate(aice_init)

         deallocate(aicen_init)
         deallocate(vicen_init)

     end subroutine dealloc_state

!=======================================================================
!BOP
!
! !IROUTINE: bound_state - bound calls for ice state variables
!
! !INTERFACE:
!
      subroutine bound_state (aicen, trcrn, &
                              vicen, vsnon, &
                              eicen, esnon)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntrcr,ncat,max_blocks), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntilyr,max_blocks),intent(inout) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntslyr,max_blocks),intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
         call ice_HaloUpdate (aicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (trcrn,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (eicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (esnon,            halo_info, &
                              field_loc_center, field_type_scalar)

      end subroutine bound_state

!=======================================================================
#ifdef GEOS


#ifdef AVG_NB
      subroutine neighborcontribute(i, j, fr, nbs, ai, vi, vs, tr, & 
                                      ei, es, iblk)

      use ice_grid,                only: tmask

       integer (kind=int_kind), intent(in) :: &
         iblk           , & ! block index
         i, j
       integer (kind=int_kind), intent(inout) :: &
         nbs
       real (kind=dbl_kind), dimension(:,:,:), intent(in):: &
         fr
       real (kind=dbl_kind), dimension(:), intent(inout):: &
          ai, vi, vs, ei, es
       real (kind=dbl_kind), dimension(:,:), intent(inout):: &
          tr 
!Local vars
       integer (kind=int_kind) :: &
            k, n

         if(fr(i,j,iblk) > 1.e-6 .and. tmask(i,j,iblk)) then
             nbs = nbs + 1 
             do n=1, ncat 
             ai(n) = ai(n) + aicen(i,j,n,iblk)
             vi(n) = vi(n) + vicen(i,j,n,iblk)
             vs(n) = vs(n) + vsnon(i,j,n,iblk)
             do k=1, ntrcr 
                 tr(k,n) = tr(k,n) + trcrn(i,j,k,n,iblk) 
             enddo
             enddo
             do k=1,ntilyr  
                 ei(k) = ei(k) + eicen(i,j,k,iblk) 
             enddo
             do k=1,ntslyr  
                 es(k) = es(k) + esnon(i,j,k,iblk) 
             enddo
         endif
      end subroutine neighborcontribute

#endif

!BOP
!
!
! !INTERFACE:
!
      subroutine set_ice_state (m_aicen, m_trcrn, &
                                m_vicen, m_vsnon, &
                                m_eicen, m_esnon, &
                                m_uvel,  m_vvel,  &
                                frocean )
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,                only: tmask, umask
      !use ice_itd,                 only: aggregate 
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=real_kind), &
         dimension(:,:,:), intent(in) :: &
         m_aicen , & ! fractional ice area
         m_vicen , & ! volume per unit area of ice          (m)
         m_vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(:,:), intent(in) :: &
         m_uvel, m_vvel

      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
         frocean   ! fractional open water 

      real (kind=real_kind), &
         dimension(:,:,:,:), &
         intent(in) :: &
         m_trcrn     ! ice tracers

      real (kind=real_kind), &
         dimension(:,:,:),intent(in) :: &
         m_eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=real_kind), &
         dimension(:,:,:),intent(in) :: &
         m_esnon     ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block

       !real (kind=real_kind) :: fro  


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
              !if(tmask(i,j,iblk) .and. fro > 1.e-6) then
              if(tmask(i,j,iblk)) then
              do n=1, ncat 
                 aicen(i,j,n,iblk) =  m_aicen(i1,j1,n)
                 vicen(i,j,n,iblk) =  m_vicen(i1,j1,n)
                 vsnon(i,j,n,iblk) =  m_vsnon(i1,j1,n)
                 do k=1, ntrcr 
                   trcrn(i,j,k,n,iblk) = m_trcrn(i1,j1,k,n)
                 enddo
              enddo
              do k=1,ntilyr  
                 eicen(i,j,k,iblk) =  m_eicen(i1,j1,k)
              enddo
              do k=1,ntslyr  
                 esnon(i,j,k,iblk) =  m_esnon(i1,j1,k)
              enddo
              else
              !*** for some reason, some grid points with tmask > 0(indicating ocean)
              !*** do not have associated tiles, such that their ice states CAN NOT be 
              !*** assembled. For these grid points, every state variable, except
              !*** Tsfc, is set to zero. Tsfc is set to ocean freezing temp. 
              do n=1, ncat 
                 aicen(i,j,n,iblk) = c0 
                 vicen(i,j,n,iblk) = c0 
                 vsnon(i,j,n,iblk) = c0 
                 trcrn(i,j,nt_Tsfc,n,iblk) = Tocnfrz 
                 do k=2, ntrcr 
                   trcrn(i,j,k,n,iblk) = c0 
                 enddo
              enddo
              do k=1,ntilyr  
                 eicen(i,j,k,iblk) = c0 
              enddo
              do k=1,ntslyr  
                 esnon(i,j,k,iblk) = c0 
              enddo
              aice0(i,j,iblk) = c1 
              endif
              if(umask(i,j,iblk))  then
                uvel(i,j,iblk)  = m_uvel(i1,j1)
                vvel(i,j,iblk)  = m_vvel(i1,j1)
              else
                uvel(i,j,iblk)  = c0 
                vvel(i,j,iblk)  = c0 
              endif
           enddo 
           enddo 
        enddo 

         call ice_HaloUpdate (aicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (trcrn,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (eicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (esnon,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (uvel,             halo_info, &
                              field_loc_NEcorner, field_type_vector)
         call ice_HaloUpdate (vvel,             halo_info, &
                              field_loc_NEcorner, field_type_vector)

      end subroutine set_ice_state

!=======================================================================
!BOP
!
! !IROUTINE: set_therm_ice_state - set only thermo state vars from 
!                                  saltwater  
!
! !INTERFACE:
      subroutine set_therm_ice_state (m_aicen, m_trcrn, &
                                m_vicen, m_vsnon, &
                                m_eicen, m_esnon)
!
! !DESCRIPTION:
!

! !REVISION HISTORY:
!
! author: Bin Zhao NASA GMAO
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,                only: tmask
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=real_kind), &
         dimension(:,:,:), intent(in) :: &
         m_aicen , & ! fractional ice area
         m_vicen , & ! volume per unit area of ice          (m)
         m_vsnon     ! volume per unit area of snow         (m)

      real (kind=real_kind), &
         dimension(:,:,:,:), &
         intent(in) :: &
         m_trcrn     ! ice tracers

      real (kind=real_kind), &
         dimension(:,:,:),intent(in) :: &
         m_eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=real_kind), &
         dimension(:,:,:),intent(in) :: &
         m_esnon     ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
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
              do n=1, ncat 
                 aicen(i,j,n,iblk) =  m_aicen(i1,j1,n)
                 vicen(i,j,n,iblk) =  m_vicen(i1,j1,n)
                 vsnon(i,j,n,iblk) =  m_vsnon(i1,j1,n)
                 do k=1, ntrcr 
                   trcrn(i,j,k,n,iblk) = m_trcrn(i1,j1,k,n)
                 enddo
              enddo
              do k=1,ntilyr  
                 eicen(i,j,k,iblk) = m_eicen(i1,j1,k)
              enddo
              do k=1,ntslyr  
                 esnon(i,j,k,iblk) = m_esnon(i1,j1,k)
              enddo
              else
              !*** for some reason, some grid points with tmask > 0(indicating ocean)
              !*** do not have associated tiles, such that their ice states CAN NOT be 
              !*** assembled. For these grid points, every state variable, except
              !*** Tsfc, is set to zero. Tsfc is set to ocean freezing temp. 
              do n=1, ncat 
                 aicen(i,j,n,iblk) = c0 
                 vicen(i,j,n,iblk) = c0 
                 vsnon(i,j,n,iblk) = c0 
                 trcrn(i,j,nt_Tsfc,n,iblk) = Tocnfrz 
                 do k=2, ntrcr 
                   trcrn(i,j,k,n,iblk) = c0 
                 enddo
              enddo
              do k=1,ntilyr  
                 eicen(i,j,k,iblk) = c0 
              enddo
              do k=1,ntslyr  
                 esnon(i,j,k,iblk) = c0 
              enddo
              aice0(i,j,iblk) = c1 
              endif
           enddo 
           enddo 
        enddo 

         call ice_HaloUpdate (aicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (trcrn,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (eicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (esnon,            halo_info, &
                              field_loc_center, field_type_scalar)

      end subroutine set_therm_ice_state

!=======================================================================
!BOP
!
! !IROUTINE: set_dyn_ice_state - set only dynamic state vars from 
!                                cicedyna 
!
! !INTERFACE:

      subroutine set_dyn_ice_state (m_uvel,  m_vvel)
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
! author: Bin Zhao NASA GMAO
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,                only: umask
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=dbl_kind), &
         dimension(:,:), intent(in) :: &
         m_uvel, m_vvel  
!
!EOP
!
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
              if(umask(i,j,iblk))  then
                uvel(i,j,iblk)     = m_uvel(i1,j1)
                vvel(i,j,iblk)     = m_vvel(i1,j1)
              else
                uvel(i,j,iblk)     = c0 
                vvel(i,j,iblk)     = c0 
              endif
           enddo 
           enddo 
        enddo 

         call ice_HaloUpdate (uvel,               halo_info, &
                              field_loc_NEcorner, field_type_vector)
         call ice_HaloUpdate (vvel,               halo_info, &
                              field_loc_NEcorner, field_type_vector)

      end subroutine set_dyn_ice_state

!=======================================================================
!BOP
!
!
! !INTERFACE:
!
      subroutine get_ice_state (m_aicen, m_trcrn, &
                                m_vicen, m_vsnon, &
                                m_eicen, m_esnon, &
                                m_uvel,  m_vvel,  &
                                m_divu,  m_shear, &
                                m_strength, frocean)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: to_tgrid, u2tgrid_vector, ANGLET, &
                                 tmask, umask
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=real_kind), &
         dimension(:,:,:), intent(out) :: &
         m_aicen , & ! fractional ice area
         m_vicen , & ! volume per unit area of ice          (m)
         m_vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(:,:), intent(out) :: &
         m_uvel, m_vvel, m_divu,  & ! dynamic fields 
         m_shear, m_strength    

      real (kind=real_kind), &
         dimension(:,:,:,:), &
         intent(out) :: &
         m_trcrn     ! ice tracers

      real (kind=real_kind), &
         dimension(:,:,:),intent(out) :: &
         m_eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=real_kind), &
         dimension(:,:,:),intent(out) :: &
         m_esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=real_kind), &
         dimension(:,:), intent(in) :: &
         frocean   ! fractional open water 

!
!EOP
!
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
              do n=1, ncat 
                 m_aicen(i1,j1,n) = aicen(i,j,n,iblk)
                 m_vicen(i1,j1,n) = vicen(i,j,n,iblk)
                 m_vsnon(i1,j1,n) = vsnon(i,j,n,iblk)
                 do k=1, ntrcr 
                   m_trcrn(i1,j1,k,n) = trcrn(i,j,k,n,iblk) 
                 enddo
              enddo
              do k=1,ntilyr  
                 m_eicen(i1,j1,k) = eicen(i,j,k,iblk) 
              enddo
              do k=1,ntslyr  
                 m_esnon(i1,j1,k) = esnon(i,j,k,iblk)
              enddo
	      m_strength(i1,j1) = strength(i,j,iblk)  
              m_divu(i1,j1)     = divu(i,j,iblk)   
              m_shear(i1,j1)    = shear(i,j,iblk)  
              endif
              if(umask(i,j,iblk)) then
                m_uvel(i1,j1)     = uvel(i,j,iblk)
                m_vvel(i1,j1)     = vvel(i,j,iblk)
              endif
           enddo 
           enddo 
        enddo 

      end subroutine get_ice_state

      subroutine get_aggregate_ice_state (m_vice, m_vsno)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: tmask
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), &
         dimension(:,:), intent(out) :: &
         m_vice , & ! volume per unit area of ice          (m)
         m_vsno     ! volume per unit area of snow         (m)

!
!EOP
!
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
                 m_vice(i1,j1) = vice(i,j,iblk)
                 m_vsno(i1,j1) = vsno(i,j,iblk)
              endif
           enddo 
           enddo 
        enddo 

      end subroutine get_aggregate_ice_state
!=======================================================================
!BOP
!
!
! !INTERFACE:
!
      subroutine get_ice_fr_thickness_tgrid ( &
                                m_hice,  m_aice, undefined_val)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: tmask
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         m_hice, m_aice

      real (kind=real_kind), intent(in) :: &
         undefined_val
         
!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, k, n
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
              if(tmask(i,j,iblk)) then
                m_hice(i-nghost,j-nghost)  =  vice(i,j,iblk)
                m_aice(i-nghost,j-nghost)  =  aice(i,j,iblk)
              else
                m_hice(i-nghost,j-nghost)  =  undefined_val
                m_aice(i-nghost,j-nghost)  =  undefined_val
              endif
           enddo 
           enddo 
        enddo 

      end subroutine get_ice_fr_thickness_tgrid

      subroutine get_snow_thickness_tgrid ( &
                                m_hsno,  undefined_val)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: tmask
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         m_hsno

      real (kind=real_kind), intent(in) :: &
         undefined_val
         
!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, k, n
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
              if(tmask(i,j,iblk)) then
                m_hsno(i-nghost,j-nghost)  =  vsno(i,j,iblk)
              else
                m_hsno(i-nghost,j-nghost)  =  undefined_val
              endif
           enddo 
           enddo 
        enddo 

      end subroutine get_snow_thickness_tgrid

!=======================================================================
!BOP
!
!
! !INTERFACE:
!
      subroutine get_ice_vel_tgrid_lon_lat ( &
                                m_uvel,  m_vvel, undef, setundef)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: to_tgrid, u2tgrid_vector, ANGLET, &
                                 tmask, umask
      use ice_work,        only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         m_uvel, m_vvel
      real (kind=real_kind), intent(in) :: &
         undef
      logical (kind=log_kind), intent(in) :: &
         setundef

!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, k, n
       type (block) :: &
         this_block           ! block information for current block

       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::&
         uvelT, vvelT
       real (kind=real_kind) :: &
          workx, worky 
          
        uvelT(:,:,:) = uvel(:,:,:)  
        vvelT(:,:,:) = vvel(:,:,:)  
        call u2tgrid_vector(uvelT)
        call u2tgrid_vector(vvelT)

        do iblk = 1, nblocks  
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
           do j = jlo, jhi
           do i = ilo, ihi
              workx = uvelT(i,j,iblk)
              worky = vvelT(i,j,iblk)
              uvelT(i,j,iblk)  =  workx*cos(ANGLET(i,j,iblk)) &
                                - worky*sin(ANGLET(i,j,iblk))  
              vvelT(i,j,iblk)  =  worky*cos(ANGLET(i,j,iblk)) &
                                + workx*sin(ANGLET(i,j,iblk))   
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
              if(tmask(i,j,iblk)) then
                     m_uvel(i-nghost,j-nghost) = uvelT(i,j,iblk)
                     m_vvel(i-nghost,j-nghost) = vvelT(i,j,iblk)
              else
                 if(setundef) then
                     m_uvel(i-nghost,j-nghost) = undef
                     m_vvel(i-nghost,j-nghost) = undef
                 else
                     m_uvel(i-nghost,j-nghost) = c0
                     m_vvel(i-nghost,j-nghost) = c0
                 endif
              endif
           enddo 
           enddo 
        enddo 

      end subroutine get_ice_vel_tgrid_lon_lat 

      subroutine get_ice_vol_transport(transix, transiy, rotatevec)

      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: tmask, u2tgrid_vector, ANGLET,  &
                           HTN, HTE, dxt, dyt
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         transix, transiy

      logical (kind=log_kind), intent(in) :: &
         rotatevec
!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block

       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::&
         uvelT, vvelT
       real (kind=dbl_kind) :: &
          workx, worky 

        !uvelT(:,:,:) = uvel(:,:,:)
        !vvelT(:,:,:) = vvel(:,:,:)
        !call u2tgrid_vector(uvelT)
        !call u2tgrid_vector(vvelT)

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
              if(tmask(i,j,iblk)) then
                workx = HTE(i,j,iblk)*p5*(uvel(i,j,iblk)+uvel(i,j-1,iblk))* &
                                  (p5*(vice(i,j,iblk)+vice(i+1,j,iblk)) + &
                                   p5*(vsno(i,j,iblk)+vsno(i+1,j,iblk)))
                worky = HTN(i,j,iblk)*p5*(vvel(i,j,iblk)+vvel(i-1,j,iblk))* &
                                  (p5*(vice(i,j,iblk)+vice(i,j+1,iblk)) + &
                                   p5*(vsno(i,j,iblk)+vsno(i,j+1,iblk)))
                !workx = dyt(i,j,iblk)*uvelT(i,j,iblk)* &
                !                      (vice(i,j,iblk) + & 
                !                       vsno(i,j,iblk))
                !worky = dxt(i,j,iblk)*vvelT(i,j,iblk)* &
                !                      (vice(i,j,iblk) + &
                !                       vsno(i,j,iblk))
                if(rotatevec) then  
                    ! rotate onto lat-lon as CMIP5 requires
                    transix(i1,j1)  =  workx*cos(ANGLET(i,j,iblk)) &
                                     - worky*sin(ANGLET(i,j,iblk))  
                    transiy(i1,j1)  =  worky*cos(ANGLET(i,j,iblk)) &
                                     + workx*sin(ANGLET(i,j,iblk))   
                else
                    transix(i1,j1)  =  workx
                    transiy(i1,j1)  =  worky
                endif
              endif
           enddo 
           enddo 
        enddo 

      end  subroutine get_ice_vol_transport

      subroutine get_ice_mass_transport(transix, transiy, rotatevec)

      use ice_boundary
      use ice_domain
      use ice_constants
      use ice_grid,        only: tmask, u2tgrid_vector, ANGLET,  &
                                 HTN, HTE, dxt, dyt
!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=real_kind), &
         dimension(:,:), intent(out) :: &
         transix, transiy

      logical (kind=log_kind), intent(in) :: &
         rotatevec
!
!EOP
!
       integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, i1, j1, k, n
       type (block) :: &
         this_block           ! block information for current block

       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::&
         uvelT, vvelT
       real (kind=dbl_kind) :: &
          workx, worky 

        !uvelT(:,:,:) = uvel(:,:,:)
        !vvelT(:,:,:) = vvel(:,:,:)
        !call u2tgrid_vector(uvelT)
        !call u2tgrid_vector(vvelT)

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
              if(tmask(i,j,iblk)) then
                workx  = HTE(i,j,iblk)*p5*(uvel(i,j,iblk)+uvel(i,j-1,iblk))* &
                                  (rhoi*p5*(vice(i,j,iblk)+vice(i+1,j,iblk)) + &
                                   rhos*p5*(vsno(i,j,iblk)+vsno(i+1,j,iblk)))
                worky  = HTN(i,j,iblk)*p5*(vvel(i,j,iblk)+vvel(i-1,j,iblk))* &
                                  (rhoi*p5*(vice(i,j,iblk)+vice(i,j+1,iblk)) + &
                                   rhos*p5*(vsno(i,j,iblk)+vsno(i,j+1,iblk)))
                !workx = dyt(i,j,iblk)*uvelT(i,j,iblk)* &
                !                      (rhoi*vice(i,j,iblk) + & 
                !                       rhos*vsno(i,j,iblk))
                !worky = dxt(i,j,iblk)*vvelT(i,j,iblk)* &
                !                      (rhoi*vice(i,j,iblk) + &
                !                       rhos*vsno(i,j,iblk))
                if(rotatevec) then  
                    ! rotate onto lat-lon as CMIP5 requires
                    transix(i1,j1)  =  workx*cos(ANGLET(i,j,iblk)) &
                                     - worky*sin(ANGLET(i,j,iblk))  
                    transiy(i1,j1)  =  worky*cos(ANGLET(i,j,iblk)) &
                                     + workx*sin(ANGLET(i,j,iblk))   
                else
                    transix(i1,j1)  =  workx
                    transiy(i1,j1)  =  worky
                endif
              endif
           enddo 
           enddo 
        enddo 

      end  subroutine get_ice_mass_transport

! =======================================================================

      subroutine init_trcr_depend(iage, pond)
       
       logical(kind=log_kind), intent(in) :: &
           iage, pond 

       trcr_depend(nt_Tsfc)  = 0   ! ice/snow surface temperature
       if (iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
       if (pond) trcr_depend(nt_volpn) = 0   ! melt pond volume

      end subroutine init_trcr_depend

#endif
      end module ice_state

!=======================================================================
