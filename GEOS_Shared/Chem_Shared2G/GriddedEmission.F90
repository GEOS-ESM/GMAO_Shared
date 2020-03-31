
!-------------------------------------------------------------------------
!
! !MODULE:  GriddedEmission.F90 --- Calculate the gridded emissions
!
! !INTERFACE:
!

   module  GriddedEmission

! !USES:
!  Uses nothing. Only instrinsic types and functions are allowed.
!   use ESMF
!   use MAPL

   implicit none

! !PUBLIC TYPES:
!
   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
   public DustEmissionGOCART2G
   public DistributePointEmission

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

!
! !DESCRIPTION:
!
!  This module implements the gridded emission calculations
!
! !REVISION HISTORY:
!
!  11Feb2020  E.Sherman, A.da Silva, T.Clune, A.Darmenov -  First attempt at refactor
!
!EOP
!-------------------------------------------------------------------------
CONTAINS


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DustEmissionGOCART2G

   subroutine DustEmissionGOCART2G(radius, fraclake, gwettop, oro, u10m, &
                                   v10m, Ch_DU, sfrac, du_src, grav, &
                                   emissions, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
!   integer, intent(in) :: i1, i2, j1, j2, nbins            ! grid dimensions
!   real, dimension(:) :: radius, sfrac                      ! particle radius [m]
!   real, pointer, dimension(:,:) :: fraclake, gwettop, oro, u10m, v10m, du_src
!   real :: Ch_DU, grav

   real, dimension(:), intent(in) :: radius, sfrac                      ! particle radius [m]
   real, pointer, dimension(:,:), intent(in) :: fraclake, gwettop, oro, u10m, v10m, du_src
   real, intent(in) :: Ch_DU, grav

! !OUTPUT PARAMETERS:
!   real  ::  emissions(i1:i2, j1:j2, nbins)    ! Local emission
   real, intent(out)  ::  emissions(:,:,:)    ! Local emission
!   real, pointer, intent(inout)  ::  emissions(:,:,:)    ! Local emission

   integer, intent(out) :: rc                   ! Error return code:


! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at refactor
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer         ::  i, j, e
   real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density  = 2650.  ! km m-3
   real            ::  diameter         ! dust effective diameter [m]
   real            ::  u_thresh0
   real            ::  u_thresh
   real            ::  w10m
   integer         ::  i1, i2, j1, j2, nbins
   integer         ::  dims(2)
!   _UNUSED_DUMMY(km)

!  Initialize local variables
!  --------------------------
   emissions(:,:,:) = 0.
   rc = 824

!  Get dimensions
!  ---------------
   nbins = size(radius)
   dims = shape(u10m)
   i1 = 1; j1 = 1
   i2 = dims(1); j2 = dims(2)

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed 
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.


  do e = 1, nbins
   diameter = 2. * radius(e)

   u_thresh0 = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
           / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)

!  Spatially dependent part of calculation
!  ---------------------------------------
   do j = j1, j2
    do i = i1, i2

     if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

     w10m = sqrt(u10m(i,j)**2.+v10m(i,j)**2.)

!    Modify the threshold depending on soil moisture as in Ginoux et al. [2001]
     if(gwettop(i,j) .lt. 0.5) then
      u_thresh = amax1(0.,u_thresh0* &
       (1.2+0.2*alog10(max(1.e-3,gwettop(i,j)))))

       if(w10m .gt. u_thresh) then     
!       Emission of dust [kg m-2 s-1]
       emissions(i,j,e) = (1.-fraclake(i,j)) * w10m**2. * (w10m-u_thresh)

       endif
      endif

     end do   ! i
    end do    ! j
    emissions(:,:,e) =  Ch_DU * sfrac(e) * du_src *  emissions(:,:,e)
 end do ! e

   rc=0

   end subroutine DustEmissionGOCART2G

!------------------------------------------------------------------------
!BOP
! !IROUTINE: DistributePointEmissions

! !INTERFACE:

   subroutine DistributePointEmission(km, delp, rhoa, z_bot, z_top, &
                                      GRAV, emissions_point, &
                                      point_column_emissions, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)             :: km
   real, dimension(:), intent(in)  :: delp, rhoa
   real, intent(in)                :: z_bot, z_top
   real,               intent(in)  :: GRAV, emissions_point


! !OUTPUT PARAMETERS:
   real, dimension(:), intent(out) ::  point_column_emissions

   integer, intent(out) :: rc                                 ! Error return code:


! !DESCRIPTION: Distributes piont emissions
!
! !REVISION HISTORY:
! ??? P. Colarco
! 16March2020 E.Sherman - Moved from original DU_GridCompMod.F90 and generalized
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
    integer :: k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_



!   find level height
    z = 0.0
    z_= 0.0

    do k = km, 1, -1
       dz(k) = delp(k)/rhoa(k)/grav
       z_    = z_ + dz(k)
       z(k)  = z_
    end do

!   find the bottom level
    do k = km, 1, -1
       if (z(k) >= z_bot) then
           k_bot = k
           exit
       end if
    end do

!   find the top level
    do k = k_bot, 1, -1
       if (z(k) >= z_top) then
           k_top = k
           exit
       end if
    end do

!   find the weights
    w_ = 0


!   if (k_top > k_bot) then
!       need to bail - something went wrong here
!   end if

    if (k_bot .eq. k_top) then
        w_(k_bot) = z_top - z_bot
    else
     do k = k_bot, k_top, -1
        if ((k < k_bot) .and. (k > k_top)) then
             w_(k) = dz(k)
        else
             if (k == k_bot) then
                 w_(k) = (z(k) - z_bot)
             end if

             if (k == k_top) then
                 w_(k) = z_top - (z(k)-dz(k))
             end if
        end if
     end do
    end if

!   distribute emissions in the vertical
    point_column_emissions(:) = (w_ / sum(w_)) * emissions_point

    rc = 0

    end subroutine DistributePointEmission





 end module GriddedEmission
