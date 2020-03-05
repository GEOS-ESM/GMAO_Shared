
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



   subroutine DustEmissionGOCART2G( i1, i2, j1, j2, n_bins, radius, &
                                  fraclake, gwettop, oro, u10m, v10m, &
                                  Ch_DU, sfrac, du_src, grav, &
                                  emissions, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, n_bins           ! grid dimensions
   real, dimension(:) :: radius, sfrac                      ! particle radius [m]
   real, pointer, dimension(:,:) :: fraclake, gwettop, oro, u10m, v10m, du_src
   real :: Ch_DU, grav

! !OUTPUT PARAMETERS:
!   real, intent(inout) ::  emissions(i1:i2, j1:j2, n_bins)    ! Local emission
   real  ::  emissions(i1:i2, j1:j2, n_bins)    ! Local emission

   integer, intent(out) :: rc                                 ! Error return code:


! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 11Feb202 E.Sherman - First attempt at refactor
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, e
   real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density  = 2650.  ! km m-3
   real            ::  diameter         ! dust effective diameter [m]
   real            ::  u_thresh0
   real            ::  u_thresh
   real            ::  w10m

!   _UNUSED_DUMMY(km)

!print*,'DustEmissionsGOCART2G BEGIN'

!  Initialize local variables
!  --------------------------
   emissions(:,:,:) = 0.

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed 
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.


!print*,'DustEmissionsGOCART2G i1 = ', i1
!print*,'DustEmissionsGOCART2G i2 = ', i2
!print*,'DustEmissionsGOCART2G j1 = ', j1
!print*,'DustEmissionsGOCART2G j2 = ', j2
!print*,'DustEmissionsGOCART2G shape(gwettop) = ',shape(gwettop)
!print*,'DustEmissionsGOCART2G shape(u10m)    = ',shape(u10m)
!print*,'DustEmissionsGOCART2G shape(v10m)    = ',shape(v10m)

!print*,'DustEmissionsGOCART2G Ch_DU  = ', Ch_DU
!print*,'DustEmissionsGOCART2G radius  = ', radius
!print*,'DustEmissionsGOCART2G n_bins  = ', n_bins


  do e = 1, n_bins
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
    emissions(:,:,e) = emissions(:,:,e) * Ch_DU * sfrac(e) * du_src
 end do ! e

!print*,'DustEmissionsGOCART2G diameter = ', diameter
!print*,'DustEmissionsGOCART2G u_thresh0 = ', u_thresh0
!print*,'DustEmissionsGOCART2G w10m = ', w10m

   rc=0
!print*,'DustEmissionsGOCART2G sum(emissions) = ',sum(emissions)
!print*,'DustEmissionsGOCART2G END'

   end subroutine DustEmissionGOCART2G




   end module GriddedEmission
