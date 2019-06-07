!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: setup_grid: setup regriding field
!
! !INTERFACE:
!
   subroutine setup_grid(im,jm,im2,jm2,dphi,dlam,lons,lats)

! !REVISION HISTORY:
!
!  01Feb2001 Cruz     Created
!
!EOP
!-------------------------------------------------------------------------

   implicit none
   integer :: i,j,im,jm,im2,jm2,loc
   real :: lats(im2*jm2), dphi(jm), dphi2, dlpin
   real :: lons(im2*jm2), dlam(im), dlam2, dlin, dpin
   real :: lon,lat
   real :: pi

   pi = 4.0*atan(1.0)
   dlin = 2*pi/im
   dpin = pi/(jm-1)
   dlam(:) = dlin
   dphi(:) = dpin

   dlam2   = 2.*pi/ im2
   dphi2   = pi/(jm2-1)

   loc=0
   do j=1,jm2
     do i=1,im2
       loc = loc + 1
       lon = -pi + (i-1)*dlam2
       lons(loc) = lon
     enddo
   enddo
   loc = 0
   do j=1,jm2
     lat = -pi/2.0 + (j-1)*dphi2
     do i=1,im2
       loc = loc + 1
       lats(loc) = lat
     enddo
   enddo

   end subroutine

!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: comput_gaus: compute guassian latitudes and
!           create output lons and lats
!
! !INTERFACE:
!
   subroutine comput_gaus(im,jm,im2,jm2,dphi,dlam,lons,lats,flip)

! !REVISION HISTORY:
!
!  01Feb2001 Cruz     Created
!  30Jul2004 Todling  Added flip option
!
!EOP
!-------------------------------------------------------------------------

   implicit none
   integer :: i,j,im,jm,im2,jm2,loc
   real :: lats(im2*jm2), dphi(jm)
   real :: lons(im2*jm2), dlam(im)
   logical :: flip
   integer :: lon,lat,n,m
   real :: dl,dp,pi
   real :: glat(jm2) ! automatic

   pi = 4.0*atan(1.0)
   dl = 2.0*pi/im
   dp = pi/( jm-1 )
   dlam(:) = dl  ! Uniform grid input
   dphi(:) = dp  ! Uniform grid input

   call gauss_lat_nmc (glat,jm2)
   dl = 2.0*pi/im2
   if (flip) glat(:) = glat(jm2:1:-1)

   loc = 0
   do j=1,jm2
     do i=1,im2
       loc = loc + 1
       lons (loc) = -pi + (i-1)*dl
     enddo
   enddo

   loc = 0
   do j=1,jm2
     do i=1,im2
       loc = loc + 1
       lats (loc) = glat(j)*pi/180.0
     enddo
   enddo

   end subroutine


