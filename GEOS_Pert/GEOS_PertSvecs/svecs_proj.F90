module svecs_proj

use ESMF
use MAPL_mod
use svecs_cf

implicit none
private

#include "mpif.h"
integer :: ierror

public proj_init, proj_start_z2xpert, proj_start_xpert2z, proj_final_xpert2xpert

contains

subroutine proj_init(SVFILE)

 implicit none

 character(len=ESMF_MAXSTR) :: SVFILE
 type (ESMF_Config)         :: cf
 integer                    :: rc

 integer                    :: i,j,l
 logical                    :: invalid

 real(8), allocatable       :: maxlons(:)
 real(8), allocatable       :: minlons(:)
 real(8), allocatable       :: maxlats(:)
 real(8), allocatable       :: minlats(:)
 real(8) :: maxlon, minlon, maxlat, minlat

 cf = ESMF_ConfigCreate(rc=rc)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = rc )

 allocate(maxlons(im))
 allocate(minlons(im))
 allocate(maxlats(jm))
 allocate(minlats(jm))

 !Initial time projection
 !-----------------------

 invalid = .false.

 projlonI(1) =  -180.0_8
 projlonI(2) =   180.0_8
 projlatI(1) =  -90.0_8
 projlatI(2) =   90.0_8
 projlevI(1) =    1
 projlevI(2) =   lm

 !west
 call ESMF_ConfigGetAttribute( cf, projlonI(1), label='start_svec_west:' , default=projlonI(1), rc = rc)

 !east
 call ESMF_ConfigGetAttribute( cf, projlonI(2), label='start_svec_east:' , default=projlonI(2), rc = rc)

 !south
 call ESMF_ConfigGetAttribute( cf, projlatI(1), label='start_svec_south:', default=projlatI(1), rc = rc)

 !north
 call ESMF_ConfigGetAttribute( cf, projlatI(2), label='start_svec_north:', default=projlatI(2), rc = rc)

 !upper
 call ESMF_ConfigGetAttribute( cf, projlevI(1), label='start_svec_upper:', default=projlevI(1), rc = rc)

 !lower
 call ESMF_ConfigGetAttribute( cf, projlevI(2), label='start_svec_lower:', default=projlevI(2), rc = rc)

 !Quick sanity check
 if ( projlatI(1) > projlatI(2) ) invalid = .true.
 if ( projlevI(1) > projlevI(2) ) invalid = .true.

 if ( invalid ) then
    projlonI(1) =  -180.0_8
    projlonI(2) =  180.0_8
    projlatI(1) =  -90._8
    projlatI(2) =   90._8
    projlevI(1) =    1
    projlevI(2) =   99
    if (MAPL_AM_I_ROOT()) &
    print*, 'Something went wrong while setting projection box. Taking default local projection (-180,180) (-90,90) (All Levels)'
 endif

 if (.not.invalid) then
    maxlons = -9999.0_8
    minlons =  9999.0_8
    maxlats = -9999.0_8
    minlats =  9999.0_8
    coefftpI = zero
    do i = 1,im
       do j = 1,jm
          do l = 1,lm
   
             if ( lons(i) >= projlonI(1) .and. lons(i) <= projlonI(2) .and. &
                  lats(j) >= projlatI(1) .and. lats(j) <= projlatI(2) .and. &
                  levs(l) >= projlevI(1) .and. levs(l) <= projlevI(2) ) then
   
                maxlons(i) = lons(i)
                minlons(i) = lons(i)
                maxlats(j) = lats(j)
                minlats(j) = lats(j)
                coefftpI(i,j,l) = one
   
             endif

          enddo
       enddo
    enddo
    maxlon = maxval(maxlons)
    minlon = minval(minlons)
    maxlat = maxval(maxlats)
    minlat = minval(minlats)
 else
    coefftpI = one
    maxlon = maxval(lons)
    minlon = minval(lons)
    maxlat = maxval(lats)
    minlat = minval(lats)
 endif

 call mpi_allreduce_max(maxlon)
 call mpi_allreduce_min(minlon)
 call mpi_allreduce_max(maxlat)
 call mpi_allreduce_min(minlat)

 projlongridI(1) = minlon
 projlongridI(2) = maxlon
 projlatgridI(1) = minlat
 projlatgridI(2) = maxlat

 if (MAPL_AM_I_ROOT()) then
    print*, 'User specified initial time projection: '
    print*, 'From Lon ', real(projlonI(1),4), ' to Lon ', real(projlonI(2),4)
    print*, 'From Lat ', real(projlatI(1),4), ' to Lat ', real(projlatI(2),4)
    print*, 'From Lev ', projlevI(1), ' to Lev ', projlevI(2)
    print*, ''
    print*, 'Grid lon: from ', real(projlongridI(1),4), ' to ', real(projlongridI(2),4)
    print*, 'Grid lat: from ', real(projlatgridI(1),4), ' to ', real(projlatgridI(2),4)
 endif


 !Final time projection
 !---------------------

 invalid = .false.

 projlonF(1) =  -180.0_8
 projlonF(2) =   180.0_8
 projlatF(1) =  -90.0_8
 projlatF(2) =   90.0_8
 projlevF(1) =    1
 projlevF(2) =   lm

 !west
 call ESMF_ConfigGetAttribute( cf, projlonF(1), label='final_svec_west:' , default=projlonF(1), rc = rc)

 !east
 call ESMF_ConfigGetAttribute( cf, projlonF(2), label='final_svec_east:' , default=projlonF(2), rc = rc)

 !south
 call ESMF_ConfigGetAttribute( cf, projlatF(1), label='final_svec_south:', default=projlatF(1), rc = rc)

 !north
 call ESMF_ConfigGetAttribute( cf, projlatF(2), label='final_svec_north:', default=projlatF(2), rc = rc)

 !upper
 call ESMF_ConfigGetAttribute( cf, projlevF(1), label='final_svec_upper:', default=projlevF(1), rc = rc)

 !lower
 call ESMF_ConfigGetAttribute( cf, projlevF(2), label='final_svec_lower:', default=projlevF(2), rc = rc)

 !Quick sanity check
 if ( projlatF(1) > projlatF(2) ) invalid = .true.
 if ( projlevF(1) > projlevF(2) ) invalid = .true.

 if ( invalid ) then
    projlonF(1) =  -180.0_8
    projlonF(2) =  180.0_8
    projlatF(1) =  -90._8
    projlatF(2) =   90._8
    projlevF(1) =    1
    projlevF(2) =   99
    if (MAPL_AM_I_ROOT()) &
    print*, 'Something went wrong while setting projection box. Taking default local projection (-180,180) (-90,90) (All Levels)'
 endif

 if (.not.invalid) then
    maxlons = -9999.0_8
    minlons =  9999.0_8
    maxlats = -9999.0_8
    minlats =  9999.0_8
    coefftpF = zero
    do i = 1,im
       do j = 1,jm
          do l = 1,lm
   
             if ( lons(i) >= projlonF(1) .and. lons(i) <= projlonF(2) .and. &
                  lats(j) >= projlatF(1) .and. lats(j) <= projlatF(2) .and. &
                  levs(l) >= projlevF(1) .and. levs(l) <= projlevF(2) ) then
   
                maxlons(i) = lons(i)
                minlons(i) = lons(i)
                maxlats(j) = lats(j)
                minlats(j) = lats(j)
                coefftpF(i,j,l) = one
   
             endif

          enddo
       enddo
    enddo
    maxlon = maxval(maxlons)
    minlon = minval(minlons)
    maxlat = maxval(maxlats)
    minlat = minval(minlats)
 else
    maxlon = maxval(lons)
    minlon = minval(lons)
    maxlat = maxval(lats)
    minlat = minval(lats)
    coefftpF = one
 endif

 call mpi_allreduce_max(maxlon)
 call mpi_allreduce_min(minlon)
 call mpi_allreduce_max(maxlat)
 call mpi_allreduce_min(minlat)

 projlongridF(1) = minlon
 projlongridF(2) = maxlon
 projlatgridF(1) = minlat
 projlatgridF(2) = maxlat

 if (MAPL_AM_I_ROOT()) then
    print*, ''
    print*, 'User specified final time projection: '
    print*, 'From Lon ', real(projlonF(1),4), ' to Lon ', real(projlonF(2),4)
    print*, 'From Lat ', real(projlatF(1),4), ' to Lat ', real(projlatF(2),4)
    print*, 'From Lev ', projlevF(1), ' to Lev ', projlevF(2)
    print*, ''
    print*, 'Grid lon: from ', real(projlongridF(1),4), ' to ', real(projlongridF(2),4)
    print*, 'Grid lat: from ', real(projlatgridF(1),4), ' to ', real(projlatgridF(2),4)
    print*, ''
 endif

 deallocate(maxlons,minlons,minlats,maxlats)

end subroutine proj_init

subroutine proj_start_z2xpert(z,xpert,nz)

!Project z into xpert at initial time

 implicit none

 integer, intent(in) :: nz
 real(8), intent(in) :: z(nz)
 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 integer :: zcount

 zcount = 0

 !Set all variables to zero
 xpert%u = 0.0_8
 xpert%v = 0.0_8
 xpert%t = 0.0_8
 xpert%delp = 0.0_8
 xpert%sphu = 0.0_8
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > 0.0) then
                zcount = zcount + 1
                xpert%u(i,j,k) = z(zcount)
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > 0.0) then
                zcount = zcount + 1
                xpert%v(i,j,k) = z(zcount)
             endif

          enddo
       enddo
    enddo

    if (svecnormI == 'te' .or. svecnormI == 'we') then

       !t dry temperature
       do k = 1,lm
          do j = 1,jm
             do i = 1,im

                if (coefftpI(i,j,k) > 0.0) then
                   zcount = zcount + 1
                   xpert%t(i,j,k) = z(zcount)
                endif

             enddo
          enddo
       enddo

       !delp pressure thickness
       do k = 1,lm
          do j = 1,jm
             do i = 1,im
   
                if (coefftpI(i,j,k) > 0.0) then
                   zcount = zcount + 1
                   xpert%delp(i,j,k) = z(zcount)
                endif

             enddo
          enddo
       enddo

       if (svecnormI == 'we') then

          !sphu specific humidity
          do k = 1,lm
             do j = 1,jm
                do i = 1,im

                   if (coefftpI(i,j,k) > 0.0) then
                      zcount = zcount + 1
                      xpert%sphu(i,j,k) = z(zcount)
                   endif

                enddo
             enddo
          enddo

       endif

    endif

 endif

end subroutine proj_start_z2xpert

subroutine proj_start_xpert2z(z,xpert,nz)

!Project xpert into z at initial time

 implicit none

 integer, intent(in) :: nz
 real(8), intent(inout) :: z(nz)
 type(pertState), intent(in) :: xpert

 integer :: i,j,k
 integer :: zcount

 z = 0.0_8 !reset Lanczos vector
 zcount = 0

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > 0.0) then
                zcount = zcount + 1
                z(zcount) = xpert%u(i,j,k)
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > 0.0) then
                zcount = zcount + 1
                z(zcount) = xpert%v(i,j,k)
             endif

          enddo
       enddo
    enddo

    if (svecnormI == 'te' .or. svecnormI == 'we') then

       !t dry temperature
       do k = 1,lm
          do j = 1,jm
             do i = 1,im

                if (coefftpI(i,j,k) > 0.0) then
                   zcount = zcount + 1
                   z(zcount) = xpert%t(i,j,k)
                endif

             enddo
          enddo
       enddo

       !delp pressure thickness
       do k = 1,lm
          do j = 1,jm
             do i = 1,im
   
                if (coefftpI(i,j,k) > 0.0) then
                   zcount = zcount + 1
                   z(zcount) = xpert%delp(i,j,k)
                endif

             enddo
          enddo
       enddo

       if (svecnormI == 'we') then

          !sphu specific humidity
          do k = 1,lm
             do j = 1,jm
                do i = 1,im

                   if (coefftpI(i,j,k) > 0.0) then
                      zcount = zcount + 1
                      z(zcount) = xpert%sphu(i,j,k)
                   endif

                enddo
             enddo
          enddo

       endif

    endif

 endif


end subroutine proj_start_xpert2z

subroutine proj_final_xpert2xpert(xpert)

 !Screen xpert using the projection

 implicit none

 type(pertState), intent(inout) :: xpert

 !u winds
 xpert%u = xpert%u * coefftpF

 !v winds
 xpert%v = xpert%v * coefftpF

 !t temperature
 xpert%t = xpert%t * coefftpF

 !delp pressure thickness
 xpert%delp = xpert%delp * coefftpF

 !sphu specific humidity
 xpert%sphu = xpert%sphu * coefftpF

 !qitot cloud liquid ice
 xpert%qitot = xpert%qitot * coefftpF

 !qltot cloud liquid water
 xpert%qltot = xpert%qltot * coefftpF

 !ozone
 xpert%ozone = xpert%ozone * coefftpF

end subroutine proj_final_xpert2xpert


subroutine mpi_allreduce_max(mymax)

 real(8), intent(INOUT)  :: mymax
 real(8) :: gmax

 call MPI_ALLREDUCE( mymax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                     MPI_COMM_WORLD, ierror )

 mymax = gmax

end subroutine mpi_allreduce_max

subroutine mpi_allreduce_min(mymin)

 real(8), intent(inout)  :: mymin
 real(8)                 :: gmin

 call MPI_ALLREDUCE( mymin, gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                     MPI_COMM_WORLD, ierror )

 mymin = gmin

end subroutine mpi_allreduce_min

end module svecs_proj
