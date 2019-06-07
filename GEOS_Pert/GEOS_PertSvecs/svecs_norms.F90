module svecs_norms

#include "MAPL_Generic.h"

use ESMF
use MAPL_mod
use svecs_cf

implicit none
private

real(8), parameter :: Tref = 300.0_8
real(8), parameter :: Pref = 10000.0_8

public svec_norm_a, svec_norm_bc1, svec_norm_bc2, svec_norm_d, &
       svec_norm_a_inv

contains

subroutine svec_norm_a(xpert)

 !Take in z from Lanczos, apply norm and put in xpert
 !xpert = S^-1/2 z

 implicit none

 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 real(8) :: Tfactor, pfactor, qfactor

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%u(i,j,k) = xpert%u(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%v(i,j,k) = xpert%v(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    Tfactor = MAPL_CP / Tref

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpI(i,j,k) > zero) then
                xpert%t(i,j,k) = xpert%t(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*Tfactor)
             endif

          enddo
       enddo
    enddo

    pfactor = MAPL_RDRY * Tref / Pref**2

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%delp(i,j,k) = xpert%delp(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*pfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%t = 0.0_8
    xpert%delp = 0.0_8

 endif

 if (svecnormI == 'we') then

    qfactor = eps_eer * MAPL_ALHL**2 / (MAPL_CP * Tref)

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%sphu(i,j,k) = xpert%sphu(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*qfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%sphu = 0.0_8

 endif

 !Not included in the norm (superfluous):
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine svec_norm_a

subroutine svec_norm_bc1(xpert)

 !First part of final time norm
 !x = EF^1/2 x

 implicit none

 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 real(8) :: Tfactor, pfactor, qfactor

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%u(i,j,k) = xpert%u(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k))
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%v(i,j,k) = xpert%v(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k))
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    Tfactor = MAPL_CP / Tref

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpF(i,j,k) > zero) then
                xpert%t(i,j,k) = xpert%t(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*Tfactor)
             endif

          enddo
       enddo
    enddo

    pfactor = MAPL_RDRY * Tref / Pref**2

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%delp(i,j,k) = xpert%delp(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*pfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%t = 0.0_8
    xpert%delp = 0.0_8

 endif

 if (svecnormI == 'we') then

    qfactor = eps_eer * MAPL_ALHL**2 / (MAPL_CP * Tref)

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%sphu(i,j,k) = xpert%sphu(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*qfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%sphu = 0.0_8

 endif

 !Not included in the norm (superfluous):
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine svec_norm_bc1


subroutine svec_norm_bc2(xpert)

 !Second part of final time norm
 !x = EF^T/2 x

 implicit none

 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 real(8) :: Tfactor, pfactor, qfactor

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%u(i,j,k) = xpert%u(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k))
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%v(i,j,k) = xpert%v(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k))
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    Tfactor = MAPL_CP / Tref

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpF(i,j,k) > zero) then
                xpert%t(i,j,k) = xpert%t(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*Tfactor)
             endif

          enddo
       enddo
    enddo

    pfactor = MAPL_RDRY * Tref / Pref**2

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%delp(i,j,k) = xpert%delp(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*pfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%t = 0.0_8
    xpert%delp = 0.0_8

 endif

 if (svecnormI == 'we') then

    qfactor = eps_eer * MAPL_ALHL**2 / (MAPL_CP * Tref)

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpF(i,j,k) > zero) then
                xpert%sphu(i,j,k) = xpert%sphu(i,j,k) * sqrt(0.5_8*gridweightF(i,j,k)*qfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%sphu = 0.0_8

 endif

 !Not included in the norm (superfluous):
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine svec_norm_bc2


subroutine svec_norm_d(xpert)

 !Take in xpert from adjoint apply norm and convert to z for Lanczos
 !z = S^-T/2 xpert

 implicit none

 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 real(8) :: Tfactor, pfactor, qfactor

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%u(i,j,k) = xpert%u(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%v(i,j,k) = xpert%v(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    Tfactor = MAPL_CP / Tref

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpI(i,j,k) > zero) then
                xpert%t(i,j,k) = xpert%t(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*Tfactor)
             endif

          enddo
       enddo
    enddo

    pfactor = MAPL_RDRY * Tref / Pref**2

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%delp(i,j,k) = xpert%delp(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*pfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%t = 0.0_8
    xpert%delp = 0.0_8

 endif

 if (svecnormI == 'we') then

    qfactor = eps_eer * MAPL_ALHL**2 / (MAPL_CP * Tref)

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%sphu(i,j,k) = xpert%sphu(i,j,k) / sqrt(0.5_8*gridweightI(i,j,k)*qfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%sphu = 0.0_8

 endif

 !Not included in the norm (superfluous):
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine svec_norm_d


subroutine svec_norm_a_inv(xpert)

 !Take in z from Lanczos, apply norm and put in xpert
 !xpert = S^-1/2 z

 implicit none

 type(pertState), intent(inout) :: xpert

 integer :: i,j,k
 real(8) :: Tfactor, pfactor, qfactor

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%u(i,j,k) = xpert%u(i,j,k) * sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%v(i,j,k) = xpert%v(i,j,k) * sqrt(0.5_8*gridweightI(i,j,k))
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    Tfactor = MAPL_CP / Tref

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpI(i,j,k) > zero) then
                xpert%t(i,j,k) = xpert%t(i,j,k) * sqrt(0.5_8*gridweightI(i,j,k)*Tfactor)
             endif

          enddo
       enddo
    enddo

    pfactor = MAPL_RDRY * Tref / Pref**2

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%delp(i,j,k) = xpert%delp(i,j,k) * sqrt(0.5_8*gridweightI(i,j,k)*pfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%t = 0.0_8
    xpert%delp = 0.0_8

 endif

 if (svecnormI == 'we') then

    qfactor = eps_eer * MAPL_ALHL**2 / (MAPL_CP * Tref)

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                xpert%sphu(i,j,k) = xpert%sphu(i,j,k) * sqrt(0.5_8*gridweightI(i,j,k)*qfactor)
             endif

          enddo
       enddo
    enddo

 else

    xpert%sphu = 0.0_8

 endif

 !Not included in the norm (superfluous):
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine svec_norm_a_inv

end module svecs_norms
