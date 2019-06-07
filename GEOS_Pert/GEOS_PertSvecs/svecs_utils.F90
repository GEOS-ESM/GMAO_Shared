module svecs_utils

#include "MAPL_Generic.h"

use ESMF
use MAPL_mod
use svecs_cf

implicit none
private

#include "mpif.h"

public numvars, norm_zvecsize, xpert_to_state, state_to_xpert, &
       t2tv, tv2t, t2tv_ad, tv2t_ad, dot_prod_xx, dot_prod_zz!, &
!       replicate_pole_points

! !DESCRIPTION: Config for singular vector calculation.
!
! !REVISION HISTORY:
!
!  28Apr2017  Holdaway  Initialized.


contains

subroutine numvars( n, m, im, jm, lm )

 implicit none

 integer, intent(in) :: im,jm,lm
 integer, intent(inout) :: n, m

 n = 0

 !u
 n = n + im*jm*lm
 !v
 n = n + im*jm*lm
 !pt
 n = n + im*jm*lm
 !delp
 n = n + im*jm*lm
 !q, qi, ql, o3
 n = n + 4*im*jm*lm

 m = n

end subroutine numvars


subroutine norm_zvecsize(nzvecsize)

 integer, intent(inout) :: nzvecsize

 integer :: i, j, k

 nzvecsize = 0

 if (svecnormI == 'ke' .or. svecnormI == 'te' .or. svecnormI == 'we') then

    !u winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                nzvecsize = nzvecsize + 1
             endif

          enddo
       enddo
    enddo

    !v winds
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                nzvecsize = nzvecsize + 1
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'te' .or. svecnormI == 'we') then

    !t dry temperature
    do k = 1,lm
       do j = 1,jm
          do i = 1,im
 
             if (coefftpI(i,j,k) > zero) then
                nzvecsize = nzvecsize + 1
             endif

          enddo
       enddo
    enddo

    !delp pressure thickness
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                nzvecsize = nzvecsize + 1
             endif

          enddo
       enddo
    enddo

 endif

 if (svecnormI == 'we') then

    !sphu specific humidity
    do k = 1,lm
       do j = 1,jm
          do i = 1,im

             if (coefftpI(i,j,k) > zero) then
                nzvecsize = nzvecsize + 1
             endif

          enddo
       enddo
    enddo

 endif

end subroutine norm_zvecsize

subroutine xpert_to_state(xpert,ESMFSTATE,ROOT)

 implicit none

 integer, intent(in) :: ROOT
 type(pertState), intent(inout) :: xpert
 type(ESMF_State), pointer, intent(inout) :: ESMFSTATE(:)

 real, pointer, dimension(:,:,:) :: fptr

 character(len=ESMF_MAXSTR) :: Iam="xpert_to_state"

 integer :: status, rc

 !u wind
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'U' , alloc=.true., rc = STATUS)
 VERIFY_(STATUS)
 fptr = real(xpert%u,4)
 nullify(fptr)

 !v wind
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'V', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%v,4)
 nullify(fptr)

 !t temperature
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'TV', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%t,4)
 nullify(fptr)

 !delp pressure
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'DP', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%delp,4)
 nullify(fptr)

 !q specifc humidity
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QV', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%sphu,4)
 nullify(fptr)

 !qitot cloud liquid ice
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QI', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%qitot,4)
 nullify(fptr)

 !qltot cloud liquid water
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QL', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%qltot,4)
 nullify(fptr)

 !ozone
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'O3', alloc=.true., rc = status)
 VERIFY_(status)
 fptr = real(xpert%ozone,4)
 nullify(fptr)

end subroutine xpert_to_state

subroutine state_to_xpert(ESMFSTATE,xpert,ROOT)

 implicit none

 integer, intent(in) :: ROOT
 type(pertState), intent(inout) :: xpert
 type(ESMF_State), pointer, intent(in) :: ESMFSTATE(:)

 real, pointer, dimension(:,:,:) :: fptr

 character(len=ESMF_MAXSTR)   :: Iam="state_to_xpert"

 integer :: status, rc

 !u wind
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'U', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%u = 0.0_8
 xpert%u = dble(fptr)
 nullify(fptr)

 !v wind
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'V', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%v = 0.0
 xpert%v = dble(fptr)
 nullify(fptr)

 !t temperature
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'TV', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%t = 0.0_8
 xpert%t = dble(fptr)
 nullify(fptr)

 !delp pressure
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'DP', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%delp = 0.0_8
 xpert%delp = dble(fptr)
 nullify(fptr)

 !q specifc humidity
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QV', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%sphu = 0.0_8
 xpert%sphu = dble(fptr)
 nullify(fptr)

 !qi cloud liquid ice
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QI', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%qitot = 0.0_8
 xpert%qitot = dble(fptr)
 nullify(fptr)

 !ql cloud liquid water
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'QL', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%qltot = 0.0_8
 xpert%qltot = dble(fptr)
 nullify(fptr)

 !ozone
 call MAPL_GetPointer(ESMFSTATE(ROOT), fptr, 'O3', alloc=.true., rc = status)
 VERIFY_(status)
 xpert%ozone = 0.0_8
 xpert%ozone = dble(fptr)
 nullify(fptr)

end subroutine state_to_xpert


subroutine t2tv(xpert,xtraj)

 implicit none

 type(pertstate), intent(inout) :: xpert
 type(trajpointer), intent(in) :: xtraj

 !Dry
 xpert%t = xpert%t * (one + dble(MAPL_VIREPS) * dble(xtraj%sphu))

 !Moist dependent
 !xpert%t = xpert%t*(one + dble(MAPL_VIREPS)*dble(xtraj%sphu)) + dble(xtraj%t)*dble(MAPL_VIREPS)*xpert%sphu

end subroutine t2tv

subroutine tv2t(xpert,xtraj)

 implicit none

 type(pertstate), intent(inout) :: xpert
 type(trajpointer), intent(in) :: xtraj

 !Dry
 xpert%t = xpert%t / (one + dble(MAPL_VIREPS) * xtraj%sphu)

 !Moist dependent
 !xpert%t = (xpert%t*(one + dble(MAPL_VIREPS) * dble(xtraj%sphu)) - dble(xtraj%t)*dble(MAPL_VIREPS)*xpert%sphu) &
  !                / (one + dble(MAPL_VIREPS) * dble(xtraj%sphu))**2

end subroutine tv2t

subroutine t2tv_ad(xpert,xtraj)

 implicit none

 type(pertstate), intent(inout) :: xpert
 type(trajpointer), intent(in) :: xtraj

 !Dry
 xpert%t = xpert%t * (one + dble(MAPL_VIREPS) * dble(xtraj%sphu))

 !Moist dependent
 !xpert%sphu = xpert%sphu + dble(xtraj%t)*dble(MAPL_VIREPS) * xpert%t
 !xpert%t = (one + dble(MAPL_VIREPS)*dble(xtraj%sphu)) * xpert%t

end subroutine t2tv_ad

subroutine tv2t_ad(xpert,xtraj)

 implicit none

 type(pertstate), intent(inout) :: xpert
 type(trajpointer), intent(in) :: xtraj

 !Dry
 xpert%t = xpert%t / (one + dble(MAPL_VIREPS) * dble(xtraj%sphu))

 !Moist dependent
 !xpert%sphu = xpert%sphu - (dble(xtraj%t)*dble(MAPL_VIREPS)/(1 + dble(MAPL_VIREPS)*dble(xtraj%sphu))**2 ) * xpert%t
 !xpert%t = ( one/(one + dble(MAPL_VIREPS)*dble(xtraj%sphu))) * xpert%t

end subroutine tv2t_ad


subroutine dot_prod_xx(xpert1,xpert2,dotp)

 implicit none

 type(pertState), intent(in) :: xpert1,xpert2
 real(8), intent(out) :: dotp
 
 type(ESMF_VM) :: vm
 integer :: i, j, k, rc, status
 real(8) :: dot_proc(1)

 dot_proc = zero

 do k = 1,npz
    do j = 1,jm
       do i = 1,im
          dot_proc = dot_proc + xpert1%u(i,j,k)*xpert2%u(i,j,k)
          dot_proc = dot_proc + xpert1%v(i,j,k)*xpert2%v(i,j,k)
          dot_proc = dot_proc + xpert1%t(i,j,k)*xpert2%t(i,j,k)
          dot_proc = dot_proc + xpert1%delp(i,j,k)*xpert2%delp(i,j,k)
          dot_proc = dot_proc + xpert1%sphu(i,j,k)*xpert2%sphu(i,j,k)
          dot_proc = dot_proc + xpert1%qitot(i,j,k)*xpert2%qitot(i,j,k)
          dot_proc = dot_proc + xpert1%qltot(i,j,k)*xpert2%qltot(i,j,k)
          dot_proc = dot_proc + xpert1%ozone(i,j,k)*xpert2%ozone(i,j,k)
       enddo
    enddo
 enddo

 !Gather
 rc=0
 call ESMF_VmGetCurrent(vm)
 call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_proc), recvbuf=dotp, cnt= 1, rc=status)

end subroutine dot_prod_xx

subroutine dot_prod_zz(z1,z2,nz,dotp)

 implicit none

 integer, intent(in) :: nz
 real(8), dimension(nz), intent(in) :: z1,z2
 real(8), intent(out) :: dotp
 
 type(ESMF_VM) :: vm
 integer :: k, rc, status
 real(8) :: dot_proc(1)

 dot_proc = zero

 do k = 1,nz
    dot_proc = dot_proc + z1(k)*z2(k)
 enddo

 !Gather
 rc=0
 call ESMF_VmGetCurrent(vm)
 call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_proc), recvbuf=dotp, cnt= 1, rc=status)

end subroutine dot_prod_zz

!subroutine replicate_pole_points(xpert)
!
! implicit none
! type(pertstate), intent(inout) :: xpert
!
! integer :: i,j,k
! logical :: send, recv
! integer :: my_rank, size
! integer :: commglobal, ierr, buddy, status
!
! integer, allocatable :: all_ijs(:,:)
!
! integer, allocatable :: spole_recv_ids_guess(:), spole_recv_ids(:)
! integer, allocatable :: npole_recv_ids_guess(:), npole_recv_ids(:)
!
! integer :: npole_send_id, spole_send_id
!
! integer :: tosend(5)
!
! call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
! call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
!
! tosend(1) = my_rank
! tosend(2) = ifirst
! tosend(3) = ilast
! tosend(4) = jfirst
! tosend(5) = jlast
!
! allocate(all_ijs(0:size-1,5))
! all_ijs = 0
!
! if (my_rank.eq.0) then
!
!    all_ijs(0,1) = tosend(1)
!    all_ijs(0,2) = tosend(2)
!    all_ijs(0,3) = tosend(3)
!    all_ijs(0,4) = tosend(4)
!    all_ijs(0,5) = tosend(5)
!
!    do buddy=1,size-1
!
!       call mpi_recv(tosend, 5, MPI_INT, buddy, 1,  MPI_COMM_WORLD, status, ierr)
!
!       print*, tosend(1), tosend(2), tosend(3), tosend(4), tosend(5)
!
!       all_ijs(buddy,1) = tosend(1)
!       all_ijs(buddy,2) = tosend(2)
!       all_ijs(buddy,3) = tosend(3)
!       all_ijs(buddy,4) = tosend(4)
!       all_ijs(buddy,5) = tosend(5)
!
!    end do
!
! else
!
!     call mpi_send(tosend, 5, MPI_INT, 0, 1, MPI_COMM_WORLD, ierr)
!
! end if
!
! !Make sure first communication is complete
! call mpi_barrier(MPI_COMM_WORLD, ierr)
!
! if (my_rank.eq.0) then
!
!    print*, all_ijs(:,1)
!
! endif
!
! if (my_rank.eq.0) then
!
!    do buddy=1,size-1
!       call mpi_send(all_ijs, size*5, MPI_INT, buddy, 2, MPI_COMM_WORLD, ierr)
!    end do
!
! else
!
!    call mpi_recv(all_ijs, size*5, MPI_INT, 0, 2,  MPI_COMM_WORLD, status, ierr)
!
! endif
!
!
! do i = 0,size-1
!
!    if (all_ijs(i,2) == 1 .and. all_ijs(i,4) == 1) then
!       spole_send_id = all_ijs(i,1)
!    elseif (all_ijs(i,2) == 1 .and. all_ijs(i,4) == npy) then
!       npole_send_id = all_ijs(i,1)
!    endif
!
! enddo
!
!
!
!
!
!
!
!
!
!
!
! deallocate(all_ijs)
!
!! !South pole
!! !----------
!!
!! send = .false.
!! recv = .false.
!!
!! if (ifirst == 1 .and. jfirst == 1) then
!!    send = .true.
!! elseif (ifirst /= 1 .and. jfirst == 1) then
!!    recv = .true.
!! endif
!!
!! if (send) then
!!    MPI_SEND(xpert%u(1,1,:),lm,MPI_DOUBLE_PRECISION)
!! elseif (recv) then
!!
!!
!! !Replicate along the pole points
!! if (jfirst == 1) then
!!    do 2,im
!!       xpert%u(i,1,:) = xpert%u(1,1,:)
!!       xpert%v(i,1,:) = xpert%v(1,1,:)
!!       xpert%t(i,1,:) = xpert%t(1,1,:)
!!       xpert%delp(i,1,:) = xpert%delp(1,1,:)
!!       xpert%sphu(i,1,:) = xpert%sphu(1,1,:)
!!       xpert%qitot(i,1,:) = xpert%qitot(1,1,:)
!!       xpert%qltot(i,1,:) = xpert%qltot(1,1,:)
!!       xpert%ozone(i,1,:) = xpert%ozone(1,1,:)
!!    enddo
!! endif
!!
!!
!! !North pole
!! !----------
!!
!! !Replicate along the pole points
!! if (jlast == npy) then
!!    do 2,im
!!       xpert%u(i,jm,:) = xpert%u(1,jm,:)
!!       xpert%v(i,jm,:) = xpert%v(1,jm,:)
!!       xpert%t(i,jm,:) = xpert%t(1,jm,:)
!!       xpert%delp(i,jm,:) = xpert%delp(1,jm,:)
!!       xpert%sphu(i,jm,:) = xpert%sphu(1,jm,:)
!!       xpert%qitot(i,jm,:) = xpert%qitot(1,jm,:)
!!       xpert%qltot(i,jm,:) = xpert%qltot(1,jm,:)
!!       xpert%ozone(i,jm,:) = xpert%ozone(1,jm,:)
!!    enddo
!! endif
!
!
!end subroutine replicate_pole_points



end module svecs_utils
