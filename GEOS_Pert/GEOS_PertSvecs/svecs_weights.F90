module svecs_weights

#include "MAPL_Generic.h"

use ESMF
use MAPL_mod
use svecs_cf

implicit none
private

public compute_weights

contains

subroutine compute_weights(svfile,TrajPtrI,TrajPtrF)

 implicit none

 type(trajPointer), intent(in) :: TrajPtrI, TrajPtrF

 character(len=ESMF_MAXSTR) :: svfile
 type (ESMF_Config)         :: cf
 integer                    :: rc
 integer                    :: status

 character(len=ESMF_MAXSTR)   :: Iam="compute_weights"

 integer                    :: vweight
 integer                    :: i,j,k

 real(8) :: dtheta, dthetatot
 real(8) :: dlambda, dlambdatot
 real(8), allocatable, dimension(:,:,:) :: peI, peF, levlwI, levlwF

! Compute vertical weighting
! --------------------------

 cf = ESMF_ConfigCreate(rc=status)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = status )

 call ESMF_ConfigGetAttribute( cf, vweight, label='vweight:', default=1, rc = status)
 VERIFY_(status)

 allocate(peI(im,jm,lm+1), stat=status)
 VERIFY_(status)

 allocate(peF(im,jm,lm+1), stat=status)
 VERIFY_(status)

 allocate(levlwI(im,jm,lm), stat=status)
 VERIFY_(status)
 levlwI = 0.0

 allocate(levlwF(im,jm,lm), stat=status)
 VERIFY_(status)
 levlwF = 0.0

 !Compute pe
 peI = 0.0_8
 peI(:,:,1) = ptop
 do k = 2,lm+1
    peI(:,:,k) = peI(:,:,k-1) + dble(TrajPtrI%delp(:,:,k-1))
 enddo

 peF = 0.0_8
 peF(:,:,1) = ptop
 do k = 2,lm+1
    peF(:,:,k) = peF(:,:,k-1) + dble(TrajPtrF%delp(:,:,k-1))
 enddo

 
 if (vweight == 1) then
    if (MAPL_AM_I_ROOT()) print*, 'Vertical weighting with delta sigma'

    do k = projlevI(1),projlevI(2)
       levlwI(:,:,k) = dble(TrajPtrI%delp(:,:,k)) / (peI(:,:,projlevI(2)+1) - peI(:,:,projlevI(1)))
    enddo
    do k = projlevF(1),projlevF(2)
       levlwF(:,:,k) = dble(TrajPtrF%delp(:,:,k)) / (peF(:,:,projlevF(2)+1) - peF(:,:,projlevF(1)))
    enddo

 elseif (vweight == 2) then
    if (MAPL_AM_I_ROOT()) print*, 'Vertical weighting with height'

!dh later

 endif

 deallocate(peI,peF)


! Compute area weighting
! ----------------------

 dtheta  = MAPL_PI / dble(npy-1)
 dlambda = 2.0_8 * MAPL_PI / dble(npx)


 !Initial time projection weighting
 !---------------------------------

 gridweightI = 0.0_8

 dthetatot = sin(min(projlatgridI(2)*d2r+half*dtheta,90.0_8*d2r)) - sin(max(projlatgridI(1)*d2r-half*dtheta,-90.0_8*d2r))
 dlambdatot = (projlongridI(2)*d2r + half*dlambda) - (projlongridI(1)*d2r - half*dlambda) !Cyclic so no max/min required

 do i = 1,im
    do j = 1,jm
       do k = projlevI(1),projlevI(2)

          if ( lons(i) >= projlongridI(1) .and. lons(i) <= projlongridI(2) .and. &
               lats(j) >= projlatgridI(1) .and. lats(j) <= projlatgridI(2) ) then

             if (j == 1 .and. jfirst == 1) then
                gridweightI(i,j,k) = levlwI(i,j,k) * &
                                     (one/dthetatot) * (sin(theta(j)+half*dtheta) + one) * &
                                     (one/dlambdatot) * dlambda
             elseif (j == jm .and. jlast == npy) then
                gridweightI(i,j,k) = levlwI(i,j,k) * &
                                     (one/dthetatot) * (one - sin(theta(j)-half*dtheta)) * &
                                     (one/dlambdatot) * dlambda
             else
                gridweightI(i,j,k) = levlwI(i,j,k) * &
                                     (one/dthetatot) * (sin(theta(j)+half*dtheta) - sin(theta(j)-half*dtheta)) * &
                                     (one/dlambdatot) * dlambda
             endif

          endif

       enddo
    enddo
 enddo


 !Final time projection weighting
 !-------------------------------

 gridweightF = 0.0_8

 dthetatot = sin(min(projlatgridF(2)*d2r + half*dtheta, 90.0_8*d2r)) - sin(max(projlatgridF(1)*d2r-half*dtheta,-90.0_8*d2r))
 dlambdatot = (projlongridF(2)*d2r + half*dlambda) - (projlongridF(1)*d2r - half*dlambda) !Cyclic so no max/min required

 do i = 1,im
    do j = 1,jm
       do k = projlevF(1),projlevF(2)

          if ( lons(i) >= projlongridF(1) .and. lons(i) <= projlongridF(2) .and. &
               lats(j) >= projlatgridF(1) .and. lats(j) <= projlatgridF(2) ) then

             if (j == 1 .and. jfirst == 1) then
                gridweightF(i,j,k) = levlwF(i,j,k) * &
                                     (one/dthetatot) * (sin(theta(j)+half*dtheta) + one) * &
                                     (one/dlambdatot) * dlambda
             elseif (j == jm .and. jlast == npy) then
                gridweightF(i,j,k) = levlwF(i,j,k) * &
                                     (one/dthetatot) * (one - sin(theta(j)-half*dtheta)) * &
                                     (one/dlambdatot) * dlambda
             else
                gridweightF(i,j,k) = levlwF(i,j,k) * &
                                     (one/dthetatot) * (sin(theta(j)+half*dtheta) - sin(theta(j)-half*dtheta)) * &
                                     (one/dlambdatot) * dlambda
             endif

          endif

       enddo
    enddo
 enddo


!Deallocate locals
!-----------------

 deallocate(levlwI,levlwF)

end subroutine compute_weights

end module svecs_weights
