!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_date - make fvDAS date based on SSI date

!
! !INTERFACE:
!
   subroutine make_date ( idate, nymd, nhms, ndt )

   implicit none
   integer, intent(in) :: idate(4)  ! SSI date
   integer, intent(out) :: nymd     ! input date
   integer, intent(out) :: nhms     !   "   time
   integer, intent(in)  :: ndt      ! time increment

!
!EOP
!-------------------------------------------------------------------------

   ! idate(1) = hour
   ! idate(2) = month
   ! idate(3) = day
   ! idate(4) = year

   nymd = idate(4)*10000 + idate(2)*100 + idate(3)
   nhms = idate(1)*10000

   call tick (nymd, nhms, ndt*3600)

   end subroutine

