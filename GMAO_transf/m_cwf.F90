!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_cwf - cloud water fraction data file i/o
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_cwf
      implicit none
      private	! except

      public :: cwf_get
      public :: cwf_put

    interface cwf_get; module procedure get_; end interface
    interface cwf_put; module procedure put_; end interface

! !REVISION HISTORY:
! 	19Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_cwf'
#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - read data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(file,nymd,nhms,cwf)
      use m_die,only : die
      implicit none
      character(len=*),intent(in) :: file
      integer,intent(in) :: nymd
      integer,intent(in) :: nhms
      real,dimension(:,:,:),intent(out) :: cwf

! !REVISION HISTORY:
! 	19Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: ier

  call get_cwf(file,nymd,nhms,cwf,	&
  	size(cwf,1),size(cwf,2),size(cwf,3),ier)
	if(ier/=0) call die(myname_,'get_cwf('//trim(file)//')',ier)

end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: put_ - write data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine put_(file,nymd,nhms,cwf)
      use m_die,only : die
      implicit none
      character(len=*),intent(in) :: file
      integer,intent(in) :: nymd
      integer,intent(in) :: nhms
      real,dimension(:,:,:),intent(in) :: cwf

! !REVISION HISTORY:
! 	19Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::put_'

! nothing implemented so far.
  call die(myname_,'not implemented, file="'//trim(file)//'"')
end subroutine put_

   subroutine get_cwf ( fname, nymd, nhms, cwf, im, jm, km, rc )

   implicit none
   
   character(len=*) :: fname
   integer  nymd, nhms
   integer  im, jm, km
   integer  rc, ngatts

   real     cwf(im,jm,km)

!  Local variables
!  ---------------
   character(len=255)              :: title, source, contact, levunits
   character(len=255), allocatable :: vname(:), vtitle(:), vunits(:)

   real,    allocatable :: lat(:), lon(:), lev(:)
   real,    allocatable :: valid_range(:,:), packing_range(:,:)
   integer, allocatable :: kmvar(:), yyyymmdd(:), hhmmss(:)

   integer, parameter :: READ_ONLY = 1
   integer :: nymd1, nhms1
   integer :: im1, jm1, km1, lm1, nvars, timinc
   integer :: fid, ierr
   real    :: amiss
   
   character(len=*), parameter :: cwfname = 'CLOUD'
   rc = 0

!  Open the file
!  -------------
   call GFIO_Open ( fname, READ_ONLY, fid, ierr )
   if ( ierr .ne. 0 ) then
      rc = 1
      return
   end if

!  Get dimensions
!  --------------
   call GFIO_DimInquire ( fid, im1, jm1, km1, lm1, nvars, ngatts, ierr)
   if ( ierr .ne. 0 ) then
      rc = 2
      call GFIO_close ( fid, ierr )
      return
   end if
   if ( im.ne.im1 .or. jm.ne.jm1 .or. km.ne.km1 ) then
      rc = 3
      call GFIO_close ( fid, ierr )
      return
   end if

   call init_ ( ierr )
   if ( ierr .ne. 0 ) then
      rc = 4
      call GFIO_close ( fid, ierr )
      return
   endif

!  Get file attributes
!  -------------------
   call GFIO_Inquire ( fid, im1, jm1, km1, lm1, nvars,     &
                       title, source, contact, amiss,  &
                       lon, lat, lev, levunits,        &
                       yyyymmdd, hhmmss, timinc,       &
                       vname, vtitle, vunits, kmvar,   &
                       valid_range , packing_range, ierr )
   if ( ierr .ne. 0 ) then
      rc = 5
      call GFIO_close ( fid, ierr )
      call clean_
      return
   end if

!  Set 1st time in file as alternative time to extract data
!  in case data not found at current time
!  --------------------------------------------------------
   nymd1 = yyyymmdd(1)
   nhms1 =   hhmmss(1)

!  Try to extract data at current time
!  -----------------------------------
   call GFIO_GetVar ( fid, cwfname,  nymd, nhms,   &
                      im, jm, 1, km, cwf,     ierr )
    if( ierr/=0 ) then  ! if fails, try getting data for 1st time in file
        call GFIO_GetVar ( fid, cwfname,  nymd1, nhms1, &
                           im, jm, 1, km, cwf,     ierr )
          if(ierr/=0) then  ! if it fails again, give up!!
             rc = 6
             call GFIO_close ( fid, ierr )
             call clean_
             return
          else
             print *, ' Retrieved Cloud-Water Fraction '
          endif
    endif

!  Close GFIO file
!  ---------------
   call GFIO_close ( fid, ierr )
   call clean_()
   return

   CONTAINS

     subroutine init_ ( err )       ! allocates local memory
     integer err
     allocate ( lat(jm1), lon(im1), lev(km1), yyyymmdd(lm1), hhmmss(lm1),    &
              vname(nvars), vunits(nvars), vtitle(nvars), kmvar(nvars), &
              valid_range(2,nvars), packing_range(2,nvars),             &
              stat=err )
     end subroutine init_

     subroutine clean_()             ! de-allocates local memory
     deallocate ( lat, lon, lev, yyyymmdd, hhmmss,   &
                  vname, vunits, vtitle, kmvar,      &
                  valid_range, packing_range,        &
                  stat=ierr )
     end subroutine clean_

   end subroutine get_cwf
end module m_cwf
