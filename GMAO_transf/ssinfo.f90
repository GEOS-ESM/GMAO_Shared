!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ssinfo: print info from SSI spectral file header
!
! !INTERFACE:

   program ssinfo
 
! !USES:

   use m_ss
   use util
   implicit NONE

! !DESCRIPTION: Reads header from Spectral files and echoes information.
!
! !REVISION HISTORY:
!
!  2002.10.30  C. Cruz:   Initial Code
!
!-------------------------------------------------------------------------
!EOP

   character(len=*), parameter :: myname = 'ssinfo'

   type(ss_vect)  :: ss     ! analysis vector in spectral space

!  Grid

!  Locals

   integer :: n,iargc,rc
   integer :: nymd, nhms, ndt
   integer :: im,jm,km,lm,tt

!  File names

   character(len=255) :: ssifile

! start

   call init ( ssifile )

   call myout ( 6, myname, 'read meta data')
   call ss_init ( ssifile, ss, rc )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot read from file', rc )
   end if

! clean up

   call myout ( 6, myname, 'Clean up' )
   call ss_clean ( ss )
   
   call myout ( 6, myname, '-- ssinfo.x has successfully ended --' )
   call exit(0)

   CONTAINS

   subroutine init ( ssifile )

   implicit NONE

   character*255, intent(out) :: ssifile

   character*4, parameter :: myname = 'init'
   character*255 :: res
   integer :: n,iargc,nargs
   character*255 :: argv
   character*255, allocatable :: arg(:)

   nargs =  iargc()
   if( nargs.lt.1 ) then 
     print *,' Usage: ssinfo.x ssfile'
     stop
   end if
   allocate ( arg(nargs) )

! process options

   do n=1,nargs
     call getarg(n,arg(n))
   enddo
   ssifile=trim(arg(1))

   deallocate(arg)
   rc = 0

   end subroutine init

   end program ssinfo
