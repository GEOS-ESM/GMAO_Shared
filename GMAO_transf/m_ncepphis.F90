!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_NCEPphis - Reader of NCEP phis on a lat-lon grid
!
! !DESCRIPTION:
!
!   Read file _fname_, which is a NCEP defined surface geopotential.
!
! !INTERFACE:

    module m_NCEPphis
      implicit none
      private	! except

      public :: NCEPphis_read	! input an object

    interface NCEPphis_read; module procedure mpread_; end interface

! !REVISION HISTORY:
! 	16Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_NCEPphis'
#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mpread_ - a multiple PE read() interface
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mpread_(fname,intLev1,phis,comm)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_die
      implicit none
      character(len=*),intent(in) :: fname	! where to read
      type(interleavedObject),intent(in) :: intLev1	! distributation
      real,dimension(:,:,:)  ,intent(out):: phis	! on S-A-L
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	16Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mpread_'
  integer :: myPE,ier
  integer :: im,jm
  integer :: mintLev1,nintLev1

  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  im=size(phis,1)
  jm=size(phis,2)

  mintLev1=localsize(intLev1)
  nintLev1=totalsize(intLev1)

	ASSERT(nintLev1==1)
	ASSERT(mintLev1==size(phis,3))

  if(mintLev1==1) then
    call get_NCEP_phis(fname,phis(:,:,1),im,jm)
  endif

end subroutine mpread_
end module m_NCEPphis
