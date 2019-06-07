!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_swapij - swap the order of the first two indices.
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_swapij
      implicit none
      private	! except

      public :: swapij		! data structure

    interface swapij; module procedure	&
      swapij2_,	&
      swapij3_; end interface

! !REVISION HISTORY:
! 	21Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- moved out m_fgInterp
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_swapij'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: swapij3_ - i-j swap
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine swapij3_(aij,aji)
      use m_die,only : assert_
      implicit none
      real,dimension(:,:,:),intent(in ) :: aij
      real,dimension(:,:,:),intent(out) :: aji

! !REVISION HISTORY:
! 	16Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::swapij3_'
  integer :: i,j,k,isz,jsz

  ASSERT(size(aij,1)==size(aji,2))
  ASSERT(size(aij,2)==size(aji,1))
  ASSERT(size(aij,3)==size(aji,3))

  isz=size(aij,1)
  jsz=size(aij,2)

  do k=1,size(aij,3)
    do i=1,isz		! loop ordered is tuned to optimize aji(:,:)
      aji(1:jsz,i,k)=aij(i,1:jsz,k)
    end do
  end do
end subroutine swapij3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: swapij2_ - i-j swap
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine swapij2_(aij,aji)
      use m_die,only : assert_
      implicit none
      real,dimension(:,:),intent(in ) :: aij
      real,dimension(:,:),intent(out) :: aji

! !REVISION HISTORY:
! 	16Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::swapij2_'
  integer :: i,j,isz,jsz

  ASSERT(size(aij,1)==size(aji,2))
  ASSERT(size(aij,2)==size(aji,1))

  isz=size(aij,1)
  jsz=size(aij,2)

  do i=1,isz		! loop ordered is tuned to optimize aji(:,:)
    aji(1:jsz,i)=aij(i,1:jsz)
  end do
end subroutine swapij2_
end module m_swapij
