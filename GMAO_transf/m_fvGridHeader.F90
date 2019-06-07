!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fvGridHeader - read and distribute FV-grid header
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_fvGridHeader
      implicit none
      private	! except

      public :: fvGridHeader_allread
      public :: fvGridHeader_read

    interface fvGridHeader_allread ; module procedure	&
      allread_ ; end interface
    interface fvGridHeader_read ; module procedure	&
      read_ ; end interface

! !REVISION HISTORY:
! 	01Apr05	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fvGridHeader'
#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allread_ - read and distribute the "header".
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allread_(dynfile,ob,comm,root,stat,gridverify,ignorediff)
      use m_dyn,only : dyn_vect
      use m_dyn,only : dyn_get
      use m_dyn,only : dyn_clean

      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_init
      use m_fvGrid,only : fvGrid_get
      use m_fvGrid,only : fvGrid_bcast
      use m_fvGridThickness,only : fvGridThickness_verify

      use m_mpout ,only : mpout_log
      use m_mpif90,only : MP_comm_rank,MP_type
      use m_die,only : assert_,die,MP_perr,perr
      implicit none
      character(len=*)    ,intent(in ) :: dynfile	! dyn input file
      type(fvGrid)        ,intent(out) :: ob	! all header information
      integer,intent(in) :: comm,root		! communicator and root
      integer,optional,intent(out) :: stat	! status at return
      logical,optional,intent(in ) :: gridverify	! verify eta?
      logical,optional,intent(in ) :: ignorediff	! ignore verify?

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allread_'
  integer :: nstep,freq
  integer :: nymd,nhms
  type(dyn_vect) :: w_fv
  integer ::  im,jm,km,lm
  integer :: myPE, ier
_ALLENTRY_

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(myPE==root) call read_(dynfile,ob,stat=ier,	&
	gridverify=gridverify,ignorediff=ignorediff	)

  call MPI_bcast(ier,1,MP_type(ier),root,comm,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(ier)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call fvGrid_bcast(ob,root,comm)
_ALLEXIT_
end subroutine allread_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read_ - read and distribute the "header".
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read_(dynfile,ob,stat,gridverify,ignorediff)
      use m_dyn,only : dyn_vect
      use m_dyn,only : dyn_get
      use m_dyn,only : dyn_clean

      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_init
      use m_fvGrid,only : fvGrid_get
      use m_fvGridThickness,only : fvGridThickness_verify

      use m_mpout ,only : mpout_log
      use m_die,only : assert_,die,perr
      implicit none
      character(len=*)    ,intent(in ) :: dynfile	! dyn input file
      type(fvGrid)        ,intent(out) :: ob	! all header information
      integer,optional,intent(out) :: stat	! status at return
      logical,optional,intent(in ) :: gridverify	! verify eta?
      logical,optional,intent(in ) :: ignorediff	! ignore verify?

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read_'
  integer :: nstep,freq
  integer :: nymd,nhms
  type(dyn_vect) :: w_fv
  integer ::  im,jm,km,lm
  integer :: myPE, ier
  logical :: gridverify_
  logical :: ignorediff_
_ENTRY_

  if(present(stat)) stat=0

  gridverify_=.false.
  if(present(gridverify)) gridverify_=gridverify
  ignorediff_=.false.
  if(present(ignorediff)) ignorediff_=ignorediff

	! It is assumed that the input file is always a single time
	! data.  Therefore, the time specification is used only for
	! information only.  In particular, timidx=1 returns the first
	! time record in the input file, with (nymd,nhms) reporting the
	! date-time of the record.

    call dyn_get(dynfile, nymd, nhms, w_fv, ier, freq=freq,	&
    					nstep=nstep,timidx=1)
    	if(ier/=0) then
	  call perr(myname_,'dyn_get('//trim(dynfile)//')',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    call fvGrid_init(ob,w_fv,nymd=nymd,nhms=nhms,nstep=nstep,freq=freq)
    call fvGrid_get(ob,im=im,jm=jm,km=km,lm=lm)

    	ASSERT(im==size(w_fv%delp,1))
    	ASSERT(jm==size(w_fv%delp,2))
    	ASSERT(km==size(w_fv%delp,3))
	ASSERT(lm>=2)

    if(gridverify_) then
      call fvGridThickness_verify(w_fv%delp,ob,w_fv%ps,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'failed, fvGridThickness_verify()',ier)
	  if(ignorediff_) then
	    call perr(myname_,'ignored, fvGridThickness_verify()')
	    ier=0
	  endif
	endif

	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    endif

    call dyn_clean(w_fv)
_EXIT_
end subroutine read_
end module m_fvGridHeader
