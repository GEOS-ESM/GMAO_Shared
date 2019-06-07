!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fvGriddedState - Distributed FV-grid state variables
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_fvGriddedState
      implicit none
      private	! except

      public :: fvGriddedState		! data structure
      public :: fvGriddedState_read
      public :: fvGriddedState_write
      public :: fvGriddedState_bdump	! binary dump
      public :: fvGriddedState_incr
      public :: fvGriddedState_diff
      public :: fvGriddedState_show

      public :: fvGriddedState_init
      public :: clean
      public :: fvGriddedState_asis	! Copy components

      public :: ptr_phis
      public :: ptr_hs_stdv
      public :: ptr_ts
      public :: ptr_lwi
      public :: ptr_ps

      public :: ptr_delp
      public :: ptr_uw
      public :: ptr_vw
      public :: ptr_cw
      public :: ptr_qall
      public :: ptr_q
      public :: ptr_oz
      public :: ptr_pt

    type fvGriddedState
      private
      	! all variables are defined to be 3+ dimensions.  For
	! single level 2-d variables, it simply means the third
	! dimension has either size 1 or 0.

	! 2-d variables

      real,pointer,dimension(:,:,:) :: phis
      real,pointer,dimension(:,:,:) :: hs_stdv
      real,pointer,dimension(:,:,:) :: ts
      real,pointer,dimension(:,:,:) :: lwi
      real,pointer,dimension(:,:,:) :: ps

	! 3-d variables

      real,pointer,dimension(:,:,:) :: delp
      real,pointer,dimension(:,:,:) :: u
      real,pointer,dimension(:,:,:) :: v
      real,pointer,dimension(:,:,:,:) :: q
      real,pointer,dimension(:,:,:) :: pt

      integer :: vtype = 0
    end type fvGriddedState

    interface fvGriddedState_read ; module procedure	&
      read_ ; end interface
    interface fvGriddedState_write; module procedure	&
      write_; end interface
    interface fvGriddedState_bdump; module procedure	&
      bdump_; end interface
    interface fvGriddedState_incr; module procedure	&
      incr_; end interface
    interface fvGriddedState_diff; module procedure	&
      diff_; end interface

    interface fvGriddedState_init ; module procedure	&
      init_ ; end interface
    interface fvGriddedState_clean; module procedure	&
      clean_; end interface
    interface clean; module procedure	&
      clean_; end interface
    interface fvGriddedState_asis ; module procedure	&
      asis_ ; end interface

    interface fvGriddedState_show ; module procedure	&
      show_ ; end interface

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       21Jun05 - Todling - changed precision of output to 32-bit
!       09Jan06  -Banglin Zhang - changed for analysis of cloud water variable 
!       15May06  -Todling - back to handling qi and ql for cloud water analysis
!       31Jan07  -Todling - add bruce force opt to apply dry mass-adjustment 
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fvGriddedState'
  integer,parameter :: VTYPE_FULL = 0
  integer,parameter :: VTYPE_INCR = 1
  integer,parameter :: VTYPE_ADJV =-1

  integer,parameter :: IPREC=0  ! 32 bits
  integer,parameter :: MINLM=4	! number of %q slots

  logical, parameter :: massadj_def =.false.  ! default: apply mass adjustment
  integer, save      :: nqmadj      = 1       ! default: mass adjustment for wv-q only (no other tracer)

!_#define SHOW_INPUTFIELDS
!_#define SHOW_OUTPUTFIELDS

#ifndef SHOW_INPUTFIELDS
#ifdef DEBUG_CHECKSUMS
#define SHOW_INPUTFIELDS
#endif
#endif

#ifndef SHOW_OUTPUTFIELDS
#ifdef DEBUG_CHECKSUMS
#define SHOW_OUTPUTFIELDS
#endif
#endif

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read_ - read a state as distributed interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read_(dynfile,ob,h_fv,intLev1,intLevs,		&
	comm,root,stat,gridverify,ignorediff,vtype,hour,massadj)
      use m_dyn,only : dyn_vect
      use m_dyn,only : dyn_get
      use m_dyn,only : dyn_clean

      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_init
      use m_fvGrid,only : fvGrid_get
      use m_fvGrid,only : fvGrid_bcast
      use m_fvGridThickness,only : fvGridThickness_verify
      use m_fvGridThickness,only : fvGridThickness_define

      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init
      use m_InterleaveScatterer,only : clean
      use m_InterleaveScatterer,only : get,localsize
      use m_InterleaveScattererComm,only : alloc_scatterv
      use m_InterleaveScattererComm,only : scatterv

      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : interleavedObject_init

      use m_massadj, only : madj_init
      use m_massadj, only : madj_pre

      use m_mpout,only : mpout_log
      use m_stdio,only : stdout
      use m_mpif90,only : MP_comm_rank,MP_type
      use m_die,only : assert_,die,MP_perr,perr,MP_die
      implicit none
      character(len=*)    ,intent(in ) :: dynfile	! dyn input
      type(fvGriddedState),intent(out) :: ob	!
      type(fvGrid)        ,intent(out) :: h_fv	! all header information
      type(interleavedObject), intent(out) :: intLev1	! distribution
      type(interleavedObject), intent(out) :: intLevs	! distribution
      integer,intent(in) :: comm,root		! communicator and root
      integer,optional,intent(out) :: stat
      logical,optional,intent(in ) :: gridverify
      logical,optional,intent(in ) :: ignorediff
      integer,optional,intent(in ) :: vtype !  0 : full
                                            ! +1 : increment
					    ! -1 : adjoint vector
      integer,optional,intent(out) :: hour
      logical,optional,intent(in)  :: massadj

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       26Jun06 - Jing Guo - bug fix in implementatin of qi/ql and qt
!	29Jan07 - Todling  - add dry mass adjustment opt(to be moved to GSI)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read_'
  integer :: nstep
  integer :: nymd,nhms,freq
  type(InterleaveScatterer) :: scatLev1,scatLevs
  type(dyn_vect) :: w_fv			! for all PEs
  				! When not on ROOT, components of w_fv
  				! are referenced in their zero-size
				! forms, since they are not expected to
				! be accessed.  See dynVect_setdummy_().
  integer ::  im,jm,km,lm,mlm
  integer :: myPE, ier
  real,allocatable,dimension(:) :: href
  integer :: mlev1,mlevs,k,m
  logical :: gridverify_
  logical :: ignorediff_
  logical :: massadj_
  integer :: vtype_

_ALLENTRY_
  if(present(stat)) stat=0

  gridverify_=.false.
  if(present(gridverify)) gridverify_=gridverify
  ignorediff_=.false.
  if(present(ignorediff)) ignorediff_=ignorediff
  vtype_=VTYPE_FULL
  if(present(vtype)) vtype_=vtype
  massadj_=massadj_def
  if(present(massadj)) massadj_=massadj

  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(myPE==root) then

	! It is assumed that the input file is always a single time
	! data.  Therefore, the time specification is used for
	! information only.  In particular, timidx=1 returns the first
	! time record in the input file, with (nymd,nhms) reporting the
	! date-time of the record.
    nymd=0;nhms=0

    call dyn_get(dynfile, nymd, nhms, w_fv, ier,	&
	nstep=nstep,timidx=1,freq=freq)
    	if(ier/=0) then
	  call perr(myname_,'dyn_get('//trim(dynfile)//')',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	! If so, calculate factors for mass adjustment (only at synoptic time)
    if ( massadj_ .and. mod(nhms,060000)==0 ) then
      call madj_init ( w_fv%grid%im, w_fv%grid%jm, w_fv%grid%km, nqmadj )
      call madj_pre  ( w_fv%grid%im, w_fv%grid%jm, w_fv%grid%km, nqmadj, w_fv%grid%ak, w_fv%grid%bk, &
                       w_fv%ps, w_fv%delp, w_fv%q, ier )
    	if(ier/=0) then
	  call perr(myname_,'mfix_pre(w_fv)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    endif

	! For now, GSI dealing w/ total could condensate so store sum of 
	! ice and liquid into ice component of state vector
    if ( size(w_fv%q,4) >= 4 ) then
	! At the input, both q3>=0. and q4>=0. are assumed for this
	! block of operations.

      w_fv%q(:,:,:,3) = w_fv%q(:,:,:,3) + w_fv%q(:,:,:,4) 
      w_fv%q(:,:,:,4) = 0.

	! Or to make separation of qcliq and qcice possible at write_(),
	! one can choose to keep the ratio as Q4=q4/(q3+q4) here:
	! where(w_fv%q(:,:,:,3)/=0.)
	!   w_fv%q(:,:,:,4) = w_fv%q(:,:,:,4)/w_fv%q(:,:,:,3) 
	! endwhere
    endif

    call fvGrid_init(h_fv,w_fv,nymd=nymd,nhms=nhms,	&
	nstep=nstep,freq=freq)
    call fvGrid_get(h_fv,im=im,jm=jm,km=km,lm=lm)

    	ASSERT(im==size(w_fv%delp,1))
    	ASSERT(jm==size(w_fv%delp,2))
    	ASSERT(km==size(w_fv%delp,3))
    	ASSERT(lm==size(w_fv%q,4))
	ASSERT(lm>=2)

    	! verify %delp against %ps and %ak,%bk

    if(gridverify_) then
      ier=0
      select case(vtype_)
      case(VTYPE_FULL)
        call fvGridThickness_verify(w_fv%delp,h_fv,w_fv%ps,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'failed, fvGridThickness_verify()',ier)
	  if(ignorediff_) then
	    call perr(myname_,'ignored, fvGridThickness_verify()')
	    ier=0
	  endif
	endif
      case(VTYPE_INCR)
      case(VTYPE_ADJV)
      case default
      end select

	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    endif

  else
    	! Although no content is needed for some variables not defined
	! on the root PE, they are defined anyway to make sure their
	! applications as actual arguments are always legally defined.
	! They are also deallocated later when they are no longer
	! referenced in the code.

    call dynVect_setdummy_(w_fv)
  endif

  call fvGrid_bcast(h_fv,root,comm)
  call fvGrid_get(h_fv,im=im,jm=jm,km=km,lm=lm,hour=hour)

  call MPI_bcast(vtype_,1,MP_type(vtype_),root,comm,ier)
    if(ier/=0) call MP_die(myname_,'MPI_bcast(vtype_)',ier)
  ob%vtype = vtype_
!________________________________________
  call InterleaveScatterer_init(scatLevs,km,comm,root)
  call InterleaveScatterer_init(scatLev1, 1,comm,root)

  	! 2-d components of a single level

  call alloc_scatterv(im,jm,scatLev1,w_fv%phis   ,ob%phis   ,	&
  						root,comm,myname)
  call alloc_scatterv(im,jm,scatLev1,w_fv%hs_stdv,ob%hs_stdv,	&
  						root,comm,myname)
  call alloc_scatterv(im,jm,scatLev1,w_fv%ts     ,ob%ts     ,	&
  						root,comm,myname)
  call alloc_scatterv(im,jm,scatLev1,w_fv%lwi    ,ob%lwi    ,	&
  						root,comm,myname)
  call alloc_scatterv(im,jm,scatLev1,w_fv%ps     ,ob%ps     ,	&
  						root,comm,myname)
  mlev1=localsize(scatLev1)

  	! 3-d components of km levels

  mlevs=localsize(scatLevs)

  call alloc_scatterv(im,jm,scatLevs,w_fv%delp,ob%delp,root,comm,myname)

  call alloc_scatterv(im,jm,scatLevs,w_fv%u   ,ob%u   ,root,comm,myname)
  call alloc_scatterv(im,jm,scatLevs,w_fv%v   ,ob%v   ,root,comm,myname)

    mlm=max(lm,MINLM)
    allocate(ob%q(im,jm,mlevs,mlm))
  ob%q(:,:,:,:)=0.
  call scatterv(scatLevs,w_fv%q,ob%q(:,:,:,1:lm),root,comm)

  call alloc_scatterv(im,jm,scatLevs,w_fv%pt,ob%pt,root,comm,myname)

	ASSERT(mlevs==size(ob%pt,3))

  	! distribution controllers

  call interleavedObject_init(intLevs,km,comm,root)
  call interleavedObject_init(intLev1, 1,comm,root)

  call dyn_clean(w_fv)

#ifdef SHOW_INPUTFIELDS
  call show_(myname_,ob)
#endif
_ALLEXIT_
end subroutine read_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dynVect_setdummy_ - create a dummy dyn_vect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dynVect_setdummy_(w)
      use m_dyn,only : dyn_vect
      use m_mpout,only : mpout_log
      implicit none
      type(dyn_vect),intent(inout) :: w

! !REVISION HISTORY:
! 	15Dec04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dynVect_setdummy_'
_ENTRY_

  	! for all 2-d fields

  allocate(w%phis   (0,0))
  allocate(w%hs_stdv(0,0))
  allocate(w%ts     (0,0))
  allocate(w%lwi    (0,0))
  allocate(w%ps     (0,0))

  	! for all 3-d fields

  allocate(w%delp(0,0,0))
  allocate(w%u   (0,0,0))
  allocate(w%v   (0,0,0))
  allocate(w%pt  (0,0,0))
  allocate(w%q (0,0,0,0))

  	! for other dynamic variables than might be referenced, e.g.
	! in dyn_clean().

  allocate(w%qm(0))
  allocate(w%grid%ak(0))
  allocate(w%grid%bk(0))
_EXIT_
end subroutine dynVect_setdummy_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: write_ - output an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine write_(dynfile,ob,h_fv,intLev1,intLevs,comm,root,massadj)
      use m_dyn,only : dyn_vect
      use m_dyn,only : dyn_put
      use m_dyn,only : dyn_clean

      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_fvGridThickness,only : fvGridThickness_verify

      use m_interleavedObject,only : InterleavedObject
      use m_interleavedObject,only : totalSize

      use m_massadj, only : madj_post
      use m_massadj, only : madj_clean

      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init
      use m_InterleaveScatterer,only : get,localSize
      use m_InterleaveScatterer,only : clean
      use m_InterleaveScattererComm,only : alloc_gatherv
      use m_mpout,only : mpout_log
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_die,die,perr
      implicit none

      character(len=*),intent(in) :: dynfile
      type(fvGriddedState),intent(in) :: ob
      type(fvGrid),intent(in) :: h_fv
      type(interleavedObject), intent(in) :: intLev1
      type(interleavedObject), intent(in) :: intLevs
      integer,intent(in) :: comm
      integer,intent(in) :: root
      logical,optional,intent(in) :: massadj

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	15May06 - Todling, adjusted assert check on q(lm)
!	29Jan07 - Todling  - add dry mass adjustment opt(to be moved to GSI)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::write_'
  type(dyn_vect) :: w_fv	! no content touching if not on root.
  integer :: im,jm,km,lm
  integer :: myPE,ier
  type(InterleaveScatterer) :: scatLev1,scatLevs
  integer :: nymd,nhms,nstep
  integer :: k,mlevs,m
  real,allocatable,dimension(:) :: href
  logical :: massadj_
_ALLENTRY_
 
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call fvGrid_get(h_fv,im=im,jm=jm,km=km,lm=lm)

  	ASSERT(lm>=2)
		! Expect w%q(...,:) contains at least q and oz
	call mpout_log(myname_,'lm=size(w_fv%q,4)',lm)

	ASSERT(km==totalSize(intLevs))
	ASSERT( 1==totalSize(intLev1))

#ifdef SHOW_OUTPUTFIELDS
  call show_(myname_,ob)
#endif
!________________________________________

  call InterleaveScatterer_init(scatLevs,km,comm,root)
  call InterleaveScatterer_init(scatLev1, 1,comm,root)

  	! 2-d components of a single level

  call alloc_gatherv(im,jm,scatLev1,ob%phis   ,w_fv%phis   ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%hs_stdv,w_fv%hs_stdv,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%ts     ,w_fv%ts     ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%lwi    ,w_fv%lwi    ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%ps     ,w_fv%ps     ,	&
  						root,comm,myname)

  	! 3-d components of km levels

  mlevs=localsize(scatLevs)
	ASSERT(mlevs==size(ob%pt,3))

  call alloc_gatherv(im,jm,scatLevs,ob%delp,w_fv%delp,root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,ob%u   ,w_fv%u   ,root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,ob%v   ,w_fv%v   ,root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,ob%q(:,:,:,1:lm),w_fv%q, &
    root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,ob%pt  ,w_fv%pt  ,root,comm,myname)

	! Upto this point, it is expected that no pointer in w_fv has
	! been touched if the process is not on the root PE.


  massadj_=massadj_def
  if(present(massadj)) massadj_=massadj

  if(myPE==root) then

    	! verify %delp against %ps and %ak,%bk

    ier=0
    select case(ob%vtype)
    case(VTYPE_FULL)
      call fvGridThickness_verify(w_fv%delp,h_fv,w_fv%ps,stat=ier)
    case(VTYPE_INCR)
    case(VTYPE_ADJV)
    case default
    end select

	if(ier/=0) then
	  call perr(myname_,'failed, fvGridThickness_verify()',ier)
	  call perr(myname_,	&
		'ignored for output, fvGridThickness_verify()')
	  ier=0
	endif

	! For now, GSI dealing w/ total could condensate so store sum of 
	! ice and liquid into ice component of state vector.  However, one
	! may choose to separate the two components with the same ratio as
	! the background.  See read_() for a comparison.
!    if ( size(w_fv%q,4) >= 4 ) then
!	w_fv%q(:,:,:,4) = w_fv%q(:,:,:,3) * w_fv%q(:,:,:,4) 
!	w_fv%q(:,:,:,3) = w_fv%q(:,:,:,3) - w_fv%q(:,:,:,4) 
!    endif

    call fvGrid_get(h_fv,nymd=nymd,nhms=nhms,nstep=nstep,dyn=w_fv)

	! If so, apply mass adjustment to analyzed field
    if ( massadj_  .and. mod(nhms,060000)==0 ) then
      call madj_post ( w_fv%grid%im, w_fv%grid%jm, w_fv%grid%km, nqmadj, 0, w_fv%grid%ptop, &
                       w_fv%ps, w_fv%delp, w_fv%q, ier )
    	if(ier/=0) call die(myname_,'madj_post()',ier)
      call madj_clean ()
    endif
    call dyn_put(dynfile,nymd,nhms,IPREC,w_fv,ier,nstep=nstep)
    	if(ier/=0) call die(myname_,'dyn_put()',ier)
    call dyn_clean(w_fv)

  endif
	! Note there is nothing but the names defined for the components
	! of w_fv not on ROOT.
_ALLEXIT_
end subroutine write_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bdump_ - binary-dump output an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine bdump_(binfile,ob,h_fv,intLev1,intLevs,comm,root)
      use m_dyn,only : dyn_vect

      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_fvGridThickness,only : fvGridThickness_verify

      use m_interleavedObject,only : InterleavedObject
      use m_interleavedObject,only : totalSize

      use m_daInterp,only : daInterp_vdtoa
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init
      use m_InterleaveScatterer,only : get,localSize
      use m_InterleaveScatterer,only : clean
      use m_InterleaveScattererComm,only : alloc_gatherv
      use m_InterleaveScattererComm,only : gatherv
      use m_mpout,only : mpout_log
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_die,die,perr
      use m_ioutil,only : luavail
      use m_realkinds,only : R4 => kind_r4
      implicit none

      character(len=*),intent(in) :: binfile
      type(fvGriddedState),intent(in) :: ob
      type(fvGrid),intent(in) :: h_fv
      type(interleavedObject), intent(in) :: intLev1
      type(interleavedObject), intent(in) :: intLevs
      integer,intent(in) :: comm
      integer,intent(in) :: root

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bdump_'
  type(dyn_vect) :: w_fv	! no content touching if not on root.
  integer :: im,jm,km,lm
  integer :: myPE,ier,lu
  type(InterleaveScatterer) :: scatLev1,scatLevs
  integer :: nymd,nhms,nstep
  integer :: k,mlevs,m
  real,allocatable,dimension(:) :: href
  real,allocatable,dimension(:,:) :: tmp_ps
  real,allocatable,dimension(:,:,:) :: tmp_u,tmp_v
  logical :: hflip_=.false. ! in case we need something different
  logical :: vdtoa_=.false. ! in case we need something different
_ALLENTRY_
 
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call fvGrid_get(h_fv,im=im,jm=jm,km=km,lm=lm)

  	ASSERT(lm>=2)
		! Expect w%q(...,:) contains at least q and oz
	call mpout_log(myname_,'lm=size(w_fv%q,4)',lm)

	ASSERT(km==totalSize(intLevs))
	ASSERT( 1==totalSize(intLev1))

#ifdef SHOW_OUTPUTFIELDS
  call show_(myname_,ob)
#endif
!________________________________________

  call InterleaveScatterer_init(scatLevs,km,comm,root)
  call InterleaveScatterer_init(scatLev1, 1,comm,root)

  	! 2-d components of a single level
  call alloc_gatherv(im,jm,scatLev1,ob%phis   ,w_fv%phis   ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%hs_stdv,w_fv%hs_stdv,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%ts     ,w_fv%ts     ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%lwi    ,w_fv%lwi    ,	&
  						root,comm,myname)
  call alloc_gatherv(im,jm,scatLev1,ob%ps     ,w_fv%ps     ,	&
  						root,comm,myname)

  	! 3-d components of km levels

  mlevs=localsize(scatLevs)
	ASSERT(mlevs==size(ob%pt,3))

  call alloc_gatherv(im,jm,scatLevs,ob%delp,w_fv%delp,root,comm,myname)

    allocate(tmp_u(size(ob%v,1),size(ob%v,2),size(ob%v,3)))
    allocate(tmp_v(size(ob%v,1),size(ob%v,2),size(ob%v,3)))
    if(vdtoa_) then
      call daInterp_vdtoa(ob%u,ob%v,tmp_u,tmp_v)
    else
       tmp_u=ob%u
       tmp_v=ob%v
    endif

  call alloc_gatherv(im,jm,scatLevs,tmp_u   ,w_fv%u   ,root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,tmp_v   ,w_fv%v   ,root,comm,myname)

    deallocate(tmp_u)
    deallocate(tmp_v)

  call alloc_gatherv(im,jm,scatLevs,ob%q   ,w_fv%q   ,root,comm,myname)
  call alloc_gatherv(im,jm,scatLevs,ob%pt  ,w_fv%pt  ,root,comm,myname)

	! Upto this point, it is expected that no pointer in w_fv has
	! been touched if the process is not on the root PE.


  if(myPE==root) then

    allocate(tmp_ps(size(w_fv%ps,1),size(w_fv%ps,2)))

    	! verify %delp against %ps and %ak,%bk

    ier=0
    select case(ob%vtype)
    case(VTYPE_FULL)
      call fvGridThickness_verify(w_fv%delp,h_fv,w_fv%ps,stat=ier)
    case(VTYPE_INCR)
        ! compute
      tmp_ps(:,:)=0.
      do k=1,size(w_fv%delp,3)
        tmp_ps(:,:)=tmp_ps(:,:)+w_fv%delp(:,:,k)
      end do
      w_fv%delp(:,:,km)=tmp_ps(:,:)-w_fv%delp(:,:,km)
      do k=km-1,1,-1
        w_fv%delp(:,:,k)=w_fv%delp(:,:,k+1)-w_fv%delp(:,:,k)
      end do
    case(VTYPE_ADJV)
    case default
    end select

	if(ier/=0) then
	  call perr(myname_,'failed, fvGridThickness_verify()',ier)
	  call perr(myname_,	&
		'ignored for output, fvGridThickness_verify()')
	  ier=0
	endif

  !!  call fvGrid_get(h_fv,nymd=nymd,nhms=nhms,nstep=nstep,dyn=w_fv)

    lu=luavail()
    open(lu,file=binfile,status='unknown', &
      form='unformatted',access='sequential',iostat=ier)
      	if(ier/=0) call die(myname_,'open("'//trim(binfile)//'")',ier)

    if(ier/=0) write(lu,iostat=ier) tmp_ps

    if(hflip_) then
      ASSERT(mod(im,2)==0)

      call hflip2d_(tmp_ps)
      call hflip3d_(w_fv%delp)
      call hflip3d_(w_fv%u)
      call hflip3d_(w_fv%v)
      call hflip3d_(w_fv%pt)
      call hflip3d_(w_fv%q(:,:,:,1))
    endif

    if(ier/=0) write(lu,iostat=ier) real(w_fv%u,kind=R4)
    if(ier/=0) write(lu,iostat=ier) real(w_fv%v,kind=R4)
    if(ier/=0) write(lu,iostat=ier) real(w_fv%pt,kind=R4)
    if(ier/=0) write(lu,iostat=ier) real(w_fv%delp,kind=R4)
    if(ier/=0) write(lu,iostat=ier) real(tmp_ps,kind=R4)
    if(ier/=0) write(lu,iostat=ier) real(w_fv%q(:,:,:,1),kind=R4)
      deallocate(tmp_ps)

    close(lu,iostat=ier)
      	if(ier/=0) call die(myname_,'close("'//trim(binfile)//'")',ier)

  endif
	! Note there is nothing but the names defined for the components
	! of w_fv not on ROOT.
_ALLEXIT_
contains
subroutine hflip2d_(f)
  implicit none
  real,dimension(:,:),intent(inout) :: f
  integer :: im,ih,j
  real,dimension(size(f,1)/2) :: d
  im=size(f,1);ih=im/2
  do j=1,size(f,2)
    d(1:ih)=f(1:ih,j)
    f(1:ih,j)=f(ih+1:im,j)
    f(ih+1:im,j)=d(1:ih)
  end do
end subroutine hflip2d_
subroutine hflip3d_(f)
  implicit none
  real,dimension(:,:,:),intent(inout) :: f
  integer :: im,ih,j
  real,dimension(size(f,1)/2) :: d
  im=size(f,1);ih=im/2
  do k=1,size(f,3)
  do j=1,size(f,2)
    d(1:ih)=f(1:ih,j,k)
    f(1:ih,j,k)=f(ih+1:im,j,k)
    f(ih+1:im,j,k)=d(1:ih)
  end do
  end do
end subroutine hflip3d_
end subroutine bdump_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: incr_ - compute analysis through increment
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine incr_(a,b,c)
      use m_mpout,only : mpout_log
      implicit none
      type(fvGriddedState),intent(inout) :: a ! a=a+b-c ! analysis
      type(fvGriddedState),intent(inout) :: b ! b=b-c   ! increment
      type(fvGriddedState),intent(in) :: c

! !REVISION HISTORY:
! 	19Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::incr_'
_ENTRY_

#ifdef DEBUG_CHECKSUMS
  call show_(myname_//':a.in',a)
  call show_(myname_//':b.in',b)
  call show_(myname_//':c.in',c)
#endif

    ! Compute the increment variables.  Note b%ps is not changed.  It
    ! is kept as the full variable for grid definition purposes.

    ! d:=b-c; b:=a; a:=a+d
  call plusdiff_sv_(a%ps     ,b%ps     ,c%ps     )
  call plusdiff_sv_(a%phis   ,b%phis   ,c%phis   )
  call plusdiff_sv_(a%hs_stdv,b%hs_stdv,c%hs_stdv)
  call plusdiff_sv_(a%lwi    ,b%lwi    ,c%lwi    )

    ! b:=b-c; a:=a+b
  call plusdiff_  (a%ts  ,b%ts  ,c%ts  )
  call plusdiff_  (a%delp,b%delp,c%delp)
  call plusdiff_  (a%u   ,b%u   ,c%u   )
  call plusdiff_  (a%v   ,b%v   ,c%v   )
  call plusdiff4d_(a%q   ,b%q   ,c%q   )
  a%q = max(0.0, a%q)
  call plusdiff_  (a%pt  ,b%pt  ,c%pt  )

  b%vtype = VTYPE_INCR

#ifdef DEBUG_CHECKSUMS
  call show_(myname_//':a.out',a)
#endif

_EXIT_
contains
subroutine plusdiff_sv_(a,b,c)
  ! -- hope to have a better cache performance than array expressions.
  implicit none
  real,dimension(:,:,:),intent(inout) :: a
  real,dimension(:,:,:),intent(inout) :: b
  real,dimension(:,:,:),intent(in) :: c
  integer :: i,j,k
  real :: d
    ! d:=b-c; b:=a; a:=a+d
  do k=lbound(a,3),ubound(a,3)
  do j=lbound(a,2),ubound(a,2)
  do i=lbound(a,1),ubound(a,1)
    d=b(i,j,k)-c(i,j,k)
    b(i,j,k)=a(i,j,k)
    a(i,j,k)=a(i,j,k)+d
  end do
  end do
  end do
end subroutine plusdiff_sv_
subroutine plusdiff_(a,b,c)
  ! -- hope to have a better cache performance than array expressions.
  implicit none
  real,dimension(:,:,:),intent(inout) :: a
  real,dimension(:,:,:),intent(inout) :: b
  real,dimension(:,:,:),intent(in) :: c
  integer :: i,j,k
  do k=lbound(a,3),ubound(a,3)
  do j=lbound(a,2),ubound(a,2)
    b(:,j,k)=b(:,j,k)-c(:,j,k)
    a(:,j,k)=a(:,j,k)+b(:,j,k)
  end do
  end do
end subroutine plusdiff_
subroutine plusdiff4d_(a,b,c)
  ! -- hope to have a better cache performance than array expressions.
  implicit none
  real,dimension(:,:,:,:),intent(inout) :: a
  real,dimension(:,:,:,:),intent(inout) :: b
  real,dimension(:,:,:,:),intent(in) :: c
  integer :: i,j,k,l
  do l=lbound(a,4),ubound(a,4)
  do k=lbound(a,3),ubound(a,3)
  do j=lbound(a,2),ubound(a,2)
    b(:,j,k,l)=b(:,j,k,l)-c(:,j,k,l)
    a(:,j,k,l)=a(:,j,k,l)+b(:,j,k,l)
  end do
  end do
  end do
end subroutine plusdiff4d_
end subroutine incr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diff_ - compute difference fields
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine diff_(b,c)
      use m_mpout,only : mpout_log
      implicit none
      type(fvGriddedState),intent(inout) :: b ! b=b-c   ! increment
      type(fvGriddedState),intent(in) :: c

! !REVISION HISTORY:
! 	19Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diff_'
_ENTRY_

#ifdef DEBUG_CHECKSUMS
  call show_(myname_//':b.in',b)
  call show_(myname_//':c.in',c)
#endif

    ! Compute difference variables.  Note b%ps is not changed.  It
    ! is kept as the full variable for grid definition purposes.

    ! b:=b-c
  b%ts  =b%ts  -c%ts
  b%delp=b%delp-c%delp
  b%u   =b%u   -c%u
  b%v   =b%v   -c%v
  b%q   =b%q   -c%q
  b%pt  =b%pt  -c%pt

  b%vtype = VTYPE_INCR

#ifdef DEBUG_CHECKSUMS
  call show_(myname_//':b.out',b)
#endif

_EXIT_
end subroutine diff_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialized an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,h_fv,intLev1,intLevs)
      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : localSize
      use m_mpout,only : mpout_log
      implicit none
      type(fvGriddedState),intent(out) :: ob
      type(fvGrid),intent(in) :: h_fv
      type(interleavedObject),intent(in) :: intLev1	! single-level
      type(interleavedObject),intent(in) :: intLevs	! multi-levels

! !REVISION HISTORY:
! 	22Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: im,jm,lm
  integer :: mlev1,mlevs
_ENTRY_

  call fvGrid_get(h_fv,im=im,jm=jm,lm=lm)
  mlev1=localSize(intLev1)		! 0 or 1
  mlevs=localSize(intLevs)

	! 2-d variables

	allocate(ob%phis   (im,jm,mlev1))
	allocate(ob%hs_stdv(im,jm,mlev1))
	allocate(ob%ts	   (im,jm,mlev1))
	allocate(ob%lwi	   (im,jm,mlev1))
	allocate(ob%ps	   (im,jm,mlev1))

	! 3-d variables

	allocate(ob%delp(im,jm,mlevs))
	allocate(ob%u	(im,jm,mlevs))
	allocate(ob%v	(im,jm,mlevs))
	allocate(ob%q	(im,jm,mlevs,lm))
	allocate(ob%pt	(im,jm,mlevs))
_EXIT_
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(fvGriddedState),intent(inout) :: ob

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_
	deallocate(ob%phis   )
	deallocate(ob%hs_stdv)
	deallocate(ob%ts     )
	deallocate(ob%lwi    )
	deallocate(ob%ps     )

	deallocate(ob%delp)
	deallocate(ob%u   )
	deallocate(ob%v   )
	deallocate(ob%pt  )
	deallocate(ob%q   )
_EXIT_
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: asis_ - copy components of the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine asis_(ob,ref)
      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : localSize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fvGriddedState),intent(inout) :: ob
      type(fvGriddedState),intent(in ) :: ref

! !REVISION HISTORY:
! 	22Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::asis_'
_ENTRY_
  ASSERT(all(shape(ob%phis)==shape(ref%phis)))
!TORM  ASSERT(all(shape(ob%hs_stdv)==shape(ref%hs_stdv)))
!TORM  ASSERT(all(shape(ob%ts  )==shape(ref%ts  )))
!TORM  ASSERT(all(shape(ob%lwi )==shape(ref%lwi )))
!TORM  ASSERT(all(shape(ob%ps  )==shape(ref%ps  )))

  ASSERT(all(shape(ob%delp)==shape(ref%delp)))
!TORM  ASSERT(all(shape(ob%u   )==shape(ref%u   )))
!TORM  ASSERT(all(shape(ob%v   )==shape(ref%v   )))
!TORM  ASSERT(all(shape(ob%q   )==shape(ref%q   )))
!TORM  ASSERT(all(shape(ob%pt  )==shape(ref%pt  )))

	ob%phis=ref%phis
	ob%hs_stdv=ref%hs_stdv
	ob%ts  =ref%ts
	ob%lwi =ref%lwi
	ob%ps  =ref%ps

	ob%delp=ref%delp
	ob%u   =ref%u
	ob%v   =ref%v
	ob%q   =ref%q
	ob%pt  =ref%pt
_EXIT_
end subroutine asis_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: show_ - show local checksums of all state variables
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine show_(where,ob)
      use m_checksums,only : checksums_show
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in) :: where
      type(fvGriddedState),intent(in) :: ob	!

! !REVISION HISTORY:
! 	09Jun05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::show_'
  integer :: ier,jm,l,m
  character(len=*),parameter :: lq='123456789*'
_ENTRY_

call mpout_log(myname_,where)

  call checksums_show(ob%phis,'UNKW','fv:%phis')
  call checksums_show(ob%ps  ,'UNKW','fv:%ps'  )
  call checksums_show(ob%ts  ,'UNKW','fv:%ts'  )

  call checksums_show(ob%delp,'UNKW','fv:%delp')
  jm=size(ob%u,2)
  call checksums_show(ob%u(:,2:jm,:),'UNKW','fv:%u@d-grid')
  call checksums_show(ob%v,'UNKW','fv:%v@d-grid' )

  call checksums_show(ob%pt,'UNKW','fv:%pt')
  do l=1,size(ob%q,4)
    m=min(l,10)
    call checksums_show(ob%q(:,:,:,l),'UNKW','fv:%q('//lq(m:m)//')')
  end do

_EXIT_
end subroutine show_

function ptr_phis (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_phis
  ptr_phis => ob%phis
end function ptr_phis
function ptr_hs_stdv (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_hs_stdv
  ptr_hs_stdv => ob%hs_stdv
end function ptr_hs_stdv
function ptr_ts (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_ts
  ptr_ts => ob%ts
end function ptr_ts
function ptr_lwi (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_lwi
  ptr_lwi => ob%lwi
end function ptr_lwi
function ptr_ps (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_ps
  ptr_ps => ob%ps
end function ptr_ps
function ptr_delp (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_delp
  ptr_delp => ob%delp
end function ptr_delp
function ptr_uw (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_uw
  ptr_uw => ob%u
end function ptr_uw
function ptr_vw (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_vw
  ptr_vw => ob%v
end function ptr_vw
function ptr_cw (ob)
  use m_die,only : assert_
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_cw
  character(len=*),parameter :: myname_=myname//"::ptr_cw"
  ASSERT(size(ob%q,4)>=3)
  ptr_cw => ob%q(:,:,:,3)
end function ptr_cw
function ptr_qall (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:,:) :: ptr_qall
  ptr_qall => ob%q
end function ptr_qall
function ptr_q (ob)
  use m_die,only : assert_
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_q
  character(len=*),parameter :: myname_=myname//"::ptr_q"
  ASSERT(size(ob%q,4)>=1)
  ptr_q => ob%q(:,:,:,1)
end function ptr_q
function ptr_oz (ob)
  use m_die,only : assert_
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_oz
  character(len=*),parameter :: myname_=myname//"::ptr_oz"
  ASSERT(size(ob%q,4)>=2)
  ptr_oz => ob%q(:,:,:,2)
end function ptr_oz
function ptr_pt (ob)
  implicit none
  type(fvGriddedState),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_pt
  ptr_pt => ob%pt
end function ptr_pt
end module m_fvGriddedState
