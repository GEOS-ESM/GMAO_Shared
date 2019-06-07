!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_InterleaveScatterer - Scatterer of an InterleavedObject
!
! !DESCRIPTION:
!	An InterleaveScatterer defines communication operations (see
!   m_InterleaveScatterComm also), such as scatterv() from a given
!   root-only array to an interleaved distributed array, or gatterv()
!   from an interleaved distributed array to a root-only array.
!
! !INTERFACE:
!#include "regime.H"

    module m_InterleaveScatterer
      use m_interleavedObject,only : interleavedObject
      implicit none
      private	! except

      public :: InterleaveScatterer		! data structure
      public :: InterleaveScatterer_init,init	! initialize an object
      public :: InterleaveScatterer_clean,clean	! clean an object

      public :: get			! get() index of a local element
      public :: deepcopy		! copy a whole object
      public :: locate			! locate() any element
      public :: inquire			! for communication parameters
      public :: totalSize		! global object size
      public :: localSize		! local object size
      public :: maximSize		! the maximum size of all PEs.

      public :: ptr_interleave

    type InterleaveScatterer
      private
      type(interleavedObject),pointer :: interleave	! the base type
      integer :: nPEs
      integer :: myPE
      integer :: root		! root of the scatterer
      integer :: comm		! the communicator
    end type InterleaveScatterer

    interface InterleaveScatterer_init; module procedure	&
    	init_; end interface
    interface init; module procedure init_; end interface
    interface InterleaveScatterer_clean; module procedure	&
    	clean_; end interface
    interface clean; module procedure clean_; end interface

    interface get; module procedure get_; end interface
    interface deepcopy; module procedure deepcopy_; end interface
    interface locate; module procedure locate_; end interface
    interface inquire; module procedure inquire_; end interface
    interface totalSize; module procedure totalSize_; end interface
    interface localSize; module procedure localSize_; end interface
    interface maximSize; module procedure maximSize_; end interface

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_InterleaveScatterer'

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,necount,comm,root)
      use m_interleavedObject,only : interleavedObject_init
      use m_mpout,only : mpout_log
      use m_mpif90,only : MP_comm_size,MP_comm_rank
      use m_die,only : MP_die
      implicit none
      type(InterleaveScatterer),intent(out) :: ob
      integer,intent(in) :: necount	! no. of ecount
      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where necount is defined

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier
_ALLENTRY_

  	allocate(ob%interleave)
  call interleavedObject_init(ob%interleave,necount,comm,root)

  ob%comm=comm
  ob%root=root
  	call MP_comm_size(comm,ob%nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,ob%myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

_ALLEXIT_
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_InterleavedObject,only : clean
      use m_mpif90,only : MP_COMM_NULL
      implicit none
      type(InterleaveScatterer),intent(inout) :: ob

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call clean(ob%interleave)
  	deallocate(ob%interleave)
  ob%root=-1
  ob%comm=MP_COMM_NULL
  ob%myPE=-1
  ob%nPEs=-1
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: inquire_ - inquire communicator info.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine inquire_(ob,comm,root)
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer,optional,intent(out) :: root
      integer,optional,intent(out) :: comm

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inquire_'
  if(present(comm)) comm=ob%comm
  if(present(root)) root=ob%root
end subroutine inquire_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get the global index of a local level
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(ob,lindex,gindex,dstrcount,		&
    	rootcounts,rootdispls,rootlbound,rootubound)
      use m_InterleavedObject,only : get
      use m_InterleavedObject,only : maximSize
      use m_InterleavedObject,only : localSize
      use m_InterleavedObject,only : totalSize
      use m_interleavedObject,only : ptr_where
      use m_interleavedObject,only : ptr_which
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer,intent(in) :: lindex	! given a local index

      		! get its corresponding "root" information
      integer,optional,intent(out) :: gindex	! its global index
      integer,optional,intent(out) :: dstrcount	! scattered count

      		! counts/displs pair is only significant on root.
      integer,optional,dimension(0:),intent(out) :: rootcounts
      integer,optional,dimension(0:),intent(out) :: rootdispls

      		! global lbound:ubound on all PEs of the current lindex
		! they are significant only on root.  However, for
		! convenience, they are also defined as for a zero-sized
		! array if not on root.
      integer,optional,intent(out) :: rootlbound
      integer,optional,intent(out) :: rootubound

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: ubnd,lbnd,i
  integer,pointer,dimension(:) :: iptr
  integer :: msize,lsize,gsize
_ENTRY_

  msize=maximSize(ob%interleave)
  lsize=localSize(ob%interleave)
  gsize=totalSize(ob%interleave)

	ASSERT(lindex>=1.and.lindex<=msize)

  	! get the global index of the given local element
  if(present(gindex)) then
    gindex=0
    if(lindex<=lsize) call get(ob%interleave,lindex,gindex)
  endif

  	! get the local (scattered) element count.
  if(present(dstrcount)) then
    dstrcount=0
    if(lindex<=lsize) dstrcount=1
  endif

	! set counts/displs on root
    		! set default counts/displs to skip empty PEs.

  if(ob%myPE==ob%root) then

    call get(ob%interleave,lindex,gindex=i,lbound=lbnd,ubound=ubnd)
	! (lbound:ubound) is the range of levels with the given
	! local index (lindex) on different processors, in actual
	! level index (gindex).

#ifndef NDEBUG
    ASSERT(i>=lbnd.and.i<=ubnd)

    	! lindex of all gindex in [lbnd:ubnd] should be the same.
    iptr=>ptr_which(ob%interleave)		! all iPEs(:)
    ASSERT(all(iptr(lbnd:ubnd)==lindex))
	nullify(iptr)
#endif

    if(present(rootcounts).or.present(rootdispls)) then
      iptr=>ptr_where(ob%interleave)		! all iPEs(:)
    	! *gindex*'s of the given lindex are range [lbnd:ubnd].
	! the locations of the range are %where(lbnd:ubnd).

      if(present(rootcounts)) then
        rootcounts(:)=0
        rootcounts(iptr(lbnd:ubnd))=1
      endif

      if(present(rootdispls)) then
        rootdispls(:)=0
        do i=lbnd,ubnd
          rootdispls(iptr(i))=i-lbnd
        end do
      endif
      nullify(iptr)
    endif


    if(present(rootlbound)) rootlbound=lbnd
    if(present(rootubound)) rootubound=ubnd
    	! when on root, rootlbound:rootubound represents a real array
	! segment, bufr(lbnd:ubnd).
  else
    if(present(rootlbound).or.present(rootubound)) then
  	! when not on root, rootlbound:rootubound represents a zero-
	! size array segment, bufr(lbnd:lbnd-1).
      call get(ob%interleave,lindex,gindex=lbnd)
      if(present(rootlbound)) rootlbound=lbnd
      if(present(rootubound)) rootubound=lbnd-1
    endif
  endif
_EXIT_
end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: locate_ - locate any given level by its global index
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine locate_(ob,gindex,iPE,lindex)
      use m_interleavedObject,only : locate
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer,intent(in) :: gindex	! a level by its global index
      integer,intent(out) :: iPE	! where it is in processor ID
      integer,intent(out) :: lindex	! where it is on that processor

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::locate_'
  call locate(ob%interleave,gindex,iPE,lindex)
end subroutine locate_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: deepcopy_ - deep copy an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine deepcopy_(a,b)
      use m_interleavedObject,only : deepcopy
      implicit none
      type(InterleaveScatterer),intent(in) :: a
      type(InterleaveScatterer),intent(out) :: b

! !REVISION HISTORY:
! 	27Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::deepcopy_'

  b=a
  allocate(b%interleave)
  call deepcopy(a%interleave,b%interleave)
end subroutine deepcopy_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: totalSize_ - total size in elements
!
! !DESCRIPTION:
!
! !INTERFACE:

    function totalSize_(ob)
      use m_InterleavedObject,only : totalSize
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer :: totalSize_

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::totalSize_'
  totalSize_=totalSize(ob%interleave)
end function totalSize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: localSize_ - local size in elements
!
! !DESCRIPTION:
!
! !INTERFACE:

    function localSize_(ob)
      use m_InterleavedObject,only : localSize
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer :: localSize_

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::localSize_'
  localSize_=localSize(ob%interleave)
end function localSize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: maximSize_ - the maximum local size of all processors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function maximSize_(ob)
      use m_interleavedObject,only : maximSize
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      integer :: maximSize_

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::maximSize_'
  maximSize_=maximSize(ob%interleave)
end function maximSize_

function ptr_interleave(ob)
  use m_interleavedObject,only : interleavedObject
  implicit none
  type(InterleaveScatterer),target,intent(in) :: ob
  type(interleavedObject),pointer :: ptr_interleave

  ptr_interleave => ob%interleave
end function ptr_interleave
end module m_InterleaveScatterer
