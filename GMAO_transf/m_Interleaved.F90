!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Interleaved - object of interleaving distribution
!
! !DESCRIPTION:
!   This object describe an interleave distribution pattern of
!
! !INTERFACE:
!#include "regime.H"

    module m_Interleaved
      implicit none
      private	! except

      public :: Interleaved		! data structure
      public :: Interleaved_init, init	! initialize an object
      public :: Interleaved_clean,clean	! clean an object

      public :: get			! get() index of a local level
      public :: locate			! locate() any level
      public :: inquire			! for communication parameters
      public :: totalSize		! total size in elements
      public :: localSize		! local size in elements
      public :: maximSize		! the maximum localSize of all

      public :: deepcopy		! copy a whole object

      public :: ptr_indices		! acc/upd. of %indices(:)
      public :: ptr_where		! acc/upd. of %where(:)
      public :: ptr_which		! acc/upd. of %which(:)

    type Interleaved
      private
      integer :: root = -1
      integer :: myPE
      integer :: nPEs
      integer :: next ! the PE where the "next" element may go
      integer :: mecount	! "my" element count (on this PE)
      integer :: mxcount	! maximum element count (on all PEs)
      integer :: necount	! total element count (on all PEs)

      				! local-to-total index
      integer,pointer,dimension(:) :: indices	! (mecount)
      				! total-to-local index
      integer,pointer,dimension(:) :: where	! (necount)
      integer,pointer,dimension(:) :: which	! (necount)
    end type Interleaved

    interface Interleaved_init
      module procedure init_
    end interface
    interface init; module procedure init_; end interface
    interface Interleaved_clean
      module procedure clean_
    end interface
    interface clean; module procedure clean_; end interface

    interface get      ; module procedure get_      ; end interface
    interface locate   ; module procedure locate_   ; end interface
    interface inquire; module procedure inquire_; end interface
    interface totalSize; module procedure totalSize_; end interface
    interface localSize; module procedure localSize_; end interface
    interface maximSize; module procedure maximSize_; end interface

    interface deepcopy; module procedure deepcopy_; end interface

    interface ptr_indices; module procedure ptr_indices_; end interface
    interface ptr_where; module procedure ptr_where_; end interface
    interface ptr_which; module procedure ptr_which_; end interface

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Interleaved'

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

    subroutine init_(ob,necount,comm,root,next)
      use m_mpif90,only : MP_comm_size,MP_comm_rank,MP_type
      use m_mpout,only : mpout_log
      use m_die   ,only : MP_die,die
      implicit none
      type(Interleaved),intent(out) :: ob
      integer,intent(in) :: necount	! no. of ecount
      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where necount is defined
      integer,optional,intent(out) :: next ! next returns the PE ID
      					   ! suggesting the location of
					   ! the %necount+1-th element.

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: rem
  integer :: myPE,nPEs
  integer :: ier,i,j,iPE
_ALLENTRY_

  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

  ob%root=root
  ob%myPE=myPE
  ob%nPEs=nPEs
  iPE=mod(myPE+nPEs-ob%root,nPEs)	! count offset from root

  ob%necount=necount
  call MPI_bcast(ob%necount,1,MP_type(ob%necount), root,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast()',ier)

  	! count on this PE

  rem=mod(ob%necount,nPEs)		! remainder
  ob%mecount=ob%necount/nPEs
  if(iPE<rem) ob%mecount=ob%mecount+1	! me-count
  ob%next=mod(ob%root+rem,nPEs)
  if(present(next)) next=ob%next

  	! the maximum count of all PEs
  ob%mxcount=(ob%necount+nPEs-1)/nPEs	! mx-count

  allocate( ob%indices(ob%mecount),	&
  	    ob%where(ob%necount),	&
	    ob%which(ob%necount), stat=ier)
  	if(ier/=0) call die(myname_,	&
		'allocate(%indices,%where,%which)',ier)

  do j=1,ob%mecount		! define indices of local elements
    ob%indices(j)=1+iPE+(j-1)*nPEs
  end do

  j=1		! which/j is-in [1,mecount]
  iPE=0		! where/iPE is-in [0,nPEs-1]
  do i=1,ob%necount		! define mapping of distributed elements
    ob%where(i)=mod(iPE+root,nPEs) ! starting from root
    ob%which(i)=j

    iPE=iPE+1
    if(iPE==nPEs) then
      iPE=0
      j=j+1
    endif
  end do

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
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none
      type(Interleaved),intent(inout) :: ob

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier
_ENTRY_

  ob%root=-1
  ob%myPE=-1
  ob%nPEs=0
  ob%necount=0
  ob%mecount=0
  ob%mxcount=0

  deallocate(ob%indices,ob%where,ob%which,stat=ier)
  	if(ier/=0) call die(myname_,'deallocate()',ier)

_EXIT_
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

    subroutine inquire_(ob,root,size,rank,next)
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,optional,intent(out) :: root
      integer,optional,intent(out) :: size
      integer,optional,intent(out) :: rank
      integer,optional,intent(out) :: next

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inquire_'
  if(present(size)) size=ob%nPEs
  if(present(root)) root=ob%root
  if(present(rank)) rank=ob%myPE
  if(present(next)) next=ob%next
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

    subroutine get_(ob,lindex,gindex,lcount,lbound,ubound, &
      root,dstrcount,rootcounts,rootdispls,rootlbound,rootubound)
      use m_mpout,only : mpout_log
      use m_die,only : die,assert_
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,intent(in) :: lindex	! a given local level reference
      integer,optional,intent(out):: gindex	! its global index
      integer,optional,intent(out):: lcount	! 0|1, whether present
      integer,optional,intent(out):: lbound	! lbound of this index
      integer,optional,intent(out):: ubound	! ubound of this index

        ! Output arguments below are defined for scatter or gather.
	! Their values depend on if the processor is root.
      integer,optional,intent(in) :: root
      integer,optional,intent(out) :: dstrcount ! ==lcount, 1 or 0
        ! for scatterv, the recv count
	! for gatherv, the send count
      integer,dimension(0:),optional,intent(out) :: rootcounts
      integer,dimension(0:),optional,intent(out) :: rootdispls
        ! for scatterv, (scounts, sdispls)
        ! for gatherv, (rcounts, rdispls)
      integer,optional,intent(out) :: rootlbound
      integer,optional,intent(out) :: rootubound
        ! for scatterv, range of the send buffer
	! for gatherv, range of the recv buffer


! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: lcount_,root_
  integer :: lbnd,ubnd
  integer :: i,l
_ENTRY_

  	ASSERT(lindex>0.and.lindex<=ob%mxcount)

  lcount_=0
  if(lindex>=1.and.lindex<=ob%mecount) lcount_=1

  if(present(lcount)) lcount=lcount_
  if(present(gindex)) then
    gindex=0 ! in case there is no element (lcount_==0)
    if(lcount_==1) gindex=ob%indices(lindex)
  endif

  	! The element range [lbound:ubound] at this index.

  lbnd=(lindex-1)*ob%nPEs+1
  ubnd=min(lindex*ob%nPEs,ob%necount)
  if(present(lbound)) lbound=lbnd
  if(present(ubound)) ubound=ubnd

  	! Output arguments specified below are defined differently
	! depending which processor the process is on.
  root_=ob%root
  if(present(root)) root_=root

  if(present(dstrcount)) dstrcount=lcount_

  if(ob%myPE==root_) then
    if(present(rootlbound)) rootlbound=lbnd
    if(present(rootubound)) rootubound=ubnd

    if(present(rootcounts)) then
      rootcounts(:)=0
        ! rootcounts(ob%where(lbnd:ubnd))=1
      do i=lbnd,ubnd
	l=ob%where(i)
        rootcounts(l)=1
      end do
    endif
    if(present(rootdispls)) then
      rootdispls(:)=0
        ! rootdispls(ob%where(lbnd:ubnd))=(/(i,i=0,ubnd-lbnd)/)
      do i=lbnd,ubnd
        l=ob%where(i)
        rootdispls(l)=i-lbnd
      end do
    endif

  else
    if(present(rootlbound)) rootlbound=lbnd
    if(present(rootubound)) rootubound=lbnd-1
    if(present(rootcounts)) rootcounts=0
    if(present(rootdispls)) rootdispls=0
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

    subroutine locate_(ob,gindex,where,which)
      use m_mpout,only : mpout_log
      use m_die,only : die,assert_
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,intent(in) :: gindex	! an element by its global index
      integer,intent(out) :: where	! where it is in processor ID
      integer,intent(out) :: which	! where it is on that processor

! !REVISION HISTORY:
! 	13Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::locate_'
_ENTRY_

  	ASSERT(gindex>0.and.gindex<=ob%necount)

  where=ob%where(gindex)
  which=ob%which(gindex)
_EXIT_
end subroutine locate_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: totalsize_ - total size of the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    function totalsize_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer :: totalsize_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::totalsize_'
_ENTRY_
  totalsize_=ob%necount
_EXIT_
end function totalsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: localsize_ - local size of the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    function localsize_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer :: localsize_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::localsize_'
_ENTRY_
  localsize_=ob%mecount
_EXIT_
end function localsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: maximsize_ - local size of the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    function maximsize_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer :: maximsize_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::maximsize_'
_ENTRY_
  maximsize_=ob%mxcount
_EXIT_
end function maximsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: deepcopy_ - copy the contents of dynamic components.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine deepcopy_(a,b)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in ) :: a
      type(Interleaved),intent(out) :: b

! !REVISION HISTORY:
! 	16Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::deepcopy_'
_ENTRY_

  b=a

  allocate(b%indices(size(a%indices)))
  allocate(b%where(size(a%where)))
  allocate(b%which(size(a%which)))

  b%indices=a%indices
  b%where  =a%where
  b%which  =a%which

_EXIT_
end subroutine deepcopy_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_indices_ - accessor/updater of %indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_indices_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,pointer,dimension(:) :: ptr_indices_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_indices_'
_ENTRY_
  ptr_indices_ => ob%indices(:)
_EXIT_
end function ptr_indices_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_where_ - accessor/updater of %where
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_where_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,pointer,dimension(:) :: ptr_where_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_where_'
_ENTRY_
  ptr_where_ => ob%where(:)
_EXIT_
end function ptr_where_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_which_ - accessor/updater of %which
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_which_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(Interleaved),intent(in) :: ob
      integer,pointer,dimension(:) :: ptr_which_

! !REVISION HISTORY:
! 	15Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_which_'
_ENTRY_
  ptr_which_ => ob%which(:)
_EXIT_
end function ptr_which_
end module m_Interleaved
