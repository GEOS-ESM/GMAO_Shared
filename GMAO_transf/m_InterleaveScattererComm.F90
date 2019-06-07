!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_InterleaveScattererComm - communcation operations
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_InterleaveScattererComm
      implicit none
      private	! except

      public :: alloc_scatterv
      public :: alloc_gatherv
      public :: scatterv
      public :: gatherv

      interface alloc_scatterv; module procedure	&
      	alloc_scatterv2_,	&
      	alloc_scatterv3_,	&
      	alloc_scatterv4_; end interface
      interface alloc_gatherv ; module procedure	&
      	alloc_gatherv2_ ,	&
      	alloc_gatherv3_ ,	&
      	alloc_gatherv4_ ; end interface
      interface scatterv; module procedure	&
      	scatterv2dr_,	&
      	scatterv3dr_,	&
      	scatterv4dr_; end interface
      interface gatherv ; module procedure	&
      	gatherv2dr_ ,	&
      	gatherv3dr_ ,	&
      	gatherv4dr_ ; end interface

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_InterleaveScattererComm'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_scatterv4_ - alloc before scatterv of rank 4 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_scatterv4_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank,MP_type
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:,:,:),intent(in) :: send ! from root
      real,pointer,dimension(:,:,:,:)            :: recv ! to all
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_scatterv4_'
  integer :: ier
  integer :: km,lm
  integer :: myPE

  	! This code avoids to ever touch the content of send if the
	! PE is not ROOT, although its name can be referenced as an
	! argument of a procedure, as long as the procedure will not
	! touch the content of it if the PE is not ROOT.

  if(present(stat)) stat=0

  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MP_comm_rank()', ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  if(myPE==root) then
    km=totalSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
  	ASSERT(km==size(send,3))
    lm=size(send,4)
  endif
	call MPI_bcast(lm,1,MP_type(lm),root,comm,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MPI_bcast()', ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  km=localSize(ob)
	allocate(recv(im,jm,km,lm))
		! for alloc() tracking, do
		!	call mall_mci(recv,name)

  call scatterv4dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'scatterv4dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine alloc_scatterv4_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_scatterv3_ - alloc before scatterv of rank 3 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_scatterv3_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:,:),intent(in) :: send	! from root
      real,pointer,dimension(:,:,:)            :: recv	! to all
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_scatterv3_'
  integer :: ier
  integer :: km
  integer :: myPE

  	! This code avoids to ever touch the content of send if the
	! PE is not ROOT, although its name can be referenced as an
	! argument of a procedure, as long as the procedure will not
	! touch the content of it if the PE is not ROOT.

  if(present(stat)) stat=0

  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MP_comm_rank()', ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  if(myPE==root) then
    km=totalSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
  	ASSERT(km==size(send,3))
  endif

  km=localSize(ob)
	allocate(recv(im,jm,km))
		! for alloc() tracking, do
		!	call mall_mci(recv,name)

  call scatterv3dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'scatterv3dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine alloc_scatterv3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_scatterv2_ - alloc before scatterv of rank 2 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_scatterv2_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:  ),intent(in) :: send	! from root
      real,pointer,dimension(:,:,:)            :: recv	! to all
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_scatterv2_'
  integer :: ier
  integer :: km
  integer :: myPE

  	! This code avoids to ever touch the content of send if the
	! PE is not ROOT, although its name can be referenced as an
	! argument of a procedure, as long as the procedure will not
	! touch the content of it if the PE is not ROOT.

  if(present(stat)) stat=0

  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MP_comm_rank()', ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  if(myPE==root) then
    km=totalSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
	ASSERT(km==1)
  endif

  km=localSize(ob)
	allocate(recv(im,jm,km))
		! for alloc() tracking, CALL MALL_MCI(RECV,NAME) here.

  call scatterv2dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'scatterv2dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine alloc_scatterv2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatterv4dr_ - scatterv of a rank 3 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatterv4dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,die,perr,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:,:,:),intent(in) :: send
      real,dimension(:,:,:,:),intent(out):: recv
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatterv4dr_'
  real,dimension(0,0,0) :: sdum
  integer :: ier,l,myPE
!_______________________________________________________________________

  if(present(stat)) stat=0
	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MP_comm_rank()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  if(myPE==root) then
	ASSERT(size(send,4)==size(recv,4))
  endif

  do l=1,size(recv,4)
    if(myPE==root) then
      call scatterv3dr_(ob,send(:,:,:,l),	&
			   recv(:,:,:,l),root,comm,stat=ier)
    else
      call scatterv3dr_(ob,sdum,		&
			   recv(:,:,:,l),root,comm,stat=ier)
    endif

    if(ier/=0) then
      call perr(myname_,'scatterv3dr_()',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
  end do
end subroutine scatterv4dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatterv3dr_ - scatterv of a rank 3 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatterv3dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : get
      use m_InterleaveScatterer,only : maximSize
      use m_InterleaveScatterer,only : totalSize
      use m_InterleaveScatterer,only : localSize
      use m_InterleaveScatterer,only : inquire
      use m_mpif90,only : MP_type,MP_comm_rank,MP_comm_size
      use m_die,only : assert_,die,perr,MP_die,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: send
      real,dimension(:,:,:),intent(out):: recv
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatterv3dr_'
  integer,parameter :: NDIM=3
  integer :: ier
  integer :: myPE,nPEs
  integer :: ob_comm,ob_root
  integer :: nelem,melem,lelem
  integer :: isend,lsend,lblock_send
  integer :: irecv,lrecv,lblock_recv
  integer,allocatable,dimension(:) :: scounts,sdispls
  integer :: rcount
  integer :: mtype
!_______________________________________________________________________

  if(present(stat)) stat=0
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!________________________________________

  nelem=totalSize(ob)
  melem=localSize(ob)

  if(myPE==root) then
	ASSERT(nelem==size(send,NDIM))
  endif

	ASSERT(melem==size(recv,NDIM))

! What are not checked, include
!	size(send,1)==size(recv,1)
!	size(send,2)==size(recv,2)
! It is because one can only check them on ROOT, unless additional
! MPI_bcast() calls are invoked.  On the other hand, checking them
! only on root, may left a false safe impression if there is no error
! message.
!________________________________________

  call inquire(ob,comm=ob_comm,root=ob_root)
	ASSERT(comm==ob_comm)
	ASSERT(root==ob_root)
!________________________________________

	! The significance of lblock_recv and lblock_send is
	! conditional, depending on the counts of the receiving and
	! sending buffers.  Therefore, local lblock_send and lblock_recv
	! and all-PE consistancy of ! lblock_recv and lblock_send.

	! All-PE consistency of lblock_recv is assumed where actual data
	! is to be received.  Actually checking it would require some
	! special and non-trivial handling where receiving buffer may
	! not be required.

lblock_recv=size(recv)/max(1,size(recv,NDIM))
if(myPE==root) then
  lblock_send=size(send)/max(1,size(send,NDIM))
  allocate(scounts(0:nPEs-1),sdispls(0:nPEs-1))

  if(melem>0 .and. lblock_recv/=lblock_send)		&
  	call die(myname_,'lblock_recv',lblock_recv,	&
  			 'lblock_send',lblock_send	)
else
  lblock_send=0
  allocate(scounts(0),sdispls(0))
endif
!_______________________________________________________________________

do irecv=1,maximSize(ob)
  	! To scatter messages from root to local PEs, get related
	! messaging information.

  call get(ob,irecv,dstrcount=rcount,		&
  	rootcounts=scounts,rootdispls=sdispls,	&
	rootlbound=isend  ,rootubound=lsend	)

  lrecv=irecv+rcount-1
  mtype=MP_type(send)

  	! Scatter the messages . ..

  call MPI_scatterv(	&
      send(:,:,isend:lsend),scounts*lblock_send,	&
    			    sdispls*lblock_send,mtype,	&
      recv(:,:,irecv:lrecv),rcount *lblock_recv,mtype,root,comm,ier)
     	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
end do
!_______________________________________________________________________

deallocate(scounts,sdispls)
end subroutine scatterv3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatterv2dr_ - scatterv of a rank 2 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatterv2dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : get
      use m_InterleaveScatterer,only : maximSize
      use m_InterleaveScatterer,only : totalSize
      use m_InterleaveScatterer,only : localSize
      use m_InterleaveScatterer,only : inquire
      use m_mpif90,only : MP_type,MP_comm_rank,MP_comm_size
      use m_die,only : assert_,die,perr,MP_die,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:  ),intent(in) :: send
      real,dimension(:,:,:),intent(out):: recv
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatterv2dr_'
  integer,parameter :: NDIM=3
  integer :: ier
  integer :: myPE,nPEs
  integer :: ob_comm,ob_root
  integer :: nelem,melem,lelem
  integer :: isend,lsend,lblock_send
  integer :: irecv,lrecv,lblock_recv
  integer,allocatable,dimension(:) :: scounts,sdispls
  integer :: rcount
  integer :: mtype
!_______________________________________________________________________

  if(present(stat)) stat=0
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!________________________________________

  nelem=totalSize(ob)
  	ASSERT(nelem==1)
  melem=maximSize(ob)
  	ASSERT(melem==1)
  melem=localSize(ob)

	ASSERT(melem==size(recv,NDIM))

! See comment about size(.,1) and size(.,2) checking in scatterv3dr_().
!________________________________________

  call inquire(ob,comm=ob_comm,root=ob_root)
	ASSERT(comm==ob_comm)
	ASSERT(root==ob_root)
!________________________________________

	! The significance of lblock_recv and lblock_send is
	! conditional, depending on the counts of the receiving and
	! sending buffers.  Therefore, local lblock_send and lblock_recv
	! and all-PE consistancy of ! lblock_recv and lblock_send.

	! All-PE consistency of lblock_recv is assumed where actual data
	! is to be received.  Actually checking it would require some
	! special and non-trivial handling where receiving buffer may
	! not be required.

lblock_recv=size(recv)/max(1,size(recv,NDIM))
if(myPE==root) then
  lblock_send=size(send)
  allocate(scounts(0:nPEs-1),sdispls(0:nPEs-1))

  if(melem>0 .and. lblock_recv/=lblock_send)		&
  	call die(myname_,'lblock_recv',lblock_recv,	&
  			 'lblock_send',lblock_send	)
else
  lblock_send=0
  allocate(scounts(0),sdispls(0))
endif
!_______________________________________________________________________

irecv=1
  	! To scatter messages from root to local PEs, get related
	! messaging information.

  call get(ob,irecv,dstrcount=rcount,		&
  	rootcounts=scounts,rootdispls=sdispls,	&
	rootlbound=isend  ,rootubound=lsend	)

  lrecv=irecv+rcount-1
  mtype=MP_type(send)

  	! Scatter the messages . ..

  call MPI_scatterv(	&
    send		 ,scounts*lblock_send,		&
    			  sdispls*lblock_send,mtype,	&
    recv(:,:,irecv:lrecv),rcount *lblock_recv,mtype,root,comm,ier)
     	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!_______________________________________________________________________

deallocate(scounts,sdispls)
end subroutine scatterv2dr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_gatherv4_ - alloc before gatherv of rank 4 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_gatherv4_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:,:,:),intent(in) :: send ! from all
      real,pointer,dimension(:,:,:,:)            :: recv ! to root
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_gatherv4_'
  integer :: ier
  integer :: km,lm
  integer :: myPE
  real,allocatable,dimension(:,:,:,:) :: rdum

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myPE,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  km=localSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
  	ASSERT(km==size(send,3))
  lm=size(send,4)

  km=totalSize(ob)
  if(myPE==root) then
	allocate(recv(im,jm,km,lm))
	! for alloc() tracking, call mall_mci(recv,name) here

    call gatherv4dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv4dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  else
	! This structure is there to make sure recv is not touched
	! in anyway unless it is on root.

	allocate(rdum(0,0,0,0))
    call gatherv4dr_(ob,send,rdum,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv4dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	deallocate(rdum)
  endif

end subroutine alloc_gatherv4_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_gatherv3_ - alloc before gatherv of rank 3 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_gatherv3_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:,:),intent(in) :: send	! from all
      real,pointer,dimension(:,:,:)            :: recv	! to root
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_gatherv3_'
  integer :: ier
  integer :: km
  integer :: myPE
  real,allocatable,dimension(:,:,:) :: rdum

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myPE,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  km=localSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
  	ASSERT(km==size(send,3))

  km=totalSize(ob)
  if(myPE==root) then
	allocate(recv(im,jm,km))
	! for alloc() tracking, call mall_mci(recv,name) here

    call gatherv3dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv3dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  else
	allocate(rdum(0,0,0))
    call gatherv3dr_(ob,send,rdum,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv3dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	deallocate(rdum)
  endif

end subroutine alloc_gatherv3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_gatherv2_ - alloc before gatherv of rank 2 REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_gatherv2_(im,jm,ob,send,recv,root,comm,name,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : localSize,totalSize
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,MP_perr,perr,die
      implicit none
      integer,intent(in) :: im,jm
      type(InterleaveScatterer),intent(in) :: Ob
      real,        dimension(:,:,:),intent(in) :: send	! from all
      real,pointer,dimension(:,:  )            :: recv	! to root
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),intent(in) :: name
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_gatherv2_'
  integer :: ier
  integer :: km
  integer :: myPE
  real,allocatable,dimension(:,:) :: rdum

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myPE,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  km=localSize(ob)
  	ASSERT(im==size(send,1))
  	ASSERT(jm==size(send,2))
  	ASSERT(km==size(send,3))

  km=totalSize(ob)
  	ASSERT(km==1)

  if(myPE==root) then
	allocate(recv(im,jm))
	! for alloc() tracking, call mall_mci(recv,name) here

    call gatherv2dr_(ob,send,recv,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv2dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  else
	allocate(rdum(0,0))
    call gatherv2dr_(ob,send,rdum,root,comm,stat=ier)
  	if(ier/=0) then
	  call perr(myname_,'gatherv2dr_('//trim(name)//')', ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	deallocate(rdum)
  endif

end subroutine alloc_gatherv2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv4dr_ - gatherv of a rank 3 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv4dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_mpif90,only : MP_comm_rank
      use m_die,only : assert_,die,perr,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:,:,:),intent(in) :: send
      real,dimension(:,:,:,:),intent(out):: recv
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Jun05	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv4dr_'
  real,dimension(0,0,0) :: rdum
  integer :: ier,l,myPE
!_______________________________________________________________________

  if(present(stat)) stat=0
	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MP_comm_rank()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
!________________________________________

  do l=1,size(send,4)
    if(myPE==root) then
      call gatherv3dr_(ob,send(:,:,:,l),	&
			  recv(:,:,:,l),root,comm,stat=ier)
    else
      call gatherv3dr_(ob,send(:,:,:,l),	&
			  rdum,         root,comm,stat=ier)
    endif
    if(ier/=0) then
      call perr(myname_,'scatterv3dr_()',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
  end do

end subroutine gatherv4dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv3dr_ - gatherv of a rank 3 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv3dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : get
      use m_InterleaveScatterer,only : maximSize
      use m_InterleaveScatterer,only : totalSize
      use m_InterleaveScatterer,only : localSize
      use m_InterleaveScatterer,only : inquire
      use m_mpif90,only : MP_type,MP_comm_rank,MP_comm_size
      use m_die,only : assert_,die,perr,MP_die,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: send
      real,dimension(:,:,:),intent(out):: recv
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv3dr_'
  integer,parameter :: NDIM=3
  integer :: ier
  integer :: myPE,nPEs
  integer :: scat_comm,scat_root
  integer :: nelem,melem,lelem
  integer :: isend,lsend,lblock_send
  integer :: irecv,lrecv,lblock_recv
  integer :: scount
  integer,allocatable,dimension(:) :: rcounts,rdispls
  integer :: mtype
!_______________________________________________________________________

  if(present(stat)) stat=0
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!________________________________________

  nelem=totalSize(ob)
  melem=localSize(ob)
	ASSERT(melem==size(send,NDIM))

  if(myPE==root) then
    	ASSERT(nelem==size(recv,NDIM))
  endif

! See comment about size(.,1) and size(.,2) checking in scatterv3dr_().
!________________________________________

  call inquire(ob,comm=scat_comm,root=scat_root)
	ASSERT(comm==scat_comm)
	ASSERT(root==scat_root)
!________________________________________

lblock_send=size(send)/max(1,size(send,NDIM))
if(myPE==root) then
  lblock_recv=size(recv)/max(1,size(recv,NDIM))
  allocate(rcounts(0:nPEs-1),rdispls(0:nPEs-1))
  if(lblock_send/=lblock_recv)	&
  	call die(myname_,'lblock_send',lblock_send,	&
			 'lblock_recv',lblock_recv	)
else
  lblock_recv=0
  allocate(rcounts(0),rdispls(0))
endif
!________________________________________

do isend=1,maximSize(ob)
  	! This is a message passing from local PE to root.

  call get(ob,isend,dstrcount=scount,		&
  	rootcounts=rcounts,rootdispls=rdispls,	&
	rootlbound=irecv  ,rootubound=lrecv	)

  lsend=isend+scount-1

  mtype=MP_type(send)
  call MPI_gatherv(	&
    send(:,:,isend:lsend),scount *lblock_send,mtype,	&
    recv(:,:,irecv:lrecv),rcounts*lblock_recv,		&
			  rdispls*lblock_recv,mtype, root,comm,ier)
     	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
end do

deallocate(rcounts,rdispls)
end subroutine gatherv3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv2dr_ - gatherv of a rank 2 REAL array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv2dr_(ob,send,recv,root,comm,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : get
      use m_InterleaveScatterer,only : maximSize
      use m_InterleaveScatterer,only : totalSize
      use m_InterleaveScatterer,only : localSize
      use m_InterleaveScatterer,only : inquire
      use m_mpif90,only : MP_type,MP_comm_rank,MP_comm_size
      use m_die,only : assert_,die,perr,MP_die,MP_perr
      implicit none
      type(InterleaveScatterer),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: send	! from all
      real,dimension(:,:  ),intent(out):: recv	! to root
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	18Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv2dr_'
  integer,parameter :: NDIM=3
  integer :: ier
  integer :: myPE,nPEs
  integer :: scat_comm,scat_root
  integer :: nelem,melem,lelem
  integer :: isend,lsend,lblock_send
  integer :: irecv,lrecv,lblock_recv
  integer :: scount
  integer,allocatable,dimension(:) :: rcounts,rdispls
  integer :: mtype
!_______________________________________________________________________

  if(present(stat)) stat=0
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!________________________________________

  nelem=totalSize(ob)
  	ASSERT(nelem==1)
  melem=maximSize(ob)
  	ASSERT(melem==1)
  melem=localSize(ob)
  	ASSERT(melem==size(send,NDIM))

! See comment about size(.,1) and size(.,2) checking in scatterv3dr_().
!________________________________________

  call inquire(ob,comm=scat_comm,root=scat_root)
	ASSERT(comm==scat_comm)
	ASSERT(root==scat_root)
!________________________________________

lblock_send=size(send)/max(1,size(send,NDIM))
if(myPE==root) then
  lblock_recv=size(recv)
  allocate(rcounts(0:nPEs-1),rdispls(0:nPEs-1))
  if(lblock_send/=lblock_recv)	&
  	call die(myname_,'lblock_send',lblock_send,	&
			 'lblock_recv',lblock_recv	)
else
  lblock_recv=0
  allocate(rcounts(0),rdispls(0))
endif
!________________________________________

isend=1
  	! This is a message passing from local PE to root.

  call get(ob,isend,dstrcount=scount,		&
  	rootcounts=rcounts,rootdispls=rdispls,	&
	rootlbound=irecv  ,rootubound=lrecv	)

  lsend=isend+scount-1

  mtype=MP_type(send)
  call MPI_gatherv(	&
    send(:,:,isend:lsend),scount *lblock_send,mtype,	&
    recv		 ,rcounts*lblock_recv,		&
			  rdispls*lblock_recv,mtype, root,comm,ier)
     	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
!________________________________________

deallocate(rcounts,rdispls)
end subroutine gatherv2dr_

end module m_InterleaveScattererComm
