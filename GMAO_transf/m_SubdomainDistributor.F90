!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SubdomainDistributor - distributor for GSI subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_SubdomainDistributor
      implicit none
      private	! except

      public :: SubdomainDistributor		! data structure
      public :: SubdomainDistributor_init,init	! initialize an object
      public :: SubdomainDistributor_clean,clean ! finalize an object
      public :: get		! get values
      public :: show		! show values for verification purposes
      public :: bufferSize	! buffer size
      public :: localSize	! local subdomain size

      public :: Subdomain_pack		! preparing a message buffer
      public :: Subdomain_unpack		! unpack a message buffer

      				! message configurations.
      public :: ptr_buffercounts
      public :: ptr_bufferdispls
      public :: ptr_subdomcounts	! be careful with this one
      public :: ptr_subdomdispls	! be careful with this one too

    type SubdomainDistributor
      private
      integer :: ni		! (global) i-size of all core subdomains
      integer :: nj		! (global) j-size of all core subdomains
      integer :: iloc,ilen	! i-deminsion of local core subdomain
      integer :: jloc,jlen	! j-dimension of local core subdomain
      integer :: lsize		! size of the local subdomain only
      integer :: bsize		! size of all subdomains
      integer,pointer,dimension(:) :: imap_bufr	! (bsize)
      integer,pointer,dimension(:) :: jmap_bufr	! (bsize)
      integer,pointer,dimension(:) :: imap_core	! (bsize)
      integer,pointer,dimension(:) :: jmap_core	! (bsize)
      integer,pointer,dimension(:) :: nodata_counts
      integer,pointer,dimension(:) :: nodata_displs
      integer,pointer,dimension(:) :: buffer_counts
      integer,pointer,dimension(:) :: buffer_displs
      integer,pointer,dimension(:) :: subdom_counts
      integer,pointer,dimension(:) :: subdom_displs
    end type SubdomainDistributor

    interface SubdomainDistributor_init; module procedure	&
      allgather_init_,	&
      init_; end interface
    interface init; module procedure	&
      allgather_init_,	&
      init_; end interface
    interface SubdomainDistributor_clean; module procedure	&
      clean_; end interface
    interface clean; module procedure clean_; end interface

    interface get; module procedure get_; end interface
    interface show; module procedure show_; end interface
    interface bufferSize; module procedure bufferSize_; end interface
    interface localSize; module procedure localSize_; end interface

    interface Subdomain_pack; module procedure pack3dr_; end interface
    interface Subdomain_unpack; module procedure unpack3dr_; end interface

    interface ptr_buffercounts; module procedure	&
    	ptr_buffercounts_; end interface
    interface ptr_bufferdispls; module procedure	&
    	ptr_bufferdispls_; end interface
    interface ptr_subdomcounts; module procedure	&
    	ptr_subdomcounts_; end interface
    interface ptr_subdomdispls; module procedure	&
    	ptr_subdomdispls_; end interface

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SubdomainDistributor'

  integer,parameter :: IHALO=1
  integer,parameter :: JHALO=1
#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgather_init_ - allgather subdomain info. before init_()
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgather_init_(ob, ni,iloc,ilen,iperiodic,	&
				   nj,jloc,jlen,jperiodic, comm	)
      use m_mpif90,only : MP_comm_size,MP_type
      use m_mpout ,only : mpout_log
      use m_die,only : MP_die
      implicit none
      type(SubdomainDistributor),intent(out) :: ob

      		! i refers to the first index of a level.

      integer,intent(in) :: ni
      integer,intent(in) :: iloc,ilen
      logical,intent(in) :: iperiodic

      		! j refers to the second index of a level.

      integer,intent(in) :: nj
      integer,intent(in) :: jloc,jlen
      logical,intent(in) :: jperiodic

      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgather_init_'
  integer,allocatable,dimension(:) :: ilocs,ilens
  integer,allocatable,dimension(:) :: jlocs,jlens
  integer :: ier,mtype,nPEs
_ALLENTRY_

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

	allocate( ilocs(0:nPEs-1),ilens(0:nPEs-1),	&
		  jlocs(0:nPEs-1),jlens(0:nPEs-1))

  mType=MP_type(iloc)
  call MPI_allgather(iloc,1,mtype,ilocs,1,mtype,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_allgather(iloc)',ier)
  call MPI_allgather(ilen,1,mtype,ilens,1,mtype,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_allgather(ilen)',ier)
  call MPI_allgather(jloc,1,mtype,jlocs,1,mtype,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_allgather(jloc)',ier)
  call MPI_allgather(jlen,1,mtype,jlens,1,mtype,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_allgather(jlen)',ier)

  call init_(ob,ni,ilocs,ilens,iperiodic,	&
  		nj,jlocs,jlens,jperiodic, comm	)

	deallocate(ilocs,ilens,jlocs,jlens)
_ALLEXIT_
end subroutine allgather_init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - intialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,ni,iloc,ilen,iperiodic,		&
			nj,jloc,jlen,jperiodic,comm	)
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_comm_rank
      use m_mpout ,only : mpout_log
      use m_die,only : assert_,MP_die
      implicit none
      type(SubdomainDistributor),intent(out) :: ob

      integer,intent(in) :: ni
      integer,dimension(0:),intent(in) :: iloc,ilen
      logical,intent(in) :: iperiodic

      integer,intent(in) :: nj
      integer,dimension(0:),intent(in) :: jloc,jlen
      logical,intent(in) :: jperiodic

      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nPEs,myPE
  integer :: iPE
  integer :: ier
  integer :: i,j,l
  integer :: im,jm
_ALLENTRY_

  call MP_comm_size(comm,nPEs,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  	ASSERT(size(iloc)==nPEs)
  	ASSERT(size(jloc)==nPEs)
  	ASSERT(size(ilen)==nPEs)
  	ASSERT(size(jlen)==nPEs)

  ob%ni=ni
  ob%nj=nj
  ob%iloc=iloc(myPE)
  ob%ilen=ilen(myPE)
  ob%jloc=jloc(myPE)
  ob%jlen=jlen(myPE)
  ob%lsize=(ob%ilen+2*IHALO)*(ob%jlen+2*JHALO)

	allocate( ob%buffer_counts(0:nPEs-1),	&
		  ob%buffer_displs(0:nPEs-1),	&
		  ob%subdom_counts(0:nPEs-1),	&
		  ob%subdom_displs(0:nPEs-1),	&
		  ob%nodata_counts(0:nPEs-1),	&
		  ob%nodata_displs(0:nPEs-1)	)

  ob%nodata_counts(:)=0
  ob%nodata_displs(:)=0
  ob%subdom_counts(:)=-1
  ob%subdom_displs(:)=-1

  ob%buffer_counts(0)=(ilen(0)+2*IHALO) * (jlen(0)+2*JHALO)
  ob%bsize=ob%buffer_counts(0)
  ob%buffer_displs(0)=0
  do iPE=1,nPEs-1
    ob%buffer_displs(iPE)=ob%bsize
    ob%buffer_counts(iPE)=(ilen(iPE)+2*IHALO) * (jlen(iPE)+2*JHALO)

    		! Subdomains in the buffer are ordered according to
		! their process IDs, thus the calculation of the
		! %buffer_displs as a simple orderly accumulation.

    ob%bsize=ob%bsize+ob%buffer_counts(iPE)
  end do

	allocate( ob%imap_bufr(1:ob%bsize),	&
		  ob%jmap_bufr(1:ob%bsize),	&
		  ob%imap_core(1:ob%bsize),	&
		  ob%jmap_core(1:ob%bsize)	)

  do iPE=0,nPEs-1
    l=ob%buffer_displs(iPE)

    do jm=jloc(iPE)-JHALO,jloc(iPE)-1
      j=jm
      if(j<=0) then
	j=j+nj
	if(.not.jperiodic) j=1
      endif

      do im=iloc(iPE)-IHALO,iloc(iPE)-1
        i=im
	if(i<=0) then
	  i=i+ni
	  if(.not.iperiodic) i=1
	endif
	l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      end do

      do i=iloc(iPE),iloc(iPE)+ilen(iPE)-1
        l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      end do

      do im=iloc(iPE)+ilen(iPE)-1+1,iloc(iPE)+ilen(iPE)-1+IHALO
        i=im
	if(i>ni) then
	  i=i-ni
	  if(.not.iperiodic) i=ni
	endif
	l=l+1
	ob%imap_bufr(l)=i
	ob%jmap_bufr(l)=j
	ob%imap_core(l)=-1
	ob%jmap_core(l)=-1
      enddo
    enddo

    do j=jloc(iPE),jloc(iPE)+jlen(iPE)-1
      do im=iloc(iPE)-IHALO,iloc(iPE)-1
        i=im
	if(i<=0) then
	  i=i+ni
	  if(.not.iperiodic) i=1
	endif
	l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      enddo

      do i=iloc(iPE),iloc(iPE)+ilen(iPE)-1
        l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=i
        ob%jmap_core(l)=j
      end do

      do im=iloc(iPE)+ilen(iPE)-1+1,iloc(iPE)+ilen(iPE)-1+IHALO
        i=im
	if(i>ni) then
	  i=i-ni
	  if(.not.iperiodic) i=ni
	endif
        l=l+1
	ob%imap_bufr(l)=i
	ob%jmap_bufr(l)=j
	ob%imap_core(l)=-1
	ob%jmap_core(l)=-1
      enddo
    end do

    do jm=jloc(iPE)+jlen(iPE)-1+1,jloc(iPE)+jlen(iPE)-1+JHALO
      j=jm
      if(j>nj) then
        j=j-nj
        if(.not.jperiodic) j=nj
      endif

      do im=iloc(iPE)-IHALO,iloc(iPE)-1
        i=im
	if(i<=0) then
	  i=i+ni
	  if(.not.iperiodic) i=1
	endif
	l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      enddo

      do i=iloc(iPE),iloc(iPE)+ilen(iPE)-1
        l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      end do

      do im=iloc(iPE)+ilen(iPE)-1+1,iloc(iPE)+ilen(iPE)-1+IHALO
        i=im
	if(i>ni) then
	  i=i-ni
	  if(.not.iperiodic) i=ni
	endif
        l=l+1
        ob%imap_bufr(l)=i
        ob%jmap_bufr(l)=j
        ob%imap_core(l)=-1
        ob%jmap_core(l)=-1
      enddo
    end do
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
      implicit none
      type(SubdomainDistributor),intent(inout) :: ob

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_

  ob%ni=0
  ob%nj=0
  ob%iloc=0
  ob%ilen=0
  ob%jloc=0
  ob%jlen=0
  ob%lsize=0
  ob%bsize=0

  deallocate(ob%imap_bufr,ob%imap_core)
  deallocate(ob%jmap_bufr,ob%jmap_core)
  deallocate(ob%buffer_counts,ob%buffer_displs)
  deallocate(ob%nodata_counts,ob%nodata_displs)
  deallocate(ob%subdom_counts,ob%subdom_displs)
_EXIT_
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get configurations
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(ob,itotalSize, jtotalSize,	&
		       ilocalSize, jlocalSize,	&
		       iLoc, iLen, jLoc, jLen	)
      use m_mpout,only : mpout_log
      implicit none
      type(SubdomainDistributor),intent(in) :: ob

      integer,optional,intent(out) :: itotalSize
      integer,optional,intent(out) :: jtotalSize
      integer,optional,intent(out) :: ilocalSize
      integer,optional,intent(out) :: jlocalSize
      integer,optional,intent(out) :: iLoc, iLen
      integer,optional,intent(out) :: jLoc, jLen

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
_ENTRY_

  if(present(itotalSize)) itotalSize=ob%ni
  if(present(jtotalSize)) jtotalSize=ob%nj

  if(present(ilocalSize)) ilocalSize=ob%ilen+2*IHALO
  if(present(jlocalSize)) jlocalSize=ob%jlen+2*JHALO

  if(present(iLen)) iLen=ob%ilen
  if(present(jLen)) jLen=ob%jlen

  if(present(iLoc)) iloc=ob%iloc
  if(present(jLoc)) jloc=ob%jloc
_EXIT_
end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: pack3dr_ - pack3dr a 2d level to a 1d subdomain buffer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine pack3dr_(ob,grid,bufr)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(SubdomainDistributor),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: grid
      real,dimension(:,:),intent(out) :: bufr

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::pack3dr_'
  integer :: i,j,k,l
_ENTRY_

	ASSERT(size(grid,1)==ob%ni)
	ASSERT(size(grid,2)==ob%nj)
	ASSERT(size(bufr,1)==ob%bsize)
  	ASSERT(size(grid,3)==size(bufr,2))

  do k=1,size(grid,3)
    do l=1,ob%bsize
      i=ob%imap_bufr(l)
      j=ob%jmap_bufr(l)
      bufr(l,k)=grid(i,j,k)
    end do
  end do
_EXIT_
end subroutine pack3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpack3dr_ - unpack a buffer to interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpack3dr_(ob,bufr,grid)
      use m_mpout ,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(SubdomainDistributor),intent(in) :: ob
      real,dimension(:,:),intent(in) :: bufr
      real,dimension(:,:,:),intent(out) :: grid

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpack3dr_'
  integer :: i,j,k,l
_ENTRY_

  	ASSERT(size(grid,3)==size(bufr,2))

  do k=1,size(grid,3)
    do l=1,ob%bsize
      i=ob%imap_core(l)
      j=ob%jmap_core(l)
      if(i/=-1.and.j/=-1) grid(i,j,k)=bufr(l,k)
    end do
  end do

_EXIT_
end subroutine unpack3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_buffercounts_ - conditional reference to buffer_counts
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_buffercounts_(ob,intlev,lindex)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : get
      use m_mpout ,only : mpout_log
      implicit none
      type(subdomainDistributor),intent(in) :: ob
      type(interleavedObject)   ,intent(in) :: intlev
      integer,intent(in) :: lindex
      integer,pointer,dimension(:) :: ptr_buffercounts_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_buffercounts_'
  integer :: lcount
_ENTRY_

  call get(intlev,lindex,lcount=lcount)

  select case(lcount)
  case(0)
    ptr_buffercounts_ => ob%nodata_counts(:)
  case default
    ptr_buffercounts_ => ob%buffer_counts(:)
  end select
_EXIT_
end function ptr_buffercounts_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_bufferdispls_ - conditional reference to buffer_displs
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_bufferdispls_(ob,intlev,lindex)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : get
      use m_mpout ,only : mpout_log
      implicit none
      type(subdomainDistributor),intent(in) :: ob
      type(interleavedObject)   ,intent(in) :: intlev
      integer,intent(in) :: lindex
      integer,pointer,dimension(:) :: ptr_bufferdispls_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_bufferdispls_'
  integer :: lcount
_ENTRY_

  call get(intlev,lindex,lcount=lcount)

  select case(lcount)
  case(0)
    ptr_bufferdispls_ => ob%nodata_displs(:)
  case default
    ptr_bufferdispls_ => ob%buffer_displs(:)
  end select
_EXIT_
end function ptr_bufferdispls_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_subdomcounts_ - conditional reference to subdom_counts
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_subdomcounts_(ob,intlev,lindex)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : get,ptr_where
      use m_mpout ,only : mpout_log
      implicit none
      type(subdomainDistributor),intent(in) :: ob
      type(interleavedObject)   ,intent(in) :: intlev
      integer,intent(in) :: lindex
      integer,pointer,dimension(:) :: ptr_subdomcounts_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_subdomcounts_'
  integer :: lbnd,ubnd
  integer :: lsize
  integer :: i
  integer,pointer,dimension(:) :: iPEs
_ENTRY_

  call get(intlev,lindex,lbound=lbnd,ubound=ubnd)
  iPEs => ptr_where(intlev)

  ptr_subdomcounts_ => ob%subdom_counts(:)

  lsize=ob%lsize

	! statements below seem inconsistent with the intent(in)
	! attribute of ob, since ob%subdomcounts is referenced by
	! ptr_subdomcounts_.  However, let's see what is going to
	! happen.
  ptr_subdomcounts_(:)=0
  do i=lbnd,ubnd
    ptr_subdomcounts_(iPEs(i)+1) = lsize
  end do
  nullify(iPEs)
_EXIT_
end function ptr_subdomcounts_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_subdomdispls_ - conditional reference to subdom_displs
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_subdomdispls_(ob,intlev,lindex)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : get,ptr_where
      use m_mpout ,only : mpout_log
      implicit none
      type(subdomainDistributor),intent(in) :: ob
      type(interleavedObject)   ,intent(in) :: intlev
      integer,intent(in) :: lindex
      integer,pointer,dimension(:) :: ptr_subdomdispls_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_subdomdispls_'
  integer :: lbnd,ubnd
  integer :: lsize
  integer,pointer,dimension(:) :: iPEs
  integer :: i
_ENTRY_

  call get(intlev,lindex,lbound=lbnd,ubound=ubnd)
  iPEs => ptr_where(intlev)

  ptr_subdomdispls_ => ob%subdom_displs(:)

  lsize=ob%lsize

	! statements below seem inconsistent with the intent(in)
	! attribute of ob, since ob%subdomdispls is referenced by
	! ptr_subdomdispls_.  However, let's see what is going to
	! happen.
  ptr_subdomdispls_(:)=0
  do i=lbnd,ubnd
    ptr_subdomdispls_(iPEs(i)+1) = (i-lbnd)*lsize
  end do
_EXIT_
end function ptr_subdomdispls_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bufferSize_ - size of the buffer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function bufferSize_(ob)
      use m_mpout ,only : mpout_log
      implicit none
      type(SubdomainDistributor),intent(in) :: ob
      integer :: bufferSize_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bufferSize_'
_ENTRY_
  bufferSize_=ob%bsize
_EXIT_
end function bufferSize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: localSize_ - size of the local subdomain
!
! !DESCRIPTION:
!
! !INTERFACE:

    function localSize_(ob)
      use m_mpout ,only : mpout_log
      implicit none
      type(SubdomainDistributor),intent(in) :: ob
      integer :: localSize_

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::localSize_'
_ENTRY_
  localSize_=ob%lsize
_EXIT_
end function localSize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: show_ - show size information
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine show_(ob,where,lu)
      use m_mpout,only : mpout	! default output unit
      use m_mpout,only : mpout_log
      implicit none
      type(SubdomainDistributor),intent(in) :: ob
      character(len=*),intent(in) :: where
      integer,optional,intent(in) :: lu

! !REVISION HISTORY:
! 	16Feb05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::show_'
  integer :: lu_

  integer :: itotalSize,ilocalSize,iLoc,iLen
  integer :: jtotalSize,jlocalSize,jLoc,jLen
_ENTRY_

  lu_=mpout
  if(present(lu)) lu_=lu

  call get(ob,	&
  	itotalsize=itotalsize, jtotalsize =jtotalsize,	&
  	ilocalsize=ilocalsize, jlocalsize =jlocalsize,	&
	iLoc=iLoc, iLen=iLen,  jLoc=jLoc, jLen=jLen	)

  write(lu_,'(2a,i5)') trim(where),': itotalSize =',itotalSize
  write(lu_,'(2a,i5)') trim(where),': jtotalSize =',jtotalSize
  write(lu_,'(2a,i5)') trim(where),': ilocalSize =',ilocalSize
  write(lu_,'(2a,i5)') trim(where),': jlocalSize =',jlocalSize
  write(lu_,'(a,2(a,i3))') trim(where),': iLoc, iLen =',iLoc,',',iLen
  write(lu_,'(a,2(a,i3))') trim(where),': jLoc, jLen =',jLoc,',',jLen
_EXIT_
end subroutine show_

end module m_SubdomainDistributor
