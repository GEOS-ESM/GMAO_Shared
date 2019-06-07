
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SubdomainDistributorComm - message communication ops.
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_SubdomainDistributorComm
      implicit none
      private	! except

      public :: distribute		! interleaved -> subdomains
      public :: undistribute		! subdomains -> interleaved

    interface distribute; module procedure	&
    	distr2dr_,	&
    	distr3dr_; end interface
    interface undistribute; module procedure	&
    	undistr2dr_,	&
    	undistr3dr_; end interface

! !REVISION HISTORY:
! 	20Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SubdomainDistributorComm'
#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: distr3dr_ - distribute an interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine distr3dr_(intlev,subdom,send,recv,comm)
      use m_interleavedObject,only : interleavedObject,get
      use m_interleavedObject,only : totalSize,localSize,maximSize
      use m_SubdomainDistributor,only : subdomainDistributor
      use m_SubdomainDistributor,only : get,Subdomain_pack
      use m_SubdomainDistributor,only : bufferSize,localSize
      use m_SubdomainDistributor,only : ptr_buffercounts
      use m_SubdomainDistributor,only : ptr_bufferdispls
      use m_SubdomainDistributor,only : ptr_subdomcounts
      use m_SubdomainDistributor,only : ptr_subdomdispls
      use m_mpif90,only : MP_type
      use m_mpout,only : mpout_log
      use m_die,only : die,MP_die,assert_
      implicit none
      type(interleavedObject)   ,intent(in) :: intlev
      type(SubdomainDistributor),intent(in) :: subdom
      real,dimension(:,:,:) ,intent(in) :: send		! interleaved
      real,dimension(:,:,:) ,intent(out) :: recv	! subdomains
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::distr3dr_'
  integer :: ni,nj,mlev,nlev,bsize,lsize
  real,allocatable,dimension(:,:) :: bufr
  integer :: isend,lsend
  integer :: lcount,lbnd,ubnd
  integer :: mtype
  integer :: ier
  integer :: nrij
_ALLENTRY_

  bsize=bufferSize(subdom)	! send buffer size
  lsize=localSize(subdom)	! receive buffer block size
  nlev=totalSize(intlev)

  nrij=size(recv,1)*size(recv,2)	! actual receive buffer
  ASSERT(nrij==lsize)
  ASSERT(size(recv,3)==nlev)

  call get(subdom,itotalSize=ni,jtotalSize=nj)
  mlev=localSize(intlev)

  ASSERT(size(send,1)==ni)
  ASSERT(size(send,2)==nj)
  ASSERT(size(send,3)==mlev)

  	allocate(bufr(bsize,1))

  do isend=1,maximSize(intlev)	! for every local level
    call get(intlev,isend,lcount=lcount,lbound=lbnd,ubound=ubnd)

	! Subdomain_pack() may operates on zero-size arrays (lcount=0).
    lsend=isend+lcount-1
    call Subdomain_pack(subdom,send(:,:,isend:lsend),bufr(:,1:lcount))

    mtype=MP_type(bufr)
    call MPI_alltoallv(			&
    		bufr(:,1:lcount),	&
    		ptr_buffercounts(subdom,intlev,isend),	&
    		ptr_bufferdispls(subdom,intlev,isend), mtype,	&
    		recv(:,:,lbnd:ubnd),	&
    		ptr_subdomcounts(subdom,intlev,isend),	&
    		ptr_subdomdispls(subdom,intlev,isend), mtype,	&
		comm,ier)
	if(ier/=0) call MP_die(myname_,'MP_alltoallv()',ier)
  end do

  	deallocate(bufr)

_ALLEXIT_
end subroutine distr3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: distr2dr_ - distribute an interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine distr2dr_(intlev,subdom,send,recv,comm)
      use m_interleavedObject,only : interleavedObject,get
      use m_interleavedObject,only : totalSize,localSize,maximSize
      use m_SubdomainDistributor,only : subdomainDistributor
      use m_SubdomainDistributor,only : get,Subdomain_pack
      use m_SubdomainDistributor,only : bufferSize,localSize
      use m_SubdomainDistributor,only : ptr_buffercounts
      use m_SubdomainDistributor,only : ptr_bufferdispls
      use m_SubdomainDistributor,only : ptr_subdomcounts
      use m_SubdomainDistributor,only : ptr_subdomdispls
      use m_mpif90,only : MP_type
      use m_mpout,only : mpout_log
      use m_die,only : die,MP_die,assert_
      implicit none
      type(interleavedObject)   ,intent(in) :: intlev
      type(SubdomainDistributor),intent(in) :: subdom
      real,dimension(:,:,:) ,intent(in) :: send		! interleaved
      real,dimension(:,:  ) ,intent(out) :: recv	! subdomains
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::distr2dr_'
  integer :: ni,nj,mlev,nlev,bsize,lsize
  real,allocatable,dimension(:,:) :: bufr
  integer :: isend,lsend
  integer :: lcount,lbnd,ubnd
  integer :: mtype
  integer :: ier
  integer :: nrij
_ALLENTRY_

  bsize=bufferSize(subdom)	! send buffer size
  lsize=localSize(subdom)	! receive buffer block size
  nlev=totalSize(intlev)
  	ASSERT(nlev==1)

  nrij=size(recv,1)*size(recv,2)
  ASSERT(nrij==lsize)

  call get(subdom,itotalSize=ni,jtotalSize=nj)
  mlev=localSize(intlev)

  ASSERT(size(send,1)==ni)
  ASSERT(size(send,2)==nj)
  ASSERT(size(send,3)==mlev)

  	allocate(bufr(bsize,1))

  	! there is only one (or zero) local level.

  isend=1
  call get(intlev,isend,lcount=lcount,lbound=lbnd,ubound=ubnd)
  	ASSERT(lbnd==1)
  	ASSERT(ubnd==1)

	! Subdomain_pack() may operates on zero-size arrays (lcount=0).
  lsend=isend+lcount-1
  call Subdomain_pack(subdom,send(:,:,isend:lsend),bufr(:,1:lcount))

  mtype=MP_type(bufr)
  call MPI_alltoallv(			&
    		bufr(:,1:lcount),	&
    		ptr_buffercounts(subdom,intlev,isend),	&
    		ptr_bufferdispls(subdom,intlev,isend), mtype,	&
    		recv(:,:),	&
    		ptr_subdomcounts(subdom,intlev,isend),	&
    		ptr_subdomdispls(subdom,intlev,isend), mtype,	&
		comm,ier)
	if(ier/=0) call MP_die(myname_,'MP_alltoallv()',ier)

  	deallocate(bufr)

_ALLEXIT_
end subroutine distr2dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: undistr3dr_ - un-distribute subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine undistr3dr_(intlev,subdom,send,recv,comm)
      use m_interleavedObject,only : interleavedObject,get
      use m_interleavedObject,only : totalSize,localSize,maximSize
      use m_SubdomainDistributor,only : subdomainDistributor
      use m_SubdomainDistributor,only : get,Subdomain_unpack
      use m_SubdomainDistributor,only : bufferSize,localSIze
      use m_SubdomainDistributor,only : ptr_buffercounts
      use m_SubdomainDistributor,only : ptr_bufferdispls
      use m_SubdomainDistributor,only : ptr_subdomcounts
      use m_SubdomainDistributor,only : ptr_subdomdispls
      use m_mpif90,only : MP_type
      use m_mpout,only : mpout_log
      use m_die,only : die,MP_die,assert_
      implicit none
      type(interleavedObject)   ,intent(in) :: intlev
      type(SubdomainDistributor),intent(in) :: subdom
      real,dimension(:,:,:) ,intent(in) :: send		! subdomains
      real,dimension(:,:,:) ,intent(out) :: recv	! interleaved
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::undistr3dr_'
  integer :: ni,nj,mlev,nlev,bsize,lsize
  real,allocatable,dimension(:,:) :: bufr
  integer :: isend,lsend
  integer :: lcount,lbnd,ubnd
  integer :: mtype
  integer :: ier
  integer :: nsij
_ALLENTRY_

  bsize=bufferSize(subdom)
  lsize=localSize(subdom)
  nlev=totalSize(intlev)

  nsij=size(send,1)*size(send,2)
  ASSERT(nsij==lsize)
  ASSERT(size(send,3)==nlev)

  call get(subdom,itotalSize=ni,jtotalSize=nj)
  mlev=localSize(intlev)

  ASSERT(size(recv,1)==ni)
  ASSERT(size(recv,2)==nj)
  ASSERT(size(recv,3)==mlev)

  	allocate(bufr(bsize,1))

  do isend=1,maximSize(intlev)	! for every local level
    call get(intlev,isend,lcount=lcount,lbound=lbnd,ubound=ubnd)

	! Subdomain_unpack() may operates on zero-size arrays (lcount=0).
    lsend=isend+lcount-1

    mtype=MP_type(bufr)
    call MPI_alltoallv(			&
    		send(:,:,lbnd:ubnd),	&
    		ptr_subdomcounts(subdom,intlev,isend),	&
    		ptr_subdomdispls(subdom,intlev,isend), mtype,	&
    		bufr(:,1:lcount),	&
    		ptr_buffercounts(subdom,intlev,isend),	&
    		ptr_bufferdispls(subdom,intlev,isend), mtype,	&
		comm,ier)
	if(ier/=0) call MP_die(myname_,'MP_alltoallv()',ier)

    call Subdomain_unpack(subdom,bufr(:,1:lcount),recv(:,:,isend:lsend))
  end do

  	deallocate(bufr)

_ALLEXIT_
end subroutine undistr3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: undistr2dr_ - un-distribute subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine undistr2dr_(intlev,subdom,send,recv,comm)
      use m_interleavedObject,only : interleavedObject,get
      use m_interleavedObject,only : totalSize,localSize,maximSize
      use m_SubdomainDistributor,only : subdomainDistributor
      use m_SubdomainDistributor,only : get,Subdomain_unpack
      use m_SubdomainDistributor,only : bufferSize,localSize
      use m_SubdomainDistributor,only : ptr_buffercounts
      use m_SubdomainDistributor,only : ptr_bufferdispls
      use m_SubdomainDistributor,only : ptr_subdomcounts
      use m_SubdomainDistributor,only : ptr_subdomdispls
      use m_mpif90,only : MP_type
      use m_mpout,only : mpout_log
      use m_die,only : die,MP_die,assert_
      implicit none
      type(interleavedObject)   ,intent(in) :: intlev
      type(SubdomainDistributor),intent(in) :: subdom
      real,dimension(:,:  ) ,intent(in) :: send		! subdomains
      real,dimension(:,:,:) ,intent(out) :: recv	! interleaved
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::undistr2dr_'
  integer :: ni,nj,mlev,nlev,bsize,lsize
  real,allocatable,dimension(:,:) :: bufr
  integer :: isend,lsend
  integer :: lcount,lbnd,ubnd
  integer :: mtype
  integer :: ier
  integer :: nsij
_ALLENTRY_

  bsize=bufferSize(subdom)
  lsize=localSize(subdom)
  nlev=totalSize(intlev)
  	ASSERT(nlev==1)

  nsij=size(send,1)*size(send,2)
  ASSERT(nsij==lsize)

  call get(subdom,itotalSize=ni,jtotalSize=nj)
  mlev=localSize(intlev)

  ASSERT(size(recv,1)==ni)
  ASSERT(size(recv,2)==nj)
  ASSERT(size(recv,3)==mlev)

  	allocate(bufr(bsize,1))

  	! There is only one level to send

  isend=1
  call get(intlev,isend,lcount=lcount,lbound=lbnd,ubound=ubnd)

	! Subdomain_unpack() may operates on zero-size arrays (lcount=0).
  lsend=isend+lcount-1

  mtype=MP_type(bufr)
  call MPI_alltoallv(			&
    		send(:,:),	&
    		ptr_subdomcounts(subdom,intlev,isend),	&
    		ptr_subdomdispls(subdom,intlev,isend), mtype,	&
    		bufr(:,1:lcount),	&
    		ptr_buffercounts(subdom,intlev,isend),	&
    		ptr_bufferdispls(subdom,intlev,isend), mtype,	&
		comm,ier)
	if(ier/=0) call MP_die(myname_,'MP_alltoallv()',ier)

  call Subdomain_unpack(subdom,bufr(:,1:lcount),recv(:,:,isend:lsend))

  	deallocate(bufr)

_ALLEXIT_
end subroutine undistr2dr_
end module m_SubdomainDistributorComm
