!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_simpleGrid - simple partitioned grid
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_simpleGrid
      implicit none
      private	! except

      public :: simpleGrid	! data structure
      public :: simpleGrid_set	! constructor
      public :: simpleGrid_get	! inquirer
      public :: simpleGrid_clean, clean ! destructor

    type simpleGrid
      private
      integer :: ndim=-1
      integer,pointer,dimension(:) :: counts
      integer,pointer,dimension(:) :: displs
    end type simpleGrid

    interface simpleGrid_set; module procedure &
      init1_, &
      init2_, &
      init0_, &
      initn_; end interface

    interface simpleGrid_get  ; module procedure   get_; end interface
    interface simpleGrid_clean; module procedure clean_; end interface
    interface            clean; module procedure clean_; end interface

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_simpleGrid'
!!
!! 1-d, init1_()
! call simpleGrid_set(nlon,comm,icount,idispl)
! iloc=idispl+1
!!
!! 2-d, init2_()
! call simpleGrid_set(nlon,nlat,comm,icount,idispl,jcount,jdispl)
! iloc=idispl+1
! jloc=jdispl+1
!!
!! 2-d, initn_()
! call simpleGrid_set(ob,(/nlon,nlat/),comm)
! call simpleGrid_get(ob,1,count=icount,displ=idispl)
! iloc=idispl+1
! call simpleGrid_get(ob,2,count=jcount,displ=jdispl)
! jloc=jdispl+1
!!

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init1_ - partition a rank 1 grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init1_(lsize,comm,count,displ)
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_die,only : MP_die
      implicit none
      integer,intent(in) :: lsize
      integer,intent(in) :: comm
      integer,intent(out) :: count
      integer,intent(out) :: displ

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init1_'
  integer :: nPEs,myPE,ier

  call MP_comm_size(comm,nPEs,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  call MP_comm_rank(comm,myPE,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
  call SimplePartition_(lsize,nPEs,myPE,count=count,displ=displ)
end subroutine init1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2_ - partition a rank 2 grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2_(isize,jsize,comm,icount,jcount,idispl,jdispl)
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_die,only : MP_die
      implicit none
      integer,intent(in) :: isize,jsize
      integer,intent(in) :: comm
      integer,intent(out) :: icount,jcount
      integer,intent(out) :: idispl,jdispl

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2_'
  integer :: nPEs,myPE,ier
  integer,dimension(2) :: lsizes,counts,displs

  lsizes(1)=isize
  lsizes(2)=jsize
  call init0_(lsizes,comm,counts,displs)
  icount=counts(1);idispl=displs(1)
  jcount=counts(2);jdispl=displs(2)
end subroutine init2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initn_ - partition a rank n grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initn_(ob,lsizes,comm)
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_die,only : MP_die
      implicit none
      type(simpleGrid),intent(out) :: ob
      integer,dimension(:),intent(in) :: lsizes
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initn_'

  ob%ndim=size(lsizes)
  allocate(ob%counts(ob%ndim))
  allocate(ob%displs(ob%ndim))
  call init0_(lsizes,comm,ob%counts,ob%displs)
end subroutine initn_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_die,only : MP_die
      implicit none
      type(simpleGrid),intent(out) :: ob

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  deallocate(ob%counts)
  deallocate(ob%displs)
  ob%ndim=-1
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get partition a rank n grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(ob,idim,count,displ,ndim)
      use m_die,only : assert_
      implicit none
      type(simpleGrid),intent(out) :: ob
      integer,intent(in) :: idim
      integer,intent(out) :: count
      integer,intent(out) :: displ
      integer,optional,intent(out) :: ndim

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: ndim_

  ndim_=ob%ndim
    ASSERT(idim> 0)
    ASSERT(idim<=ob%ndim)

  count=ob%counts(idim)
  displ=ob%displs(idim)
  if(present(ndim)) ndim=ndim_
end subroutine get_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init0_ - partition a rank-n grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init0_(sizes,comm,counts,displs)
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_die,only : MP_die,die,assert_
      implicit none
      integer,dimension(:),intent(in) :: sizes
      integer,intent(in) :: comm
      integer,dimension(:),intent(out) :: counts
      integer,dimension(:),intent(out) :: displs

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init0_'
  integer :: nPEs,myPE,ier
  integer :: lPi,nPi,mPi ! total, size, and rank of the i-th subspace.
  integer,dimension(0:size(sizes)-1) :: ldims
  integer :: ndim,i
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  call MP_comm_rank(comm,myPE,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  ndim=size(ldims)
    ASSERT(ndim==size(counts))
    ASSERT(ndim==size(displs))

  ldims(:)=0
  call MP_dims_create(nPEs,ndim,ldims,ier)
    if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  lPi=myPE ! total PE counts in remaining dimensions
  do i=0,ndim-1
    nPi=ldims(i)

    mPi=mod(lPi,nPi) ! 
    lPi=    lPi/nPi
    call SimplePartition_(sizes(i),nPi,mPi,counts(i),displs(i))
  end do
end subroutine init0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: simplePartition_ - simple partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine simplePartition_(ngrid,nproc,iproc,count,displ)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: ngrid	! number of total grid points
      integer,intent(in) :: nproc	! number of PEs
      integer,intent(in) :: iproc	! PE ID
      integer,intent(out) :: count	! number of local grid points
      integer,intent(out) :: displ	! displacement of the local grid

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::simplePartition_'
  integer :: resid

        ASSERT(nproc> 0)
        ASSERT(iproc>=0 .and. iproc< nproc)

  resid=mod(ngrid,nproc)

  count=ngrid/nproc
  if(iproc< resid) count=count+1

  displ=count*iproc
  if(iproc>=resid) displ=displ+resid
end subroutine simplePartition_
end module m_simpleGrid
