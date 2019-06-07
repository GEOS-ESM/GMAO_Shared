!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_parGrADS - a parallel GrADS file reader or writer.
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_parGrADS
      use m_GrADS, only : GrADS_reader
      use m_GrADS, only : GrADS_writer
      implicit none
      private	! except

      public :: parGrADS		! data structure
      public :: parGrADS_create         ! open a file for write
      public :: parGrADS_write          ! output a gridded variable
      public :: parGrADS_open           ! open a file for read
      public :: parGrADS_read           ! input a gridded variable
      public :: parGrADS_close,close    ! close a file in either form

    type parGrADS
      private
      integer :: root = -1
      integer :: myPE = -1
      logical :: isroot = .false.

      type(GrADS_reader),pointer :: rf => null()
      type(GrADS_writer),pointer :: wf => null()
    end type parGrADS

    interface parGrADS_create; module procedure &
      wopen_; end interface

    interface parGrADS_write; module procedure &
      !write2dlev_, &
      !write3dlev_, &
      write2dsub_, &
      write3dsub_; end interface

    interface parGrADS_open; module procedure &
      ropen_; end interface

    interface parGrADS_read; module procedure &
      !read2dlev_, &
      !read3dlev_, &
      read2dsub_, &
      read3dsub_; end interface

    interface parGrADS_close; module procedure close_; end interface
    interface          close; module procedure close_; end interface

! !REVISION HISTORY:
! 	31May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_parGrADS'
#include "assert.H"

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: wopen_ - open for output
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine wopen_(fd,ctrl,nlon,nlat,nlev,nvar,comm, &
      root, nymd,nhms,incr,dset,udef,unit,xdef,ydef,zdef)
      use m_mpif90,only : MP_comm_rank
      use m_GrADS,only : GrADS_create
      use m_GrADS,only : bcast
      use m_die,only : MP_die,die,assert_
      implicit none
      type(parGrADS),intent(inout) :: fd
      character(len=*),intent(in) :: ctrl
      integer,intent(in) :: nlon
      integer,intent(in) :: nlat
      integer,intent(in) :: nlev
      integer,intent(in) :: nvar
      integer,intent(in) :: comm
      integer,optional,intent(in) :: root
      integer,optional,intent(in) :: nymd
      integer,optional,intent(in) :: nhms
      integer,optional,intent(in) :: incr
      character(len=*),optional,intent(in) :: dset
      real   ,optional,intent(in) :: udef
      integer,optional,intent(in) :: unit ! speicify the logical unit
      real,dimension(:),optional,intent(in) :: xdef
      real,dimension(:),optional,intent(in) :: ydef
      real,dimension(:),optional,intent(in) :: zdef

! !REVISION HISTORY:
! 	31May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::wopen_'
  real,allocatable,dimension(:) :: xdef_,ydef_,zdef_
  integer :: ier,k

  fd%root=0
  if(present(root)) fd%root=root

  if(present(xdef)) call die(myname_,'yet to be implemented "xdef="')
  if(present(ydef)) call die(myname_,'yet to be implemented "ydef="')
  
  call MP_comm_rank(comm,fd%myPE,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  fd%isroot = fd%myPE==fd%root

  allocate(fd%wf)
  if(fd%isroot) then
    allocate(zdef_(nlev))
    if(present(zdef)) then
      ASSERT(nlev==size(zdef))
      zdef_=zdef
    else
      zdef_=(/(1.*k,k=1,nlev)/)
    endif
    call GrADS_create(fd%wf,ctrl,nlon,nlat,nlev,zdef_,nvar,    &
      nymd=nymd,nhms=nhms,incr=incr,dset=dset,udef=udef, &
      unit=unit,stat=ier   )
!      unit=unit,xdef=xdef,ydef=ydef,zdef=zdef,stat=ier   )
      if(ier/=0) call die(myname_, &
        'GrADS_create("'//trim(ctrl)//'")',ier)
    deallocate(zdef_)
  endif
  call bcast(fd%wf,fd%root,comm)
end subroutine wopen_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ropen_ - open for input
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine ropen_(fd,ctrl,comm,root,unit)
      use m_GrADS,only : GrADS_open
      use m_GrADS,only : bcast
      use m_mpif90,only : MP_comm_rank
      use m_die,only : MP_die,die
      implicit none
      type(parGrADS),intent(out) :: fd    ! the object
      character(len=*),intent(in) :: ctrl ! file name
      integer,intent(in) :: comm          ! communicator
      integer,optional,intent(in) :: root ! where I/O will take place.
      integer,optional,intent(in) :: unit ! explicitly required unit
        ! Arguments ctrl and unit are significant only on root

! !REVISION HISTORY:
! 	31May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ropen_'
  integer :: ier

  call MP_comm_rank(comm,fd%myPE,ier)
    if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  fd%root=0
  if(present(root)) fd%root=root
  fd%isroot = fd%myPE==fd%root

  allocate(fd%rf)
  if(fd%isroot) then
    call GrADS_open(fd%rf,ctrl,stat=ier,unit=unit)
      if(ier/=0) call die(myname_,'GrADS_open("'//trim(ctrl)//'")',ier)
  endif
  call bcast(fd%rf,fd%root,comm)
  end subroutine ropen_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: close_ - close an file
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine close_(fd)
      !use m_GrADS,only : GrADS_close
      use m_rGrADS,only : rGrADS_close
      use m_wGrADS,only : wGrADS_close
      use m_die  ,only : die
      implicit none
      type(parGrADS),intent(inout) :: fd

! !REVISION HISTORY:
! 	31May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::close_'
  integer :: ier

  if(associated(fd%wf)) then
    if(fd%isroot) then
      call wGrADS_close(fd%wf,stat=ier)
        if(ier/=0) call die(myname_,'GrADS_close(%wf)',ier)
    endif
    deallocate(fd%wf)
  endif
  if(associated(fd%rf)) then
    if(fd%isroot) then
      call rGrADS_close(fd%rf,stat=ier)
        if(ier/=0) call die(myname_,'GrADS_close(%rf)',ier)
    endif
    deallocate(fd%rf)
  endif
  fd%isroot=.false.
  fd%root  = -1
  fd%myPE  = -1
end subroutine close_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read2dsub_ - read 2-d grid in Subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read2dsub_(fd,vname,grid,v2d,comm, itim,swapij)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : interleavedObject_init
      use m_interleavedObject,only : clean
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_swapij,only : swap => swapij
      use m_GrADS,only : inquire
      use m_GrADS,only : GrADS_read
      use m_die,only : assert_,die

      implicit none
      type(parGrADS),intent(inout) :: fd   ! the input file object
      character(len=*),intent(in) :: vname ! name of the var. to read
      type(SubdomainDistributor),intent(in) :: grid   ! distribution
      real,dimension(:,:),target,intent(out) :: v2d   ! output grid
      integer,intent(in) :: comm ! communicator
      integer,optional,intent(in) :: itim   ! a sepecified time?
      logical,optional,intent(in) :: swapij ! swap i and j?
        ! arguments vname and itim are siginificant only on root

! !REVISION HISTORY:
! 	01Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read2dsub_'

  type(InterleavedObject) :: scat
  integer :: ier
  integer :: im,il,jm,jl,kl,itim_
  integer :: nlon,nlat
  logical :: swapij_
  real,allocatable,target,dimension(:,:,:) :: vinp
  real,pointer    ,dimension(:,:,:) :: send
  real,pointer    ,dimension(:,:  ) :: recv

  swapij_=.false.
  if(present(swapij)) swapij_=swapij

    ! Get dimensions, according to the arguments defining the output
  call get(grid,itotalsize=im,ilocalsize=il, &
		jtotalsize=jm,jlocalsize=jl  )

    ! Get dimensions, according to the source of input
  call inquire(fd%rf,nlon=nlon,nlat=nlat)

  ASSERT(im==nlon)
  ASSERT(jm==nlat)
  if(swapij_) then
    ASSERT(il==size(v2d,2))
    ASSERT(jl==size(v2d,1))
  else
    ASSERT(il==size(v2d,1))
    ASSERT(jl==size(v2d,2))
  endif

    ! Define a root PE input buffer
  kl=0
  if(fd%isroot) kl=1
  allocate(vinp(im,jm,kl))

  if(fd%isroot) then
    itim_=1
    if(present(itim)) itim_=itim
    call GrADS_read(fd%rf,vname,itim_,vinp(:,:,1),stat=ier)
      if(ier/=0) call die(myname_,'GrADS_read("'//trim(vname)//'")',ier)
  endif

    ! Define a single PE input buffer
  call InterleavedObject_init(scat,1,comm,fd%root)
  send => vinp(:,:,:)

  if(swapij_) then
    ! Define the final output through a swapping buffer
    allocate(recv(il,jl))
    call distribute(scat,grid,send,recv,comm)
    deallocate(vinp)
    nullify(send)
    call swap(recv,v2d) ! swap i and j
    deallocate(recv)

  else
    ! Define the final output buffer if no swapping is required.
    recv => v2d
    call distribute(scat,grid,send,recv,comm)
    deallocate(vinp)
    nullify(send)
    nullify(recv)
  endif
end subroutine read2dsub_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read3dsub_ - read 3-d grid in Subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read3dsub_(fd,vname,grid,v3d,comm, itim,swapij)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : interleavedObject_init
      use m_interleavedObject,only : clean
      use m_interleavedComm  ,only : alloc_scatterv
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_swapij,only : swap => swapij
      use m_GrADS,only : inquire
      use m_GrADS,only : GrADS_read
      use m_die,only : assert_,die
      use m_mpout,only : mpout_log

      implicit none
      type(parGrADS),intent(inout) :: fd   ! the input file object
      character(len=*),intent(in) :: vname ! name of the var. to read
      type(SubdomainDistributor),intent(in) :: grid   ! distribution
      real,dimension(:,:,:),target,intent(out) :: v3d ! output grid
      integer,intent(in) :: comm ! communicator
      integer,optional,intent(in) :: itim   ! a sepecified time?
      logical,optional,intent(in) :: swapij ! swap i and j?
        ! arguments vname and itim are siginificant only on root

! !REVISION HISTORY:
! 	01Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read3dsub_'

  type(InterleavedObject) :: scat
  integer :: ier
  integer :: im,il,jm,jl,km,kl,itim_
  integer :: nlon,nlat,nlev
  logical :: swapij_
  real,allocatable,dimension(:,:,:) :: vinp
  real,pointer,dimension(:,:,:) :: send
  real,pointer,dimension(:,:,:) :: recv

  swapij_=.false.
  if(present(swapij)) swapij_=swapij

    ! Get dimensions, according to the arguments defining the output
  call get(grid,itotalsize=im,ilocalsize=il, &
		jtotalsize=jm,jlocalsize=jl  )
  km=size(v3d,3)

    ! Get dimensions, according to the source of input
  call inquire(fd%rf,nlon=nlon,nlat=nlat,nlev=nlev)

  ASSERT(im==nlon)
  ASSERT(jm==nlat)
  call mpout_log(myname_,'km',km)
  call mpout_log(myname_,'nlev',nlev)
  ASSERT(km==nlev)
  if(swapij_) then
    ASSERT(il==size(v3d,2))
    ASSERT(jl==size(v3d,1))
  else
    ASSERT(il==size(v3d,1))
    ASSERT(jl==size(v3d,2))
  endif

    ! Define a single PE input buffer (vinp)
  kl=0
  if(fd%isroot) kl=km
  allocate(vinp(im,jm,kl))
  if(fd%isroot) then
    itim_=1
    if(present(itim)) itim_=itim
    call GrADS_read(fd%rf,vname,itim_,vinp(:,:,:),stat=ier)
      if(ier/=0) call die(myname_,'GrADS_read("'//trim(vname)//'")',ier)
  endif

    ! Define a multiple PEs buffer for interleaved levels (send)
  call InterleavedObject_init(scat,nlev,comm,fd%root)
  call alloc_scatterv(im,jm,scat,vinp,send,fd%root,comm,vname)
  deallocate(vinp)

  if(swapij_) then
    ! Define the final otuput through a swapping buffer as required.
    allocate(recv(il,jl,km))
    call distribute(scat,grid,send,recv,comm)
    deallocate(send)
    call swap(recv,v3d) ! swap i and j
    deallocate(recv)

  else
    ! Define the final output buffer if no swapping is required.
    recv => v3d
    call distribute(scat,grid,send,recv,comm)
    deallocate(send)
    nullify(recv)
  endif
end subroutine read3dsub_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read2dlev_ - read a 2-d grid as interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read2dlev_(fd,vname,scat,v2d,comm, itim)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : inquire
      use m_interleavedObject,only : localsize
      use m_interleavedObject,only : clean
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_GrADS,only : inquire
      use m_GrADS,only : GrADS_read
      use m_die,only : assert_,die

      implicit none
      type(parGrADS),intent(inout) :: fd   ! the input file object
      character(len=*),intent(in) :: vname ! name of the var. to read
      type(interleavedObject),intent(in) :: scat ! distribution
      real,dimension(:,:,:),target,intent(out) :: v2d ! output grid
      integer,intent(in) :: comm ! communicator
      integer,optional,intent(in) :: itim   ! a sepecified time?
        ! arguments vname and itim are siginificant only on root

! !REVISION HISTORY:
! 	01Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read2dlev_'
  integer :: ier
  integer :: im,jm,kl,root
  integer :: nlon,nlat
  integer :: itim_
  
    ! Check the consistancy between the communicators defined by the
    ! output arguments and by the input data file.
  call inquire(scat,root=root)
  ASSERT(root==fd%root)

  itim_=1
  if(present(itim)) itim_=itim

  if(fd%isroot) then
      ! Get the 3rd dimenstion defined by the input arguments
    kl=localsize(scat) ! 0 or 1
    ASSERT(kl==1)
    ASSERT(kl==size(v2d,3))

      ! Get the other dimensions defined by the input data file
    call inquire(fd%rf,nlon=nlon,nlat=nlat)

      ! Get the other dimensions defined by the input arguments
    im=size(v2d,1)
    jm=size(v2d,2)
    ASSERT(im==nlon)
    ASSERT(jm==nlat)

    call GrADS_read(fd%rf,vname,itim_,v2d(:,:,1),stat=ier)
      if(ier/=0) call die(myname_,'GrADS_read("'//trim(vname)//'")',ier)

  else
    ! I don't really care what is going on on other processors.  Well,
    ! almost.  Let us just make sure there is NOTHING to be returned.

    kl=localsize(scat) ! 0 or 1
    ASSERT(kl==0)
    ASSERT(kl==size(v2d,3))
  endif
end subroutine read2dlev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read3dlev_ - read a 3-d grid as interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read3dlev_(fd,vname,scat,v3d,comm, itim)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : inquire
      use m_interleavedObject,only : localsize
      use m_interleavedObject,only : totalsize
      use m_interleavedObject,only : clean
      use m_InterleavedComm,only : scatterv
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_GrADS,only : inquire
      use m_GrADS,only : GrADS_read
      use m_die,only : assert_,die

      implicit none
      type(parGrADS),intent(inout) :: fd   ! the input file object
      character(len=*),intent(in) :: vname ! name of the var. to read
      type(interleavedObject),intent(in) :: scat ! distribution
      real,dimension(:,:,:),target,intent(out) :: v3d ! output grid
      integer,intent(in) :: comm ! communicator
      integer,optional,intent(in) :: itim   ! a sepecified time?
        ! arguments vname and itim are siginificant only on root

! !REVISION HISTORY:
! 	01Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read3dlev_'
  integer :: ier
  integer :: im,jm,km,kl,root
  integer :: nlon,nlat,nlev
  integer :: itim_
  real,allocatable,dimension(:,:,:) :: vinp
  
    ! Check the consistancy between the communicators defined by the
    ! output arguments and by the input data file.
  call inquire(scat,root=root)
  ASSERT(root==fd%root)

    ! Get the 3rd dimenstion defined by the input arguments
  kl=localsize(scat) ! 0 or 1
  ASSERT(kl==size(v3d,3))
  km=totalsize(scat)

    ! Get the other dimensions defined by the input data file
  call inquire(fd%rf,nlon=nlon,nlat=nlat,nlev=nlev)

      ! Get the other dimensions defined by the input arguments
    im=size(v3d,1)
    jm=size(v3d,2)
    ASSERT(im==nlon)
    ASSERT(jm==nlat)
    ASSERT(km==nlev)

  itim_=1
  if(present(itim)) itim_=itim

  kl=0
  if(fd%isroot) kl=km
  allocate(vinp(im,jm,kl))

    ! Read from the input file
  if(fd%isroot) then
    call GrADS_read(fd%rf,vname,itim_,vinp(:,:,:),stat=ier)
      if(ier/=0) call die(myname_,'GrADS_read("'//trim(vname)//'")',ier)
  endif

    ! Send levels to their final places.
  call scatterv(scat,vinp,v3d,fd%root,comm)
  deallocate(vinp)
end subroutine read3dlev_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: write2dsub_ - write out a 2-d grid in Subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

subroutine write2dsub_(fd,vname,grid,v2d,comm,swapij)
  !-- write out a 2-d subdomain distribution
  use m_InterleavedObject,only : InterleavedObject
  use m_InterleavedObject,only : InterleavedObject_init
  use m_InterleavedObject,only : localsize
  use m_InterleavedObject,only : clean
  use m_SubdomainDistributor,only : SubdomainDistributor
  use m_SubdomainDistributor,only : get
  use m_SubdomainDistributorComm,only : undistribute
  use m_swapij,only : swap => swapij
  use m_GrADS,only : GrADS_write
  use m_die,only : die
      implicit none
  type(parGrADS),intent(inout) :: fd
  character(len=*),intent(in) :: vname
  type(SubdomainDistributor),intent(in) :: grid
  real,dimension(:,:),target,intent(in) :: v2d
  integer,intent(in) :: comm
  logical,optional,intent(in) :: swapij

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::write2dsub_'

  type(InterleavedObject) :: scat
  real,pointer,dimension(:,:,:) :: recv
  real,pointer,dimension(:,:) :: send
  logical :: swapij_
  integer :: im,il,jm,jl,mlev,ier

  swapij_=.false.
  if(present(swapij)) swapij_=swapij

  call get(grid,itotalsize=im,ilocalsize=il, &
		jtotalsize=jm,jlocalsize=jl  )
  if(swapij_) then
    allocate(send(il,jl))
    call swap(v2d,send)
  else
    send => v2d
  endif

  call InterleavedObject_init(scat,1,comm,fd%root)

  mlev=localsize(scat)
  allocate(recv(im,jm,mlev))
  call undistribute(scat,grid,send,recv,comm)

  if(swapij_) then
    deallocate(send)
  else
    nullify(send)
  endif

  if(fd%isroot) then
    call GrADS_write(fd%wf,vname,recv(:,:,1),stat=ier)
      if(ier/=0) call die(myname_, &
        'GrADS_write("'//trim(vname)//'")',ier)
  endif
  deallocate(recv)
  call clean(scat)
end subroutine write2dsub_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: write3dsub_ - write out a 3-d subdomain distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

subroutine write3dsub_(fd,vname,grid,v3d,comm,swapij)
  !-- write out a 3-d subdomain distribution
  use m_interleavedObject,only : interleavedObject
  use m_interleavedObject,only : interleavedObject_init
  use m_InterleavedObject,only : localsize
  use m_interleavedObject,only : clean
  use m_InterleavedComm,only : alloc_gatherv
  use m_SubdomainDistributor,only : SubdomainDistributor
  use m_SubdomainDistributor,only : get
  use m_SubdomainDistributorComm,only : undistribute
  use m_swapij,only : swap => swapij
  use m_GrADS,only : GrADS_write
  use m_die,only : die
  use m_mpout,only : mpout_log
  implicit none
  type(parGrADS),intent(inout) :: fd
  character(len=*),intent(in) :: vname
  type(SubdomainDistributor),intent(in) :: grid
  real,dimension(:,:,:),target,intent(in) :: v3d
  integer,intent(in) :: comm
  logical,optional,intent(in) :: swapij

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::write3dsub_'
  type(InterleavedObject) :: scat
  real,pointer,dimension(:,:,:) :: recv
  real,pointer,dimension(:,:,:) :: send
  real,allocatable,dimension(:,:,:) :: bufr
  logical :: swapij_
  integer :: im,il,jm,jl,km,mlev,ier

  swapij_=.false.
  if(present(swapij)) swapij_=swapij

  call get(grid,itotalsize=im,ilocalsize=il, &
		jtotalsize=jm,jlocalsize=jl  )
		
  km=size(v3d,3)
  if(swapij_) then
    allocate(send(il,jl,km))
    call swap(v3d,send)
  else
    send => v3d
  endif

  call InterleavedObject_init(scat,km,comm,fd%root)
  mlev=localsize(scat)
  allocate(bufr(im,jm,mlev))

  call undistribute(scat,grid,send,bufr,comm)
  if(swapij_) then
    deallocate(send)
  else
    nullify(send)
  endif

    ! Note that recv is only allocated on root at return.
  call alloc_gatherv(im,jm,scat,bufr,recv,fd%root,comm,myname_)
  deallocate(bufr)

  if(fd%isroot) then
    call GrADS_write(fd%wf,vname,recv(:,:,:),stat=ier)
      if(ier/=0) call die(myname_, &
        'GrADS_write("'//trim(vname)//'")',ier)
    deallocate(recv)
  endif
  call clean(scat)
end subroutine write3dsub_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: write2dlev_ - write a 2-d grid in interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

subroutine write2dlev_(fd,vname,scat,v2d,comm)
!-- write a 2d interleaved level
  use m_InterleavedObject,only : InterleavedObject
  use m_InterleavedComm,only : alloc_gatherv
  use m_GrADS,only : GrADS_write
  use m_die,only : die
  use m_die,only : assert_
  implicit none
  type(parGrADS),intent(inout) :: fd
  character(len=*),intent(in) :: vname
  type(InterleavedObject),intent(in) :: scat
  real,target,dimension(:,:,:),intent(in) :: v2d
  integer,intent(in) :: comm

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::write2dlev_'
  integer :: im,jm,mx,ij,ier
  real,pointer,dimension(:,:,:) :: send
  real,pointer,dimension(:,:) :: recv

  mx=size(v2d,3) ! mx == 0 or 1
    ASSERT(mx/2==0)

  im=size(v2d,1)
  jm=size(v2d,2)
  send => v2d
  call alloc_gatherv(im,jm,scat,send,recv,fd%root,comm,myname_)
  nullify(send)

  if(fd%isroot) then
    call GrADS_write(fd%wf,vname,recv(:,:),stat=ier)
      if(ier/=0) call die(myname_, &
        'GrADS_write("'//trim(vname)//'")',ier)
    deallocate(recv)
  endif
end subroutine write2dlev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: write3dlev_ - write a 3-d grid in interleaved levels
!
! !DESCRIPTION:
!
! !INTERFACE:

subroutine write3dlev_(fd,vname,scat,v3d,comm)
  use m_InterleavedObject,only : InterleavedObject
  use m_InterleavedComm,only : alloc_gatherv
  use m_GrADS,only : GrADS_write
  use m_die,only : die
  use m_die,only : assert_
  implicit none
  type(parGrADS),intent(inout) :: fd
  character(len=*),intent(in) :: vname
  type(InterleavedObject),intent(in) :: scat
  real,target,dimension(:,:,:),intent(in) :: v3d
  integer,intent(in) :: comm

! !REVISION HISTORY:
! 	02Jun06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::write3dlev_'
  integer :: im,jm,mx,ij,ier
  real,pointer,dimension(:,:,:) :: send
  real,pointer,dimension(:,:,:) :: recv

  mx=size(v3d,3) ! mx == 0 or 1
    ASSERT(mx/2==0)

  im=size(v3d,1)
  jm=size(v3d,2)
  send => v3d
  call alloc_gatherv(im,jm,scat,send,recv,fd%root,comm,myname_)
  nullify(send)

  if(fd%isroot) then
    call GrADS_write(fd%wf,vname,recv(:,:,:),stat=ier)
      if(ier/=0) call die(myname_, &
        'GrADS_write("'//trim(vname)//'")',ier)
    deallocate(recv)
  endif
end subroutine write3dlev_
end module m_parGrADS
