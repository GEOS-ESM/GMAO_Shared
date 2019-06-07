!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ggGradient - parallel gradient operator on Gaussian grid
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_ggGradient
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_interleavedObject   ,only : interleavedObject
      use m_ggGradientSP        ,only : ggGradientSP
      implicit none
      private	! except

      public :: ggGradient		! data structure
      public :: ggGradient_init,init
      public :: clean
      public :: ggGrad		! (x,y) = grad(f)
      public :: ggDivo		! d=div(x,y) & z=curl(x,y)

    type ggGradient
      private
      type(SubdomainDistributor) :: gGrid
      type(interleavedObject) :: intLev1
      type(interleavedObject) :: intLevs
      type(ggGradientSP) :: grad
    end type ggGradient

    interface ggGradient_init; module procedure	&
      gginit_,	&
      init_; end interface
    interface init ; module procedure	&
      gginit_,	&
      init_; end interface
    interface clean; module procedure clean_; end interface

    interface ggGrad; module procedure	&
      grad2_,	&
      grad3_; end interface

    interface ggDivo; module procedure	&
      divvor2_,	&
      divvor3_; end interface

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	01Dec05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- Appended pbr_ggDivo(), a non-module pass-by-reference
!		  interface of ggDivo(), to this file.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ggGradient'
  integer,parameter :: ROOT=0

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gginit_ - initialize this object on a Gaussian grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gginit_(ob,nGlon,nGlat,	&
    	lcLon,lnLon,lcLat,lnLat,nlevs,comm, radius)
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_interleavedObject   ,only : interleavedObject_init
      use m_ggGradientSP        ,only : ggGradientSP_init
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradient),intent(out) :: ob
      integer,intent(in) :: nGlon,nGlat		! nglon, nglat
      integer,intent(in) :: lcLon,lnLon		! local lon.-subdomain
      integer,intent(in) :: lcLat,lnLat		! local lat.-subdomain
      integer,intent(in) :: nlevs		! size of 3-d grid
      integer,intent(in) :: comm		! communicator
      real,optional,intent(in) :: radius	! default to Earth

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gginit_'
  integer :: ier
_ALLENTRY_

  	! Define a ggGradientgrid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
	nGlon,lcLon,lnLon,.true.,		&
	nGlat,lcLat,lnLat,.false., comm		)

	! Define two vertical level distribution controllers for
	! 2-d and 3-d grids.

  call interleavedObject_init(ob%intLev1,    1,comm,root=ROOT)
  call interleavedObject_init(ob%intLevs,nlevs,comm,root=ROOT)

  	! Define a single processor gradient operator
  call ggGradientSP_init(ob%grad,nGlon,nGlat,radius=radius)
_ALLEXIT_
end subroutine gginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize this object with given latitudes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,nGlon,rlats,	&
    	lcLon,lnLon,lcLat,lnLat,nlevs,comm, radius)
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_interleavedObject   ,only : interleavedObject_init
      use m_ggGradientSP        ,only : ggGradientSP_init
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradient),intent(out) :: ob
      integer,intent(in) :: nGlon		! nglon
      real,dimension(:),intent(in) :: rlats	! a given list of lat.
      integer,intent(in) :: lcLon,lnLon		! local lon.-subdomain
      integer,intent(in) :: lcLat,lnLat		! local lat.-subdomain
      integer,intent(in) :: nlevs		! size of 3-d grid
      integer,intent(in) :: comm		! communicator
      real,optional,intent(in) :: radius	! default to Earth

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nlat_
  integer :: ier
_ALLENTRY_
  nlat_=size(rlats)	! known global latitudes

  	! Define a ggGradientgrid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
	nGlon,lcLon,lnLon,.true.,		&
	nlat_,lcLat,lnLat,.false., comm		)

	! Define two vertical level distribution controllers for
	! 2-d and 3-d grids.

  call interleavedObject_init(ob%intLev1,    1,comm,root=ROOT)
  call interleavedObject_init(ob%intLevs,nlevs,comm,root=ROOT)

  	! Define a single processor gradient operator
  call ggGradientSP_init(ob%grad,nGlon,rlats,radius=radius)
_ALLEXIT_
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean this object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_SubdomainDistributor,only : clean
      use m_interleavedObject   ,only : clean
      use m_ggGradientSP        ,only : clean
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradient),intent(inout) :: ob

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_

  call clean(ob%gGrid)
  call clean(ob%intLev1)
  call clean(ob%intLevs)
  call clean(ob%grad)
_EXIT_
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divvor2_ - get vor and div for a rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divvor2_(ob,u,v,div,vor,comm)
      use m_SubdomainDistributorComm,only : undistribute
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : totalSize
      use m_swapij,only : swapij
      use m_ggGradientSP,only : ggDivo
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradient),intent(in) :: ob		! ggGradient operator
      real,dimension(:,:),intent(in ) :: u  ,v	! a vector field
      real,dimension(:,:),intent(out) :: vor,div ! curl(u,v) & div(u,v)
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divvor2_'
  integer :: mlon,mlat		! sizes of arguments are explicitly
				! named for clearity.
  integer :: nGLon,lnLon	! total and "my" longitude counts
  integer :: nGLat,lnLat	! total and "my" latitude counts
  integer :: nLev1,lnLev	! total and "my" level counts

  real,allocatable,dimension(:,:  ) :: u_GS,v_GS,div_GS,vor_GS
  real,allocatable,dimension(:,:,:) :: u_GL,v_GL,div_GL,vor_GL
  integer :: ier
!________________________________________
		! sizes of subdomain storage
_ALLENTRY_

	mlon=size(u,2)
	mlat=size(u,1)

  call get(ob%gGrid,ilocalSize=lnLon,jlocalSize=lnLat)
  nLev1=totalSize(ob%intLev1)

		! consistance check

  	ASSERT(mlon==lnLon)
  	ASSERT(mlat==lnLat)
  	ASSERT(   1==nLev1)
  
  		! sizes of interleaved storage

  call get(ob%gGrid,itotalSize=nGLon,jtotalSize=nGLat)
  lnLev=localSize(ob%intLev1)

  	ASSERT(lnLev==0.or.lnLev==1)
!________________________________________
	allocate( u_GS(mlon,mlat),	&
		  v_GS(mlon,mlat))

  	! swap the storage from (lat,lon) to (lon,lat)
  call swapij(u,u_GS)
  call swapij(v,v_GS)

  	allocate( u_GL(nGlon,nGlat,lnLev),	&
		  v_GL(nGlon,nGlat,lnLev))

  call undistribute(ob%intLev1,ob%gGrid, u_GS,u_GL, comm)
  call undistribute(ob%intLev1,ob%gGrid, v_GS,v_GL, comm)

  	deallocate( u_GS, v_GS)
	allocate( div_GL(nGlon,nGlat,lnLev),	&
		  vor_GL(nGlon,nGlat,lnLev)	)

  call ggDivo(ob%grad,u_GL,v_GL,div_GL,vor_GL)

  	deallocate( u_GL, v_GL)
	allocate( div_GS(mlon,mlat),	&
		  vor_GS(mlon,mlat)	)

  call distribute(ob%intLev1,ob%gGrid, div_GL,div_GS, comm)
  call distribute(ob%intLev1,ob%gGrid, vor_GL,vor_GS, comm)

  	deallocate( div_GL, vor_GL)

  	! swap the storage back to (lat,lon)

  call swapij(div_GS,div)
  call swapij(vor_GS,vor)

  	deallocate( div_GS, vor_GS)
_ALLEXIT_
end subroutine divvor2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divvor3_ - get vor and div for a rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divvor3_(ob,u,v,div,vor,comm)
      use m_SubdomainDistributorComm,only : undistribute
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : totalSize
      use m_ggGradientSP,only : ggDivo
      use m_swapij,only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradient),intent(in) :: ob		! ggGradient operator
      real,dimension(:,:,:),intent(in ) :: u  ,v	! a vector field
      real,dimension(:,:,:),intent(out) :: vor,div ! curl(u,v) & div(u,v)
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divvor3_'
  integer :: mlon,mlat,mlev	! sizes of arguments are explicitly
				! named for clearity.
  integer :: nGLon,lnLon	! total and "my" longitude counts
  integer :: nGLat,lnLat	! total and "my" latitude counts
  integer :: nLevs,lnLev	! total and "my" level counts
  integer :: ier
!________________________________________

  	! variables stored in subdomains
  real,allocatable,dimension(:,:,:) :: u_GS,v_GS,div_GS,vor_GS

  	! variables stored in interleaved levels
  real,allocatable,dimension(:,:,:) :: u_GL,v_GL,div_GL,vor_GL
!________________________________________
		! sizes of subdomain storage
_ALLENTRY_

	mlon=size(u,2)
	mlat=size(u,1)
	mlev=size(u,3)

  call get(ob%gGrid,ilocalSize=lnLon,jlocalSize=lnLat)
  nLevs=totalSize(ob%intLevs)

		! consistance check

  	ASSERT(mlon==lnLon)
  	ASSERT(mlat==lnLat)
  	ASSERT(mlev==nLevs)
  
  		! sizes of interleaved storage

  call get(ob%gGrid,itotalSize=nGLon,jtotalSize=nGLat)
  lnLev=localSize(ob%intLevs)

	allocate( u_GS(mlon,mlat,mlev),	&
		  v_GS(mlon,mlat,mlev))

  	! swap the storage from (lat,lon,lev) to (lon,lat,lev)
  call swapij(u,u_GS)
  call swapij(v,v_GS)

  	allocate( u_GL(nGLon,nGLat,lnLev),	&
		  v_GL(nGLon,nGLat,lnLev))

	! move data from subdomains to interleaved levels
  call undistribute(ob%intLevs,ob%gGrid, u_GS,u_GL, comm)
  call undistribute(ob%intLevs,ob%gGrid, v_GS,v_GL, comm)

  	deallocate( u_GS, v_GS)
	allocate( div_GL(nGLon,nGLat,lnLev),	&
		  vor_GL(nGLon,nGLat,lnLev)	)

  	! compute divergence and vortisity

  call ggDivo(ob%grad,u_GL,v_GL,div_GL,vor_GL)

  	deallocate( u_GL, v_GL)
	allocate( div_GS(mlon,mlat,mlev),	&
		  vor_GS(mlon,mlat,mlev)	)

  	! move data from interleaved levels back to subdomains
  call distribute(ob%intLevs,ob%gGrid, div_GL,div_GS, comm)
  call distribute(ob%intLevs,ob%gGrid, vor_GL,vor_GS, comm)

  	deallocate( div_GL, vor_GL)

  	! swap the storage back to (lat,lon,lev)

  call swapij(div_GS,div)
  call swapij(vor_GS,vor)

  	deallocate( div_GS, vor_GS)
_ALLEXIT_
end subroutine divvor3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: grad2_ - get vor and div for a rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine grad2_(ob,f,x,y,comm)
      use m_SubdomainDistributorComm,only : undistribute
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : totalSize
      use m_swapij,only : swapij
      use m_ggGradientSP,only : ggGrad
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradient),intent(in) :: ob		! ggGradient operator
      real,dimension(:,:),intent(in ) :: f	! a scalar field
      real,dimension(:,:),intent(out) :: x,y	! (x,y)=grad(f)
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::grad2_'
  integer :: mlon,mlat		! sizes of arguments are explicitly
				! named for clearity.
  integer :: nGLon,lnLon	! total and "my" longitude counts
  integer :: nGLat,lnLat	! total and "my" latitude counts
  integer :: nLev1,lnLev	! total and "my" level counts

  real,allocatable,dimension(:,:  ) :: f_GS,x_GS,y_GS
  real,allocatable,dimension(:,:,:) :: f_GL,x_GL,y_GL
  integer :: ier
!________________________________________
		! sizes of subdomain storage
_ALLENTRY_

	mlon=size(f,2)
	mlat=size(f,1)

  call get(ob%gGrid,ilocalSize=lnLon,jlocalSize=lnLat)
  nLev1=totalSize(ob%intLev1)

		! consistance check

  	ASSERT(mlon==lnLon)
  	ASSERT(mlat==lnLat)
  	ASSERT(   1==nLev1)
  
  		! sizes of interleaved storage

  call get(ob%gGrid,itotalSize=nGLon,jtotalSize=nGLat)
  lnLev=localSize(ob%intLev1)

  	ASSERT(lnLev==0 .or. lnLev==1)
!________________________________________
	allocate(f_GS(mlon,mlat))

  	! swap the storage from (lat,lon) to (lon,lat)
  call swapij(f,f_GS)

  	allocate(f_GL(nGlon,nGlat,lnLev))

  call undistribute(ob%intLev1,ob%gGrid, f_GS,f_GL, comm)

  	deallocate(f_GS)
	allocate( x_GL(nGlon,nGlat,lnLev),	&
		  y_GL(nGlon,nGlat,lnLev)	)

  call ggGrad(ob%grad,f_GL,x_GL,y_GL)

  	deallocate(f_GL)
	allocate( x_GS(mlon,mlat),	&
		  y_GS(mlon,mlat)	)

  call distribute(ob%intLev1,ob%gGrid, x_GL,x_GS, comm)
  call distribute(ob%intLev1,ob%gGrid, y_GL,y_GS, comm)

  	deallocate( x_GL, y_GL)

  	! swap the storage back to (lat,lon)

  call swapij(x_GS,x)
  call swapij(y_GS,y)

  	deallocate( x_GS, y_GS)
_ALLEXIT_
end subroutine grad2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: grad3_ - get vor and div for a rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine grad3_(ob,f,x,y,comm)
      use m_SubdomainDistributorComm,only : undistribute
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : totalSize
      use m_swapij,only : swapij
      use m_ggGradientSP,only : ggGrad
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradient),intent(in) :: ob		! ggGradient operator
      real,dimension(:,:,:),intent(in ) :: f	! a scalar field
      real,dimension(:,:,:),intent(out) :: x,y	! (x,y)=grad(f)
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::grad3_'
  integer :: mlon,mlat,mlev	! sizes of arguments are explicitly
				! named for clearity.
  integer :: nGLon,lnLon	! total and "my" longitude counts
  integer :: nGLat,lnLat	! total and "my" latitude counts
  integer :: nLevs,lnLev	! total and "my" level counts

  real,allocatable,dimension(:,:,:) :: f_GS,x_GS,y_GS
  real,allocatable,dimension(:,:,:) :: f_GL,x_GL,y_GL
  integer :: ier
!________________________________________
		! sizes of subdomain storage
_ALLENTRY_

	mlon=size(f,2)
	mlat=size(f,1)
	mlev=size(f,3)

  call get(ob%gGrid,ilocalSize=lnLon,jlocalSize=lnLat)
  nLevs=totalSize(ob%intLevs)

		! consistance check

  	ASSERT(mlon==lnLon)
  	ASSERT(mlat==lnLat)
  	ASSERT(mlev==nLevs)
  
  		! sizes of interleaved storage

  call get(ob%gGrid,itotalSize=nGLon,jtotalSize=nGLat)
  lnLev=localSize(ob%intLevs)
!________________________________________
	allocate(f_GS(mlon,mlat,mlev))

  	! swap the storage from (lat,lon,lev) to (lon,lat,lev)
  call swapij(f,f_GS)

  	allocate(f_GL(nGlon,nGlat,lnLev))

  call undistribute(ob%intLevs,ob%gGrid, f_GS,f_GL, comm)

  	deallocate(f_GS)
	allocate( x_GL(nGlon,nGlat,lnLev),	&
		  y_GL(nGlon,nGlat,lnLev)	)

  call ggGrad(ob%grad,f_GL,x_GL,y_GL)

  	deallocate(f_GL)
	allocate( x_GS(mlon,mlat,mlev),	&
		  y_GS(mlon,mlat,mlev)	)

  call distribute(ob%intLevs,ob%gGrid, x_GL,x_GS, comm)
  call distribute(ob%intLevs,ob%gGrid, y_GL,y_GS, comm)

  	deallocate( x_GL, y_GL)

  	! swap the storage back to (lat,lon,lev)

  call swapij(x_GS,x)
  call swapij(y_GS,y)

  	deallocate( x_GS, y_GS)
_ALLEXIT_
end subroutine grad3_

end module m_ggGradient
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: pbr_ggDivo - a pass-by-reference interface of ::ggDivo()
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine pbr_ggDivo(ob,lat2,lon2,nsig,u,v,div,vor,comm)
      use m_ggGradient,only : ggGradient,ggDivo
      implicit none
      type(ggGradient),intent(in) :: ob		! ggGradient operator
      integer,intent(in) :: lat2,lon2,nsig
      real,dimension(lat2,lon2,nsig),intent(in ) :: u  ,v
      real,dimension(lat2,lon2,nsig),intent(out) :: vor,div
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	01Dec05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  call ggDivo(ob,u,v,div,vor,comm)
end subroutine pbr_ggDivo
