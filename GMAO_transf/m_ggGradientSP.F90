!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ggGradientSP - low-order 2-d gradient on an any-lat grid.
!
! !DESCRIPTION:
!
!     ggGradientSP computes low-order numerical gradient on a 2-d
!   global grid without reference or correction to the vertical grid.
!   The ordering of the 2-d grid points may be either longitude-major
!   (lon,lat) or latitude-major (lat,lon), as long as it is specified.
!
!     For notation purposes, operators supported in this module are
!   defined as such, for scalar field f:
!
!	grad(f) := ( fx, fy)	: 2d-gradient of f
!	strm(f) := (-fy, fx)	: -curl(f*ek), or k.cross.grad(f)
!
!   where fx := pf/px and fy := pf/py are _partial_ derivatives; and
!   for vector field (u,v):
!
!	dive(u,v) := pu/px + pu/py	: divergence of (u,v)
!	vort(u,v) :=-pu/py + pu/px	: vortisity is 2-d curl(u,v)
!
!     There are also the following second derivative operations:
!
!	lapla(f) := dive(grad(f))
!	psolv(z) := solve lapla(f)=z for f	! const.dlat for now.
!
!   For example, from (d,z) to (u,v), it includes the following
!   two steps:
!
!	1. x=psolv(d); s=psolv(z)	: from (d,z) to (x,s)
!	2. (u,v)=grad(x)+strm(s)	! derive (u,v) fro (x,s)
!
!   For vector fields:
!
!         /--------v
!       xs -> uv -> dz
!         ^--------/
!
!.	uv -> dz			: d=div(u,v), z=vor(u,v)
!.	      dz -> xs			: x=solve(d), s=solve(z)
!.	            xs -> uv		: (u,v)=grad(x)+strm(s)
!?	uv -------> xs			: uv->dz; dz->xs  (uv2xs)
!.	      dz -------> uv		: dz->xs; xs->uv  (dz2uv)
!.		    xs -------> dz	: d=lapla(x), z=lapla(s)
!
! !INTERFACE:
!#include "regime.H"

    module m_ggGradientSP
      use m_elliptic, only: PWSSSP
      implicit none
      private	! except

      public :: ggGradientSP

      	! _TODO_: Need a new name to redefine the init_() interface.
	!	  However, it seems difficult.  I might just go ahead
	!	  using the same name with a new definition.

      public :: ggGradientSP_init,init
      public :: ggGradientSP_gginit
      public :: ggGradientSP_llinit

      public :: clean
      public :: ggGrad	! grad(f):=(d/dx,d/dy)(f); f is a scalar field
      public :: ggStrm	! strm(f):=(-d/dy,d/dx)(f); f is a scalar field
      public :: ggdive  ! d = dive(fx,fy)
      public :: ggvort  ! z = vort(fx,fy) = dive(-fy,fx)
      public :: ggLapla ! solve d=lapla(f) for d
      public :: ggPsolv ! Poisson solver of lapla(f)=d for f

      public :: ggDivo	! div() and vor() of a vector field (u,v)
      			! solve {d,z} = {div(u,v),vor(u,v)} for {d,z}

      type ggGradientSP
	private
        real,pointer,dimension(:) :: sinlam
        real,pointer,dimension(:) :: coslam

	real :: radius
        real,pointer,dimension(:) :: r2dx
        real,pointer,dimension(:) :: rdyp
        real,pointer,dimension(:) :: rdym
        real,pointer,dimension(:) :: cosphi
	logical :: evenlat
      end type ggGradientSP

    interface ggGradientSP_init; module procedure	&
    	init_,	&
	gginit_; end interface
    interface init; module procedure init_, gginit_; end interface
    interface ggGradientSP_gginit; module procedure	&
	gginit_; end interface
    interface ggGradientSP_llinit; module procedure	&
	llinit_; end interface
    interface clean; module procedure clean_; end interface

    interface ggGrad; module procedure		&
      gradr2_,	&
      gradr3_; end interface
    interface ggStrm; module procedure		&
      strmr2_,	&
      strmr3_; end interface
    interface ggdive; module procedure		&
      diver2_,	&
      diver3_; end interface
    interface ggvort; module procedure		&
      vortr2_,	&
      vortr3_; end interface
    interface ggLapla; module procedure &
      lapla2_, &
      lapla3_; end interface
    interface ggPsolv; module procedure &
      ellip2_, &
      ellip3_; end interface

    interface ggDivo; module procedure		&
      divvor2_,	&
      divvor3_; end interface

!    interface xs2uv; module procedure &
!      xs2uv2_, &
!      xs2uv3_; end interface
!    interface uv2xs; module procedure &
!      uv2xs2_, &
!      uv2xs3_; end interface
!    interface uv2dz; module procedure &
!      uv2dz2_, &
!      uv2dz3_; end interface
    interface dz2uv; module procedure &
      dz2uv2_, &
      dz2uv3_; end interface

!    interface xs2dv; module procedure &
!      xs2dv2_, &
!      xs2dv3_; end interface
!    interface dv2xs; module procedure &
!      dv2xs2_, &
!      dv2xs3_; end interface

! !REVISION HISTORY:
! 	03Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ggGradientSP'

  logical, save :: verbose_ = .false.

! Usecases:
!
!   use m_ggGradientSP,only : ggGradientSP,ggGradientSP_init,clean
!   use m_ggGradientSP,only : ggDivo,ggGrad
!
!   real,dimension(:,:,:) :: f3,gx3,gy3,div3,vor3
!   real,dimension(:,:) :: f2,gx2,gy2,div2,vor2
!
!
!	type(ggGradientSP) :: gr
!
!	call ggGradientSP_init(gr,im,glats)
!
!	call ggGrad(gr,f2,gx2,gy2)		! (gx2,gy2)=grad(f2)
!	call ggGrad(gr,f3,gx3,gy3)		! (gx3,gx3)=grad(f3)
!	call ggDivo(gr,gx2,gy2,div2,vor2)	! div2=div(gx2,gy2)
!	call ggDivo(gr,gx3,gy3,div3,vor3)	! vor2=curl(gx2,gy2)
!
!	call clean(gr)

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gginit_ - initialize assuming a global gaussian grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gginit_(ob,nglon,nglat,radius,verbose)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(out) :: ob	! this operator
      integer,intent(in) :: nglon	! No. of longitudes in 360.
      integer,intent(in) :: nglat	! No. of Gaussian latitudes
      real,optional,intent(in ) :: radius	! radius of the sphere
      logical,optional,intent(in) :: verbose    ! sets verbose level

! !REVISION HISTORY:
! 	21Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	17Aug07 - Todling, implemented verbose
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gginit_'
  real,dimension(nglat) :: glats
_ENTRY_

  	ASSERT(nglon>1)
  	ASSERT(nglat>2)

  	! Set Gaussian latitude grid points in degrees without poles.

  call gauss_lat_nmc(glats(2:nglat-1),nglat-2)

  	! Set Gaussian latitude grid points at the poles.

  glats(1) = -90.		! at the South Pole
  glats(nglat) = +90.		! at the North Pole

  call init_(ob,nglon,glats,radius=radius,verbose=verbose)
_EXIT_
end subroutine gginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: llinit_ - initialize assuming a global lat-lon grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine llinit_(ob,nglon,nglat,radius,verbose)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(out) :: ob	! this operator
      integer,intent(in) :: nglon	! No. of longitudes in 360.
      integer,intent(in) :: nglat	! No. of Gaussian latitudes
      real,optional,intent(in ) :: radius	! radius of the sphere
      logical,optional,intent(in) :: verbose    ! sets verbose level

! !REVISION HISTORY:
! 	21Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	17Aug07 - Todling, implemented verbose
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::llinit_'
  real,dimension(nglat) :: glats
  real :: dlat
  integer :: j
_ENTRY_

  	ASSERT(nglon>1)
  	ASSERT(nglat>2)

  	! Set Gaussian latitude grid points in degrees without poles.

  dlat=180./max(1,nglat-1)
  do j=1,nglat
    glats(j)=-90.+(j-1)*dlat
  end do

  call init_(ob,nglon,glats,radius=radius,verbose=verbose)
_EXIT_
end subroutine llinit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,nlon,glats,radius,verbose)
      use m_geosapi,only : getcon
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(out) :: ob	! this operator
      integer         ,intent(in ) :: nlon	! longitudes in 360.
      real,dimension(:),intent(in) :: glats	! latitudes, SP. to NP.
      real,optional  ,intent(in ) :: radius	! radius of the sphere
      logical,optional,intent(in) :: verbose    ! sets verbose level

      		! Assumptions about the grid include:
		!	. an equal-interval longitude grid
		!	. an any-interval monotonic latitude grid
		!	. two pole points are included in glats(:).

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	17Aug07 - Todling, implemented verbose
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  real :: dlon,x,dlatm,dlatp
  real :: r2dy,cosphi
  real :: RAD
  integer :: i,j,nlat
  real :: radius_
_ENTRY_

  nlat=size(glats)

		! check for minimum arithmatic requirement

  	ASSERT(nlon>=1)
  	ASSERT(mod(nlon,2)==0)
	ASSERT(nlat>=2)

		! It is also expected, but not checked, that:
		!   glats(   1) == S.P. (== -90.)
		!   glats(nlat) == N.P. (== +90.)
		!   are_monotopic(glats)

  RAD = 4.*atan(1.)/180.	! atan(1.)==pi/4

  	! override the default, "EARTH RADIUS"

  radius_ = getcon('EARTH RADIUS')
  if(present(radius)) radius_=radius
  if(present(verbose)) verbose_=verbose
!________________________________________

	allocate( ob%sinlam(nlon),ob%coslam(nlon),	&
		  ob%  r2dx(nlat),ob%cosphi(nlat),	&
		  ob%  rdyp(nlat),ob%  rdym(nlat)	)
!________________________________________
  dlon=360./nlon
  do i=1,nlon
    x=(i-1)*dlon*RAD
    ob%sinlam(i)=sin(x)
    ob%coslam(i)=cos(x)
  end do
!________________________________________

  dlatp = glats(2)-glats(1)

  r2dy = 1./(2.*radius_*dlatp*RAD)

  ob%radius  = radius_
  ob%r2dx(1) = r2dy
  ob%rdyp(1) = r2dy
  ob%rdym(1) = r2dy
  ob%cosphi(1) = cos(glats(1)*RAD)
	!________________________________________

  do j=2,nlat-1
    		! r2dx = 1./( 2 * ER * dlon * cos(phi) )

    x = glats(j)*RAD
    ob%cosphi(j)=cos(x)
    ob%r2dx(j)=1./(2.*radius_*ob%cosphi(j)*dlon*RAD)

		! r2dy = 1./( ER * (phi(j+1)-phi(j-1)) )
    		! rdyp = (phi(j+1)-phi(j))/(phi(j)-phi(j-1)) * r2dy
    		! rdym = (phi(j)-phi(j-1))/(phi(j+1)-phi(j)) * r2dy

    dlatm = glats(j)-glats(j-1)
    dlatp = glats(j+1)-glats(j)

    		! dlatp is used as phi(j+1)-phi(j) in this iteration,
		! and as phi(j)-phi(j-1) in the next iteration.

    r2dy = 1./( radius_*(glats(j+1)-glats(j-1))*RAD )

    ob%rdyp(j) = dlatp/dlatm * r2dy
    ob%rdym(j) = dlatm/dlatp * r2dy
  end do
	!________________________________________

  dlatm= glats(nlat)-glats(nlat-1)
  r2dy = 1./( 2.*radius_*dlatm*RAD )
  ob%r2dx(nlat) = r2dy
  ob%rdyp(nlat) = r2dy
  ob%rdym(nlat) = r2dy
  x = glats(nlat)*RAD
  ob%cosphi(nlat) = cos(x)

  	! Determine if the given latitide grid has a constant interval.

  dlatp=maxval( glats(2:nlat)-glats(1:nlat-1) )		! maximum value
  dlatm=minval( glats(2:nlat)-glats(1:nlat-1) )		! minimum value
  x = (glats(nlat)-glats(1))/(nlat-1)			! mean value

	! For <a=n*b>, the <truncation error of variable a> is no more
	! than <n> times the <truncation error of variable b>.  Taking
	! <n> at its maximum, one gets an estimated <maximum error of
	! variable a> from an estimated <maximum error of variable b>.
	! Thus, the condition to consider a given list of values to be
	! practically the same, is defined as
	!
	!   maxval-minmin <= maxerr
	!
	! where, for the given latitude grid value problem,
	!
	!   maxerr/mean := nlat*epsilon()

  ob%evenlat = dlatp-dlatm <= abs(x)*nlat*epsilon(x)

	! A real case testing suggests that
	!
	!   maxval-minval ~= 1/2 maxerr

_EXIT_
end subroutine init_
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
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(inout) :: ob

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	17Aug07 - Todling, implemented verbose
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_

	deallocate( ob%sinlam,ob%coslam,	&
  		    ob%  r2dx,ob%cosphi,	&
		    ob%  rdyp,ob%  rdym	)

        verbose_=.false. ! reset verbose
_EXIT_
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divrlp_ - 2-D divergence of a vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divrlp_(u,v,div,r2dx,cosphi,r2dyp,r2dym)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none

      real,dimension(:,:),intent(in) :: u,v
      			! (u,v) is a vector wind field

      real,dimension(:,:),intent(out) :: div
      			! div(u,v) := du/dx+d(v*cosphi)/(cosphi*dy)

      real,dimension(:),intent(in) :: r2dx	! := 1/(a*cosphi*2*dlam)
      real,dimension(:),intent(in) :: cosphi	! := cos(phi)
      real,dimension(:),intent(in) :: r2dyp
      			! := dphi(j  )/[a*dphi(j-1)*(dphi(j-1)+dphi(j)]

      real,dimension(:),intent(in) :: r2dym
      			! := dphi(j-1)/[a*dphi(j  )*(dphi(j-1)+dphi(j)]
			! where dphi(j) := phi(j+1)-phi(j-1)

! !REVISION HISTORY:
! 	11Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- fixed div(u,v) calculation bugs
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divrlp_'
  integer :: nx,ny,i,j
_ENTRY_

  nx=size(div,1)
  ny=size(div,2)

	ASSERT(nx>0)
	ASSERT(ny>1)
	ASSERT(nx==size(u,1))
	ASSERT(ny==size(u,2))
	ASSERT(nx==size(v,1))
	ASSERT(ny==size(v,2))


  	! at the South Pole, j=1

  div(:,1)=+sum(v(:,2))*r2dx(1)*4./nx

  	! at all grid points excluding the poles.

  	! Cache v(2)*cosphi(2)-v(1)*cosphi(1)

  div(:,2)=v(:,2)*cosphi(2)-v(:,1)*cosphi(1)
  do j=2,ny-1
    	! cache v(j+1)*cosphi(j+1)-v(j)*cosphi(j)
    div(:,j+1)=v(:,j+1)*cosphi(j+1)-v(:,j)*cosphi(j)

    	! at the same time, du(j) and dv(j) are available from vor(j)
	! and div(j), respectively.  Therefore,
	! du/dy=du(j)*r2dyp+du(j+1)*r2dy

    	! r2dyp(j) := dy(j)/[a*dy(j-1)*(dy(j)+dy(j-1))]
    	! r2dym(j) := dy(j-1)/[a*dy(j)*(dy(j)+dy(j-1))]

    	! d(v*cosphi)/(cosphi*dy)

    div(:,j)=(div(:,j)*r2dyp(j)+div(:,j+1)*r2dym(j))/cosphi(j)

    	! div = du/dx + d(v*cosphi)/(cosphi*dy)

    div( 1     ,j)=(u(2   ,j)-u(nx    ,j))*r2dx(j) + div(1     ,j)
    div( 2:nx-1,j)=(u(3:nx,j)-u(1:nx-2,j))*r2dx(j) + div(2:nx-1,j)
    div(nx     ,j)=(u(1   ,j)-u(nx-1  ,j))*r2dx(j) + div(nx    ,j)
  end do

  	! at the North Pole, j=ny

  div(:,ny)=-sum(v(:,ny-1))*r2dx(ny)*4./nx
_EXIT_
end subroutine divrlp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divrpl_ - 2-D divergence of a latitude-major vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divrpl_(u,v,div,r2dx,cosphi,r2dyp,r2dym)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none

      real,dimension(:,:),intent(in) :: u,v
      			! (u,v) is a vector wind field

      real,dimension(:,:),intent(out) :: div
      			! div(u,v) := du/dx+d(v*cosphi)/(cosphi*dy)

      real,dimension(:),intent(in) :: r2dx	! := 1/(a*cosphi*2*dlam)
      real,dimension(:),intent(in) :: cosphi	! := cos(phi)
      real,dimension(:),intent(in) :: r2dyp
      			! := dphi(j  )/[a*dphi(j-1)*(dphi(j-1)+dphi(j)]

      real,dimension(:),intent(in) :: r2dym
      			! := dphi(j-1)/[a*dphi(j  )*(dphi(j-1)+dphi(j)]
			! where dphi(j) := phi(j+1)-phi(j-1)

! !REVISION HISTORY:
! 	11Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- fixed div calculation bugs
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divrpl_'
  integer :: nx,ny,i,j
_ENTRY_

  nx=size(div,2)
  ny=size(div,1)
!!  call shownorm2_(myname_,'u',u)
!!  call shownorm2_(myname_,'v',v)

	ASSERT(nx>0)
	ASSERT(ny>1)
	ASSERT(nx==size(u,2))
	ASSERT(ny==size(u,1))
	ASSERT(nx==size(v,2))
	ASSERT(ny==size(v,1))
!________________________________________
  	! at the South Pole, j=1

  div(1,:)=+sum(v(2,:))*r2dx(1)*4./nx

  	! at all grid points excluding the poles.

  	! Cache v(2)*cosphi(2)-v(1)*cosphi(1)

  div(2,:)=v(2,:)*cosphi(2)-v(1,:)*cosphi(1)
!!  call shownorm2_(myname_,'div.1',div(1:1,:))
  do j=2,ny-1
    	! cache v(j+1)*cosphi(j+1)-v(j)*cosphi(j)
    div(j+1,:)=v(j+1,:)*cosphi(j+1)-v(j,:)*cosphi(j)

    	! at the same time, du(j) and dv(j) are available from vor(j)
	! and div(j), respectively.  Therefore,
	! du/dy=du(j)*r2dyp+du(j+1)*r2dy

    	! r2dyp(j) := dy(j)/[a*dy(j-1)*(dy(j)+dy(j-1))]
    	! r2dym(j) := dy(j-1)/[a*dy(j)*(dy(j)+dy(j-1))]

    	! d(v*cosphi)/(cosphi*dy)

    div(j,:)=(div(j,:)*r2dyp(j)+div(j+1,:)*r2dym(j))/cosphi(j)

    	! div = du/dx + d(v*cosphi)/(cosphi*dy)

    div(j, 1     )=(u(j,2   )-u(j,nx    ))*r2dx(j) + div(j,1     )
    div(j, 2:nx-1)=(u(j,3:nx)-u(j,1:nx-2))*r2dx(j) + div(j,2:nx-1)
    div(j,nx     )=(u(j,1   )-u(j,nx-1  ))*r2dx(j) + div(j,nx    )
!!  call shownorm2_(myname_,'div.j',div(1:j,:))
  end do

  	! at the North Pole, j=ny

  div(ny,:)=-sum(v(ny-1,:))*r2dx(ny)*4./nx
!!  call shownorm2_(myname_,'div',div)
_EXIT_
end subroutine divrpl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divvor2_ - divergence and vortisity of rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divvor2_(ob,u,v,div,vor,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: u,v
      real,dimension(:,:),optional,intent(out) :: div,vor
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divvor2_'
  logical :: latlon_
_ENTRY_

  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  if(latlon_) then
  	! curl(u,v) == div(v,-u)
    if(present(vor)) &
      call divrpl_(v,-u,vor,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    if(present(div)) &
      call divrpl_(u, v,div,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
  else
  	! curl(u,v) == div(v,-u)
    if(present(vor)) &
      call divrlp_(v,-u,vor,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    if(present(div)) &
      call divrlp_(u, v,div,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
  endif
_EXIT_
end subroutine divvor2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: divvor3_ - divergence and vortisity of rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine divvor3_(ob,u,v,div,vor,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: u,v
      real,dimension(:,:,:),optional,intent(out) :: div,vor
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::divvor3_'
  integer :: l,n
  logical :: latlon_
_ENTRY_

  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  n=size(u,3)
  	ASSERT(n==size(v,3))
  if(present(div)) then
  	ASSERT(n==size(div,3))
  endif
  if(present(vor)) then
  	ASSERT(n==size(vor,3))
  endif
!!  call shownorm3_(myname_,'u',u)
!!  call shownorm3_(myname_,'v',v)

  if(latlon_) then
    do l=1,n
  	! curl(u,v) == div(v,-u)
      if(present(vor)) &
        call divrpl_(v(:,:,l),-u(:,:,l),vor(:,:,l),		&
      			ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
      if(present(div)) &
        call divrpl_(u(:,:,l), v(:,:,l),div(:,:,l),		&
      			ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    end do
  else
    do l=1,n
  	! curl(u,v) == div(v,-u)
      if(present(vor)) &
        call divrlp_(v(:,:,l),-u(:,:,l),vor(:,:,l),		&
      			ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
      if(present(div)) &
        call divrlp_(u(:,:,l), v(:,:,l),div(:,:,l),		&
      			ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    end do
  endif
!!  if(present(div)) call shownorm3_(myname_,'div',div)
!!  if(present(vor)) call shownorm3_(myname_,'vor',vor)
_EXIT_
end subroutine divvor3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gradlp_ - gradiant at poles
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gradlp_(slam,clam,r2dx,rdyp,rdym,f,xk,yk)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:),intent(in) :: slam,clam
      real,dimension(:),intent(in) :: r2dx,rdyp,rdym
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: xk,yk

! !REVISION HISTORY:
! 	07Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gradlp_'
  real :: ys,yc, xr,yr, d, rd, sf,cf
  integer :: i,m,j
  integer :: nlon,nlat
_ENTRY_

  nlon=size(f,1)
  nlat=size(f,2)

  	ASSERT(nlon==size(slam))
  	ASSERT(nlon==size(clam))
  	ASSERT(nlat==size(r2dx))
  	ASSERT(nlat==size(rdyp))
  	ASSERT(nlat==size(rdym))

  	ASSERT(nlon==size(xk,1))
  	ASSERT(nlat==size(xk,2))
  	ASSERT(nlon==size(yk,1))
  	ASSERT(nlat==size(yk,2))
!________________________________________
! Compute gradient at the South Pole

  ys=0.
  yc=0.

  do i=1,nlon
    m=mod(i-1+nlon/2,nlon)+1
    d = f(m,2)-f(i,2)		! cross-pole difference
    ys = ys + d * slam(i)
    yc = yc + d * clam(i)
  end do

  xr = -2.*ys*r2dx(1)/nlon
  yr = -2.*yc*r2dx(1)/nlon
!________________________________________
! or

  sf = sum(f(:,2)*slam(:))
  cf = sum(f(:,2)*clam(:))

  rd = 4.*r2dx(1)/nlon
  xr = rd*sf
  yr = rd*cf
!________________________________________
  do i=1,nlon
    xk(i,1) = +clam(i)*xr -slam(i)*yr
    yk(i,1) = +slam(i)*xr +clam(i)*yr
  end do
!________________________________________
! Compute gradient at grid points away from poles

  	! Precompute a row of del f at different latitude lines

  yk(:,2)=f(:,2)-f(:,1)

  	! Compute
  do j=2,nlat-1
    xk(1,j) = (f(2,j)-f(nlon,j)) * r2dx(j)
    do i=2,nlon-1
      xk(i,j) = (f(i+1,j)-f(i-1,j)) * r2dx(j)
    end do
    xk(nlon,j) = (f(1,j)-f(nlon-1,j)) * r2dx(j)

    yk(:,j+1) = f(:,j+1)-f(:,j)		! next-gridpoint difference
    yk(:,j) = yk(:,j)*rdyp(j) + yk(:,j+1)*rdym(j)
  end do
!________________________________________
! Compute gradient at the North Pole.

  ys=0.
  yc=0.

  do i=1,nlon
    m=mod(i-1+nlon/2,nlon)+1
    d = f(m,nlat-1)-f(i,nlat-1)		! cross-pole difference
    ys = ys + d * slam(i)
    yc = yc + d * clam(i)
  end do

  xr = -2.*ys*r2dx(nlat)/nlon
  yr =  2.*yc*r2dx(nlat)/nlon
!________________________________________
! or

  sf = sum(f(:,nlat-1)*slam(:))
  cf = sum(f(:,nlat-1)*clam(:))
  rd = 4.*r2dx(nlat)/nlon
  xr = rd*sf
  yr =-rd*cf
!________________________________________

  do i=1,nlon
    xk(i,nlat) = +clam(i)*xr +slam(i)*yr
    yk(i,nlat) = -slam(i)*xr +clam(i)*yr
  end do

  !! g := grad(f) = xk*ik + yk*jk = xr*ir + yr+jr
  !!
  !! xr = (ir,g)
  !!    = xk*(ir,ik) + yk*(ir,jk)
  !!    = xk*cos(ai) + yk*cos(aj)
  !!    = xk*cos(ai) + yk*cos(ai+pi/2)
  !!    = xk*cos(ai) - yk*sin(ai)
  !!
  !! yr = (jr,g)
  !!    = xk*(jr,ik) + yk*(jr,jk)
  !!	= xk*sin(ai) + yk*sin(aj)
  !!	= xk*sin(ai) + yk*sin(ai+pi/2)
  !!	= xk*sin(ai) + yk*cos(ai)
  !!
  !! /xr\   / cos(ai) -sin(ai)\ /xk\
  !! \yr/ = \ sin(ai)  cos(ai)/ \yk/
  !!
  !! /xk\   / cos(ai)  sin(ai)\ /xr\
  !! \yk/ = \-sin(ai)  cos(ai)/ \yr/

_EXIT_
end subroutine gradlp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gradpl_ - gradiant for latitude-major storage
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gradpl_(slam,clam,r2dx,rdyp,rdym,f,xk,yk)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:),intent(in) :: slam,clam
      real,dimension(:),intent(in) :: r2dx,rdyp,rdym
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: xk,yk

! !REVISION HISTORY:
! 	07Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gradpl_'
  real :: ys,yc, xr,yr, d, rd, sf,cf
  integer :: i,m,j
  integer :: nlon,nlat
_ENTRY_

  nlon=size(f,2)
  nlat=size(f,1)

  	ASSERT(nlon==size(slam))
  	ASSERT(nlon==size(clam))
  	ASSERT(nlat==size(r2dx))
  	ASSERT(nlat==size(rdyp))
  	ASSERT(nlat==size(rdym))

  	ASSERT(nlon==size(xk,2))
  	ASSERT(nlat==size(xk,1))
  	ASSERT(nlon==size(yk,2))
  	ASSERT(nlat==size(yk,1))
!________________________________________
! Compute gradient at the South Pole

  ys=0.
  yc=0.

!!  call shownorm2_(myname_,'f',f)
  do i=1,nlon
    m=mod(i-1+nlon/2,nlon)+1
    d = f(2,m)-f(2,i)		! cross-pole difference
    ys = ys + d * slam(i)
    yc = yc + d * clam(i)
  end do

  xr = -2.*ys*r2dx(1)/nlon
  yr = -2.*yc*r2dx(1)/nlon
!________________________________________
! or

  sf = sum(f(2,:)*slam(:))
  cf = sum(f(2,:)*clam(:))

  rd = 4.*r2dx(1)/nlon
  xr = rd*sf
  yr = rd*cf
!________________________________________
  do i=1,nlon
    xk(1,i) = +clam(i)*xr -slam(i)*yr
    yk(1,i) = +slam(i)*xr +clam(i)*yr
  end do

!!  call shownorm2_(myname_,'xk',xk)
!!  call shownorm2_(myname_,'yk',yk)
!________________________________________
! Compute gradient at grid points away from poles

  	! Precompute a row of del f at different latitude lines

  yk(2,:)=f(2,:)-f(1,:)

  	! Compute
  do j=2,nlat-1
    xk(j,1) = (f(j,2)-f(j,nlon)) * r2dx(j)
    do i=2,nlon-1
      xk(j,i) = (f(j,i+1)-f(j,i-1)) * r2dx(j)
    end do
    xk(j,nlon) = (f(j,1)-f(j,nlon-1)) * r2dx(j)

    yk(j+1,:) = f(j+1,:)-f(j,:)		! next-gridpoint difference
    yk(j  ,:) = yk(j,:)*rdyp(j) + yk(j+1,:)*rdym(j)
  end do
!________________________________________
! Compute gradient at the North Pole.

  ys=0.
  yc=0.

  do i=1,nlon
    m=mod(i-1+nlon/2,nlon)+1
    d = f(nlat-1,m)-f(nlat-1,i)		! cross-pole difference
    ys = ys + d * slam(i)
    yc = yc + d * clam(i)
  end do

  xr = -2.*ys*r2dx(nlat)/nlon
  yr =  2.*yc*r2dx(nlat)/nlon
!________________________________________
! or

  sf = sum(f(nlat-1,:)*slam(:))
  cf = sum(f(nlat-1,:)*clam(:))
  rd = 4.*r2dx(nlat)/nlon
  xr = rd*sf
  yr =-rd*cf
!________________________________________

  do i=1,nlon
    xk(nlat,i) = +clam(i)*xr +slam(i)*yr
    yk(nlat,i) = -slam(i)*xr +clam(i)*yr
  end do

  !! g := grad(f) = xk*ik + yk*jk = xr*ir + yr+jr
  !!
  !! xr = (ir,g)
  !!    = xk*(ir,ik) + yk*(ir,jk)
  !!    = xk*cos(ai) + yk*cos(aj)
  !!    = xk*cos(ai) + yk*cos(ai+pi/2)
  !!    = xk*cos(ai) - yk*sin(ai)
  !!
  !! yr = (jr,g)
  !!    = xk*(jr,ik) + yk*(jr,jk)
  !!	= xk*sin(ai) + yk*sin(aj)
  !!	= xk*sin(ai) + yk*sin(ai+pi/2)
  !!	= xk*sin(ai) + yk*cos(ai)
  !!
  !! /xr\   / cos(ai) -sin(ai)\ /xk\
  !! \yr/ = \ sin(ai)  cos(ai)/ \yk/
  !!
  !! /xk\   / cos(ai)  sin(ai)\ /xr\
  !! \yk/ = \-sin(ai)  cos(ai)/ \yr/

_EXIT_
end subroutine gradpl_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gradr2_ - gradient of a rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gradr2_(ob,f,xk,yk,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: xk,yk
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gradr2_'
  logical :: latlon_
_ENTRY_

  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  if(latlon_) then
    call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,xk,yk)
  else
    call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,xk,yk)
  endif
end subroutine gradr2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: strmr2_ - "stream vector" of a rank-2 array
!
! !DESCRIPTION:
!	(xk,yk) := strm(f) := ek.cross.grad(f) = (-fy,fx)
!
! !INTERFACE:

    subroutine strmr2_(ob,f,xk,yk,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: xk,yk
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::strmr2_'
  logical :: latlon_
_ENTRY_
  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  if(latlon_) then
    call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,yk,xk)
  else
    call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,yk,xk)
  endif
  xk(:,:)=-xk(:,:)
_EXIT_
end subroutine strmr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gradr3_ - gradient of a rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gradr3_(ob,f,xk,yk,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: f
      real,dimension(:,:,:),intent(out) :: xk,yk
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gradr3_'
  integer :: l,n
  logical :: latlon_
_ENTRY_

  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  n=size(f,3)
  	ASSERT(n==size(xk,3))
  	ASSERT(n==size(yk,3))

!!  call shownorm3_(myname_,'f=',f)
  if(latlon_) then
    do l=1,n
      call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),xk(:,:,l),yk(:,:,l)	)
    end do
  else
    do l=1,n
      call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),xk(:,:,l),yk(:,:,l)	)
    end do
  endif
!!  call shownorm3_(myname_,'xk=',xk)
!!  call shownorm3_(myname_,'yk=',yk)
_EXIT_
end subroutine gradr3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: strmr3_ - "stream vector" of a rank-3 array
!
! !DESCRIPTION:
!	(xk,yk) := strm(f) := ek.cross.grad(f) = (-fy,fx)
!
! !INTERFACE:

    subroutine strmr3_(ob,f,xk,yk,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: f
      real,dimension(:,:,:),intent(out) :: xk,yk
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::strmr3_'
  integer :: l,n
  logical :: latlon_
_ENTRY_

  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  n=size(f,3)
  	ASSERT(n==size(xk,3))
  	ASSERT(n==size(yk,3))

  if(latlon_) then
    do l=1,n
      call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),yk(:,:,l),xk(:,:,l)	)
      xk(:,:,l)=-xk(:,:,l)
    end do
  else
    do l=1,n
      call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),yk(:,:,l),xk(:,:,l)	)
      xk(:,:,l)=-xk(:,:,l)
    end do
  endif
_EXIT_
end subroutine strmr3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diver2_ - divergence of rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine diver2_(ob,u,v,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: u,v
      real,dimension(:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diver2_'
_ENTRY_
  call divvor2_(ob,u,v,div=d,latlon=latlon)
_EXIT_
end subroutine diver2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: diver3_ - divergence of rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine diver3_(ob,u,v,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: u,v
      real,dimension(:,:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::diver3_'
_ENTRY_
!!  call shownorm3_(myname_,'u',u)
!!  call shownorm3_(myname_,'v',v)
  call divvor3_(ob,u,v,div=d,latlon=latlon)
!!  call shownorm3_(myname_,'d',d)
_EXIT_
end subroutine diver3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vortr2_ - vortisity of rank-2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vortr2_(ob,u,v,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: u,v
      real,dimension(:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vortr2_'
_ENTRY_
  call divvor2_(ob,u,v,vor=d,latlon=latlon)
_EXIT_
end subroutine vortr2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vortr3_ - vortisity of rank-3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vortr3_(ob,u,v,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: u,v
      real,dimension(:,:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vortr3_'
_ENTRY_
  call divvor3_(ob,u,v,vor=d,latlon=latlon)
_EXIT_
end subroutine vortr3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lapla2_ - solve d = div(grad(f)) for d
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lapla2_(ob,f,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	26May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lapla2_'
  real,dimension(size(f,1),size(f,2)) :: fx,fy
_ENTRY_
  call gradr2_(ob,f,fx,fy,latlon=latlon)
  call diver2_(ob,fx,fy,d,latlon=latlon)
_EXIT_
end subroutine lapla2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lapla3_ - solve d = div(grad(f)) for d
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lapla3_(ob,f,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: f
      real,dimension(:,:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	26May06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lapla3_'
  real,dimension(size(f,1),size(f,2),size(f,3)) :: fx,fy
_ENTRY_
  call gradr3_(ob,f,fx,fy,latlon=latlon)
  call diver3_(ob,fx,fy,d,latlon=latlon)
_EXIT_
end subroutine lapla3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dz2uv2_ - solve d=div(u,v),z=vor(u,v) for u,v, rank-2
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dz2uv2_(ob,div,vor,u,v,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : die,assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: div,vor
      real,dimension(:,:),intent(out) :: u,v
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dz2uv2_'
  integer :: ier,imx,jmx
  real,allocatable,dimension(:,:) :: pot,uwk,vwk
_ENTRY_
  imx=size(div,1)
  jmx=size(div,2)
  	ASSERT(imx==size(vor,1))
  	ASSERT(imx==size(  u,1))
  	ASSERT(imx==size(  v,1))
  	ASSERT(jmx==size(vor,2))
  	ASSERT(jmx==size(  u,2))
  	ASSERT(jmx==size(  v,2))

  	allocate(pot(imx,jmx))
  	allocate(uwk(imx,jmx),vwk(imx,jmx))

  	! invert div(u,v) by solving div=lapla(vpot) then grad(vpot).

  call ellip2_(ob,div,pot,stat=ier,latlon=latlon)
  	if(ier/=0) call die(myname_,'ellip2_(div)',ier)
  call gradr2_(ob,pot,u,v,latlon=latlon)

  	! invert vor(u,v) by solving vor=lapla(sfun) then strm(sfun).

  call ellip2_(ob,vor,pot,stat=ier,latlon=latlon)
  	if(ier/=0) call die(myname_,'ellip2_(vor)',ier)
  call strmr2_(ob,pot,uwk,vwk,latlon=latlon)
  u=u+uwk; v=v+vwk

  	deallocate(pot)
	deallocate(uwk,vwk)
_EXIT_
end subroutine dz2uv2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dz2uv3_ - solve d=div(u,v),z=vor(u,v) for u,v, rank-3
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dz2uv3_(ob,div,vor,u,v,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : die,assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: div,vor
      real,dimension(:,:,:),intent(out) :: u,v
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dz2uv3_'
  integer :: ier,imx,jmx,kmx,k
  real,allocatable,dimension(:,:) :: pot,uwk,vwk
_ENTRY_
  imx=size(div,1)
  jmx=size(div,2)
  kmx=size(div,3)
  	ASSERT(imx==size(vor,1))
  	ASSERT(imx==size(  u,1))
  	ASSERT(imx==size(  v,1))
  	ASSERT(jmx==size(vor,2))
  	ASSERT(jmx==size(  u,2))
  	ASSERT(jmx==size(  v,2))
  	ASSERT(jmx==size(vor,3))
  	ASSERT(jmx==size(  u,3))
  	ASSERT(jmx==size(  v,3))

  	allocate(pot(imx,jmx))
  	allocate(uwk(imx,jmx),vwk(imx,jmx))

  do k=1,kmx

  	! invert div(u,v) by solving div=lapla(vpot) then grad(vpot).

    call ellip2_(ob,div(:,:,k),pot,stat=ier,latlon=latlon)
  	if(ier/=0) call die(myname_,'ellip2_(div)',ier)
    call gradr2_(ob,pot,u(:,:,k),v(:,:,k),latlon=latlon)

  	! invert vor(u,v) by solving vor=lapla(sfun) then strm(sfun).

    call ellip2_(ob,vor(:,:,k),pot,stat=ier,latlon=latlon)
  	if(ier/=0) call die(myname_,'ellip2_(vor)',ier)
    call strmr2_(ob,pot,uwk,vwk,latlon=latlon)
    u(:,:,k)=u(:,:,k)+uwk(:,:)
    v(:,:,k)=v(:,:,k)+vwk(:,:)
  end do

  	deallocate(pot)
	deallocate(uwk,vwk)
_EXIT_
end subroutine dz2uv3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ellip2_ - 2-D elliptic equation solver
!
! !DESCRIPTION:
!	solve Lapla(f) := div(grad(f)) = d for f
!
! !INTERFACE:

    subroutine ellip2_(ob,d,f,stat,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none

      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: d
      real,dimension(:,:),intent(out) :: f
      integer,optional,intent(out) :: stat
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- fixed div(u,v) calculation bugs
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ellip2_'
  integer :: imx,jmx,imp,iwk,m,n,idimf,ier
  logical :: latlon_
  real :: ONEPI,TWOPI, radsq, pertrb,fmean

  real,allocatable,dimension(:) :: vp,w
  real,allocatable,dimension(:) :: bdts,bdtf
  real,allocatable,dimension(:) :: bdps,bdpf
_ENTRY_
	ASSERT(ob%evenlat)

  if(present(stat)) stat=0
  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  ONEPI=4.*atan(1.)
  TWOPI=8.*atan(1.)
  radsq=ob%radius*ob%radius

  imx=size(d,1)
  jmx=size(d,2)
  	ASSERT(imx==size(f,1))
  	ASSERT(jmx==size(f,2))

  imp=imx+1
  iwk=11*jmx+6*imp
  m=jmx-1
  n=imx
  idimf=jmx

  allocate(vp(jmx*imp), w(iwk))
  allocate(bdts(imp),bdtf(imp))
  allocate(bdps(jmx),bdpf(jmx))

  if(latlon_) then
    call packpl_(d,vp)		! d(lat,lon) => vp(lat*(lon+1))
  else
    call packlp_(d,vp)		! d(lon,lat) => vp(lat*(lon+1))
  endif

  	! boundary conditions
  bdts(:)=0.
  bdtf(:)=0.
  bdps(:)=0.
  bdpf(:)=0.

  call PWSSSP(  0, 0.,ONEPI, jmx-1,9,bdts,bdtf,	&
		   0.,TWOPI, imx  ,0,bdps,bdpf,	&
	           0.,   vp, jmx  ,pertrb,ier , w)
	if(ier/=0) then
	  call myperr_(myname_,'PWSSSP()',ier,pertrb)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(latlon_) then
    call unpackpl_(vp,f)		! vp(lat*(lon+1)) -> f(lat,lon)
    fmean=areameanpl_(f,ob%cosphi)
  else
    call unpacklp_(vp,f)		! vp(lat*(lon+1)) -> f(lon,lat)
    fmean=areameanlp_(f,ob%cosphi)
  endif
  	! remove the mean and rescale to the physical dimension of the
	! sphere.
  f(:,:)=(f(:,:)-fmean)*radsq
_EXIT_
end subroutine ellip2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ellip3_ - 3-D elliptic equation solver
!
! !DESCRIPTION:
!	solve Lapla(f) := div(grad(f)) = d for f
!
! !INTERFACE:

    subroutine ellip3_(ob,d,f,stat,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none

      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: d
      real,dimension(:,:,:),intent(out) :: f
      integer,optional,intent(out) :: stat
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- fixed div(u,v) calculation bugs
! 	20Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ellip3_'
  integer :: imx,jmx,kmx,k,imp,iwk,m,n,idimf,ier
  logical :: latlon_
  real :: ONEPI,TWOPI, radsq, pertrb, fmean

  real,allocatable,dimension(:) :: vp,w
  real,allocatable,dimension(:) :: bdts,bdtf
  real,allocatable,dimension(:) :: bdps,bdpf
  integer :: INTL
_ENTRY_
  if(present(stat)) stat=0
  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  ONEPI=4.*atan(1.)
  TWOPI=8.*atan(1.)
  radsq=ob%radius*ob%radius

  imx=size(d,1)
  jmx=size(d,2)
  kmx=size(d,3)
  	ASSERT(imx==size(f,1))
  	ASSERT(jmx==size(f,2))
  	ASSERT(kmx==size(f,3))

  if(latlon_) then	! swap imx and jmx
    iwk=imx
    imx=jmx
    jmx=iwk
  endif
  	! at this point, imx is for lon and jmx is for lat

  imp=imx+1
  iwk=11*jmx+6*imp
  m=jmx-1
  n=imx
  idimf=jmx

  allocate(vp(jmx*imp), w(iwk))
  allocate(bdts(imp),bdtf(imp))
  allocate(bdps(jmx),bdpf(jmx))

  	! boundary conditions
  bdts(:)=0.
  bdtf(:)=0.
  bdps(:)=0.
  bdpf(:)=0.

  	! loop over all levels
  INTL=0
  do k=1,kmx
    if(latlon_) then
      call packpl_(d(:,:,k),vp)		! d(lat,lon) => vp(lat*(lon+1))
    else
      call packlp_(d(:,:,k),vp)		! d(lon,lat) => vp(lat*(lon+1))
    endif

    call PWSSSP(INTL,0.,ONEPI, jmx-1,9,bdts,bdtf,	&
		     0.,TWOPI, imx  ,0,bdps,bdpf,	&
	             0.,   vp, jmx  ,pertrb,ier , w)
	if(ier/=0) then
	  call myperr_(myname_,'PWSSSP()',ier,pertrb,level=k)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    INTL=1

    if(latlon_) then
      call unpackpl_(vp,f(:,:,k))	! vp(lat*(lon+1)) -> f(lat,lon)
      fmean=areameanpl_(f(:,:,k),ob%cosphi)
    else
      call unpacklp_(vp,f(:,:,k))	! vp(lat*(lon+1)) -> f(lon,lat)
      fmean=areameanlp_(f(:,:,k),ob%cosphi)
    endif

  	! remove the mean and rescale to the physical dimension of the
	! sphere.
    f(:,:,k)=(f(:,:,k)-fmean)*radsq
  end do
_EXIT_
end subroutine ellip3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: packlp_ - pack 2-d d(lon,lat) -> 1-d w(lat*(lon+1))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine packlp_(d,w)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:,:),intent(in) :: d
      real,dimension(:),intent(out) :: w

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::packlp_'
  integer :: ij,i
  integer :: imx,jmx
_ENTRY_
  imx=size(d,1)
  jmx=size(d,2)

  	ASSERT((imx+1)*jmx==size(w))

    	! Copy from a (lon,lat) array to a lat*(lon+1) flat storage with
	! an extra band of imx+1 for periodic boundary condition.

  ij=0
  do i=1,imx	! for each longitude, copy a slice of all latitudes
    w(ij+1:ij+jmx)=d(i,1:jmx)
    ij=ij+jmx
  end do
  w(ij+1:ij+jmx)=d(1,1:jmx)	! periodic boundary condition
_EXIT_
end subroutine packlp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: packpl - pack 2-d d(lat,lon) -> 1-d w(lat*(lon+1))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine packpl_(d,w)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:,:),intent(in) :: d
      real,dimension(:),intent(out) :: w

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::packpl_'
  integer :: ij,i
  integer :: imx,jmx
_ENTRY_
  imx=size(d,2)		! lon
  jmx=size(d,1)		! lat

  	ASSERT((imx+1)*jmx==size(w))

    	! Copy from a (lat,lon) array to a lat*(lon+1) flat storage with
	! an extra band of imx+1 for periodic boundary condition.

  ij=0
  do i=1,imx	! for each longitude, copy a slice of all latitudes
    w(ij+1:ij+jmx)=d(1:jmx,i)
    ij=ij+jmx
  end do
  w(ij+1:ij+jmx)=d(1:jmx,1)	! periodic boundary condition
_EXIT_
end subroutine packpl_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpacklp_ - unpack 2-d d(lon,lat) -> 1-d w(lat*(lon+1))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpacklp_(w,d)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:),intent(in) :: w
      real,dimension(:,:),intent(out) :: d

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpacklp_'
  integer :: ij,j
  integer :: imx,jmx
_ENTRY_
  imx=size(d,1)		! lon
  jmx=size(d,2)		! lat

  	ASSERT((imx+1)*jmx==size(w))

    	! Copy to a (lon,lat) array from a lat*(lon+1) flat storage with
	! an extra band of imx+1 for periodic boundary condition.

  ij=(imx-1)*jmx
  do j=1,jmx	! for each latitude, copy a slice of all longitudes
    d(1:imx,j)=w(j:j+ij:jmx)
  end do
_EXIT_
end subroutine unpacklp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpackpl - unpack 2-d d(lat,lon) -> 1-d w(lat*(lon+1))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpackpl_(w,d)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real,dimension(:),intent(in) :: w
      real,dimension(:,:),intent(out) :: d

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpackpl_'
  integer :: ij,j
  integer :: imx,jmx
_ENTRY_
  imx=size(d,2)		! lon
  jmx=size(d,1)		! lat

  	ASSERT((imx+1)*jmx==size(w))

    	! Copy to a (lat,lon) array from a lat*(lon+1) flat storage with
	! an extra band of imx+1 for periodic boundary condition.

  ij=(imx-1)*jmx
  do j=1,jmx	! for each latitude, copy a slice of all longitudes
    d(j,1:imx)=w(j:j+ij:jmx)
  end do
_EXIT_
end subroutine unpackpl_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: areameanlp_ - area-mean of lon-lat array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function areameanlp_(f,clats)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real :: areameanlp_
      real,dimension(:,:),intent(in) :: f	! f(lon,lat)
      real,dimension(:),intent(in) :: clats	! cos(lats)

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::areameanlp_'
  integer :: imx,jmx,j
  real :: dxdy,vmean,amean
_ENTRY_
  imx=size(f,1)
  jmx=size(f,2)
  	ASSERT(jmx>2)
	ASSERT(jmx==size(clats))

  	! This area-weighted mean is based on a simple implementation
	! of an integration form of
	!
	!   intg(f(lam,phi)cos(phi)d(lam)d(phi); lam=0,360; phi=-90,+90)
	!
	! Any numerical error in the same order of the grid mesh
	! resolution, thus, is ignored, which leads to the fact that
	! the values the two poles are not relevent to the integration
	! result.  However, to make the function more general, the
	! values of the latitude grid are used directly, regardless if
	! at any given point, contributions to the total are vanished.
	!
	! Therefore, it is expected by this function, jmx>>2 and
	! grid (1:imx) represents a periodic grid.

  vmean=0.
  amean=0.
  do j=1,jmx
    dxdy=clats(j)		! one at all other points
    if(j==1.or.j==jmx) dxdy=.5*dxdy ! half at the boundaries.
    vmean=vmean+sum(f(1:imx,j))*dxdy
    amean=amean+imx*dxdy
  end do
  	! amean coule be zero if jmx <= 2
  areameanlp_=vmean/amean
_EXIT_
end function areameanlp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: areameanpl_ - area-mean of lat-lon array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function areameanpl_(f,clats)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      real :: areameanpl_
      real,dimension(:,:),intent(in) :: f	! f(lat,lon)
      real,dimension(:),intent(in) :: clats	! cos(lats)

! !REVISION HISTORY:
! 	03Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::areameanpl_'
  integer :: imx,jmx,j
  real :: dxdy,vmean,amean
_ENTRY_
  imx=size(f,2)
  jmx=size(f,1)
  	ASSERT(jmx>2)
	ASSERT(jmx==size(clats))

  	! This area-weighted mean is based on a simple implementation
	! of an integration form of
	!
	!   intg(f(lam,phi)cos(phi)d(lam)d(phi); lam=0,360; phi=-90,+90)
	!
	! Any numerical error in the same order of the grid mesh
	! resolution, thus, is ignored, which leads to the fact that
	! the values the two poles are not relevent to the integration
	! result.  However, to make the function more general, the
	! values of the latitude grid are used directly, regardless if
	! at any given point, contributions to the total are vanished.
	!
	! Therefore, it is expected by this function, jmx>>2 and
	! grid (1:imx) represents a periodic grid.

  vmean=0.
  amean=0.
  do j=1,jmx
    dxdy=clats(j)		! one at all other points
    if(j==1.or.j==jmx) dxdy=.5*dxdy ! half at the boundaries
    vmean=vmean+sum(f(j,1:imx))*dxdy
    amean=amean+imx*dxdy
  end do
  	! amean coule be zero if jmx <= 2
  areameanpl_=vmean/amean
_EXIT_
end function areameanpl_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: laplar2_ - Laplasian of a rank-2 array, div(grad(f))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine laplar2_(ob,f,d,latlon)
      use m_mpout,only : mpout_log
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:),intent(in) :: f
      real,dimension(:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::laplar2_'
  logical :: latlon_
  real,dimension(size(f,1),size(f,2)) :: xk,yk
_ENTRY_
  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  if(latlon_) then
    call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,xk,yk)
    call divrpl_(xk,yk,d,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
  else
    call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,f,xk,yk)
    call divrlp_(xk,yk,d,ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
  endif
_EXIT_
end subroutine laplar2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: laplar3_ - Laplasian of a rank-3 array, div(grad(f))
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine laplar3_(ob,f,d,latlon)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(ggGradientSP),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: f
      real,dimension(:,:,:),intent(out) :: d
      logical,optional,intent(in) :: latlon

! !REVISION HISTORY:
! 	11Jan05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::laplar3_'
  integer :: l,n
  logical :: latlon_
  real,dimension(size(f,1),size(f,2)) :: xk,yk
_ENTRY_
  latlon_=.false.
  if(present(latlon)) latlon_=latlon

  n=size(f,3)
  	ASSERT(n==size(d,3))

  if(latlon_) then
    do l=1,n
      call gradpl_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),xk,yk	)
      call divrpl_(xk,yk,d(:,:,l),ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    end do
  else
    do l=1,n
      call gradlp_(ob%sinlam,ob%coslam,ob%r2dx,ob%rdyp,ob%rdym,	&
      				f(:,:,l),xk,yk	)
      call divrlp_(xk,yk,d(:,:,l),ob%r2dx,ob%cosphi,ob%rdyp,ob%rdym)
    end do
  endif
_EXIT_
end subroutine laplar3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: myperr_ - generate module specific error messages
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine myperr_(name,mesg,ier,err,level)
      use m_mpout,only : mpout_log
      use m_die,only : perr
      implicit none
      character(len=*),intent(in) :: name
      character(len=*),intent(in) :: mesg
      integer,intent(in) :: ier
      real,intent(in) :: err
      integer,optional,intent(in) :: level


! !REVISION HISTORY:
! 	05Jan06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	17Aug2007 - Todling, implemented verbose
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::myperr_'
  character(len=16) :: cerr
  call perr(name,mesg,ier)
  write(cerr,'(1p,e16.8)') err
  call perr(name,'error size ='//trim(adjustl(cerr)))
  if(present(level)) call perr(name,'level =',level)
end subroutine myperr_
subroutine shownorm2_(name,w,f)
  implicit none
  character(len=*),intent(in) :: name
  character(len=*),intent(in) :: w
  real,dimension(:,:),intent(in) :: f
  real :: s
  integer :: j
  if(.not.verbose_)return
  s=0.
  do j=1,size(f,2)
    s=s+dot_product(f(:,j),f(:,j))
  end do
  print*,name,': norm('//w//') =',sqrt(s)
end subroutine shownorm2_
subroutine shownorm3_(name,w,f)
  implicit none
  character(len=*),intent(in) :: name
  character(len=*),intent(in) :: w
  real,dimension(:,:,:),intent(in) :: f
  real :: s
  integer :: j,k
  if(.not.verbose_)return
  s=0.
  do k=1,size(f,3)
  do j=1,size(f,2)
    s=s+dot_product(f(:,j,k),f(:,j,k))
  end do
  end do
  print*,name,': norm('//w//') =',sqrt(s)
end subroutine shownorm3_

end module m_ggGradientSP
