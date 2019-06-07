!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_daInterp - interface of d-a/a-d grid transformation
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_daInterp
      implicit none
      private	! except

      public :: daInterp_vdtoa
      public :: daInterp_vatod
      public :: daInterp_vatod_ad
      public :: daInterp_vdtoa_ad

    interface daInterp_vdtoa; module procedure	&
	vdtoa2d_,	&
  	vdtoa3d_; end interface daInterp_vdtoa

    interface daInterp_vatod; module procedure	&
	vatod2d_,	&
  	vatod3d_; end interface daInterp_vatod

    interface daInterp_vatod_ad; module procedure	&
	vatod2d_ad_,	&
  	vatod3d_ad_; end interface daInterp_vatod_ad

    interface daInterp_vdtoa_ad; module procedure	&
	vdtoa2d_ad_,	&
  	vdtoa3d_ad_; end interface daInterp_vdtoa_ad	

! !REVISION HISTORY:
! 	09Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       01Mar06 - Yanqiu Zhu   add adjoint codes
!	14Jun07 - Todling  add adjoint of dtoa
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_daInterp'

  integer,parameter :: UTYPE=2
  integer,parameter :: VTYPE=1

#include "assert.H"
#ifndef _VECATOD_
#define _VECATOD_
#define _VECDTOA_
#endif
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtoa3d_ - d-grid to a-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtoa3d_(ud,vd,ua,va)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:,:),intent(in ) :: ud,vd	! d-grid wind
      real,dimension(:,:,:),intent(out) :: ua,va	! a-grid wind

! !REVISION HISTORY:
! 	10Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtoa3d_'
  integer :: im,jm,lm,l

  im=size(ud,1)
  jm=size(ud,2)
  lm=size(ud,3)

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))
  	ASSERT(lm==size(vd,3))

  	ASSERT(im==size(ua,1))
  	ASSERT(jm==size(ua,2))
  	ASSERT(lm==size(ua,3))

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))
  	ASSERT(lm==size(va,3))

#ifdef _VECDTOA_
  do l=1,lm
    call vecdtoa_(ud(:,:,l),vd(:,:,l),ua(:,:,l),va(:,:,l))
  end do
#else
  call die(myname_,'not in use for this version, dtoa()')
  call dtoa(ud,ua,im,jm,lm,UTYPE)
  call dtoa(vd,va,im,jm,lm,VTYPE)
#endif

end subroutine vdtoa3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtoa3d_ad_ - adjoint of d-grid to a-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtoa3d_ad_(ud,vd,ua,va)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:,:),intent(inout) :: ud,vd	! d-grid wind
      real,dimension(:,:,:),intent(inout) :: ua,va	! a-grid wind

! !REVISION HISTORY:
! 	14Jun2007 Todling Initial code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtoa3d_ad_'
  integer :: im,jm,lm,l

  im=size(ud,1)
  jm=size(ud,2)
  lm=size(ud,3)

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))
  	ASSERT(lm==size(vd,3))

  	ASSERT(im==size(ua,1))
  	ASSERT(jm==size(ua,2))
  	ASSERT(lm==size(ua,3))

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))
  	ASSERT(lm==size(va,3))

#ifdef _VECDTOA_
  do l=1,lm
    call vecdtoa_ad_(ud(:,:,l),vd(:,:,l),ua(:,:,l),va(:,:,l))
  end do
#else
  call die(myname_,'not in use for this version, dtoa_ad_()')
  call dtoa_ad(ud,ua,im,jm,lm,UTYPE)
  call dtoa_ad(vd,va,im,jm,lm,VTYPE)
#endif

end subroutine vdtoa3d_ad_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtoa2d_ - d-grid to a-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtoa2d_(ud,vd,ua,va)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:),intent(in ) :: ud,vd	! d-grid wind
      real,dimension(:,:),intent(out) :: ua,va	! a-grid wind

! !REVISION HISTORY:
! 	10Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtoa2d_'
  integer :: im,jm

  im=size(ud,1)
  jm=size(ud,2)

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))

  	ASSERT(im==size(ua,1))
  	ASSERT(jm==size(ua,2))

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))

#ifdef _VECDTOA_
  call vecdtoa_(ud,vd,ua,va)
#else
  call die(myname_,'not in use for this version, dtoa()')
  call dtoa(ud,ua,im,jm,1,UTYPE)
  call dtoa(vd,va,im,jm,1,VTYPE)
#endif

end subroutine vdtoa2d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtoa2d_ad_ - adjoint of d-grid to a-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtoa2d_ad_(ud,vd,ua,va)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:),intent(inout) :: ud,vd	! d-grid wind
      real,dimension(:,:),intent(inout) :: ua,va	! a-grid wind

! !REVISION HISTORY:
! 	14Jun2007  Todling  Initial code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtoa2d_ad_'
  integer :: im,jm

  im=size(ud,1)
  jm=size(ud,2)

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))

  	ASSERT(im==size(ua,1))
  	ASSERT(jm==size(ua,2))

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))

#ifdef _VECDTOA_
  call vecdtoa_ad_(ud,vd,ua,va)
#else
  call die(myname_,'not in use for this version, dtoa_ad()')
  call dtoa_ad(ud,ua,im,jm,1,UTYPE)
  call dtoa_ad(vd,va,im,jm,1,VTYPE)
#endif

end subroutine vdtoa2d_ad_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vatod3d_ - a-grid to d-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vatod3d_(ua,va,ud,vd,umiss,uskip)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:,:),intent(in ) :: ua,va	! a-grid wind
      real,dimension(:,:,:),intent(out) :: ud,vd	! d-grid wind
      real,optional,intent(in) :: umiss ! "missing-value" of u at j=1
      logical,optional,intent(in) :: uskip ! leave u at j=1 alone

! !REVISION HISTORY:
! 	10Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vatod3d_'
  integer :: im,jm,lm,l

  im=size(ua,1)
  jm=size(ua,2)
  lm=size(ua,3)

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))
  	ASSERT(lm==size(va,3))

  	ASSERT(im==size(ud,1))
  	ASSERT(jm==size(ud,2))
  	ASSERT(lm==size(ud,3))

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))
  	ASSERT(lm==size(vd,3))

#ifdef _VECATOD_
  do l=1,lm
    call vecatod_(ua(:,:,l),va(:,:,l),ud(:,:,l),vd(:,:,l), &
    	umiss=umiss,uskip=uskip)
  end do
#else
  call die(myname_,'not in use for this version, atod()')
  call atod(ua,ud,im,jm,lm,UTYPE)
  call atod(va,vd,im,jm,lm,VTYPE)
#endif

end subroutine vatod3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vatod3d_ad_ - adjoint of a-grid to d-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vatod3d_ad_(ua,va,ud,vd)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:,:),intent(inout ) :: ua,va        ! a-grid wind
      real,dimension(:,:,:),intent(in) :: ud,vd            ! d-grid wind

! !REVISION HISTORY:
!       02Mar06 - Yanqiu Zhu
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vatod3d_ad_'
  integer :: im,jm,lm,l

  im=size(ua,1)
  jm=size(ua,2)
  lm=size(ua,3)

        ASSERT(im==size(va,1))
        ASSERT(jm==size(va,2))
        ASSERT(lm==size(va,3))

        ASSERT(im==size(ud,1))
        ASSERT(jm==size(ud,2))
        ASSERT(lm==size(ud,3))

        ASSERT(im==size(vd,1))
        ASSERT(jm==size(vd,2))
        ASSERT(lm==size(vd,3))

#ifdef _VECATOD_
  do l=1,lm
    call vecatod_ad_(ua(:,:,l),va(:,:,l),ud(:,:,l),vd(:,:,l))
  end do
#else
  call die(myname_,'not in use for this version, atod_ad()')
  call atod_ad(ua,ud,im,jm,lm,UTYPE)
  call atod_ad(va,vd,im,jm,lm,VTYPE)
#endif

end subroutine vatod3d_ad_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vatod2d_ - d-grid to a-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vatod2d_(ua,va,ud,vd,umiss,uskip)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:),intent(in ) :: ua,va	! a-grid wind
      real,dimension(:,:),intent(out) :: ud,vd	! d-grid wind
      real,optional,intent(in) :: umiss ! "missing-value" of u at j=1
      logical,optional,intent(in) :: uskip ! leave u at j=1 alone

! !REVISION HISTORY:
! 	10Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vatod2d_'
  integer :: im,jm

  im=size(ua,1)
  jm=size(ua,2)

  	ASSERT(im==size(va,1))
  	ASSERT(jm==size(va,2))

  	ASSERT(im==size(ud,1))
  	ASSERT(jm==size(ud,2))

  	ASSERT(im==size(vd,1))
  	ASSERT(jm==size(vd,2))

#ifdef _VECATOD_
  call vecatod_(ua,va,ud,vd,umiss=umiss,uskip=uskip)
#else
  call die(myname_,'not in use for this version, atod()')
  call atod(ua,ud,im,jm,1,UTYPE)
  call atod(va,vd,im,jm,1,VTYPE)
#endif
end subroutine vatod2d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vatod2d_ad_ - adjoint of a-grid to d-grid vector interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vatod2d_ad_(ua,va,ud,vd)
      use m_die,only : assert_,die
      implicit none
      real,dimension(:,:),intent(inout ) :: ua,va  ! a-grid wind
      real,dimension(:,:),intent(in) :: ud,vd  ! d-grid wind

! !REVISION HISTORY:
!       02Mar06 - Yanqiu Zhu
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vatod2d_ad_'
  integer :: im,jm

  im=size(ua,1)
  jm=size(ua,2)

        ASSERT(im==size(va,1))
        ASSERT(jm==size(va,2))

        ASSERT(im==size(ud,1))
        ASSERT(jm==size(ud,2))

        ASSERT(im==size(vd,1))
        ASSERT(jm==size(vd,2))

#ifdef _VECATOD_
  call vecatod_ad_(ua,va,ud,vd)
#else
  call die(myname_,'not in use for this version, atod_ad()')
  call atod_ad(ua,ud,im,jm,1,UTYPE)
  call atod_ad(va,vd,im,jm,1,VTYPE)
#endif
end subroutine vatod2d_ad_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vecatod_ - interpolate a-grid wind to d-grid wind
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vecatod_(ua,va,ud,vd,umiss,uskip)
      implicit none
      real,dimension(:,:),intent(in) :: ua,va
      real,dimension(:,:),intent(out) :: ud,vd
      real,optional,intent(in) :: umiss ! "missing-value" of u at j=1
      logical,optional,intent(in) :: uskip	! leave u at j=1 alone

! !REVISION HISTORY:
! 	15Mar06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	16Feb07	- Larry Takacs <takacs@gmao.gsfc.nasa.gov>
!		- Fixed the storage of vd.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vecatod_'
  integer :: nlon,j,nlat
  real :: wx,wy
  logical :: uskip_
#ifdef ADJOINT_CHECK_vecatod_
  real :: sum1,sum2
  real,dimension(size(ua,1),size(ua,2)) :: a_ua,a_va
#endif

  uskip_=.false.
  if(present(uskip)) uskip_=uskip

  nlon=size(ua,1)
  nlat=size(ua,2)

  ! (ua,va) -> vd at the South Pole.  ud is not touched by setuv_()
  call getxy_(wx,wy,ua(:,1),va(:,1),stag=0.,pole=-1.)
  	! vd(i) stores data at a loccation a half interval to the left
	! of the i-th grid point (-0.5 dlon).
  call setuv_(wx,wy,ud(:,1),vd(:,1),stag=-.5,pole=-1.,skipuv=-1)
  if(.not.uskip_) then
    ud(:,1)=0.
    if(present(umiss)) ud(:,1)=umiss
  endif

  do j=2,nlat-1
    ud(:,j)=.5*(ua(:,j-1)+ua(:,j))
  end do
  do j=2,nlat-1
	! The vd indices agree with the upper corners of va box.
    vd(1       ,j)=.5*(va(  nlon  ,j)+va(1     ,j))
    vd(2:nlon  ,j)=.5*(va(1:nlon-1,j)+va(2:nlon,j))
  end do

  ! (ua,va) -> vd at the North Pole.  ud is not touched by setuv_()
  call getxy_(wx,wy,ua(:,nlat),va(:,nlat),stag=0.,pole=+1.)
  	! vd(i) stores data at a loccation a half interval to the left
	! of the i-th grid point (-0.5 dlon).
  call setuv_(wx,wy,ud(:,nlat),vd(:,nlat),stag=-.5,pole=+1.,skipuv=-1)
  ud(:,nlat)=.5*(ua(:,nlat-1)+ua(:,nlat))

#ifdef ADJOINT_CHECK_vecatod_
  a_ua(:,:)=0.
  a_va(:,:)=0.
  call vecatod_ad_(a_ua,a_va,ud,vd)
  sum1=0.
  sum2=0.
  do j=1,nlat
    sum1=sum1+dot_product(a_ua(:,j),ua(:,j))
    sum1=sum1+dot_product(a_va(:,j),va(:,j))
    if(j>1) sum2=sum2+dot_product(ud(:,j),ud(:,j))
    sum2=sum2+dot_product(vd(:,j),vd(:,j))
  enddo
  print*,myname_,': sum1=',sum1,'  sum2=',sum2
#endif

end subroutine vecatod_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vecatod_ad_ - adjoint operator for interpolating a-grid wind to d-grid wind
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vecatod_ad_(ua,va,ud,vd)
      implicit none
      real,dimension(:,:),intent(inout) :: ua,va
      real,dimension(:,:) :: ud,vd

! !REVISION HISTORY:
! 	17Mar06	- Yanqiu Zhu
!	
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vecatod_ad_'
  integer :: nlon,j,nlat
  real :: wx,wy

  nlon=size(ua,1)
  nlat=size(ua,2)

  ! (ua,va) -> vd at the North Pole.  ud is not touched by setuv_()
  ua(:,nlat-1) = ua(:,nlat-1) + 0.5*ud(:,nlat)
  ua(:,nlat)   = ua(:,nlat) + 0.5*ud(:,nlat)
  wx = 0.0
  wy = 0.0
  call setuv_ad_(wx,wy,ud(:,nlat),vd(:,nlat),stag=.5,pole=+1.,skipuv=-1)
  call getxy_ad_(wx,wy,ua(:,nlat),va(:,nlat),stag=0.,pole=+1.)

  do j=nlat-1,2,-1
     va(1   ,j) = va(1   ,j) + 0.5*vd(nlon,j)
     va(nlon,j) = va(nlon,j) + 0.5*vd(nlon,j)
     va(2:nlon,j)   = va(2:nlon,j) + 0.5*vd(1:nlon-1,j)
     va(1:nlon-1,j) = va(1:nlon-1,j) + 0.5*vd(1:nlon-1,j)
  end do
  do j=nlat-1,2,-1
     ua(:,j-1) = ua(:,j-1) + 0.5*ud(:,j)
     ua(:,j)   = ua(:,j) + 0.5*ud(:,j)
  end do

  ! (ua,va) -> vd at the South Pole.  ud is not touched by setuv_()
  ud(:,1)=0.
  wx = 0.0
  wy = 0.0
  call setuv_ad_(wx,wy,ud(:,1),vd(:,1),stag=.5,pole=-1.,skipuv=-1)
  call getxy_ad_(wx,wy,ua(:,1),va(:,1),stag=0.,pole=-1.)

end subroutine vecatod_ad_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vecdtoa_ - interpolate d-grid wind to a-grid wind
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vecdtoa_(ud,vd,ua,va)
      implicit none
      real,dimension(:,:),intent(in ) :: ud,vd
      real,dimension(:,:),intent(out) :: ua,va

! !REVISION HISTORY:
! 	15Mar06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	16Feb07	- Larry Takacs <takacs@gmao.gsfc.nasa.gov>
!		- Fixed the storage of vd.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vecdtoa_'
  integer :: nlon,j,nlat
  real :: wx,wy

  nlon=size(ud,1)
  nlat=size(ud,2)

  ! vd -> (ua,va) at the South Pole.  ud is not used by getxy_()
  	! vd(i) stores data at a loccation a half interval to the left
	! of the i-th grid point (-0.5 dlon).
  call getxy_(wx,wy,ud(:,1),vd(:,1),stag=-.5,pole=-1.,skipuv=-1)
  call setuv_(wx,wy,ua(:,1),va(:,1),stag=0.,pole=-1.)

  do j=2,nlat-1
    ua(:,j)=.5*(ud(:,j)+ud(:,j+1))
  end do

  do j=2,nlat-1
	! The va locations agree with the lower corners of vd box.
    va(1:nlon-1,j)=.5*(vd(1:nlon-1,j)+vd(2:nlon,j))
    va(  nlon  ,j)=.5*(vd(  nlon  ,j)+vd(1     ,j))
  end do

  ! vd -> (ua,va) at the North Pole.  ud is not used by getxy_()
  	! vd(i) stores data at a loccation a half interval to the left
	! of the i-th grid point (-0.5 dlon).
  call getxy_(wx,wy,ud(:,nlat),vd(:,nlat),stag=-.5,pole=+1.,skipuv=-1)
  call setuv_(wx,wy,ua(:,nlat),va(:,nlat),stag=0.,pole=+1.)
end subroutine vecdtoa_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vecdtoa_ad_ - adjoint of interpolate d-grid wind to a-grid wind
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vecdtoa_ad_(ud,vd,ua,va)
      implicit none
      real,dimension(:,:),intent(inout) :: ud,vd
      real,dimension(:,:),intent(inout) :: ua,va

! !REVISION HISTORY:
! 	14Jun2007 Todling  Initial code; followed YZhu's prescription
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vecdtoa_ad_'
  integer :: nlon,j,nlat
  real :: wx,wy

  nlon=size(ud,1)
  nlat=size(ud,2)

  wx = 0.0
  wy = 0.0
  call setuv_ad_(wx,wy,ua(:,nlat),va(:,nlat),stag=0.,pole=+1.)
  call getxy_ad_(wx,wy,ud(:,nlat),vd(:,nlat),stag=-.5,pole=+1.,skipuv=-1)
  
  do j=nlat-1,2,-1
    vd(    nlon,j)= vd(  nlon  ,j) + .5*va(     1,j)
    vd(       1,j)= vd(     1  ,j) + .5*va(     1,j)
    vd(1:nlon-1,j)= vd(1:nlon-1,j) + .5*va(2:nlon,j)
    vd(2:nlon  ,j)= vd(2:nlon  ,j) + .5*va(2:nlon,j)
  end do
  
  do j=nlat-1,2,-1
    ud(:,j+1) = ud(:,j+1) + .5*ua(:,j)
    ud(:,j)   = ud(:,j)   + .5*ua(:,j)
  end do
  
  call setuv_ad_(wx,wy,ua(:,1),va(:,1),stag=0.,pole=-1.)
  call getxy_ad_(wx,wy,ud(:,1),vd(:,1),stag=-.5,pole=-1.,skipuv=-1)
  
end subroutine vecdtoa_ad_

subroutine getxy_(wx,wy,u,v,stag,pole,reflon,skipuv)
  implicit none
  real,intent(out) :: wx,wy
  real,dimension(:),intent(in) :: u,v
  real,intent(in) :: stag
  real,intent(in) :: pole
  real   ,optional,intent(in) :: reflon ! lon (in rad) at i=1
  integer,optional,intent(in) :: skipuv ! -1 if u(:) is to be ignored
                                        ! +1 if v(:) is to be ignored
  integer :: i,nlon
  real :: slat,dlon,rlon,slon,clon
  real :: elx,ely,emx,emy
  real :: emxx,emxy,emyy, emxd,emyd
  real :: reflon_
  integer :: skipuv_
  reflon_=0.
  if(present(reflon)) reflon_=reflon
  skipuv_=0
  if(present(skipuv)) skipuv_=skipuv

slat=pole     ! slat == -1 for SP (rlat=-pi/2.) and
              ! slat == +1 for NP (rlat=+pi/2.)
nlon=size(v)
dlon=4.*asin(1.)/nlon
emxx=0.
emxy=0.
emyy=0.
emxd=0.
emyd=0.
do i=1,nlon
  rlon=(i-1.+stag)*dlon + reflon_
  slon=sin(rlon)
  clon=cos(rlon)

  if(skipuv_<=0) then  ! .not. (skipuv_>0) .eqv. .not.skipv
    emx=-clon*slat
    emy=-slon*slat
    emxx=emxx+emx*emx
    emxy=emxy+emx*emy
    emyy=emyy+emy*emy
    emxd=emxd+emx*v(i)
    emyd=emyd+emy*v(i)
  endif

  if(skipuv_>=0) then  ! .not.(skipuv_ <0) .eqv. .not.skipu
    elx=-slon
    ely=+clon
    emxx=emxx+elx*elx
    emxy=emxy+elx*ely
    emyy=emyy+ely*ely
    emxd=emxd+elx*u(i)
    emyd=emyd+ely*u(i)
  endif
end do
wy=1./(emxx*emyy-emxy*emxy)
wx=(+emyy*emxd-emxy*emyd)*wy
wy=(-emxy*emxd+emxx*emyd)*wy
end subroutine getxy_

subroutine getxy_ad_(wx,wy,u,v,stag,pole,reflon,skipuv)
  implicit none
  real,intent(in) :: wx,wy
  real,dimension(:),intent(inout) :: u,v
  real,intent(in) :: stag
  real,intent(in) :: pole
  real   ,optional,intent(in) :: reflon ! lon (in rad) at i=1
  integer,optional,intent(in) :: skipuv ! -1 if u(:) is to be ignored
                                        ! +1 if v(:) is to be ignored
  integer :: i,nlon
  real :: slat,dlon,rlon,slon,clon
  real :: elx,ely,emx,emy,ww
  real :: emxx,emxy,emyy, emxd,emyd
  real :: reflon_
  integer :: skipuv_
  reflon_=0.
  if(present(reflon)) reflon_=reflon
  skipuv_=0
  if(present(skipuv)) skipuv_=skipuv

slat=pole     ! slat == -1 for SP (rlat=-pi/2.) and
              ! slat == +1 for NP (rlat=+pi/2.)
nlon=size(v)
dlon=4.*asin(1.)/nlon
emxx=0.
emxy=0.
emyy=0.
emxd=0.
emyd=0.
do i=1,nlon
  rlon=(i-1.+stag)*dlon + reflon_
  slon=sin(rlon)
  clon=cos(rlon)

  if(skipuv_<=0) then  ! .not. (skipuv_>0) .eqv. .not.skipv
    emx=-clon*slat
    emy=-slon*slat
    emxx=emxx+emx*emx
    emxy=emxy+emx*emy
    emyy=emyy+emy*emy
  endif

  if(skipuv_>=0) then  ! .not.(skipuv_ <0) .eqv. .not.skipu
    elx=-slon
    ely=+clon
    emxx=emxx+elx*elx
    emxy=emxy+elx*ely
    emyy=emyy+ely*ely
  endif
end do
ww=1./(emxx*emyy-emxy*emxy)
emxd = (emyy*wx - emxy*wy)*ww
emyd = (emxx*wy - emxy*wx)*ww

do i=1,nlon
  rlon=(i-1.+stag)*dlon + reflon_
  slon=sin(rlon)
  clon=cos(rlon)

  if(skipuv_<=0) then  ! .not. (skipuv_>0) .eqv. .not.skipv
    emx=-clon*slat
    emy=-slon*slat
    v(i) = v(i) + emx*emxd + emy*emyd
  endif

  if(skipuv_>=0) then  ! .not.(skipuv_ <0) .eqv. .not.skipu
    elx=-slon
    ely=+clon
    u(i) = u(i) + elx*emxd + ely*emyd
  endif
end do
end subroutine getxy_ad_

subroutine setuv_(wx,wy,u,v,stag,pole,reflon,skipuv)
  implicit none
  real,intent(in) :: wx,wy ! x and y components of the polar wind
  real,dimension(:),intent(inout) :: u
  real,dimension(:),intent(inout) :: v
  real,intent(in) :: stag  ! in grid index
  real,intent(in) :: pole  ! in sin(lat)
  real,optional,intent(in) :: reflon ! in radiance
  integer,optional,intent(in) :: skipuv ! -1/=1 to skip u/v

  integer :: i,nlon
  real :: dlon,rlon,slon,clon,slat
  real :: elx,ely,emx,emy
  real :: reflon_
  integer :: skipuv_

  reflon_=0.
  if(present(reflon)) reflon_=reflon
  skipuv_=0
  if(present(skipuv)) skipuv_=skipuv

nlon=size(v)
dlon=4.*asin(1.)/nlon
slat=pole
do i=1,nlon
  rlon=(i-1.+stag)*dlon + reflon_
  slon=sin(rlon)
  clon=cos(rlon)

  if(.not. skipuv_>0) then
    emx=-clon*slat
    emy=-slon*slat
    v(i)=wx*emx+wy*emy
  endif

  if(.not. skipuv_<0) then
    elx=-slon
    ely=+clon
    u(i)=wx*elx+wy*ely
  endif
end do
end subroutine setuv_

subroutine setuv_ad_(wx,wy,u,v,stag,pole,reflon,skipuv)
  implicit none
  real,intent(inout) :: wx,wy ! x and y components of the polar wind
  real,dimension(:) :: u
  real,dimension(:) :: v
  real,intent(in) :: stag  ! in grid index
  real,intent(in) :: pole  ! in sin(lat)
  real,optional,intent(in) :: reflon ! in radiance
  integer,optional,intent(in) :: skipuv ! -1/=1 to skip u/v

  integer :: i,nlon
  real :: dlon,rlon,slon,clon,slat
  real :: elx,ely,emx,emy
  real :: reflon_
  integer :: skipuv_

  reflon_=0.
  if(present(reflon)) reflon_=reflon
  skipuv_=0
  if(present(skipuv)) skipuv_=skipuv

nlon=size(v)
dlon=4.*asin(1.)/nlon
slat=pole
do i=1,nlon
  rlon=(i-1.+stag)*dlon + reflon_
  slon=sin(rlon)
  clon=cos(rlon)

  if(.not. skipuv_>0) then
    emx=-clon*slat
    emy=-slon*slat
    wx = wx + v(i)*emx
    wy = wy + v(i)*emy
  endif

  if(.not. skipuv_<0) then
    elx=-slon
    ely=+clon
    wx = wx + u(i)*elx
    wy = wy + u(i)*ely
  endif
end do
end subroutine setuv_ad_

end module m_daInterp
