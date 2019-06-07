#define _SPECSOLV_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_gsiGriddedState - Convert GEOS FV to GSI state
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_gsiGriddedState
      implicit none
      private	! except

      public :: gsiGriddedState_convert		! v = f(w)
      public :: gsiGriddedState_invconv		! w = h(v)
      public :: gsiGriddedState_increto		! w += h(v; w), or
      						! w += H(v-f(w); w)
      public :: fvGvec

    interface gsiGriddedState_convert; module procedure &
      fwdconvert_; end interface
    interface gsiGriddedState_invconv; module procedure &
      invconvert_; end interface
    interface gsiGriddedState_increto; module procedure &
      incrUpdate_; end interface

! !REVISION HISTORY:
! 	27Sep06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       22Jun09 - Todling - define default fv-vector as GEOS-5
!       22Jun09 - Jing Guo
!		- fixed a couple of missing deallocate().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_gsiGriddedState'
  real,parameter :: kPa_per_Pa = 1./1000.
  real,parameter :: Pa_per_kPa = 1000./1.
  integer,save:: fvGvec = 5 ! defined to be new GEOS-5 grid

#include "assert.H"
#include "mytrace.H"
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fwdconvert_ - forward convertion to gsiGriddedState
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine fwdconvert_(im,rlats,ak,bk, &
      w_phis,w_psfc,w_pt,w_qw,w_uw,w_vw,   &
      zsfc,lnps,tv,qw,di,ze, uw,vw, mwave )

      use m_daInterp,only : daInterp_vdtoa ! assume lat-lon
      use m_ggGradientSP,only : ggGradientSP
      use m_ggGradientSP,only : ggGradientSP_init
      use m_ggGradientSP,only : ggDIVo
      use m_ggGradientSP,only : clean
      use m_geosapi,only : getcon
      use m_die,only : die,assert_
      use m_mpout, only: mpout_log

      implicit none
      					! grid information
      integer,intent(in) :: im
      real,dimension(:),intent(in) :: rlats
      real,dimension(:),intent(in) :: ak,bk

      					! FV state
      real,dimension(:,:),intent(in) :: w_phis ! surface geo-p. height
      real,dimension(:,:),intent(in) :: w_psfc ! surface pressure
      real,dimension(:,:,:),intent(in) :: w_pt ! potential temperatur
      real,dimension(:,:,:),intent(in) :: w_qw ! water vapor mix. ratio
      real,dimension(:,:,:),intent(in) :: w_uw ! d-grid u-wind
      real,dimension(:,:,:),intent(in) :: w_vw ! d-grid v-wind

      					! GSI state on the same grid
      real,dimension(:,:),intent(out) :: zsfc ! elevation
      real,dimension(:,:),intent(out) :: lnps ! log(p_sfc:kPa)
      real,dimension(:,:,:),intent(out) :: tv ! virtual temperature
      real,dimension(:,:,:),intent(out) :: qw ! water vapor mixing ratio
      real,dimension(:,:,:),intent(out) :: di ! divergence
      real,dimension(:,:,:),intent(out) :: ze ! vorticity

      real,target,dimension(:,:,:),optional,intent(out) :: uw ! u-wind
      real,target,dimension(:,:,:),optional,intent(out) :: vw ! v-wind

      integer,optional,intent(in) :: mwave ! number of wave defined only
      					   ! for wind spectral solver

! !REVISION HISTORY:
! 	27Sep06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       22Jun09 - Todling - handle for GEOS-5 vector
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fwdconvert_'
  integer :: jm,km
  real,allocatable,dimension(:,:) :: sm,rm

#ifndef _SPECSOLV_
  type(ggGradientSP) :: gg
#else
  integer :: MWAVE_,MAXWV,MAXWV2
  real,allocatable,dimension(:,  :) :: s_di,s_ze
#endif
  real,pointer    ,dimension(:,:,:) :: a_uw,a_vw

  integer :: i,j,k,kp
  real :: p,s,r
  real :: KAPPA
  real :: GRAV
_ENTRY_
  	! Verify the consistency of some arguments
  jm=size(rlats)
  km=size(ak)-1

    ASSERT(km==size(bk)-1)

    ASSERT(im==size(w_phis,1))
    ASSERT(jm==size(w_phis,2))

    ASSERT(im==size(w_psfc,1))
    ASSERT(jm==size(w_psfc,2))

    ASSERT(im==size(w_pt,1))
    ASSERT(jm==size(w_pt,2))
    ASSERT(km==size(w_pt,3))

    ASSERT(im==size(w_qw,1))
    ASSERT(jm==size(w_qw,2))
    ASSERT(km==size(w_qw,3))

    ASSERT(im==size(w_uw,1))
    ASSERT(jm==size(w_uw,2))
    ASSERT(km==size(w_uw,3))

    ASSERT(im==size(w_vw,1))
    ASSERT(jm==size(w_vw,2))
    ASSERT(km==size(w_vw,3))

    ASSERT(im==size(zsfc,1))
    ASSERT(jm==size(zsfc,2))

    ASSERT(im==size(lnps,1))
    ASSERT(jm==size(lnps,2))

    ASSERT(im==size(tv,1))
    ASSERT(jm==size(tv,2))
    ASSERT(km==size(tv,3))

    ASSERT(im==size(qw,1))
    ASSERT(jm==size(qw,2))
    ASSERT(km==size(qw,3))

    ASSERT(im==size(di,1))
    ASSERT(jm==size(di,2))
    ASSERT(km==size(di,3))

    ASSERT(im==size(ze,1))
    ASSERT(jm==size(ze,2))
    ASSERT(km==size(ze,3))

  if(present(uw)) then
    ASSERT(im==size(uw,1))
    ASSERT(jm==size(uw,2))
    ASSERT(km==size(uw,3))
  endif

  if(present(vw)) then
    ASSERT(im==size(vw,1))
    ASSERT(jm==size(vw,2))
    ASSERT(km==size(vw,3))
  endif

    GRAV=getcon('GRAVITY')
  zsfc(:,:) = w_phis(:,:)/GRAV
  lnps(:,:) = log(w_psfc(:,:)*kPa_per_Pa)
  
 if(fvGvec==4)then
    KAPPA=getcon('KAPPA')
    allocate(sm(im,jm),rm(im,jm))
  do j=1,jm; do i=1,im
    p=ak(1)+bk(1)*w_psfc(i,j)
    sm(i,j)=p**KAPPA
    rm(i,j)=log(sm(i,j))
  end do; end do
  do k=1,km
    kp=k+1
    do j=1,jm; do i=1,im
      p=ak(kp)+bk(kp)*w_psfc(i,j)
      s=p**KAPPA
      r=log(s)
      tv(i,j,k)=w_pt(i,j,k)*(s-sm(i,j))/(r-rm(i,j))
      sm(i,j)=s
      rm(i,j)=r
    end do; end do
  end do
    deallocate(sm,rm)
 else  ! fvGvec
    tv=w_pt ! already tv
 endif ! fvGvec

  	! qw: w.m.r.
  qw(:,:,:) = w_qw(:,:,:)

  	! uw,vw,di,ze
    if(present(uw)) then
      a_uw => uw
    else
      allocate(a_uw(im,jm,km))
    endif
    if(present(vw)) then
      a_vw => vw
    else
      allocate(a_vw(im,jm,km))
    endif

  if (fvGvec==4) then 
    call daInterp_vdtoa(w_uw,w_vw, a_uw,a_vw)
  else
    a_uw=w_uw; a_vw=w_vw
  endif

#ifndef _SPECSOLV_
    call ggGradientSP_init(gg,im,rlats)
    call ggDivo(gg,a_uw,a_vw,di,ze)
    call clean(gg)

#else

!!  nlat=(jcap+1)*3/2
!!  nlat=nlat+2           ! to include the polse
!!  nlon=nlat*2
!!  jcap ?= (nlat-2)*2/3-1

!!  From nlat and nlon as functions of jcap shown above, an expression
!!  of MWAVE_ (:=jcap) is derived below.

  if(mod(im,2)/=0) call die(myname_,'an even im is required',im)

  MWAVE_ = ((jm-3+1)/2*2+2)/3  ! MWAVE for IDRT=0 and IROMB=0
  MWAVE_ = (jm+1)*2/3-1  ! however, this result looks better
  MWAVE_ = (jm-1)*2/3    ! or this?
  if(present(mwave)) MWAVE_=mwave
  	print'(2a,3i5)',myname_,': IM,JM = ',im,jm
  	print'(2a,1i5)',myname_,': MWAVE = ',MWAVE_
  MAXWV = (MWAVE_+1)*(MWAVE_+2)/2
  MAXWV2= MAXWV*2

	allocate(s_di(MAXWV2,km),s_ze(MAXWV2,km))

    ! Derive spectral (s_di,s_e) from a-grid (a_uw,a_vw).  The
    ! latitude grid of (uw,vw) is swapped before used by SPTEZMV() as
    ! the input.
  call SPTEZMV(0,MWAVE_,0,im,jm,km,s_di,s_ze, &
    a_uw(:,jm:1:-1,:),a_vw(:,jm:1:-1,:), -1)

    ! Derive a-grid (di,ze) from spectral (s_di,s_ze).
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_di,di,+1)
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_ze,ze,+1)

    di(:,1:jm,:)=di(:,jm:1:-1,:)
    ze(:,1:jm,:)=ze(:,jm:1:-1,:)

  	deallocate(s_di,s_ze)
#endif

  if(.not.present(uw)) deallocate(a_uw); nullify(a_uw)
  if(.not.present(vw)) deallocate(a_vw); nullify(a_vw)
_EXIT_
end subroutine fwdconvert_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: invconvert_ - update through inverse conversion
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine invconvert_(im,rlats,ak,bk, &
        w_psfc,w_dp,w_pt,w_qw,w_uw,w_vw, lnps,tv,qw,di,ze, mwave)

      use m_daInterp,only : daInterp_vatod ! assume lat-lon
      use m_ggGradientSP,only : ggGradientSP
      use m_ggGradientSP,only : ggGradientSP_init
      use m_ggGradientSP,only : ggPsolv
      use m_ggGradientSP,only : ggGrad,ggStrm
      use m_ggGradientSP,only : clean
      use m_geosapi,only : getcon
      use m_die,only : die,assert_
      use m_mpout, only: mpout_log

      implicit none
      				! Grid imformation
      integer,intent(in) :: im
      real,dimension(:),intent(in) :: rlats
      real,dimension(:),intent(in) :: ak,bk

				! wf is a state to be Updated.
      real,dimension(:,:  ),intent(inout) :: w_psfc
      real,dimension(:,:,:),intent(inout) :: w_dp
      real,dimension(:,:,:),intent(inout) :: w_pt
      real,dimension(:,:,:),intent(inout) :: w_qw
      real,dimension(:,:,:),intent(inout) :: w_uw
      real,dimension(:,:,:),intent(inout) :: w_vw

				! Variables below is the state set as
				! the target of updating wf.  The
				! physical differences from wf are
				! expected to be small, such that
				! increasing tangent-linearly is still
				! valid.
      real,dimension(:,:  ),intent(in) :: lnps
      real,dimension(:,:,:),intent(in) :: tv
      real,dimension(:,:,:),intent(in) :: qw
      real,dimension(:,:,:),intent(in) :: di ! 2-d divergence
      real,dimension(:,:,:),intent(in) :: ze ! 2-d vorticity

      integer,optional,intent(in) :: mwave ! number of wave defined only
      					   ! for wind spectral solver

! !REVISION HISTORY:
! 	27Sep06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       22Jun09 - Todling - handle for GEOS-5 vector
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::invconvert_'
  integer :: jm,km
  real,allocatable,dimension(:,:) :: sm,rm
  real,allocatable,dimension(:,:,:) :: a_uw,a_vw,u_uw
  real,allocatable,dimension(:,:,:) :: chi,psi

  integer :: i,j,k,kp,ier
  real :: p,s,r
  real :: KAPPA
#ifndef _SPECSOLV_
  type(ggGradientSP) :: gg
#else
  integer :: MWAVE_,MAXWV,MAXWV2
  real,allocatable,dimension(:,  :) :: s_di,s_ze
#endif

_ENTRY_
  	! Verify the consistency of some arguments
  jm=size(rlats)
  km=size(ak)-1

    ASSERT(km==size(bk)-1)

    ASSERT(im==size(w_psfc,1))
    ASSERT(jm==size(w_psfc,2))

    ASSERT(im==size(w_dp,1))
    ASSERT(jm==size(w_dp,2))
    ASSERT(km==size(w_dp,3))

    ASSERT(im==size(w_pt,1))
    ASSERT(jm==size(w_pt,2))
    ASSERT(km==size(w_pt,3))

    ASSERT(im==size(w_qw,1))
    ASSERT(jm==size(w_qw,2))
    ASSERT(km==size(w_qw,3))

    ASSERT(im==size(w_uw,1))
    ASSERT(jm==size(w_uw,2))
    ASSERT(km==size(w_uw,3))

    ASSERT(im==size(w_vw,1))
    ASSERT(jm==size(w_vw,2))
    ASSERT(km==size(w_vw,3))

    ASSERT(im==size(lnps,1))
    ASSERT(jm==size(lnps,2))

    ASSERT(im==size(tv,1))
    ASSERT(jm==size(tv,2))
    ASSERT(km==size(tv,3))

    ASSERT(im==size(qw,1))
    ASSERT(jm==size(qw,2))
    ASSERT(km==size(qw,3))

    ASSERT(im==size(di,1))
    ASSERT(jm==size(di,2))
    ASSERT(km==size(di,3))

    ASSERT(im==size(ze,1))
    ASSERT(jm==size(ze,2))
    ASSERT(km==size(ze,3))

	! Compute ps := Pa_per_kPa*exp(lnps)
  w_psfc(:,:) = Pa_per_kPa*exp(lnps(:,:))

  	! Determine the overall sign of p(:,:,1:km)-p(:,:,2:km+1).
	! s > 0, i.e. p(k+1)>p(k) is expected.  i.e. 
  s = 1.
  p = minval(w_psfc(:,:))
  kp=km+1
  ALWAYS_ASSERT(ak(1)+bk(1)*p < ak(kp)+bk(kp)*p)

  	! Compute delp(k) := sign(p(k+1)-p(k),1.)
  do k=1,km
    kp=k+1

      !<< delp(k) := s*(p(kp)-p(k))	( >=0. )
      !<< p(k) := a(k)+b(k)*ps
      !>> delp(k) = s*((a(kp)-a(k)) + (b(kp)-b(k))*ps)

    w_dp(:,:,k) = (ak(kp)-ak(k)) + (bk(kp)-bk(k))*w_psfc(:,:)
  end do
  ASSERT(all(w_dp>0.))

    !<< pt(k) = tv(k)*(log(pke(kp))-log(pke(k)))/(pke(kp)-pke(k))
    !         = tv(k)*(r(kp)-r(k))/(s(kp)-s(k))
    !         = tv(k)*(R/S); for k=1,...,km and kp:=k+1
    !:: R := r(kp)-r(k)
    !:: S := s(kp)-s(k)
    !:: s(k) := pke(k)
    !:: r(k) := log(pke(k))=log(s(k))
    !:: pke(k) := p(k)**kap
    !:: p(k) = a(k)+b(k)*ps

      ! a bunch of workspace arrays to store variables computed at the
      ! current layer for the next layer.

 if (fvGvec==4) then
    KAPPA=getcon('KAPPA')
    allocate(sm(im,jm),rm(im,jm))

  do j=1,jm; do i=1,im
    p = ak(1)+bk(1)*w_psfc(i,j)
    sm(i,j) = p**KAPPA
    rm(i,j) = log(sm(i,j))
  end do; end do

  do k=1,km
    kp=k+1
    do j=1,jm; do i=1,im
      p = ak(kp)+bk(kp)*w_psfc(i,j)
      s = p**KAPPA
      r = log(s)

      w_pt(i,j,k) = tv(i,j,k)*(r-rm(i,j))/(s-sm(i,j))

		! saved for the next layer.
      rm(i,j)=r
      sm(i,j)=s
    end do; end do
  end do
    deallocate(sm,rm)

 else ! <fvGvec>
    w_pt = tv ! need tv and already tv
 endif ! < fvGvec>

  	! qw: water vapor mixing ratio, in gram/gram
  w_qw(:,:,:) = qw(:,:,:)

  	! uw,vw,di,ze

    allocate(a_uw(im,jm,km),a_vw(im,jm,km))

#ifndef _SPECSOLV_
   ! The limitation of requiring equal-spacing in rlats in ggPsolv()
   ! is verified inside ggPsolv() as an ASSERT().

    allocate( chi(im,jm,km), psi(im,jm,km))
    allocate(u_uw(im,jm,km))

  call ggGradientSP_init(gg,im,rlats)
  call ggPsolv(gg,di(:,:,:),chi(:,:,:),stat=ier)
  	if(ier/=0) call die(myname_,'ggPsolv(gg,di,chi)',ier)
  call ggPsolv(gg,ze(:,:,:),psi(:,:,:),stat=ier)
  	if(ier/=0) call die(myname_,'ggPsolv(gg,ze,psi)',ier)

  call ggGrad(gg,chi,a_uw,a_vw)
  call ggStrm(gg,psi,u_uw,w_vw)
    deallocate(chi,psi)

  a_uw=a_uw+u_uw
  a_vw=a_vw+w_vw
  call clean(gg)

  deallocate(u_uw)

#else
  if(mod(im,2)/=0) call die(myname_,'an even im is required',im)

  MWAVE_ = ((jm-3+1)/2*2+2)/3  ! MWAVE for IDRT=0 and IROMB=0
  MWAVE_ = (jm+1)*2/3-1  ! however, this result looks better
  MWAVE_ = (jm-1)*2/3    ! or this?
  if(present(mwave)) MWAVE_=mwave
  	print'(2a,3i5)',myname_,': IM,JM = ',im,jm
  	print'(2a,1i5)',myname_,': MWAVE = ',MWAVE_
  MAXWV = (MWAVE_+1)*(MWAVE_+2)/2
  MAXWV2= MAXWV*2

	allocate(s_di(MAXWV2,km),s_ze(MAXWV2,km))

    ! Derive spectral (s_di,s_ze) from a-grid (di,ze).
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_di,di(:,jm:1:-1,:),-1)
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_ze,ze(:,jm:1:-1,:),-1)

    ! Derive a-grid (a_uw,a_vw) from spectral (s_di,s_ze)
  call SPTEZMV(0,MWAVE_,0,im,jm,km,s_di,s_ze,a_uw,a_vw,+1)

  a_uw(:,1:jm,:)=a_uw(:,jm:1:-1,:)
  a_vw(:,1:jm,:)=a_vw(:,jm:1:-1,:)

  	deallocate(s_di,s_ze)
#endif

  if (fvGvec==4) then
  	! Flag uskip will leave slab w_uw(:,j=1,:) untouched.
    call daInterp_vatod(a_uw,a_vw,w_uw,w_vw,uskip=.true.)
  else  ! No need to convert (A-grid needed)
    w_uw=a_uw
    w_vw=a_vw
  endif ! <fvGvec>
      deallocate(a_uw,a_vw)
_EXIT_
end subroutine invconvert_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: incrUpdate_ - incrementally update through TL conversion
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine incrUpdate_(im,rlats,ak,bk, &
        w_psfc,w_dp,w_pt,w_qw,w_uw,w_vw, lnps,tv,qw,di,ze, mwave)

      use m_daInterp,only : daInterp_vatod ! assume lat-lon
      use m_ggGradientSP,only : ggGradientSP
      use m_ggGradientSP,only : ggGradientSP_init
      use m_ggGradientSP,only : ggPsolv
      use m_ggGradientSP,only : ggGrad,ggStrm
      use m_ggGradientSP,only : clean
      use m_geosapi,only : getcon
      use m_die,only : die,assert_
      use m_mpout, only: mpout_log

      implicit none
      				! Grid imformation
      integer,intent(in) :: im
      real,dimension(:),intent(in) :: rlats
      real,dimension(:),intent(in) :: ak,bk

				! wf is a state to be updated.
      real,dimension(:,:  ),intent(inout) :: w_psfc
      real,dimension(:,:,:),intent(inout) :: w_dp
      real,dimension(:,:,:),intent(inout) :: w_pt
      real,dimension(:,:,:),intent(inout) :: w_qw
      real,dimension(:,:,:),intent(inout) :: w_uw
      real,dimension(:,:,:),intent(inout) :: w_vw

				! Variables below is the state set as
				! the target of updating wf.  The
				! physical differences from wf are
				! expected to be small, such that
				! increasing tangent-linearly is still
				! valid.
      real,dimension(:,:  ),intent(in) :: lnps
      real,dimension(:,:,:),intent(in) :: tv
      real,dimension(:,:,:),intent(in) :: qw
      real,dimension(:,:,:),intent(in) :: di ! 2-d divergence
      real,dimension(:,:,:),intent(in) :: ze ! 2-d vorticity

      integer,optional,intent(in) :: mwave ! number of wave defined only
      					   ! for wind spectral solver

! !REVISION HISTORY:
! 	27Sep06	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       22Jun09 - Todling - handle for GEOS-5 vector
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::incrUpdate_'
  integer :: jm,km
  real,allocatable,dimension(:,:) :: b_lnps
  real,allocatable,dimension(:,:,:) :: b_tv,b_qw,b_di,b_ze
  real,allocatable,dimension(:,:) :: d_ps
  real,allocatable,dimension(:,:) :: sm,dsm,rm,drm
  real,allocatable,dimension(:,:,:) :: a_uw,a_vw,d_uw,d_vw
  real,allocatable,dimension(:,:,:) :: chi,psi

  integer :: i,j,k,kp,ier
  real :: p,s,ds,sd,dsd
  real ::   r,dr,rd,drd
  real :: dtv,btv,dpt
  real :: KAPPA
#ifndef _SPECSOLV_
  type(ggGradientSP) :: gg
#else
  integer :: MWAVE_,MAXWV,MAXWV2
  real,allocatable,dimension(:,  :) :: s_di,s_ze
#endif
_ENTRY_
  	! Verify the consistency of some arguments
  jm=size(rlats)
  km=size(ak)-1

    ASSERT(km==size(bk)-1)

    ASSERT(im==size(w_psfc,1))
    ASSERT(jm==size(w_psfc,2))

    ASSERT(im==size(w_dp,1))
    ASSERT(jm==size(w_dp,2))
    ASSERT(km==size(w_dp,3))

    ASSERT(im==size(w_pt,1))
    ASSERT(jm==size(w_pt,2))
    ASSERT(km==size(w_pt,3))

    ASSERT(im==size(w_qw,1))
    ASSERT(jm==size(w_qw,2))
    ASSERT(km==size(w_qw,3))

    ASSERT(im==size(w_uw,1))
    ASSERT(jm==size(w_uw,2))
    ASSERT(km==size(w_uw,3))

    ASSERT(im==size(w_vw,1))
    ASSERT(jm==size(w_vw,2))
    ASSERT(km==size(w_vw,3))

    ASSERT(im==size(lnps,1))
    ASSERT(jm==size(lnps,2))

    ASSERT(im==size(tv,1))
    ASSERT(jm==size(tv,2))
    ASSERT(km==size(tv,3))

    ASSERT(im==size(qw,1))
    ASSERT(jm==size(qw,2))
    ASSERT(km==size(qw,3))

    ASSERT(im==size(di,1))
    ASSERT(jm==size(di,2))
    ASSERT(km==size(di,3))

    ASSERT(im==size(ze,1))
    ASSERT(jm==size(ze,2))
    ASSERT(km==size(ze,3))

    allocate(b_lnps(im,jm),b_tv(im,jm,km),b_qw(im,jm,km), &
             d_ps  (im,jm),b_di(im,jm,km),b_ze(im,jm,km)  )

  call fwdconvert_(im,rlats,ak,bk,(w_psfc),w_psfc,w_pt,w_qw,w_uw,w_vw, &
  	d_ps,b_lnps,b_tv,b_qw,b_di,b_ze)
	! (w_psfc) and d_ps are only referenced above as scratch arrays
	! for w_phis and b_zsfc, both are not used in this procedure.

      !<< ps   = Pa_per_kPa*exp(lnps)
      !>> d_ps = Pa_per_kPa*exp(lnps)*d_lnps = ps*d_lnps

  d_ps(:,:) = w_psfc(:,:)*(lnps(:,:)-b_lnps(:,:))
    deallocate(b_lnps)

  	! Determine the overall sign of p(:,:,2:km+1)-p(:,:,1:km).
  s = 1.
  p = minval(w_psfc(:,:))
  kp=km+1
  ALWAYS_ASSERT(ak(1)+bk(1)*p < ak(kp)+bk(kp)*p)

  do k=1,km
    kp=k+1

      !<< delp(k) := p(kp)-p(k)
      !<< p(k) := a(k)+b(k)*ps
      !>> delp(k) = (a(kp)-a(k)) + (b(kp)-b(k))*ps
      !>> d_delp(k) = (b(kp)-b(k))*d_ps
      !>> d_delp(k) = (b(kp)-b(k))*ps*d_lnps, << d_ps=ps*d_lnps
      !! d_delpk(:,:) = (bk(kp)-bk(k))*d_ps_(:,:)
      !! w_delp(:,:,k) = delp(:,:,k) + d_delpk(:,:)

    w_dp(:,:,k) = w_dp(:,:,k) + (bk(kp)-bk(k))*d_ps(:,:)
  end do
  ASSERT(all(w_dp>0.))

    !<< pt(k) = tv(k)*(log(pke(kp))-log(pke(k)))/(pke(kp)-pke(k))
    !         = tv(k)*(r(kp)-r(k))/(s(kp)-s(k))
    !         = tv(k)*(R/S)
    !:: R := r(kp)-r(k)
    !:: S := s(kp)-s(k)
    !:: s(k) := pke(k)
    !:: r(k) := log(pke(k))=log(s(k))
    !:: pke(k) := p(k)**kap
    !
    !>> d_pt(k) = d_tv(k)*(R/S) + tv(k)*d_[R/S]
    !
    !<< d_[R/S] = d_R/S - R*d_S/(S*S)
    !           = d_R/R*(R/S) - d_S/S*(R/S)
    !           = (d_R/R-d_S/S)*(R/S)
    !:: d_R = d_r(kp)-d_r(k)
    !:: d_S = d_s(kp)-d_s(k)
    !
    !>> d_pt(k) = (d_tv(k) + tv(k)*(d_R/R-d_S/S))*(R/S)
    !
    !<< d_pke(k) = kap*p(k)**(kap-1)*d_p(k)
    !            = kap*pke(k)*d_p(k)/p(k)
    !>> d_s(k) = kap*s(k)*d_p(k)/p(k)
    !<< d_r(k) = d_[log(s(k))]
    !          = d_s(k)/s(k)
    !          = kap*d_p(k)/p(k)
    !:: d_p(k)/p(k) = b(k)*d_ps/(a(k)+b(k)*ps)
    !               = b(k)*d_ps/p(k)
    !>> d_s(k) = d_r(k)*s(k)
    !
    ! Reversely, given a k,
    !   p(k) = a(k)+b(k)*ps
    !   d_p(k) = b(k)*d_ps
    !   s(k) = p(k)**kappa
    !   r(k) = log(s(k))
    !   d_r(k) = kap*d_p(k)/p(k)
    !   d_s(k) = d_r(k)*s(k)
    !
    ! Then, at kp=k+1
    !   R = r(kp)-r(k)
    !   S = s(kp)-s(k)
    !   d_R = d_r(kp)-d_r(k)
    !   d_S = d_s(kp)-d_s(k)
    !
    ! Finally, for the k-th layer between the k-th and kp=k+1 levels,
    ! it is solved,
    !>> d_pt(k) = (d_tv(k) + tv(k)*(d_R/R-d_S/S))*(R/S)
    !
    ! This approach seems use less calculations than other more forms.

      ! a bunch of workspace arrays to store variables computed at the
      ! current layer for the next layer.

 if (fvGvec==4) then
    KAPPA=getcon('KAPPA')
    allocate(sm(im,jm),dsm(im,jm),rm(im,jm),drm(im,jm))

  do j=1,jm; do i=1,im
    p = ak(1)+bk(1)*w_psfc(i,j)
    sm(i,j) = p**KAPPA
    rm(i,j) = log(sm(i,j))
    drm(i,j) = KAPPA*bk(1)*d_ps(i,j)/p
    dsm(i,j) = sm(i,j)*drm(i,j)
  end do; end do

  do k=1,km
    kp=k+1
    do j=1,jm; do i=1,im
      btv = b_tv(i,j,k)
      dtv = tv(i,j,k)-btv ! d_ps(i,j) was computed earlier

      p = ak(kp)+bk(kp)*w_psfc(i,j)
      s = p**KAPPA
      r = log(s)
      dr = KAPPA*bk(kp)*d_ps(i,j)/p
      ds = s*dr

      rd = r-rm(i,j)
      drd = dr-drm(i,j)
      sd = s-sm(i,j)
      dsd = ds-dsm(i,j)

      dpt = rd/sd*(dtv+btv*(drd/rd-dsd/sd))
      w_pt(i,j,k) = w_pt(i,j,k) + dpt

		! saved for the next layer.
      rm(i,j)=r
      sm(i,j)=s
      drm(i,j)=dr
      dsm(i,j)=ds
    end do; end do
  end do
    deallocate(b_tv)
    deallocate(sm,dsm,rm,drm)

 else  ! <fvGvec>
    deallocate(b_tv)
    w_pt = tv	! == "w_pt = w_pt + (tv - b_tv)"
 endif ! <fvGvec>

  	! qw: w.m.r.
  w_qw(:,:,:) = w_qw(:,:,:) + qw(:,:,:)-b_qw(:,:,:)
    deallocate(b_qw)

  	! uw,vw,di,ze
    allocate(a_uw(im,jm,km),a_vw(im,jm,km))
    allocate(d_uw(im,jm,km),d_vw(im,jm,km))

#ifndef _SPECSOLV_
    allocate( chi(im,jm,km), psi(im,jm,km))

   ! The limitation of requiring equal-spacing in rlats in ggPsolv()
   ! is verified inside ggPsolv() as an ASSERT().

  call ggGradientSP_init(gg,im,rlats)
  call ggPsolv(gg,di(:,:,:)-b_di(:,:,:),chi(:,:,:),stat=ier)
  	if(ier/=0) call die(myname_,'ggPsolv(gg,di,chi)',ier)
  call ggPsolv(gg,ze(:,:,:)-b_ze(:,:,:),psi(:,:,:),stat=ier)
  	if(ier/=0) call die(myname_,'ggPsolv(gg,ze,psi)',ier)
    deallocate(b_di,b_ze)

  call ggGrad(gg,chi,a_uw,a_vw)
  call ggStrm(gg,psi,d_uw,d_vw)
    deallocate(chi,psi)

  a_uw=a_uw+d_uw
  a_vw=a_vw+d_vw
  call clean(gg)

#else
  if(mod(im,2)/=0) call die(myname_,'an even im is required',im)

  MWAVE_ = ((jm-3+1)/2*2+2)/3  ! MWAVE for IDRT=0 and IROMB=0
  MWAVE_ = (jm+1)*2/3-1  ! however, this result looks better
  MWAVE_ = (jm-1)*2/3    ! or this?
  if(present(mwave)) MWAVE_=mwave
  	print'(2a,3i5)',myname_,': IM,JM = ',im,jm
  	print'(2a,1i5)',myname_,': MWAVE = ',MWAVE_
  MAXWV = (MWAVE_+1)*(MWAVE_+2)/2
  MAXWV2= MAXWV*2

	allocate(s_di(MAXWV2,km),s_ze(MAXWV2,km))

    ! Derive spectral (s_di,s_ze) from a-grid (di,ze).
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_di, &
    di(:,jm:1:-1,:)-b_di(:,jm:1:-1,:), -1)
  call SPTEZM(0,MWAVE_,0,im,jm,km,s_ze, &
    ze(:,jm:1:-1,:)-b_ze(:,jm:1:-1,:),-1)

    ! Derive a-grid (a_uw,a_vw) from spectral (s_di,s_ze)
  call SPTEZMV(0,MWAVE_,0,im,jm,km,s_di,s_ze,a_uw,a_vw,+1)
   
  a_uw(:,1:jm,:)=a_uw(:,jm:1:-1,:)
  a_vw(:,1:jm,:)=a_vw(:,jm:1:-1,:)

  	deallocate(s_di,s_ze)
#endif

 if (fvGvec==4) then

  call daInterp_vatod(a_uw,a_vw,d_uw,d_vw)
    deallocate(a_uw,a_vw)

  	! values of w_uw(:,j=1,:) are not touched.
  w_uw(:,2:,:) = w_uw(:,2:,:) + d_uw(:,2:,:)
  w_vw(:,1:,:) = w_vw(:,1:,:) + d_vw(:,1:,:)
    deallocate(d_uw,d_vw)

  else  ! <fvGvec>
    deallocate(d_uw,d_vw)
    w_uw = w_uw + a_uw
    w_vw = w_vw + a_vw
    deallocate(a_uw,a_vw)
  endif ! <fvGvec>

  w_psfc(:,:) = w_psfc(:,:)+d_ps(:,:)
    deallocate(d_ps)
_EXIT_
end subroutine incrUpdate_
end module m_gsiGriddedState
