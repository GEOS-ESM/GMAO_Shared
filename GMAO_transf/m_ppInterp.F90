!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ppInterp - A different modularized mappm() in gmap.F
!
! !DESCRIPTION:
!	ppInterp() interpolates grid in vertical.  An input grid or an
!   output grid can be either in a *top-first* (height descending) or in
!   a *sfc-first* (height ascending) order.  In particular, the fv-dyn
!   grid is *top-first*; the gsi-guess grid is *sfc-first*.  Procedures
!   in this module always check the order of a given grid first, before
!   derive *top-first* temporary variables from the given grid.  The
!   temporary result of interpolation procedure is then reordered as
!   specified by the given output grid if necessary.
!
! !INTERFACE:
!#include "regime.H"

    module m_ppInterp
      implicit none
      private	! except

      public :: ppInterp		! data structure

      public :: ppInterp_init
      public :: ppInterp_clean, clean

      public :: ppInterp_intp

      public :: THTA	! vertical interpolation on p**kappa
      public :: MASS	! vertical interpolation on p
      public :: WIND	! vertical vector interpolation on p

      		! Type of the temperature variable at the RHS of the
		! forward interpolation (or at the LHS of the reversed
		! interpolation.

      public :: VIRTUAL_TEMPERATURE
      public :: POTENTIAL_TEMPERATURE   ! virtual-potential-temperature

      type ppInterp
        private
			! Source grid
	real,pointer,dimension(:,:  ) :: ps_srcS
	real,pointer,dimension(    :) :: ak_src	! Same unit as ps_src
	real,pointer,dimension(    :) :: bk_src	!

			! Target grid
	real,pointer,dimension(:,:  ) :: ps_tgtS
	real,pointer,dimension(    :) :: ak_tgt	! Same unit as ps_tgt
	real,pointer,dimension(    :) :: bk_tgt
      end type ppInterp

    interface ppInterp_init; module procedure	&
      init_; end interface

    interface ppInterp_clean; module procedure	&
      clean_; end interface
    interface clean; module procedure clean_; end interface

    interface ppInterp_intp; module procedure	&
      intp_; end interface


! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ppInterp'

  integer,parameter :: THTA= 1
  integer,parameter :: MASS= 0
  integer,parameter :: WIND=-1

  integer,parameter :: VIRTUAL_TEMPERATURE   = 1
  integer,parameter :: POTENTIAL_TEMPERATURE = 2

  real,parameter :: FORCED_PTOP = 1.	! in Pascal

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,aGrid,					&
	srcLev1,ps_srcL,phis_srcL,ak_src,bk_src,srcLevs,th_srcL,&
	tgtLev1,ps_tgtL,phis_tgtL,ak_tgt,bk_tgt,tgtLevs, comm,	&
	Temp,delp)

      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : totalSize
      use m_interleavedObject,only : localSize
      use m_interleavedObject,only : get
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributorComm,only : undistribute
      use m_checksums,only : checksums_show

      use m_mpout,only : mpout_log,mpout_ison,mpout
      use m_die,only : assert_,MP_die,die

      implicit none
      type(ppInterp),intent(out) :: ob

      type(SubdomainDistributor),intent(in) :: aGrid

		! The grid specifications of the source field.

        type(interleavedObject),intent(in) :: srcLev1
        real,dimension(:,:,:),intent(in) ::   ps_srcL
        real,dimension(:,:,:),intent(in) :: phis_srcL
        real,dimension(:    ),intent(in) ::   ak_src
        real,dimension(:    ),intent(in) ::   bk_src
        type(interleavedObject),intent(in) :: srcLevs
        real,dimension(:,:,:),intent(in) ::   th_srcL

		! The grid specifications of the target field.

        type(interleavedObject),intent(in) :: tgtLev1
        real,dimension(:,:,:),intent(out)::   ps_tgtL
        real,dimension(:,:,:),intent(in) :: phis_tgtL
        real,dimension(:)    ,intent(in) ::   ak_tgt
        real,dimension(:)    ,intent(in) ::   bk_tgt
        type(interleavedObject),intent(in) :: tgtLevs

      integer,intent(in) :: comm

      integer,optional,intent(in) :: Temp	! type of th_srcL
      real,optional,dimension(:,:,:),intent(out) :: delp ! of tgt

! !REVISION HISTORY:
! 	10Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Init_'

  integer :: iaLlen,iaGlen,jaLlen,jaGlen
  integer :: nSrcLev1,nSrcLevs,mSrcLev1,mSrcLevs
  integer :: ntgtLev1,ntgtLevs,mtgtLev1,mtgtLevs
  integer :: ier
  integer :: k,m,kb,ke
  integer :: srcT

  		! _[src/tgt]S : src or tgt grid, in subdomains

  real,allocatable,dimension(:,:,:) ::   th_srcS
  real,allocatable,dimension(:,:  ) :: phis_srcS
  real,allocatable,dimension(:,:  ) :: phis_tgtS
  real,allocatable,dimension(:,:,:) ::   dp_tgtS
!________________________________________
_ALLENTRY_

  call get(aGrid,itotalSize=iaGlen,jtotalSize=jaGlen,	&
  		 ilocalSize=iaLlen,jlocalSize=jaLlen)

  nSrcLev1=totalSize(srcLev1)
  mSrcLev1=localSize(srcLev1)

	ASSERT(nSrcLev1==1)

  	ASSERT(iaGlen  ==size(  ps_srcL,1))
  	ASSERT(jaGlen  ==size(  ps_srcL,2))
  	ASSERT(mSrcLev1==size(  ps_srcL,3))

  	ASSERT(iaGlen  ==size(phis_srcL,1))
  	ASSERT(jaGlen  ==size(phis_srcL,2))
  	ASSERT(mSrcLev1==size(phis_srcL,3))

  nSrcLevs=totalSize(srcLevs)
  mSrcLevs=localSize(srcLevs)

  	ASSERT(nSrcLevs==size(ak_src)-1)
  	ASSERT(nSrcLevs==size(bk_src)-1)

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Surface pressures on interleaved source grid levels")
  call checksums_show(ps_srcL,'UNKW',"ps_srcL")
  call showakbk_(ak_src/100.,bk_src,myname_,'source grid')
  call showakbk_(ak_tgt/100.,bk_tgt,myname_,'target grid')
#endif

  	ASSERT(iaGlen  ==size(  th_srcL,1))
  	ASSERT(jaGlen  ==size(  th_srcL,2))
  	ASSERT(mSrcLevs==size(  th_srcL,3))

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Potential Temp. on interleaved source grid levels")
  call checksums_show(th_srcL,'THTA',"th_srcL")
#endif

  ntgtLev1=totalSize(tgtLev1)
  mtgtLev1=localSize(tgtLev1)

	ASSERT(ntgtLev1==1)

  	ASSERT(iaGlen  ==size(  ps_tgtL,1))
  	ASSERT(jaGlen  ==size(  ps_tgtL,2))
	ASSERT(mtgtLev1==size(  ps_tgtL,3))		! 0 or 1

  	ASSERT(iaGlen  ==size(phis_tgtL,1))
  	ASSERT(jaGlen  ==size(phis_tgtL,2))
	ASSERT(mtgtLev1==size(phis_tgtL,3))		! 0 or 1

  ntgtLevs=totalSize(tgtLevs)
  mtgtLevs=localSize(tgtLevs)

  	ASSERT(ntgtLevs==size(ak_tgt)-1)
  	ASSERT(ntgtLevs==size(bk_tgt)-1)
  if(present(delp)) then
  	ASSERT(iaGlen  ==size(delp,1))
  	ASSERT(jaGlen  ==size(delp,2))
	ASSERT(mtgtLevs==size(delp,3))
  endif
!________________________________________

  	allocate(ob%ps_srcS(iaLlen,jaLlen))
  	allocate(ob%ak_src(nSrcLevs+1))
  	allocate(ob%bk_src(nSrcLevs+1))

  	allocate(ob%ps_tgtS(iaLlen,jaLlen))
  	allocate(ob%ak_tgt(ntgtLevs+1))
  	allocate(ob%bk_tgt(ntgtLevs+1))

  	allocate( phis_srcS(iaLlen,jaLlen))
  	allocate(   th_srcS(iaLlen,jaLlen,nSrcLevs))
  	allocate( phis_tgtS(iaLlen,jaLlen))

  call distribute(srcLev1,aGrid,  ps_srcL,ob%ps_srcS,comm)
  call distribute(srcLev1,aGrid,phis_srcL, phis_srcS,comm)

  ob%ak_src(:)=ak_src	! ak_src and bk_src are already replicated
  ob%bk_src(:)=bk_src	! on all PEs.

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Surface pressures on source grid subdomains")
  call checksums_show(ob%ps_srcS,'UNKW',"ob%ps_srcS")
#endif

  call distribute(srcLevs,aGrid,  th_srcL,   th_srcS,comm)
  call distribute(tgtLev1,aGrid,phis_tgtL, phis_tgtS,comm)

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Potential Temp. on source grid subdomains")
  call checksums_show(phis_srcS,'UNKW',"vt:phis_srcS")
  call checksums_show(ob%ps_srcS,'UNKW',"vt:ps_srcS")
  call checksums_show(th_srcS,'THTA',"vt:th_srcS")
  call checksums_show(phis_tgtS,'UNKW',"vt:phis_tgtS")
#endif

  srcT=POTENTIAL_TEMPERATURE
  if(present(Temp)) srcT=Temp
  select case(srcT)
  case(VIRTUAL_TEMPERATURE)
    call psVTXtrp_(phis_srcS,ob%ps_srcS,ob%ak_src,ob%bk_src,th_srcS, &
  		   phis_tgtS,ob%ps_tgtS)

  case(POTENTIAL_TEMPERATURE)
    call psPTXtrp_(phis_srcS,ob%ps_srcS,ob%ak_src,ob%bk_src,th_srcS, &
  		   phis_tgtS,ob%ps_tgtS)

  case default
    call die(myname_,'undefined tag ("[Temp=]")',srcT)
  end select

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Surface pressures on target grid subdomains")
  call checksums_show(ob%ps_tgtS,'UNKW',"ob%ps_tgtS")
#endif

	deallocate(phis_srcS)
	deallocate(  th_srcS)
	deallocate(phis_tgtS)

  call undistribute(tgtLev1,aGrid,ob%ps_tgtS,ps_tgtL,comm)

#ifdef DEBUG_CHECKSUMS
_ALWAYS_ALLTRACE_("Surface pressure on target grid")
  call checksums_show(ps_tgtL,'UNKW',"ps_tgtL")
#endif

  ob%ak_tgt(:)=ak_tgt(1:ntgtLevs+1)
  ob%bk_tgt(:)=bk_tgt(1:ntgtLevs+1)

  if(present(delp)) then
  	allocate(dp_tgtS(iaLlen,jaLlen,ntgtLevs))

    do k=1,ntgtLevs
      dp_tgts(:,:,k)=(ak_tgt(k+1)-ak_tgt(k)) +	&
		     (bk_tgt(k+1)-bk_tgt(k))*ob%ps_tgtS(:,:)
    end do

	! If ak(1)+bk(1)*ps is the surface, then ak(2)+bk(2)*ps would
	! be smaller, then dp would be negtive and should be reset to
	! -dp.

    if(sfcFirst_(ak_tgt,bk_tgt,maxval(ob%ps_tgtS)))	&
	dp_tgts(:,:,:)=-dp_tgtS(:,:,:)

    call undistribute(tgtLevs,aGrid,dp_tgts,delp,comm)
	deallocate(dp_tgts)
  endif

_ALLEXIT_
end subroutine init_
subroutine showakbk_(ak,bk,where,tag)
  use m_mpout, only: mpout_ison	! control
  use m_mpout, only: mpout	! logical unit
  use m_mpout, only: mpout_log	! message
  use m_die,only : assert_
  implicit none
  real,dimension(:),intent(in) :: ak,bk
  character(len=*),intent(in) :: where,tag

  integer :: nsig,k,kb,ke

  nsig=size(ak)-1
  	ASSERT(nsig==size(bk)-1)

  call mpout_log(where,tag//':nsig =',nsig,flush=.true.)
  if(.not.mpout_ison()) return

  write(mpout,'(2a)') where,':: -- nsig+1 (a,1000*b) levels in hPa--'
  do kb=1,nsig+1,4
    ke=min(kb+3,nsig+1)
	! expecting in hPa ...
    write(mpout,'(f9.4,7(f10.4))') (ak(k),1000.*bk(k),k=kb,ke)
  end do
end subroutine showakbk_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: psPTXtrp_ - calculate p_sfc at given topography
!
! !DESCRIPTION:
!   This version of p_sfc interpolation uses v.p.t. as the temperature
!   variable.  The geopotential height integration is based on equation
!
!	d(phi) = - cp*pt* d(p^k)
!
! !INTERFACE:

    subroutine psPTXtrp_(phis,ps,ak,bk,thv,phisx,psx)
      use m_geosapi,only : getcon
      use m_die,only : assert_
      use m_checksums,only : checksums_show
      implicit none
      real,dimension(:,:  ),intent(in) :: ps	! surface pressure
      real,dimension(:,:  ),intent(in) :: phis	! surface geo.port.hght.
      real,dimension(    :),intent(in) :: ak,bk	! vertical grid ref.
      real,dimension(:,:,:),intent(in) :: thv	! virt.port.temp.

      real,dimension(:,:  ),intent(in) :: phisx	! target geo.port.hght.
      real,dimension(:,:  ),intent(out) :: psx	! output sfc.pres.

! !REVISION HISTORY:
! 	14Apr05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::psPTXtrp_'
  integer :: msfc,mtop,minc,kmm,ksfc,ier
  integer :: km,k,m
  integer :: ilen,i
  integer :: jlen,j
  real,allocatable,dimension(:,:) :: phx	! a workspace
  real :: cp,kappa,akap,pk_,ph_,psmax

  ilen=size(thv,1)
  jlen=size(thv,2)
  	ASSERT(ilen==size(phis ,1))
  	ASSERT(jlen==size(phis ,2))
  	ASSERT(ilen==size(ps   ,1))
  	ASSERT(jlen==size(ps   ,2))
  	ASSERT(ilen==size(phisx,1))
  	ASSERT(jlen==size(phisx,2))
  	ASSERT(ilen==size(psx  ,1))
  	ASSERT(jlen==size(psx  ,2))

  km=size(thv,3)
  	ASSERT(km==size(ak)-1)
  	ASSERT(km==size(bk)-1)

  psmax=maxval(ps)
  if(sfcFirst_(ak,bk,psmax)) then
    msfc= 1
    mtop=km
    minc=+1
     kmm= 1
    ksfc=msfc-minc+kmm
	ASSERT(ksfc==1)
  else
    msfc=km
    mtop= 1
    minc=-1
     kmm= 0
    ksfc=msfc-minc+kmm
	ASSERT(ksfc==km+1)
  endif
#ifdef DEBUG_CHECKSUMS
call checksums_show(phis,'UNKW','y:phis:0')
call checksums_show(phisx,'UNKW','y:phisx:0')
call checksums_show(phis-phisx,'UNKW','y:phis-phisx:0')
#endif

  	allocate(phx(ilen,jlen))

  cp   =getcon('CP')
  kappa=getcon('KAPPA')
  akap =1./kappa

  phx(:,:)=phis(:,:)
  psx(:,:)=ak(ksfc)+bk(ksfc)*ps(:,:)
  psx(:,:)= psx(:,:)**kappa
#ifdef DEBUG_CHECKSUMS
call checksums_show(ps,'UNKW','y:ps:0')
call checksums_show(psx**akap,'UNKW','y:ps:00')
#endif

  	! If anywhere the target topography phisx is below the source
	! topography phis, the target surface pressure psx is determined
	! by extrapolating, from the source surface pressure ps at the
	! source topography phis.

  if( any( phx(:,:) > phisx(:,:) ) ) then
    do j=1,jlen
      do i=1,ilen
        if(phx(i,j) > phisx(i,j)) then

		! integrate psx _downward_ to psx at phisx
	  psx(i,j)=psx(i,j)-(phisx(i,j)-phx(i,j))/(cp*thv(i,j,msfc))

		! integrate phx _downward_ to phisx for consistancy
	  phx(i,j) = phisx(i,j)
	endif
      end do
    end do
  endif
#ifdef DEBUG_CHECKSUMS
call checksums_show(psx**akap,'UNKW','y:ps:1')
#endif

  do m=msfc,mtop,minc
    if( all( phx(:,:) >= phisx(:,:) ) ) exit

    	! m is for thv(_middle); k is for [a,b](_edge)
    k=m+kmm		! i.e. k = m + (k-m); where k-m=0 or 1
#ifdef DEBUG_CHECKSUMS
call checksums_show(psx**akap,'UNKW','y:ps:2')
#endif

    do j=1,jlen
      do i=1,ilen
        if(phx(i,j) >= phisx(i,j)) continue	! skip this point

        pk_ = (ak(k)+bk(k)*ps(i,j))**kappa

		! d(phi) = - cp*pt*d(p^k), compute ph at k

        ph_ = phx(i,j) - cp*thv(i,j,m)*(pk_-psx(i,j))

	if(ph_ < phisx(i,j)) then
		! increase p^k by a full grid interval
	  psx(i,j)=pk_
	  phx(i,j)=ph_

	else
		! increase p^k by a fractioanl portion of the
		! grid interval in phi, until phi is as high as
		! the target topography phisx, based on equation
		!
		!	d(p^k) = - d(phi)/(cp*pt)

          psx(i,j)=psx(i,j)-(phisx(i,j)-phx(i,j))/(cp*thv(i,j,m))
	  phx(i,j)=phisx(i,j)
	endif

      end do	! i=1,ilen
    end do	! j=1,jlen
  end do	! k=msfc,mtop
#ifdef DEBUG_CHECKSUMS
call checksums_show(psx**akap,'UNKW','y:ps:9')
#endif

  	! convert the working variable from p^k to p

  psx(:,:) = psx(:,:)**akap
#ifdef DEBUG_CHECKSUMS
call checksums_show(psx,'UNKW','y:ps:99')
#endif

  	deallocate(phx)
end subroutine psPTXtrp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: psVTXtrp_ - calculate p_sfc at given topography
!
! !DESCRIPTION:
!	This version of p_sfc interpolation uses virtual temperature
!   field as the temperature variable.  Increment of pressure vs.
!   geopotential height is determined from equation
!
!	d(phi) = - cp*vt*d(log(p^k))
!
! !INTERFACE:

    subroutine psVTXtrp_(phis,ps,ak,bk,tv,phisx,psx)
      use m_geosapi,only : getcon
      use m_die,only : assert_
      use m_checksums,only : checksums_show
      implicit none
      real,dimension(:,:  ),intent(in) :: ps	! surface pressure
      real,dimension(:,:  ),intent(in) :: phis	! surface geo.port.hght.
      real,dimension(    :),intent(in) :: ak,bk	! vertical grid ref.
      real,dimension(:,:,:),intent(in) :: tv	! virt.temp.

      real,dimension(:,:  ),intent(in) :: phisx	! target geo.port.hght.
      real,dimension(:,:  ),intent(out) :: psx	! output sfc.pres.

! !REVISION HISTORY:
! 	14Apr05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::psVTXtrp_'
  integer :: msfc,mtop,minc,kmm,ksfc,ier
  integer :: km,k,m
  integer :: ilen,i
  integer :: jlen,j
  real,allocatable,dimension(:,:) :: phx	! a workspace
  real :: cp,kappa,akap,pk_,ph_,psm

  ilen=size(tv,1)
  jlen=size(tv,2)
  	ASSERT(ilen==size(phis ,1))
  	ASSERT(jlen==size(phis ,2))
  	ASSERT(ilen==size(ps   ,1))
  	ASSERT(jlen==size(ps   ,2))
  	ASSERT(ilen==size(phisx,1))
  	ASSERT(jlen==size(phisx,2))
  	ASSERT(ilen==size(psx  ,1))
  	ASSERT(jlen==size(psx  ,2))

  km=size(tv,3)
  	ASSERT(km==size(ak)-1)
  	ASSERT(km==size(bk)-1)

  psm=maxval(ps)
  if(sfcFirst_(ak,bk,psm)) then
    msfc= 1
    mtop=km
    minc=+1
     kmm= 1
    ksfc=msfc-minc+kmm
     ASSERT(ksfc==1)
  else
    msfc=km
    mtop= 1
    minc=-1
     kmm= 0
    ksfc=msfc-minc+kmm
     ASSERT(ksfc==km+1)
  endif
#ifdef DEBUG_CHECKSUMS
call checksums_show(phis,'UNKW','x:phis:0')
call checksums_show(phisx,'UNKW','x:phisx:0')
call checksums_show(phis-phisx,'UNKW','x:phis-phisx:0')
#endif

  	allocate(phx(ilen,jlen))

  cp   =getcon('CP')
  kappa=getcon('KAPPA')
  akap =1./kappa

  phx(:,:)=phis(:,:)
  psx(:,:)=ak(ksfc)+bk(ksfc)*ps(:,:)
  psx(:,:)=kappa*log(psx(:,:))
#ifdef DEBUG_CHECKSUMS
call checksums_show(ps,'UNKW','x:ps:0')
call checksums_show(exp(psx*akap),'UNKW','x:ps:00')
#endif

  	! If anywhere the target topography phisx is below the source
	! topography phis, the target surface pressure psx is determined
	! by extrapolating, from the source surface pressure ps at the
	! source topography phis.

  if( any( phx(:,:) > phisx(:,:) ) ) then
    do j=1,jlen
      do i=1,ilen
        if(phx(i,j) > phisx(i,j)) then

		! integrate psx _downward_ to psx at phisx
	  psx(i,j)=psx(i,j)-(phisx(i,j)-phx(i,j))/(cp*tv(i,j,msfc))

		! integrate phx _downward_ to phisx for consistancy
	  phx(i,j)=phisx(i,j)
	endif
      end do
    end do
  endif
#ifdef DEBUG_CHECKSUMS
call checksums_show(exp(psx*akap),'UNKW','x:ps:1')
#endif

  do m=msfc,mtop,minc
    if( all( phx(:,:) >= phisx(:,:) ) ) exit

    	! m is for tv(_middle); k is for [a,b](_edge)
    k=m+kmm		! i.e. k = m + (k-m); where k-m=0 or 1
#ifdef DEBUG_CHECKSUMS
call checksums_show(exp(psx*akap),'UNKW','x:ps:2')
#endif

    do j=1,jlen
      do i=1,ilen
        if(phx(i,j) >= phisx(i,j)) continue

        pk_ = kappa*log(ak(k)+bk(k)*ps(i,j))

		! d(phi) = - cp*vt*d(log(p^k)), compute ph at k

        ph_ = phx(i,j) - cp*tv(i,j,m)*(pk_-psx(i,j))

	if(ph_ < phisx(i,j)) then
		! increase log(p^k) by a full grid interval
	  psx(i,j)=pk_
	  phx(i,j)=ph_

	else
		! increase log(p^k) by a fractioanl portion of the
		! grid interval in phi, until phi is as high as
		! the target topography phisx, based on equation
		!
		!	d(log(p^k)) = - d(phi)/(cp*vt)

          psx(i,j)=psx(i,j)-(phisx(i,j)-phx(i,j))/(cp*tv(i,j,m))
	  phx(i,j)=phisx(i,j)
	endif

      end do	! i=1,ilen
    end do	! j=1,jlen
  end do	! k=msfc,mtop
#ifdef DEBUG_CHECKSUMS
call checksums_show(exp(psx*akap),'UNKW','x:ps:9')
#endif

  	! convert the working variable from log(p^k) to p

  psx(:,:) = exp(psx(:,:)*akap)
#ifdef DEBUG_CHECKSUMS
call checksums_show(psx,'UNKW','x:ps:99')
#endif

  	deallocate(phx)
end subroutine psVTXtrp_

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
      type(ppInterp),intent(inout) :: ob

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_

  deallocate(ob%ps_srcS)
  deallocate(ob%ak_src)
  deallocate(ob%bk_src)
  deallocate(ob%ps_tgtS)
  deallocate(ob%ak_tgt)
  deallocate(ob%bk_tgt)

_EXIT_
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_ - eta to sigma interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_(ob,fm,fn,key,SrcTemp,TgtTemp)
      use m_mpout,only : mpout_log
      use m_die,only : assert_,perr,die
      implicit none
      type(ppInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in)  :: fm
      real,dimension(:,:,:),intent(out) :: fn
      integer,intent(in) :: key
      integer,optional,intent(in) :: SrcTemp	! type of the source th
      integer,optional,intent(in) :: TgtTemp	! type of the target th

! !REVISION HISTORY:

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_'
  integer :: srcT_
  integer :: tgtT_
  integer :: isz,jsz,ksz
_ENTRY_

  isz=size(fm,1)
  jsz=size(fm,2)
  ksz=size(fm,3)

   	ASSERT(isz==size(ob%ps_srcS,1))
   	ASSERT(jsz==size(ob%ps_srcS,2))
   	ASSERT(ksz==size(ob%ak_src)-1)
   	ASSERT(ksz==size(ob%bk_src)-1)

  isz=size(fn,1)
  jsz=size(fn,2)
  ksz=size(fn,3)

   	ASSERT(isz==size(ob%ps_tgtS,1))
   	ASSERT(jsz==size(ob%ps_tgtS,2))
   	ASSERT(ksz==size(ob%ak_tgt)-1)
   	ASSERT(ksz==size(ob%bk_tgt)-1)
!________________________________________

  select case(key)
  case(THTA)

    srcT_=POTENTIAL_TEMPERATURE
    if(present(SrcTemp)) srcT_=SrcTemp
    if( srcT_/=POTENTIAL_TEMPERATURE .and.	&
        srcT_/=VIRTUAL_TEMPERATURE )		&
	call die(myname_,'invalid SrcTemp',srcT_)

    tgtT_=POTENTIAL_TEMPERATURE
    if(present(TgtTemp)) tgtT_=TgtTemp
    if( tgtT_/=POTENTIAL_TEMPERATURE .and.	&
        tgtT_/=VIRTUAL_TEMPERATURE )		&
	call die(myname_,'invalid TgtTemp',tgtT_)

    call ppIntp_( ob%ps_srcS,ob%ak_src,ob%bk_src,fm,	&
		  ob%ps_tgtS,ob%ak_tgt,ob%bk_tgt,fn,	&
  	key,SrcTemp=srcT_,TgtTemp=tgtT_)

  case default
    !if(present(SrcTemp)) call perr(myname_,'SrcTemp ignored (key)',key)
    !if(present(TgtTemp)) call perr(myname_,'TgtTemp ignored (key)',key)

  	! Other variables
    call ppIntp_( ob%ps_srcS,ob%ak_src,ob%bk_src,fm,	&
		  ob%ps_tgtS,ob%ak_tgt,ob%bk_tgt,fn,	&
  	key)
  end select
_EXIT_
end subroutine intp_

function sfcFirst_(ak,bk,ps)
  implicit none
  real,dimension(:),intent(in) :: ak,bk
  real,intent(in) :: ps
  logical :: sfcFirst_
  integer :: m
  m=size(ak)
	!   p_sfc=pres(1)   >  p_top=pres(m)
  sfcFirst_= ak(1)+bk(1)*ps > ak(m)+bk(m)*ps
end function sfcFirst_

function dscend_(ak,bk)
  implicit none
  real,dimension(:),intent(in) :: ak,bk
  logical :: dscend_
  integer :: m
  m=size(ak)	! == size(bk)
	!   p_top=pres(1)   <  p_sfc=pres(m)
  dscend_= ak(1)+1.e5*bk(1) < ak(m)+1.e5*bk(m)
end function dscend_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ppIntp_ - pressure-to-pressure grid "mapping"
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine ppIntp_( ps_m,ak_m,bk_m,fm,		&
			ps_n,ak_n,bk_n,fn,	key,	&
			SrcTemp,TgtTemp)
      use m_geosapi,only : getcon
      use m_mpout,only : mpout_log
      use m_die,only : die,assert_
      implicit none
      real,dimension(:,:  ),intent(in)  :: ps_m
      real,dimension(    :),intent(in)  :: ak_m	! starting from top
      real,dimension(    :),intent(in)  :: bk_m	! starting from top
      real,dimension(:,:,:),intent(in)  :: fm	! top to sfc.
      real,dimension(:,:  ),intent(in)  :: ps_n
      real,dimension(    :),intent(in)  :: ak_n	! starting from sfc.
      real,dimension(    :),intent(in)  :: bk_n	! starting from sfc.
      real,dimension(:,:,:),intent(out) :: fn	! sfc. to top
      integer,intent(in) :: key
      integer,optional,intent(in) :: SrcTemp
      integer,optional,intent(in) :: TgtTemp

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ppIntp_'
  integer :: srcT_,tgtT_
  integer :: isz,jsz,km,kn,i,j,k,m
  real,allocatable,dimension(:,:) :: prm_,dpm_,fm_
  real,allocatable,dimension(:,:) :: prn_,dpn_,fn_
  logical :: sfcFirst_m,sfcFirst_n
  real :: psmax

  real :: kappa
_ENTRY_

  isz=size(fm,1)
  jsz=size(fm,2)
  km =size(fm,3)
  kn =size(fn,3)

  	ASSERT(isz==size(ps_m,1))
  	ASSERT(isz==size(ps_n,1))
  	ASSERT(isz==size(  fn,1))

  	ASSERT(jsz==size(ps_m,2))
  	ASSERT(jsz==size(ps_n,2))
  	ASSERT(jsz==size(  fn,2))

  	ASSERT( km==size(ak_m)-1)
  	ASSERT( km==size(bk_m)-1)

  	ASSERT( kn==size(ak_n)-1)
  	ASSERT( kn==size(bk_n)-1)
!________________________________________
  	allocate(prm_(isz,km+1),dpm_(isz,km),fm_(isz,km))
  	allocate(prn_(isz,kn+1),dpn_(isz,kn),fn_(isz,kn))

  srcT_=POTENTIAL_TEMPERATURE
  if(present(SrcTemp)) srcT_=SrcTemp
 	select case(srcT_)
	case(POTENTIAL_TEMPERATURE,VIRTUAL_TEMPERATURE)
	case default
	  call die(myname_,'invalid SrcTemp',srcT_)
	end select

  tgtT_=POTENTIAL_TEMPERATURE
  if(present(TgtTemp)) tgtT_=TgtTemp
 	select case(tgtT_)
	case(POTENTIAL_TEMPERATURE,VIRTUAL_TEMPERATURE)
	case default
	  call die(myname_,'invalid TgtTemp',tgtT_)
	end select

  psmax=maxval(ps_m)
  sfcFirst_m=sfcFirst_(ak_m,bk_m,psmax)

  psmax=maxval(ps_n)
  sfcFirst_n=sfcFirst_(ak_n,bk_n,psmax)

  do j=1,jsz
    select case(key)
    case(THTA)
      kappa = getcon('KAPPA')

      call setprdp_(prm_,dpm_,ak_m,bk_m,ps_m(:,j),	&
	myname_,j,swapab=sfcFirst_m,kappa=kappa)

      call setprdp_(prn_,dpn_,ak_n,bk_n,ps_n(:,j),	&
	myname_,j,swapab=sfcFirst_n,kappa=kappa)

    case(MASS,WIND)

      call setprdp_(prm_,dpm_,ak_m,bk_m,ps_m(:,j),	&
	myname_,j,swapab=sfcFirst_m)

      call setprdp_(prn_,dpn_,ak_n,bk_n,ps_n(:,j),	&
	myname_,j,swapab=sfcFirst_n)

    case default

      call die(myname_,'invalid key',key)
    end select
!________________________________________
! Get slab -j from fm(:,:,:), swap -k if necessary.

    call copyslab_(fm(:,j,:),fm_(:,:),swapk=sfcFirst_m)
!________________________________________
! Do a grid "mapping"
    select case(key)
    case(THTA)		! theta_v from/to t_v
	!________________________________________
	! convert the input, if required, before mapping takes place.
      select case(srcT_)
      case(VIRTUAL_TEMPERATURE)	! then convert it to p.t.

      		!  T_v(k)*[log(pke(k+1)-log(pke(k))] ==	&
		! TH_v(k)*[    pke(k+1)-    pke(k) ]
                                            
	do k=1,km
          fm_(:,k)=fm_(:,k) *				&
	  	(log(prm_(:,k+1))-log(prm_(:,k))) /	&
	  	(    prm_(:,k+1) -    prm_(:,k) )
	end do

      end select

	!________________________________________
	! map a slab of p.t. from one grid to another

      call mappm( km,prm_,dpm_,fm_(:,:),	&
		  kn,prn_,dpn_,fn_(:,:), isz,key,7)

	!________________________________________
	! convert the output, if required, after mapping takes place.
      select case(tgtT_)
      case(VIRTUAL_TEMPERATURE)	! then convert it from p.t.

      		!  T_v(k)*[log(pke(k+1)-log(pke(k))] ==	&
		! TH_v(k)*[    pke(k+1)-    pke(k) ]
                                            
        do k=1,kn
          fn_(:,k)=fn_(:,k) *				&
	  	(    prn_(:,k+1) -    prn_(:,k) ) /	&
		(log(prn_(:,k+1))-log(prn_(:,k)))
	end do
      end select

    case default	! other variables

      call mappm(km,prm_,dpm_,fm_(:,:),	&
      		 kn,prn_,dpn_,fn_(:,:), isz,key,7)

    end select
!________________________________________
! Set a slab -j to fn, swap -k if necessary.

    call copyslab_(fn_(:,:),fn(:,j,:),swapk=sfcFirst_n)
  end do	! j

  	deallocate(prm_,dpm_,fm_)
  	deallocate(prn_,dpn_,fn_)
_EXIT_
end subroutine ppIntp_

subroutine copyslab_(fi,fo,swapk)
  implicit none
  real,dimension(:,:),intent(in ) :: fi
  real,dimension(:,:),intent(out) :: fo
  logical,intent(in) :: swapk
  integer :: m
  if(swapk) then
    m=size(fi,2)
    fo(:,1:m:+1)=fi(:,m:1:-1)
  else
    fo(:,:)=fi(:,:)
  endif
end subroutine copyslab_

subroutine setprdp_(pr,dp,ak,bk,ps,myname_,j,swapab,kappa)
  use m_die,only : assert_
  implicit none
  real,dimension(:,:),intent(out) :: pr,dp
  real,dimension(  :),intent(in ) :: ak,bk
  real,dimension(:  ),intent(in ) :: ps
  character(len=*),intent(in) :: myname_
  integer,intent(in) :: j
  logical,intent(in) :: swapab
  real,optional,intent(in) :: kappa

  integer :: i,k,m,km
  km=size(ak)	! == size(bk) == size(pr,3), == size(dp,3)+1
  km=km-1

  if(swapab) then
    m=1
    do k=km+1,1,-1
      pr(:,k) = ak(m) + bk(m)*ps(:)
      m=m+1
    enddo

  else
    do k=km+1,1,-1
      pr(:,k) = ak(k) + bk(k)*ps(:)
    enddo
  endif

! Force a minimum ptop > 0., when necessary
  do i=1,size(pr,1)	! == size(ps,1)
    if(pr(i,1)<=FORCED_PTOP) pr(i,1)=FORCED_PTOP
  end do

#ifdef DEBUG_COLUMNCHECKING
  call show_minval_(j,ps(:),pr(:,:),myname_,"pr")
#endif

  ASSERT(all(ps(:)>0.))
  ASSERT(all(pr(:,2)>FORCED_PTOP))

  if(present(kappa)) pr(:,1:km+1) = pr(:,1:km+1)**kappa

  dp(:,1:km) = pr(:,2:km+1) - pr(:,1:km)
  ASSERT(all(dp(:,:)>0.))

end subroutine setprdp_

subroutine show_minval_(j,ps,pr,where,tag)
use m_mpout,only : mpout,mpout_ison
implicit none
integer,intent(in) :: j
real,dimension(:),intent(in) :: ps
real,dimension(:,:),intent(in) :: pr
character(len=*),intent(in) :: where
character(len=*),intent(in) :: tag

integer :: m,k,kk,mdim(2),km

if(.not.mpout_ison()) return

if(minval(pr(:,:))<FORCED_PTOP) then
  mdim=minloc(pr(:,:))
  m=mdim(1)
  k=mdim(2)
  write(mpout,'(1x,2a,1p,e15.7,0p,a,2i5)') trim(where),	&
	': .01*'//trim(tag)//'(i,j,:) with .01*ps =',	&
	.01*ps(m),' at (i,j)=',m,j
  km=size(pr,2)
  do k=1,km,5
    write(mpout,'(1x,1p,5e15.7)') (.01*pr(m,km-kk+1),kk=k,min(k+4,km))
  end do
endif
end subroutine show_minval_

end module m_ppInterp
