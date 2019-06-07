!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fvGridThickness - pressure thickness of a fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_fvGridThickness
      implicit none
      private	! except

      public :: fvGridThickness_define	! define pressure thickness data
      public :: fvGridThickness_verify	! verify given thickness data

    interface fvGridThickness_define; module procedure	&
	define_; end interface
    interface fvGridThickness_verify; module procedure	&
	verify_; end interface

! !REVISION HISTORY:
! 	30Mar05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fvGridThickness'

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: define_ - define vertical grid in pressures thickness
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine define_(delp,ob,ps,show)
      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_die,only : assert_,die
      use m_mpout,only : mpout_log
      implicit none
      real,dimension(:,:,:),intent(out) :: delp	! pres. thickness in Pa
      type(fvGrid)       ,intent(in) :: ob
      real,dimension(:,:),intent(in) :: ps	! in Pa
      logical,optional   ,intent(in) :: show

! !REVISION HISTORY:
! 	30Mar05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::define_'
  real,allocatable,dimension(    :) :: ak,bk
  integer :: im,jm,km,k
  logical :: show_
  logical :: topdown
  real :: pr
_ENTRY_
!________________________________________
  	! Get the dimensions of the local block, either as lat-lon or
	! as lon-lat.
  im=size(ps,1)
  jm=size(ps,2)
  call fvGrid_get(ob,km=km)

  	ASSERT(im==size(delp,1))
  	ASSERT(jm==size(delp,2))
  	ASSERT(km==size(delp,3))
!________________________________________
	! Get the vertical grid references
	allocate( ak(km+1), bk(km+1) )

  call fvGrid_get(ob,ak=ak,bk=bk)	! ak is also in Pa.

  pr=sum(ps)/size(ps)
  topdown=ak(1)+bk(1)*pr < ak(km+1)+bk(km+1)*pr

  do k=1,km
    delp(:,:,k)=(ak(k+1)-ak(k)) + (bk(k+1)-bk(k))*ps(:,:)
  end do
  if(.not.topdown) delp(:,:,:)=-delp(:,:,:)
!________________________________________
  show_=.false.
  if(present(show)) show_=show
  if(show_) then
    if(topdown) then
      call showdelp_(ak(km+1)+bk(km+1)*ps(:,:),delp,topdown)
    else
      call showdelp_(ak(1)+bk(1)*ps(:,:),delp,topdown)
    endif
    call showakbk_(ak,bk,ps,myname_)
  endif

  	deallocate(ak,bk)
_EXIT_
end subroutine define_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: verify_ - verify vertical grid in pressures thickness
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine verify_(delp,ob,ps,prec,stat)
      use m_checksums,only : checksums_show
      use m_fvGrid,only : fvGrid
      use m_fvGrid,only : fvGrid_get
      use m_die,only : assert_,die,perr
      use m_mpout,only : mpout_log
      implicit none
      real,dimension(:,:,:),intent(in) :: delp
      type(fvGrid)       ,intent(in) :: ob
      real,dimension(:,:),intent(in) :: ps
      real,optional      ,intent(in) :: prec ! relative precision
      integer,optional   ,intent(out) :: stat

! !REVISION HISTORY:
! 	30Mar05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::verify_'
  real,allocatable,dimension(:,:) :: pdp,pab
  real,allocatable,dimension(:) :: ak,bk
  logical :: failedAny,failedLev,topdown
  integer :: im,jm,km,k
  real :: prec_
  real :: pr
_ENTRY_

  if(present(stat)) stat=0

  prec_=max(1.e-5,tiny(prec_))		! default precision value
  if(present(prec)) prec_=prec
  	ASSERT(prec_>=tiny(prec_))

  	! Get the dimensions of the local block, either as lat-lon or
	! as lon-lat.
  im=size(ps,1)
  jm=size(ps,2)
  call fvGrid_get(ob,km=km)

  	ASSERT(im==size(delp,1))
  	ASSERT(jm==size(delp,2))
  	ASSERT(km==size(delp,3))

  	allocate( ak(km+1), bk(km+1) )
  	allocate( pdp(im,jm), pab(im,jm) )

  call fvGrid_get(ob,ak=ak,bk=bk)

  pr=sum(ps)/size(ps)
  topdown=ak(1)+bk(1)*pr < ak(km+1)+bk(km+1)*pr
  if(.not.topdown) call die(myname_,'unexpected bottom-up levels')

  failedAny=.false.

  pab(:,:) = ak(km+1)+bk(km+1)*ps(:,:)
  pdp(:,:) = ps(:,:)
  call verifyLev_(failedLev,pab,pdp,prec_)
  failedAny=failedAny.or.failedLev

  do k=km,1,-1
    pab(:,:)=ak(k)+bk(k)*ps(:,:)
    pdp(:,:)=pdp(:,:)-delp(:,:,k)
		! check if there is anything wrong on this level.
    call verifyLev_(failedLev,pab,pdp,prec_)
    failedAny=failedAny.or.failedLev
  end do

  if(failedAny) then	! if there is any level wrong, show it all.
    call perr(myname_,'inconsistent vertical pressure grids')
    call showdelp_(ps(:,:),delp,topdown)
    call showakbk_(ak,bk,ps,myname_)
    if(.not.present(stat)) call die(myname_)
  	deallocate(ak,bk)
  	deallocate(pdp,pab)
    stat=1
    return
  endif
  	deallocate(ak,bk)
  	deallocate(pdp,pab)

_EXIT_
end subroutine verify_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: verifyLev_ - verify two grids of pressures
!
! !DESCRIPTION:
!   If any(abs(pab-pdp)>prec) show checksums of pab and pdp.
!
! !INTERFACE:

    subroutine verifyLev_(failed,pab,pdp,prec)

      implicit none
      logical,intent(out) :: failed
      real,dimension(:,:),intent(in) :: pab,pdp
      real,intent(in) :: prec

! !REVISION HISTORY:
! 	01Apr05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::verifyLev_'

  real :: dmax,d
  integer :: i,j

  failed=.false.
  dmax=0.
outerloop:	&
    do j=1,size(pab,2)
      do i=1,size(pab,1)
        d=prec
        if(pab(i,j)>0. .and. pdp(i,j)>0.) then
          d=abs(pab(i,j)-pdp(i,j))
          d=d/sqrt(pab(i,j)*pdp(i,j))
        endif
        if(dmax<d) dmax=d
        if(dmax>=prec) exit outerloop
      end do
    end do outerloop

  failed=dmax>=prec
end subroutine verifyLev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: showdelp_ - show pressure levels from thickness
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine showdelp_(ps,delp,topdown)
      use m_checksums,only : checksums_show
      implicit none
      real,dimension(:,:),intent(in) :: ps
      real,dimension(:,:,:),intent(in) :: delp
      logical,intent(in) :: topdown

! !REVISION HISTORY:
! 	01Apr05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::showdelp_'
  real,dimension(size(delp,1),size(delp,2),size(delp,3)+1) :: pres
  integer :: km,k

!!  call checksums_show(ps,'UNKW','ps')
!!  call checksums_show(delp,'UNKW','delp')

  km=size(delp,3)

  if(topdown) then
			! top level first, and surface last
    pres(:,:,km+1)=ps(:,:)
    do k=km,1,-1
      pres(:,:,k)=pres(:,:,k+1)-delp(:,:,k)
    end do

  else
			! surface level first, top level last
    pres(:,:,1)=ps(:,:)
    do k=1,km
      pres(:,:,k+1)=pres(:,:,k)-delp(:,:,k)
    end do
  endif
  call checksums_show(pres,'UNKW','pk=ps-sum_k(dp), in ps unit')
end subroutine showdelp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: showakbk_ - show pressure levels from (ak,bk)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine showakbk_(ak,bk,ps,where)
      use m_checksums,only : checksums_show
      use m_mpout,only : mpout_ison,mpout
      use m_die,only : assert_
      implicit none
      real,dimension(:),intent(in) :: ak	! (km+1)
      real,dimension(:),intent(in) :: bk	! (km+1)
      real,dimension(:,:),intent(in) :: ps	! (im,jm)
      character(len=*),intent(in) :: where

! !REVISION HISTORY:
! 	01Apr05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::showakbk_'
  real,dimension(size(ps,1),size(ps,2),size(ak)) :: pres
  integer :: km,k
  integer :: im,jm,mp
  real :: pk,pr

  km=size(ak)-1
  	ASSERT(km==size(bk)-1)

	! It is assumed in this procedure, the unit of ak is the same
	! as the unit of ps, while bk is dimensionless.  Therefore, a
	! constant reference ps is estimated from ps(:,:) to generate
	! a table of the (ak,bk) values.

  if(mpout_ison()) then
    pk=sum(ps)/size(ps)
    mp=nint(log(pk)/log(10.))
    pr=1.
    do k=1,mp
      pr=10.*pr
    end do

	! It is also assumed, the order of the (ak,bk) for fvGrid is
	! top-down, i.e. the top level is the first in the (ak,bk).
	! Therefore, while the order of table output is made in its
	! natural order, an error message is created if the order is
	! not as expected.

    write(mpout,'(2a,i4,1a)') where,':: -- nlev = ',km,' --'
    write(mpout,'(2a,2(i4,1a))') where,':: -- shape(ps) = ',	&
	size(ps,1),' X ',size(ps,2),' --'
    write(mpout,'(2a,f15.4,1a)') where,':: -- pref = ',pr,' --'
    write(mpout,'(2a)')	where,	&
	':: -- (pk,ak,pref*bk) at nlev+1 levels --'

    do k=1,km+1
      pk=ak(k)+pr*bk(k)
      write(mpout,'(i6,3f15.4)') k,pk,ak(k),pr*bk(k)
    end do
    write(mpout,'(2a)') where,'::'
  endif

  do k=1,km+1
    pres(:,:,k)=ak(k)+ps(:,:)*bk(k)
  end do
  call checksums_show(pres,'UNKW','pk=ak+bk*ps, in ps unit')
end subroutine showakbk_
end module m_fvGridThickness
