!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fvGrid - FV-Grid and communications
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_fvGrid
      use m_dyn,only : dyn_grid
      use m_dyn,only : dyn_meta
      implicit none
      private	! except

      public :: fvGrid			! data structure
      public :: fvGrid_init
      public :: fvGrid_clean, clean
      public :: fvGrid_bcast
      public :: fvGrid_get
      public :: fvGrid_show

      type fvGrid
	private
	type(dyn_grid) :: grid

	type(dyn_meta) :: phism
  	type(dyn_meta) :: hs_stdvm
  	type(dyn_meta) :: psm
  	type(dyn_meta) :: tsm
  	type(dyn_meta) :: lwim
  	type(dyn_meta) :: delpm
  	type(dyn_meta) :: um
  	type(dyn_meta) :: vm
  	type(dyn_meta) :: ptm
  	type(dyn_meta),pointer,dimension(:) :: qm

	real :: missing_value
	integer :: nymd
	integer :: nhms
	integer :: nstep
	integer :: freq
      end type fvGrid

    interface fvGrid_init ; module procedure  init_; end interface
    interface fvGrid_clean; module procedure clean_; end interface
    interface clean; module procedure clean_; end interface
    interface fvGrid_bcast; module procedure bcast_; end interface
    interface fvGrid_get  ; module procedure   get_; end interface
    interface fvGrid_show ; module procedure	&
	showdyngrid_,	&
	showwfvgrid_; end interface


! !REVISION HISTORY:
! 	25Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fvGrid'

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an object from a dyn_vect, etc.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,dyn,nymd,nhms,nstep,freq)
      use m_dyn,only : dyn_vect
      use m_mpout,only : mpout_log
      implicit none
      type(fvGrid) ,intent(out) :: ob
      type(dyn_vect),intent(in) :: dyn
      integer	    ,intent(in) :: nymd
      integer	    ,intent(in) :: nhms
      integer       ,intent(in) :: nstep
      integer       ,intent(in) :: freq

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: m,i
_ENTRY_


  				! deep-copy %grid component
  call dynGrid_dup_(dyn%grid,ob%grid)

				! copy all %meta componentS
  ob%phism   = dyn%phism
  ob%hs_stdvm= dyn%hs_stdvm
  ob%psm     = dyn%psm
  ob%tsm     = dyn%tsm
  ob%lwim    = dyn%lwim
  ob%delpm   = dyn%delpm
  ob%um      = dyn%um
  ob%vm      = dyn%vm
  ob%ptm     = dyn%ptm

	m=size(dyn%qm)
	allocate(ob%qm(m))

  do i=1,m
    ob%qm(i) = dyn%qm(i)
  end do

  				! copy additional components
  ob%missing_value=dyn%missing_value
  ob%nymd =nymd
  ob%nhms =nhms
  ob%nstep=nstep
  ob%freq =freq

_EXIT_
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dynGrid_dup_ - duplicate %grid of a dyn_vect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dynGrid_dup_(dyn,ob)
      use m_dyn,only : dyn_grid
      use m_mpout,only : mpout_log
      implicit none
      type(dyn_grid),intent(in ) :: dyn
      type(dyn_grid),intent(out) :: ob

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dynGrid_dup_'
  integer :: na,nb
_ENTRY_

  	! shallow copy
  ob=dyn

  	! deeper copy.

  na=size(ob%lat)
  nb=size(ob%lon)
	allocate(ob%lat(na))
	allocate(ob%lon(nb))
  ob%lat=dyn%lat
  ob%lon=dyn%lon

  na=size(ob%ak)
  nb=size(ob%bk)
	allocate(ob%ak(na))
	allocate(ob%bk(nb))
  ob%ak=dyn%ak
  ob%bk=dyn%bk
	
_EXIT_
end subroutine dynGrid_dup_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get local dimensions
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(ob,im,jm,km,lm,nymd,nhms,nstep,ak,bk,dyn,	&
	year,month,day,hour,freq)
      use m_dyn,only : dyn_vect
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      type(fvGrid),intent(in) :: ob
      integer,optional,intent(out) :: im
      integer,optional,intent(out) :: jm
      integer,optional,intent(out) :: km
      integer,optional,intent(out) :: lm
      integer,optional,intent(out) :: nymd
      integer,optional,intent(out) :: nhms
      integer,optional,intent(out) :: nstep
      real,dimension(:),optional,intent(out) :: ak,bk
      type(dyn_vect),optional,intent(inout) :: dyn
      integer,optional,intent(out) :: year,month,day,hour
      integer,optional,intent(out) :: freq

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: m,i
_ENTRY_

  if(present(im)) im=ob%grid%im
  if(present(jm)) jm=ob%grid%jm
  if(present(km)) km=ob%grid%km
  if(present(lm)) lm=ob%grid%lm
  if(present(nymd)) nymd=ob%nymd
  if(present(nhms)) nhms=ob%nhms
  if(present(nstep)) nstep=ob%nstep

  if(present(ak)) then
    ASSERT(size(ak)==size(ob%grid%ak))
    ak=ob%grid%ak
  endif

  if(present(bk)) then
    ASSERT(size(bk)==size(ob%grid%bk))
    bk=ob%grid%bk
  endif

  if(present(dyn)) then
    call dynGrid_dup_(ob%grid    ,dyn%grid )	! deep-copy %grid

    dyn%phism   = ob%phism	! deep-copy %meta data
    dyn%hs_stdvm= ob%hs_stdvm
    dyn%psm     = ob%psm
    dyn%tsm     = ob%tsm
    dyn%lwim    = ob%lwim
    dyn%delpm   = ob%delpm
    dyn%um      = ob%um
    dyn%vm      = ob%vm
    dyn%ptm     = ob%ptm

	m=size(ob%qm)
	allocate(dyn%qm(m))

    do i=1,m
      dyn%qm(i) = ob%qm(i)
    end do
  endif

  if(present(day  )) day  =mod(ob%nymd      ,100)
  if(present(month)) month=mod(ob%nymd/100  ,100)
  if(present(year )) year =    ob%nymd/10000
  if(present(hour )) hour =    ob%nhms/10000
  if(present(freq )) freq=ob%freq
_EXIT_
end subroutine get_
subroutine showdyngrid_(where,tag,dyn)
use m_dyn,only : dyn_vect
implicit none
character(len=*),intent(in) :: where,tag
type(dyn_vect),intent(in) :: dyn
print*,where,':: -- ',tag,' --'
print*,size(dyn%grid%lon)
print*,dyn%grid%lon_min,dyn%grid%lon_max,dyn%grid%lon_del
print*,dyn%grid%lon
print*,size(dyn%grid%lat)
print*,dyn%grid%lat_min,dyn%grid%lat_max,dyn%grid%lat_del
print*,dyn%grid%lat
print*,where,'::'
end subroutine showdyngrid_
subroutine showwfvgrid_(where,tag,wfv)
implicit none
character(len=*),intent(in) :: where,tag
type(fvGrid),intent(in) :: wfv
print*,where,':: -- ',tag,' --'
print*,size(wfv%grid%lon)
print*,wfv%grid%lon_min,wfv%grid%lon_max,wfv%grid%lon_del
print*,wfv%grid%lon
print*,size(wfv%grid%lat)
print*,wfv%grid%lat_min,wfv%grid%lat_max,wfv%grid%lat_del
print*,wfv%grid%lat
print*,where,'::'
end subroutine showwfvgrid_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: isize_ - size of integer variables in a fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function isize_(ob)
      implicit none
      type(fvGrid),intent(in) :: ob
      integer :: isize_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::isize_'

  isize_=isize_grid_(ob%grid) + 5	! +1 for size(ob%qm),
  					! +2 for ob%ymd, %hms,
					! +1 for %nstep
					! +1 for %freq
end function isize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: isize_grid_ - size of integer variables in a dyn_grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function isize_grid_(ob)
      use m_dyn,only : dyn_grid
      implicit none
      type(dyn_grid),intent(in) :: ob
      integer :: isize_grid_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::isize_grid_'

  isize_grid_=13
end function isize_grid_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rsize_ - size of real variables in a dyn_grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function rsize_(ob)
      implicit none
      type(fvGrid),intent(in) :: ob
      integer :: rsize_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rsize_'

  rsize_ = rsize_grid_(ob%grid) + 1	! +1 for ob%missing_value
end function rsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rsize_grid_ - size of real variables in a dyn_grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function rsize_grid_(ob)
      use m_dyn,only : dyn_grid
      implicit none
      type(dyn_grid),intent(in) :: ob
      integer :: rsize_grid_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rsize_grid_'

  rsize_grid_=8 + 2*(ob%km+1)
end function rsize_grid_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: csize_ - size of character variables in a fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function csize_(ob)
      implicit none
      type(fvGrid),intent(in) :: ob
      integer :: csize_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::csize_'
  integer :: i,m

  csize_ = csize_meta_(ob%phism)
  csize_ = csize_ + csize_meta_(ob%hs_stdvm)
  csize_ = csize_ + csize_meta_(ob%psm)
  csize_ = csize_ + csize_meta_(ob%tsm)
  csize_ = csize_ + csize_meta_(ob%lwim)
  csize_ = csize_ + csize_meta_(ob%delpm)
  csize_ = csize_ + csize_meta_(ob%um)
  csize_ = csize_ + csize_meta_(ob%vm)
  csize_ = csize_ + csize_meta_(ob%ptm)

  m=size(ob%qm)
  do i=1,m
    csize_ = csize_ + csize_meta_(ob%qm(i))
  end do
end function csize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: csize_meta_ - size of character variables in a dyn_meta
!
! !DESCRIPTION:
!
! !INTERFACE:

    function csize_meta_(ob)
      use m_dyn,only : dyn_meta
      implicit none
      type(dyn_meta),intent(in) :: ob
      integer :: csize_meta_

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::csize_meta_'

  csize_meta_ = len(ob%name)
  csize_meta_ = csize_meta_ + len(ob%long_name)
  csize_meta_ = csize_meta_ + len(ob%units)
end function csize_meta_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: pack_ - pack graibles in a fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine pack_(ob,ibuf,rbuf,cbuf)
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      type(fvGrid),intent(in) :: ob
      integer,dimension(:),intent(out) :: ibuf
      real   ,dimension(:),intent(out) :: rbuf
      character(len=1),dimension(:),intent(out) :: cbuf

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::pack_'
  integer :: ilen,rlen,clen
  integer :: ni,nr,nc
  integer :: i,m
_ENTRY_

  ilen=isize_(ob)
  rlen=rsize_(ob)
  clen=csize_(ob)
  	ASSERT(ilen==size(ibuf))
  	ASSERT(rlen==size(rbuf))
  	ASSERT(clen==size(cbuf))
!________________________________________

  ni=0
  nr=0
  nc=0
!________________________________________
! pack %grid

  call pack_grid_(ob%grid,ibuf,ni,rbuf,nr)
!________________________________________
! pack all %meta

  call pack_meta_(ob%phism   ,cbuf,nc)
  call pack_meta_(ob%hs_stdvm,cbuf,nc)
  call pack_meta_(ob%psm     ,cbuf,nc)
  call pack_meta_(ob%tsm     ,cbuf,nc)
  call pack_meta_(ob%lwim    ,cbuf,nc)
  call pack_meta_(ob%delpm   ,cbuf,nc)
  call pack_meta_(ob%um      ,cbuf,nc)
  call pack_meta_(ob%vm      ,cbuf,nc)
  call pack_meta_(ob%ptm     ,cbuf,nc)

  m=size(ob%qm)
  ni=ni+1
  ibuf(ni)=m

  do i=1,m
    call pack_meta_(ob%qm(i) ,cbuf,nc)
  end do
!________________________________________
! pack rest variables

  ni=ni+1
  ibuf(ni)=ob%nymd
  ni=ni+1
  ibuf(ni)=ob%nhms
  ni=ni+1
  ibuf(ni)=ob%nstep
  ni=ni+1
  ibuf(ni)=ob%freq

  nr=nr+1
  rbuf(nr)=ob%missing_value
!________________________________________

  	ASSERT(ni==ilen)
  	ASSERT(nr==rlen)
  	ASSERT(nc==clen)
_EXIT_
end subroutine pack_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: pack_grid_ - pack variables in a dyn_grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine pack_grid_(ob,ibuf,ni,rbuf,nr)
      use m_dyn,only : dyn_grid
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      type(dyn_grid),intent(in) :: ob
      integer,dimension(:),intent(inout) :: ibuf
      integer,intent(inout) :: ni
      real   ,dimension(:),intent(inout) :: rbuf
      integer,intent(inout) :: nr

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::pack_grid_'
  integer :: lr
_ENTRY_

    ibuf(ni+ 1:ni+ 3)=(/ob%ib,ob%ie,ob%im/)
    ibuf(ni+ 4:ni+ 6)=(/ob%jb,ob%je,ob%jm/)
    ibuf(ni+ 7:ni+10)=(/ob%kb,ob%ke,ob%km,ob%ks/)
    ibuf(ni+11:ni+13)=(/ob%lb,ob%le,ob%lm/)

ni=ni+13

    	ASSERT(13==isize_grid_(ob))
!________________________________________
    rbuf(nr+ 1:nr+ 3)=(/ob%lon_min,ob%lon_max,ob%lon_del/)
    rbuf(nr+ 4:nr+ 6)=(/ob%lat_min,ob%lat_max,ob%lat_del/)
    rbuf(nr+ 7:nr+ 8)=(/ob%pint,ob%ptop/)
    lr=8

    rbuf(nr+lr+1:nr+lr+ob%km+1)=ob%ak(1:ob%km+1)
    lr=lr+ob%km+1
    rbuf(nr+lr+1:nr+lr+ob%km+1)=ob%bk(1:ob%km+1)
    lr=lr+ob%km+1

nr=nr+lr

    	ASSERT(lr==rsize_grid_(ob))
_EXIT_
end subroutine pack_grid_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: pack_meta_ - pack characters in dyn_meta
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine pack_meta_(ob,cbuf,nc)
      use m_dyn,only : dyn_meta
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      type(dyn_meta),intent(in) :: ob
      character(len=1),dimension(:),intent(inout) :: cbuf
      integer,intent(inout) :: nc

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::pack_meta_'
  integer :: i
_ENTRY_

  do i=1,len(ob%name)
    nc=nc+1
    cbuf(nc)=ob%name(i:i)
  end do
  do i=1,len(ob%long_name)
    nc=nc+1
    cbuf(nc)=ob%long_name(i:i)
  end do
  do i=1,len(ob%units)
    nc=nc+1
    cbuf(nc)=ob%units(i:i)
  end do

_EXIT_
end subroutine pack_meta_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpack_ - unpack to a fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpack_(ibuf,rbuf,cbuf,ob)
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      integer,dimension(:),intent(in) :: ibuf
      real   ,dimension(:),intent(in) :: rbuf
      character(len=1),dimension(:),intent(in) :: cbuf
      type(fvGrid),intent(out) :: ob

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpack_'
  integer :: ilen,rlen,clen
  integer :: ni,nr,nc
  integer :: i,m
_ENTRY_

  ilen=size(ibuf)
  rlen=size(rbuf)
  clen=size(cbuf)

  ni=0
  nr=0
  nc=0
!________________________________________
! unpack %grid

  call unpack_grid_(ibuf,ni,rbuf,nr,ob%grid)
!________________________________________
! unpack all %meta

  call unpack_meta_(cbuf,nc,ob%phism   )
  call unpack_meta_(cbuf,nc,ob%hs_stdvm)
  call unpack_meta_(cbuf,nc,ob%psm     )
  call unpack_meta_(cbuf,nc,ob%tsm     )
  call unpack_meta_(cbuf,nc,ob%lwim    )
  call unpack_meta_(cbuf,nc,ob%delpm   )
  call unpack_meta_(cbuf,nc,ob%um      )
  call unpack_meta_(cbuf,nc,ob%vm      )
  call unpack_meta_(cbuf,nc,ob%ptm     )

  ni=ni+1
  m=ibuf(ni)
  	allocate(ob%qm(m))
  do i=1,m
    call unpack_meta_(cbuf,nc,ob%qm(i) )
  end do
!________________________________________
! unpack additional variables with checkings

  ni=ni+1
  ob%nymd =ibuf(ni)
  ni=ni+1
  ob%nhms =ibuf(ni)
  ni=ni+1
  ob%nstep=ibuf(ni)
  ni=ni+1
  ob%freq =ibuf(ni)

  nr=nr+1
  ob%missing_value=rbuf(nr)
!________________________________________

  	ASSERT(ilen==ni)
  	ASSERT(rlen==nr)
  	ASSERT(clen==nc)

	ASSERT(ilen==isize_(ob))
	ASSERT(rlen==rsize_(ob))
	ASSERT(clen==csize_(ob))
_EXIT_
end subroutine unpack_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpack_grid_ - unpack variables in a dyn_grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpack_grid_(ibuf,ni,rbuf,nr,ob)
      use m_dyn,only : dyn_grid
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      integer,dimension(:),intent(in) :: ibuf
      integer,intent(inout) :: ni
      real   ,dimension(:),intent(in) :: rbuf
      integer,intent(inout) :: nr
      type(dyn_grid),intent(out) :: ob

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpack_grid_'
  integer :: lr
_ENTRY_

ob%ib=ibuf(ni+ 1)
ob%ie=ibuf(ni+ 2)
ob%im=ibuf(ni+ 3)

ob%jb=ibuf(ni+ 4)
ob%je=ibuf(ni+ 5)
ob%jm=ibuf(ni+ 6)

ob%kb=ibuf(ni+ 7)
ob%ke=ibuf(ni+ 8)
ob%km=ibuf(ni+ 9)
ob%ks=ibuf(ni+10)

ob%lb=ibuf(ni+11)
ob%le=ibuf(ni+12)
ob%lm=ibuf(ni+13)

ni=ni+13

	ASSERT(13==isize_grid_(ob))
!________________________________________

	allocate(ob%ak(ob%km+1))
	allocate(ob%bk(ob%km+1))

ob%lon_min=rbuf(nr+ 1)
ob%lon_max=rbuf(nr+ 2)
ob%lon_del=rbuf(nr+ 3)

ob%lat_min=rbuf(nr+ 4)
ob%lat_max=rbuf(nr+ 5)
ob%lat_del=rbuf(nr+ 6)

ob%pint   =rbuf(nr+ 7)
ob%ptop   =rbuf(nr+ 8)

	lr=8
ob%ak(1:ob%km+1)=rbuf(nr+lr+1:nr+lr+ob%km+1)
	lr=lr+ob%km+1
ob%bk(1:ob%km+1)=rbuf(nr+lr+1:nr+lr+ob%km+1)
	lr=lr+ob%km+1

nr=nr+lr

	ASSERT(lr==rsize_grid_(ob))
_EXIT_
end subroutine unpack_grid_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpack_meta_ - unpack characters into dyn_meta
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpack_meta_(cbuf,nc,ob)
      use m_dyn,only : dyn_meta
      use m_die,only : assert_
      use m_mpout,only : mpout_log
      implicit none
      character(len=1),dimension(:),intent(in) :: cbuf
      integer,intent(inout) :: nc
      type(dyn_meta),intent(out) :: ob

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpack_meta_'
  integer :: i
_ENTRY_

  do i=1,len(ob%name)
    nc=nc+1
    ob%name(i:i)=cbuf(nc)
  end do
  do i=1,len(ob%long_name)
    nc=nc+1
    ob%long_name(i:i)=cbuf(nc)
  end do
  do i=1,len(ob%units)
    nc=nc+1
    ob%units(i:i)=cbuf(nc)
  end do

_EXIT_
end subroutine unpack_meta_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast fvGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine bcast_(mesg,root,comm)
      use m_mpif90,only : MP_type,MP_comm_rank
      use m_die,only : MP_die
      use m_mpout,only : mpout_log
      implicit none
      type(fvGrid),intent(inout) :: mesg	! message buffer
      integer,intent(in) :: root	! where is sent from
      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	25Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'

  integer :: ni,nr,nc
  integer,dimension(3) :: nbuf
  integer,allocatable,dimension(:) :: ibuf
  real   ,allocatable,dimension(:) :: rbuf
  character(len=1),allocatable,dimension(:) :: cbuf
  integer :: myPE
  integer :: ier
_ALLENTRY_

  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myPE==root) then
	nbuf(1)=isize_(mesg)
	nbuf(2)=rsize_(mesg)
	nbuf(3)=csize_(mesg)
  endif

  call MPI_bcast(nbuf,size(nbuf),MP_type(nbuf),root,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast(ni,nr)',ier)

	ni=nbuf(1)
	nr=nbuf(2)
	nc=nbuf(3)

  	allocate(ibuf(ni))
  	allocate(rbuf(nr))
  	allocate(cbuf(nc))

  if(myPE==root) then
    call pack_(mesg,ibuf,rbuf,cbuf)
  endif

  call MPI_bcast(ibuf,size(ibuf),MP_type(ibuf),root,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast(ibuf)',ier)

  call MPI_bcast(rbuf,size(rbuf),MP_type(rbuf),root,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast(rbuf)',ier)

  call MPI_bcast(cbuf,size(cbuf),MP_type(cbuf),root,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast(rbuf)',ier)

  if(myPE/=root) then
    call unpack_(ibuf,rbuf,cbuf,mesg)
  endif

	deallocate(ibuf)
	deallocate(rbuf)
	deallocate(cbuf)

_ALLEXIT_
end subroutine bcast_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object to a undefined state.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_mpout,only : mpout_log
      implicit none
      type(fvGrid),intent(inout) :: ob

! !REVISION HISTORY:
! 	18Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  type(fvGrid) :: undefined_fvGrid	! whatever that value is
_ENTRY_

	! first, get rid of everything dynamic

	deallocate(ob%grid%ak)
	deallocate(ob%grid%bk)
  	deallocate(ob%qm)

  ob = undefined_fvGrid
_EXIT_
end subroutine clean_
end module m_fvGrid
