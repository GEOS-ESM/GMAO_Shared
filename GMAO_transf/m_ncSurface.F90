!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ncSurface - Object of NCEP surface variables
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_ncSurface
      use m_interleavedObject,only : interleavedObject
      implicit none
      private	! except

      public :: ncSurface		! data structure
      public :: ncSurface_getheader	! read sfcio header
      public :: ncSurface_allgetheader	! read sfcio header
      public :: ncSurface_read		! read sfcio head+data
      public :: clean

      public :: get

      public :: ptr_soil_type		! soil type
      public :: ptr_veg_type		! vegetation type
      public :: ptr_veg_frac		! vegetation fraction

    integer,parameter :: ncSurface_UNDEFINED=-1	! undefined
    integer,parameter :: ncSurface_HEADONLY = 0	! header only
    integer,parameter :: ncSurface_DEFINED  = 1	! defined with data

    type ncSurface
      private
      integer :: state = ncSurface_UNDEFINED

      integer :: idim,jdim
      		! idim and jdim are global dimensions of the data
		! components

      integer :: nymd,nhms,freq
      		! date, time, and interval
      		! in yyyymmdd, hhmmss, and hhmmss

      type(interleavedObject),pointer :: scat
      		! scat contains the information on how
		! the data arrays are distributed.  In
		! this implementation, an interleaved
		! levels distribution is used.

      real,pointer,dimension(:,:,:) :: soil_type
      real,pointer,dimension(:,:,:) ::  veg_type
      real,pointer,dimension(:,:,:) ::  veg_frac
    end type ncSurface

    interface ncSurface_read; module procedure read_; end interface
    interface ncSurface_getheader; module procedure head_; end interface
    interface ncSurface_allgetheader; module procedure allhead_; end interface
    interface clean; module procedure clean_; end interface
    interface   get; module procedure   get_; end interface

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ncSurface'

!#define SHOW_INPUTFIELDS

#ifndef SHOW_INPUTFIELDS
#ifdef DEBUG_CHECKSUMS
#define SHOW_INPUTFIELDS
#endif
#endif

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: head_ - read a sfcio header
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine head_(fname,ob,stat)
      use sfcio_module,only : sfcio_head
      use sfcio_module,only : sfcio_sropen
      use sfcio_module,only : sfcio_srhead
      use sfcio_module,only : sfcio_sclose
      use sfcio_module,only : sfcio_axhead

      use m_ioutil,only : luavail
      use m_die,only : perr,die,assert_
      use m_mpout,only : mpout_log
      use m_mpout,only : mpout_ison,mpout
      implicit none

      character(len=*),intent(in) :: fname
      type(ncSurface),intent(out) :: ob
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::head_'
  integer :: ier,lu
  integer :: im,jm
  integer :: nymd,nhms,freq
  integer,dimension(0:5) :: ibufr
  type(sfcio_head) :: sfchead
!________________________________________

  if(present(stat)) stat=0
  ob%state = ncSurface_UNDEFINED
!________________________________________
! Read in NCEP SFC file (file number wired-in still)
! ---------------------

    lu=luavail()
    call sfcio_sropen(lu,fname,ier)
	if(ier/=0) then
	  call perr(myname_,'sfcio_sropen("'//trim(fname)//	&
	  	'") opening NCEP SFC file',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    call sfcio_srhead(lu,sfchead,ier)
	if(ier/=0) then
	  call perr(myname_,'sfcio_srhead("'//trim(fname)//	&
	  	'") reading NCEP SFC header',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
    call sfcio_sclose(lu,ier)
      	if(ier/=0) then
	  call perr(myname_,'sfcio_sclose("'//trim(fname)//	&
	  	'") closing NCEP SFC file',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    im=sfchead%lonb
    jm=sfchead%latb+2		! add points at two poles

    nymd=sfchead%idate(4)		! year
    nymd=sfchead%idate(2) + nymd*100	! append month
    nymd=sfchead%idate(3) + nymd*100	! append day

    nhms=sfchead%idate(1)		! hours
    nhms=00+nhms*100			! append minutes
    nhms=00+nhms*100			! append seconds

    freq=nint(sfchead%fhour)		! truncated intervals in hours
    freq=00+freq*100			! append minutes
    freq=00+freq*100			! append seconds

    ob%idim=im
    ob%jdim=jm
    ob%nymd=nymd
    ob%nhms=nhms
    ob%freq=freq

  nullify(ob%scat)
  nullify(ob%soil_type)
  nullify(ob%veg_type)
  nullify(ob%veg_frac)
  ob%state=ncSurface_HEADONLY

  call sfcio_axhead(sfchead,ier)
	if(ier/=0) then
	  call perr(myname_,'sfcio_axhead()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
end subroutine head_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allhead_ - read a sfcio header
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allhead_(fname,ob,comm,root,stat)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : interleavedObject_init
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init,clean
      use m_InterleaveScattererComm,only : alloc_scatterv

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_die,only : perr,die,assert_
      use m_die,only : MP_perr
      use m_die,only : MP_die
      use m_mpout,only : mpout_log
      use m_mpout,only : mpout_ison,mpout
      implicit none

      character(len=*),intent(in) :: fname
      type(ncSurface),intent(out) :: ob
      integer,intent(in) :: comm
      integer,intent(in) :: root
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allhead_'
  integer :: ier,lu
  integer :: myPE
  integer :: im,jm
  integer :: nymd,nhms,freq
  integer,dimension(0:5) :: ibufr
!________________________________________

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_ranl()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  ob%state = ncSurface_UNDEFINED
!________________________________________
! Read in NCEP SFC file (file number wired-in still)
! ---------------------
  if(myPE==root) then
    call head_(fname,ob,stat=ier)
    call pack_(ibufr,ier,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
    call clean(ob)
  endif

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call unpack_(ibufr,ier,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
	if(ier/=0) then
	  call perr(myname_,'sfcio',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(ob%scat)
  nullify(ob%soil_type)
  nullify(ob%veg_type)
  nullify(ob%veg_frac)
  ob%state=ncSurface_HEADONLY
end subroutine allhead_

subroutine pack_(ibufr,ier,idim,jdim,nymd,nhms,freq)
use m_die,only : assert_
implicit none
integer,dimension(0:),intent(out) :: ibufr
integer,intent(in) :: ier,idim,jdim,nymd,nhms,freq
	ASSERT(size(ibufr)==6)
  ibufr(0)=ier
  ibufr(1)=idim
  ibufr(2)=jdim
  ibufr(3)=nymd
  ibufr(4)=nhms
  ibufr(5)=freq
end subroutine pack_

subroutine unpack_(ibufr,ier,idim,jdim,nymd,nhms,freq)
use m_die,only : assert_
implicit none
integer,dimension(0:),intent(in) :: ibufr
integer,intent(out) :: ier,idim,jdim,nymd,nhms,freq
	ASSERT(size(ibufr)==6)
  ier =ibufr(0)
  idim=ibufr(1)
  jdim=ibufr(2)
  nymd=ibufr(3)
  nhms=ibufr(4)
  freq=ibufr(5)
end subroutine unpack_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: read_ - read an object from a file
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine read_(fname,ob,comm,root,stat)
      use sfcio_module,only : sfcio_data
      use sfcio_module,only : sfcio_head
      use sfcio_module,only : sfcio_srohdc
      use sfcio_module,only : sfcio_axdata
      use sfcio_module,only : sfcio_axhead

      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : interleavedObject_init
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init,clean
      use m_InterleaveScattererComm,only : alloc_scatterv

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_ioutil,only : luavail
      use m_die,only : perr,die,assert_
      use m_die,only : MP_perr
      use m_die,only : MP_die
      use m_mpout,only : mpout_log
      use m_mpout,only : mpout_ison,mpout
      implicit none

      character(len=*),intent(in) :: fname
      type(ncSurface),intent(out) :: ob
      integer,intent(in) :: comm
      integer,intent(in) :: root
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read_'
  integer :: ier,lu
  integer :: myPE
  integer :: im,jm,km
  integer :: nymd,nhms,freq
  type(sfcio_data) :: sfcdata
  type(sfcio_head) :: sfchead
  integer,dimension(0:5) :: ibufr
  real,allocatable,dimension(:,:,:) :: rbufr
  type(InterleaveScatterer) :: scatr
!________________________________________

  if(present(stat)) stat=0
  call MP_comm_rank(comm,myPE,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_ranl()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  ob%state = ncSurface_UNDEFINED
!________________________________________
! Read in NCEP SFC file (file number wired-in still)
! ---------------------

  if(myPE==root) then
    call head_(fname,ob,stat=ier)
	if(ier/=0) call perr(myname_,'head_("'//trim(fname)//'")',ier)

    if(ier==0) then
      lu=luavail()
      call sfcio_srohdc(lu,fname,sfchead,sfcdata,ier)
	if(ier/=0) call perr(myname_,'sfcio_srohdc("'//	&
		trim(fname)//'") reading NCEP SFC file',ier)
	  
      ASSERT(ob%idim==sfchead%lonb)
      ASSERT(ob%jdim==sfchead%latb+2)

      if(ier==0) then
	call sfcio_axhead(sfchead,ier)
	if(ier/=0) call perr(myname_,'sfcio_axhead()',ier)
      endif

      call pack_(ibufr,ier,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)

    	ASSERT(ob%idim==size(sfcdata%vtype,1))
    	ASSERT(ob%jdim==size(sfcdata%vtype,2)+2)
    	ASSERT(ob%idim==size(sfcdata%vfrac,1))
    	ASSERT(ob%jdim==size(sfcdata%vfrac,2)+2)
    	ASSERT(ob%idim==size(sfcdata%stype,1))
    	ASSERT(ob%jdim==size(sfcdata%stype,2)+2)

      call mpout_log(myname_,	&
	'original NCEP surface data file with lon= ',ob%idim)
      call mpout_log(myname_,	&
	'original NCEP surface data file with lat= ',ob%jdim)

      call clean(ob)
    endif
  endif

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(ibufr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call unpack_(ibufr,ier,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
	if(ier/=0) then
	  call perr(myname_,'sfcio',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  im=ob%idim
  jm=ob%jdim
  km=0
  if(myPE==root) km=1
  allocate(rbufr(im,jm,km))

  allocate(ob%scat)
  call interleavedObject_init(ob%scat,1,comm,root)
  call InterleaveScatterer_init(scatr,1,comm,root)

!________________________________________
! Note the fliping of the source North-South latitude grid direction
! to South-North, before scattering the input data.

  if(myPE==root) call setSP2NP_(rbufr(:,:,km),sfcdata%stype(:,:))
  call alloc_scatterv(im,jm,scatr,rbufr,ob%soil_type,	&
   						root,comm,myname)
  if(myPE==root) call setSP2NP_(rbufr(:,:,km),sfcdata%vtype(:,:))
  call alloc_scatterv(im,jm,scatr,rbufr,ob%veg_type,	&
   						root,comm,myname)
  if(myPE==root) call setSP2NP_(rbufr(:,:,km),sfcdata%vfrac(:,:))
  call alloc_scatterv(im,jm,scatr,rbufr,ob%veg_frac,	&
   						root,comm,myname)
  call clean(scatr)
  deallocate(rbufr)

  if(myPE==root) then
    call sfcio_axdata(sfcdata,ier)
	if(ier/=0) call perr(myname_,'sfcio_axdata()',ier)
  endif

  call MPI_bcast(ier,1,MP_type(ier),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(ier)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  ob%state=ncSurface_DEFINED

#ifdef SHOW_INPUTFIELDS
  call show_(myname_,ob)
#endif
contains
subroutine setSP2NP_(rbufr,sfcio)
  use sfcio_module, only : sfcio_realkind
  use m_die,only : assert_
  implicit none
  real,dimension(:,:),intent(out) :: rbufr
  real(kind=sfcio_realkind),dimension(:,:),intent(in ) :: sfcio

  integer :: im,jm,j

  im=size(rbufr,1)
  jm=size(rbufr,2)
  ASSERT(im  ==size(sfcio,1))
  ASSERT(im  > 0)
  ASSERT(jm-2==size(sfcio,2))

  rbufr(:, 1)=sum(sfcio(:,jm-2))/im	! add South Pole points
  do j=2,jm-1
				! for j :=    2,   3, ..., jm-2,jm-1
				!  jm-j == jm-2,jm-3, ...,    2,   1
    rbufr(:,j)=sfcio(:,jm-j) 
  end do
  rbufr(:,jm)=sum(sfcio(:,   1))/im	! add North Pole points
end subroutine setSP2NP_
end subroutine read_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ 
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ob)
      use m_interleavedObject,only : clean
      implicit none
      type(ncSurface),intent(inout) :: ob

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(ob%state==ncSurface_DEFINED) then
    deallocate(ob%soil_type)
    deallocate(ob%veg_type)
    deallocate(ob%veg_frac)
    call clean(ob%scat)
    deallocate(ob%scat)
  endif

  ob%idim=-1
  ob%jdim=-1
  ob%nymd=-1
  ob%nhms=-1
  ob%freq=-1

  ob%state = ncSurface_UNDEFINED
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get parameters
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(ob,idim,jdim,scat,nymd,nhms,freq,	&
    	year,month,day,hour)
      use m_interleavedObject,only : interleavedObject
      use m_interleavedObject,only : deepcopy
      use m_die,only : die
      implicit none
      type(ncSurface),intent(in) :: ob
      integer,optional,intent(out) :: idim,jdim
      				! idim and jdim are global dimensions
      type(interleavedObject),optional,intent(out) :: scat
      				! scat contains the information on how
				! the data arrays are distributed.  In
				! this implementation, an interleaved
				! levels distribution is used.
      integer,optional,intent(out) :: nymd,nhms,freq
		! in yyyymmdd, hhmmss, and hhmmss
      integer,optional,intent(out) :: year,month,day,hour

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'

  if(present(idim)) idim=ob%idim
  if(present(jdim)) jdim=ob%jdim
  if(present(scat)) then
    if(ob%state/=ncSurface_DEFINED)	&
    	call die(myname_,'%state',ob%state)
    call deepcopy(ob%scat,scat)
  endif

  if(present(nymd)) nymd=ob%nymd
  if(present(nhms)) nhms=ob%nhms

  if(present(day  )) day  =mod(ob%nymd      ,100)
  if(present(month)) month=mod(ob%nymd/100  ,100)
  if(present(year )) year =    ob%nymd/10000
  if(present(hour )) hour =    ob%nhms/10000
  if(present(freq )) freq=ob%freq

end subroutine get_
!-----------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: show_ - show local checksums of all state variables
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine show_(where,ob)
      use m_checksums,only : checksums_show
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in) :: where
      type(ncSurface),intent(in) :: ob	!

! !REVISION HISTORY:
! 	09Jun05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::show_'
_ENTRY_
call mpout_log(myname_,where)

  call checksums_show(ob%veg_type,'UNKW','nc:%veg_type')
  call checksums_show(ob%veg_frac,'UNKW','nc:%veg_frac')
  call checksums_show(ob%soil_type,'UNKW','nc:%soil_type')
_EXIT_
end subroutine show_
!-----------------------------------------------------------------------
function ptr_soil_type (ob)
  use m_die,only : die
  implicit none
  type(ncSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_soil_type
  character(len=*),parameter :: myname_=myname//'::ptr_soil_type'
  if(ob%state/=ncSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_soil_type => ob%soil_type
end function ptr_soil_type
function ptr_veg_type (ob)
  use m_die,only : die
  implicit none
  type(ncSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_veg_type
  character(len=*),parameter :: myname_=myname//'::ptr_veg_type'
  if(ob%state/=ncSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_veg_type => ob%veg_type
end function ptr_veg_type
function ptr_veg_frac (ob)
  use m_die,only : die
  implicit none
  type(ncSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_veg_frac
  character(len=*),parameter :: myname_=myname//'::ptr_veg_frac'
  if(ob%state/=ncSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_veg_frac => ob%veg_frac
end function ptr_veg_frac
end module m_ncSurface
