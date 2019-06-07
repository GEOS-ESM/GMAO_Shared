!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fvSurface - Object of FV surface variables
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_fvSurface
      use m_interleavedObject,only : interleavedObject
      implicit none
      private	! except

      public :: fvSurface		! data structure
      public :: fvSurface_getheader
      public :: fvSurface_allgetheader
      public :: fvSurface_read
      public :: fvSurface_alloc_getvar
      public :: clean

      public :: get

      public :: ptr_u10m,ptr_v10m	! 10m winds
      public :: ptr_tskn		! t_skin
      public :: ptr_snow_dep		! snow-depth
      public :: ptr_soil_mois		! top layer soil wetness (frac.)
      public :: ptr_soil_temp		! top layer soil temperature
      public :: ptr_isli_mask		! island-land-ice index
      public :: ptr_sfc_rough		! surface roughness

    integer,parameter :: fvSurface_UNDEFINED=-1
    integer,parameter :: fvSurface_HEADONLY = 0
    integer,parameter :: fvSurface_DEFINED  = 1

    real,parameter :: UNDEF_SOIL_TEMP = 1.E+15
    real,parameter :: UNDEF_SNOW_DEP  = 1.E+12

    type fvSurface
      private
      integer :: state=fvSurface_UNDEFINED
      integer :: idim,jdim
      integer :: nymd,nhms
      integer :: freq
      		! idim and jdim are global dimensions of the data
		! components

      type(interleavedObject),pointer :: scat
      		! scat contains the information on how
		! the data arrays are distributed.  In
		! this implementation, an interleaved
		! levels distribution is used.

      real,pointer,dimension(:,:,:) :: u10m,v10m	! wind at 10m
      real,pointer,dimension(:,:,:) :: tskn		! t_skin
      real,pointer,dimension(:,:,:) :: snow_dep		! snow-depth
      real,pointer,dimension(:,:,:) :: soil_temp
      					! top layer soil temperature
      real,pointer,dimension(:,:,:) :: soil_mois
      					! top layer soil wetness (frac)
      real,pointer,dimension(:,:,:) :: isli_mask
      					! island/land/ice index
      real,pointer,dimension(:,:,:) :: sfc_rough
      					! surface roughness
    end type fvSurface

    interface fvSurface_getheader; module procedure head_; end interface
    interface fvSurface_allgetheader; module procedure allhead_; end interface
    interface fvSurface_read; module procedure read_; end interface
    interface fvSurface_alloc_getvar; module procedure	&
    	alloc_getvar_; end interface
    interface clean; module procedure clean_; end interface
    interface   get; module procedure   get_; end interface

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       04Feb09 - Todling - Takacs add get oro alternative
!                         - Guo add z0
!                         - Nadeau change size of char in file
!                         - Why should I be adding comments for others?
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fvSurface'

!_#define SHOW_INPUTFIELDS

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
! !IROUTINE: alloc_getvar_ - single PE read a variable from a file
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_getvar_(fname,var,tag,stat)
      use m_undef,only : undef_fv2ssi
      use m_undef,only : undef_2ssi
      use m_undef,only : undef_ssi

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_die,only : perr,die,assert_
      use m_die,only : MP_perr
      use m_die,only : MP_die
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in) :: fname
      real,pointer,dimension(:,:) :: var	! <= It is a pointer !
      character(len=*),intent(in) :: tag
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_getvar_'
  integer :: iid,rc
  integer :: im,jm,km,lm,nvars,ngatt,timinc
  integer :: myPE
!________________________________________
   logical :: fliplon
   integer :: nhmss, nymds
   integer, allocatable :: yymmdd(:)
   integer, allocatable :: hhmmss(:)
   integer, allocatable ::  kmvar(:)
   
   real :: undef_in
   real, allocatable, dimension(:,:) :: vrange, prange
   real, allocatable, dimension(:)   :: lon, lat, lev
   real, allocatable, dimension(:,:) :: tskn

   character(len = 256) :: title
   character(len = 256) :: source
   character(len = 256) :: contact
   character(len = 256) :: levunits
   character(len = 256), allocatable :: vname(:)
   character(len = 256), allocatable :: vtitle(:)
   character(len = 256), allocatable :: vunits(:)
!________________________________________

  if(present(stat)) stat=0
!________________________________________
! Open fv surface file

    call Gfio_Open ( trim(fname), 1, iid, rc )
    	if(rc/=0) then
	  call perr(myname_,	&
		'GFIO_open("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    call GFIO_DimInquire ( iid,im,jm,km,lm,nvars,ngatt,rc)
    	if(rc/=0) then
	  call perr(myname_,	&
    		'GFIO_DimInquire("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

   allocate ( lon(im) )
   allocate ( lat(jm) )
   allocate ( lev(km) )
   allocate ( yymmdd(lm) )
   allocate ( hhmmss(lm) )
   allocate ( vname(nvars) )
   allocate ( vtitle(nvars) )
   allocate ( vunits(nvars) )
   allocate (  kmvar(nvars) )
   allocate ( vrange(2,nvars) )
   allocate ( prange(2,nvars) )

    call GFIO_INQUIRE ( iid,im,jm,km,lm,nvars, &
                       title,source,contact,undef_in, &
                       lon,lat,lev,levunits, &
                       yymmdd,hhmmss,timinc, &
                       vname,vtitle,vunits,kmvar, &
                       vrange,prange,rc )
    	if(rc/=0) then
	  call perr(myname_,	&
		'GFIO_Inquire("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    fliplon=fliplon_(lon(1))
    nymds=yymmdd(1)
    nhmss=hhmmss(1)

   deallocate ( lon )
   deallocate ( lat )
   deallocate ( lev )
   deallocate ( yymmdd )
   deallocate ( hhmmss )
   deallocate ( vname )
   deallocate ( vtitle )
   deallocate ( vunits )
   deallocate (  kmvar )
   deallocate ( vrange )
   deallocate ( prange )
!________________________________________
! Define variable names.  "*" marks the variables used in output
!
! "U10M"    ! 10 meter u wind (m/s)
! "V10M"    ! 10 meter V wind (m/s)
! "TSKIN"   ! Surface skin temp (K)
! "SNOWDP"  ! Snow depth (2D)
! "ORO"     ! Surface type flag (flag)
! "GWETTOP" ! Top Soil Layer Wetness (fraction)
! "TSOIL1"  ! Top layer soil temperature (K)
! "ALDIF"   ! Albedo: longwave, diffuse
! "ALDIR"   ! Albedo: longwave, direct
! "ASDIF"   ! Albedo: shortwave, diffuse
! "ASDIR"   ! Albedo: shortwave, direct
! ** Read veg fraction, soil type and veg type from FIX file

!________________
     call mpout_log( myname_, '---get "'//trim(tag)//'"' )

     allocate(var(im,jm))

     if( trim(tag).eq.'ORO' ) then
         call get_oro ( iid,nymds,nhmss,im,jm,var,rc )
     else
         call Gfio_GetVar ( iid, trim(tag), nymds, nhmss, im, jm, 1, 1, var, rc )
     endif
     if (rc/=0) call die(myname_,'GFIO_getvar("'//trim(tag)//'")',rc)

     call undef_2ssi(var,undef_in,myname_,verb=.true.,vname=trim(tag))
!________________
! Try again for a different UNDEF value.

! When this code section was implemented, two different UNDEF values
! are used in FV surface data files.  They are 1E+12 for SNOWDP and
! 1E+15 for TSOIL1, as they are shown in plots.  Since GSI does not
! appear to give any special treatment of UNDEF values of SNOWDP (as
! snow_dep in gsi/) and TSOIL1 (as soil_temp in gsi/), it is not clear
! what may happen if all these UNDEF grid points are set to undef_ssi().

     select case(tag)
     case("TSOIl1")	! Top layer soil temperature (K)
       call undef_2ssi(var,UNDEF_SOIL_TEMP,myname_,	&
		verb=.true.,vname=trim(tag))
	if(any(var==undef_ssi())) then
     	  allocate(tskn(im,jm))
          call Gfio_GetVar ( iid, 'TSKIN', nymds, nhmss, im, jm, 1, &
                    1, tskn, rc )
	    if (rc/=0) call die(myname_,'GFIO_getvar("TSKIN")',rc)
	  where(var==undef_ssi())
	    var=tskn
	  endwhere
	  deallocate(tskn)
	endif

     case("SNOWDP")	! Snow depth (2D)
       call undef_2ssi(var,UNDEF_SNOW_DEP,myname_,	&
		verb=.true.,vname=trim(tag))
       where(var==undef_ssi())
	 var=0.
       endwhere
     end select
     if(fliplon) call hflip2_(var)
!________________
     call Gfio_Close (iid, rc)
	if (rc/=0) call die(myname_,'GFIO_close("'//trim(tag)//'")',rc)
end subroutine alloc_getvar_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: head_ - read an object header
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine head_(fname,ob,stat)
      use m_die,only : perr,die,assert_
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in ) :: fname
      type(fvSurface) ,intent(out) :: ob
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::head_'
  integer :: iid,rc,ier
  integer :: im,jm,km,lm,nvars,ngatt,timinc
!________________________________________
   integer :: nhmss, nymds, freqs
   integer, allocatable :: yymmdd(:)
   integer, allocatable :: hhmmss(:)
   integer, allocatable ::  kmvar(:)
   
   real :: undef_in
   real, allocatable, dimension(:,:) :: vrange, prange
   real, allocatable, dimension(:)   :: lon, lat, lev

   character(len = 256) :: title
   character(len = 256) :: source
   character(len = 256) :: contact
   character(len = 256) :: levunits
   character(len = 256), allocatable :: vname(:)
   character(len = 256), allocatable :: vtitle(:)
   character(len = 256), allocatable :: vunits(:)
!________________________________________
_ENTRY_

  if(present(stat)) stat=0

  ob%state=fvSurface_UNDEFINED
!________________________________________
! Open fv surface file

    call Gfio_Open ( trim(fname), 1, iid, rc )
    	if(rc/=0) then
	  call perr(myname_,	&
		'GFIO_open("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    call GFIO_DimInquire ( iid,im,jm,km,lm,nvars,ngatt,rc)
	if(rc/=0) then
	  call perr(myname_,	&
    		'GFIO_DimInquire("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    allocate ( lon(im) )
    allocate ( lat(jm) )
    allocate ( lev(km) )
    allocate ( yymmdd(lm) )
    allocate ( hhmmss(lm) )
    allocate ( vname(nvars) )
    allocate ( vtitle(nvars) )
    allocate ( vunits(nvars) )
    allocate (  kmvar(nvars) )
    allocate ( vrange(2,nvars) )
    allocate ( prange(2,nvars) )

    call GFIO_INQUIRE ( iid,im,jm,km,lm,nvars, &
                       title,source,contact,undef_in, &
                       lon,lat,lev,levunits, &
                       yymmdd,hhmmss,timinc, &
                       vname,vtitle,vunits,kmvar, &
                       vrange,prange,rc )

    	if(rc/=0) then
	  call perr(myname_,	&
		'GFIO_DimInquire("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    call Gfio_Close (iid, rc)
	if(rc/=0) then
	  call perr(myname_,	&
		'GFIO_Close("'//trim(fname)//'")',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

    ob%idim=im
    ob%jdim=jm
    ob%nymd=yymmdd(1)
    ob%nhms=hhmmss(1)
    ob%freq=timinc

    deallocate ( lon )
    deallocate ( lat )
    deallocate ( lev )
    deallocate ( yymmdd )
    deallocate ( hhmmss )
    deallocate ( vname )
    deallocate ( vtitle )
    deallocate ( vunits )
    deallocate (  kmvar )
    deallocate ( vrange )
    deallocate ( prange )
!________________________________________

  nullify(ob%scat)
  nullify(ob%u10m,ob%v10m)
  nullify(ob%tskn)
  nullify(ob%snow_dep)
  nullify(ob%soil_mois)
  nullify(ob%soil_temp)
  nullify(ob%isli_mask)
  nullify(ob%sfc_rough)

  ob%state=fvSurface_HEADONLY
_EXIT_
end subroutine head_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allhead_ - all read in a header object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allhead_(fname,ob,comm,root,stat)
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init
      use m_InterleaveScatterer,only : clean
      use m_InterleaveScattererComm,only : alloc_scatterv
      use m_interleavedObject,only : interleavedObject_init

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_die,only : perr,die,assert_
      use m_die,only : MP_perr
      use m_die,only : MP_die
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in) :: fname
      type(fvSurface),intent(out) :: ob
      integer,intent(in) :: comm
      integer,intent(in) :: root
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allhead_'
  integer :: rc,ier
  integer :: myPE
  integer :: ibufr(0:5)
!________________________________________
_ALLENTRY_

  if(present(stat)) stat=0
  call MP_comm_rank(comm,myPE,rc)
  	if(rc/=0) then
	  call MP_perr(myname_,'MP_comm_ranl()',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif
!________________________________________
! Open fv surface file on root
  if(myPE==root) then
    call head_(fname,ob,stat=rc)
    call pack_(ibufr,rc,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
    call clean_(ob)
  endif
!________________________________________

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,rc)
  	if(rc/=0) then
	  call MP_perr(myname_,'MPI_bcast()',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

  call unpack_(ibufr,rc,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
  if(rc/=0) then
    call perr(myname_,'head_()',rc)
    if(.not.present(stat)) call die(myname_)
    stat=rc
    return
  endif

  nullify(ob%scat)
  nullify(ob%u10m,ob%v10m)
  nullify(ob%tskn)
  nullify(ob%snow_dep)
  nullify(ob%soil_mois)
  nullify(ob%soil_temp)
  nullify(ob%isli_mask)
  nullify(ob%sfc_rough)

  ob%state=fvSurface_HEADONLY
_ALLEXIT_
end subroutine allhead_

subroutine pack_(ibufr,ier,idim,jdim,nymd,nhms,freq)
use m_die,only : assert_
implicit none
integer,dimension(0:),intent(out) :: ibufr
integer,intent(in) :: ier,idim,jdim,nymd,nhms,freq
	ASSERT(size(ibufr)==6)
ibufr(0:5)=(/ier,idim,jdim,nymd,nhms,freq/)
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
      use m_InterleaveScatterer,only : InterleaveScatterer
      use m_InterleaveScatterer,only : InterleaveScatterer_init
      use m_InterleaveScatterer,only : clean
      use m_InterleaveScattererComm,only : alloc_scatterv
      use m_interleavedObject,only : interleavedObject_init
      use m_undef,only : undef_fv2ssi
      use m_undef,only : undef_2ssi
      use m_undef,only : undef_ssi
!!TORM      use m_undef,only : unlike,like

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type
      use m_die,only : perr,die,assert_
      use m_die,only : MP_perr
      use m_die,only : MP_die,assert_
      use m_mpout,only : mpout_log
      implicit none
      character(len=*),intent(in) :: fname
      type(fvSurface),intent(out) :: ob
      integer,intent(in) :: comm
      integer,intent(in) :: root
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	20Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::read_'
  integer :: iid,rc
  integer :: im,jm,km,lm,nvars,ngatt,timinc,n
  integer :: myPE
  integer :: ibufr(0:5)
  real,allocatable,dimension(:,:,:) :: rbufr
  type(InterleaveScatterer) :: scatr
!________________________________________
   logical :: fliplon
   integer :: nhmss, nymds, freqs
   integer, allocatable :: yymmdd(:)
   integer, allocatable :: hhmmss(:)
   integer, allocatable ::  kmvar(:)
   
   real :: undef_in
   real, allocatable, dimension(:,:) :: vrange, prange
   real, allocatable, dimension(:)   :: lon, lat, lev

   character(len = 256) :: title
   character(len = 256) :: source
   character(len = 256) :: contact
   character(len = 256) :: levunits
   character(len = 256), allocatable :: vname(:)
   character(len = 256), allocatable :: vtitle(:)
   character(len = 256), allocatable :: vunits(:)
!________________________________________

  if(present(stat)) stat=0
  call MP_comm_rank(comm,myPE,rc)
  	if(rc/=0) then
	  call MP_perr(myname_,'MP_comm_ranl()',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

  ob%state=fvSurface_UNDEFINED

	allocate(ob%scat)
  call interleavedObject_init(ob%scat,1,comm,root)
  call InterleaveScatterer_init(scatr,1,comm,root)

  if(myPE==root) then
!________________________________________
! Open fv surface file

    call Gfio_Open ( trim(fname), 1, iid, rc )
    	if(rc/=0) call perr(myname_,	&
		'GFIO_open("'//trim(fname)//'")',rc)

    call GFIO_DimInquire ( iid,im,jm,km,lm,nvars,ngatt,rc)
    if(rc/=0) call perr(myname_,	&
    	'GFIO_DimInquire("'//trim(fname)//'")',rc)

   allocate ( lon(im) )
   allocate ( lat(jm) )
   allocate ( lev(km) )
   allocate ( yymmdd(lm) )
   allocate ( hhmmss(lm) )
   allocate ( vname(nvars) )
   allocate ( vtitle(nvars) )
   allocate ( vunits(nvars) )
   allocate (  kmvar(nvars) )
   allocate ( vrange(2,nvars) )
   allocate ( prange(2,nvars) )

    call GFIO_INQUIRE ( iid,im,jm,km,lm,nvars, &
                       title,source,contact,undef_in, &
                       lon,lat,lev,levunits, &
                       yymmdd,hhmmss,timinc, &
                       vname,vtitle,vunits,kmvar, &
                       vrange,prange,rc )

    	if(rc/=0) call perr(myname_,	&
		'GFIO_DimInquire("'//trim(fname)//'")',rc)

    call pack_(ibufr,rc,im,jm,yymmdd(1),hhmmss(1),timinc)
    fliplon=fliplon_(lon(1))

   deallocate ( lon )
   deallocate ( lat )
   deallocate ( lev )
   deallocate ( yymmdd )
   deallocate ( hhmmss )
   deallocate ( vname )
   deallocate ( vtitle )
   deallocate ( vunits )
   deallocate (  kmvar )
   deallocate ( vrange )
   deallocate ( prange )
  endif
!________________________________________

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,rc)
  	if(rc/=0) then
	  call MP_perr(myname_,'MPI_bcast()',rc)
	  if(.not.present(stat)) call die(myname_)
	  stat=rc
	  return
	endif

  call unpack_(ibufr,rc,ob%idim,ob%jdim,ob%nymd,ob%nhms,ob%freq)
  if(rc/=0) then
    call perr(myname_,'GFIO',rc)
    if(.not.present(stat)) call die(myname_)
    stat=rc
    return
  endif

  im=ob%idim
  jm=ob%jdim
  nymds=ob%nymd
  nhmss=ob%nhms

  km=0
  if(myPE==root) km=1
  allocate(rbufr(im,jm,km))

!________________________________________
! Define variable names.  "*" marks the variables used in output
!
!*"U10M"    ! 10 meter u wind (m/s)
!*"V10M"    ! 10 meter V wind (m/s)
!*"TSKIN"   ! Surface skin temp (K)
!*"SNOWDP"  ! Snow depth (2D)
!*"ORO"     ! Surface type flag (flag)
!*"GWETTOP" ! Top Soil Layer Wetness (fraction)
!*"TSOIL1"  ! Top layer soil temperature (K)
! "ALDIF"   ! Albedo: longwave, diffuse
! "ALDIR"   ! Albedo: longwave, direct
! "ASDIF"   ! Albedo: shortwave, diffuse
! "ASDIR"   ! Albedo: shortwave, direct
! ** Read veg fraction, soil type and veg type from FIX file

   	! There is actually no scattering if only one level presents.
	! However, the call is made anyway in case there will be some
	! future changes to the distribution pattern
!________________
   if(myPE==root) then
     call mpout_log( myname_, '---get U10M' )
     call Gfio_GetVar ( iid, "U10M", nymds, nhmss, im, jm, 1, &
                    1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading U10M', rc )
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='U10M')
     if(fliplon) call hflip3_(rbufr)
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%u10m,root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get V10M' )
     call Gfio_GetVar ( iid, "V10M", nymds, nhmss, im, jm, 1, &
                     1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading V10M', rc )
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='V10M')
     if(fliplon) call hflip3_(rbufr)
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%v10m,root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get TSKIN' )
     call Gfio_GetVar ( iid, "TSKIN", nymds, nhmss, im, jm, 1, &
                     1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading TSKIN', rc )
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='TSKIN')
     if(fliplon) call hflip3_(rbufr)
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%tskn,root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get ORO' )
     call get_oro ( iid,nymds,nhmss,im,jm,rbufr,rc )
     if (rc/=0) call die(myname_,'reading ORO', rc )
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='ORO')
     if(fliplon) call hflip3_(rbufr)

     	! rbufr=max(0.,min(nint(rbufr),2.))	! set to (0., 1., 2.)
     where(rbufr == undef_ssi())	! in case there is any
     elsewhere(rbufr <.5)
       rbufr=0.
     elsewhere(rbufr <1.5)
       rbufr=1.
     elsewhere
       rbufr=2.
     endwhere
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%isli_mask,root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get SNOWDP' )
     call Gfio_GetVar ( iid, "SNOWDP", nymds, nhmss, im, jm, 1,&
                     1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading SNOWDP', rc )

	! Not sure which UNDEF value is actually used in the input
	! surface data file.  Both _fv2ssi() and _2ssi() interfaces
	! were used here to ensure a standard output.
     call undef_2ssi(rbufr(:,:,1),UNDEF_SNOW_DEP,myname_,	&
	verb=.true.,vname='SNOWDP')
     call undef_2ssi(rbufr(:,:,1),undef_in      ,myname_,	&
	verb=.true.,vname='SNOWDP')
     if(fliplon) call hflip3_(rbufr)

     where(rbufr == undef_ssi())	! in case of any
       rbufr=0.
     elsewhere(rbufr <0.)
       rbufr=0.
     endwhere

   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%snow_dep,		&
   						root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get GWETTOP' )
     call Gfio_GetVar ( iid, "GWETTOP", nymds, nhmss, im, jm, 1,&
                     1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading GWETTOP', rc )

	! Expect no undef value, but checked anyway.
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='GWETTOP')
     if(fliplon) call hflip3_(rbufr)

     where(rbufr == undef_ssi())	! in case of any
     elsewhere(rbufr <=0.)		! why not <.05?
		! This block would undo any :=undef_ssi() == -9.99e33
       rbufr=.05
     elsewhere(rbufr > 1.)
		! This block would undo any :=undef_in == 1.e15
       rbufr=1.
     endwhere
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%soil_mois,	&
   						root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get TSOIL1' )
     call Gfio_GetVar ( iid, "TSOIL1", nymds, nhmss, im, jm, 1,&
                     1, rbufr, rc )
     if (rc/=0) call die(myname_,'reading TSOIL1', rc )
	! Not sure which UNDEF value is actually used in the input
	! surface data file.  Both _fv2ssi() and _2ssi() interfaces
	! were used here to ensure a standard output.
     call undef_2ssi(rbufr(:,:,1),UNDEF_SOIL_TEMP,myname_,	&
	verb=.true.,vname='TSOIL1')
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,		&
	verb=.true.,vname='TSOIL1')
     if(fliplon) call hflip3_(rbufr)

! A hack solution to workarround possible undesired undef values.
     where(rbufr == undef_ssi())
       rbufr=ob%tskn
     endwhere
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%soil_temp,	&
   						root,comm,myname_)
!________________
   if(myPE==root) then
     call mpout_log ( myname_, '---get Z0' )
     call Gfio_GetVar ( iid, "Z0", nymds, nhmss, im, jm, 1,&
                     1, rbufr, rc )
     !! if (rc/=0) call die(myname_,'reading Z0', rc )
     if (rc/=0) then
       call perr(myname_,'reading Z0', rc )
       rbufr(:,:,1)=0.
     endif

	! Expect no undef value, but checked anyway.
     call undef_2ssi(rbufr(:,:,1),undef_in,myname_,	&
	verb=.true.,vname='Z0')
     if(fliplon) call hflip3_(rbufr)
   endif

   call alloc_scatterv(im,jm,scatr,rbufr,ob%sfc_rough,	&
   						root,comm,myname_)
!________________
   call clean(scatr)
   deallocate(rbufr)
   if(myPE==root) then
     call Gfio_Close (iid, rc)
   endif

  ob%state=fvSurface_DEFINED

#ifdef SHOW_INPUTFIELDS
  call show_(myname_,ob)
#ifdef _TORM_
  n=size(ob%isli_mask)
if(n/=0.and..false.) then
  print'(2a,i6))', myname_,': size(isli_mask) = ',n
  print'(2a,i6))', myname_,				&
    ': count( sea (isli_mask == 0)) = ',count(nint(ob%isli_mask)==0)
  print'(2a,i6))', myname_,				&
    ': count(land (isli_mask == 1)) = ',count(nint(ob%isli_mask)==1)
  print'(2a,i6))', myname_,				&
    ': count( ice (isli_mask == 2)) = ',count(nint(ob%isli_mask)==2)

  print'(2a,2i6))',myname_,				&
    ':  sea points (0) with/without soil_temp ',	&
    count(nint(ob%isli_mask)==0.and.unlike(ob%soil_temp,undef_ssi())),&
    count(nint(ob%isli_mask)==0.and.  like(ob%soil_temp,undef_ssi()))
  print'(2a,2i6))',myname_,				&
    ': land points (1) with/without soil_temp ',	&
    count(nint(ob%isli_mask)==1.and.unlike(ob%soil_temp,undef_ssi())),&
    count(nint(ob%isli_mask)==1.and.  like(ob%soil_temp,undef_ssi()))
  print'(2a,2i6))',myname_,				&
    ':  ice points (2) with/without soil_temp ',	&
    count(nint(ob%isli_mask)==2.and.unlike(ob%soil_temp,undef_ssi())),&
    count(nint(ob%isli_mask)==2.and.  like(ob%soil_temp,undef_ssi()))

  print'(2a,2i6))',myname_,				&
    ':  sea points (0) with/without snow_dep ',		&
    count(nint(ob%isli_mask)==0.and.	&
		unlike(ob%snow_dep,undef_ssi())),	&
    count(nint(ob%isli_mask)==0.and.	&
		  like(ob%snow_dep,undef_ssi()))
  print'(2a,2i6))',myname_,				&
    ': land points (1) with/without snow_dep ',		&
    count(nint(ob%isli_mask)==1.and.	&
		unlike(ob%snow_dep,undef_ssi())),	&
    count(nint(ob%isli_mask)==1.and.	&
		  like(ob%snow_dep,undef_ssi()))
  print'(2a,2i6))',myname_,				&
    ':  ice points (2) with/without snow_dep ',		&
    count(nint(ob%isli_mask)==2.and.	&
		unlike(ob%snow_dep,undef_ssi())),	&
    count(nint(ob%isli_mask)==2.and.	&
		  like(ob%snow_dep,undef_ssi()))
endif
#endif
#endif
end subroutine read_

subroutine get_oro ( id,nymd,nhms,im,jm,oro,rc )
! Add by L. Takacs (comment by Todling)
implicit none
integer  nymd,nhms
integer  im,jm,id,rc
real                 oro      (im,jm)
real, allocatable :: TS       (:,:)
real, allocatable :: FRLAKE   (:,:)
real, allocatable :: FROCEAN  (:,:)
real, allocatable :: FRSEAICE (:,:)

      allocate( TS       (im,jm) )
      allocate( FRLAKE   (im,jm) )
      allocate( FROCEAN  (im,jm) )
      allocate( FRSEAICE (im,jm) )

! Create ORO Using Model Land Fractions
! -------------------------------------
                call GFIO_GetVar ( id, 'TSKIN'    , nymd, nhms, im, jm, 0, 1, TS,        rc )
      if(rc==0) call GFIO_GetVar ( id, 'FRLAKE'   , nymd, nhms, im, jm, 0, 1, FRLAKE,    rc )
      if(rc==0) call GFIO_GetVar ( id, 'FROCEAN'  , nymd, nhms, im, jm, 0, 1, FROCEAN,   rc )
      if(rc==0) call GFIO_GetVar ( id, 'FRSEAICE' , nymd, nhms, im, jm, 0, 1, FRSEAICE,  rc )
      if(rc==0) then
                                             oro = 1.0  ! Land
      where (  FROCEAN+FRLAKE >= 0.6       ) oro = 0.0  ! Water
      where (  oro==0 .and. FRSEAICE > 0.5 ) oro = 2.0  ! Ice
      where (  oro==0 .and.     TS < 271.4 ) oro = 2.0  ! Ice

      else

! READ ORO Directly from Dataset
! ------------------------------
      call Gfio_GetVar ( id, "ORO", nymd, nhms, im, jm, 0, 1, oro, rc )

      endif

      deallocate( TS        )
      deallocate( FRLAKE    )
      deallocate( FROCEAN   )
      deallocate( FRSEAICE  )

end subroutine get_oro

subroutine hflipi_(q)
    use m_die,only : assert_
    implicit none
    real,dimension(:),intent(inout) :: q
    integer :: im
    real,dimension(size(q,1)/2) :: d
    im=size(q,1)
	ASSERT(mod(im,2)==0)
    d(     1:im/2) = q(     1:im/2)
    q(     1:im/2) = q(im/2+1:im  )
    q(im/2+1:im  ) = d(     1:im/2)
end subroutine hflipi_
subroutine hflip2_(q)
    implicit none
    real,dimension(:,:),intent(inout) :: q
    integer :: j
    do j=1,size(q,2)
      call hflipi_(q(:,j))
    end do
end subroutine hflip2_
subroutine hflip3_(q)
    implicit none
    real,dimension(:,:,:),intent(inout) :: q
    integer :: k
    do k=1,size(q,3)
      call hflip2_(q(:,:,k))
    end do
end subroutine hflip3_

function fliplon_(reflon)
  use m_die,only : die
  use m_mpout,only : mpout_log
  implicit none
  real,intent(in) :: reflon
  logical :: fliplon_
  character(len=*),parameter :: myname_=myname//'::fliplon_'

  if(abs(reflon)<1.e-5) then
    fliplon_=.false.

  else if(abs(reflon+180.)<1.e-5) then
    call mpout_log(myname_,'Fliping Longitudes of Surface File')
    fliplon_=.true.

  else
    call die(myname_,'longitudes not compatible w/ neither GEOS-4 nor 5')
  endif
end function fliplon_
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
      type(fvSurface),intent(inout) :: ob

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(ob%state==fvSurface_DEFINED) then
    deallocate(ob%u10m,ob%v10m)
    deallocate(ob%tskn)
    deallocate(ob%snow_dep)
    deallocate(ob%soil_temp)
    deallocate(ob%soil_mois)
    deallocate(ob%isli_mask)
    call clean(ob%scat)
    deallocate(ob%scat)
  endif

  ob%idim=-1
  ob%jdim=-1
  ob%nymd=-1
  ob%nhms=-1
  ob%freq=-1
  ob%state=fvSurface_UNDEFINED
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
      implicit none
      type(fvSurface),intent(in) :: ob
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
  if(present(scat)) call deepcopy(ob%scat,scat)

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
      type(fvSurface),intent(in) :: ob	!

! !REVISION HISTORY:
! 	09Jun05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::show_'
_ENTRY_

call mpout_log(myname_,where)

  call checksums_show(ob%u10m,'UNKW','fv:%u10m')
  call checksums_show(ob%v10m,'UNKW','fv:%v10m')
  call checksums_show(sqrt(ob%u10m*ob%u10m+ob%v10m*ob%v10m),	&
			      'UNKW','fv:%ws10')
  call checksums_show(ob%tskn,'UNKW','fv:%tskn')
  call checksums_show(ob%snow_dep,'UNKW','fv:%snow_dep')
  call checksums_show(ob%soil_temp,'UNKW','fv:%soil_temp')
  call checksums_show(ob%soil_mois,'UNKW','fv:%soil_mois')
  call checksums_show(ob%isli_mask,'UNKW','fv:%isli_mask')
  call checksums_show(ob%sfc_rough,'UNKW','fv:%sfc_rough')
_EXIT_
end subroutine show_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ptr_u10m (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_u10m
  character(len=*),parameter :: myname_=myname//'::ptr_u10m'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_u10m => ob%u10m
end function ptr_u10m
function ptr_v10m (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_v10m
  character(len=*),parameter :: myname_=myname//'::ptr_v10m'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_v10m => ob%v10m
end function ptr_v10m
function ptr_tskn (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_tskn
  character(len=*),parameter :: myname_=myname//'::ptr_tskn'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_tskn => ob%tskn
end function ptr_tskn
function ptr_snow_dep (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_snow_dep
  character(len=*),parameter :: myname_=myname//'::ptr_snow_dep'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_snow_dep => ob%snow_dep
end function ptr_snow_dep
function ptr_soil_temp (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_soil_temp
  character(len=*),parameter :: myname_=myname//'::ptr_soil_temp'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_soil_temp => ob%soil_temp
end function ptr_soil_temp
function ptr_soil_mois (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_soil_mois
  character(len=*),parameter :: myname_=myname//'::ptr_soil_mois'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_soil_mois => ob%soil_mois
end function ptr_soil_mois
function ptr_isli_mask (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_isli_mask
  character(len=*),parameter :: myname_=myname//'::ptr_isli_mask'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_isli_mask => ob%isli_mask
end function ptr_isli_mask
function ptr_sfc_rough (ob)
  use m_die,only : die
  implicit none
  type(fvSurface),target,intent(in) :: ob
  real,pointer,dimension(:,:,:) :: ptr_sfc_rough
  character(len=*),parameter :: myname_=myname//'::ptr_sfc_rough'
  if(ob%state/=fvSurface_DEFINED) call die(myname_,'%state',ob%state)
  ptr_sfc_rough => ob%sfc_rough
end function ptr_sfc_rough
end module m_fvSurface
