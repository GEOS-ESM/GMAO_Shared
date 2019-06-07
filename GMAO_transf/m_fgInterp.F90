!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_fgInterp - fv-grid to GSI sigma-Gaussian grid interpolator
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_fgInterp
      use m_interleavedObject   ,only : interleavedObject
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_llInterp		,only : llInterp
      use m_ppInterp		,only : MASS,THTA
      use m_ppInterp,only : VIRTUAL_TEMPERATURE
      use m_ppInterp,only : POTENTIAL_TEMPERATURE

      implicit none
      private	! except

      public :: fgInterp		! data structure

      public :: MASS,THTA
      public :: VIRTUAL_TEMPERATURE
      public :: POTENTIAL_TEMPERATURE

#ifdef _TORM_
      public :: fgInterp_init,init	! initialize an object
#endif
      public :: fgInterp_l2ginit	! for fv dyn data
      public :: fgInterp_l2linit

      public :: fgInterp_g2ginit	! not used directly
      public :: fgInterp_g2linit

      public :: fgInterp_surface_l2ginit	! for fv surface data
      public :: fgInterp_surface_l2linit
      public :: fgInterp_surface_g2ginit	! for ncep surface data
      public :: fgInterp_surface_g2linit

      public :: fgInterp_clean,clean	! clean an object

      public :: psInterp_ftog		! interpolate ps from fv to s-g
      public :: psInterp_gtof		! interpolate ps from s-g to fv

      public :: fgInterp_dtog	! interpolate from fv d-grid to s-g
      public :: fgInterp_gtod	! interpolate from s-g to fv d-grid

      public :: fgInterp_ftog	! interpolate 2/3-d grid from fv to s-g
      public :: fgInterp_gtof	! interpolate 2/3-d grid from s-g to fv

      public :: fgInterp_lh2rh	! interpolate 2-d grid from lhs to rhs
      public :: fgInterp_rh2lh	! interpolate 2-d grid from rhs to lhs

      public :: getFgrid		! get fv-grid dimensions, or
      public :: getGgrid		! sigma-Gaussian-grid dimensions

      public :: populate	! a hack to get subdomains to all_global

    type fgInterp
      private
      type(interleavedObject) :: etaLev1
      type(interleavedObject) :: etaLevs
      type(SubdomainDistributor) :: aGrid
      type(interleavedObject) :: sigLev1
      type(interleavedObject) :: sigLevs
      type(SubdomainDistributor) :: gGrid

      type(llInterp) :: llIntp
      character(len=256) :: ncepphis
    end type fgInterp

    interface fgInterp_l2ginit; module procedure	&
    	l2ginit_; end interface
    interface fgInterp_l2linit; module procedure	&
    	l2linit_; end interface
    interface fgInterp_g2ginit; module procedure	&
    	g2ginit_; end interface
    interface fgInterp_g2linit; module procedure	&
    	g2linit_; end interface

    interface fgInterp_surface_l2ginit; module procedure	&
    	surface_l2ginit_; end interface
    interface fgInterp_surface_l2linit; module procedure	&
    	surface_l2linit_; end interface
    interface fgInterp_surface_g2ginit; module procedure	&
    	surface_g2ginit_; end interface
    interface fgInterp_surface_g2linit; module procedure	&
    	surface_g2linit_; end interface

#ifdef _TORM_
    interface fgInterp_init; module procedure	&
    	init_; end interface
#endif
    interface fgInterp_clean; module procedure	&
    	clean_; end interface

#ifdef _TORM_
    interface  init; module procedure  init_; end interface
#endif
    interface clean; module procedure clean_; end interface

    interface psInterp_ftog; module procedure	&
    	psIntp_ftog_; end interface
    interface psInterp_gtof; module procedure	&
    	psIntp_gtof_; end interface

    interface fgInterp_dtog; module procedure	&
	vdtog2d_,	&
	vdtog3d_; end interface

    interface fgInterp_gtod; module procedure	&
	vgtod2d_,	&
	vgtod3d_; end interface

    interface fgInterp_ftog; module procedure	&
	satog2d_,	&
	satog3d_; end interface

    interface fgInterp_gtof; module procedure	&
	sgtoa2d_,	&
	sgtoa3d_; end interface

    interface fgInterp_lh2rh; module procedure	&
	vatog2d_,	&
	satog2d_; end interface

    interface fgInterp_rh2lh; module procedure	&
	vgtoa2d_,	&
	sgtoa2d_; end interface

    interface getFgrid; module procedure getFGrid_; end interface
    interface getGgrid; module procedure getGGrid_; end interface

    interface populate; module procedure	&
    	populate2d_; end interface

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_fgInterp'

  integer,parameter :: ROOT=0

#include "assert.H"
#include "mytrace.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: surface_l2ginit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine surface_l2ginit_(ob,iAdim,jAdim,etaLev1,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_mpout ,only : mpout_log
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::surface_l2ginit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  call l2ginit_(ob,'.none.',iAdim,jAdim,etaLev1,etaLev1,	&
  	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,1,comm)
_ALLEXIT_
end subroutine surface_l2ginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: surface_g2ginit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine surface_g2ginit_(ob,iAdim,jAdim,etaLev1,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_mpout ,only : mpout_log
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::surface_g2ginit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  call g2ginit_(ob,'.none.',iAdim,jAdim,etaLev1,etaLev1,	&
  	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,1,comm)
_ALLEXIT_
end subroutine surface_g2ginit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: surface_l2linit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine surface_l2linit_(ob,iAdim,jAdim,etaLev1,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_mpout ,only : mpout_log
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::surface_l2linit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  call l2linit_(ob,'.none.',iAdim,jAdim,etaLev1,etaLev1,	&
  	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,1,comm)
_ALLEXIT_
end subroutine surface_l2linit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: surface_g2linit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine surface_g2linit_(ob,iAdim,jAdim,etaLev1,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_mpout ,only : mpout_log
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::surface_g2linit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  call g2linit_(ob,'.none.',iAdim,jAdim,etaLev1,etaLev1,	&
  	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,1,comm)
_ALLEXIT_
end subroutine surface_g2linit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2ginit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2ginit_(ob,ncepphis,	&
    	iAdim,jAdim,etaLev1,etaLevs,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,nSigm, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_interleavedObject   ,only : deepcopy
      use m_interleavedObject   ,only : interleavedObject_init
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_SubdomainDistributor,only : show

      use m_llInterp,only : llInterp_l2ginit

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_mpout ,only : mpout_log
      use m_die   ,only : MP_die,assert_
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      character(len=*),intent(in) :: ncepphis	! filename
      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1
      type(interleavedObject),intent(in) :: etaLevs

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution
      integer,intent(in) :: nSigm	! number of sigma-levels.

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2ginit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  	! Save the filename for NCEP phis input on fv grid.

  	ASSERT(len(ob%ncepphis)>=len_trim(ncepphis))

  ob%ncepphis=ncepphis
!________________________________________
  	! define vertical fv-eta grid distribution controllers.

  call deepcopy(etaLev1,ob%etaLev1)
  call deepcopy(etaLevs,ob%etaLevs)
!________________________________________

  	call MP_comm_size(comm,nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

	ldims(1:2)=0
	call MP_dims_create(nPEs,2,ldims,ier)
  		if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  	! ldims(1:2) could be as bad as niPE X njPE = 1 X nPEs

  niPE=ldims(1)
  njPE=ldims(2)

	! partition for a 2-d working subdomain distribution.

  jPE=    myPE/niPE	! latitudinal
  iPE=mod(myPE,niPE)	! longitudinal

  call SimplePartition_(iAdim,niPE,iPE,count=iAlen,displ=iAloc)
  call SimplePartition_(jAdim,njPE,jPE,count=jAlen,displ=jAloc)

  iAloc=iAloc+1
  jAloc=jAloc+1

  	! Define an A-grid distributation operator

  call SubdomainDistributor_init(ob%aGrid,	&
  	iAdim,iAloc,iAlen,.true. ,		&	! longitude
	jAdim,jAloc,jAlen,.false., comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%aGrid")
  call show(ob%aGrid,myname_)
#endif
!________________________________________
  	! Define an horizontal lat-lon interpolation operator.

  call llInterp_l2ginit(ob%llIntp,iAdim,jAdim,iGdim,jGdim)
!________________________________________
	! define a G-grid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
  	iGdim,iGloc,iGlen,.true.,		&	! longitude
  	jGdim,jGloc,jGlen,.false. , comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%gGrid")
  call show(ob%gGrid,myname_)
#endif
!________________________________________
	! define vertical sigma-grid distribution controllers

  call interleavedobject_init(ob%sigLev1,    1,comm,root=ROOT)
  call interleavedobject_init(ob%sigLevs,nSigm,comm,root=ROOT)
_ALLEXIT_
end subroutine l2ginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2linit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2linit_(ob,ncepphis,	&
    	iAdim,jAdim,etaLev1,etaLevs,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,nSigm, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_interleavedObject   ,only : deepcopy
      use m_interleavedObject   ,only : interleavedObject_init
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_SubdomainDistributor,only : show

      use m_llInterp,only : llInterp_l2linit

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_mpout ,only : mpout_log
      use m_die   ,only : MP_die,assert_
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      character(len=*),intent(in) :: ncepphis	! filename
      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1
      type(interleavedObject),intent(in) :: etaLevs

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution
      integer,intent(in) :: nSigm	! number of sigma-levels.

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2linit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  	! Save the filename for NCEP phis input on fv grid.

  	ASSERT(len(ob%ncepphis)>=len_trim(ncepphis))

  ob%ncepphis=ncepphis
!________________________________________
  	! define vertical fv-eta grid distribution controllers.

  call deepcopy(etaLev1,ob%etaLev1)
  call deepcopy(etaLevs,ob%etaLevs)
!________________________________________

  	call MP_comm_size(comm,nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

	ldims(1:2)=0
	call MP_dims_create(nPEs,2,ldims,ier)
  		if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  	! ldims(1:2) could be as bad as niPE X njPE = 1 X nPEs

  niPE=ldims(1)
  njPE=ldims(2)

	! partition for a 2-d working subdomain distribution.

  jPE=    myPE/niPE	! latitudinal
  iPE=mod(myPE,niPE)	! longitudinal

  call SimplePartition_(iAdim,niPE,iPE,count=iAlen,displ=iAloc)
  call SimplePartition_(jAdim,njPE,jPE,count=jAlen,displ=jAloc)

  iAloc=iAloc+1
  jAloc=jAloc+1

  	! Define an A-grid distributation operator

  call SubdomainDistributor_init(ob%aGrid,	&
  	iAdim,iAloc,iAlen,.true. ,		&	! longitude
	jAdim,jAloc,jAlen,.false., comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%aGrid")
  call show(ob%aGrid,myname_)
#endif
!________________________________________
  	! Define an horizontal lat-lon interpolation operator.

  call llInterp_l2linit(ob%llIntp,iAdim,jAdim,iGdim,jGdim)
!________________________________________
	! define a G-grid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
  	iGdim,iGloc,iGlen,.true.,		&	! longitude
  	jGdim,jGloc,jGlen,.false. , comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%gGrid")
  call show(ob%gGrid,myname_)
#endif
!________________________________________
	! define vertical sigma-grid distribution controllers

  call interleavedobject_init(ob%sigLev1,    1,comm,root=ROOT)
  call interleavedobject_init(ob%sigLevs,nSigm,comm,root=ROOT)
_ALLEXIT_
end subroutine l2linit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: g2ginit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine g2ginit_(ob,ncepphis,	&
    	iAdim,jAdim,etaLev1,etaLevs,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,nSigm, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_interleavedObject   ,only : deepcopy
      use m_interleavedObject   ,only : interleavedObject_init
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_SubdomainDistributor,only : show

      use m_llInterp,only : llInterp_g2ginit

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_mpout ,only : mpout_log
      use m_die   ,only : MP_die,assert_
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      character(len=*),intent(in) :: ncepphis	! filename
      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1
      type(interleavedObject),intent(in) :: etaLevs

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution
      integer,intent(in) :: nSigm	! number of sigma-levels.

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::g2ginit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  	! Save the filename for NCEP phis input on fv grid.

  	ASSERT(len(ob%ncepphis)>=len_trim(ncepphis))

  ob%ncepphis=ncepphis
!________________________________________
  	! define vertical fv-eta grid distribution controllers.

  call deepcopy(etaLev1,ob%etaLev1)
  call deepcopy(etaLevs,ob%etaLevs)
!________________________________________

  	call MP_comm_size(comm,nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

	ldims(1:2)=0
	call MP_dims_create(nPEs,2,ldims,ier)
  		if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  	! ldims(1:2) could be as bad as niPE X njPE = 1 X nPEs

  niPE=ldims(1)
  njPE=ldims(2)

	! partition for a 2-d working subdomain distribution.

  jPE=    myPE/niPE	! latitudinal
  iPE=mod(myPE,niPE)	! longitudinal

  call SimplePartition_(iAdim,niPE,iPE,count=iAlen,displ=iAloc)
  call SimplePartition_(jAdim,njPE,jPE,count=jAlen,displ=jAloc)

  iAloc=iAloc+1
  jAloc=jAloc+1

  	! Define an A-grid distributation operator

  call SubdomainDistributor_init(ob%aGrid,	&
  	iAdim,iAloc,iAlen,.true. ,		&	! longitude
	jAdim,jAloc,jAlen,.false., comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%aGrid")
  call show(ob%aGrid,myname_)
#endif
!________________________________________
  	! Define an horizontal lat-lon interpolation operator.

  call llInterp_g2ginit(ob%llIntp,iAdim,jAdim,iGdim,jGdim)
!________________________________________
	! define a G-grid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
  	iGdim,iGloc,iGlen,.true.,		&	! longitude
  	jGdim,jGloc,jGlen,.false. , comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%gGrid")
  call show(ob%gGrid,myname_)
#endif
!________________________________________
	! define vertical sigma-grid distribution controllers

  call interleavedobject_init(ob%sigLev1,    1,comm,root=ROOT)
  call interleavedobject_init(ob%sigLevs,nSigm,comm,root=ROOT)
_ALLEXIT_
end subroutine g2ginit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: g2linit_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine g2linit_(ob,ncepphis,	&
    	iAdim,jAdim,etaLev1,etaLevs,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,nSigm, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_interleavedObject   ,only : deepcopy
      use m_interleavedObject   ,only : interleavedObject_init
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_SubdomainDistributor,only : show

      use m_llInterp,only : llInterp_g2linit

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_mpout ,only : mpout_log
      use m_die   ,only : MP_die,assert_
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      character(len=*),intent(in) :: ncepphis	! filename
      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1
      type(interleavedObject),intent(in) :: etaLevs

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution
      integer,intent(in) :: nSigm	! number of sigma-levels.

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::g2linit_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  	! Save the filename for NCEP phis input on fv grid.

  	ASSERT(len(ob%ncepphis)>=len_trim(ncepphis))

  ob%ncepphis=ncepphis
!________________________________________
  	! define vertical fv-eta grid distribution controllers.

  call deepcopy(etaLev1,ob%etaLev1)
  call deepcopy(etaLevs,ob%etaLevs)
!________________________________________

  	call MP_comm_size(comm,nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

	ldims(1:2)=0
	call MP_dims_create(nPEs,2,ldims,ier)
  		if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  	! ldims(1:2) could be as bad as niPE X njPE = 1 X nPEs

  niPE=ldims(1)
  njPE=ldims(2)

	! partition for a 2-d working subdomain distribution.

  jPE=    myPE/niPE	! latitudinal
  iPE=mod(myPE,niPE)	! longitudinal

  call SimplePartition_(iAdim,niPE,iPE,count=iAlen,displ=iAloc)
  call SimplePartition_(jAdim,njPE,jPE,count=jAlen,displ=jAloc)

  iAloc=iAloc+1
  jAloc=jAloc+1

  	! Define an A-grid distributation operator

  call SubdomainDistributor_init(ob%aGrid,	&
  	iAdim,iAloc,iAlen,.true. ,		&	! longitude
	jAdim,jAloc,jAlen,.false., comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%aGrid")
  call show(ob%aGrid,myname_)
#endif
!________________________________________
  	! Define an horizontal lat-lon interpolation operator.

  call llInterp_g2linit(ob%llIntp,iAdim,jAdim,iGdim,jGdim)
!________________________________________
	! define a G-grid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
  	iGdim,iGloc,iGlen,.true.,		&	! longitude
  	jGdim,jGloc,jGlen,.false. , comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%gGrid")
  call show(ob%gGrid,myname_)
#endif
!________________________________________
	! define vertical sigma-grid distribution controllers

  call interleavedobject_init(ob%sigLev1,    1,comm,root=ROOT)
  call interleavedobject_init(ob%sigLevs,nSigm,comm,root=ROOT)
_ALLEXIT_
end subroutine g2linit_

#ifdef _TORM_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ob,ncepphis,	&
    	iAdim,jAdim,etaLev1,etaLevs,	&
	iGdim,jGdim,iGloc,iGlen,jGloc,jGlen,nSigm, comm)

      use m_interleavedObject   ,only : interleavedObject
      use m_interleavedObject   ,only : deepcopy
      use m_interleavedObject   ,only : interleavedObject_init
      use m_SubdomainDistributor,only : SubdomainDistributor
      use m_SubdomainDistributor,only : SubdomainDistributor_init
      use m_SubdomainDistributor,only : show

      use m_llInterp,only : llInterp_init

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_dims_create
      use m_mpout ,only : mpout_log
      use m_die   ,only : MP_die,assert_
      implicit none

      type(fgInterp),intent(out) :: ob

      		! i is for longitude and j is for latitude!

      character(len=*),intent(in) :: ncepphis	! filename
      integer,intent(in) :: iAdim,jAdim	! global A-grid dimensions
      type(interleavedObject),intent(in) :: etaLev1
      type(interleavedObject),intent(in) :: etaLevs

      		! i is for longitude and j is for latitude!

      integer,intent(in) :: iGdim,jGdim	! global G-grid dimensions
      integer,intent(in) :: iGloc,iGlen	! local longitude distribution
      integer,intent(in) :: jGloc,jGlen	! local latitude distribution
      integer,intent(in) :: nSigm	! number of sigma-levels.

      integer,intent(in) :: comm	! communicator

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nPEs,myPE,ier
  integer :: iPE,niPE
  integer :: jPE,njPE
  integer :: iAloc,iAlen
  integer :: jAloc,jAlen
  integer,dimension(2) :: ldims
!________________________________________
_ALLENTRY_
  	! Save the filename for NCEP phis input on fv grid.

  	ASSERT(len(ob%ncepphis)>=len_trim(ncepphis))

  ob%ncepphis=ncepphis
!________________________________________
  	! define vertical fv-eta grid distribution controllers.

  call deepcopy(etaLev1,ob%etaLev1)
  call deepcopy(etaLevs,ob%etaLevs)
!________________________________________

  	call MP_comm_size(comm,nPEs,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
  	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

	ldims(1:2)=0
	call MP_dims_create(nPEs,2,ldims,ier)
  		if(ier/=0) call MP_die(myname_,'MP_dims_create()',ier)

  	! ldims(1:2) could be as bad as niPE X njPE = 1 X nPEs

  niPE=ldims(1)
  njPE=ldims(2)

	! partition for a 2-d working subdomain distribution.

  jPE=    myPE/niPE	! latitudinal
  iPE=mod(myPE,niPE)	! longitudinal

  call SimplePartition_(iAdim,niPE,iPE,count=iAlen,displ=iAloc)
  call SimplePartition_(jAdim,njPE,jPE,count=jAlen,displ=jAloc)

  iAloc=iAloc+1
  jAloc=jAloc+1

  	! Define an A-grid distributation operator

  call SubdomainDistributor_init(ob%aGrid,	&
  	iAdim,iAloc,iAlen,.true. ,		&	! longitude
	jAdim,jAloc,jAlen,.false., comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%aGrid")
  call show(ob%aGrid,myname_)
#endif
!________________________________________
  	! Define an horizontal lat-lon interpolation operator.

  call llInterp_init(ob%llIntp,iAdim,jAdim,iGdim,jGdim)
!________________________________________
	! define a G-grid distributation operator

  call SubdomainDistributor_init(ob%gGrid,	&
  	iGdim,iGloc,iGlen,.true.,		&	! longitude
  	jGdim,jGloc,jGlen,.false. , comm)		! latitude

#ifdef DEBUG_CHECKSUMS
  _ALLTRACE_("%gGrid")
  call show(ob%gGrid,myname_)
#endif
!________________________________________
	! define vertical sigma-grid distribution controllers

  call interleavedobject_init(ob%sigLev1,    1,comm,root=ROOT)
  call interleavedobject_init(ob%sigLevs,nSigm,comm,root=ROOT)
_ALLEXIT_
end subroutine init_
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: simplePartition_ - simple partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine simplePartition_(ngrid,nproc,iproc,count,displ)
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: ngrid	! number of total grid points
      integer,intent(in) :: nproc	! number of PEs
      integer,intent(in) :: iproc	! PE ID
      integer,intent(out) :: count	! number of local grid points
      integer,intent(out) :: displ	! displacement of the local grid

! !REVISION HISTORY:
! 	26Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::simplePartition_'
  integer :: resid
_ENTRY_

        ASSERT(nproc> 0)
        ASSERT(iproc>=0 .and. iproc< nproc)

  resid=mod(ngrid,nproc)

  count=ngrid/nproc
  if(iproc< resid) count=count+1

  displ=count*iproc
  if(iproc>=resid) displ=displ+resid

_EXIT_
end subroutine simplePartition_
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
      use m_interleavedObject	,only : clean
      use m_SubdomainDistributor,only : clean
      use m_llInterp		,only : clean
      use m_mpout,only : mpout_log
      implicit none
      type(fgInterp),intent(inout) :: ob

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
_ENTRY_

  call clean(ob%etaLev1)
  call clean(ob%etaLevs)
  call clean(ob%aGrid)
  call clean(ob%sigLev1)
  call clean(ob%sigLevs)
  call clean(ob%gGrid)
  call clean(ob%llIntp)
_EXIT_
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getFgrid_ - get dimensional parameters of the A-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getFgrid_(ob, itotalSize, jtotalSize,	&
    			     ilocalSize, jlocalSize,	&
			     iLoc, iLen, jLoc, jLen	)
      use m_SubdomainDistributor,only : get
      use m_mpout,only : mpout_log
      implicit none
      type(fgInterp),intent(in) :: ob
      integer,optional,intent(out) :: itotalsize
      integer,optional,intent(out) :: jtotalsize
      integer,optional,intent(out) :: ilocalsize
      integer,optional,intent(out) :: jlocalsize
      integer,optional,intent(out) :: iLoc, iLen
      integer,optional,intent(out) :: jLoc, jLen

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getFgrid_'
_ENTRY_

  call get(ob%aGrid,	&
  	itotalsize=itotalsize, jtotalsize =jtotalsize,	&
  	ilocalsize=ilocalsize, jlocalsize =jlocalsize,	&
	iLoc=iLoc, iLen=iLen,  jLoc=jLoc, jLen=jLen	)

_EXIT_
end subroutine getFgrid_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getGGrid_ - get dimensional parameters of the A-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getGGrid_(ob, itotalSize ,jtotalSize ,	&
    			     ilocalSize, jlocalSize,	&
			     iLoc, iLen, jLoc, jLen	)
      use m_SubdomainDistributor,only : get
      use m_mpout,only : mpout_log
      implicit none
      type(fgInterp),intent(in) :: ob
      integer,optional,intent(out) :: itotalsize
      integer,optional,intent(out) :: jtotalsize
      integer,optional,intent(out) :: ilocalsize
      integer,optional,intent(out) :: jlocalsize
      integer,optional,intent(out) :: iLoc, iLen
      integer,optional,intent(out) :: jLoc, jLen

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getGGrid_'
_ENTRY_

  call get(ob%gGrid,	&
  	itotalsize=itotalsize, jtotalsize =jtotalsize,	&
  	ilocalsize=ilocalsize, jlocalsize =jlocalsize,	&
	iLoc=iLoc, iLen=iLen,  jLoc=jLoc, jLen=jLen	)

_EXIT_
end subroutine getGGrid_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: psIntp_ftog_ - ps interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine psIntp_ftog_(ob,				&
    		ps_fvL,phis_fvL,ak_fv,bk_fv,th_fvL,	&
    	ppIntp, ps_ggS,phis_ggS,ak_gg,bk_gg, comm)

      use m_ppInterp, only : POTENTIAL_TEMPERATURE
      use m_ppInterp, only : ppInterp
      use m_ppInterp, only : ppInterp_Init
      use m_llInterp, only : llInterp_lh2rh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_interleavedObject,only : totalsize,localsize
      use m_ncepphis,only : NCEPphis_read
      use m_swapij,only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_,MP_die

      implicit none

      type(fgInterp),intent(in) :: ob

      real,dimension(:,:,:),intent(in) ::   ps_fvL
      real,dimension(:,:,:),intent(in) :: phis_fvL
      real,dimension(    :),intent(in) ::   ak_fv
      real,dimension(    :),intent(in) ::   bk_fv
      real,dimension(:,:,:),intent(in) ::   th_fvL

      type(ppInterp),intent(out) :: ppIntp
      real,dimension(:,:  ),intent(out)::   ps_ggS
      real,dimension(:,:  ),intent(out):: phis_ggS
      real,dimension(    :),intent(in) ::   ak_gg
      real,dimension(    :),intent(in) ::   bk_gg

      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::psIntp_ftog_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1,nsigLev1
  integer :: metaLev1,netaLev1
  integer :: msigLevs,nsigLevs
  integer :: metaLevs,netaLevs
  integer :: ier
  real,allocatable,dimension(:,:,:) ::  ps_SAL, phis_SAL
  real,allocatable,dimension(:,:,:) ::  ps_SGL, phis_SGL
  real,allocatable,dimension(:,:  ) :: zps_SGS,zphis_SGS
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1
	metaLevs=localsize(ob%etaLevs)
	netaLevs=totalsize(ob%etaLevs)

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1
	msigLevs=localsize(ob%sigLevs)
	nsigLevs=totalsize(ob%sigLevs)
!________________________________________
! Check for dimension consistancy of the input arguments

	ASSERT(netaLev1==1)
	ASSERT(   iAdim==size(  ps_fvL,1))
	ASSERT(   jAdim==size(  ps_fvL,2))
	ASSERT(metaLev1==size(  ps_fvL,3))
	ASSERT(   iAdim==size(phis_fvL,1))
	ASSERT(   jAdim==size(phis_fvL,2))
	ASSERT(metaLev1==size(phis_fvL,3))
	ASSERT(netaLevs==size(  ak_fv )-1)
	ASSERT(netaLevs==size(  bk_fv )-1)
	ASSERT(   iAdim==size(  th_fvL,1))
	ASSERT(   jAdim==size(  th_fvL,2))
	ASSERT(metaLevs==size(  th_fvL,3))
!________________________________________
! Check for dimension consistancy of the output arguments

	ASSERT(nsigLev1==1)
	ASSERT(   iGlen==size(  ps_ggS,2))
	ASSERT(   jGlen==size(  ps_ggS,1))
	ASSERT(   iGlen==size(phis_ggS,2))
	ASSERT(   jGlen==size(phis_ggS,1))
	ASSERT(nsigLevs==size(  ak_gg )-1)
	ASSERT(nsigLevs==size(  bk_gg )-1)

!________________________________________
  	! Get topography field, from either the given state, or
	! from a file.

  	allocate(phis_SAL(iAdim,jAdim,msigLev1))
  select case(ob%ncepphis)
  case(".none.")
    phis_SAL(:,:,:)=phis_fvL(:,:,:)
  case default
    call NCEPphis_read(ob%ncepphis,ob%sigLev1,phis_SAL,comm)
  endselect
!________________________________________

  	allocate(ps_SAL(iAdim,jAdim,msigLev1))

  call ppInterp_init(ppIntp,ob%aGrid,	&
  	ob%etaLev1,ps_fvL,phis_fvL,ak_fv,bk_fv,ob%etaLevs,th_fvL, &
	ob%sigLev1,ps_SAL,phis_SAL,ak_gg,bk_gg,ob%sigLevs,	  &
	comm,Temp=POTENTIAL_TEMPERATURE)

	! interpolate topography phis_SAL to Gaussian grid.

  	allocate(phis_SGL(iGdim,jGdim,msigLev1))
  call llInterp_lh2rh(ob%llIntp,phis_SAL,phis_SGL)
	deallocate(phis_SAL)

	allocate(zphis_SGS(iGlen,jGlen))
  call distribute(ob%sigLev1,ob%gGrid,phis_SGL,zphis_SGS,comm)
  	deallocate(phis_SGL)

  call swapij(zphis_SGS,phis_ggS)
	deallocate(zphis_SGS)

  	allocate(ps_SGL(iGdim,jGdim,msigLev1))
  call llInterp_lh2rh(ob%llIntp,ps_SAL,ps_SGL)
  	deallocate(ps_SAL)

	allocate(zps_SGS(iGlen,jGlen))
  call distribute(ob%sigLev1,ob%gGrid,ps_SGL,zps_SGS,comm)
  	deallocate(ps_SGL)

  call swapij(zps_SGS,ps_ggS)
	deallocate(zps_SGS)

_ALLEXIT_
end subroutine psIntp_ftog_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: psIntp_gtof_ - ps interpolation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine psIntp_gtof_(ob,				&
    		ps_ggS,         ak_gg,bk_gg,tv_ggS,	&
    	ppIntp, ps_fvL,phis_fvL,ak_fv,bk_fv, comm, delp )

      use m_ppInterp, only : VIRTUAL_TEMPERATURE
      use m_ppInterp, only : ppInterp
      use m_ppInterp, only : ppInterp_Init
      use m_llInterp, only : llInterp_rh2lh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : undistribute
      use m_interleavedObject,only : totalsize,localsize
      use m_ncepphis,only : NCEPphis_read
      use m_swapij,only : swapij
      use m_mpout ,only : mpout_log
      use m_die,only : assert_,MP_die

      implicit none

      type(fgInterp),intent(in) :: ob

      real,dimension(:,:),intent(in) :: ps_ggS
      real,dimension(:)  ,intent(in) :: ak_gg,bk_gg
      real,dimension(:,:,:),intent(in) :: tv_ggS

      type(ppInterp),intent(out) :: ppIntp
      real,dimension(:,:,:),intent(out) ::   ps_fvL
      real,dimension(:,:,:),intent(in ) :: phis_fvL
      real,dimension(:)  ,intent(in) :: ak_fv,bk_fv

      integer,intent(in) :: comm
      real,optional,dimension(:,:,:),intent(out) :: delp

! !REVISION HISTORY:
! 	15Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::psIntp_gtof_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1,nsigLev1
  integer :: metaLev1,netaLev1
  integer :: msigLevs,nsigLevs
  integer :: metaLevs,netaLevs
  real,allocatable,dimension(:,:  ) :: zps_SGS
  real,allocatable,dimension(:,:,:) :: zth_SGS
  real,allocatable,dimension(:,:,:) :: ps_SGL,th_SGL
  real,allocatable,dimension(:,:,:) :: ps_SAL,th_SAL,phis_SAL
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1
	metaLevs=localsize(ob%etaLevs)
	netaLevs=totalsize(ob%etaLevs)

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1
	msigLevs=localsize(ob%sigLevs)
	nsigLevs=totalsize(ob%sigLevs)
!________________________________________
! Check for dimension consistancy of the input arguments

	ASSERT(nsigLev1==1)
	ASSERT(   iGlen==size(ps_ggS,2))
	ASSERT(   jGlen==size(ps_ggS,1))
	ASSERT(nsigLevs==size(ak_gg )-1)
	ASSERT(nsigLevs==size(bk_gg )-1)
	ASSERT(   iGlen==size(tv_ggS,2))
	ASSERT(   jGlen==size(tv_ggS,1))
	ASSERT(nsigLevs==size(tv_ggS,3))
!________________________________________
! Check for dimension consistancy of the output arguments

	ASSERT(netaLev1==1)			! of 2-d surface data
	ASSERT(   iAdim==size(  ps_fvL,1))
	ASSERT(   jAdim==size(  ps_fvL,2))
	ASSERT(metaLev1==size(  ps_fvL,3))
	ASSERT(   iAdim==size(phis_fvL,1))
	ASSERT(   jAdim==size(phis_fvL,2))
	ASSERT(metaLev1==size(phis_fvL,3))
	ASSERT(netaLevs==size(  ak_fv )-1)
	ASSERT(netaLevs==size(  bk_fv )-1)
!________________________________________
	allocate(zps_SGS(iGlen,jGlen))
  call swapij(ps_ggS,zps_SGS)

  	allocate(ps_SGL(iGdim,jGdim,msigLev1))
  call undistribute(ob%sigLev1,ob%gGrid, zps_SGS,ps_SGL, comm)
	deallocate(zps_SGS)

  	allocate(ps_SAL(iAdim,jAdim,msigLev1))
  call llInterp_rh2lh(ob%llIntp,ps_SGL,ps_SAL)
  	deallocate(ps_SGL)

	allocate(zth_SGS(iGlen,jGlen,nsigLevs))
  call swapij(tv_ggS,zth_SGS)

  	allocate(th_SGL(iGdim,jGdim,msigLevs))
  call undistribute(ob%sigLevs,ob%gGrid, zth_SGS,th_SGL, comm)
	deallocate(zth_SGS)

  	allocate(th_SAL(iAdim,jAdim,msigLevs))
  call llInterp_rh2lh(ob%llIntp,th_SGL,th_SAL)
  	deallocate(th_SGL)
!________________________________________
  	! Get topography field, from either the given state, or
	! from a file.

  	allocate(phis_SAL(iAdim,jAdim,msigLev1))
  select case(ob%ncepphis)
  case(".none.")
    phis_SAL(:,:,:)=phis_fvL(:,:,:)
  case default
    call NCEPphis_read(ob%ncepphis,ob%sigLev1,phis_SAL,comm)
  endselect
!________________________________________
  call ppInterp_init(ppIntp,ob%aGrid,	&
    ob%sigLev1,ps_SAL,phis_SAL,ak_gg,bk_gg,ob%sigLevs,th_SAL,	&
    ob%etaLev1,ps_fvL,phis_fvL,ak_fv,bk_fv,ob%etaLevs,		&
    comm, Temp=VIRTUAL_TEMPERATURE,delp=delp)

	deallocate(phis_SAL)
	deallocate(  th_SAL)
	deallocate(  ps_SAL)

_ALLEXIT_
end subroutine psIntp_gtof_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: atog3d_ - from a-intlevs to sigma-Gaussian-subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine atog3d_(ob,ppIntp,w_EAL,w_SGS, key,comm,	&
    	TgtTemp, norder,nocheck,undef)
      use m_ppInterp,only : POTENTIAL_TEMPERATURE
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : ppInterp_intp
      use m_ppInterp,only : WIND,MASS,THTA
      use m_llInterp,only : llInterp_lh2rh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributorComm,only : undistribute
      use m_interleavedObject,only : totalsize,localsize
      use m_swapij, only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,dimension(:,:,:),intent(in) :: w_EAL
      real,dimension(:,:,:),intent(out):: w_SGS
      integer,intent(in) :: key		! variable class
      integer,intent(in) :: comm	! communicator
      integer,optional,intent(in) :: TgtTemp	! type of the target T
      integer,optional,intent(in) :: norder
      logical,optional,intent(in) :: nocheck
      real   ,optional,intent(in) :: undef

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::atog3d_'
  integer :: iAdim,jAdim,metaLev1,metaLevs	! EAL
  integer :: iAlen,jAlen,netaLev1,netaLevs	! EAS
  integer :: iGdim,jGdim,msigLev1,msigLevs	! SGL
  integer :: iGlen,jGlen,nsigLev1,nsigLevs	! SGS

  real,allocatable,dimension(:,:,:) :: w_EAS
  real,allocatable,dimension(:,:,:) :: w_SAS
  real,allocatable,dimension(:,:,:) :: w_SAL
  real,allocatable,dimension(:,:,:) :: w_SGL
  real,allocatable,dimension(:,:,:) :: z_SGS
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1
	metaLevs=localsize(ob%etaLevs)
	netaLevs=totalsize(ob%etaLevs)

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1
	msigLevs=localsize(ob%sigLevs)
	nsigLevs=totalsize(ob%sigLevs)
  	call get(ob%aGrid,itotalsize=iAdim,jtotalsize=jAdim)
	metalevs=localsize(ob%etaLevs)
!________________________________________
! Consistancy checking of input arguments

	ASSERT(   iAdim==size(w_EAL,1))
	ASSERT(   jAdim==size(w_EAL,2))
	ASSERT(metaLevs==size(w_EAL,3))
!________________________________________
! Consistancy checking of output arguments

	ASSERT(   iGlen==size(w_SGS,2))		! swapped i-j indices
	ASSERT(   jGlen==size(w_SGS,1))
	ASSERT(nsigLevs==size(w_SGS,3))
!________________________________________
  	! dstr. to eta-a-subdomain

  	allocate(w_EAS(iAlen,jAlen,netaLevs))
  call distribute(ob%etaLevs,ob%aGrid, w_EAL,w_EAS, comm)
!________________________________________
  	! v-intp. to sigma-a-subdomain.  (%kmsig,%si) defines the sigma
	! level grid, (%ilen,%jlen) defines the local a-subdomain.

  	allocate(w_SAS(iAlen,jAlen,nsigLevs))
  call ppInterp_intp(ppIntp,w_EAS,w_SAS, key,	&
  	SrcTemp=POTENTIAL_TEMPERATURE, TgtTemp=TgtTemp)
  	deallocate(w_EAS)
!________________________________________
	! dstr. to sigma-a-intlevs

  	allocate(w_SAL(iAdim,jAdim,msigLevs))
  call undistribute(ob%sigLevs,ob%aGrid, w_SAS,w_SAL, comm)
  	deallocate(w_SAS)
!________________________________________
  	! h-intp. to sigma-g-intlevs

	allocate(w_SGL(iGdim,jGdim,msigLevs))
  select case(key)
  case(WIND)
    call llInterp_lh2rh(ob%llIntp,w_SAL,w_SGL,	&
  	vector=.true. ,norder=norder,nocheck=nocheck,undef=undef)
  case(MASS,THTA)
    call llInterp_lh2rh(ob%llIntp,w_SAL,w_SGL,	&
  	vector=.false.,norder=norder,nocheck=nocheck,undef=undef)
  case default
    call die(myname_,'unknown variable key',key)
  end select
  	deallocate(w_SAL)
!________________________________________
	! dstr. to sigma-g-subdomain

  	allocate(z_SGS(iGlen,jGlen,nsigLevs))
  call distribute(ob%sigLevs,ob%gGrid, w_SGL,z_SGS, comm)
  	deallocate(w_SGL)
	
	! reorder the subdomain storage from (longitude,latitude) to
	! (latitude,longitude) ordering.

  call swapij(z_SGS,w_SGS)		! Swap i and j
	deallocate(z_SGS)

_ALLEXIT_
end subroutine atog3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: atog2d_ - from a-intlevs to sigma-Gaussian-subdomains
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine atog2d_(ob,w_AL,w_GS, key,comm,	&
    	norder,nocheck,undef)
      use m_ppInterp,only : WIND,MASS,THTA
      use m_llInterp,only : llInterp_lh2rh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_interleavedObject,only : totalsize,localsize
      use m_swapij, only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: w_AL
      real,dimension(:,:  ),intent(out):: w_GS
      integer,intent(in) :: key		! variable class
      integer,intent(in) :: comm	! communicator
      integer,optional,intent(in) :: norder
      logical,optional,intent(in) :: nocheck
      real   ,optional,intent(in) :: undef

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::atog2d_'
  integer :: iAdim,jAdim,metaLev1	! EAL
  integer :: iAlen,jAlen,netaLev1	! EAS
  integer :: iGdim,jGdim,msigLev1	! SGL
  integer :: iGlen,jGlen,nsigLev1	! SGS

  real,allocatable,dimension(:,:,:) :: w_GL
  real,allocatable,dimension(:,:  ) :: z_GS
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1

	ASSERT(metaLev1==msigLev1)
!________________________________________
! Consistancy checking of input arguments

	ASSERT(netaLev1==1)
	ASSERT(   iAdim==size(w_AL,1))
	ASSERT(   jAdim==size(w_AL,2))
	ASSERT(metaLev1==size(w_AL,3))
!________________________________________
! Consistancy checking of output arguments

	ASSERT(nsigLev1==1)
	ASSERT(   iGlen==size(w_GS,2))		! indices are swapped.
	ASSERT(   jGlen==size(w_GS,1))
!________________________________________
  	! h-intp. to g-intlevs

	allocate(w_GL(iGdim,jGdim,metaLev1))
  select case(key)
  case(WIND)
    call llInterp_lh2rh(ob%llIntp,w_AL,w_GL,	&
  	vector=.true. ,norder=norder,nocheck=nocheck,undef=undef)
  case(MASS)
    call llInterp_lh2rh(ob%llIntp,w_AL,w_GL,	&
  	vector=.false.,norder=norder,nocheck=nocheck,undef=undef)
  case(THTA)
    call llInterp_lh2rh(ob%llIntp,w_AL,w_GL,	&
  	vector=.false.,norder=norder,nocheck=nocheck,undef=undef)
  case default
    call die(myname_,'unknown variable key',key)
  end select
!________________________________________
	! dstr. to g-subdomain

  	allocate(z_GS(iGlen,jGlen))
  call distribute(ob%etaLev1,ob%gGrid, w_GL,z_GS, comm)
  	deallocate(w_GL)

	! reorder the subdomain storage from (longitude,latitude) to
	! (latitude,longitude) ordering.

  call swapij(z_GS,w_GS)		! Swap the index i and j
  	deallocate(z_GS)

_ALLEXIT_
end subroutine atog2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtog3d_ - vector interpolation/distribution from d-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtog3d_(ob,ppIntp,u_EDL,v_EDL,u_SGS,v_SGS, comm)
      use m_daInterp,only : daInterp_vdtoa
      use m_daInterp,only : daInterp_vatod
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,target,dimension(:,:,:),intent(in) :: u_EDL
      real,target,dimension(:,:,:),intent(in) :: v_EDL
      real,target,dimension(:,:,:),intent(out) :: u_SGS
      real,target,dimension(:,:,:),intent(out) :: v_SGS
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtog3d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLevs,nsigLevs
  integer :: metaLevs,netaLevs

  real,target,allocatable,dimension(:,:,:) :: u_EAL,v_EAL
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	netaLevs=totalsize(ob%etaLevs)
	metaLevs=localsize(ob%etaLevs)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	nsigLevs=totalsize(ob%sigLevs)
	msigLevs=localsize(ob%sigLevs)

  	ASSERT(   iAdim==size(u_EDL,1))
  	ASSERT(   jAdim==size(u_EDL,2))
  	ASSERT(metaLevs==size(u_EDL,3))
  	ASSERT(   iAdim==size(v_EDL,1))
  	ASSERT(   jAdim==size(v_EDL,2))
  	ASSERT(metaLevs==size(v_EDL,3))

  	ASSERT(   iGlen==size(u_SGS,2))
  	ASSERT(   jGlen==size(u_SGS,1))
  	ASSERT(nsigLevs==size(u_SGS,3))
  	ASSERT(   iGlen==size(v_SGS,2))
  	ASSERT(   jGlen==size(v_SGS,1))
  	ASSERT(nsigLevs==size(v_SGS,3))

  	allocate(u_EAL(iAdim,jAdim,metaLevs))
  	allocate(v_EAL(iAdim,jAdim,metaLevs))

  call daInterp_vdtoa(u_EDL,v_EDL,u_EAL,v_EAL)
  call atog3d_(ob,ppIntp,u_EAL,u_SGS, WIND,comm)
  call atog3d_(ob,ppIntp,v_EAL,v_SGS, WIND,comm)

  	deallocate(u_EAL)
  	deallocate(v_EAL)
_ALLEXIT_
end subroutine vdtog3d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vatog2d_ - vector interpolation/distribution from a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vatog2d_(ob,u_AL,v_AL,u_GS,v_GS, comm)
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: u_AL
      real,dimension(:,:,:),intent(in) :: v_AL
      real,dimension(:,:),intent(out) :: u_GS
      real,dimension(:,:),intent(out) :: v_GS
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vatog2d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1,nsigLev1
  integer :: metaLev1,netaLev1
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	netaLev1=totalsize(ob%etaLev1)
	metaLev1=localsize(ob%etaLev1)

	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	nsigLev1=totalsize(ob%sigLev1)

  	ASSERT(netaLev1==1)
  	ASSERT(   iAdim==size(u_AL,1))
  	ASSERT(   jAdim==size(u_AL,2))
  	ASSERT(metaLev1==size(u_AL,3))
  	ASSERT(   iAdim==size(v_AL,1))
  	ASSERT(   jAdim==size(v_AL,2))
  	ASSERT(metaLev1==size(v_AL,3))

  	ASSERT(nsigLev1==1)
  	ASSERT(   iGlen==size(u_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(u_GS,1))
  	ASSERT(   iGlen==size(v_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(v_GS,1))

  call atog2d_(ob,u_AL,u_GS, WIND,comm)
  call atog2d_(ob,v_AL,v_GS, WIND,comm)
_ALLEXIT_
end subroutine vatog2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdtog2d_ - vector interpolation/distribution from d-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vdtog2d_(ob,u_DL,v_DL,u_GS,v_GS, comm)
      use m_daInterp,only : daInterp_vdtoa
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: u_DL
      real,dimension(:,:,:),intent(in) :: v_DL
      real,dimension(:,:),intent(out) :: u_GS
      real,dimension(:,:),intent(out) :: v_GS
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdtog2d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1,nsigLev1
  integer :: metaLev1,netaLev1

  real,allocatable,dimension(:,:,:) :: u_AL,v_AL
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	netaLev1=totalsize(ob%etaLev1)
	metaLev1=localsize(ob%etaLev1)

	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	nsigLev1=totalsize(ob%sigLev1)

  	ASSERT(netaLev1==1)
  	ASSERT(   iAdim==size(u_DL,1))
  	ASSERT(   jAdim==size(u_DL,2))
  	ASSERT(metaLev1==size(u_DL,3))
  	ASSERT(   iAdim==size(v_DL,1))
  	ASSERT(   jAdim==size(v_DL,2))
  	ASSERT(metaLev1==size(v_DL,3))

  	ASSERT(nsigLev1==1)
  	ASSERT(   iGlen==size(u_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(u_GS,1))
  	ASSERT(   iGlen==size(v_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(v_GS,1))

  	allocate(u_AL(iAdim,jAdim,metaLev1))
  	allocate(v_AL(iAdim,jAdim,metaLev1))

  call daInterp_vdtoa(u_DL,v_DL,u_AL,v_AL)
  call atog2d_(ob,u_AL,u_GS, WIND,comm)
  call atog2d_(ob,v_AL,v_GS, WIND,comm)

  	deallocate(u_AL)
  	deallocate(v_AL)
_ALLEXIT_
end subroutine vdtog2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: satog3d_ - scalar interpolation/distribution from a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine satog3d_(ob,ppIntp,w_EAL,w_SGS, key,comm, TgtTemp)
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : MASS,THTA,WIND
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,dimension(:,:,:),intent(in ) :: w_EAL
      real,dimension(:,:,:),intent(out) :: w_SGS
      integer,intent(in) :: key
      integer,intent(in) :: comm
      integer,optional,intent(in) :: TgtTemp	! type of the target th

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::satog3d_'
  integer :: ier
_ALLENTRY_

  select case(key)
  case(MASS)
    call atog3d_(ob,ppIntp,w_EAL,w_SGS, key,comm, norder=1)
  case(THTA)
    call atog3d_(ob,ppIntp,w_EAL,w_SGS, key,comm, TgtTemp=TgtTemp)
  case(WIND)
    call die(myname_,'use vdtog3d_() for WIND')
  case default
    call die(myname_,'unknown key value',key)
  end select
_ALLEXIT_
end subroutine satog3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: satog2d_ - scalar interpolation/distribution from a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine satog2d_(ob,w_AL,w_GS, key,comm)
      use m_ppInterp,only : MASS,THTA,WIND
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in ) :: w_AL
      real,dimension(:,:  ),intent(out) :: w_GS
      integer,intent(in) :: key
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::satog2d_'
  integer :: ier
_ALLENTRY_

  select case(key)
  case(MASS)
    call atog2d_(ob,w_AL,w_GS, key,comm,norder=1)
  case(THTA)
    call atog2d_(ob,w_AL,w_GS, key,comm)
  case(WIND)
    call die(myname_,'use vdtog2d_() for WIND')
  case default
    call die(myname_,'unknown key value',key)
  end select
_ALLEXIT_
end subroutine satog2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gtoa3d_ - inverse of atog3d_()
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gtoa3d_(ob,ppIntp, w_SGS,w_EAL, key,comm,	&
    	SrcTemp, norder,nocheck,undef)
      use m_ppInterp,only : POTENTIAL_TEMPERATURe
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : ppInterp_intp
      use m_ppInterp,only : WIND,MASS,THTA
      use m_llInterp,only : llInterp_rh2lh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : distribute
      use m_SubdomainDistributorComm,only : undistribute
      use m_interleavedObject,only : totalsize,localsize
      use m_swapij, only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,dimension(:,:,:),intent(in ) :: w_SGS
      real,dimension(:,:,:),intent(out) :: w_EAL
      integer,intent(in) :: key
      integer,intent(in) :: comm
      integer,optional,intent(in) :: SrcTemp	! type of the source T
      integer,optional,intent(in) :: norder
      logical,optional,intent(in) :: nocheck
      real   ,optional,intent(in) :: undef

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gtoa3d_'
  integer :: iAdim,jAdim,metaLev1,metaLevs	! EAL
  integer :: iAlen,jAlen,netaLev1,netaLevs	! EAS
  integer :: iGdim,jGdim,msigLev1,msigLevs	! SGL
  integer :: iGlen,jGlen,nsigLev1,nsigLevs	! SGS

  real,allocatable,dimension(:,:,:) :: w_EAS
  real,allocatable,dimension(:,:,:) :: w_SAS
  real,allocatable,dimension(:,:,:) :: w_SAL
  real,allocatable,dimension(:,:,:) :: w_SGL
  real,allocatable,dimension(:,:,:) :: z_SGS
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1
	metaLevs=localsize(ob%etaLevs)
	netaLevs=totalsize(ob%etaLevs)

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1
	msigLevs=localsize(ob%sigLevs)
	nsigLevs=totalsize(ob%sigLevs)
!________________________________________
! Consistancy checking of input arguments

	ASSERT(   iGlen==size(w_SGS,2))		! indices are swapped.
	ASSERT(   jGlen==size(w_SGS,1))
	ASSERT(nsigLevs==size(w_SGS,3))
!________________________________________
! Consistancy checking of output arguments

	ASSERT(   iAdim==size(w_EAL,1))
	ASSERT(   jAdim==size(w_EAL,2))
	ASSERT(metaLevs==size(w_EAL,3))
!________________________________________
	! un-dstr. to sigma-g-intlevs

	! reorder the subdomain storage from (latitude,longitude) to
	! (longitude,latitude) ordering.

	allocate(z_SGS(iGlen,jGlen,nsigLevs))
  call swapij(w_SGS,z_SGS)

	allocate(w_SGL(iGdim,jGdim,msigLevs))
  call undistribute(ob%sigLevs,ob%gGrid, z_SGS,w_SGL, comm)
  	deallocate(z_SGS)
!________________________________________
  	! h-intp. to sigma-a-intlevs

	allocate(w_SAL(iAdim,jAdim,msigLevs))
  select case(key)
  case(WIND)
    call llInterp_rh2lh(ob%llIntp,w_SGL,w_SAL,	&
  	vector=.true. ,norder=norder,nocheck=nocheck,undef=undef)
  case(MASS,THTA)
    call llInterp_rh2lh(ob%llIntp,w_SGL,w_SAL,	&
  	vector=.false.,norder=norder,nocheck=nocheck,undef=undef)
  case default
    call die(myname_,'unknown variable key',key)
  end select
  	deallocate(w_SGL)
!________________________________________
	! dstr. to sigma-a-subdomains

  	allocate(w_SAS(iAlen,jAlen,nsigLevs))
  call distribute(ob%sigLevs,ob%aGrid, w_SAL,w_SAS, comm)
  	deallocate(w_SAL)
!________________________________________
  	! v-intp. to eta-a-subdomain.  (%km_eta,%ps,%ds) defines the
	! eta level grid.

  	allocate(w_EAS(iAlen,jAlen,netaLevs))
  call ppInterp_intp(ppIntp,w_SAS,w_EAS,key,	&
  	SrcTemp=SrcTemp,TgtTemp=POTENTIAL_TEMPERATURE)
  	deallocate(w_SAS)
!________________________________________
  	! un-dstr. to eta-a-intlevs

  call undistribute(ob%etaLevs,ob%aGrid, w_EAS,w_EAL, comm)
  	deallocate(w_EAS)

_ALLEXIT_
end subroutine gtoa3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gtoa2d_ - inverse of atog3d_()
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gtoa2d_(ob, w_GS,w_AL, key,comm,	&
    	norder,nocheck,undef)
      use m_ppInterp,only : WIND,MASS,THTA
      use m_llInterp,only : llInterp_rh2lh
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : undistribute
      use m_interleavedObject,only : totalsize,localsize
      use m_swapij, only : swapij
      use m_mpout,only : mpout_log
      use m_die,only : assert_,die
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:  ),intent(in ) :: w_GS
      real,dimension(:,:,:),intent(out) :: w_AL
      integer,intent(in) :: key
      integer,intent(in) :: comm
      integer,optional,intent(in) :: norder
      logical,optional,intent(in) :: nocheck
      real   ,optional,intent(in) :: undef

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gtoa2d_'
  integer :: iAdim,jAdim,metaLev1,metaLevs	! EAL
  integer :: iAlen,jAlen,netaLev1,netaLevs	! EAS
  integer :: iGdim,jGdim,msigLev1,msigLevs	! SGL
  integer :: iGlen,jGlen,nsigLev1,nsigLevs	! SGS

  real,allocatable,dimension(:,:,:) :: w_GL
  real,allocatable,dimension(:,:) :: z_GS
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1

	ASSERT(metaLev1==msigLev1)
!________________________________________
! Consistancy checking of input arguments

	ASSERT(nsigLev1==1)
	ASSERT(   iGlen==size(w_GS,2))		! indices are swapped.
	ASSERT(   jGlen==size(w_GS,1))
!________________________________________
! Consistancy checking of output arguments

	ASSERT(   iAdim==size(w_AL,1))
	ASSERT(   jAdim==size(w_AL,2))
	ASSERT(metaLev1==size(w_AL,3))

!________________________________________
	! un-dstr. to g-intlevs

	! reorder the subdomain storage from (latitude,longitude) to
	! (longitude,latitude) ordering.

	allocate(z_GS(iGlen,jGlen))
  call swapij(w_GS,z_GS)

	allocate(w_GL(iGdim,jGdim,msigLev1))
  call undistribute(ob%sigLev1,ob%gGrid, z_GS,w_GL, comm)
  	deallocate(z_GS)
!________________________________________
  	! h-intp. to a-intlevs

  select case(key)
  case(WIND)
    call llInterp_rh2lh(ob%llIntp,w_GL,w_AL,	&
  	vector=.true. ,norder=norder,nocheck=nocheck,undef=undef)
  case(MASS,THTA)
    call llInterp_rh2lh(ob%llIntp,w_GL,w_AL,	&
  	vector=.false.,norder=norder,nocheck=nocheck,undef=undef)
  case default
    call die(myname_,'unknown variable key',key)
  end select
  	deallocate(w_GL)

_ALLEXIT_
end subroutine gtoa2d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: populate2d_ - populate subdomains to all
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine populate2d_(ob, w_GS,w_all, comm)
      use m_SubdomainDistributor,only : get
      use m_SubdomainDistributorComm,only : undistribute
      use m_interleavedObject,only : totalsize,localsize
      use m_swapij, only : swapij
      use m_mpout,only : mpout_log
      use m_mpif90,only : MP_comm_rank,MP_type
      use m_die,only : assert_,die,MP_die
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:),intent(in ) :: w_GS
      real,dimension(:,:),intent(out) :: w_all
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::populate2d_'
  integer :: iAdim,jAdim,metaLev1,metaLevs	! EAL
  integer :: iAlen,jAlen,netaLev1,netaLevs	! EAS
  integer :: iGdim,jGdim,msigLev1,msigLevs	! SGL
  integer :: iGlen,jGlen,nsigLev1,nsigLevs	! SGS

  real,allocatable,dimension(:,:,:) :: w_GL
  real,allocatable,dimension(:,:) :: z_GS
  integer :: ier,myPE
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	call get(ob%aGrid,ilocalSize=iAlen,jlocalSize=jAlen)
	metaLev1=localsize(ob%etaLev1)	! msigLev1 = 0 or 1
	netaLev1=totalsize(ob%etaLev1)	! nsigLev1 ==1

	call get(ob%gGrid,itotalSize=iGdim,jtotalSize=jGdim)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	msigLev1=localsize(ob%sigLev1)	! msigLev1 = 0 or 1
	nsigLev1=totalsize(ob%sigLev1)	! nsigLev1 ==1

	ASSERT(metaLev1==msigLev1)
!________________________________________
	call MP_comm_rank(comm,myPE,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
!________________________________________
! Consistancy checking of input arguments

	ASSERT(nsigLev1==1)
	ASSERT(   iGlen==size(w_GS,2))		! indices are swapped.
	ASSERT(   jGlen==size(w_GS,1))

	if(myPE==ROOT) then
	  ASSERT(msigLev1==1)
	endif
!________________________________________
! Consistancy checking of output arguments

	ASSERT(   iGdim==size(w_all,2))		! also swapped
	ASSERT(   jGdim==size(w_all,1))
!________________________________________
	! un-dstr. to g-intlevs

	! reorder the subdomain storage from (latitude,longitude) to
	! (longitude,latitude) ordering.

	allocate(z_GS(iGlen,jGlen))
  call swapij(w_GS,z_GS)

	allocate(w_GL(iGdim,jGdim,msigLev1))
  call undistribute(ob%sigLev1,ob%gGrid, z_GS,w_GL, comm)
  	deallocate(z_GS)

  if(myPE==ROOT) call swapij(w_GL(:,:,1),w_all)
  	deallocate(w_GL)

  call MPI_bcast(w_all,size(w_all),MP_type(w_all),ROOT,comm,ier)
  	if(ier/=0) call MP_die(myname_,'MPI_bcast()',ier)
_ALLEXIT_
end subroutine populate2d_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vgtod3d_ - vector interpolation/distribution to d-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vgtod3d_(ob,ppIntp,u_SGS,v_SGS,u_EDL,v_EDL, comm)
      use m_daInterp,only : daInterp_vatod
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,dimension(:,:,:),intent(in) :: u_SGS
      real,dimension(:,:,:),intent(in) :: v_SGS
      real,dimension(:,:,:),intent(out) :: u_EDL
      real,dimension(:,:,:),intent(out) :: v_EDL
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vgtod3d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLevs,nsigLevs
  integer :: metaLevs,netaLevs
  integer :: ier

  real,allocatable,dimension(:,:,:) :: u_EAL,v_EAL
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	metaLevs=localsize(ob%etaLevs)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)
	nsigLevs=totalsize(ob%sigLevs)

  	ASSERT(   iAdim==size(u_EDL,1))
  	ASSERT(   jAdim==size(u_EDL,2))
  	ASSERT(metaLevs==size(u_EDL,3))
  	ASSERT(   iAdim==size(v_EDL,1))
  	ASSERT(   jAdim==size(v_EDL,2))
  	ASSERT(metaLevs==size(v_EDL,3))

  	ASSERT(   iGlen==size(u_SGS,2))
  	ASSERT(   jGlen==size(u_SGS,1))
  	ASSERT(nsigLevs==size(u_SGS,3))
  	ASSERT(   iGlen==size(v_SGS,2))
  	ASSERT(   jGlen==size(v_SGS,1))
  	ASSERT(nsigLevs==size(v_SGS,3))

  	allocate(u_EAL(iAdim,jAdim,metaLevs))
  	allocate(v_EAL(iAdim,jAdim,metaLevs))

  call gtoa3d_(ob,ppIntp,u_SGS,u_EAL, WIND,comm)
  call gtoa3d_(ob,ppIntp,v_SGS,v_EAL, WIND,comm)

  call daInterp_vatod(u_EAL,v_EAL,u_EDL,v_EDL)
  	deallocate(u_EAL)
  	deallocate(v_EAL)
_ALLEXIT_
end subroutine vgtod3d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vgtoa2d_ - vector interpolation/distribution to a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vgtoa2d_(ob,u_GS,v_GS,u_AL,v_AL, comm)
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:  ),intent(in) :: u_GS
      real,dimension(:,:  ),intent(in) :: v_GS
      real,dimension(:,:,:),intent(out) :: u_AL
      real,dimension(:,:,:),intent(out) :: v_AL
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vgtoa2d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1
  integer :: metaLev1
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	metaLev1=localsize(ob%etaLev1)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)

  	ASSERT(   iGlen==size(u_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(u_GS,1))
  	ASSERT(   iGlen==size(v_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(v_GS,1))

  	ASSERT(   iAdim==size(u_AL,1))
  	ASSERT(   jAdim==size(u_AL,2))
  	ASSERT(metaLev1==size(u_AL,3))
  	ASSERT(   iAdim==size(v_AL,1))
  	ASSERT(   jAdim==size(v_AL,2))
  	ASSERT(metaLev1==size(v_AL,3))

  call gtoa2d_(ob,u_GS,u_AL, WIND,comm)
  call gtoa2d_(ob,v_GS,v_AL, WIND,comm)
_ALLEXIT_
end subroutine vgtoa2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vgtod2d_ - vector interpolation/distribution to d-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine vgtod2d_(ob,u_GS,v_GS,u_DL,v_DL, comm)
      use m_daInterp,only : daInterp_vatod
      use m_ppInterp,only : WIND
      use m_SubdomainDistributor,only : get
      use m_interleavedObject,only : localsize,totalsize
      use m_mpout,only : mpout_log
      use m_die,only : assert_
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:  ),intent(in) :: u_GS
      real,dimension(:,:  ),intent(in) :: v_GS
      real,dimension(:,:,:),intent(out) :: u_DL
      real,dimension(:,:,:),intent(out) :: v_DL
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vgtod2d_'
  integer :: iAdim,jAdim,iAlen,jAlen
  integer :: iGdim,jGdim,iGlen,jGlen
  integer :: msigLev1
  integer :: metaLev1

  real,allocatable,dimension(:,:,:) :: u_AL,v_AL
  integer :: ier
!________________________________________
! Get all dimension parameters
_ALLENTRY_

	call get(ob%aGrid,itotalSize=iAdim,jtotalSize=jAdim)
	metaLev1=localsize(ob%etaLev1)
	call get(ob%gGrid,ilocalSize=iGlen,jlocalSize=jGlen)

  	ASSERT(   iGlen==size(u_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(u_GS,1))
  	ASSERT(   iGlen==size(v_GS,2))		! swapped i-j
  	ASSERT(   jGlen==size(v_GS,1))

  	ASSERT(   iAdim==size(u_DL,1))
  	ASSERT(   jAdim==size(u_DL,2))
  	ASSERT(metaLev1==size(u_DL,3))
  	ASSERT(   iAdim==size(v_DL,1))
  	ASSERT(   jAdim==size(v_DL,2))
  	ASSERT(metaLev1==size(v_DL,3))

  	allocate(u_AL(iAdim,jAdim,metaLev1))
  	allocate(v_AL(iAdim,jAdim,metaLev1))

  call gtoa2d_(ob,u_GS,u_AL, WIND,comm)
  call gtoa2d_(ob,v_GS,v_AL, WIND,comm)

  call daInterp_vatod(u_AL,v_AL,u_DL,v_DL)
  	deallocate(u_AL)
  	deallocate(v_AL)
_ALLEXIT_
end subroutine vgtod2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sgtoa3d_ - scalar interpolation/distribution to a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sgtoa3d_(ob,ppIntp,w_SGS,w_EAL, key,comm, SrcTemp)
      use m_ppInterp,only : ppInterp
      use m_ppInterp,only : MASS,THTA,WIND
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none
      type(fgInterp),intent(in) :: ob
      type(ppInterp),intent(in) :: ppIntp
      real,dimension(:,:,:),intent(in ) :: w_SGS
      real,dimension(:,:,:),intent(out) :: w_EAL
      integer,intent(in) :: key
      integer,intent(in) :: comm
      integer,optional,intent(in) :: SrcTemp	! type of the source th

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sgtoa3d_'
  integer :: ier
_ALLENTRY_

  select case(key)
  case(MASS)
    call gtoa3d_(ob,ppIntp,w_SGS,w_EAL, key,comm, norder=1)
  case(THTA)
    call gtoa3d_(ob,ppIntp,w_SGS,w_EAL, key,comm,SrcTemp=SrcTemp)
  case(WIND)
    call die(myname_,'try vdtog3d_() for WIND')
  case default
    call die(myname_,'unknown key value',key)
  end select
_ALLEXIT_
end subroutine sgtoa3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sgtoa2d_ - scalar interpolation/distribution to a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sgtoa2d_(ob,w_GS,w_AL, key,comm)
      use m_ppInterp,only : MASS,THTA,WIND
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none
      type(fgInterp),intent(in) :: ob
      real,dimension(:,:  ),intent(in ) :: w_GS
      real,dimension(:,:,:),intent(out) :: w_AL
      integer,intent(in) :: key
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sgtoa2d_'
  integer :: ier
_ALLENTRY_

  select case(key)
  case(MASS)
    call gtoa2d_(ob,w_GS,w_AL, key,comm,norder=1)
  case(THTA)
    call gtoa2d_(ob,w_GS,w_AL, key,comm)
  case(WIND)
    call die(myname_,'try vdtog2d_() for WIND')
  case default
    call die(myname_,'unknown key value',key)
  end select
_ALLEXIT_
end subroutine sgtoa2d_
end module m_fgInterp
