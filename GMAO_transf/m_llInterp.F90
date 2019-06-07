!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_llInterp - lat-lon interpolations through interp_h()
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_llInterp
      implicit none
      private	! except

      public :: llInterp		! data structure

      public :: llInterp_g2ginit	! Gaussian to Gaussian
      public :: llInterp_l2ginit	! Lat-long to Gaussian
      public :: llInterp_g2linit	! Gaussian to lat-long
      public :: llInterp_l2linit	! Lat-long to lat-long

      public :: llInterp_clean,clean

      public :: llInterp_lh2rh
      public :: llInterp_rh2lh

      public :: llInterp_init ,init	! replaced by _ltoginit()
      public :: llInterp_atog		! replaced by _lh2rh()
!TORM      public :: llInterp_gtoa	! replaced by _rh2lh()

    type llInterp
      private
      		! a-grid, longitude by latitude, including poles.
      			! ticks of the two separable grid axses
      real,pointer,dimension(:) :: adlam	! d-lon in rad.
      real,pointer,dimension(:) :: adphi	! d-lat in rad.
      			! grid values of all grid points.
      real,pointer,dimension(:) :: alons	! lon. in rad.
      real,pointer,dimension(:) :: alats	! lat. in rad.

      		! Gaussian-grid, longitude by latitude, including poles.
      			! ticks of the two separable grid axses
      real,pointer,dimension(:) :: gdlam	! d-lon in rad.
      real,pointer,dimension(:) :: gdphi	! d-lat in rad.
      			! grid values of all grid points.
      real,pointer,dimension(:) :: glons	! lon. in rad.
      real,pointer,dimension(:) :: glats	! lat. in rad.
    end type llInterp

    interface llInterp_g2ginit; module procedure	&
    	g2ginit_; end interface
    interface llInterp_l2ginit; module procedure	&
    	l2ginit_; end interface
    interface llInterp_g2linit; module procedure	&
    	g2linit_; end interface
    interface llInterp_l2linit; module procedure	&
    	l2linit_; end interface

    	! These two are here for backward compatibility.
    interface llInterp_init; module procedure	&
    	l2ginit_; end interface
    interface  init; module procedure l2ginit_; end interface

    interface llInterp_clean; module procedure	&
    	clean_; end interface
    interface clean; module procedure clean_; end interface

    interface llInterp_lh2rh; module procedure	&
	atog2dr_,	&
    	atog3dr_; end interface

    interface llInterp_rh2lh; module procedure	&
	gtoa2dr_,	&
    	gtoa3dr_; end interface

    	! Remain here for backward compatibility
    interface llInterp_atog; module procedure	&
	atog2dr_,	&
    	atog3dr_; end interface

    interface llInterp_gtoa; module procedure	&
	gtoa2dr_,	&
    	gtoa3dr_; end interface

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_interpH'

  real,parameter :: NOROT = 0.
  real,parameter :: NOTLT = 90.
  real,parameter :: NOPRE = 0.

  integer,parameter :: DEFAULT_ORDER=3

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: g2ginit_ - initialize an object for a Gaussian-to-Gaussian
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine g2ginit_(ob,im_a,jm_a,im_g,jm_g)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! Dimensions of Gaussian grid A
      integer,intent(in) :: im_g,jm_g	! Dimensions of Guassian grid B

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::g2ginit_'

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(jm_g>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLG_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

	call setLG_(im_g,jm_g,dlam=ob%gdlam(1:im_g-1   ),	&
			      dphi=ob%gdphi(1:jm_g-1   ),	&
			      lons=ob%glons(1:im_g*jm_g),	&
			      lats=ob%glats(1:im_g*jm_g)	)

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine g2ginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2ginit_ - initialize an object for a lat-long to Gaussian
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2ginit_(ob,im_a,jm_a,im_g,jm_g)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! a-grid dimensions
      integer,intent(in) :: im_g,jm_g	! Gaussian-grid dimensions

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2ginit_'

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(jm_g>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLL_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

	call setLG_(im_g,jm_g,dlam=ob%gdlam(1:im_g-1   ),	&
			      dphi=ob%gdphi(1:jm_g-1   ),	&
			      lons=ob%glons(1:im_g*jm_g),	&
			      lats=ob%glats(1:im_g*jm_g)	)

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine l2ginit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: g2linit_ - initialize an object for Gaussian to lat-long
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine g2linit_(ob,im_a,jm_a,im_g,jm_g)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! a-grid dimensions
      integer,intent(in) :: im_g,jm_g	! Gaussian-grid dimensions

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::g2linit_'

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(jm_g>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLG_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

	call setLL_(im_g,jm_g,dlam=ob%gdlam(1:im_g-1   ),	&
			      dphi=ob%gdphi(1:jm_g-1   ),	&
			      lons=ob%glons(1:im_g*jm_g),	&
			      lats=ob%glats(1:im_g*jm_g)	)

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine g2linit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2linit_ - initialize an object for a lat-long to lat-long
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2linit_(ob,im_a,jm_a,im_g,jm_g)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! a-grid dimensions
      integer,intent(in) :: im_g,jm_g	! Gaussian-grid dimensions

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2linit_'

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(jm_g>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLL_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

	call setLL_(im_g,jm_g,dlam=ob%gdlam(1:im_g-1   ),	&
			      dphi=ob%gdphi(1:jm_g-1   ),	&
			      lons=ob%glons(1:im_g*jm_g),	&
			      lats=ob%glats(1:im_g*jm_g)	)

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine l2linit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: g2ainit_ - initialize an object for Gaussian to any
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine g2ainit_(ob,im_a,jm_a,im_g,rlats,inrad)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! Longitude/Gaussian dimensions
      integer,intent(in) :: im_g	! other grid dimension in long.
      real,dimension(:),intent(in) :: rlats	! lat grid in deg.
      logical,optional,intent(in) :: inrad	! rlats is in radiance

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::g2ainit_'
  real :: onepi,deg
  integer :: jm_g
  logical :: inrad_
  inrad_=.false.
  if(present(inrad)) inrad_=inrad

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(size(rlats)>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	jm_g=size(rlats)
	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLG_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

      if(.not.inrad_) then
		! If values of rlats(:) are in degree, they are
		! consistent with setLA_() assumption
	call setLA_(im_g,rlats,dlam=ob%gdlam(1:im_g-1   ),	&
			       dphi=ob%gdphi(1:jm_g-1   ),	&
			       lons=ob%glons(1:im_g*jm_g),	&
			       lats=ob%glats(1:im_g*jm_g)	)
      else
		! If values of rlats(:) are in radiance, they need to
		! be converted to degrees.

	onepi=4.*atan(1.); deg=180./onepi
	call setLA_(im_g,rlats*deg,				&
			       dlam=ob%gdlam(1:im_g-1   ),	&
			       dphi=ob%gdphi(1:jm_g-1   ),	&
			       lons=ob%glons(1:im_g*jm_g),	&
			       lats=ob%glats(1:im_g*jm_g)	)
      endif

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine g2ainit_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2ainit_ - initialize an object for Lat-long to Any-long
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2ainit_(ob,im_a,jm_a,im_g,rlats,inrad)
      use m_die,only : assert_
      implicit none
      type(llInterp),intent(out) :: ob
      integer,intent(in) :: im_a,jm_a	! a-grid dimensions
      integer,intent(in) :: im_g	! other grid dimension in long.
      real,dimension(:),intent(in) :: rlats	! lat grid in deg.
      logical,optional,intent(in) :: inrad	! rlats is in radiance

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2ainit_'
  real :: onepi,deg
  integer :: jm_g
  logical :: inrad_
  inrad_=.false.
  if(present(inrad)) inrad_=inrad

  	ASSERT(im_a>1)
  	ASSERT(jm_a>2)
  	ASSERT(im_g>1)
  	ASSERT(size(rlats)>2)

	allocate(ob%adlam(im_a))
	allocate(ob%adphi(jm_a))
	allocate(ob%alons(im_a*jm_a))
	allocate(ob%alats(im_a*jm_a))

	jm_g=size(rlats)
	allocate(ob%gdlam(im_g))
	allocate(ob%gdphi(jm_g))
	allocate(ob%glons(im_g*jm_g))
	allocate(ob%glats(im_g*jm_g))

	call setLL_(im_a,jm_a,dlam=ob%adlam(1:im_a-1   ),	&
			      dphi=ob%adphi(1:jm_a-1   ),	&
			      lons=ob%alons(1:im_a*jm_a),	&
			      lats=ob%alats(1:im_a*jm_a)	)

      if(.not.inrad_) then
		! If values of rlats(:) are in degree, they are
		! consistent with setLA_() assumption
	call setLA_(im_g,rlats,dlam=ob%gdlam(1:im_g-1   ),	&
			       dphi=ob%gdphi(1:jm_g-1   ),	&
			       lons=ob%glons(1:im_g*jm_g),	&
			       lats=ob%glats(1:im_g*jm_g)	)
      else
		! If values of rlats(:) are in radiance, they need to
		! be converted to degrees.

	onepi=4.*atan(1.); deg=180./onepi
	call setLA_(im_g,rlats*deg,				&
			       dlam=ob%gdlam(1:im_g-1   ),	&
			       dphi=ob%gdphi(1:jm_g-1   ),	&
			       lons=ob%glons(1:im_g*jm_g),	&
			       lats=ob%glats(1:im_g*jm_g)	)
      endif

  	! Although last elements in both dlam(:) and dphi(:) are not
	! used, they are passed to interp_h() anyway, since the
	! interface expects that.

end subroutine l2ainit_

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
      implicit none
      type(llInterp),intent(inout) :: ob

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  deallocate(ob%adlam)
  deallocate(ob%adphi)
  deallocate(ob%alons)
  deallocate(ob%alats)
  deallocate(ob%gdlam)
  deallocate(ob%gdphi)
  deallocate(ob%glons)
  deallocate(ob%glats)

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: atog3dr_ - interpolate a 3-d field on a-grid to gaussian
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine atog3dr_(ob,wa,wg, norder,vector,nocheck,undef)
      use m_die,only : assert_
      use m_undef,only : undef_ssi
      implicit none
      type(llInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: wa	! a field on a-grid
      real,dimension(:,:,:),intent(out) :: wg	! a field on g-grid

      integer,optional,intent(in) :: norder	! order of interp.
      logical,optional,intent(in) :: vector	! is w a vector
      logical,optional,intent(in) :: nocheck	! for undef in wa
      real   ,optional,intent(in) :: undef	! "undefined" value

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::atog3dr_'
  integer :: im_a,jm_a
  integer :: im_g,jm_g
  integer :: lm
  logical :: chk_
  logical :: vec_
  integer :: sgn_
  integer :: ord_
  real    :: undef_
  
  	! the field type (although I don't see how it works)
  vec_=.false.
  if(present(vector)) vec_=vector
  sgn_=1
  if(vec_) sgn_=-1

  	! order of interpolation scheme
  ord_=DEFAULT_ORDER
  if(present(norder)) ord_=norder

  	! checking for undefined values.

  undef_=undef_ssi()
  chk_=.true.
  if(present(nocheck)) chk_=.not.nocheck
  if(chk_.and.present(undef)) undef_=undef

  	! dimension checking
  im_a=size(ob%adlam)
  jm_a=size(ob%adphi)
  	ASSERT(im_a==size(wa,1))
  	ASSERT(jm_a==size(wa,2))
	ASSERT(im_a*jm_a==size(ob%alons))
	ASSERT(im_a*jm_a==size(ob%alats))
  im_g=size(ob%gdlam)
  jm_g=size(ob%gdphi)
  	ASSERT(im_g==size(wg,1))
  	ASSERT(jm_g==size(wg,2))
	ASSERT(im_g*jm_g==size(ob%glons))
	ASSERT(im_g*jm_g==size(ob%glats))

  lm=size(wa,3)
  	ASSERT(lm==size(wg,3))

  	! interpolation
  call interp_h( wa,im_a,jm_a,lm, ob%adlam,ob%adphi,		&
  					NOROT,NOTLT,NOPRE,	&
  		 wg,im_g*jm_g,    ob%glons,ob%glats,		&
		 			sgn_,ord_,chk_,undef_)
end subroutine atog3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gtoa3dr_ - interpolate a 3-d field on gaussian to a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gtoa3dr_(ob,wg,wa,norder,vector,nocheck,undef)
      use m_die,only : assert_
      use m_undef,only : undef_ssi
      implicit none
      type(llInterp),intent(in) :: ob
      real,dimension(:,:,:),intent(in) :: wg	! a field on Gaussian
      real,dimension(:,:,:),intent(out) :: wa	! a field on a-grid

      integer,optional,intent(in) :: norder	! order of interp.
      logical,optional,intent(in) :: vector	! is w a vector
      logical,optional,intent(in) :: nocheck	! for undef in wa
      real   ,optional,intent(in) :: undef	! value of "undefined"

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gtoa3dr_'
  integer :: im_a,jm_a
  integer :: im_g,jm_g
  integer :: lm
  logical :: chk_
  logical :: vec_
  integer :: sgn_
  integer :: ord_
  real    :: undef_
  
  	! the field type (although I don't see how it works)
  vec_=.false.
  if(present(vector)) vec_=vector
  sgn_=1
  if(vec_) sgn_=-1

  	! order of interpolation scheme
  ord_=DEFAULT_ORDER
  if(present(norder)) ord_=norder

  	! checking for undefined values.

  undef_=undef_ssi()
  chk_=.true.
  if(present(nocheck)) chk_=.not.nocheck
  if(chk_.and.present(undef)) undef_=undef

  	! dimension checking
  im_a=size(ob%adlam)
  jm_a=size(ob%adphi)
  	ASSERT(im_a==size(wa,1))
  	ASSERT(jm_a==size(wa,2))
	ASSERT(im_a*jm_a==size(ob%alons))
	ASSERT(im_a*jm_a==size(ob%alats))
  im_g=size(ob%gdlam)
  jm_g=size(ob%gdphi)
  	ASSERT(im_g==size(wg,1))
  	ASSERT(jm_g==size(wg,2))
	ASSERT(im_g*jm_g==size(ob%glons))
	ASSERT(im_g*jm_g==size(ob%glats))

  lm=size(wa,3)
  	ASSERT(lm==size(wg,3))

  	! interpolation
  call interp_h( wg,im_g,jm_g,lm, ob%gdlam,ob%gdphi,		&
  					NOROT,NOTLT,NOPRE,	&
  		 wa,im_a*jm_a,    ob%alons,ob%alats,		&
		 			sgn_,ord_,chk_,undef_)
end subroutine gtoa3dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: atog2dr_ - interpolate a 2-d field on a-grid to gaussian
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine atog2dr_(ob,wa,wg, norder,vector,nocheck,undef)
      use m_die,only : assert_
      use m_undef,only : undef_ssi
      implicit none
      type(llInterp),intent(in) :: ob
      real,dimension(:,:),intent(in) :: wa	! a field on a-grid
      real,dimension(:,:),intent(out) :: wg	! a field on g-grid

      integer,optional,intent(in) :: norder	! order of interp.
      logical,optional,intent(in) :: vector	! is w a vector
      logical,optional,intent(in) :: nocheck	! for undef in wa
      real   ,optional,intent(in) :: undef	! "undefined" value

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::atog2dr_'
  integer :: im_a,jm_a
  integer :: im_g,jm_g
  logical :: chk_
  logical :: vec_
  integer :: sgn_
  integer :: ord_
  real    :: undef_
  
  	! the field type (although I don't see how it works)
  vec_=.false.
  if(present(vector)) vec_=vector
  sgn_=1
  if(vec_) sgn_=-1

  	! order of interpolation scheme
  ord_=DEFAULT_ORDER
  if(present(norder)) ord_=norder

  	! checking for undefined values.

  undef_=undef_ssi()
  chk_=.true.
  if(present(nocheck)) chk_=.not.nocheck
  if(chk_.and.present(undef)) undef_=undef

  	! dimension checking
  im_a=size(ob%adlam)
  jm_a=size(ob%adphi)
  	ASSERT(im_a==size(wa,1))
  	ASSERT(jm_a==size(wa,2))
	ASSERT(im_a*jm_a==size(ob%alons))
	ASSERT(im_a*jm_a==size(ob%alats))
  im_g=size(ob%gdlam)
  jm_g=size(ob%gdphi)
  	ASSERT(im_g==size(wg,1))
  	ASSERT(jm_g==size(wg,2))
	ASSERT(im_g*jm_g==size(ob%glons))
	ASSERT(im_g*jm_g==size(ob%glats))

  	! interpolation
  call interp_h( wa,im_a,jm_a, 1, ob%adlam,ob%adphi,		&
  					NOROT,NOTLT,NOPRE,	&
  		 wg,im_g*jm_g,    ob%glons,ob%glats,		&
		 			sgn_,ord_,chk_,undef_)
end subroutine atog2dr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gtoa2dr_ - interpolate a 2-d field on gaussian to a-grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gtoa2dr_(ob,wg,wa,norder,vector,nocheck,undef)
      use m_die,only : assert_
      use m_undef,only : undef_ssi
      implicit none
      type(llInterp),intent(in) :: ob
      real,dimension(:,:),intent(in) :: wg	! a field on Gaussian
      real,dimension(:,:),intent(out) :: wa	! a field on a-grid

      integer,optional,intent(in) :: norder	! order of interp.
      logical,optional,intent(in) :: vector	! is w a vector
      logical,optional,intent(in) :: nocheck	! for undef in wa
      real   ,optional,intent(in) :: undef	! value of "undefined"

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gtoa2dr_'
  integer :: im_a,jm_a
  integer :: im_g,jm_g
  logical :: chk_
  logical :: vec_
  integer :: sgn_
  integer :: ord_
  real    :: undef_
  
  	! the field type (although I don't see how it works)
  vec_=.false.
  if(present(vector)) vec_=vector
  sgn_=1
  if(vec_) sgn_=-1

  	! order of interpolation scheme
  ord_=DEFAULT_ORDER
  if(present(norder)) ord_=norder

  	! checking for undefined values.

  undef_=undef_ssi()
  chk_=.true.
  if(present(nocheck)) chk_=.not.nocheck
  if(chk_.and.present(undef)) undef_=undef

  	! dimension checking
  im_a=size(ob%adlam)
  jm_a=size(ob%adphi)
  	ASSERT(im_a==size(wa,1))
  	ASSERT(jm_a==size(wa,2))
	ASSERT(im_a*jm_a==size(ob%alons))
	ASSERT(im_a*jm_a==size(ob%alats))
  im_g=size(ob%gdlam)
  jm_g=size(ob%gdphi)
  	ASSERT(im_g==size(wg,1))
  	ASSERT(jm_g==size(wg,2))
	ASSERT(im_g*jm_g==size(ob%glons))
	ASSERT(im_g*jm_g==size(ob%glats))

  	! interpolation
  call interp_h( wg,im_g,jm_g, 1, ob%gdlam,ob%gdphi,		&
  					NOROT,NOTLT,NOPRE,	&
  		 wa,im_a*jm_a,    ob%alons,ob%alats,		&
		 			sgn_,ord_,chk_,undef_)
end subroutine gtoa2dr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setLL_ - set a-grid lat-lon grid data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setLL_(im,jm,dlam,dphi,lons,lats)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: im,jm
      real,dimension(:),optional,intent(out) :: dlam,dphi
      real,dimension(:),optional,intent(out) :: lons,lats

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- aglorithm adapted from set_regrid()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setLL_'
  real :: halfpi,onepi,twopi
  real :: dlon,dlat
  integer :: i,j,ij

  halfpi= 2.*atan(1.)
  onepi = 4.*atan(1.)
  twopi = 8.*atan(1.)

  if(present(dlam)) then
  	ASSERT(im-1==size(dlam))

    dlam(:) = twopi/im
  endif

  if(present(dphi)) then
  	ASSERT(jm-1==size(dphi))

    dphi(:) = onepi/(jm-1)
  endif

  if(present(lons)) then
  	ASSERT(im*jm==size(lons))

    ij=0
    dlon=twopi/im
    do j=1,jm
      do i=1,im
        ij=ij+1
        lons(ij) = -onepi + (i-1)*dlon
      enddo
    enddo
  endif

  if(present(lats)) then
  	ASSERT(im*jm==size(lats))
    ij=0
    dlat=onepi/(jm-1)
    do j=1,jm
      lats(ij+1:ij+im) = (j-1)*dlat - halfpi
      ij=ij+im
   enddo
  endif

end subroutine setLL_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setLA_ - set a-grid lat-lon grid data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setLA_(im,rlats,dlam,dphi,lons,lats)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: im
      real,dimension(:),intent(in) :: rlats	! in degree-angle
      real,dimension(:),optional,intent(out) :: dlam,dphi
      real,dimension(:),optional,intent(out) :: lons,lats

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- aglorithm adapted from set_regrid()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setLA_'
  real :: halfpi,onepi,twopi
  real :: dlon,dlat
  integer :: i,j,ij
  integer :: jm
  real :: RAD

  halfpi= 2.*atan(1.)
  onepi = 4.*atan(1.)
  twopi = 8.*atan(1.)
  RAD   = onepi/180.

  if(present(dlam)) then
  	ASSERT(im-1==size(dlam))

    dlam(:) = twopi/im
  endif

  jm=size(rlats)
  if(present(dphi)) then
  	ASSERT(jm-1==size(dphi))

    dphi(:) = (rlats(2:jm)-rlats(1:jm-1))*RAD
  endif

  if(present(lons)) then
  	ASSERT(im*jm==size(lons))

    ij=0
    dlon=twopi/im
    do j=1,jm
      do i=1,im
        ij=ij+1
        lons(ij) = -onepi + (i-1)*dlon
      enddo
    enddo
  endif

  if(present(lats)) then
  	ASSERT(im*jm==size(lats))
    ij=0
    do j=1,jm
      lats(ij+1:ij+im) = rlats(j)*RAD
      ij=ij+im
   enddo
  endif

end subroutine setLA_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setLG_ - set Gaussian-grid lat-lon grid data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setLG_(im,jm,dlam,dphi,lons,lats)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: im,jm	! numbers of grid points
      real,dimension(:),optional,	&
      	intent(out) :: dlam,dphi	! grid intervals, im-1|jm-1
      real,dimension(:),optional,	&
      	intent(out) :: lons,lats	! grid values, im*jm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- aglorithm adapted from set_regrid()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setLG_'
  real :: halfpi,onepi,twopi,rad
  real :: dlon,dlat
  real,dimension(jm-2) :: glats
  integer :: i,j,ij

  	ASSERT(im>1)
  	ASSERT(jm>2)

  halfpi= 2.*atan(1.)
  onepi = 4.*atan(1.)
  twopi = 8.*atan(1.)
  rad   = onepi/180.

  	! Get Gaussian latitude grid points in degrees without poles.

  call gauss_lat_nmc(glats(1:jm-2),jm-2)

  if(present(dlam)) then
  	ASSERT(im-1==size(dlam))

    dlam(1:im-1) = twopi/im
  endif

  if(present(dphi)) then
  	ASSERT(jm-1==size(dphi))

    	! Compute latitude grid intervals in radiance with poles.

    dphi(1) = glats(1)*rad + halfpi
    do j=2,jm-2
      dphi(j) = (glats(j)-glats(j-1))*rad
    end do
    dphi(jm-1) = halfpi - glats(jm-2)*rad
  endif

  if(present(lons)) then
  	ASSERT(im*jm==size(lons))

	! Compute longitude grid points in radiance with poles

    ij=0
    dlon=twopi/im
    do j=1,jm
      do i=1,im
        ij=ij+1
        lons(ij) = -onepi + (i-1)*dlon
      enddo
    enddo
  endif

  if(present(lats)) then
  	ASSERT(im*jm==size(lats))

	! Compute latitude grid points in radiance with poles
    ij=0
    lats(ij+1:ij+im) = -halfpi
    ij=ij+im
    do j=2,jm-1
      lats(ij+1:ij+im) = glats(j-1)*rad
      ij=ij+im
    enddo
    lats(ij+1:ij+im) = halfpi
  endif

end subroutine setLG_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setXX_ - set lat-lon grid data with a given latitude grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setXX_(im,rlats,dlam,dphi,lons,lats)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: im		! numbers of grid points
      real,dimension(:),intent(in) :: rlats
      real,dimension(:),optional,	&
      	intent(out) :: dlam,dphi	! grid intervals, im-1|jm-1
      real,dimension(:),optional,	&
      	intent(out) :: lons,lats	! grid values, im*jm

! !REVISION HISTORY:
! 	27Oct04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- aglorithm adapted from set_regrid()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setXX_'
  real :: halfpi,onepi,twopi,rad
  real :: dlon,dlat
  integer :: i,j,ij,jm

  	ASSERT(im>1)
  	ASSERT(size(rlats)>2)
  jm=size(rlats)

  halfpi= 2.*atan(1.)
  onepi = 4.*atan(1.)
  twopi = 8.*atan(1.)
  rad   = onepi/180.

  	! Get latitude grid points in degrees without poles.

  if(present(dlam)) then
  	ASSERT(im-1==size(dlam))

    dlam(1:im-1) = twopi/im
  endif

  if(present(dphi)) then
  	ASSERT(jm-1==size(dphi))

    	! Compute latitude grid intervals in radiance with poles.

    dphi(1) = rlats(2)*rad + halfpi
    do j=2,jm-2
      dphi(j) = (rlats(j)-rlats(j-1))*rad
    end do
    dphi(jm-1) = halfpi - rlats(jm-2)*rad
  endif

  if(present(lons)) then
  	ASSERT(im*jm==size(lons))

	! Compute longitude grid points in radiance with poles

    ij=0
    dlon=twopi/im
    do j=1,jm
      do i=1,im
        ij=ij+1
        lons(ij) = -onepi + (i-1)*dlon
      enddo
    enddo
  endif

  if(present(lats)) then
  	ASSERT(im*jm==size(lats))

	! Compute latitude grid points in radiance with poles
    ij=0
    lats(ij+1:ij+im) = -halfpi
    ij=ij+im
    do j=2,jm-1
      lats(ij+1:ij+im) = rlats(j-1)*rad
      ij=ij+im
    enddo
    lats(ij+1:ij+im) = halfpi
  endif

end subroutine setXX_
end module m_llInterp
