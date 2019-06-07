!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_undef - reset UNDEF values from one system to another
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_undef
      implicit none
      private	! except

      public :: undef_ssi		! NCEP/SSI UNDEF value
      public :: undef_fvgcm		! GMAO/fvGCM UNDEF value
      public :: undef_geos3		! DAO/GEOS3 UNDEF value

      public :: undef_fv2ssi
      public :: undef_2ssi
      public :: unlike,like

    interface undef_fv2ssi; module procedure	&
    	fv2ssi1d_,	&
	fv2ssi2d_; end interface
    interface undef_2ssi; module procedure	&
    	udf2ssi1d_,	&
	udf2ssi2d_; end interface
    interface unlike; module procedure	&
    	unlike1ds_,	&
    	unlike2ds_,	&
	unlike3ds_; end interface
    interface like; module procedure	&
    	like1ds_,	&
    	like2ds_,	&
	like3ds_; end interface

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_undef'

  real,parameter :: undef_ssi_   = -9.99e33		! SSI UNDEF
  real,parameter :: undef_fvgcm_ = +1.e25		! fvGCM UNDEF
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: undef_ssi - NCEP/SSI UNDEF
!
! !DESCRIPTION:
!
! !INTERFACE:

    function undef_ssi()
      implicit none
      real :: undef_ssi

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*),parameter :: myname_=myname//'::undef_ssi'
  undef_ssi=undef_ssi_
end function undef_ssi

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: undef_fvgcm - GMAO/fvGCM UNDEF
!
! !DESCRIPTION:
!
! !INTERFACE:

    function undef_fvgcm()
      implicit none
      real :: undef_fvgcm

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*),parameter :: myname_=myname//'::undef_fvgcm'
  undef_fvgcm=undef_fvgcm_
end function undef_fvgcm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: undef_geos3 - GMAO/fvGCM UNDEF
!
! !DESCRIPTION:
!
! !INTERFACE:

    function undef_geos3()
      use m_geosapi,only : getcon
      implicit none
      real :: undef_geos3

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*),parameter :: myname_=myname//'::undef_geos3'
  undef_geos3=getcon('UNDEF')
end function undef_geos3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fv2ssi1d_ - Replace undef FV-GCM values by SSI undef
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine fv2ssi1d_ (var, where, reltol, verb, vname )
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none

      real,dimension(:),intent(inout) ::     var
      character(len=*) ,intent(in) :: where
      real,optional    ,intent(in) :: reltol
      logical,optional ,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*), parameter :: myname_ = myname//'::fv2ssi1d_'

  call udf2udf1d_(var, undef_fvgcm(),undef_ssi(),where,	&
			reltol=reltol, verb=verb, vname=vname )

end subroutine fv2ssi1d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fv2ssi2d_ - Replace undef FV-GCM values by SSI undef
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine fv2ssi2d_ (var, where, reltol, verb, vname )
      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none

      real,dimension(:,:),intent(inout) ::     var
      character(len=*) ,intent(in) :: where
      real,optional    ,intent(in) :: reltol
      logical,optional ,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*), parameter :: myname_ = myname//'::fv2ssi2d_'

  call udf2udf2d_(var, undef_fvgcm(),undef_ssi(),where,	&
			reltol=reltol, verb=verb, vname=vname )

end subroutine fv2ssi2d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: udf2ssi1d_ - Substitute a given undef value by SSI undef
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine udf2ssi1d_ (var, undef_in, where, reltol, verb, vname )
      use m_mpout,only : mpout_ison,mpout
      use m_die,only : die
      implicit none

      real,dimension(:),intent(inout) ::     var
      real	       ,intent(in) :: undef_in
      character(len=*) ,intent(in) :: where
      real,optional    ,intent(in) :: reltol
      logical,optional ,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*), parameter :: myname_ = myname//'::udf2ssi1d_'
  real :: reltol_

  reltol_=.01
  if(present(reltol)) reltol_=reltol

  if(mpout_ison()) then
    call lookfor_(undef_ssi()  ,reltol_,var,mpout,where,vname)	! -.999e33
    call lookfor_(undef_fvgcm(),reltol_,var,mpout,where,vname)	! 1.e25
    call lookfor_(undef_geos3(),reltol_,var,mpout,where,vname)	! 1.e15
    call lookfor_(undef_in     ,reltol_,var,mpout,where,vname)	! given
    call lookfor_(1.e12        ,reltol_,var,mpout,where,vname)	! 1.e12
  endif

  call udf2udf1d_(var, undef_in,undef_ssi(),where,	&
			reltol=reltol, verb=verb, vname=vname )
contains
subroutine lookfor_(spval,reltol,v,lu,where,vname)
  implicit none
  real,intent(in) :: spval
  real,intent(in) :: reltol
  real,dimension(:),intent(in) :: v
  integer,intent(in) :: lu
  character(len=*),intent(in) :: where,vname

  integer :: n,m
  n=count(abs(v-spval)<=abs(reltol*spval))
  m=size(v)
  if(n/=0) write(lu,'(2a,e12.6,3a,i6,a,i6)') trim(where),	&
	': find ',spval,' in variable "',trim(vname),'", ',	&
	n,' out of ',m
end subroutine lookfor_
end subroutine udf2ssi1d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: udf2ssi2d_ - Substitute a given undef value by SSI undef
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine udf2ssi2d_ (var, undef_in, where, reltol, verb, vname )
      use m_mpout,only : mpout_ison,mpout
      use m_die,only : die
      implicit none

      real,dimension(:,:),intent(inout) ::     var
      real             ,intent(in) :: undef_in
      character(len=*) ,intent(in) :: where
      real,optional    ,intent(in) :: reltol
      logical,optional ,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !REVISION HISTORY:
! 	26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________
  character(len=*), parameter :: myname_ = myname//'::udf2ssi2d_'
  real :: reltol_

  reltol_=.01
  if(present(reltol)) reltol_=reltol

  if(mpout_ison()) then
    call lookfor_(undef_ssi()  ,reltol_,var,mpout,where,vname)	! -.999e33
    call lookfor_(undef_fvgcm(),reltol_,var,mpout,where,vname)	! 1.e25
    call lookfor_(undef_geos3(),reltol_,var,mpout,where,vname)	! 1.e15
    call lookfor_(undef_in     ,reltol_,var,mpout,where,vname)	! given
    call lookfor_(1.e12        ,reltol_,var,mpout,where,vname)	! 1.e12
  endif

  call udf2udf2d_(var, undef_in,undef_ssi(),where,	&
			reltol=reltol, verb=verb, vname=vname )

contains
subroutine lookfor_(spval,reltol,v,lu,where,vname)
  implicit none
  real,intent(in) :: spval
  real,intent(in) :: reltol
  real,dimension(:,:),intent(in) :: v
  integer,intent(in) :: lu
  character(len=*),intent(in) :: where,vname

  integer :: n,m
  n=count(abs(v-spval)<=abs(reltol*spval))
  m=size(v)
  if(n/=0) write(lu,'(2a,e12.6,3a,i6,a,i6)') trim(where),	&
	': find ',spval,' in variable "',trim(vname),'", ',	&
	n,' out of ',m
end subroutine lookfor_
end subroutine udf2ssi2d_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  udf2udf1d_ - Replace an UNDEF in the input to another
! 
! !INTERFACE:
!
      subroutine udf2udf1d_(var,undef_in,undef_out,where,	&
      	reltol,verb,vname)
!
! !USES:


      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none

! !INPUT PARAMETERS: 

	real,intent(in) :: undef_in
	real,intent(in) :: undef_out
	character(len=*),intent(in) :: where

! !INPUT/OUTPUT PARAMETERS: 
 
      real,dimension(:),intent(inout) ::     var

! ! OPTIONAL PARAMETERS:

      real,optional,intent(in) :: reltol
      logical,optional,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  13Apr00   Todling    - Initial code.
!  26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- modified/extended for multi-dimensional interfaces
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::udf2udf1d_'

      integer  i
      integer  nonmiss
      real     tol
      real     val_max_in, val_max_out
      logical :: verbose

      verbose=.false.
      if(present(verb)) verbose=verb

      tol   = 0.01
      if(present(reltol)) tol=reltol
      tol   = tol*undef_in

!     Make sure there is no precision problem with undef's
!     ----------------------------------------------------
      nonmiss = 0
      val_max_in  = maxval(var)
      do i = 1, size(var)
         if ( abs(var(i)-undef_in) .lt. tol ) then
              nonmiss = nonmiss + 1
              var(i) = undef_out
         end if
      end do
      val_max_out = maxval(abs(var))
      if (nonmiss.ne.0) then
         if(verbose) then
	   if(present(vname)) then
	     call mpout_log(where,	&
	     	'No. of UNDEF_in fixed to UNDEF_out for "'//	&
		trim(vname)//'"',nonmiss)
	   else
	     call mpout_log(where,	&
	     	'No. of UNDEF_in fixed to UNDEF_out',nonmiss)
	   endif
	 endif

         if ( val_max_out .ne. abs(undef_out) ) then
	   call mpout_log(where,	&
	   	' Largest  value on  input: ',  val_max_in)
           call mpout_log(where,	&
	   	' Largest  value on output: ',  val_max_out)
           call mpout_log(where,	&
		' Undef_in     value spec.: ',  undef_in)
           call mpout_log(where,	&
		' Undef_out    value spec.: ',  undef_out)
         end if
      end if

      return
      end subroutine udf2udf1d_
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  udf2udf2d_ - Replace undef FV values on input by SSI undef
! 
! !INTERFACE:
!
      subroutine udf2udf2d_ (var,undef_in,undef_out,where,	&
      	reltol,verb,vname)
!
! !USES:

      use m_mpout,only : mpout_log
      use m_die,only : die
      implicit none

! !INPUT/OUTPUT PARAMETERS: 

      real,intent(in) :: undef_in
      real,intent(in) :: undef_out
      character(len=*),intent(in) :: where

! !INPUT/OUTPUT PARAMETERS: 
 
      real,dimension(:,:),intent(inout) ::     var

! ! OPTIONAL PARAMETERS:

      real,optional,intent(in) :: reltol
      logical,optional,intent(in) :: verb
      character(len=*),optional,intent(in) :: vname

! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  13Apr00   Todling    - Initial code.
!  26Jul05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- modified/extended for multi-dimensional interfaces
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::udf2udf2d_'

      integer  i,j
      integer  nonmiss
      real     tol
      real     val_max_in, val_max_out
      logical :: verbose

      verbose=.false.
      if(present(verb)) verbose=verb

      tol   = 0.01
      if(present(reltol)) tol=reltol
      tol   = tol*undef_in

!     Make sure there is no precision problem with undef's
!     ----------------------------------------------------
      nonmiss = 0
      val_max_in  = maxval(var)
      do j = 1, size(var,2)
      do i = 1, size(var,1)
         if ( abs(var(i,j)-undef_in) .lt. tol ) then
              nonmiss = nonmiss + 1
              var(i,j) = undef_out
         end if
      end do
      end do
      val_max_out = maxval(abs(var))

      if (nonmiss.ne.0) then
         if(verbose) then
	   if(present(vname)) then
	     call mpout_log(where,	&
	     	'No. of UNDEF_in fixed to UNDEF_out for "'//	&
		trim(vname)//'"',nonmiss)
	   else
	     call mpout_log(where,	&
	     	'No. of UNDEF_in fixed to UNDEF_out',nonmiss)
	   endif
	 endif

         if ( val_max_out .ne. abs(undef_out) ) then
	   call mpout_log(where,	&
	   	' Largest  value on  input: ',  val_max_in)
           call mpout_log(where,	&
	   	' Largest  value on output: ',  val_max_out)
           call mpout_log(where,	&
		' Undef_in     value spec.: ',  undef_in)
           call mpout_log(where,	&
		' Undef_out    value spec.: ',  undef_out)
         end if
      end if

      return
      end subroutine udf2udf2d_

function unlike1ds_(a,v,reltol)
  implicit none
  real,dimension(:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a,1)) :: unlike1ds_

  unlike1ds_(:)=.not.like1ds_(a(:),v,reltol=reltol)
end function unlike1ds_
function unlike2ds_(a,v,reltol)
  implicit none
  real,dimension(:,:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a,1),size(a,2)) :: unlike2ds_

  unlike2ds_(:,:)=.not.like2ds_(a(:,:),v,reltol=reltol)
end function unlike2ds_
function unlike3ds_(a,v,reltol)
  implicit none
  real,dimension(:,:,:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a,1),size(a,2),size(a,3)) :: unlike3ds_

  unlike3ds_(:,:,:)=.not.like3ds_(a(:,:,:),v,reltol=reltol)
end function unlike3ds_
function like1ds_(a,v,reltol)
  implicit none
  real,dimension(:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a)) :: like1ds_

  real :: tol
  tol=.01
  if(present(reltol)) tol=reltol
  tol=abs(v*tol)
  like1ds_(:)=abs(a(:)-v)<=tol
end function like1ds_
function like2ds_(a,v,reltol)
  implicit none
  real,dimension(:,:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a,1),size(a,2)) :: like2ds_

  real :: tol
  tol=.01
  if(present(reltol)) tol=reltol
  tol=abs(v*tol)
  like2ds_(:,:)=abs(a(:,:)-v)<=tol
end function like2ds_
function like3ds_(a,v,reltol)
  implicit none
  real,dimension(:,:,:),intent(in) :: a
  real,intent(in) :: v
  real,optional,intent(in) :: reltol
  logical,dimension(size(a,1),size(a,2),size(a,3)) :: like3ds_

  real :: tol
  tol=.01
  if(present(reltol)) tol=reltol
  tol=abs(v*tol)
  like3ds_(:,:,:)=abs(a(:,:,:)-v)<=tol
end function like3ds_

end module m_undef
