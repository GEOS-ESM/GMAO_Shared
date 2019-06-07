!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_gg --- FVSSI gaussian grid state class
!
! !INTERFACE:
!
   module m_gg
!
!USES:
!
   use m_ss
   use m_ioutil, only : luavail, opnieee, opntext, clsieee, clstext
   use util
   implicit NONE
!
! !PUBLIC TYPES:
!
   Private
   Public gg_vect      ! gaussian fields type

! !PUBLIC MEMBER FUNCTIONS:
!
   Public gg_init      ! initializes a gaussian grid state vector
   Public gg_clean     ! deallocates memory
   Public gg_null      ! nullify pointers
   Public gg_put       ! writes gaussian grid file to a grads file
   Public gg_get       ! reads  gaussian grid file from a grads file

   interface gg_init 
      module procedure gg_init_
   end interface
   interface gg_clean
      module procedure gg_clean_
   end interface
   interface gg_put
      module procedure gg_put_
   end interface
   interface gg_get
      module procedure gg_get_
   end interface  
   interface gg_null
      module procedure gg_null_
   end interface 
!
! !DESCRIPTION: This module defines data types and methods for dealing
!               with gaussian grid state vectors derived from NCEP's
!               spectral space state variables.
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  02Aug2003  Todling  Initilized rc=0 in a number of routines.
!  15Mar2004  Todling  - Added slot for could water fraction; revised init interf.
!                      - cleaned up init_
!  28Sep2004  Todling  Nullifying pointers at initialization
!  10May2004  Todling  Renamed var cwf to cmwr for consistency w/ GSI.
!
!EOP
!-------------------------------------------------------------------------

!BOC

!  Public Types

!  SSI gaussian grid vector
!  ------------------------

   type gg_vect
     real,pointer :: hs(:,:)   =>null()  ! surface height in m
     real,pointer :: ps(:,:)   =>null() ! surface pressure in Pa
     real,pointer :: p(:,:,:)  =>null() ! pressure in Pa
     real,pointer :: dp(:,:,:) =>null() ! delta pressure in Pa
     real,pointer :: t(:,:,:)  =>null() ! dry bulb temperature in K
     real,pointer :: q(:,:,:)  =>null() ! specific humidity (kg/kg)
     real,pointer :: u(:,:,:)  =>null() ! zonal wind (m/s)
     real,pointer :: v(:,:,:)  =>null() ! meridional wind 
     real,pointer :: o(:,:,:)  =>null() ! ozone
     real,pointer :: cwmr(:,:,:)=>null()! cloud water condensate mixing ratio
   end type

!EOC
!---------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  gg_init_ --- Initializes a gaussian grid state 
!
! !INTERFACE:
!
   subroutine  gg_init_ ( nsig, glon, glat, gg, rc )
!
! !USES:
!
   implicit NONE
!
! !INPUT PARAMETERS:
!
   integer, intent(in)   :: nsig         ! number of sigma levels
   integer, intent(in)   :: glon, glat   ! gaussian grid dimensions
!
! !OUTPUT PARAMETERS:
!
   type(gg_vect),intent(inout)  :: gg  ! Gaussian grid vector 
   integer,intent(out)          :: rc  ! error return code
!
! !DESCRIPTION: Initializes the gg state vector. 
!               array sizes are not those specified in the ss 
!               metadata but specified in the arguments 
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  15Mar2004  Todling  Added cloud-water fraction.
!  15Dec2004  Kokron   Added OMP paralelizm.  Change to .F90 from .f90
!  30Mar2005  Todling  Somehow omp implementation was commented out - uncommented
!
!EOP
!-------------------------------------------------------------------------

   integer :: err,i

! allocate data
! -------------

   rc = 0
#if   (openmp)
!$omp parallel do default (shared), private (i)
#endif
do i=1,10
if (i==1) call myalloc(gg%hs,glon,glat)
if (i==2) call myalloc(gg%ps,glon,glat)
if (i==3) call myalloc(gg%u,glon,glat,nsig)
if (i==4) call myalloc(gg%v,glon,glat,nsig)
if (i==5) call myalloc(gg%q,glon,glat,nsig)
if (i==6) call myalloc(gg%p,glon,glat,nsig)
if (i==7) call myalloc(gg%dp,glon,glat,nsig)
if (i==8) call myalloc(gg%t,glon,glat,nsig)
if (i==9) call myalloc(gg%o,glon,glat,nsig)
if (i==10) call myalloc(gg%cwmr,glon,glat,nsig)
enddo

   end subroutine gg_init_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  gg_clean_ --- Deallocates memory used by gaussian grid vector
!
! !INTERFACE:
!
   subroutine  gg_clean_ ( gg ) 
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(gg_vect),intent(inout) :: gg ! gaussian grid vector

! !DESCRIPTION: Deallocates memory used by gg state vector
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  15Mar2004  Todling  Added cloud-water fraction.
!
!EOP
!-------------------------------------------------------------------------

   if(associated(gg%hs)) deallocate(gg%hs)
   if(associated(gg%ps)) deallocate(gg%ps)
   if(associated(gg%u))  deallocate(gg%u)
   if(associated(gg%v))  deallocate(gg%v)
   if(associated(gg%q))  deallocate(gg%q)
   if(associated(gg%p))  deallocate(gg%p)
   if(associated(gg%dp)) deallocate(gg%dp)
   if(associated(gg%t))  deallocate(gg%t)
   if(associated(gg%o))  deallocate(gg%o)
   if(associated(gg%cwmr))deallocate(gg%cwmr)
   
   end subroutine gg_clean_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  gg_null_ --- Nullify pointers used by gaussian grid vector
!
! !INTERFACE:
!
   subroutine  gg_null_ ( gg ) 
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(gg_vect),intent(inout) :: gg ! gaussian grid vector

! !DESCRIPTION: Deallocates memory used by gg state vector
!
! !REVISION HISTORY:
!
!  28Sep2004  Todling  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   if(associated(gg%hs)) nullify(gg%hs)
   if(associated(gg%ps)) nullify(gg%ps)
   if(associated(gg%u))  nullify(gg%u)
   if(associated(gg%v))  nullify(gg%v)
   if(associated(gg%q))  nullify(gg%q)
   if(associated(gg%p))  nullify(gg%p)
   if(associated(gg%dp)) nullify(gg%dp)
   if(associated(gg%t))  nullify(gg%t)
   if(associated(gg%o))  nullify(gg%o)
   if(associated(gg%cwmr))nullify(gg%cwmr)
   
   end subroutine gg_null_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  gg_put_ --- write out a gg vector
!
! !INTERFACE:
!
   subroutine  gg_put_ ( ss, nsig, glon, glat, gg, fname, rc)
!
! !USES:
!
  implicit NONE

!
! !INPUT PARAMETERS:
!
   type(ss_vect),intent(in) :: ss           ! spectral space vector
   type(gg_vect),intent(in) :: gg           ! gaussian grid vector
   integer, intent(in)      :: nsig         ! number of sigma levels
   integer, intent(in)      :: glon, glat   ! gaussian grid dimensions
   character*(*)            :: fname        ! grads file name
!
! !OUTPUT PARAMETERS:
!
   integer, intent(out)     :: rc

!
! !DESCRIPTION: write a gg vector to a grads file with its ctl 
!               counterpart
!
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!
!EOP
!-------------------------------------------------------------------------

   integer :: k,n,jhr,err, lugr,luctl
   character(len=32) :: gradsfile, ctlfile
   integer:: idat(8)
   character(len=10)  :: cdat(8)
   real(4), allocatable, dimension(:,:) :: s1
   real(4) :: xx(5)
   integer, allocatable, dimension(:) :: nlev

! start

   rc = 0
   gradsfile=fname//'.dat'
   ctlfile=fname//'.ctl'
   allocate(s1(glon,glat),stat=err)
   if(err/=0) then
     print *,' no more memory in m_gg'
     stop
   end if
   allocate(nlev(nsig))
   do k=1,nsig
     nlev(k)=k
   end do

   lugr=luavail()
   call opnieee(lugr,gradsfile,'unknown',err)
   if(err.ne.0) call myexit ('m_gg', 'error opening grads file',2 )

   s1=gg%hs; write(lugr)s1/9.81
   s1=gg%ps; write(lugr)s1
   do k=1,nsig
     s1=gg%dp(:,:,k); write(lugr)s1
   enddo
   do k=1,nsig
     s1=gg%t(:,:,k); write(lugr)s1
   enddo
   do k=1,nsig
     s1=gg%u(:,:,k); write(lugr)s1
   enddo
   do k=1,nsig
     s1=gg%v(:,:,k); write(lugr)s1
   enddo
   do k=1,nsig
     s1=gg%q(:,:,k); write(lugr)s1
   enddo

   call clsieee(lugr,err)

   luctl=luavail()
   call opntext(luctl,ctlfile,'unknown',err)
   if(err.ne.0) call myexit ( 'm_gg', 'error opening grads ctl file',6 )
   xx = 0. ; xx(2) = ss%meta%fhour
   call w3movdat((/xx(1), xx(2), xx(3), xx(4), xx(5)/),&
                 (/ss%meta%idate(4),ss%meta%idate(2),ss%meta%idate(3),0,&
                   ss%meta%idate(1),0,0,0/),idat)
   call w3pradat(idat,cdat)
   jhr=12

   write(luctl,'("dset ^",a)') trim(gradsfile)
   write(luctl,'("options sequential")')
   write(luctl,'("undef -9.99E+33")')
   write(luctl,'("title debug state")')
   write(luctl,'("xdef",i6," linear",2f12.6)') glon,0.d0,360.d0/glon
   write(luctl,'("ydef",i6," linear",2f12.6)') glat,-90.d0,180.d0/glat
   write(luctl,'("zdef",i6," levels")') nsig
   write(luctl,'(7i4,1x)') nlev
   write(luctl,'("tdef",i6," linear ",i2.2,"Z",i2.2,a3,i4.4,1x,i6,"hr")')&
                  1,idat(5),idat(3),cdat(2)(1:3),idat(1),jhr
   write(luctl,'("vars",i6)') 7
   write(luctl,'("HS  ",i3," 99 surface orography (m)")') 1
   write(luctl,'("PS  ",i3," 99 surface pressure (Pa)")') 1
   write(luctl,'("DP  ",i3," 99 delta pressure (Pa)")') nsig
   write(luctl,'("T   ",i3," 99 temperature (K)")') nsig
   write(luctl,'("U   ",i3," 99 zonal wind (m/s)")') nsig
   write(luctl,'("V   ",i3," 99 meridional wind (m/s)")') nsig
   write(luctl,'("Q   ",i3," 99 specific humidity (kg/kg)")') nsig
   write(luctl,'("endvars")')

   call clstext(luctl,err)

   deallocate(s1,nlev)

   end subroutine gg_put_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  gg_get_ --- read gaussian grid vector from a file
!
! !INTERFACE:
!
   subroutine  gg_get_ ( ss, gg, fname, rc )
!
! !USES:
!
   implicit NONE
!
! !INPUT PARAMETERS:
!
   type(ss_vect),intent(inout) :: ss      ! spectral space vector
   character*(*)               :: fname   ! grads file name
!
! !OUTPUT PARAMETERS:
!
   type(gg_vect),intent(out)   :: gg  ! gaussian grid vector
   integer,intent(out)         :: rc  ! return code

!
! !DESCRIPTION: read gaussian grid vector from a grads file
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!
!EOP
!-------------------------------------------------------------------------

   integer :: luin
   integer :: k, nsig, glon, glat, ios
   real(4), allocatable, dimension(:,:) :: s1

! start

   nsig=ss%meta%nsig
   glon=ss%meta%lonf
   glat=ss%meta%latf
   allocate(s1(glon,glat),stat=ios)
   if(ios/=0) then
     print *,' no more memory in m_gg'
     stop
   end if

   rc=1

   luin = luavail()
   call opnieee(luin,fname,'old',ios)
   read(luin) s1; gg%hs=s1 
   read(luin) s1; gg%ps=s1 
   do k=1,nsig
     read(luin) s1; gg%dp(:,:,k)=s1
   enddo
   do k=1,nsig
     read(luin) s1; gg%t(:,:,k)=s1
   enddo
   do k=1,nsig
     read(luin) s1; gg%u(:,:,k)=s1
   enddo
   do k=1,nsig
     read(luin) s1; gg%v(:,:,k)=s1
   enddo
   do k=1,nsig
     read(luin) s1; gg%q(:,:,k)=s1
   enddo

   call clsieee(luin, ios)
   deallocate(s1) 
   rc=0

   end subroutine gg_get_

!-------------------------------------------------------------------------

   end module m_gg
