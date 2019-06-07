!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_ss --- FVSSI spectral space state class
!
! !INTERFACE:
!
   module m_ss
!
!USES:
!
   use m_ioutil,     only : luavail, opnieee, clsieee
   use m_ioutil,     only : opntext, clstext
   use m_die,        only : die
   use util
   implicit NONE
!
! !PUBLIC TYPES:
!
   Private
   Public ss_meta     ! spectral space vector metadata
   Public ss_vect     ! spectral space vector
!
! !PUBLIC MEMBER FUNCTIONS:
!
   Public ss_init      ! initializes a spectral space state
   Public ss_clean     ! deallocates memory
   Public ss_null      ! nullify pointers
   Public ss_put       ! writes spectral file 
   Public ss_put_sigio ! writes spectral file using sigio_module
   Public ss_get       ! reads  spectral file 
   Public ss_get_meta  ! reads  spectral file metadata

!
! !DESCRIPTION: This module defines data types and methods for dealing
!               with spectral space fields from NCEP's SSI.
!
! !REMARKS:
!  29Oct2004  Todling  Commented out pointers initialized with =>null()
!                      since code crashes on the Compaq when this is
!                      done; this is possibly indication of a bug in this
!                      module.
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  28Set2004  Todling  Nullified initialization of pointers; 
!                      added interfaces and changed names accordingly.
!  15Feb2005  Dee      Optionally handle eta coordinate produced by fv2ss
!  24Sep2009  Sienkiewicz   add output option via sigio_module
!
!EOP
!-------------------------------------------------------------------------

!BOC

   interface ss_init
      module procedure ss_init_
   end interface
   interface ss_clean
      module procedure ss_clean_
   end interface
   interface ss_null
      module procedure ss_null_
   end interface
   interface ss_put
      module procedure ss_put_
   end interface
   interface ss_put_sigio
      module procedure ss_put_sigio_
   end interface
   interface ss_get
      module procedure ss_get_
   end interface
   interface ss_get_meta
      module procedure ss_get_meta_
   end interface
   interface ss_put_meta
      module procedure ss_put_meta_
   end interface
    interface ss_create_meta
      module procedure ss_create_meta_
   end interface
   
   integer,parameter:: ss_lhead1=32   ! length of first header record
   integer,parameter:: ss_levmax=100  ! maximum allowed number of levels
   integer,parameter:: ss_nwext=44    ! word length of header extension

!  Public Types

!  analysis file metadata
!  ----------------------

   type ss_meta
     character(ss_lhead1):: clabsig  ! Character(sigio_lmeta1) ON85 label
     real*4              :: fhour    ! forecast hour
     integer :: idate(4)             ! initial date 
                                     ! (hour, month, day, 4-digit year)
     logical :: eta                  ! true if eta vertical coordinate			     
     real*4, pointer     :: si(:)    ! sigma interfaces if eta=.false.
     real*4, pointer     :: sl(:)    ! sigma levels if eta=.false.
     real*4, pointer     :: ak(:), bk(:)  ! vertical grid coefficients if eta=.true.
     integer :: jcap             ! spectral truncation
     integer :: nsig             ! number of levels
     integer :: itrun            ! truncation flag (=1 for triangular)
     integer :: iorder           ! coefficient order flag (=2 for ibm order)
     integer :: irealf           ! floating point flag (=1 for ibm)
     integer :: igen             ! model generating flag
     integer :: latf             ! number of latitudes in dynamics (=(jcap+1)*3/2)
     integer :: lonf             ! number of longitudes in dynamics 
                                             ! (>=(jcap+1)*3 appropriate for fft)
     integer :: latb             ! number of latitudes in physics
     integer :: lonb             ! number of longitudes in physics
     integer :: latr             ! number of latitudes in radiation 
     integer :: lonr             ! number of longitudes in radiation
     integer :: ntrac            ! number of tracers
     integer :: icen2            ! subcenter id
     integer :: iens(2)          ! ensemble ids
     integer :: idpp             ! processing id
     integer :: idsl             ! semi-lagrangian id
     integer :: idvc             ! vertical coordinate id
     integer :: idvm             ! mass variable id
     integer :: idvt             ! tracer variable id
     integer :: idrun            ! run id
     integer :: idusr            ! user-defined id
     real*4  :: pdryini          ! global mean dry air pressure (kPa)
     integer :: ncldt            !  number of cloud types
   end type

!  Sigma file second header record
!  -------------------------------

   type ss_meta2
     sequence
     real*4   :: fhour               ! forecast hour
     integer  :: idate(4)            ! initial date
     real*4   :: sisl(2*ss_levmax+1) ! sigma information
     real*4   :: akbk(2*ss_levmax+2) ! eta grid information
     real*4   :: ext(ss_nwext)       ! model identifiers
   end type

!  SSI vector 
!  ----------

   type ss_vect
     type(ss_meta) :: meta
!    real,pointer :: hs(:)   =>null()   ! coefficients of surface height in m
!    real,pointer :: ps(:)   =>null()   ! coefficients of log of surface pressure over 1 kPa 
!    real,pointer :: t(:,:)  =>null()   ! coefficients of virtual temperature by level in K
!    real,pointer :: d(:,:)  =>null()   ! coefficients of divergence by level in 1/second
!    real,pointer :: z(:,:)  =>null()   ! coefficients of vorticity by level in 1/second
!    real,pointer :: q(:,:,:)=>null()   ! coefficients of tracers by level and tracer number
!                                       ! (mixing ratio, ozone, cloud water cond mix ratio)
     real,pointer :: hs(:)              ! coefficients of surface height in m
     real,pointer :: ps(:)              ! coefficients of log of surface pressure over 1 kPa 
     real,pointer :: t(:,:)             ! coefficients of virtual temperature by level in K
     real,pointer :: d(:,:)             ! coefficients of divergence by level in 1/second
     real,pointer :: z(:,:)             ! coefficients of vorticity by level in 1/second
     real,pointer :: q(:,:,:)           ! coefficients of tracers by level and tracer number
                                        ! (mixing ratio, ozone, cloud water cond mix ratio)
   end type

!
!EOC
!---------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_init_ --- Initializes a spectral space states
!
! !INTERFACE:
!
  subroutine  ss_init_ ( fname, ss, rc, &
                         jcap, nsig, nymd, nhms, fhour, ak, bk )  ! optional

!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS:
!
    character(len=255),intent(in)  :: fname  ! input file name
    integer, intent(in), optional  :: jcap   ! triangular truncation
    integer, intent(in), optional  :: nsig   ! number of levels
    integer, intent(in), optional  :: nhms   ! Time: hour-min-sec
    integer, intent(in), optional  :: nymd   ! Date: year-month-day
    integer, intent(in), optional  :: fhour  ! Forecast time
    real   , intent(in), optional  :: ak(:), bk(:) ! eta grid coefficients
!
! !OUTPUT PARAMETERS:
!
    type(ss_vect),intent(out)   :: ss  ! SSI vector
    integer,intent(out)         :: rc  ! error return code:
!
! !DESCRIPTION: Initializes the spectral space state.
!               
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  09Dec2003  Todling  Added forecast time as opt arg
!  02Sep2004  Todling  Zeroing out init arrays
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname='ss_init_'
   integer :: myfhour
   integer :: nc,dim1,dim2,dim3q,err
   integer :: id
  
! start

   rc = 0
   myfhour = 6       ! default is a 6 hour forecast 
   if (present(fhour)) myfhour = fhour

!  Either get metadata from an existing file (ss2fv)
!  or create metadata for an output file (fv2ss)
!  ------------------------------------------------
   if (.not.present(jcap)) then
      call ss_get_meta ( fname, ss, rc )
      dim2=ss%meta%nsig
   else
      dim2 = nsig
      if (.not.present(ak) .or. .not.present(bk)) then
      
         allocate(ss%meta%si(nsig+1), stat=err)
         allocate(ss%meta%sl(nsig  ), stat=err)
      
         if(err.ne.0) then
            call mywarn(myname,'problems allocating memory')
            rc = 2
            return
         end if
      
         call ss_create_meta ( ss, jcap, nsig, nymd, nhms, myfhour, rc )
         
      else ! inherit fv eta coordinate defined by ak, bk
      
         allocate(ss%meta%ak(nsig+1), stat=err)
         allocate(ss%meta%bk(nsig+1), stat=err)
      
         if(err.ne.0) then
            call mywarn(myname,'problems allocating memory')
            rc = 2
            return
         end if
	 
         call ss_create_meta ( ss, jcap, nsig, nymd, nhms, myfhour, rc, ak, bk )   
	 
      end if

   end if
   
!  allocate data
!  ------------- 

   if(associated(ss%hs)  .or. &
      associated(ss%ps)  .or. &
      associated(ss%t )  .or. &
      associated(ss%d )  .or. &
      associated(ss%z )  .or. &
      associated(ss%q )) then
      rc = 1
      return
   endif

   nc=(ss%meta%jcap+1)*(ss%meta%jcap+2)
   dim1=nc
   dim3q=ss%meta%ntrac
   allocate(ss%hs(dim1),    stat=err)
   allocate(ss%ps(dim1),    stat=err)
   allocate(ss%t(dim1,dim2), stat=err)
   allocate(ss%d(dim1,dim2), stat=err)
   allocate(ss%z(dim1,dim2), stat=err)
   allocate(ss%q(dim1,dim2,dim3q),stat=err)
   if(err.ne.0) then
     call mywarn(myname,'problems allocating memory')
      rc = 2
     return
   end if

!  Zero out to avoid garbage
!  -------------------------
   ss%hs  = 0.0
   ss%ps  = 0.0
   ss%t   = 0.0
   ss%d   = 0.0
   ss%z   = 0.0
   ss%q   = 0.0

   end subroutine ss_init_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_clean_ --- Deallocates memory used by spectral space state
!
! !INTERFACE:
!
   subroutine  ss_clean_ ( ss ) 
!
! !USES:
!
   implicit NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(ss_vect),intent(inout) :: ss  ! spectral space vector


! !DESCRIPTION: Deallocates memory used by spectral space state.
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!
!EOP
!-------------------------------------------------------------------------

   if(associated(ss%hs)) deallocate(ss%hs)
   if(associated(ss%ps)) deallocate(ss%ps)
   if(associated(ss%t))  deallocate(ss%t ) 
   if(associated(ss%d))  deallocate(ss%d )
   if(associated(ss%z))  deallocate(ss%z )
   if(associated(ss%q))  deallocate(ss%q )

   end subroutine ss_clean_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_null_ --- Nullify pointers used by spectral space state
!
! !INTERFACE:
!
   subroutine  ss_null_ ( ss ) 
!
! !USES:
!
   implicit NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(ss_vect),intent(inout) :: ss  ! spectral space vector


! !DESCRIPTION: Nullify pointers used by spectral space state.
!
! !REVISION HISTORY:
!
!  28Sep2004  Todling  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   if(associated(ss%hs)) nullify(ss%hs)
   if(associated(ss%ps)) nullify(ss%ps)
   if(associated(ss%t))  nullify(ss%t ) 
   if(associated(ss%d))  nullify(ss%d )
   if(associated(ss%z))  nullify(ss%z )
   if(associated(ss%q))  nullify(ss%q )

   end subroutine ss_null_
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_put_sigio_ --- writes out SSI vector using sigio_module
!
! !INTERFACE:
!
   subroutine  ss_put_sigio_ ( fname, ss, rc )
!
! !USES:
!
  use sigio_module
  implicit NONE

!
! !INPUT PARAMETERS:
!
    character(len=*), intent(in)   :: fname  ! output file name
    type(ss_vect), intent(inout)   :: ss     ! spectral space vector
   
!
! !OUTPUT PARAMETERS:
!
    integer, intent(out)   :: rc    ! error return code:
                                  !  0 - all is well
                                  !  >0 - errors
!
! !DESCRIPTION: writes a file containing fields in spectral coordinates 
!                   via sigio\_module, to be used by the SSI/GSI analysis
!
!
! !REVISION HISTORY:
!
!  24Sep2009  Sienkiewicz  Initial code.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter :: myname='ss_put_sigio'
  
    integer   :: jcap, nc, err, izero, i, ii, ii1, k, kk, l, m, n
    real      :: zero, zero1, one
  
    integer(sigio_intkind)      :: luout, irets
    type(sigio_head)            :: head
    type(sigio_data)            :: data
  

    real*4,  allocatable, dimension(:) :: s1
!  factsml  - factor to ensure proper scalar coefficients are zero
!  factvml  - factor to ensure proper vector coefficients are zero
    real, allocatable, dimension(:) :: factsml,factvml

    rc = 1

    jcap = ss%meta%jcap
    nc=(jcap+1)*(jcap+2)   
    allocate(s1(nc), stat=err)
    if(err.ne.0) then
       call mywarn(myname,'problems allocating memory')
       return
    end if

    head%ivs = 198410
    head%clabsig = ss%meta%clabsig
    head%fhour = ss%meta%fhour
    head%idate = ss%meta%idate
    head%jcap  = ss%meta%jcap
    head%levs  = ss%meta%nsig
    head%itrun = ss%meta%itrun
    head%iorder= ss%meta%iorder
    head%irealf= ss%meta%irealf
    head%igen  = ss%meta%igen
    head%latf  = ss%meta%latf
    head%lonf  = ss%meta%lonf
    head%latb  = ss%meta%latb
    head%lonb  = ss%meta%lonb
    head%latr  = ss%meta%latr
    head%lonr  = ss%meta%lonr
    head%ntrac = ss%meta%ntrac
    head%icen2 = ss%meta%icen2
    head%iens  = ss%meta%iens
    head%idpp  = ss%meta%idpp
    head%idsl  = ss%meta%idsl
    head%idvm  = ss%meta%idvm
    head%idvt  = ss%meta%idvt
    head%idrun = ss%meta%idrun
    head%idusr = ss%meta%idusr
    head%pdryini = ss%meta%pdryini
    head%ncldt   = ss%meta%ncldt
    if (ss%meta%eta) then
       head%idvc = 2
       head%ak(1:ss%meta%nsig+1) = ss%meta%ak(ss%meta%nsig+1:1:-1)
       head%bk(1:ss%meta%nsig+1) = ss%meta%bk(ss%meta%nsig+1:1:-1)
       head%nvcoord = 2
    else
       head%idvc = 1
       head%si   = ss%meta%si
       head%sl   = ss%meta%sl
       head%nvcoord = 1
    endif
    head%ixgr = 0
    
!   not yet sure about setting the vcoord array
    call sigio_alhead(head,irets)
    if (irets .ne. 0) then
       call mywarn(myname,'problem setting vcoord in head')
       return
    end if

    if (ss%meta%eta) then
       head%vcoord(1:head%levs+1,1)=head%ak(1:head%levs+1)
       head%vcoord(1:head%levs+1,2)=head%bk(1:head%levs+1)
    else
       head%vcoord(1:head%levs+1,1)=head%si(1:head%levs+1)
    end if

    call sigio_aldata(head,data,irets)
    if (irets .ne. 0) then
       call mywarn(myname,'problems allocating data structure')
       return
    end if

   ! setup factsml,factvml
   ! these factors ensure the complex parts are zeroed out
    allocate(factsml(nc),factvml(nc), stat=err)
    if(err.ne.0) then
       call mywarn(myname,'problems allocating memory')
       return
    end if
    izero=0; zero=0.0; one=1.0
    ii=-1; ii1=izero
    do l=izero,jcap
       zero1=float(min(1,l))
       do m=izero,jcap-l
          ii=ii+2; ii1=ii1+2
          factsml(ii)=one; factsml(ii1)=zero1
          factvml(ii)=one; factvml(ii1)=zero1
       end do
    end do
    factvml(1)=zero
    
    do i = 1,nc
       data%hs(i) = ss%hs(i)*factsml(i)
       data%ps(i) = ss%ps(i)*factsml(i)
    end do
    
    kk = 0
    do k=ss%meta%nsig,1,-1   !  note levels are wriiten in reverse vertical order
       kk = kk + 1
       do i=1,nc
          data%t(i,kk) = ss%t(i,k)*factsml(i)
          data%d(i,kk) = ss%d(i,k)*factvml(i)
          data%z(i,kk) = ss%z(i,k)*factvml(i)
          do n = 1,3
             data%q(i,kk,n) = ss%q(i,k,n)*factsml(i)
          end do
       end do
    end do
   
    luout = luavail()
    
    call sigio_swohdc(luout,fname,head,data,irets)
    
    if (irets .eq. 0) rc = 0
    
    return

  end subroutine ss_put_sigio_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_put_ --- writes out single instance of SSI vector
!
! !INTERFACE:
!
   subroutine  ss_put_ ( fname, ss, rc, &
                        grads )        ! optionals
!
! !USES:
!
  implicit NONE

!
! !INPUT PARAMETERS:
!
    character(len=*), intent(in)   :: fname  ! output file name
    type(ss_vect), intent(inout)   :: ss     ! spectral space vector

    logical, intent(in), optional  :: grads  ! controls grads file output
   
!
! !OUTPUT PARAMETERS:
!
  integer, intent(out)   :: rc    ! error return code:
                                  !  0 - all is well
                                  !  >0 - errors
!
! !DESCRIPTION: writes a file containing the SSI vector to be used
!               by the SSI analysis
!
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  02Apr2004  Todling  Added grads output option
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname='ss_put'

   character(len=255) :: ctlfile, grdfile
   integer :: nc, k, n, i, nvars, luctl, jhr
   integer :: id, id2, oid, ios,jcap,l,m,ii,ii1,izero,err
   integer :: idat(8)
   character(len=10)  :: cdat(8)
   real :: zero1,one,zero
   real(4) :: xx(5)
   logical :: gradsout

   integer, allocatable, dimension(:) :: nlev
   real*4,  allocatable, dimension(:) :: s1
   real*4,  allocatable, dimension(:,:) :: s2
!  factsml  - factor to ensure proper scalar coefficients are zero
!  factvml  - factor to ensure proper vector coefficients are zero
   real, allocatable, dimension(:) :: factsml,factvml
   logical :: fex
 
! start

   rc = 1
   jcap = ss%meta%jcap
   nc=(jcap+1)*(jcap+2)   
   allocate(s1(nc), stat=err)
   if(err.ne.0) then
     call mywarn(myname,'problems allocating memory')
     return
   end if

   gradsout = .false.
   if (present(grads)) then
       if(grads) gradsout = .true.
   endif

   if (gradsout) then
       allocate(s2(jcap+1,jcap+2),stat=err)
       if(err.ne.0) then
          call mywarn(myname,'problems allocating memory for s2')
         return
       end if
   endif

!  Write metadata
!  --------------

   call ss_put_meta ( fname, ss, rc )

!  Open the file
!  -------------

   rc = 1
   id=luavail() 
   call opnieee(id,fname,'append',ios)
   if(ios.ne.0) then
     call mywarn(myname,'cannot open ss file')
     return
   end if

   if (gradsout) then
       grdfile = trim(fname) // '.grd'
       id2=luavail() 
       call opnieee(id2,grdfile,'append',ios)
       if(ios.ne.0) then
         call mywarn(myname,'cannot open grads file')
         return
       end if
   end if

   ! setup factsml,factvml
   ! these factors ensure the complex parts are zeroed out
   allocate(factsml(nc),factvml(nc), stat=err)
   if(err.ne.0) then
     call mywarn(myname,'problems allocating memory')
     return
   end if
   izero=0; zero=0.0; one=1.0
   ii=-1; ii1=izero
   do l=izero,jcap
     zero1=float(min(1,l))
     do m=izero,jcap-l
       ii=ii+2; ii1=ii1+2
       factsml(ii)=one; factsml(ii1)=zero1
       factvml(ii)=one; factvml(ii1)=zero1
     end do
   end do
   factvml(1)=zero

!  write topography and pressure
!  -----------------------------
   
   do i=1,nc
     s1(i) = ss%hs(i)*factsml(i)
   end do
   write(id,iostat=ios)s1 
   if (gradsout) then
       s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
       write(id2,iostat=ios)s2
   endif
   if(ios.ne.0) then
     call mywarn(myname,'cannot write topography')
     return
   end if
   do i=1,nc
     s1(i) = ss%ps(i)*factsml(i)
   end do
   write(id,iostat=ios)s1 
   if (gradsout) then
       s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
       write(id2,iostat=ios)s2
   endif
   if(ios.ne.0) then
     call mywarn(myname,'cannot write pressure')
     return
   end if

!  write  temperature
!  ------------------
   
   do k=ss%meta%nsig,1,-1   !  note levels are wriiten in reverse vertical order
     do i=1,nc
       s1(i) = ss%t(i,k)*factsml(i)
     end do
     write(id,iostat=ios)s1
     if (gradsout) then
         s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
         write(id2,iostat=ios)s2
     endif
     if(ios.ne.0) then
       call mywarn(myname,'cannot write temperature')
       return
     end if
   enddo
 
!  write  div and vor
!  ------------------
    
   do k=ss%meta%nsig,1,-1
     do i=1,nc
       s1(i) = ss%d(i,k)*factvml(i)
     end do
     write(id,iostat=ios)s1 
     if (gradsout) then
         s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
        write(id2,iostat=ios)s2
     endif
     if(ios.ne.0) then
       call mywarn(myname,'cannot write divergence')
       return
     end if
     do i=1,nc
       s1(i) = ss%z(i,k)*factvml(i)
     end do
     write(id,iostat=ios)s1
     if (gradsout) then
         s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
         write(id2,iostat=ios)s2
     endif
     if(ios.ne.0) then
       call mywarn(myname,'cannot write vorticity')
       return
     end if
   enddo
 
!  write  tracers
!  --------------

! ozone coefs

   write(*,*) 'm_ss : writing ',ss%meta%ntrac,' tracers'   
   ! first write out q
   do n=1,ss%meta%ntrac      ! n=1 contains q, n=2 contains ozone, n=3 cloud cond mix ratio
     do k=ss%meta%nsig,1,-1
       do i=1,nc
         s1(i) = ss%q(i,k,n)*factsml(i)
       end do
       write(id,iostat=ios)s1
       if (gradsout) then
           s2(:,:) = reshape(s1,(/jcap+1,jcap+2/))
           write(id2,iostat=ios)s2
       endif
       if(ios.ne.0) then
         call mywarn(myname,'cannot write q')
         return
       end if
     enddo
   enddo

!  Close file
!  ----------
   call clsieee(id,ios)
   if(ios.ne.0) then
     call mywarn(myname,'cannot close ss file')
     return
   end if

   if (gradsout) then
       call clsieee(id2,ios)
       if(ios.ne.0) then
         call mywarn(myname,'cannot close grads file')
         return
       end if
   end if

   deallocate(s1,factsml,factvml)

   rc = 0
   if(.not.gradsout) return

   deallocate(s2)

!  Write out grads table for this file
!  -----------------------------------
   luctl=luavail()
   ctlfile = trim(fname) // '.ctl'
   call opntext(luctl,ctlfile,'unknown',err)
   if (err.ne.0) then 
       call mywarn ( myname, 'error opening grads ctl file' )
       rc = 1
       return
   endif
   xx = 0. ; xx(2) = ss%meta%fhour
   call w3movdat((/xx(1), xx(2), xx(3), xx(4), xx(5)/),&
                 (/ss%meta%idate(4),ss%meta%idate(2),ss%meta%idate(3),0,&
                   ss%meta%idate(1),0,0,0/),idat)
   call w3pradat(idat,cdat)
   jhr=12

   allocate(nlev(ss%meta%nsig))
   do k=ss%meta%nsig,1,-1
      nlev(k)=k
   end do

   write(luctl,'("dset ^",a)') trim(fname)
   write(luctl,'("options sequential")')
   write(luctl,'("Undef -9.99E+33")')
   write(luctl,'("title debug state")')
   write(luctl,'("xdef",i6," linear",3f12.6)') jcap+1,float(1),float(jcap+1)
   write(luctl,'("ydef",i6," linear",2f12.6)') jcap+2,float(1),float(jcap+2)
   write(luctl,'("zdef",i6," levels")') ss%meta%nsig
   write(luctl,'(7i4,1x)') nlev
   write(luctl,'("tdef",i6," linear ",i2.2,"Z",i2.2,a3,i4.4,1x,i6,"hr")')&
                  1,idat(5),idat(3),cdat(2)(1:3),idat(1),jhr

   deallocate(nlev)

   nvars = 7
   if(ss%meta%ntrac>=2) nvars = nvars + 1
   if(ss%meta%ntrac>=3) nvars = nvars + 1
   write(luctl,'("vars",i6)') nvars
   write(luctl,'("HS  ",i3," 99 surface orography (m)")') 1
   write(luctl,'("PS  ",i3," 99 surface pressure (Pa)")') 1
   write(luctl,'("DP  ",i3," 99 delta pressure (Pa)")') ss%meta%nsig
   write(luctl,'("T   ",i3," 99 temperature (K)")') ss%meta%nsig
   write(luctl,'("DV  ",i3," 99 divergency")') ss%meta%nsig
   write(luctl,'("VOR ",i3," 99 vorticity")') ss%meta%nsig
   write(luctl,'("Q   ",i3," 99 specific humidity (kg/kg)")') ss%meta%nsig
   if(ss%meta%ntrac>=2) write(luctl,'("O3  ",i3," 99 ozone mixing ratio (ppmv)")') ss%meta%nsig
   if(ss%meta%ntrac>=3) write(luctl,'("CW  ",i3," 99 cloud water cond mix ratio (kg/kg)")') ss%meta%nsig
   write(luctl,'("endvars")')

   call clstext(luctl,err)
   if (err/=0) then
       call mywarn ( myname, 'error closing grads ctl file' )
       rc = 1
       return
   endif
   print *, 'wrote grads table file to ', trim(ctlfile)

   rc = 0

   end subroutine ss_put_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_put_meta_ --- put metadata to SSI analysis file
!
! !INTERFACE:
!
   subroutine ss_put_meta_ ( fname, ss, rc )

! !USES

   implicit none
!
! !INPUT PARAMETERS
!
   character(len=*), intent(in)   :: fname  ! output file name
   type(ss_vect),intent(inout)    :: ss     ! spectral space vector
!
! !OUTPUT PARAMETERS
!
    integer,intent(out) :: rc  ! error return code
!
! !DESCRIPTION: this routine writes the metadata ( header ) to an
!               SSI analysis file
!
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname='ss_put_meta_'
   type(ss_meta2) :: meta2
   integer :: ios, id, nsig
   logical :: eta

! start

   rc=0
   nsig = ss%meta%nsig
   eta  = ss%meta%eta

   meta2%fhour   = ss%meta%fhour
   meta2%idate   = ss%meta%idate
   if (.not. eta) then
     meta2%sisl(1:nsig+1)=ss%meta%si
     meta2%sisl(nsig+2:2*nsig+1)=ss%meta%sl
     meta2%sisl(2*nsig+2:)=0
   else
     meta2%akbk(1:nsig+1) = ss%meta%ak
     meta2%akbk(nsig+2:2*nsig+2) = ss%meta%bk
     meta2%akbk(2*nsig+3:)=0
   end if
   meta2%ext(1)  = ss%meta%jcap
   meta2%ext(2)  = ss%meta%nsig
   meta2%ext(3)  = ss%meta%itrun
   meta2%ext(4)  = ss%meta%iorder
   meta2%ext(5)  = ss%meta%irealf
   meta2%ext(6)  = ss%meta%igen
   meta2%ext(7)  = ss%meta%lonf
   meta2%ext(8)  = ss%meta%latf
   meta2%ext(9)  = ss%meta%lonb
   meta2%ext(10) = ss%meta%latb
   meta2%ext(11) = ss%meta%lonr
   meta2%ext(12) = ss%meta%latr
   meta2%ext(13) = ss%meta%ntrac
   meta2%ext(14) = ss%meta%icen2
   meta2%ext(15:16)=ss%meta%iens
   meta2%ext(17) = ss%meta%idpp
   meta2%ext(18) = ss%meta%idsl
   meta2%ext(19) = ss%meta%idvc
   meta2%ext(20) = ss%meta%idvm
   meta2%ext(21) = ss%meta%idvt
   meta2%ext(22) = ss%meta%idrun
   meta2%ext(23) = ss%meta%idusr
   meta2%ext(24) = ss%meta%pdryini
   meta2%ext(25) = ss%meta%ncldt
   meta2%ext(26:ss_nwext)=0
! 
!  Open file
!  ---------
  
   rc=1
   id = luavail() 
   call opnieee(id,fname,'unknown',ios)
   if(ios.ne.0) then
     return
   end if

!  write header
!  ------------ 
   
   write(id,iostat=ios) ss%meta%clabsig
   if(ios.ne.0) return
   if (eta) then
      call myout( 6, myname, 'writing header for spectral-eta file' )
      write(id,iostat=ios) meta2%fhour,meta2%idate,meta2%akbk,meta2%ext
   else
      call myout( 6, myname, 'writing header for spectral-sigma file' )
      write(id,iostat=ios) meta2%fhour,meta2%idate,meta2%sisl,meta2%ext
   end if
   if(ios.ne.0) return

!  close file
!  ----------
   call clsieee(id,ios)
   if(ios.ne.0) call die('ss_put_meta','cannot close file')

   rc=0
 
   end subroutine ss_put_meta_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_get_ --- reads single instance of SSI vector
!
! !INTERFACE:
!
   subroutine  ss_get_ ( fname, ss, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS:
!
   character(len=*), intent(in) :: fname  ! input file name
!
! !OUTPUT PARAMETERS:
!
    type(ss_vect),intent(inout) :: ss ! spectral space vector 
    integer,intent(out)         :: rc ! return code

!
! !DESCRIPTION: reads a file containing the SSI vector 
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  18Aug2003  Todling  Turned ss intent into inout for Halem
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname='ss_get_'
   integer, parameter :: READ_ONLY = 1
   integer :: i,nc,n,k,ios,id,jcap,l,m,ii,ii1,izero,err
   real :: zero1,one,zero
   real*4, allocatable, dimension(:) :: s1
!  factsml  - factor to ensure proper scalar coefficients are zero
!  factvml  - factor to ensure proper vector coefficients are zero
   real, allocatable,dimension(:)::factsml,factvml

! start

   rc = 1

   jcap = ss%meta%jcap
   nc = (jcap+1)*(jcap+2)   
   allocate(s1(nc), stat=err)  
   if (err.ne. 0) then
     call mywarn(myname,'problems allocating memory')
     rc = 2
     return
   end if

!  Open file
!  ---------

   id=luavail()
   call opnieee(id,fname,'old',ios)
   if(ios.ne.0) then
     call mywarn(myname,'cannot open file')
     return
   end if

   ! setup factsml,factvml
   allocate(factsml(nc),factvml(nc), stat=err)
   if (err.ne. 0) then
     call mywarn(myname,'problems allocating memory')
     rc = 2
     return
   end if
   izero=0; zero=0.0; one=1.0
   ii=-1; ii1=izero
   do l=izero,jcap
     zero1=float(min(1,l))
     do m=izero,jcap-l
       ii=ii+2; ii1=ii1+2
       factsml(ii)=one; factsml(ii1)=zero1
       factvml(ii)=one; factvml(ii1)=zero1
     end do
   end do
   factvml(1)=zero
 
   read(id,iostat=ios)
   if(ios.ne.0) then 
     call mywarn(myname,'cannot read first blank')
     return
   end if
   read(id,iostat=ios)
   if(ios.ne.0) then 
     call mywarn(myname,'cannot read second blank')
     return
   end if

!  read topography and pressure
!  ----------------------------

   read(id,iostat=ios)s1 
   if(ios.ne.0) then 
     call mywarn(myname,'cannot read topography')
     return
   end if
   ss%hs=s1
   do i=1,nc
     ss%hs(i) = ss%hs(i)*factsml(i)
   end do

   read(id,iostat=ios)s1 
   if(ios.ne.0) then 
     call mywarn(myname,'cannot read pressure')
     return
   end if
   ss%ps=s1
   do i=1,nc
     ss%ps(i) = ss%ps(i)*factsml(i)
   end do

!  read temperature
!  ----------------

   do k=ss%meta%nsig,1,-1    ! note levels are in reverse order
     read(id,iostat=ios)s1   ! SSI's top is FV's bottom
     if(ios.ne.0) then 
       call mywarn(myname,'cannot read temperature')
       return
     end if
     ss%t(:,k)=s1            ! this is taken care of in regrid routines 
     do i=1,nc
       ss%t(i,k) = ss%t(i,k)*factsml(i)
     end do
   enddo
 
!  read div and vor
!  ----------------

   do k=ss%meta%nsig,1,-1    
     read(id,iostat=ios)s1
     if(ios.ne.0) then 
       call mywarn(myname,'cannot read divergence')
       return
     end if
     ss%d(:,k)=s1
     do i=1,nc
       ss%d(i,k) = ss%d(i,k)*factvml(i)
     end do
     read(id,iostat=ios)s1
     if(ios.ne.0) then 
       call mywarn(myname,'cannot read vorticity')
       return
     end if
     ss%z(:,k)=s1
     do i=1,nc
       ss%z(i,k) = ss%z(i,k)*factvml(i)
     end do
   enddo

!  read tracers
!  ------------

   do n=1,ss%meta%ntrac
     do k=ss%meta%nsig,1,-1  
       read(id,iostat=ios)s1
       if(ios.ne.0) then 
         call mywarn(myname,'cannot read tracers')
         return
       end if
       ss%q(:,k,n)=s1
       do i=1,nc
         ss%q(i,k,n) = ss%q(i,k,n)*factsml(i)
       end do
     enddo
   enddo

!  close file
!  ----------
   call clsieee(id, ios)
   if(ios.ne.0) then
     call mywarn(myname,'cannot close file')
     return
   end if

   deallocate(s1,factsml,factvml)

   rc=0

   end subroutine ss_get_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ss_get_meta_ --- get metadata from SSI analysis file
!
! !INTERFACE:
!
   subroutine  ss_get_meta_ ( fname, ss, rc )
!
! !USES:
!
   implicit NONE
!
! !INPUT PARAMETERS:
!
   character(len=*), intent(in)   :: fname  ! input file name
!
! !OUTPUT PARAMETERS:
!
   type(ss_vect),intent(out)   :: ss    ! spectral space vector
   integer, intent(out)        :: rc    ! error return code
!
! !DESCRIPTION: this routine reads the metadata ( header ) from an
!               input  SSI analysis file
!
! !REVISION HISTORY:
!
!  05Nov2002  Cruz     Initial code.
!  25Mar2004  Ravi/RT  Error check to read old NCEP spectral files.
!  08Mar2005  Todling  Added 382 configuration.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter  :: myname = 'ss_get_meta_'
   real*4, parameter  :: one4 = 1.
   type(ss_meta2)     :: meta2
   integer            :: ios, id, err, nsig

! start

   rc=0

!  Open the file
!  -------------
  
   id = luavail() 
   call opnieee(id,fname,'old',ios)
   if(ios.ne.0) then
     call mywarn(myname,'cannot open spectral file')
     rc = 3
     return
   end if

!  retrieve the metadata fields
!  into the second header record
!  ----------------------------   
 
   read(id,iostat=ios) ss%meta%clabsig
   if(ios.ne.0) then
     call mywarn(myname,'cannot metadata (label) from spectral file')
     rc = 3
     return
   end if

   read(id,iostat=ios) meta2%fhour, &
                       meta2%idate, &
                       meta2%akbk,  & 
                       meta2%ext
   if (ios==0) then
      ss%meta%eta = .true.
      call myout( 6, myname, 'ss metadata read for eta vertical coordinate' )
   else
      backspace id
      read(id,iostat=ios) meta2%fhour, &
                          meta2%idate, &
                          meta2%sisl,  &
                          meta2%ext
      if (ios==0) then
         ss%meta%eta = .false.
         call myout( 6, myname, 'ss metadata read for sigma vertical coordinate' )
      else
                  ! allow reading of old NCEP spec files
         call mywarn(myname,'this appears to be an OLD NCEP spectral file')
         meta2%ext = 0
         backspace id
         read(id,iostat=ios) meta2%fhour, &
                             meta2%idate, &
                             meta2%sisl,  &
                             meta2%ext(1:6)
         if (ios==0) then
            ss%meta%eta = .false.
            call myout( 6, myname, 'ss metadata read for sigma vertical coordinate, old style' )
         else
            call mywarn(myname,'cannot read metadata from spectral file')
            rc = 3
            return
         end if
      end if
   end if
   
!  initialize metadata with input
!  ------------------------------
 
   nsig=meta2%ext(2)
   ss%meta%fhour = meta2%fhour
   ss%meta%idate = meta2%idate
   ss%meta%jcap  = meta2%ext(1)
   ss%meta%nsig  = meta2%ext(2)
   if (.not. ss%meta%eta) then
     allocate(ss%meta%si(nsig+1), stat=err);   if (err.ne. 0) rc = 2
     allocate(ss%meta%sl(nsig  ), stat=err);   if (err.ne. 0) rc = 2
     if (rc .ne. 0) then
       call mywarn(myname,'problems allocating memory')
       return
     end if
     ss%meta%si=meta2%sisl(1:nsig+1)
     ss%meta%sl=meta2%sisl(nsig+2:2*nsig+1)
   else
     allocate(ss%meta%ak(nsig+1), stat=err);   if (err.ne. 0) rc = 2
     allocate(ss%meta%bk(nsig+1), stat=err);   if (err.ne. 0) rc = 2
     if (rc .ne. 0) then
       call mywarn(myname,'problems allocating memory')
       return
     end if
     ss%meta%ak = meta2%akbk(1:nsig+1)  
     ss%meta%bk = meta2%akbk(nsig+2:2*nsig+2)  
   end if
   ss%meta%itrun = meta2%ext(3)
   ss%meta%iorder= meta2%ext(4)
   ss%meta%irealf= meta2%ext(5)
   ss%meta%igen  = meta2%ext(6)
   ss%meta%lonf  = meta2%ext(7)
   ss%meta%latf  = meta2%ext(8)
   ss%meta%lonb  = meta2%ext(9)
   ss%meta%latb  = meta2%ext(10)
   ss%meta%lonr  = meta2%ext(11)
   ss%meta%latr  = meta2%ext(12)
   ss%meta%ntrac = max(meta2%ext(13),one4)
   ss%meta%icen2 = meta2%ext(14)
   ss%meta%iens  = meta2%ext(15:16)
   ss%meta%idpp  = meta2%ext(17)
   ss%meta%idsl  = meta2%ext(18)
   ss%meta%idvc  = meta2%ext(19)
   ss%meta%idvm  = meta2%ext(20)
   ss%meta%idvt  = meta2%ext(21)
   ss%meta%idrun = meta2%ext(22)
   ss%meta%idusr = meta2%ext(23)
   ss%meta%pdryini=meta2%ext(24)
   ss%meta%ncldt = meta2%ext(25)

   call myout( 6, myname, 'ss vector info' )
   write(*,*) '   get date  = ',ss%meta%idate
   write(*,*) '   get hour  = ',ss%meta%fhour
   write(*,*) '   get jcap  = ',ss%meta%jcap
   write(*,*) '   get nsig  = ',ss%meta%nsig
   write(*,*) '   get lonf  = ',ss%meta%lonf
   write(*,*) '   get latf  = ',ss%meta%latf
   write(*,*) '   get ntrac = ',ss%meta%ntrac
   if (ss%meta%eta) then
   write(*,*) ' vertical coordinate is eta:'
   write(*,*) '   get ak(top)  = ',ss%meta%ak(1)
   write(*,*) '   get ak(bot)  = ',ss%meta%ak(nsig+1)
   write(*,*) '   get bk(top)  = ',ss%meta%bk(1)
   write(*,*) '   get bk(bot)  = ',ss%meta%bk(nsig+1)
   else
   write(*,*) ' vertical coordinate is sigma:'
   write(*,*) '   get sigi(bot)  = ',ss%meta%si(1)
   write(*,*) '   get sigi(top)  = ',ss%meta%si(nsig)
   write(*,*) '   get sigl(bot)  = ',ss%meta%sl(1)
   write(*,*) '   get sigl(top)  = ',ss%meta%sl(nsig)
   end if

! check if jcap info is missing or zero

   if(ss%meta%jcap==0) then
     call myexit(myname,'JCAP is zero', -1)
   end if

! The SSI does not write ss%meta%lonf and ss%meta%latf so ...

   if(ss%meta%lonf==0 .or. ss%meta%latf==0) then
     select case (ss%meta%jcap)
       case (62)
         ss%meta%lonf=192
         ss%meta%latf=96
       case (126)
         ss%meta%lonf=384
         ss%meta%latf=192
       case (170)
         ss%meta%lonf=384
         ss%meta%latf=192
       case (254)
         ss%meta%lonf=512
         ss%meta%latf=258
       case (382)
         ss%meta%lonf=768
         ss%meta%latf=386
       case default  ! T170
         ss%meta%lonf=384
         ss%meta%latf=192
     end select
   end if

!  close file
!  ----------

   call clsieee(id,ios)

   rc = 0

   end subroutine ss_get_meta_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ss_create_meta_ create SSI analysis file metadata

!
! !INTERFACE:
!
   subroutine ss_create_meta_ ( ss, jcap, nsig, nymd, nhms, fhour, rc, &
                                                      ak, bk ) ! optional
!
! !USES
!
   implicit none
!
! !INPUT PARAMETERS
!
   integer, intent(in) :: jcap       ! triangular truncation
   integer, intent(in) :: nsig       ! number of layers 
   integer, intent(in) :: nhms       ! eta hour
   integer, intent(in) :: nymd       ! eta date
   integer, intent(in) :: fhour      ! forecast hour
   real   , intent(in), optional  :: ak(:), bk(:) ! eta grid coefficients
!
! !INPUT/OUTPUT PARAMETERS
!
   type(ss_vect), intent(out) :: ss  ! spectral space vector
!
! !OUTPUT PARAMETERS
!
   integer,intent(out) :: rc   ! return code
!
! !DESCRIPTION: create variables that make up the metadata of NCEP's
!               SSI analysis file
!
!
! !REVISION HISTORY:
!
!  05Nov2002  Cruz     Initial code.
!  02Apr2004  Todling  Fixed the definition of analysis time.
!  23Jun2004  Todling  Added 32 level case; generated from using
!                      lowpass filter on the corresponding 64 pressure levs
!                      of 64 level case [matlab: decimate(pe,2,'FIR')].
!  08Mar2005  Todling  Added 382 configuration.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter  :: myname = 'ss_create_meta_'
   integer :: nymda, nhmsa
   integer :: err

! header elements
! ---------------

   integer :: itrun,iorder,irealf,igen, &
              latf,lonf,latb,lonb,latr,lonr, &
              ntrac,icen2,iens(2),idpp,idsl, &
              idvc,idvm,idvt,idrun,idusr,ncldt
   real*4 :: pdryini
   real*4, allocatable, dimension(:) :: sigi,sigl

! start

   call myout( 6, myname, 'create a metadata header' )

   rc=0
   
   select case (jcap)
     case (62)
       lonf=192
       latf=96
     case (126)
       lonf=384  
       latf=192
     case (170)
       lonf=384  
       latf=192
     case (254)
       lonf=512
       latf=258
     case (382)
       lonf=768
       latf=386
     case default  ! T62 
       lonf=192
       latf=96
   end select

   if (present(ak) .and. present(bk)) then
   
      ss%meta%eta = .true.
      ss%meta%ak  = ak
      ss%meta%bk  = bk
      ss%meta%idvc = 2
      
   else
   
      ss%meta%eta = .false.
      ss%meta%idvc = 1
      allocate( sigi(nsig+1), sigl(nsig), stat=err)
      if (err.ne. 0) rc = 2
      if (rc .ne. 0) return     ! problems allocating memory
  
      select case (nsig)
      case(28)
       sigl(28) = 0.002726
       sigl(27) = 0.010056
       sigl(26) = 0.018336
       sigl(25) = 0.028750
       sigl(24) = 0.041795
       sigl(23) = 0.058048
       sigl(22) = 0.078152
       sigl(21) = 0.102782
       sigl(20) = 0.132611
       sigl(19) = 0.168230
       sigl(18) = 0.210057
       sigl(17) = 0.258231
       sigl(16) = 0.312479
       sigl(15) = 0.372048
       sigl(14) = 0.435679
       sigl(13) = 0.501681
       sigl(12) = 0.568091
       sigl(11) = 0.632896
       sigl(10) = 0.694265
       sigl(9)  = 0.750756
       sigl(8)  = 0.801416
       sigl(7)  = 0.845785
       sigl(6)  = 0.883844
       sigl(5)  = 0.915917
       sigl(4)  = 0.942547
       sigl(3)  = 0.964373
       sigl(2)  = 0.982082
       sigl(1)  = 0.994997
       sigi(29) = 0.000000
       sigi(28) = 0.006570
       sigi(27) = 0.013860
       sigi(26) = 0.023090
       sigi(25) = 0.034690
       sigi(24) = 0.049200
       sigi(23) = 0.067230
       sigi(22) = 0.089450
       sigi(21) = 0.116540
       sigi(20) = 0.149160
       sigi(19) = 0.187830
       sigi(18) = 0.232860
       sigi(17) = 0.284210
       sigi(16) = 0.341370
       sigi(15) = 0.403340
       sigi(14) = 0.468600
       sigi(13) = 0.535290
       sigi(12) = 0.601350
       sigi(11) = 0.664820
       sigi(10) = 0.724010
       sigi(9)  = 0.777730
       sigi(8)  = 0.825270
       sigi(7)  = 0.866420
       sigi(6)  = 0.901350
       sigi(5)  = 0.930540
       sigi(4)  = 0.954590
       sigi(3)  = 0.974180
       sigi(2)  = 0.990000
       sigi(1)  = 1.000000
      case(32)
       sigi(1)  = 1.00000000000000
       sigi(2)  = 0.98872534607179
       sigi(3)  = 0.97419476058937
       sigi(4)  = 0.95579303249680
       sigi(5)  = 0.93257762063833
       sigi(6)  = 0.90373411346272
       sigi(7)  = 0.86840256952308
       sigi(8)  = 0.82594488796729
       sigi(9)  = 0.77604539486662
       sigi(10) = 0.71892139925941
       sigi(11) = 0.65545262507372
       sigi(12) = 0.58723094893478
       sigi(13) = 0.51646259559839
       sigi(14) = 0.44570663998292
       sigi(15) = 0.37752168467858
       sigi(16) = 0.31410709914775
       sigi(17) = 0.25704949441056
       sigi(18) = 0.20722067727026
       sigi(19) = 0.16482699621217
       sigi(20) = 0.12955418065981
       sigi(21) = 0.10074674794219
       sigi(22) = 0.07757556044339
       sigi(23) = 0.05916436309299
       sigi(24) = 0.04467802721353
       sigi(25) = 0.03336620491047
       sigi(26) = 0.02458555073375
       sigi(27) = 0.01780274961459
       sigi(28) = 0.01258026965530
       sigi(29) = 0.00857329557478
       sigi(30) = 0.00550074957890
       sigi(31) = 0.00315766281337
       sigi(32) = 0.00136097899391
       sigi(33) = 0.00000000000000
       sigl(1)  = 0.99436267303590
       sigl(2)  = 0.98146005333058   
       sigl(3)  = 0.96499389654308   
       sigl(4)  = 0.94418532656756
       sigl(5)  = 0.91815586705053  
       sigl(6)  = 0.88606834149290  
       sigl(7)  = 0.84717372874518   
       sigl(8)  = 0.80099514141695
       sigl(9)  = 0.74748339706302   
       sigl(10) = 0.68718701216656   
       sigl(11) = 0.62134178700425   
       sigl(12) = 0.55184677226659
       sigl(13) = 0.48108461779066   
       sigl(14) = 0.41161416233075   
       sigl(15) = 0.34581439191316   
       sigl(16) = 0.28557829677915
       sigl(17) = 0.23213508584041   
       sigl(18) = 0.18602383674122   
       sigl(19) = 0.14719058843599   
       sigl(20) = 0.11515046430100
       sigl(21) = 0.08916115419279   
       sigl(22) = 0.06836996176819   
       sigl(23) = 0.05192119515326   
       sigl(24) = 0.03902211606200
       sigl(25) = 0.02897587782211   
       sigl(26) = 0.02119415017417   
       sigl(27) = 0.01519150963494   
       sigl(28) = 0.01057678261504
       sigl(29) = 0.00703702257684   
       sigl(30) = 0.00432920619614   
       sigl(31) = 0.00225932090364   
       sigl(32) = 0.00068048949696
      case(42)
       sigi(1)  = 1.00000000
       sigi(2)  = 0.99197000
       sigi(3)  = 0.98273998
       sigi(4)  = 0.97215998
       sigi(5)  = 0.96006000
       sigi(6)  = 0.94625998
       sigi(7)  = 0.93061000
       sigi(8)  = 0.91293001
       sigi(9)  = 0.89306003
       sigi(10) = 0.87085998
       sigi(11) = 0.84619999
       sigi(12) = 0.81902999
       sigi(13) = 0.78930998
       sigi(14) = 0.75708002
       sigi(15) = 0.72245997
       sigi(16) = 0.68564999
       sigi(17) = 0.64691001
       sigi(18) = 0.60661000
       sigi(19) = 0.56515998
       sigi(20) = 0.52305001
       sigi(21) = 0.48076999
       sigi(22) = 0.43886000
       sigi(23) = 0.39780000
       sigi(24) = 0.35804999
       sigi(25) = 0.32001001
       sigi(26) = 0.28400999
       sigi(27) = 0.25029001
       sigi(28) = 0.21901000
       sigi(29) = 0.19025999
       sigi(30) = 0.16406000
       sigi(31) = 0.14036000
       sigi(32) = 0.11906000
       sigi(33) = 0.10005000
       sigi(34) = 0.08316000
       sigi(35) = 0.06824000
       sigi(36) = 0.05512000
       sigi(37) = 0.04362000
       sigi(38) = 0.03357000
       sigi(39) = 0.02482000
       sigi(40) = 0.01722000
       sigi(41) = 0.01063000
       sigi(42) = 0.00492000
       sigi(43) = 0.00000000
       sigl(1)  = 0.99597311
       sigl(2)  = 0.987361193
       sigl(3)  = 0.977443695
       sigl(4)  = 0.966113031
       sigl(5)  = 0.953153133
       sigl(6)  = 0.93842876
       sigl(7)  = 0.921754241
       sigl(8)  = 0.9029845
       sigl(9)  = 0.88193953
       sigl(10) = 0.858512819
       sigl(11) = 0.832589447
       sigl(12) = 0.804136455
       sigl(13) = 0.773154199
       sigl(14) = 0.739722669
       sigl(15) = 0.703998089
       sigl(16) = 0.666211843
       sigl(17) = 0.62668252
       sigl(18) = 0.585798085
       sigl(19) = 0.544008136
       sigl(20) = 0.501803696
       sigl(21) = 0.459700972
       sigl(22) = 0.418210596
       sigl(23) = 0.377800435
       sigl(24) = 0.338902682
       sigl(25) = 0.301882356
       sigl(26) = 0.267023146
       sigl(27) = 0.234525785
       sigl(28) = 0.204514682
       sigl(29) = 0.177044451
       sigl(30) = 0.152100176
       sigl(31) = 0.129605815
       sigl(32) = 0.109456532
       sigl(33) = 9.151221067E-2
       sigl(34) = 7.561223209E-2
       sigl(35) = 6.159678474E-2
       sigl(36) = 4.929006845E-2
       sigl(37) = 3.851688281E-2
       sigl(38) = 2.911660634E-2
       sigl(30) = 2.093769237E-2
       sigl(40) = 1.383117214E-2
       sigl(41) = 7.646807469E-3
       sigl(42) = 2.041562228E-3
      case(64)
       sigi(1)  = 1.
       sigi(2)  = 0.994670987
       sigi(3)  = 0.988632023
       sigi(4)  = 0.98180002
       sigi(5)  = 0.974083006
       sigi(6)  = 0.96538502
       sigi(7)  = 0.955603004
       sigi(8)  = 0.94463098
       sigi(9)  = 0.932359993
       sigi(10) = 0.918677986
       sigi(11) = 0.903479993
       sigi(12) = 0.88666302
       sigi(13) = 0.868139029
       sigi(14) = 0.847829998
       sigi(15) = 0.825685024
       sigi(16) = 0.801676989
       sigi(17) = 0.775811017
       sigi(18) = 0.748133004
       sigi(19) = 0.718729019
       sigi(20) = 0.687731028
       sigi(21) = 0.655315995
       sigi(22) = 0.621704996
       sigi(23) = 0.587159991
       sigi(24) = 0.551973999
       sigi(25) = 0.516462982
       sigi(26) = 0.480955005
       sigi(27) = 0.445778012
       sigi(28) = 0.411249012
       sigi(29) = 0.377658993
       sigi(30) = 0.345268995
       sigi(31) = 0.314300001
       sigi(32) = 0.284927994
       sigi(33) = 0.257283986
       sigi(34) = 0.231454
       sigi(35) = 0.207481995
       sigi(36) = 0.185371995
       sigi(37) = 0.165098995
       sigi(38) = 0.146607995
       sigi(39) = 0.129822999
       sigi(40) = 0.114655003
       sigi(41) = 0.101002
       sigi(42) = 8.875600249E-2
       sigi(43) = 7.780800015E-2
       sigi(44) = 6.804899871E-2
       sigi(45) = 5.936999992E-2
       sigi(46) = 5.167099833E-2
       sigi(47) = 4.485499859E-2
       sigi(48) = 3.883099928E-2
       sigi(49) = 3.351499885E-2
       sigi(50) = 2.882999927E-2
       sigi(51) = 2.470799908E-2
       sigi(52) = 2.108399943E-2
       sigi(53) = 1.790099964E-2
       sigi(54) = 1.510700025E-2
       sigi(55) = 1.265799999E-2
       sigi(56) = 1.051099971E-2
       sigi(57) = 8.631000295E-3
       sigi(58) = 6.984999869E-3
       sigi(59) = 5.54399984E-3
       sigi(60) = 4.284000024E-3
       sigi(61) = 3.183000023E-3
       sigi(62) = 2.219999908E-3
       sigi(63) = 1.378000015E-3
       sigi(64) = 6.419999991E-4
       sigi(65) = 0.E+0
       sigl(1)  = 0.997334659
       sigl(2)  = 0.991650403
       sigl(3)  = 0.985214591
       sigl(4)  = 0.977939725
       sigl(5)  = 0.969731688
       sigl(6)  = 0.960491061
       sigl(7)  = 0.950113237
       sigl(8)  = 0.938490689
       sigl(9)  = 0.925512969
       sigl(10) = 0.91107142
       sigl(11) = 0.895062089
       sigl(12) = 0.877389371
       sigl(13) = 0.857970178
       sigl(14) = 0.836740077
       sigl(15) = 0.813659906
       sigl(16) = 0.78871876
       sigl(17) = 0.761942089
       sigl(18) = 0.733395934
       sigl(19) = 0.703189373
       sigl(20) = 0.67147696
       sigl(21) = 0.638457835
       sigl(22) = 0.604373693
       sigl(23) = 0.569502294
       sigl(24) = 0.534148216
       sigl(25) = 0.498633713
       sigl(26) = 0.463286996
       sigl(27) = 0.428430676
       sigl(28) = 0.394368827
       sigl(29) = 0.361377567
       sigl(30) = 0.329697907
       sigl(31) = 0.299528241
       sigl(32) = 0.271022052
       sigl(33) = 0.244287685
       sigl(34) = 0.21939002
       sigl(35) = 0.196352869
       sigl(36) = 0.175165638
       sigl(37) = 0.155788153
       sigl(38) = 0.13815479
       sigl(39) = 0.122182935
       sigl(40) = 0.107777007
       sigl(41) = 9.483192116E-2
       sigl(42) = 8.323913068E-2
       sigl(43) = 7.288959622E-2
       sigl(44) = 6.367427856E-2
       sigl(45) = 5.548869073E-2
       sigl(46) = 4.823431745E-2
       sigl(47) = 4.181715846E-2
       sigl(48) = 3.614972159E-2
       sigl(49) = 3.115151823E-2
       sigl(50) = 2.675008588E-2
       sigl(51) = 2.287890576E-2
       sigl(52) = 1.947700977E-2
       sigl(53) = 1.648990251E-2
       sigl(54) = 1.386962179E-2
       sigl(55) = 1.157263666E-2
       sigl(56) = 9.559988044E-3
       sigl(57) = 7.797650062E-3
       sigl(58) = 6.254609209E-3
       sigl(59) = 4.904353991E-3
       sigl(60) = 3.723796224E-3
       sigl(61) = 2.691220725E-3
       sigl(62) = 1.787146204E-3
       sigl(63) = 9.936115239E-4
       sigl(64) = 2.663990017E-4
     end select
   
     ss%meta%si=sigi
     ss%meta%sl=sigl
   
     deallocate(sigi,sigl)
  
   end if
 
!  create idate using nhms, nymd
!  -----------------------------

   write(*,*) '   fvDAS nymd, nhms : ',nymd, nhms
   nymda = nymd ; nhmsa = nhms

   call tick (nymda, nhmsa, -fhour*3600)
   ss%meta%fhour = real(fhour)
   write(*,*) '   fvDAS nymda, nhmda : ',nymda, nhmsa, ss%meta%fhour

   ss%meta%idate(1)=nhmsa/10000
   ss%meta%idate(3)=nymda-100*(nymda/100)
   ss%meta%idate(2)=nymda-10000*(nymda/10000)-ss%meta%idate(3)
   ss%meta%idate(2)=ss%meta%idate(2)/100
   ss%meta%idate(4)=nymda/10000

   ss%meta%jcap=jcap
   ss%meta%nsig=nsig
   ss%meta%itrun=1
   ss%meta%iorder=2
   ss%meta%irealf=1
   ss%meta%igen=82
   ss%meta%lonf=lonf
   ss%meta%latf=latf
   ss%meta%lonb=lonf
   ss%meta%latb=latf
   ss%meta%lonr=lonf
   ss%meta%latr=latf
   ss%meta%ntrac=3
   ss%meta%icen2=0
   ss%meta%iens=0
   ss%meta%idpp=0
   ss%meta%idsl=0
!   ss%meta%idvc=0   !!! variable set above
   ss%meta%idvm=0
   ss%meta%idvt=21
   ss%meta%idrun=0
   ss%meta%idusr=0
   ss%meta%pdryini=98.3161392
   ss%meta%ncldt=1

   call myout( 6, myname, 'ss vector info' )
   write(*,*) '   set date = ',ss%meta%idate
   write(*,*) '   set hour = ',ss%meta%fhour
   write(*,*) '   set jcap = ',ss%meta%jcap
   write(*,*) '   set lonf = ',ss%meta%lonf
   write(*,*) '   set latf = ',ss%meta%latf
   write(*,*) '   set nsig = ',ss%meta%nsig
   if (ss%meta%eta) then
   write(*,*) ' vertical coordinate is eta:'
   write(*,*) '   set ak(top)  = ',ss%meta%ak(1)
   write(*,*) '   set ak(bot)  = ',ss%meta%ak(nsig+1)
   write(*,*) '   set bk(top)  = ',ss%meta%bk(1)
   write(*,*) '   set bk(bot)  = ',ss%meta%bk(nsig+1)
   else
   write(*,*) ' vertical coordinate is sigma:'
   write(*,*) '   set sigi(bot)  = ',ss%meta%si(1)
   write(*,*) '   set sigi(top)  = ',ss%meta%si(nsig)
   write(*,*) '   set sigl(bot)  = ',ss%meta%sl(1)
   write(*,*) '   set sigl(top)  = ',ss%meta%sl(nsig)
   end if

   end subroutine ss_create_meta_

  end module m_ss
