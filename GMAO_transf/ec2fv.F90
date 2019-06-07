   Program ec2fv_main

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ec2fv --- Reading the ECMWF ERA40 data.
!
!
! !USAGE: see the routine usage() below

  use m_dyn
  use m_set_eta, only  : set_eta
  use m_mapz
  use m_ecdyn
  use m_ec_set_eta
  use m_ec2fv, only : ec2fv
  use m_const, only : kappa
  use m_const, only : rgas
  use m_const, only : rvap
  use m_const, only : grav
  use m_const, only : cpm
  use m_const, only: undef


   Implicit NONE

! !DESCRIPTION: Remaping ECMWF Nature 91 layers to 72 Layers GEOS-5 grid.
!
! !REVISION HISTORY:
!
!  August 29 2006 Ravi C. Govindaraju Initial design.
!  July   22 2009 Ravi C. Govindaraju Added an option to use GEOS-5.4.0 analysis eta file.
!                         
!-----------------------------------------------------------------------------------------
!EOP


   character(len=*), parameter :: myname = 'ec2fv'

  
!                              -----------------------
!                               Hardwired Parameters
!                              -----------------------

      integer, parameter :: mFiles = 256       ! Max. number of input files
      integer, parameter :: mVars  = 256       ! Max. number of variables
      integer, parameter :: mLevs  = 256       ! Max. number of levels    

      type(dyn_ecvect)          :: w_ec   ! dynamics ECMWF state vector
      type(dyn_vect)            :: w_f    ! dynamics GMAO  state vector
      type(dyn_vect)            :: w_d    ! dynamics GMAO  state vector

!                              -----------------------
!                              User Defined Parameters
!                              -----------------------


      integer            :: nFiles             ! Actual number of input files
      integer            :: kount,err          ! Variable counter
      character(len=256) :: inFile             ! Input file names
      character(len=256) :: phisFile           ! Input surface Orography
      character(len=256) :: uaFile             ! Input diag_dyn File     
      character(len=256) :: anaFile            ! Input Analysis File    
      character(len=256) :: outFile            ! Output file name
      integer  :: freq                         ! increment hours specified from command line
      integer  :: prec                         ! precission 0= 32bit,1= 64 bit
      integer  :: nstepec,nstep                ! 

      integer  :: iflag                        ! Initial flag
      integer  :: irflag                       ! Monthly Mean flag

      real big
      real, pointer     :: lon(:)              ! longitudes in deg (im)
      real, pointer     :: lat(:)              ! latitudes in deg (jm)
      real, pointer     :: lev(:)              ! levels in hPa (km)
      real              :: rair,cpair,akap,pr
      real, pointer     :: pe0(:,:,:)
      real, pointer     :: dp0(:,:,:)
      real, pointer     :: pke0(:,:,:)
      real, pointer     :: pk0(:,:,:)
      real, pointer     :: pe1(:,:,:)
      real, pointer     :: ps1(:,:)
      real, pointer     :: pke1(:,:,:)
      real, pointer     :: phi1(:,:,:)

      real, pointer     :: pe2(:,:,:)
      real, pointer     :: dp2(:,:,:)
      real, pointer     :: pke2(:,:,:)
      real, pointer     :: thv2(:,:,:)

      real, pointer     :: u0(:,:,:)
      real, pointer     :: v0(:,:,:)
      real, pointer     :: thv0(:,:,:)
      real, pointer     :: q0(:,:,:,:)
      real, pointer     :: array(:,:)


 
      integer           :: nLevs = 0           ! total number of levels
      real, pointer     :: Levs(:)             ! vertical levels
      character(len=25) :: cLevs(mLevs)        ! Character reprsentation of levels

      integer           :: nVars,nq            ! Actual number of variables
      character(len=64) :: inVars(mVars)       ! Input  variable names (nVars)
      character(len=64) :: outVars(mVars)      ! output variable names (nVars)
      character(len=64) :: outUnits(mVars)     ! Units of output variables (nVars)
      character(len=64) :: uname               ! Uwind name
      character (len=1) resolution
    
      integer          :: outPrec              ! Output file precision:
                                               ! 0 = 32 bits,  1 = 64bits

!                              -----------------------
!                                Variable Work Space
!                              -----------------------

      real, pointer :: outField(:,:,:)         ! Ouput variable
      real, pointer ::  ofn(:,:,:)              ! Ouput variable


!                                  Local Work Space
!                              -----------------------

      integer count, iff, it, iv, itest, ii, i, j, k,ind2
      integer ig,is
      real dx,dx1, dx2,dy, dy1, dy2,ptop,pint,cp
      double precision pi
      logical ptop_found,fex,wsh,blend

 
      integer    :: vectype ! GEOS-4 or GEOS-5 dyn vect
      logical    :: dgrid   ! GEOS-4 or GEOS-5 switch for winds grid


!                              -----------------------
!                                  Output Meta Data
!                              -----------------------

      character(len=255) :: title              ! meta data title
      character(len=255) :: source             ! data source
      character(len=255) :: contact            ! contact org.   
      character(len=255) :: levunits           ! Vertical levels
      character(len=25)  :: append             ! im*jm
      real               :: missing_val

      integer, pointer :: yyyymmdd(:)          ! Date
      integer, pointer :: hhmmss(:)            !
      integer          :: ndate                ! Date
      integer          :: ndate_ec             ! Date
      integer          :: ndate_old            ! Date
      integer          :: yyyymmddp,hhmmssp    ! previous Date & time
      integer          :: ntimep               ! counter for total number of times previously accumulated.
      integer          :: ntime 
      integer          :: ntime_ec

      integer          :: fid                  ! input file ID
      integer          :: out_fid              ! output file ID
      integer          :: nkount
      integer          :: rc, rc1              ! return error code

!                              -----------------------
!                                  eta information 
!                              -----------------------

      integer           :: im_e,im             ! input zonal dimension       
      integer           :: jm_e,jm             ! input meridional dimension       
      integer           :: km_e,km             ! input vertical dimension    
      integer           :: lm_e,lm             ! input time dimension    
      integer           :: in_e                ! output zonal dimension       
      integer           :: jn_e                ! output meridional dimension       
      integer           :: kn_e                ! output vertical dimension    
      integer           :: ln_e                ! output vertical dimension    
      integer           :: in                  ! output zonal dimension       
      integer           :: jn                  ! output meridional dimension       

      integer           :: nVars_e             ! input number of variables   
      integer           :: kn            ! Output number of vertical levels
      integer           :: ks55                  
      integer           :: buf(3)
      real, pointer     :: lon_e(:)            ! longitudes in deg (im)
      real, pointer     :: lat_e(:)            ! latitudes in deg (jm)
      real, pointer     :: ak55(:)             ! vertical grid a coefficien
      real, pointer     :: bk55(:)             ! vertical grid a coefficien
      real, pointer     :: dpref(:)            ! vertical grid a coefficien
      real, pointer     :: inField(:,:,:)             ! Surface Pressure
      real, pointer     :: thv1(:,:,:)    
      real, pointer     :: pk(:,:,:)    
      real*8, pointer   :: lev_e(:)            ! levels in eta (km)
      integer, pointer  :: kmVar_e(:)          ! Number of vertical levels for variables
      integer           :: nymd,nhms,nymd_phis,nhms_phis

      real,pointer      :: akec(:)
      real,pointer      :: bkec(:)
      real              :: ptopec,pb,pa,p1,p2

      character(len=255) :: vtitle_in(mVars)   ! output title
      real              :: valid_range(2, mVars)
      real              :: packing_range(2, mVars)
      real              :: p,bkh,eps
      integer           :: ngatts              ! Number of attributes for GFIO
      integer           :: imin,jmin,xmin,imax,jmax,xmax,L
      integer           :: in_f, jn_f, kn_f,ln_f
      logical              initial,file_exist
!.................................................................................

    initial = .true.
    nkount = 0
    nq = 1
!  Get user input
!  --------------


    call  Init_ ( inFile,phisFile,anaFile,uaFile,outFile,ndate_ec,ntime_ec,wsh,blend,vectype,dgrid)
 
    cp    = cpm ! used to be getcon('CP') (ie, dry, but not consistent w/ fv)
    eps   = rvap/rgas-1.0

!    -------------------------------------------------
!       Retrieve ECMWF Nature analysis fields.
!    -------------------------------------------------

    call dyn_ecget ( inFile, ndate_ec, ntime_ec, w_ec, rc)

    if ( rc /= 0 )  then
      print *,' dyn_ecget:rc ',rc
      call die (myname, 'something wrong with dyn_ecget ')
    endif

!
!     -------------------------------------
!        Retrieve ECMWF resolution indexes.
!     -------------------------------------

    call Dyn_EcGetDim ( w_ec, im_e, jm_e, km_e, lm_e )
    print *,' im_e,jm_e,km_e,lm_e ',im_e,jm_e,km_e,lm_e

!    -------------------------------------------------
!       Retrieve ECMWF Nature Orography.
!    1/28/2008 Ravi   The following change is made (reading the phis and ps
!       from the phis file) to read the input data from
!       the Nature files.
!    -------------------------------------------------

    if(km_e == 61) then
     nymd_phis = 20000601
     nhms_phis = 000000
     call get_ecphis(phisFile,nymd_phis,nhms_phis,w_ec%phism%name,w_ec%phis,im_e,jm_e,rc)
    endif

    if(km_e == 91) then
     nymd_phis = ndate_ec
     nhms_phis = ntime_ec
     print *,' nymd_phis,nhms_phis ',nymd_phis,nhms_phis 
     call get_ecphis(phisFile,nymd_phis,nhms_phis,w_ec%phism%name,w_ec%phis,im_e,jm_e,rc)
     call get_ecphis(phisFile,nymd_phis,nhms_phis,w_ec%psm%name,w_ec%ps,im_e,jm_e,rc)
      print *,' get_ecphis:rc ',rc 
     if(rc == 0) then
      print *,' get_ecphis:rc ',rc 
      w_ec%ps = exp(w_ec%ps)
      print *,' get_ecphis:rcexp ',rc 
     else
      print *,' get_ecphis:rc ',rc
      call die (myname, 'something wrong with get_ecphis ')
     endif
    endif

!   call minmaxpe(w_ec%phis,im_e,jm_e,1)

!    -------------------------------------------
!        Allocate working arrays.
!    -------------------------------------------

    print *,'im_e,jm_e,km_e ',im_e,jm_e,km_e

    allocate(akec(km_e+1))
    allocate(bkec(km_e+1))
    allocate(ps1(im_e,jm_e))
    allocate(pe1(im_e,jm_e,km_e+1))
    allocate(pke1(im_e,jm_e,km_e+1))
    allocate(phi1(im_e,jm_e,km_e+1))
    allocate(thv1(im_e,jm_e,km_e))
    allocate(pk(im_e,jm_e,km_e))

!     -------------------------------------------
!     Extract ECMWF ak and bk coefficients.
!     -------------------------------------------
!
     nstepec = 18000   ! default so that all files from this point on have it specified
    call ec_set_eta(km_e+1,ptopec,akec,bkec)
    print *,'   km_e,ptopec ',km_e,ptopec
!   do k = 1,km_e+1
!    print *,'   ',k,akec(k),bkec(k)
!   end do


!       --------------------------------------
!          Compute ECMWF Pressure Delta.
!       --------------------------------------

! ----------------------------------
! Compute ECMWF edge-level pressures
! ----------------------------------

!       TEST 2/19/2008

    do k = 1,km_e+1
       pe1(:,:,k) = akec(k) + bkec(k) * w_ec%ps
    end do

     do j = 1,jm_e
      do i = 1,im_e
       if(pe1(i,j,1) <= 0.0) pe1(i,j,1) = 0.00001
      end do
     end do

! -----------------------

!   print *,' '
!   print *,' ec2fv:ecmwf:delp '
!   print *,' '

    do k = 1,km_e
       w_ec%delp(:,:,k) = pe1(:,:,k+1)-pe1(:,:,k)
!      w_ec%delp(:,:,k) = (akec(k+1)-akec(k)) +  &
!                         (bkec(k+1)-bkec(k)) * w_ec%ps
!      call minmaxpe(w_ec%delp(:,:,k),im_e,jm_e,k)
    end do
     

! ----------------------------------
! Compute ECMWF edge-level pressures
! ----------------------------------

!       TEST 2/19/2008

!     pe1(:,:,km_e+1) = w_ec%ps(:,:)
!     do k=km_e,1,-1
!      pe1(:,:,k) = pe1(:,:,k+1)-w_ec%delp(:,:,k)
!     enddo

! -----------------------

!     print *,' kappa,cp ',kappa,cp
!     print *,' '
!     print *,' '
!     print *,' *** ec2fv: pe1 ***'
!     print *,' '

!       TEST 2/19/2008

!      do j = 1,jm_e
!       do i = 1,im_e
!        pe1(i,j,1) = 1.0
!       end do
!      end do
! -----------------------

!
!        ----------------------
!         Compute ECMWF pkappa.
!        ----------------------
!
!     print *,' '
!     print *,' *** ec2fv: pe1 ***'
!     print *,' '

!     do k=1,km_e+1
!      call minmaxpe(pe1(:,:,k),im_e,jm_e,k)
!     enddo
!     print *,' '
!     print *,' *** ec2fv: pke1 ***'
!     print *,' '

      do k=1,km_e+1
       pke1(:,:,k) = pe1(:,:,k)**kappa
!      call minmaxpe(pke1(:,:,k),im_e,jm_e,k)
      enddo

!
! ----------------------------------------------
! Construct ECMWF virtual potential temperature
! ----------------------------------------------
      do k=1,km_e
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif

      
       do j=1,jm_e
        do i=1,im_e
         pk  (i,j,k) = ( pke1(i,j,k+1)-pke1(i,j,k) )/( kappa*log(pe1(i,j,k+1)/pe1(i,j,k)) )
         thv1(i,j,k) =     w_ec%pt(i,j,k)*( 1.0+eps*max(0.0,w_ec%q(i,j,k,1)) )/pk(i,j,k)
        enddo
       enddo
     enddo

! ---------------------------------
! Construct ECMWF analysis heights
! ---------------------------------

      phi1(:,:,km_e+1) = w_ec%phis(:,:)
      do k=km_e,1,-1
       phi1(:,:,k) = phi1(:,:,k+1) + cp*thv1(:,:,k)*( pke1(:,:,k+1)-pke1(:,:,k) )
      enddo


!     print *,' '
!     print *,' *** ec2fv: thv1 ***'
!     print *,' '
!     do k = 1,km_e
!      call minmaxpe(thv1(:,:,k),im_e,jm_e,k)
!     end do

!     print *,' '
!     print *,' *** ec2fv: phi1 ***'
!     print *,' '
!     do k = 1,km_e+1
!      call minmaxpe(phi1(:,:,k),im_e,jm_e,k)
!     end do

    print *,' Getting analysis ',trim(anaFile)

     call dyn_get ( anaFile, ndate, ntime, w_f,rc,vectype=vectype) 

    call Dyn_GetDim ( w_f, in_e, jn_e, kn_e,ln_e )
!
    if(vectype == 5) then
!      ------------------------
!      Ravi: Modfied to be consistent with qltot.

      w_f%grid%lm = 4
      ln_e = w_f%grid%lm
    endif
!   ------------------------

    print *,' in_e, jn_e, kn_e,ln_e ',in_e, jn_e, kn_e,ln_e

!     print *,' '
!     print *,' *** analysis: phis ***'
!     print *,' '
!     call minmaxpe(w_f%phis(:,:),in_e,jn_e,1)


     inquire ( FILE = trim(uaFile), EXIST = fex )
     if(fex) then
      print *,' Getting uafile ',trim(uaFile)
      call dyn_get ( uaFile, ndate, ntime, w_d, rc,vectype=vectype,nstep=nstepec)
      if(rc == 0) then
       call Dyn_GetDim ( w_d, in_f, jn_f, kn_f,ln_f )

!      ------------------------
      if(vectype == 5) then
!        Ravi: Modfied to be consistent with qltot.

         w_d%grid%lm = 4
         ln_f = w_d%grid%lm
      endif
!      ------------------------
!      print *,' in_f, jn_f, kn_f,ln_f ',in_f, jn_f, kn_f,ln_f
       w_f%ps = w_d%ps
       w_f%u = w_d%u
       w_f%v = w_d%v
       w_f%pt = w_d%pt
       w_f%delp = w_d%delp
       w_f%q = w_d%q
      else
       print *,' ec2fv: *** PROBLEM reading uaFile ***',trim(uaFile)
       stop 30
      endif
     endif

     if(vectype == 5) then
      if(.not. associated(array))  allocate(array(in_e,jn_e))
       ind2 = in_e/2
       array(1:ind2,:) = w_f%ps(ind2+1:in_e,:)
       array(ind2:in_e,:) = w_f%ps(1:ind2,:)
       w_f%ps =  array

       array(1:ind2,:) = w_f%phis(ind2+1:in_e,:)
       array(ind2:in_e,:) = w_f%phis(1:ind2,:)
       w_f%phis =  array

       do k = 1,kn_e
         array(1:ind2,:)    = w_f%u(ind2+1:in_e,:,k)
         array(ind2+1:in_e,:) = w_f%u(1:ind2,:,k)
         w_f%u(:,:,k)       = array

         array(1:ind2,:)    = w_f%v(ind2+1:in_e,:,k)
         array(ind2:in_e,:) = w_f%v(1:ind2,:,k)
         w_f%v(:,:,k)       = array

         array(1:ind2,:)    = w_f%pt(ind2+1:in_e,:,k)
         array(ind2:in_e,:) = w_f%pt(1:ind2,:,k)
         w_f%pt(:,:,k)       = array
       end do

       do l = 1,ln_e
        do k = 1,kn_e
         array(1:ind2,:)    = w_f%q(ind2+1:in_e,:,k,l)
         array(ind2:in_e,:) = w_f%q(1:ind2,:,k,l)
         w_f%q(:,:,k,l)       = array
        end do
       end do
     if(associated(array))  deallocate(array)
    endif

    allocate(ak55(kn_e+1))
    allocate(bk55(kn_e+1))

    allocate(pe0(in_e,jn_e,kn_e+1))
    allocate(pk0(in_e,jn_e,kn_e))
    allocate(dp0(in_e,jn_e,kn_e+1))
    allocate(pke0(in_e,jn_e,kn_e+1))
    allocate(u0(in_e,jn_e,kn_e))
    allocate(v0(in_e,jn_e,kn_e))
    allocate(thv0(in_e,jn_e,kn_e))
    allocate(q0(in_e,jn_e,kn_e,ln_e))

    allocate(pe2(in_e,jn_e,kn_e+1))
!   allocate(dp2(in_e,jn_e,kn_e+1))
    allocate(pke2(in_e,jn_e,kn_e+1))
    allocate(thv2(in_e,jn_e,kn_e))
!
!    -----------------------------------
!     Extract GEOS ak and bk coefficients.
!    -----------------------------------
!
    call set_eta(kn_e,ks55,ptop,pint,ak55,bk55)
!   print *,'   kn_e,ptop ',kn_e,ptop
!   do k = 1,kn_e+1
!    print *,'   ',k,ak55(k),bk55(k)
!   end do
!

!Compute new surface pressure (ps1) consistent with fv topography
! ----------------------------------------------------------------
                                                                                                                     
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do j=1,jm_e
       do i=1,im_e
           L = km_e
           do while ( phi1(i,j,L).lt.w_f%phis(i,j) )
!           print *,' phi1(',i,',',j,',',L,') ', phi1(i,j,L),' phis1(',i,',',j,')  ',w_f%phis(i,j)
            L = L-1
            if (L.lt.1) then
              write(6,*) 'ec2fv: Level variable is less than 1.'
              write(6,*) 'ec2fv: i,j,l have values ',i,j,l
              write(6,*) 'ec2fv: ABORT'
              call flush(6)
              call abort()
            endif
           enddo
           ps1(i,j) = pe1(i,j,L+1)*( 1+(phi1(i,j,L+1)-w_f%phis(i,j))/(cp*thv1(i,j,L)*pke1(i,j,L+1)) )**(1.0/kappa)
       enddo
      enddo
!     print *,' *** We are here ***'

! Compute edge-level pressures
! ----------------------------
!      print *,' '
!      print *,' ec2fv,ua: delp '
!      print *,' '
!      do  k = 1,kn_e+1
!       call minmaxpe(w_f%delp(:,:,k),in_e,jn_e,k)
!      end do

! -------------------------------------------------------------------------
! Construct fv pressure variables (pe1,pke1,dp1)
! -------------------------------------------------------------------------
                                                                                                                     
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,k)
#endif


       pe2(:,:,kn_e+1) = w_f%ps(:,:)
       do k = kn_e,1,-1
        pe2(:,:,k) = pe2(:,:,k+1)-w_f%delp(:,:,k)
       end do

!    --------------------------

!      print *,' '
!      print *,' ec2fv:ua: pe2 '
!      print *,' '
!      do  k = 1,kn_e+1
!       call minmaxpe(pe2(:,:,k),in_e,jn_e,k)
!      end do
      


!      print *,' '
!      print *,' ec2fv:ua: pke2 '
!      print *,' '
      do k=1,kn_e+1
       do j=1,jn_e
        do i=1,in_e
         pke2(i,j,k) = pe2(i,j,k)**kappa
        enddo
       enddo
!      call minmaxpe(pke2(:,:,k),in_e,jn_e,k)
      enddo
!     print *,' *** passed pke2 loop ***' 

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,k)
#endif

!     do k=1,kn_e
!     do j=1,jn_e
!     do i=1,in_e
!      dp2(i,j,k) = pe2(i,j,k+1)-pe2(i,j,k)
!     enddo
!     enddo
!     enddo

      do k = 1,kn_e
       do j = 1,jn_e
        do i = 1,in_e
         thv2(i,j,k) = w_f%pt(i,j,k)
        end do
       end do
      end do

!     print *,' *** passed thv2 loop *** '
! -------------------------------------------------------------------------
! Construct fv pressure variables (pe1,pke1,dp1) using new surface pressure
! -------------------------------------------------------------------------
                                                                                                                     
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,k)
#endif

      do k=1,kn_e+1
       do j=1,jn_e
        do i=1,in_e
          pe0(i,j,k) = ak55(k) + bk55(k)*ps1(i,j)
!         if(pe0(i,j,k) <= 0.0) pe0(i,j,k) = 0.00001
          pke0(i,j,k) = pe0(i,j,k)**kappa
        enddo
       enddo
      enddo

      do k=1,kn_e
       do j=1,jn_e
        do i=1,in_e
         pk0 (i,j,k) = ( pke0(i,j,k+1)-pke0(i,j,k) )/( kappa*log(pe0(i,j,k+1)/pe0(i,j,k)) )
        enddo
       enddo
      enddo

!     print *,' *** passed construction of new pe0 and pke0 ***'

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,k)
#endif

!     print *,' ** Newly constructed delp ** '
      do k=1,kn_e
       do j=1,jn_e
        do i=1,in_e
         dp0(i,j,k) = pe0(i,j,k+1)-pe0(i,j,k)
        enddo
       enddo
!       call minmaxpe(dp0(:,:,k),in_e,jn_e,k)
      enddo

!     print *, 'passed new delp loop'

      if(vectype == 4) then
       thv0 = thv2
      else
       thv0 = thv2 / pk0
      endif

!      print *,' '
!      print *,' ec2fv:ua: ps '
!      print *,' '
!      call minmaxpe(w_f%ps(:,:),in_e,jn_e,1)

!      print *,' '
!      print *,' ec2fv:ua: ps1 '
!      print *,' '
!      call minmaxpe(ps1(:,:),in_e,jn_e,1)

!      print *,' '
!      print *,' ec2fv:ua: pe0 '
!      print *,' '
!      do  k = 1,kn_e+1
!       call minmaxpe(pe0(:,:,k),in_e,jn_e,k)
!      end do

!      print *,' '
!      print *,' ec2fv:ua: pke0 '
!      print *,' '
!      do  k = 1,kn_e+1
!       call minmaxpe(pke0(:,:,k),in_e,jn_e,k)
!      end do

!     print *,' '
!     print *,' *** ec2fv: thv0 ***'
!     print *,' '
!     do k = 1,kn_e
!      call minmaxpe(thv0(:,:,k),im_e,jm_e,k)
!     end do

      w_ec%q(:,:,:,2) = w_ec%q(:,:,:,2) * 604229.0

!     print *,' '
!     print *,' *** ec2fv: q2 ***'
!     print *,' '
!     do k = 1,km_e
!      call minmaxpe(w_ec%q(:,:,k,2),im_e,jm_e,k)
!     end do

      print *,' *** Calling ec2fv ***' 

      pa = p1*100.0 ! convert to hPa
      pb = p2*100.0 ! convert to hPa

      u0 = w_f%u
      v0 = w_f%v
      q0 = w_f%q


      if(wsh) then
       if(blend) then
        call ec2fv(im_e,jm_e,km_e,nq, pke1, pe1,                     &
                    w_ec%u, w_ec%v, thv1, w_ec%q(:,:,:,1),           &
                    w_ec%q(:,:,:,2),w_ec%q(:,:,:,3),w_ec%q(:,:,:,4), & 
                    kn_e, pke2, pe2,pke0,pe0,                        &
                    u0, v0, thv0, q0(:,:,:,1),q0(:,:,:,2),           &
                    q0(:,:,:,3),q0(:,:,:,4),dgrid,pabove = pa,       &
                    pbelow = pb,wsh=wsh,blend=blend)
       else
        call ec2fv(im_e,jm_e,km_e,nq, pke1, pe1,                     &
                    w_ec%u, w_ec%v, thv1, w_ec%q(:,:,:,1),           &
                    w_ec%q(:,:,:,2),w_ec%q(:,:,:,3),w_ec%q(:,:,:,4), & 
                    kn_e, pke2, pe2,pke0,pe0,                        &
                    u0, v0, thv0, q0(:,:,:,1),q0(:,:,:,2),           &
                    q0(:,:,:,3),q0(:,:,:,4),dgrid,pabove = pa,       &
                    pbelow = pb,wsh=wsh)
       endif
      else
       if(blend) then
        call ec2fv(im_e,jm_e,km_e,nq, pke1, pe1,                     &
                    w_ec%u, w_ec%v, thv1, w_ec%q(:,:,:,1),           &
                    w_ec%q(:,:,:,2),w_ec%q(:,:,:,3),w_ec%q(:,:,:,4), & 
                    kn_e, pke2, pe2,pke0,pe0,                        &
                    u0, v0, thv0, q0(:,:,:,1),q0(:,:,:,2),           &
                    q0(:,:,:,3),q0(:,:,:,4),dgrid,pabove = pa,       &
                    pbelow = pb,blend=blend)
       else
        call ec2fv(im_e,jm_e,km_e,nq, pke1, pe1,                     &
                    w_ec%u, w_ec%v, thv1, w_ec%q(:,:,:,1),           &
                    w_ec%q(:,:,:,2),w_ec%q(:,:,:,3),w_ec%q(:,:,:,4), & 
                    kn_e, pke2, pe2,pke0,pe0,                        &
                    u0, v0, thv0, q0(:,:,:,1),q0(:,:,:,2),           &
                    q0(:,:,:,3),q0(:,:,:,4),dgrid,pabove = pa,       &
                    pbelow = pb)


       endif
      endif

       print *,' *** Out from ec2fv **** '

!     print *,' '
!     print *,' *** ec2fv: thv0 ***'
!     print *,' '
!     do k = 1,kn_e
!      call minmaxpe(thv0(:,:,k),im_e,jm_e,k)
!     end do
     
!     For now, set unnecessary 2D fields to missing
!     ---------------------------------------------
    w_f%ps    = ps1
    w_f%u     = u0
    w_f%v     = v0

    if(vectype == 5) then

      do k=1,kn_e
       do j=1,jn_e
        do i=1,in_e
         pk0 (i,j,k) = ( pke0(i,j,k+1)-pke0(i,j,k) )/( kappa*log(pe0(i,j,k+1)/pe0(i,j,k)) )
        enddo
       enddo
      enddo

      w_f%pt    = thv0 * pk0
    else
     w_f%pt    = thv0 
    endif

    w_f%delp  = dp0

    print *,' ec2fv:size:q ',size(w_f%q,4)
    w_f%q = q0

     if(vectype == 5) then
      if(.not. associated(array))  allocate(array(in_e,jn_e))
       ind2 = in_e/2
       array(1:ind2,:) = w_f%phis(ind2+1:in_e,:)
       array(ind2+1:in_e,:) = w_f%phis(1:ind2,:)
       w_f%phis =  array

       array(1:ind2,:) = w_f%ps(ind2+1:in_e,:)
       array(ind2+1:in_e,:) = w_f%ps(1:ind2,:)
       w_f%ps =  array

       do k = 1,kn_e
         array(1:ind2,:)    = w_f%u(ind2+1:in_e,:,k)
         array(ind2+1:in_e,:) = w_f%u(1:ind2,:,k)
         w_f%u(:,:,k)       = array

         array(1:ind2,:)    = w_f%v(ind2+1:in_e,:,k)
         array(ind2+1:in_e,:) = w_f%v(1:ind2,:,k)
         w_f%v(:,:,k)       = array

         array(1:ind2,:)    = w_f%pt(ind2+1:in_e,:,k)
         array(ind2+1:in_e,:) = w_f%pt(1:ind2,:,k)
         w_f%pt(:,:,k)       = array
       end do


       print *, ' ln_e ',ln_e
       do l = 1,ln_e
        do k = 1,kn_e
         array(1:ind2,:)    = w_f%q(ind2+1:in_e,:,k,l)
         array(ind2+1:in_e,:) = w_f%q(1:ind2,:,k,l)
         w_f%q(:,:,k,l)       = array
        end do
       end do
       if(associated(array))  deallocate(array)
      endif

!   w_f%q(:,:,:,1)  =  q0(:,:,:,1)
!   if(w_f%grid%lm > 2 ) w_f%q(:,:,:,2)  =  0.0
   
    nymd = ndate_ec
    nhms = ntime_ec
    prec = 0
!   print *,' ',trim(w_f%phism%name)
!   print *,' ',trim(w_f%tsm%name)
!   print *,' ',trim(w_f%um%name)
!   print *,' ',trim(w_f%vm%name)
!   print *,' ',trim(w_f%ptm%name)

    print *,' im,jm,km,lm ',w_f%grid%im,w_f%grid%jm,w_f%grid%km,w_f%grid%lm
!   print *,' lon ',w_f%grid%lon
!   print *,' lat ',w_f%grid%lat
!   print *,' nstepec ',nstepec
!
!       print *,'*** ec2fv: out from gmap ***'
!        print *, '  '
!        print *, '  '
!        print *," ec2fv:MIN/MAX of qltot0    ECMWF ",kn_e
!        print *, '  '
!        call minmax_uv(q0(:,:,:,4),in_e,jn_e,kn_e)

    call dyn_put ( outfile, nymd, nhms, prec, w_f, rc,  &
                   nstep=nstepec,verbose=.true.,vectype=vectype)

    if ( rc /= 0 )  then
      print *,' dyn_put:rc ',rc
      call die (myname, 'something wrong with dyn_put ')
    endif


    if(associated(ak55)) deallocate(ak55)
    if(associated(bk55)) deallocate(bk55)
    if(associated(akec)) deallocate(akec)
    if(associated(bkec)) deallocate(bkec)
    if(associated(bkec)) deallocate(pe1)
    if(associated(pke1)) deallocate(pke1)
    if(associated(phi1)) deallocate(phi1)
    if(associated(thv1)) deallocate(thv1)
    if(associated(pk))   deallocate(pk)
    if(associated(ps1))   deallocate(ps1)

    if(associated(pe0))   deallocate(pe0)
    if(associated(pk0))   deallocate(pk0)
    if(associated(dp0))   deallocate(dp0)
    if(associated(pke0))  deallocate(pke0)
    if(associated(u0))    deallocate(u0)
    if(associated(v0))    deallocate(v0)
    if(associated(thv0))  deallocate(thv0)
    if(associated(q0))    deallocate(q0)

    if(associated(pe2))   deallocate(pe2)
!   if(associated(dp2))   deallocate(dp2)
    if(associated(pke2))  deallocate(pke2)
    if(associated(thv2))  deallocate(thv2)


!  All done
!  --------
   call exit(0)

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Init_ --- Parses command line 
! 
! !INTERFACE:
!

   subroutine Init_ (inFile,phisFile,anaFile,uaFile,outFile,ndate_ec,ntime_ec,wsh,blend,vectype,dgrid)

!
! !USES:
!
   Implicit NONE

!
! !INPUT PARAMETERS: 
!

!
! !OUTPUT PARAMETERS:
!

      character(len=*), intent(out) :: inFile       !  Input file names
      character(len=*), intent(out) :: phisFile     !  Input phis file name
      character(len=*), intent(out) :: uaFile       !  Output file name 
      character(len=*), intent(out) :: anaFile       !  Output file name 
      character(len=*), intent(out) :: outFile      !  Output file name 
      integer  :: freq                              !  increment hours specified
                                                    !  command line


      integer               :: ndate_ec,ntime_ec,kn_e
      

! !DESCRIPTION: This routine initializes and parses the command.
!
! !REVISION HISTORY: 
!
! Jan 2006 Ravi C. Govindaraju Initial design and prologue.
!
!EOP
!-------------------------------------------------------------------------

   integer             iarg, argc
   character(len=2048)  argv

   character(len=255)   rcfile, label

   integer, parameter :: mKm = 256  ! max no of levels

   integer i, j, n,  rc, ios
   real    xWest, p
   logical :: debug = .false.
   logical, intent(out) :: wsh,blend
 

   integer,       intent(out) :: vectype ! GEOS-4 or GEOS-5 dyn vect
   logical,       intent(out) :: dgrid   ! GEOS-4 or GEOS-5 switch for winds grid
   character(len=10) nLx, nLy

print *
print *, "-------------------------------------------------------------------"
print *, "GFIO_ecread - Reading a ECMWF ERA40 GFIO file.                      "
print *, "-------------------------------------------------------------------"
print *

   argc = command_argument_count()
   if ( argc < 1 ) call usage_()

!  Defaults
!  --------
   inFile  = '/output/ravi/era40/Y2000/M06/ecmwf.era40.r1.25x1.2002010100z.nc4'
   phisFile = '/output/ravi/era40/ecmwf.era40.phis.r1.25x1.nc'
   uaFile  = ' '
   anaFile = '/output/dao_ops/GEOS-4.0.3/c403_cer_01/ana/Y2001/M12/c403_cer_01.ana.eta.20020101_00z.nc4'
   outFile = 'GFIO_ecread.eta.nc4'

   freq       = 99999999
   ndate_ec   = 999999
   ntime_ec   = 999999
   kn_e       = 72
   p1 = 0.1
   p2 = 1.0
   wsh   = .false.
   blend = .false.


   vectype = 4          ! default: assume vector is GEOS-4-type
   dgrid   = .true.     ! default: in GEOS-4 dyn-vector winds are on D-grid

   iarg = 0
   do i = 1, 32767
      iarg = iarg + 1
      if ( iarg .gt. argc ) exit
      call GetArg ( iArg, argv )
      if(index(argv,'-o') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, outFile )


      else if(index(argv,'-phis') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, phisFile )
      else if(index(argv,'-blend') .gt. 0 ) then
         blend = .true.
      else if(index(argv,'-ws') .gt. 0 ) then
         wsh = .true.


      else if(index(argv,'-g5') .gt. 0 ) then
         vectype = 5
         dgrid = .false.
      else if(index(argv,'-anafile') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, anaFile )
      else if(index(argv,'-uafile') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, uaFile )
      else if(index(argv,'-inc') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, argv )
         read(argv,*) freq
!     else if(index(argv,'-nlev') .gt. 0 ) then
!        if ( iarg+1 .gt. argc ) call usage_()
!        iarg = iarg + 1
!        call GetArg ( iArg, argv )
!        read(argv,*) kn_e
      else if(index(argv,'-date') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, argv )
         read(argv,*) ndate_ec
      else if(index(argv,'-pa') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, argv )
         read(argv,*) p1
      else if(index(argv,'-pb') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, argv )
         read(argv,*) p2
      else if(index(argv,'-time') .gt. 0 ) then
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iArg, argv )
         read(argv,*) ntime_ec
      else
         call GetArg ( iArg, inFile )
      end if
   end do

      print *, "               ", trim(inFile)
   end subroutine Init_


    subroutine Usage_()
   
print *, "NAME"
print *, "   ec2fv  Reads the ECMWF ERA40 hdf file. "
print *
print *, "SYNOPYSIS"
print *, "   ec2fv [options] input_fname"
print *
print *, "OPTIONS"
print *
print *, "  -o         ofname    output file name (default: GFIO_ecread.nc4)"
print *, "                      (default: all variables in the input file)"
print *, "  -phis      phisFile  ECMWF phis file name "
print *, "                       (default: ecmwf.era40.phis.r1.25x1.nc)"
print *, "  -pa        pabove    Pressure above (default 1.0)"
print *, "  -pb        pbelow    Pressure below (default 0.1)"
print *, "                       (default: ecmwf.era40.phis.r1.25x1.nc)"
print *, "  -blend     blend     Blending performance (default no blending.) "
print *, "  -ws        sheer     Blend the U and V using wind sheer method. "
print *, "                         (default: Traditional way.  "
print *, "                                   becomes obselete when option "
print *, "                                   -blend is not set.           "
print *, "  -g5                  specify when using GEOS-5 vector (not all transforms work)"
print *, '                       (default: GEOS-4 vector-type)'
print *
print *, "  -anafile   etafile  GMAO eta file name "
print *, "                      (default: c403_cer_01.ana.eta.20011216_00z.nc4)"
print *, "  -uafile    etafile  GMAO diag_dyn eta file name "
print *, "                      (default: c403_cer_01.ana.eta.20011216_00z.nc4)"
print *, "  -date      ndate    Date to be intialized for the monthly mean "
print *, "                      (if not given will be computed)."
print *, "  -time      ntime    Time to be intialized for the monthly mean "
print *, "                       (if not given will be computed)."
print *
print *, "DESCRIPTION"
print *, "   ec2fv  Reads the ECMWF ERA40 hdf file. "
print *
print *, "   ex:  ec2fv.x -o outfile input_file"
print * 

    call die ( myname, 'exiting' )

    end subroutine Usage_

!
      subroutine die(myname,string)
      character(len=*) myname
      character(len=*) string
!
      print *, ' --------------------------------'
      print *, '        ',myname
      print *, '   ',string
      print *, ' --------------------------------'
      call exit(1)
      return
      end subroutine die
! ---------------------------------------------------------------
!
!BOP
!
! !ROUTINE:  minmaxpe --- Compute min and max for a given variable.
!
! !INTERFACE: minmaxpe
!
      subroutine minmaxpe(pe3d_m,im_e,jm_e,k)
!
! !USES:
!
       Implicit NONE
!
! !INPUT PARAMETERS:

      real  :: pe3d_m(:,:)
      integer :: im_e,jm_e,km_e

! !OUTPUT PARAMETERS: NONE

! WORK AREA PAMETERS:
             
      real  :: pemin,pemax
      integer :: i,j,k
      integer :: imin_loc,jmin_loc,imax_loc,jmax_loc

! !DESCRIPTION: This routine Computes the min and max for a given variable.
!
! !REVISION HISTORY: 
!
! Jan 2006 Ravi C. Govindaraju Initial design and prologue.
!
!EOP
! ---------------------------------------------------------------

       pemin = 10.e+20
       pemax = -10.e+20
       do  j = 1,jm_e
        do i = 1,im_e
         if(pe3d_m(i,j) < pemin) then
            pemin = pe3d_m(i,j)
            imin_loc = i
            jmin_loc = j
         endif
         if(pe3d_m(i,j) > pemax) then
            pemax = pe3d_m(i,j)
            imax_loc = i
            jmax_loc = j
         endif
        end do
       end do
       print *,'k, min,max ',k,pemin,pemax
      end subroutine minmaxpe
end Program ec2fv_main
