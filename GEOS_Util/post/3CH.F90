      module ThreeCornerHat_mod
         use dynamics_lattice_module
         implicit none

         type triplet_member
            character*256              ::  expid
            integer                    ::  nfiles
            integer                    ::  ndates
            character*256, allocatable ::  fnames(:)
            integer                    ::  im,jm,lm
            real                       ::  undef
            real                       ::  dlon
            real                       ::  dlat
            integer,       allocatable ::  dates(:,:)
            real,          allocatable ::   levs(:)
            real,          allocatable ::  q(:,:,:,:,:)  !  LON x LAT x LEV x DATE x Field
         endtype triplet_member

      end module ThreeCornerHat_mod

! -------------------------------------------------------------------------------

      program Main
      use dynamics_lattice_module
      use ThreeCornerHat_mod
      use iso_fortran_env, only: REAL64
      implicit none

      type ( dynamics_lattice_type ) lattice

#ifdef mpi
      include 'mpif.h'
#endif

      integer  comm,myid,npes,ierror
      integer  imglobal
      integer  jmglobal
      integer  npex,npey
      logical  root

      type(triplet_member),   target, allocatable ::  experiment(:)
      type(triplet_member),           allocatable ::    dummyexp(:)
      integer nargs, nexps, m,n,k
      integer nfiles, nexp
      integer nfield
      integer nfields
      integer nymdb, nhmsb
      integer nymde, nhmse
      integer ndate
      integer ndates
      integer ndt

      real(REAL64),  allocatable :: mse12(:,:,:)     ! Mean Square Error between Fixed Corners 1 & 2
      real(REAL64),  allocatable :: var12(:,:,:)     !   Variance        between Fixed Corners 1 & 2
      real(REAL64),  allocatable :: cov12(:,:,:)     ! CoVariance        between Fixed Corners 1 & 2
      real(REAL64),  allocatable :: rho12(:,:,:)     ! Corr. Coeff.      between Fixed Corners 1 & 2

      real(REAL64),  allocatable :: msec1(:,:,:,:)   ! Mean Square Error between Corner 1 and Corner 3 
      real(REAL64),  allocatable :: msec2(:,:,:,:)   ! Mean Square Error between Corner 2 and Corner 3
      real(REAL64),  allocatable :: varc1(:,:,:,:)   !   Variance        between Corner 1 and Corner 3
      real(REAL64),  allocatable :: varc2(:,:,:,:)   !   Variance        between Corner 2 and Corner 3
      real(REAL64),  allocatable :: covc1(:,:,:,:)   ! CoVariance        between Corner 1 and Corner 3
      real(REAL64),  allocatable :: covc2(:,:,:,:)   ! CoVariance        between Corner 2 and Corner 3
      real(REAL64),  allocatable :: rhoc1(:,:,:,:)   ! Corr. Coeff.      between Corner 1 and Corner 3
      real(REAL64),  allocatable :: rhoc2(:,:,:,:)   ! Corr. Coeff.      between Corner 2 and Corner 3

      real(REAL64),  allocatable :: mse12A(:,:,:)    ! Amplitude Mean Square Error between Fixed Corners 1 & 2
      real(REAL64),  allocatable :: mse12B(:,:,:)    ! Bias      Mean Square Error between Fixed Corners 1 & 2
      real(REAL64),  allocatable :: mse12P(:,:,:)    ! Phase     Mean Square Error between Fixed Corners 1 & 2

      real(REAL64),  allocatable :: msec1A(:,:,:,:)  ! Amplitude Mean Square Error between Corner 1 and Corner 3
      real(REAL64),  allocatable :: msec2A(:,:,:,:)  ! Amplitude Mean Square Error between Corner 2 and Corner 3
      real(REAL64),  allocatable :: msec1B(:,:,:,:)  ! Bias      Mean Square Error between Corner 1 and Corner 3
      real(REAL64),  allocatable :: msec2B(:,:,:,:)  ! Bias      Mean Square Error between Corner 2 and Corner 3
      real(REAL64),  allocatable :: msec1P(:,:,:,:)  ! Phase     Mean Square Error between Corner 1 and Corner 3
      real(REAL64),  allocatable :: msec2P(:,:,:,:)  ! Phase     Mean Square Error between Corner 2 and Corner 3

      real(REAL64),  allocatable ::  varec1(:,:,:,:)  !  3CH Corner-1 Variance Error    Estimates 
      real(REAL64),  allocatable ::  varec2(:,:,:,:)  !  3CH Corner-2 Variance Error    Estimates
      real(REAL64),  allocatable ::  varec3(:,:,:,:)  !  3CH Corner-3 Variance Error    Estimates

      real(REAL64),  allocatable ::   var(:,:,:,:)   !  Variance  Variables
      real(REAL64),  allocatable ::     q(:,:,:,:)   !  Time-Mean Variables
      integer, allocatable ::    cc(:,:,:,:)   !  Time-Mean Counters

      character*256, allocatable ::  arg(:)
      character*256              :: dummy
      character*256              :: tag
      character*256, allocatable :: field(:)
      character*1                :: carg(256)
      character*4                :: darg
      character*1                :: char1
      character*4                :: charim
      character*4                :: charjm
      character*3                :: charlm
      character*8                :: charbegdate
      character*2                :: charbeghour
      character*8                :: charenddate
      character*2                :: charendhour
      equivalence ( dummy, carg )
      equivalence ( dummy, darg )

      integer im,jm,nlev
      integer nvars
      integer imax, jmax
      integer method

      real, pointer :: zlev (:)

      integer i,j,kk,L,nymd,nhms

      real     undef
      logical  defined
      logical  valid

      integer year, month

      character*20  char
      character*256 desc(5000)
      character*20  name(5000)
      character*3   months(12)
      data months /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/

! **********************************************************************
! ****                          Interfaces                          ****
! **********************************************************************
                                                                                                          
      interface
         subroutine init_levs( experiment,num,nlev,lev )
         use ThreeCornerHat_mod
         implicit none
         integer  num, nlev
         type(triplet_member) :: experiment(num)
         real, pointer :: lev(:)
         end subroutine init_levs
      end interface

! **********************************************************************
! ****                           Usage                              ****
! **********************************************************************
                                                                                                          
          nargs = command_argument_count()
      if( nargs.eq.0 ) call usage()

! **********************************************************************
! ****                     MPI Initialization                       ****
! **********************************************************************
                                                                                                          
#ifdef mpi
      call mpi_init                ( ierror ) ; comm = mpi_comm_world
      call mpi_comm_rank ( comm,myid,ierror )
      call mpi_comm_size ( comm,npes,ierror )
      npex = nint ( sqrt( float(npes) ) )
      npey = npex
      do while ( npex*npey .ne. npes )
         npex = npex-1
         npey = nint ( float(npes)/float(npex) )
      enddo
#else
      comm = 0
      npes = 1
      npex = 1
      npey = 1
      myid = 0
      ierror = 0
#endif
      root  = myid.eq.0
      call create_dynamics_lattice ( lattice,npex,npey )
      call   init_dynamics_lattice ( lattice,comm,1,1,1 )

! **********************************************************************
! ****             Input 3CH Experiment Files and Parameters        ****
! **********************************************************************
                                                                                                          
      nlev   = -999
      imax   = -999
      jmax   = -999
      ndt    = 21600
      nymdb  = -999
      nymde  = -999
      nhmsb  = -999
      nhmse  = -999
      method =  1
      tag    = 'NULL'

                      nfields = 6
      allocate( field(nfields) )
                field(1) = 'u'
                field(2) = 'v'
                field(3) = 't'
                field(4) = 'q'
                field(5) = 'h'
                field(6) = 'p'

      allocate ( arg(nargs) )
      do n=1,nargs
      call get_command_argument(n,arg(n))
      enddo

      nexps = 0
      do n=1,nargs
             dummy = trim( arg(n) )
                                                                                                          
            if( darg.eq.'-exp' ) then
                 if( nexps .eq. 0 ) then
                     nexps = nexps + 1
                     allocate( experiment(nexps) )
                 else
                     allocate( dummyexp(nexps) )
                               dummyexp = experiment
                   deallocate( experiment )
                     nexps = nexps + 1
                     allocate( experiment(nexps) )
                               experiment(1:nexps-1) = dummyexp
                   deallocate( dummyexp )
                 endif
                 read( carg(5:6),* ) nexp
                 if( nexp.ne.nexps ) then
                     print *, 'You must enter experiments in order!'
                     stop
                 endif

                 k = n+1
                 experiment(nexps)%expid = arg(k)

                 nfiles = 1
                 read(arg(k+nfiles),fmt='(a1)') char1
                 do while (char1.ne.'-' .and. n+nfiles.ne.nargs )
                 nfiles = nfiles+1
                 read(arg(k+nfiles),fmt='(a1)') char1
                 enddo
                 if( char1.eq.'-' ) nfiles = nfiles-1

                 experiment(nexps)%nfiles = nfiles
                 allocate ( experiment(nexps)%fnames(nfiles) )

                 if(root) print *, 'EXP: ',nexps
                 do m=1,nfiles
                 experiment(nexps)%fnames(m) = arg(k+m)
                 if(root) write(*,'("   FILES: ",a)' ) trim(experiment(nexps)%fnames(m))
                 enddo
                 if(root) print *
            endif

            if( trim(arg(n)).eq.'-levs' ) then
                 nlev = 1
                 read(arg(n+nlev),fmt='(a1)') char1
                 do while (char1.ne.'-' .and. n+nlev.ne.nargs )
                 nlev = nlev+1
                 read(arg(n+nlev),fmt='(a1)') char1
                 enddo
                 if( char1.eq.'-' ) nlev = nlev-1
                 allocate ( zlev(nlev) )
                 do m=1,nlev
                 read(arg(n+m),*) zlev(m)
                 enddo
            endif

            if( trim(arg(n)).eq.'-im'      ) read(arg(n+1),*) imax
            if( trim(arg(n)).eq.'-jm'      ) read(arg(n+1),*) jmax
            if( trim(arg(n)).eq.'-nymdb'   ) read(arg(n+1),*) nymdb
            if( trim(arg(n)).eq.'-nymde'   ) read(arg(n+1),*) nymde
            if( trim(arg(n)).eq.'-nhmsb'   ) read(arg(n+1),*) nhmsb
            if( trim(arg(n)).eq.'-nhmse'   ) read(arg(n+1),*) nhmse
            if( trim(arg(n)).eq.'-ndt'     ) read(arg(n+1),*) ndt
            if( trim(arg(n)).eq.'-tag'     ) read(arg(n+1),*) tag
            if( trim(arg(n)).eq.'-method'  ) read(arg(n+1),*) method

            if( trim(arg(n)).eq.'-fields' ) then
                 nfields = 1
                 read(arg(n+nfields),fmt='(a1)') char1
                 do while (char1.ne.'-' .and. n+nfields.ne.nargs )
                 nfields = nfields+1
                 read(arg(n+nfields),fmt='(a1)') char1
                 enddo
                 if( char1.eq.'-' ) nfields = nfields-1
                 deallocate ( field )
                   allocate ( field(nfields) )
                 do m=1,nfields
                 read(arg(n+m),*) field(m)
                 enddo
            endif

      enddo ! End NARGS Do-Loop

      if( nymdb.eq.-999 .or. nymde.eq.-999 .or. &
          nhmsb.eq.-999 .or. nhmse.eq.-999 ) then
          if(root) print *, 'You must supply the beginning and ending date and time'
          call my_finalize
          stop
      endif

      if( ( method.ne.1 )  .and. &
          ( method.ne.2 )  .and. &
          ( method.ne.3 ) ) then
          if(root) print *, 'You must supply the method for horizontal interpolation'
          call my_finalize
          stop
      endif

! **********************************************************************
! ****        Initialize Dates Associated with Each Experiment      ****
! **********************************************************************

      call timebeg('main        ')
      imglobal = 0
      jmglobal = 0
      if(root) then
      print *
      print *, '3CH Triplet Calculations'
      print *, '------------------------'
      endif
      call timebeg('init_triplet')
      do n=1,nexps
         call init_triplet( experiment(n),lattice )
         if( (imglobal*jmglobal.eq.0 ) .or. (experiment(n)%im*experiment(n)%jm.lt.imglobal*jmglobal) ) then
              imglobal = experiment(n)%im
              jmglobal = experiment(n)%jm
         endif
         if(root) then
         write(*,  '(" Initializing Dates for Experiment #: ",i2,3x,"Native Resolution: (",i4,",",i4,",",i3,")  EXPID: ",a)' &
                 ) n,experiment(n)%im,experiment(n)%jm,experiment(n)%lm,trim(experiment(n)%expid)
         endif
      enddo
      call timeend('init_triplet')
      if(root) print *

      if( imax.ne.-999 ) imglobal = imax
      if( jmax.ne.-999 ) jmglobal = jmax

!     Levels Initialization
!     ---------------------
      if( nlev.eq.-999 ) then
          call timebeg('init_levs   ')
          call init_levs( experiment,nexps,nlev,zlev )
          call timeend('init_levs   ')
      endif
      if( root ) then
          write(*,'(" Output RLSV: (",i4,",",i4,",",i3,")")') imglobal,jmglobal,nlev
          print *
          print *, 'Verifying Levels:'
          print *, '----------------'
          write(6,*) (zlev(n),n=1,nlev)
          print *
          write(charim,'(i4.4)') imglobal
          write(charjm,'(i4.4)') jmglobal
          write(charlm,'(i3.3)') nlev
          write(charbegdate,'(i8.8)') nymdb
          write(charbeghour,'(i2.2)') nhmsb/10000
          write(charenddate,'(i8.8)') nymde
          write(charendhour,'(i2.2)') nhmse/10000
      endif

! **********************************************************************
! ****                Create Lattice and Experiments                ****
! **********************************************************************

      if( root ) write(*,'(" Creating ",i4," x ",i4," lattice ...")' ) imglobal,jmglobal
      call destroy_dynamics_lattice ( lattice )
      call  create_dynamics_lattice ( lattice,npex,npey )
      call    init_dynamics_lattice ( lattice,comm,imglobal,jmglobal,1 )
      if( root ) print *

      im = lattice%im( lattice%pei )
      jm = lattice%jm( lattice%pej )

           nymd = nymdb
           nhms = nhmsb
         ndates = 0
      do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
         ndates = ndates + 1
         call tick( nymd,nhms,ndt )
      enddo

      do n=1,nexps
         allocate( experiment(n)%q(im,jm,nlev,ndates,nfields) )
      enddo

! **********************************************************************
! ****                           Loop over Times                    ****
! **********************************************************************

      undef = 1.0e15
       nymd = nymdb
       nhms = nhmsb
      ndate = 0

      call timebeg('read triplet')
      do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
         if( root ) write(*,'(" Reading Triplets for: ",i8.8,1x,i6.6)' ) nymd,nhms
         ndate = ndate + 1

         do k=1,nexps
            call read_triplet( experiment(k),nymd,nhms,nlev,zlev,ndate,k,field,nfields,undef,lattice,method )
         enddo
         call tick( nymd,nhms,ndt )

         if( root ) print *
      enddo
      call timeend('read triplet')
      if( root ) print *

! **********************************************************************
! ****                 We have everything we need                   ****
! **********************************************************************

      allocate( var(im,jm,nlev,nexps) )
      allocate(   q(im,jm,nlev,nexps) )
      allocate(  cc(im,jm,nlev,nexps) )

      allocate( mse12(im,jm,nlev) )
      allocate( var12(im,jm,nlev) )
      allocate( cov12(im,jm,nlev) )
      allocate( rho12(im,jm,nlev) )

      allocate( msec1(im,jm,nlev,nexps-2) )
      allocate( msec2(im,jm,nlev,nexps-2) )
      allocate( varc1(im,jm,nlev,nexps-2) )
      allocate( varc2(im,jm,nlev,nexps-2) )
      allocate( covc1(im,jm,nlev,nexps-2) )
      allocate( covc2(im,jm,nlev,nexps-2) )
      allocate( rhoc1(im,jm,nlev,nexps-2) )
      allocate( rhoc2(im,jm,nlev,nexps-2) )

      allocate( mse12A(im,jm,nlev) )
      allocate( mse12B(im,jm,nlev) )
      allocate( mse12P(im,jm,nlev) )

      allocate( msec1A(im,jm,nlev,nexps-2) )
      allocate( msec2A(im,jm,nlev,nexps-2) )
      allocate( msec1B(im,jm,nlev,nexps-2) )
      allocate( msec2B(im,jm,nlev,nexps-2) )
      allocate( msec1P(im,jm,nlev,nexps-2) )
      allocate( msec2P(im,jm,nlev,nexps-2) )

      allocate( varec1(im,jm,nlev,nexps-2) )
      allocate( varec2(im,jm,nlev,nexps-2) )
      allocate( varec3(im,jm,nlev,nexps-2) )

! **********************************************************************
! ****           Open Grads Dataset and Loop over Fields            ****
! **********************************************************************

! Open Grads File
! ---------------
      if( root ) then
          write(*,'(" Writing Data")' )
          if( trim(tag).ne.'NULL' ) then
              open (51,file=trim(tag) // '.data',form='unformatted',access='sequential')
          else
              if( method.eq.1 ) tag = 'BL.' // trim( experiment(1)%expid )
              if( method.eq.2 ) tag = 'BC.' // trim( experiment(1)%expid )
              if( method.eq.3 ) tag = 'BA.' // trim( experiment(1)%expid )
              do n=2,nexps
              tag = trim(tag) // '.' // trim( experiment(n)%expid )
              enddo
              tag = '3CH_' // charim // 'x' // charjm // '_L' // charlm // '_' // charbegdate // '_' // charbeghour // 'z_' &
                    // charenddate // '_' // charendhour // 'z_' // trim(tag)
              open (51,file=trim(tag)//'.data',form='unformatted',access='sequential')
          endif
      endif

    ! Initialize Number of Grads Variables
    ! ------------------------------------
      nvars = 0

      do nfield = 1,nfields

       q = 0.0
      cc = 0

! **********************************************************************
! ****    Compute Time Averages for each Field in each Experiment   ****
! **********************************************************************

    ! Compute Time Averages
    ! ---------------------
      if( root ) write(*,'(/," Computing Time Averages for Variable: ",a)' ) trim(field(nfield))
      call timebeg('time mean variables')
      nymd = nymdb
      nhms = nhmsb
         n = 0
      do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
         n = n + 1

         do k=1,nexps
            do L=1,nlev
            do j=1,jm
            do i=1,im
                  valid = .true.
               do kk=1,nexps
                  valid = valid .and. defined( experiment(kk)%q(i,j,L,n,nfield),undef )
               enddo
               if( valid ) then
                    q(i,j,L,k) =  q(i,j,L,k) + experiment(k)%q(i,j,L,n,nfield)
                   cc(i,j,L,k) = cc(i,j,L,k) + 1
               endif
            enddo
            enddo
            enddo
         enddo

         call tick( nymd,nhms,ndt )
      enddo

      do k=1,nexps
         do L=1,nlev
         do j=1,jm
         do i=1,im
            if( cc(i,j,L,k).ne.0 ) then
                 q(i,j,L,k) =  q(i,j,L,k) / cc(i,j,L,k)
            else
                 q(i,j,L,k) = undef
            endif
         enddo
         enddo
         enddo
      enddo
      call timeend('time mean variables')

    ! Compute Variances
    ! -----------------
      var = 0.0
       cc = 0

      if( root ) write(*,'(" Computing Variance for Variable: ",a)' ) trim(field(nfield))
      call timebeg('variance variables')
      nymd = nymdb
      nhms = nhmsb
         n = 0
      do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
         n = n + 1

         do k=1,nexps
            do L=1,nlev
            do j=1,jm
            do i=1,im
                  valid = .true.
               do kk=1,nexps
                  valid = valid .and. defined( experiment(kk)%q(i,j,L,n,nfield),undef )
               enddo
               if( valid ) then
                  var(i,j,L,k) = var(i,j,L,k) +  ( experiment(k)%q(i,j,L,n,nfield) - q(i,j,L,k) )**2
                   cc(i,j,L,k) =  cc(i,j,L,k) + 1
               endif
            enddo
            enddo
            enddo
         enddo

         call tick( nymd,nhms,ndt )
      enddo

      do k=1,nexps
         do L=1,nlev
         do j=1,jm
         do i=1,im
            if( cc(i,j,L,k).ne.0 ) then
               var(i,j,L,k) = var(i,j,L,k) / cc(i,j,L,k)
            else
               var(i,j,L,k) = undef
            endif
         enddo
         enddo
         enddo
      enddo
      call timeend('variance variables')

! ********************************************************************************
! ****       Compute Mean Square Error and Variance between Fixed Corners     ****
! ********************************************************************************

    ! Compute Mean Square Error:  MSE(X-Y) = 1/N SUM[ X-Y ]^2
    ! -------------------------------------------------------
      if( root ) write(*,'(" Computing Mean Square Error between Corners for: ",a)' ) trim(field(nfield))
      call timebeg('X-Y variances')
      mse12 = 0.0
      cov12 = 0.0
         cc = 0.0

      nymd = nymdb
      nhms = nhmsb
         n = 0
      do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
         n = n + 1
            do L=1,nlev
            do j=1,jm
            do i=1,im
                  valid = .true.
               do kk=1,nexps
                  valid = valid .and. defined( experiment(kk)%q(i,j,L,n,nfield),undef )
               enddo
               if( valid ) then
                   mse12(i,j,L) = mse12(i,j,L) +  ( experiment(1)%q(i,j,L,n,nfield) - experiment(2)%q(i,j,L,n,nfield) )**2
                   cov12(i,j,L) = cov12(i,j,L) +  ( experiment(1)%q(i,j,L,n,nfield) - q(i,j,L,1) ) &
                                               *  ( experiment(2)%q(i,j,L,n,nfield) - q(i,j,L,2) )
                   cc(i,j,L,3) =  cc(i,j,L,3) + 1
               endif
            enddo
            enddo
            enddo
            call tick( nymd,nhms,ndt )
      enddo

    ! Compute Variance between Fixed Corners:  SIG^2( X-Y ) = MSE( X-Y ) - (Xave-Yave)^2
    ! ----------------------------------------------------------------------------------
      do L=1,nlev
      do j=1,jm
      do i=1,im
         if( cc(i,j,L,3).ne.0 ) then
             mse12(i,j,L) = mse12(i,j,L) /  cc(i,j,L,3)
             cov12(i,j,L) = cov12(i,j,L) /  cc(i,j,L,3)
             var12(i,j,L) = mse12(i,j,L) - ( q(i,j,L,1)-q(i,j,L,2) )**2
             if( ( var(i,j,L,1).eq.0 .or. var(i,j,L,2).eq.0 ) .and. cov12(i,j,L).eq.0 ) then
               rho12(i,j,L) = 1
             else
               rho12(i,j,L) = cov12(i,j,L) / sqrt( var(i,j,L,1) ) / sqrt( var(i,j,L,2) )
             endif

             mse12B(i,j,L) = ( q(i,j,L,1)-q(i,j,L,2) )**2
             mse12A(i,j,L) = ( sqrt( var(i,j,L,1) ) - sqrt( var(i,j,L,2) ) )**2
             mse12P(i,j,L) = 2 * ( 1-rho12(i,j,L) )*( sqrt( var(i,j,L,1) ) * sqrt( var(i,j,L,2) ) )
         else
             mse12(i,j,L) = undef
             cov12(i,j,L) = undef
             var12(i,j,L) = undef
             rho12(i,j,L) = undef

             mse12B(i,j,L) = undef
             mse12A(i,j,L) = undef
             mse12P(i,j,L) = undef
         endif
      enddo
      enddo
      enddo

! -------------------------------------------------------------------------------------------

      do m=1,nexps-2
         k=m+2
         msec1(:,:,:,m) = 0.0
         msec2(:,:,:,m) = 0.0
         covc1(:,:,:,m) = 0.0
         covc2(:,:,:,m) = 0.0

         cc = 0.0
         nymd = nymdb
         nhms = nhmsb
            n = 0
         do while( (nymd.lt.nymde) .or. (nymd.eq.nymde .and. nhms.le.nhmse) )
            n = n + 1
            do L=1,nlev
            do j=1,jm
            do i=1,im
                  valid = .true.
               do kk=1,nexps
                  valid = valid .and. defined( experiment(kk)%q(i,j,L,n,nfield),undef )
               enddo
               if( valid ) then
                   msec1(i,j,L,m) = msec1(i,j,L,m) +  ( experiment(1)%q(i,j,L,n,nfield) - experiment(k)%q(i,j,L,n,nfield) )**2
                   covc1(i,j,L,m) = covc1(i,j,L,m) +  ( experiment(1)%q(i,j,L,n,nfield) - q(i,j,L,1) )*( experiment(k)%q(i,j,L,n,nfield) - q(i,j,L,k) )
                      cc(i,j,L,2) =    cc(i,j,L,2) + 1
 
                   msec2(i,j,L,m) = msec2(i,j,L,m) +  ( experiment(2)%q(i,j,L,n,nfield) - experiment(k)%q(i,j,L,n,nfield) )**2
                   covc2(i,j,L,m) = covc2(i,j,L,m) +  ( experiment(2)%q(i,j,L,n,nfield) - q(i,j,L,2) )*( experiment(k)%q(i,j,L,n,nfield) - q(i,j,L,k) )
                      cc(i,j,L,1) =    cc(i,j,L,1) + 1
               endif
            enddo
            enddo
            enddo
            call tick( nymd,nhms,ndt )
         enddo

            do L=1,nlev
            do j=1,jm
            do i=1,im
               if( cc(i,j,L,2).ne.0 ) then
                   msec1(i,j,L,m) = msec1(i,j,L,m) / cc(i,j,L,2)
                   covc1(i,j,L,m) = covc1(i,j,L,m) / cc(i,j,L,2)
                   varc1(i,j,L,m) = msec1(i,j,L,m) - ( q(i,j,L,1)-q(i,j,L,k) )**2
                   if( ( var(i,j,L,1).eq.0 .or. var(i,j,L,k).eq.0 ) .and. covc1(i,j,L,m).eq.0 ) then
                     rhoc1(i,j,L,m) = 1
                   else
                     rhoc1(i,j,L,m) = covc1(i,j,L,m) / sqrt( var(i,j,L,1) ) / sqrt( var(i,j,L,k) )
                   endif

                   msec1B(i,j,L,m) = ( q(i,j,L,1)-q(i,j,L,k) )**2
                   msec1A(i,j,L,m) = ( sqrt( var(i,j,L,1) ) - sqrt( var(i,j,L,k) ) )**2
                   msec1P(i,j,L,m) = 2 * ( 1-rhoc1(i,j,L,m) )*( sqrt( var(i,j,L,1) ) * sqrt( var(i,j,L,k) ) )
               else
                   msec1(i,j,L,m) = undef
                   covc1(i,j,L,m) = undef
                   varc1(i,j,L,m) = undef
                   rhoc1(i,j,L,m) = undef

                   msec1B(i,j,L,m) = undef
                   msec1A(i,j,L,m) = undef
                   msec1P(i,j,L,m) = undef
               endif

               if( cc(i,j,L,1).ne.0 ) then
                   msec2(i,j,L,m) = msec2(i,j,L,m) / cc(i,j,L,1)
                   covc2(i,j,L,m) = covc2(i,j,L,m) / cc(i,j,L,1)
                   varc2(i,j,L,m) = msec2(i,j,L,m) - ( q(i,j,L,2)-q(i,j,L,k) )**2
                   if( ( var(i,j,L,2).eq.0 .or. var(i,j,L,k).eq.0 ) .and. covc2(i,j,L,m).eq.0 ) then
                     rhoc2(i,j,L,m) = 1
                   else
                     rhoc2(i,j,L,m) = covc2(i,j,L,m) / sqrt( var(i,j,L,2) ) / sqrt( var(i,j,L,k) )
                   endif

                   msec2B(i,j,L,m) = ( q(i,j,L,2)-q(i,j,L,k) )**2
                   msec2A(i,j,L,m) = ( sqrt( var(i,j,L,2) ) - sqrt( var(i,j,L,k) ) )**2
                   msec2P(i,j,L,m) = 2 * ( 1-rhoc2(i,j,L,m) )*( sqrt( var(i,j,L,2) ) * sqrt( var(i,j,L,k) ) )
               else
                   msec2(i,j,L,m) = undef
                   covc2(i,j,L,m) = undef
                   varc2(i,j,L,m) = undef
                   rhoc2(i,j,L,m) = undef

                   msec2B(i,j,L,m) = undef
                   msec2A(i,j,L,m) = undef
                   msec2P(i,j,L,m) = undef
               endif
            enddo
            enddo
            enddo
      enddo
      call timeend('X-Y variances')

! -------------------------------------------------------------------------------------------

    ! Compute Initial 3CH Error Variances Estimates
    ! ---------------------------------------------
          if( root ) write(*,'(" Computing Initial 3CH Variance Error Estimates for Variable: ",a)' ) trim(field(nfield))
          call timebeg('err estimate')
          do n=1,nexps-2

            do L=1,nlev
            do j=1,jm
            do i=1,im
               if( defined( real(var12(i,j,L)  ,kind=4),undef ) .and. &
                   defined( real(varc1(i,j,L,n),kind=4),undef ) .and. &
                   defined( real(varc2(i,j,L,n),kind=4),undef ) ) then
                   varec1(i,j,L,n) = 0.5*( var12(i,j,L)   + varc1(i,j,L,n) - varc2(i,j,L,n) )
                   varec2(i,j,L,n) = 0.5*( var12(i,j,L)   + varc2(i,j,L,n) - varc1(i,j,L,n) )
                   varec3(i,j,L,n) = 0.5*( varc1(i,j,L,n) + varc2(i,j,L,n) - var12(i,j,L)   )
               else
                   varec1(i,j,L,n) = real(undef,kind=4)
                   varec2(i,j,L,n) = real(undef,kind=4)
                   varec3(i,j,L,n) = real(undef,kind=4)
               endif
            enddo
            enddo
            enddo

          enddo
          call timeend('err estimate')

! Write Grads Data
! ----------------
          call timebeg('grads write ')

          ! Write 3CH VAR Error Estimates
          ! -----------------------------
            do n=1,nexps-2
               m=n+2
               write(char,'("vare1"i1)') m
               call writit ( varec1(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = '3CH_Variance_Error_Estimate_for_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_using_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("vare2"i1)') m
               call writit ( varec2(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = '3CH_Variance_Error_Estimate_for_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_using_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("vare"i1)') m
               call writit ( varec3(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = '3CH_Variance_Error_Estimate_for_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write VAR
          ! ---------
               call writit ( var12,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'var12  '
               desc(nvars) = 'Variance_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("var1"i1)') m
               call writit ( varc1(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Variance_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("var2"i1)') m
               call writit ( varc2(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Variance_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write MSE Total
          ! ---------------
               call writit ( mse12,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'mse12  '
               desc(nvars) = 'Total_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("mse1"i1)') m
               call writit ( msec1(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Total_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("mse2"i1)') m
               call writit ( msec2(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Total_Mean_Square_Error_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write MSE Bias
          ! --------------
               call writit ( mse12B,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'mse12B '
               desc(nvars) = 'BIAS_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("mse1"i1"B")') m
               call writit ( msec1b(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'BIAS_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("mse2"i1"B")') m
               call writit ( msec2b(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'BIAS_Mean_Square_Error_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write MSE Amplitude
          ! -------------------
               call writit ( mse12A,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'mse12A '
               desc(nvars) = 'AMPL_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("mse1"i1"A")') m
               call writit ( msec1a(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'AMPL_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("mse2"i1"A")') m
               call writit ( msec2a(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'AMPL_Mean_Square_Error_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write MSE Phase
          ! ---------------
               call writit ( mse12P,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'mse12P '
               desc(nvars) = 'PHAZ_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("mse1"i1"P")') m
               call writit ( msec1p(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'PHAZ_Mean_Square_Error_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("mse2"i1"P")') m
               call writit ( msec2p(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'PHAZ_Mean_Square_Error_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write RHO
          ! ---------
               call writit ( rho12,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // 'rho12  '
               desc(nvars) = 'Correlation_Coefficient_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Fixed_Corner_2_' // trim(experiment(2)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("rho1"i1)') m
               call writit ( rhoc1(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Correlation_Coefficient_between_Fixed_Corner_1_' // trim(experiment(1)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)

               write(char,'("rho2"i1)') m
               call writit ( rhoc2(1,1,1,n) ,im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Correlation_Coefficient_between_Fixed_Corner_2_' // trim(experiment(2)%expid) // '_and_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write Time Means
          ! ----------------
                                  n = 1
               write(char,'(i1)') n
               call writit ( q(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Time_Mean_for_Fixed_Corner_1_' // trim(experiment(n)%expid)
                                  n = 2
               write(char,'(i1)') n
               call writit ( q(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Time_Mean_for_Fixed_Corner_2_' // trim(experiment(n)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'(i1)') m
               call writit ( q(1,1,1,m),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Time_Mean_for_Third_Corner_' // trim(experiment(m)%expid)
            enddo

          ! Write Variances
          ! ---------------
                                       n = 1
               write(char,'("var"i1)') n
               call writit ( var(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Variance_for_Fixed_Corner_1_' // trim(experiment(n)%expid)
                                       n = 2
               write(char,'("var"i1)') n
               call writit ( var(1,1,1,n),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Variance_for_Fixed_Corner_2_' // trim(experiment(n)%expid)
            do n=1,nexps-2
               m=n+2
               write(char,'("var"i1)') m
               call writit ( var(1,1,1,m),im,jm,nlev,lattice ) ; nvars = nvars + 1 ; name(nvars) = trim(field(nfield)) // trim(char)
               desc(nvars) = 'Variance_for_Third_Corner_' // trim(experiment(m)%expid)
            enddo
            call timeend('grads write ')

! -------------------------------------------------------------------------------------------

      enddo ! ENDDO for NFIELDS Loop

! Write Grads CTL file
! --------------------
      if( root ) then
          open (86,file= trim(tag) // '.ctl',form='formatted',access='sequential')
          write(86,5001) trim(tag) // '.data'
          write(86,5002)
          write(86,5003) undef
          write(86,5004) imglobal,360.0/float(imglobal)
          write(86,5005) jmglobal,180.0/float(jmglobal-1)
          write(86,5006) nlev,(zlev(k),k=1,nlev)

           nymd = nymdb
           nhms = nhmsb
           year =     nymd/10000
          month = mod(nymd,10000)/100
          write(86,5007) ndates/ndates, months(month),year

          write(86,5008) nvars
          do n=1,nvars
             write(86,5009) name(n),nlev,trim(desc(n))
          enddo
          write(86,5010)
      endif

 5001 format('dset    ^',a)
 5002 format('title   3CH_Data',/,'options sequential')
 5003 format('undef  ',g13.6)
 5004 format('xdef  ',i4,' linear  -180 ',f8.5)
 5005 format('ydef  ',i4,' linear   -90 ',f8.5)
 5006 format('zdef  ',i3,' levels  ',5(f7.2,1x),/, &
          50(18x,                   5(f7.2,1x),/) )
 5007 format('tdef  ',i3,' linear  ',   '00:00',   'Z01'   ,a3,i4,'  6hr')
 5008 format('vars  ',i3)
 5009 format(a,2x,i3,' 0 ',a)
 5010 format('endvars')

      call timeend('main        ')
      if( root ) call timepri (6)

      call my_finalize
      stop
      end

      subroutine init_levs( experiment,num,nl,lev )
      use ThreeCornerHat_mod
      implicit none
      integer  num, nl
      type(triplet_member) :: experiment(num)
      real, pointer        :: lev(:)
      real, allocatable    :: dum(:)
      real                   zlev
      integer                nlev

      integer L,k,n

      nl = 0
      do L=1,experiment(1)%lm
         zlev = experiment(1)%levs(L)
         nlev = 1
         do n=2,num
            do k=1,experiment(n)%lm
               if( experiment(n)%levs(k).eq.zlev ) then
                   nlev = nlev + 1
                   exit
               endif
            enddo
         enddo
         if( nlev.eq.num ) then
             if( nl.eq.0 ) then
                 nl = nl + 1
                 allocate( lev(nl) )
                 lev(nl) = zlev
             else
                 allocate( dum(nl) )
                 dum = lev
                 deallocate( lev )
                 nl = nl + 1
                 allocate( lev(nl) )
                 lev(1:nl-1) = dum
                 lev(nl)     = zlev
                 deallocate( dum )
             endif
         endif
      enddo

      return
      end

      subroutine init_triplet ( experiment,lattice )
      use dynamics_lattice_module
      use ThreeCornerHat_mod
      implicit none

      type ( dynamics_lattice_type ) lattice

#ifdef mpi
      include 'mpif.h'
#endif

      integer  comm,myid,npes,ierror
      integer  npex,npey
      logical  root

      type(triplet_member) :: experiment

      integer       id,nvars,rc
      integer       ntime,ngatts,timinc
      integer       im,jm,lm
      integer       ndates,num
      real          undef

      character*256  title
      character*256  source
      character*256  contact
      character*256  levunits
      character*256, allocatable ::  vname(:)
      character*256, allocatable :: vtitle(:)
      character*256, allocatable :: vunits(:)

      real,    allocatable ::    lat(:)
      real,    allocatable ::    lon(:)
      real,    allocatable ::    lev(:)
      real,    allocatable :: vrange(:,:)
      real,    allocatable :: prange(:,:)
      integer, allocatable :: yymmdd(:)
      integer, allocatable :: hhmmss(:)
      integer, allocatable ::  kmvar(:)
      integer, allocatable ::  dates(:,:)
      integer, allocatable ::  datez(:,:)

      root = lattice%myid.eq.0

! Read All Experiment Files to Gather Date and Time Information
! -------------------------------------------------------------
!     print *, 'Examining Triplet Files ...'
!     print *
      ndates = 0

      do num=1,experiment%nfiles

         if( root ) then
         call gfio_open       ( experiment%fnames(num),1,id,rc )
         call gfio_diminquire ( id,im,jm,lm,ntime,nvars,ngatts,rc )
         if( rc.ne.0 ) then
             print *, 'Failure to Open File: ',trim( experiment%fnames(num) )
         endif
         endif

#ifdef mpi
          call mpi_bcast (    rc,1,mpi_integer,0,lattice%comm,ierror )
          call mpi_bcast (    im,1,mpi_integer,0,lattice%comm,ierror )
          call mpi_bcast (    jm,1,mpi_integer,0,lattice%comm,ierror )
          call mpi_bcast (    lm,1,mpi_integer,0,lattice%comm,ierror )
          call mpi_bcast ( ntime,1,mpi_integer,0,lattice%comm,ierror )
          call mpi_bcast ( nvars,1,mpi_integer,0,lattice%comm,ierror )
#endif

          if( rc.ne.0 ) then
              call my_finalize
              stop
          endif

          allocate ( lon(im) )
          allocate ( lat(jm) )
          allocate ( lev(lm) )
          allocate ( yymmdd(ntime) )
          allocate ( hhmmss(ntime) )
          allocate (  vname(nvars) )
          allocate ( vtitle(nvars) )
          allocate ( vunits(nvars) )
          allocate (  kmvar(nvars) )
          allocate ( vrange(2,nvars) )
          allocate ( prange(2,nvars) )
                                                                                                          
          if( root ) then
          call gfio_inquire ( id,im,jm,lm,ntime,nvars,    &
                              title,source,contact,undef, &
                              lon,lat,lev,levunits,       &
                              yymmdd,hhmmss,timinc,       &
                              vname,vtitle,vunits,kmvar,  &
                              vrange,prange,rc )
          endif

#ifdef mpi
      call mpi_bcast ( yymmdd,ntime,   mpi_integer,0,lattice%comm,ierror )
      call mpi_bcast ( hhmmss,ntime,   mpi_integer,0,lattice%comm,ierror )
      call mpi_bcast ( timinc,1,       mpi_integer,0,lattice%comm,ierror )
      call mpi_bcast ( kmvar,nvars,    mpi_integer,0,lattice%comm,ierror )
      call mpi_bcast ( lev,lm, lattice%mpi_rkind,  0,lattice%comm,ierror )
#endif

          if( ndates.eq.0 ) then
              ndates = ndates + ntime
              allocate ( dates(3,ndates) )
              allocate ( datez(3,ndates) )
              dates(1,ndates-ntime+1:ndates) = yymmdd(:)
              dates(2,ndates-ntime+1:ndates) = hhmmss(:)
              dates(3,ndates-ntime+1:ndates) = num
              datez = dates
              allocate ( experiment%levs(lm) )
                         experiment%levs = lev
          else
              deallocate(dates)
              ndates = ndates + ntime
              allocate ( dates(3,ndates) )
              dates(1,1:ndates-ntime) = datez(1,1:ndates-ntime)
              dates(2,1:ndates-ntime) = datez(2,1:ndates-ntime)
              dates(3,1:ndates-ntime) = datez(3,1:ndates-ntime)
              dates(1,ndates-ntime+1:ndates) = yymmdd(:)
              dates(2,ndates-ntime+1:ndates) = hhmmss(:)
              dates(3,ndates-ntime+1:ndates) = num
              deallocate(datez)
              allocate ( datez(3,ndates) )
              datez = dates
          endif

          if(root) call gfio_close(id,rc)
          deallocate ( lon,lat,lev,yymmdd,hhmmss,vname,vtitle,vunits,kmvar,vrange,prange )
      enddo

          allocate ( experiment%dates(3,ndates) )

                     experiment%dates  =  dates
                     experiment%im     =  im
                     experiment%jm     =  jm
                     experiment%lm     =  lm
                     experiment%ndates =  ndates
                     experiment%undef  =  undef

      ! Note:
      ! -----
      ! NDATES are the total number of Dates defined by the INPUT filenames for each Member
      ! NDATES is NOT necessarily the total number of Dates requested by the 3CH.rc

      return
      end

      subroutine read_triplet ( experiment,nymd,nhms,lmo,zlev,ndate,nexp,field,nfields,zundef,lattice,method )
      use ThreeCornerHat_mod
      use dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice
      type(triplet_member) :: experiment
#ifdef mpi
      include 'mpif.h'
#endif
      integer       nymd,nhms,ndate,nexp
      integer       lmo,nfields
      real          undef, zundef
      real          zlev(lmo)

      integer       id,nvars,rc
      integer       ntime,ngatts,timinc
      integer       L,m,n,msgn
      integer       im,jm,lm,LN
      integer       i,imglobal
      integer       j,jmglobal
      integer       iratio,jratio
      integer       method

      integer npes, ierror
      integer index(lmo*nfields)

      character*256  title
      character*256  source
      character*256  contact
      character*256  levunits
      character*256, allocatable ::  vname(:)
      character*256, allocatable :: vtitle(:)
      character*256, allocatable :: vunits(:)

      real,    allocatable ::   qglo(:,:,:)
      real,    allocatable ::      q(:,:,:)
      real,    allocatable ::    lat(:)
      real,    allocatable ::    lon(:)
      real,    allocatable ::    lev(:)
      real,    allocatable :: vrange(:,:)
      real,    allocatable :: prange(:,:)
      integer, allocatable :: yymmdd(:)
      integer, allocatable :: hhmmss(:)
      integer, allocatable ::  kmvar(:)

! Default Aliases
! ---------------
      character*256 field(nfields)
      character*256 alias

      character*256 alias_u (6)
      character*256 alias_v (6)
      character*256 alias_t (6)
      character*256 alias_q (3)
      character*256 alias_h (4)
      character*256 alias_p (4)

      data alias_u /'u','uwnd','ugrd','ugrdprs','U_velocity','ugrdes'/
      data alias_v /'v','vwnd','vgrd','vgrdprs','V_velocity','vgrdes'/
      data alias_t /'t','tmpu','tmp' ,'tmpprs','Temperature','tmpes' /
      data alias_q /'q','sphu','qv' /
      data alias_h /'h','hght','Height','hgtes' /
      data alias_p /'ps','sp','pressfc','surface_pressur' /

      logical check_names
      logical shift
      logical found_var
      logical found_date
      integer num, LL, k, loc, len
      data id      /0/
      data num     /0/
      data shift /.false./

! Read Appropriate Triplet File to Get Data
! -----------------------------------------
      imglobal = lattice%imglobal
      jmglobal = lattice%jmglobal
      allocate ( qglo(imglobal,jmglobal,lmo) )

      npes = lattice%nx * lattice%ny
      do m=1,nfields
      do L=1,lmo
            LN  = L + (m-1)*LMO
      index(LN) = mod(LN-1,npes)
      enddo
      enddo

    ! Initialize Experiment to UNDEF
    ! ------------------------------
      experiment%q(:,:,:,ndate,:) = zundef

    ! Find the Appropriate Date
    ! -------------------------
      found_date = .false.
      do num=1,experiment%ndates
      if( experiment%dates(1,num).eq.nymd .and. experiment%dates(2,num).eq.nhms ) then
          found_date = .true.

          if(lattice%myid.eq.0) then
             write(*,'("                Found: ",i8.8,1x,i6.6," for Experiment ",i3)' ) nymd,nhms,nexp
          endif
          call gfio_open       ( experiment%fnames(experiment%dates(3,num)),1,id,rc )
          if(lattice%myid.eq.0) then
          call gfio_diminquire ( id,im,jm,lm,ntime,nvars,ngatts,rc )
          endif

#ifdef mpi
          call mpi_bcast (    im,1, mpi_integer, 0, lattice%comm,ierror )
          call mpi_bcast (    jm,1, mpi_integer, 0, lattice%comm,ierror )
          call mpi_bcast (    lm,1, mpi_integer, 0, lattice%comm,ierror )
          call mpi_bcast ( ntime,1, mpi_integer, 0, lattice%comm,ierror )
          call mpi_bcast ( nvars,1, mpi_integer, 0, lattice%comm,ierror )
#endif

          allocate ( lon(im) )
          allocate ( lat(jm) )
          allocate ( lev(lm) )
          allocate ( yymmdd(ntime) )
          allocate ( hhmmss(ntime) )
          allocate (  vname(nvars) )
          allocate ( vtitle(nvars) )
          allocate ( vunits(nvars) )
          allocate (  kmvar(nvars) )
          allocate ( vrange(2,nvars) )
          allocate ( prange(2,nvars) )
                                                                                                          
          if(lattice%myid.eq.0) then
          call gfio_inquire ( id,im,jm,lm,ntime,nvars,    &
                              title,source,contact,undef, &
                              lon,lat,lev,levunits,       &
                              yymmdd,hhmmss,timinc,       &
                              vname,vtitle,vunits,kmvar,  &
                              vrange,prange,rc )
          endif

#ifdef mpi
      call mpi_bcast ( yymmdd,ntime,    mpi_integer,   0,lattice%comm,ierror )
      call mpi_bcast ( hhmmss,ntime,    mpi_integer,   0,lattice%comm,ierror )
      call mpi_bcast ( kmvar,nvars,     mpi_integer,   0,lattice%comm,ierror )
      call mpi_bcast ( vname,nvars*256, mpi_character, 0,lattice%comm,ierror )
      call mpi_bcast ( undef,1,         mpi_real,      0,lattice%comm,ierror )
      call mpi_bcast ( lon,im,          mpi_real,      0,lattice%comm,ierror )
      call mpi_bcast ( lev,lm,          mpi_real,      0,lattice%comm,ierror )
#endif

              shift = .false.
          if( lon(1).eq.0.0 ) then
              if(lattice%myid.eq.0) write(*,'("                Triplet data begins at lon: ", &
                                              g6.2," Horizontal Shift will be performed")' ) lon(1)
              shift = .true.
          endif

          allocate ( q(im,jm,lmo) )

          do m=1,nfields      ! Loop over fields of interest
             found_var = .false.  ! Initialize FOUND flag for each field
             do n=1,nvars     ! Loop over variables in Experiment File

            ! Generic Field
            ! -------------
             if( .not.found_var ) then

                 if( trim(field(m)).eq.'u' ) then ; len = size( alias_u ) ; msgn = 1 ; endif
                 if( trim(field(m)).eq.'v' ) then ; len = size( alias_v ) ; msgn = 1 ; endif
                 if( trim(field(m)).eq.'t' ) then ; len = size( alias_t ) ; msgn = 0 ; endif
                 if( trim(field(m)).eq.'q' ) then ; len = size( alias_q ) ; msgn = 0 ; endif
                 if( trim(field(m)).eq.'h' ) then ; len = size( alias_h ) ; msgn = 0 ; endif
                 if( trim(field(m)).eq.'p' ) then ; len = size( alias_p ) ; msgn = 0 ; endif

                 do k = 1,len
                    if( trim(field(m)).eq.'u' ) alias = alias_u(k)
                    if( trim(field(m)).eq.'v' ) alias = alias_v(k)
                    if( trim(field(m)).eq.'t' ) alias = alias_t(k)
                    if( trim(field(m)).eq.'q' ) alias = alias_q(k)
                    if( trim(field(m)).eq.'h' ) alias = alias_h(k)
                    if( trim(field(m)).eq.'p' ) alias = alias_p(k)

                    if( check_names( vname(n),alias ) ) then
                        found_var = .true.
                        if(lattice%myid.eq.0) then
                           write(*,'("                         Found Variable: ",a)' ) trim(field(m))
                        endif
 
                       do L=1,lmo
                          do LL=1,lm
                          if( lev(LL).eq.zlev(L) ) loc = LL
                          enddo
                          LN = L + (m-1)*LMO

                          ! Read on ROOT: Global Q(IMxJM) from Experiment
                          ! ---------------------------------------------
                          if( index(LN).eq.lattice%myid ) then

                             call timebeg('  gfio getvar ')
                             if( trim(field(m)).eq.'p' ) then
                                 call gfio_getvar ( id,vname(n),nymd,nhms,im,jm,0  ,1,q(1,1,L),rc )
                             else
                                 call gfio_getvar ( id,vname(n),nymd,nhms,im,jm,loc,1,q(1,1,L),rc )
                             endif
                             call timeend('  gfio getvar ')

                             where( q(:,:,L).eq.undef ) ; q(:,:,L) = zundef ; endwhere
                             if( shift ) call hshift ( q(1,1,L),im,jm )

                             if( im.eq.imglobal .and. jm.eq.jmglobal ) then
                                 qglo(:,:,L) = q(:,:,L)
                             else

                               ! Bi-Linear Interpolation
                               ! -----------------------
                                 if( method.eq.1 ) then
                                     call timebeg(' bi_linear')
                                     iratio =  im   / imglobal
                                     jratio = (jm-1)/(jmglobal-1)
                                     if( (  imglobal   *iratio .eq.  im    ) .and. &
                                         ( (jmglobal-1)*jratio .eq. (jm-1) )  ) then
                                         do j=1,jmglobal
                                         do i=1,imglobal
                                            qglo(i,j,L) = q( 1+(i-1)*iratio, 1+(j-1)*jratio,L )
                                         enddo
                                         enddo
                                     else
                                         call hinterp ( q(1,1,L),im,jm,qglo(1,1,L),imglobal,jmglobal,1,zundef,method )
                                     endif
                                     call timeend(' bi_linear')
                                 endif

                               ! Bi-Cubic Interpolation
                               ! ----------------------
                                 if( method.eq.2 ) then
                                     call timebeg(' bi_cubic')
                                     iratio =  im   / imglobal
                                     jratio = (jm-1)/(jmglobal-1)
                                     if( (  imglobal   *iratio .eq.  im    ) .and. &
                                         ( (jmglobal-1)*jratio .eq. (jm-1) )  ) then
                                         do j=1,jmglobal
                                         do i=1,imglobal
                                            qglo(i,j,L) = q( 1+(i-1)*iratio, 1+(j-1)*jratio,L )
                                         enddo
                                         enddo
                                     else
                                         call hinterp ( q(1,1,L),im,jm,qglo(1,1,L),imglobal,jmglobal,1,zundef,method )
                                     endif
                                     call timeend(' bi_cubic')
                                 endif

                               ! Box-Average Interpolation
                               ! -------------------------
                                 if( method.eq.3 ) then
                                     call timebeg(' box_average')
                                     call bin( q(1,1,L),im,jm,qglo(1,1,L),imglobal,jmglobal,zundef,msgn )
                                     call timeend(' box_average')
                                 endif

                             endif  ! End im.eq.imglobal .and. jm.eq.jmglobal Test

                          endif  ! End MPI-Processor Test
                       enddo   ! End Output Level Loop
                          
                    endif ! End NAME test
                 enddo
             endif ! Endif for FOUND test
             enddo ! Enddo for N=1,NVARS loop

             if( .not.found_var ) then
                  if(lattice%myid.eq.0) then
                           write(*,'("                     Not Found Variable: ",a," Setting to UNDEF")' ) trim(field(m))
                  endif
             endif

             do L=1,lmo
#ifdef mpi
                LN = L + (m-1)*LMO
                call mpi_bcast ( qglo(1,1,L),imglobal*jmglobal,lattice%mpi_rkind,index(LN),lattice%comm,ierror )
#endif
               ! Scatter Global Output QGLO(IMGLOBALxJMGLOBAL) from Experiment onto Local Array Experiment%Q(IMxJM)
               ! --------------------------------------------------------------------------------------------------
                call timebeg('      scatter')
                call scatter_2d ( qglo(1,1,L),experiment%q(1,1,L,ndate,m),lattice )
                call timeend('      scatter')
             enddo

          enddo ! Enddo for M=1,NFIELDS loop

    ! Close File and Deallocate Space
    ! -------------------------------
      call gfio_close(id,rc)
      deallocate ( q )
      deallocate ( lon )
      deallocate ( lat )
      deallocate ( lev )
      deallocate ( yymmdd )
      deallocate ( hhmmss )
      deallocate (  vname )
      deallocate ( vtitle )
      deallocate ( vunits )
      deallocate (  kmvar )
      deallocate ( vrange )
      deallocate ( prange )

      endif  ! Endif for Proper Date Test

      enddo  ! Enddo for Proper Date Loop

      if( .not.found_date ) then
          if(lattice%myid.eq.0) then
             write(*,'("       Date Not Found: ",i8.8,1x,i6.6," for Experiment ",i3," Setting Member to UNDEF")' ) nymd,nhms,nexp
          endif
      endif

      deallocate ( qglo )

      return
      end

      subroutine bin ( qin,im_in,jm_in,qout,im_out,jm_out,undef,msgn )
      implicit none
      integer   im_in ,jm_in ,msgn
      integer   im_out,jm_out
      real      undef
      real    qin(im_in ,jm_in )
      real   qout(im_out,jm_out)
      real q10x10(360*6,180*6)
 
! Parse Arbitray Field (im,jm) to 10'x10' Variable
! ------------------------------------------------
      call timebeg('    bin_10x10')
      call bin_10x10 ( qin,im_in,jm_in,q10x10 )
      call timeend('    bin_10x10')
 
! Bin 10'x10' Variable to Output Field (im_out,jm_out)
! ----------------------------------------------------
      call timebeg('    ave_10x10')
      call averaged_10x10 ( q10x10,qout,im_out,jm_out,undef,msgn )
      call timeend('    ave_10x10')

      return
      end

      subroutine bin_10x10 ( z,im,jm,z10x10 )
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Compute a (10m X 10m) array binned from an input array (im,jm)
!
!  INPUT:
!  ======
!    z .......... Input array(im,jm)
!    im ......... Longitudinal dimension of z
!    jm ......... Latitudinal  dimension of z
!
!  OUTPUT:
!  =======
!    z10x10 ..... Output array(360*6,180*6)
!
!  NOTES:
!  ======
!    Input array z(im,jm) is assumed to be on an A-grid.
!                z(i,j)   represents the value at the center of the grid-box.
!                z(1,j)   is located at lon=-180.
!                z(i,1)   is located at lat=-90.
!                z(i,jm)  is located at lat=+90.
!
!    Output array z10x10  represents values within a 10min X 10min grid-box.
!             Each box is referenced by the latitude and longitude of
!             its southwest corner, not its center point.  Thus,
!             the height associated with a coordinate actually
!             represents the heights centered to the northeast of that point.
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      implicit none
      integer im,jm
      real  z(im,jm)
      real z10x10(360*6,180*6)

      integer i,j,ii,jj,ibeg,iend,jbeg,jend
      real    zlatc,zlonc
      real    lonbeg,lonend
      real    latbeg,latend
      real    pi,dl,dp,dz 

      pi = 4.*atan(1.)
      dl = 2*pi/im
      dp = pi/(jm-1)
      dz = pi/(6.*180)
      
      do j=1,180*6
      do i=1,360*6

      zlatc = -pi/2+(j-0.5)*dz  ! Latitude  at center of 10x10 box
      zlonc = -pi  +(i-0.5)*dz  ! Longitude at center of 10x10 box

! Find bounding lat and lon on IMxJM grid
! ---------------------------------------
      iend = nint( 1.+(zlonc+pi)/dl )
      lonend = -pi + (iend-1)*dl
      if( lonend.ge.zlonc ) then
      lonbeg = -pi + (iend-2)*dl
      else
      iend = iend+1
      lonbeg = lonend
      lonend = -pi + (iend-1)*dl
      endif
      ibeg = iend-1

      jend = nint( 1.+(zlatc+pi/2)/dp )
      latend = -pi/2 + (jend-1)*dp
      if( latend.ge.zlatc ) then
      latbeg = -pi/2 + (jend-2)*dp
      else
      jend = jend+1
      latbeg = latend
      latend = -pi/2 + (jend-1)*dp
      endif
      jbeg = jend-1


      if(iend.gt.im) iend=iend-im

      if( zlonc.le.lonbeg+0.5*dl ) then
      ii = ibeg
      else
      ii = iend
      endif
      if( zlatc.le.latbeg+0.5*dp ) then
      jj = jbeg
      else
      jj = jend
      endif

      if( ii.lt.1 .or. ii.gt.im .or. &
          jj.lt.1 .or. jj.gt.jm      ) then
          write(6,1000) i,j,zlonc*180/pi,zlatc*180/pi, &
                        ii,jj,lonbeg*180/pi,lonend*180/pi,latbeg*180/pi,latend*180/pi
      endif

      z10x10(i,j) = z(ii,jj)
      
      enddo
      enddo

 1000 format(1x,'i_10x10: ',i4,'  j_10x10: ',i4,'  lon_10x10: ',f10.4,'  lat_10x10: ',f10.4, &
             4x,'i_IMxJM: ',i4,'  j_IMxJM: ',i4,'  lonbeg: ',f10.4,'  lonend: ',f10.4, &
                                                '  latbeg: ',f10.4,'  latend: ',f10.4)
      return
      end

      subroutine averaged_10x10 ( z10x10,z,im,jm,undef,msgn )
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Average a (10m X 10m) input array to an output array (im,jm)
!
!  INPUT:
!  ======
!    z10x10 ..... Input array(360*6,180*6)
!    msgn ....... Integer Flag for scalar (0) or vector (1)
!
!  OUTPUT:
!  =======
!    z .......... Output array(im,jm)
!    im ......... Longitudinal dimension of z
!    jm ......... Latitudinal  dimension of z
!
!  NOTES:
!  ======
!    Input array z10x10  represents values within a 10min X 10min grid-box.
!             Each box is referenced by the latitude and longitude of
!             its southwest corner, not its center point.  Thus,
!             the height associated with a coordinate actually
!             represents the heights centered to the northeast of that point.
!
!    Output array z(im,jm) is assumed to be on an A-grid.
!                 z(i,j)   represents the value at the center of the grid-box.
!                 z(1,j)   is located at lon=-180.
!                 z(i,1)   is located at lat=-90.
!                 z(i,jm)  is located at lat=+90.
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      implicit none
      integer im,jm,msgn
      real  z(im,jm)
      real  dlam(im), dphi(jm)
      real  z10x10(360*6,180*6)

      integer i,j,ibeg,iend,jbeg,jend
      integer ii,jj,itmp
      real    sum1,sum2
      real    zlat,zlon
      real    lon1,lon2,wx
      real    lat1,lat2,wy
      real    lonbeg,lonend,lat,coslat
      real    latbeg,latend
      real    undef 
      real    pi,dz 
      real    lon_cmp(im)
      real    lat_cmp(jm)
      logical defined

      pi   = 4.*atan(1.)
      dz   = pi/(6.*180)
      dlam = 2*pi/ im
      dphi =   pi/(jm-1)

! Compute Computational Lambda's and Phi's
! ----------------------------------------
      lon_cmp(1) = -pi
      do i=2,im
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_cmp(1) = -pi*0.5
      do j=2,jm-1
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_cmp(jm) =  pi*0.5


! Compute average away from poles
! -------------------------------
      do j=2,jm-1
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat-dphi(j-1)/2
      latend = zlat+dphi(j)  /2
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1.+(latbeg+pi/2)/dz
      jend = 1.+(latend+pi/2)/dz

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo
      enddo

! Compute average at South Pole
! -----------------------------
         j=1
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat
      latend = zlat+dphi(j)  /2
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1
      jend = 1.+(latend+pi/2)/dz

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo

! Compute average at North Pole
! -----------------------------
         j=jm
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat-dphi(j-1)/2
      latend = zlat
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1.+(latbeg+pi/2)/dz
      jend = 1080

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo

! Average Pole Values
! -------------------
      if( msgn.eq.0 ) then
      sum1 = 0
         j = 0
      do i=1,im
         if( defined(z(i,1),undef) ) then
             sum1 = sum1 + z(i,1)
                j = j + 1
         endif
      enddo
      if( j.ne.0 ) then
      z(:,1) = sum1/j
      else
      z(:,1) = undef
      endif

      sum2 = 0
         j = 0
      do i=1,im
         if( defined(z(i,jm),undef) ) then
             sum2 = sum2 + z(i,jm)
                j = j + 1
         endif
      enddo
      if( j.ne.0 ) then
      z(:,jm) = sum2/j
      else
      z(:,jm) = undef
      endif

      endif

      return
      end

      function defined ( q,undef )
      implicit none
      logical  defined
      real     q,undef
      defined = abs(q-undef).gt.0.1*abs(undef)
      return
      end

      subroutine hshift ( q,im,jm )
      real q(im,jm), dum(im,jm)
      dum(1:im/2,:)    =   q(1:im/2,:)
        q(1:im/2,:)    =   q(1+im/2:im,:)
        q(1+im/2:im,:) = dum(1:im/2,:)
      return
      end

      subroutine minmax (q,im,jm,undef)
      real   q(im,jm)
      qmin =  1e33
      qmax = -1e33
      do j=1,jm
      do i=1,im
      if(q(i,j).ne.undef) qmin = min( qmin,q(i,j) )
      if(q(i,j).ne.undef) qmax = max( qmax,q(i,j) )
      enddo
      enddo
      print *, ' qmin: ',qmin,' qmax: ',qmax
      return
      end

      subroutine writit ( q,im,jm,lm,lattice )
      use dynamics_lattice_module
      use iso_fortran_env, only: REAL64
      implicit none
      type ( dynamics_lattice_type ) lattice
      integer  im,jm,lm,L,img,jmg
      real(REAL64) q(im,jm,lm)
      real,   allocatable :: qglo(:,:)
      real,   allocatable :: qloc(:,:)
      img = lattice%imglobal
      jmg = lattice%jmglobal
      allocate ( qglo(img,jmg) )
      allocate ( qloc(im ,jm ) )
      do L=1,lm
         qloc(:,:) = q(:,:,L)
         call timebeg('      gather ')
         call gather_2d ( qglo,qloc,lattice )
         call timeend('      gather ')
         if( lattice%myid.eq.0 ) then
             write(51) qglo
         endif
      enddo
      deallocate ( qglo )
      deallocate ( qloc )
      return
      end

      subroutine tick (nymd,nhms,ndt)
!***********************************************************************
!  Purpose
!     Tick the Date (nymd) and Time (nhms) by NDT (seconds)
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      IF(NDT.NE.0) THEN
      NSEC = NSECF(NHMS) + NDT

      IF (NSEC.GT.86400)  THEN
      DO WHILE (NSEC.GT.86400)
      NSEC = NSEC - 86400
      NYMD = INCYMD (NYMD,1)
      ENDDO
      ENDIF   
               
      IF (NSEC.EQ.86400)  THEN
      NSEC = 0
      NYMD = INCYMD (NYMD,1)
      ENDIF   
               
      IF (NSEC.LT.00000)  THEN
      DO WHILE (NSEC.LT.0)
      NSEC = 86400 + NSEC
      NYMD = INCYMD (NYMD,-1)
      ENDDO
      ENDIF   
               
      NHMS = NHMSF (NSEC)
      ENDIF   

      RETURN  
      end subroutine tick

      function incymd (NYMD,M)
!***********************************************************************        
!  PURPOSE                                                                      
!     INCYMD:  NYMD CHANGED BY ONE DAY                                          
!     MODYMD:  NYMD CONVERTED TO JULIAN DATE                                    
!  DESCRIPTION OF PARAMETERS                                                    
!     NYMD     CURRENT DATE IN YYMMDD FORMAT                                    
!     M        +/- 1 (DAY ADJUSTMENT)                                           
!                                                                               
!***********************************************************************        
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *        
!***********************************************************************        

      INTEGER NDPM(12)
      DATA    NDPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      LOGICAL LEAP
      LEAP(NY) = MOD(NY,4).EQ.0 .AND. (MOD(NY,100).NE.0 .OR. MOD(NY,400).EQ.0)

!***********************************************************************        
 
      NY = NYMD / 10000
      NM = MOD(NYMD,10000) / 100
      ND = MOD(NYMD,100) + M

      IF (ND.EQ.0) THEN
      NM = NM - 1
      IF (NM.EQ.0) THEN
          NM = 12
          NY = NY - 1
      ENDIF
      ND = NDPM(NM)
      IF (NM.EQ.2 .AND. LEAP(NY))  ND = 29
      ENDIF

      IF (ND.EQ.29 .AND. NM.EQ.2 .AND. LEAP(NY))  GO TO 20

      IF (ND.GT.NDPM(NM)) THEN
      ND = 1
      NM = NM + 1
      IF (NM.GT.12) THEN
          NM = 1
          NY = NY + 1
      ENDIF
      ENDIF

   20 CONTINUE
      INCYMD = NY*10000 + NM*100 + ND
      RETURN

!***********************************************************************        
!                      E N T R Y    M O D Y M D                                 
!***********************************************************************        

      ENTRY MODYMD (NYMD)
      NY = NYMD / 10000
      NM = MOD(NYMD,10000) / 100
      ND = MOD(NYMD,100)

   40 CONTINUE
      IF (NM.LE.1)  GO TO 60
      NM = NM - 1
      ND = ND + NDPM(NM)
      IF (NM.EQ.2 .AND. LEAP(NY))  ND = ND + 1
      GO TO 40

   60 CONTINUE
      MODYMD = ND
      RETURN
      end function incymd

      function nsecf (nhms)
!***********************************************************************
!  Purpose
!     Converts NHMS format to Total Seconds
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************
      implicit none
      integer  nhms, nsecf
      nsecf =  nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)
      return
      end function nsecf

      function nhmsf (nsec)
!***********************************************************************
!  Purpose
!     Converts Total Seconds to NHMS format
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************
      implicit none
      integer  nhmsf, nsec
      nhmsf =  nsec/3600*10000 + mod(nsec,3600)/60*100 + mod(nsec,60)
      return
      end function nhmsf

      function check_names (name1,name2)
      implicit none
      logical  check_names
      character(*)   name1,name2
      integer  len,i
      character*256 uname1,uname2
      character*1   c

! Convert name1 to All UpperCase uname1
! -------------------------------------
      len = len_trim(name1)
      uname1 = ''
      do i=1,len
         c = name1(i:i)
         if( ichar(c).ge.97 .and. ichar(c).le.122 ) then
             c = achar( ichar(c)-32 )
         endif
         uname1 = trim(uname1) // c
      enddo

! Convert name2 to All UpperCase uname2
! -------------------------------------
      len = len_trim(name2)
      uname2 = ''
      do i=1,len
         c = name2(i:i)
         if( ichar(c).ge.97 .and. ichar(c).le.122 ) then
             c = achar( ichar(c)-32 )
         endif
         uname2 = trim(uname2) // c
      enddo

! Compare uname1 and uname2
! -------------------------
    ! print *, 'Check Names: ',trim(uname1),' ',trim(uname2)
      check_names = ( trim(uname1) == trim(uname2) )
      return
      end

      subroutine hinterp ( qin,iin,jin,qout,iout,jout,mlev,undef,method )
      implicit   none
      integer    iin,jin,       iout,jout, mlev
      real   qin(iin,jin,mlev), qout(iout,jout,mlev)
      real undef,pi,dlin,dpin,dlout,dpout
      real dlam(iin), lons(iout*jout), lon
      real dphi(jin), lats(iout*jout), lat
      integer i,j,loc,method

      pi = 4.0*atan(1.0)
      dlin = 2*pi/iin
      dpin = pi/(jin-1)
      dlam(:) = dlin
      dphi(:) = dpin

      dlout = 2*pi/iout
      dpout = pi/(jout-1)
      
      loc = 0
      do j=1,jout
      do i=1,iout
      loc = loc + 1
      lon = -pi + (i-1)*dlout
      lons(loc) = lon
      enddo
      enddo

      loc = 0
      do j=1,jout
      lat = -pi/2.0 + (j-1)*dpout
      do i=1,iout
      loc = loc + 1
      lats(loc) = lat
      enddo
      enddo

      call interp_h ( qin,iin,jin,mlev,dlam,dphi,qout,iout*jout,lons,lats,undef,method )

      return
      end subroutine hinterp

      subroutine interp_h ( q_cmp,im,jm,lm,dlam,dphi,q_geo,irun,lon_geo,lat_geo,undef,method )
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Performs a horizontal interpolation from a field on a computational grid
!    to arbitrary locations.
!
!  INPUT:
!  ======
!    q_cmp ...... Field q_cmp(im,jm,lm) on the computational grid
!    im ......... Longitudinal dimension of q_cmp
!    jm ......... Latitudinal  dimension of q_cmp
!    lm ......... Vertical     dimension of q_cmp
!    dlam ....... Computational Grid Delta Lambda
!    dphi ....... Computational Grid Delta Phi
!    irun ....... Number of Output Locations
!    lon_geo .... Longitude Location of Output
!    lat_geo .... Latitude  Location of Output
!
!  OUTPUT:
!  =======
!    q_geo ...... Field q_geo(irun,lm) at arbitrary locations
!
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      implicit none

! Input Variables
! ---------------
      integer im,jm,lm,irun
      integer method

      real      q_geo(irun,lm)
      real    lon_geo(irun)
      real    lat_geo(irun)

      real    q_cmp(im,jm,lm)
      real     dlam(im)
      real     dphi(jm)

! Local Variables
! ---------------
      integer  i,j,l
      integer, allocatable       :: ip1(:), ip0(:), im1(:), im2(:)
      integer, allocatable       :: jp1(:), jp0(:), jm1(:), jm2(:)

! Bi-Linear Weights
! -----------------
      real, allocatable       ::    wl_ip0jp0 (:)
      real, allocatable       ::    wl_im1jp0 (:)
      real, allocatable       ::    wl_ip0jm1 (:)
      real, allocatable       ::    wl_im1jm1 (:)

! Bi-Cubic Weights
! ----------------
      real, allocatable       ::    wc_ip1jp1 (:)
      real, allocatable       ::    wc_ip0jp1 (:)
      real, allocatable       ::    wc_im1jp1 (:)
      real, allocatable       ::    wc_im2jp1 (:)
      real, allocatable       ::    wc_ip1jp0 (:)
      real, allocatable       ::    wc_ip0jp0 (:)
      real, allocatable       ::    wc_im1jp0 (:)
      real, allocatable       ::    wc_im2jp0 (:)
      real, allocatable       ::    wc_ip1jm1 (:)
      real, allocatable       ::    wc_ip0jm1 (:)
      real, allocatable       ::    wc_im1jm1 (:)
      real, allocatable       ::    wc_im2jm1 (:)
      real, allocatable       ::    wc_ip1jm2 (:)
      real, allocatable       ::    wc_ip0jm2 (:)
      real, allocatable       ::    wc_im1jm2 (:)
      real, allocatable       ::    wc_im2jm2 (:)

      real    ap1, ap0, am1, am2
      real    bp1, bp0, bm1, bm2

      real    lon_cmp(im)
      real    lat_cmp(jm)
      real    q_tmp(irun)

      real    pi,d
      real    lam,lam_ip1,lam_ip0,lam_im1,lam_im2
      real    phi,phi_jp1,phi_jp0,phi_jm1,phi_jm2
      real    dl,dp
      real    lam_cmp
      real    phi_cmp
      real    undef
      integer im1_cmp,icmp
      integer jm1_cmp,jcmp

! Initialization
! --------------
      pi = 4.*atan(1.)
      dl = 2*pi/ im     ! Uniform Grid Delta Lambda
      dp =   pi/(jm-1)  ! Uniform Grid Delta Phi

! Allocate Memory for Weights and Index Locations
! -----------------------------------------------
      allocate ( wl_ip0jp0(irun) , wl_im1jp0(irun) )
      allocate ( wl_ip0jm1(irun) , wl_im1jm1(irun) )
      allocate ( wc_ip1jp1(irun) , wc_ip0jp1(irun) , wc_im1jp1(irun) , wc_im2jp1(irun) )
      allocate ( wc_ip1jp0(irun) , wc_ip0jp0(irun) , wc_im1jp0(irun) , wc_im2jp0(irun) )
      allocate ( wc_ip1jm1(irun) , wc_ip0jm1(irun) , wc_im1jm1(irun) , wc_im2jm1(irun) )
      allocate ( wc_ip1jm2(irun) , wc_ip0jm2(irun) , wc_im1jm2(irun) , wc_im2jm2(irun) )
      allocate (       ip1(irun) ,       ip0(irun) ,       im1(irun) ,       im2(irun) )
      allocate (       jp1(irun) ,       jp0(irun) ,       jm1(irun) ,       jm2(irun) )

! Compute Input Computational-Grid Latitude and Longitude Locations
! -----------------------------------------------------------------
      lon_cmp(1) = -pi
      do i=2,im
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_cmp(1) = -pi*0.5
      do j=2,jm-1
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_cmp(jm) =  pi*0.5

! Compute Weights for Computational to Geophysical Grid Interpolation
! -------------------------------------------------------------------
      do i=1,irun
      lam_cmp = lon_geo(i)
      phi_cmp = lat_geo(i)

! Determine Indexing Based on Computational Grid
! ----------------------------------------------
      im1_cmp = 1
      do icmp = 2,im
      if( lon_cmp(icmp).lt.lam_cmp ) im1_cmp = icmp
      enddo
      jm1_cmp = 1
      do jcmp = 2,jm
      if( lat_cmp(jcmp).lt.phi_cmp ) jm1_cmp = jcmp
      enddo

      im1(i) = im1_cmp
      ip0(i) = im1(i) + 1
      ip1(i) = ip0(i) + 1
      im2(i) = im1(i) - 1

      jm1(i) = jm1_cmp
      jp0(i) = jm1(i) + 1
      jp1(i) = jp0(i) + 1
      jm2(i) = jm1(i) - 1

! Fix Longitude Index Boundaries
! ------------------------------
      if(im1(i).eq.im) then
      ip0(i) = 1
      ip1(i) = 2
      endif
      if(im1(i).eq.1) then
      im2(i) = im
      endif
      if(ip0(i).eq.im) then
      ip1(i) = 1
      endif


! Compute Immediate Surrounding Coordinates
! -----------------------------------------
      lam     =  lam_cmp
      phi     =  phi_cmp

! Compute and Adjust Longitude Weights
! ------------------------------------
      lam_im2 =  lon_cmp(im2(i))
      lam_im1 =  lon_cmp(im1(i))
      lam_ip0 =  lon_cmp(ip0(i))
      lam_ip1 =  lon_cmp(ip1(i))

      if( lam_im2.gt.lam_im1 ) lam_im2 = lam_im2 - 2*pi
      if( lam_im1.gt.lam_ip0 ) lam_ip0 = lam_ip0 + 2*pi
      if( lam_im1.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi
      if( lam_ip0.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi


! Compute and Adjust Latitude Weights   
! Note:  Latitude Index Boundaries are Adjusted during Interpolation
! ------------------------------------------------------------------
      phi_jm2 =  lat_cmp( min(max(1,jm2(i)),jm) )
      phi_jm1 =  lat_cmp( min(max(1,jm1(i)),jm) )
      phi_jp0 =  lat_cmp( min(max(1,jp0(i)),jm) )
      phi_jp1 =  lat_cmp( min(max(1,jp1(i)),jm) )

      if( jm2(i).eq.0    ) phi_jm2 = phi_jm1 - dphi(1)
      if( jm1(i).eq.jm   ) then
                           phi_jp0 = phi_jm1 + dphi(jm-1)
                           phi_jp1 = phi_jp0 + dphi(jm-2)
      endif
      if( jp1(i).eq.jm+1 ) phi_jp1 = phi_jp0 + dphi(jm-1)


! Bi-Linear Weights
! -----------------
              d    = (lam_ip0-lam_im1)*(phi_jp0-phi_jm1)
      wl_im1jm1(i) = (lam_ip0-lam    )*(phi_jp0-phi    )/d
      wl_ip0jm1(i) = (lam    -lam_im1)*(phi_jp0-phi    )/d
      wl_im1jp0(i) = (lam_ip0-lam    )*(phi    -phi_jm1)/d
      wl_ip0jp0(i) = (lam    -lam_im1)*(phi    -phi_jm1)/d

! Bi-Cubic Weights
! ----------------
      ap1 = ( (lam    -lam_ip0)*(lam    -lam_im1)*(lam    -lam_im2) ) &
          / ( (lam_ip1-lam_ip0)*(lam_ip1-lam_im1)*(lam_ip1-lam_im2) )
      ap0 = ( (lam_ip1-lam    )*(lam    -lam_im1)*(lam    -lam_im2) ) &
          / ( (lam_ip1-lam_ip0)*(lam_ip0-lam_im1)*(lam_ip0-lam_im2) )
      am1 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam    -lam_im2) ) &
          / ( (lam_ip1-lam_im1)*(lam_ip0-lam_im1)*(lam_im1-lam_im2) )
      am2 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam_im1-lam    ) ) &
          / ( (lam_ip1-lam_im2)*(lam_ip0-lam_im2)*(lam_im1-lam_im2) )

      bp1 = ( (phi    -phi_jp0)*(phi    -phi_jm1)*(phi    -phi_jm2) ) &
          / ( (phi_jp1-phi_jp0)*(phi_jp1-phi_jm1)*(phi_jp1-phi_jm2) )
      bp0 = ( (phi_jp1-phi    )*(phi    -phi_jm1)*(phi    -phi_jm2) ) &
          / ( (phi_jp1-phi_jp0)*(phi_jp0-phi_jm1)*(phi_jp0-phi_jm2) )
      bm1 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi    -phi_jm2) ) &
          / ( (phi_jp1-phi_jm1)*(phi_jp0-phi_jm1)*(phi_jm1-phi_jm2) )
      bm2 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi_jm1-phi    ) ) &
          / ( (phi_jp1-phi_jm2)*(phi_jp0-phi_jm2)*(phi_jm1-phi_jm2) )

      wc_ip1jp1(i) = bp1*ap1
      wc_ip0jp1(i) = bp1*ap0
      wc_im1jp1(i) = bp1*am1
      wc_im2jp1(i) = bp1*am2

      wc_ip1jp0(i) = bp0*ap1
      wc_ip0jp0(i) = bp0*ap0
      wc_im1jp0(i) = bp0*am1
      wc_im2jp0(i) = bp0*am2

      wc_ip1jm1(i) = bm1*ap1
      wc_ip0jm1(i) = bm1*ap0
      wc_im1jm1(i) = bm1*am1
      wc_im2jm1(i) = bm1*am2

      wc_ip1jm2(i) = bm2*ap1
      wc_ip0jm2(i) = bm2*ap0
      wc_im1jm2(i) = bm2*am1
      wc_im2jm2(i) = bm2*am2

      enddo

! Interpolate Computational-Grid Quantities to Geophysical Grid
! -------------------------------------------------------------
      do L=1,lm
      do i=1,irun

      if( lat_geo(i).le.lat_cmp(2)     .or.  &
          lat_geo(i).ge.lat_cmp(jm-1)  .or.  &
          method.eq.1                ) then

! 1st Order Interpolation at Poles
! --------------------------------
      if( q_cmp( im1(i),jm1(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jm1(i),L ).ne.undef  .and.  &
          q_cmp( im1(i),jp0(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jp0(i),L ).ne.undef ) then

      q_tmp(i) = wl_im1jm1(i) * q_cmp( im1(i),jm1(i),L ) &
               + wl_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L ) &
               + wl_im1jp0(i) * q_cmp( im1(i),jp0(i),L ) &
               + wl_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )

      else
      q_tmp(i) = undef
      endif

      else

! Cubic Interpolation away from Poles
! -----------------------------------
      if( q_cmp( ip1(i),jp0(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jp0(i),L ).ne.undef  .and.  &
          q_cmp( im1(i),jp0(i),L ).ne.undef  .and.  &
          q_cmp( im2(i),jp0(i),L ).ne.undef  .and.  &

          q_cmp( ip1(i),jm1(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jm1(i),L ).ne.undef  .and.  &
          q_cmp( im1(i),jm1(i),L ).ne.undef  .and.  &
          q_cmp( im2(i),jm1(i),L ).ne.undef  .and.  &

          q_cmp( ip1(i),jp1(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jp1(i),L ).ne.undef  .and.  &
          q_cmp( im1(i),jp1(i),L ).ne.undef  .and.  &
          q_cmp( im2(i),jp1(i),L ).ne.undef  .and.  &

          q_cmp( ip1(i),jm2(i),L ).ne.undef  .and.  &
          q_cmp( ip0(i),jm2(i),L ).ne.undef  .and.  &
          q_cmp( im1(i),jm2(i),L ).ne.undef  .and.  &
          q_cmp( im2(i),jm2(i),L ).ne.undef ) then

      q_tmp(i) = wc_ip1jp1(i) * q_cmp( ip1(i),jp1(i),L )  &
               + wc_ip0jp1(i) * q_cmp( ip0(i),jp1(i),L )  &
               + wc_im1jp1(i) * q_cmp( im1(i),jp1(i),L )  &
               + wc_im2jp1(i) * q_cmp( im2(i),jp1(i),L )  &

               + wc_ip1jp0(i) * q_cmp( ip1(i),jp0(i),L )  &
               + wc_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )  &
               + wc_im1jp0(i) * q_cmp( im1(i),jp0(i),L )  &
               + wc_im2jp0(i) * q_cmp( im2(i),jp0(i),L )  &

               + wc_ip1jm1(i) * q_cmp( ip1(i),jm1(i),L )  &
               + wc_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )  &
               + wc_im1jm1(i) * q_cmp( im1(i),jm1(i),L )  &
               + wc_im2jm1(i) * q_cmp( im2(i),jm1(i),L )  &

               + wc_ip1jm2(i) * q_cmp( ip1(i),jm2(i),L )  &
               + wc_ip0jm2(i) * q_cmp( ip0(i),jm2(i),L )  &
               + wc_im1jm2(i) * q_cmp( im1(i),jm2(i),L )  &
               + wc_im2jm2(i) * q_cmp( im2(i),jm2(i),L )

      elseif( q_cmp( im1(i),jm1(i),L ).ne.undef  .and.    &
              q_cmp( ip0(i),jm1(i),L ).ne.undef  .and.    &
              q_cmp( im1(i),jp0(i),L ).ne.undef  .and.    &
              q_cmp( ip0(i),jp0(i),L ).ne.undef ) then

      q_tmp(i) = wl_im1jm1(i) * q_cmp( im1(i),jm1(i),L )  &
               + wl_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )  &
               + wl_im1jp0(i) * q_cmp( im1(i),jp0(i),L )  &
               + wl_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )

      else
      q_tmp(i) = undef
      endif

      endif
      enddo

! Load Temp array into Output array
! ---------------------------------
      do i=1,irun
      q_geo(i,L) = q_tmp(i)
      enddo
      enddo

      deallocate ( wl_ip0jp0 , wl_im1jp0 )
      deallocate ( wl_ip0jm1 , wl_im1jm1 )
      deallocate ( wc_ip1jp1 , wc_ip0jp1 , wc_im1jp1 , wc_im2jp1 )
      deallocate ( wc_ip1jp0 , wc_ip0jp0 , wc_im1jp0 , wc_im2jp0 )
      deallocate ( wc_ip1jm1 , wc_ip0jm1 , wc_im1jm1 , wc_im2jm1 )
      deallocate ( wc_ip1jm2 , wc_ip0jm2 , wc_im1jm2 , wc_im2jm2 )
      deallocate (       ip1 ,       ip0 ,       im1 ,       im2 )
      deallocate (       jp1 ,       jp0 ,       jm1 ,       jm2 )

      return
      end subroutine interp_h

      subroutine usage()
      write(*, '()' )
      write(*, '(" Usage:  ")' )
      write(*, '(" ------  ")' )
      write(*, '(" 3CH.x -exp1  expid plus filelist for 1st 3CH Triplet member (exp1:exp2:exp3)")' )
      write(*, '("       -exp2  expid plus filelist for 2nd 3CH Triplet member (exp1:exp2:exp3)")' )
      write(*, '("       -exp3  expid plus filelist for 3rd 3CH Triplet member (exp1:exp2:exp3)")' )
      write(*, '()' )
      write(*, '("       -nymdb Beginning Date for 3CH computation")' )
      write(*, '("       -nhmsb Beginning Time for 3CH computation")' )
      write(*, '("       -nymde    Ending Date for 3CH computation")' )
      write(*, '("       -nhmse    Ending Time for 3CH computation")' )
      write(*, '()' )
      write(*, '(" Optional Arguments:  ")' )
      write(*, '(" -------------------  ")' )
      write(*, '("      <-exp4  expid plus filelist for 3rd 3CH Triplet member of 2nd Triplet (exp1:exp2:exp4)>")' )
      write(*, '("      <-exp5  expid plus filelist for 3rd 3CH Triplet member of 3rd Triplet (exp1:exp2:exp5)>")' )
      write(*, '("         .")' )
      write(*, '("         .")' )
      write(*, '("      <-expN  expid plus filelist for 3rd 3CH Triplet member of N-2 Triplet (exp1:exp2:expN)>")' )
      write(*, '()' )
      write(*, '("      <-ndt   Time Frequency in seconds  (Default:  21600)>")' )
      write(*, '()' )
      write(*, '("      <-levs   levels used for 3CH computation (e.g., 1000 925 850 700 ) (Default: ALL coincident levels)>")' )
      write(*, '("      <-fields fields used for 3CH computation (Default: u v t q h)>")' )
      write(*, '("      <-im     Longitudinal Dimension for Output (Default: Minimum among input experiments)>")' )
      write(*, '("      <-jm     Latitudinal  Dimension for Output (Default: Minimum among input experiments)>")' )
      write(*, '()' )
      write(*, '("      <-method  1 => Bilinear    Interpolation  (Default)")' )
      write(*, '("                2 => BiCubic     Interpolation")' )
      write(*, '("                3 => Box-Average Interpolation")' )
      write(*, '()' )
      write(*, '("      <-tag    tag used for Output Filename (Default: grads)>")' )
      write(*, '()' )
      write(*, '(" creates:  tag.data")' )
      write(*, '(" --------  tag.ctl")' )
      write(*, '()' )
      write(*, '("           Contains  3-D fields of 3CH Variance Error Estimates for each Triplet Member: var(expN error),")' )
      write(*, '("           plus full 3-D fields of Member Means: ave(expN), Variance Differences: var(expM-expN)")' )
      write(*, '("                and  3-D fields of Mean Square Error Decomposition: Total, Bias, Amplitude, Phase)")' )
      write(*, '()' )
      error stop 7
      end subroutine usage

