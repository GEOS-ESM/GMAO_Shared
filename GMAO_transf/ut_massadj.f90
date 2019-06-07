      program ut_massadj

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ut_massadj: Ensure analysis does not change dry mass 
!
! !USAGE: see the routine usage() below
!
! !USES:
!
      use m_dyn,     only : dyn_get
      use m_dyn,     only : dyn_put
      use m_dyn,     only : dyn_vect
      use m_dyn,     only : dyn_clean

      use m_massadj, only : madj_init
      use m_massadj, only : madj_pre
      use m_massadj, only : madj_post
      use m_massadj, only : madj_clean

      use m_realkinds, only : r_kind => kind_r8
!     use kinds, only : r_kind

      implicit NONE

! !DESCRIPTION: Unit test to software for ensuring dry mass is 
!               kept unchanged by analysis
!
! !REVISION HISTORY:
!
!  25Jan2007   Todling  Initial unit-test code.
!
!-------------------------------------------------------------------------
!EOP

      character(len=*), parameter :: myname = 'ut_massadj'

      integer :: NQ = 1 ! force mass adjustment to be applied to water-vapor only
                      
!     File names
!     ----------
      integer, parameter :: MFILES = 2 ! max.   number of input files
      integer, parameter :: PREC   = 0 ! 32 bit output
      character(len=255) :: files(MFILES), binfiles(MFILES), etafile, st, hroot
      character(len=255) :: dynout     ! output filename
      integer            :: nfiles     ! actual no. of input files
      integer            :: nstep
      integer            :: rc
      integer            :: i,j,k
      real(r_kind)       :: delpmin

!     Dynamics/simulator vectors
!     --------------------------
      type(dyn_vect) w_b
      type(dyn_vect) w_a

!     Locals
!     ------
      integer ntimes, n, freq, nymd_b, nhms_b, nymd_a, nhms_a
      integer im, jm, km, system 

!     Initialize
!     ----------     
      call Init_ ( mfiles, files, dynout )

!     Read in backgroud fields
!     ------------------------
      n = 1 ! assumes only one file in file
      call dyn_get ( trim(files(1)), nymd_b, nhms_b, w_b, rc, timidx=n, freq=freq, nstep=nstep )
         if ( rc .ne. 0 ) then
            call die(myname,'cannot read background field')
         end if

      call dyn_get ( trim(files(2)), nymd_a, nhms_a, w_a, rc, timidx=n, freq=freq, nstep=nstep )
         if ( rc .ne. 0 ) then
            call die(myname,'cannot read dynamics vector file')
         end if

!     Check dates/time
!     ----------------
      print *, 'bkg data/time ', nymd_b, nhms_b
      print *, 'ana data/time ', nymd_a, nhms_a
      if ( nymd_a/=nymd_b .and. nhms_a/=nhms_b ) then
           call die(myname,'unmatched date and times')
      endif

!     Change analysis so it preserves dry-mass in background
!     ------------------------------------------------------

      call madj_init ( w_b%grid%im, w_b%grid%jm, w_b%grid%km, w_b%grid%lm )
      call madj_pre  ( w_b%grid%im, w_b%grid%jm, w_b%grid%km, NQ, w_b%grid%ak, w_b%grid%bk, &
                       w_b%ps, w_b%delp, w_b%q, rc )

      call madj_post ( w_a%grid%im, w_a%grid%jm, w_a%grid%km, NQ, 0, w_a%grid%ptop, &
                       w_a%ps, w_a%delp, w_a%q, rc, delpmin )
      call madj_clean ()

      print *, ' delpmin = ', delpmin


!     Write out mass-fixed analysis
!     -----------------------------
      call dyn_put ( trim(dynout), nymd_a, nhms_a, prec, w_a, rc, freq=freq, nstep=nstep )


!     Clean up mess
!     -------------
      call dyn_clean ( w_b )
      call dyn_clean ( w_a )

!     All done
!     --------
      call exit(0)

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Init_ --- Initialize dyn2dyn
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
      subroutine Init_ ( mfiles, files, dynout )

      implicit NONE

      integer,       intent(in)  :: mfiles  ! max. number of eta files
                                            ! dynamics file names (eta)
      character*255, intent(out) :: files(mfiles) 
      character*255, intent(out) :: dynout
      
! !REVISION HISTORY:
!
!       25Jan2007  Todling   Initial code (based on dyndiff.
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, iarg, argc, iargc
      character(len=255) :: etafile, argv

      print *
      print *, '     ------------------------------------------------------'
      print *, '     ut_massadj - ensure analysis keeps dry mass unchanged '
      print *, '     ------------------------------------------------------'
      print *

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 2 ) call usage()

      iarg = 0
      nfiles = 0
      dynout = 'massadj.eta.hdf'

      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iarg, argv )
         select case (argv)
           case ("-h")
             if ( iarg+1 .gt. argc ) call usage()
           case ("-o")
             iarg = iarg + 1
             call GetArg ( iarg, dynout )
           case default
             nfiles = nfiles + 1
             if ( nfiles .gt. mfiles ) call usage()
             files(nfiles) = argv
         end select
      end do

      if ( nfiles .ne. 2 ) call usage()

!     Echo the parameters
!     -------------------
      print *
      print *, '------------------------------------------------------------------'
      print *, '  Eta     Dynamics state files: '
      do i = 1, nfiles
        print *, i, ': ', trim(files(i))
      end do

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *,'Usage: '
      print *
      print *,'  ut_massadj.x [-h] [-o FNAME] bkgeta anaeta'
      print *
      print *, 'where'
      print *
      print *, '-h          Help (optional)'
      print *
      print *, '-o FNAME    filename of outout/mass corrected analysis (def: massadj.eta.hdf)'
      print *
      call exit(1)
      end subroutine usage
      
!.................................................................

      subroutine die ( myname, msg )
      character(len=*) :: myname, msg
      write(*,'(a)') trim(myname) // ': ' // trim(msg)
      call exit(1)
      end subroutine die

!.................................................................

  end program ut_massadj
