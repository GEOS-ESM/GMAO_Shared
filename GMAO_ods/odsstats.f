!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsstats:  Operates on specific attributes in ODS
!
! !INTERFACE:

      program odsstats

!
! !DESCRIPTION:
!           Reads ODS files, selects data according to criteria specified
!           via a resource file, operates on chosen attribute and summarizes
!           results by writing them to an ASC file.
!

      use m_odsmeta, only : KTMAX, KXMAX
      use m_ods
      use m_inpak90
      use m_chars,   only: lowercase
      use m_stdio,   only: stdout, stderr

      implicit NONE

!
! !REVISION HISTORY:
!
!     23Jan2008 Todling  Initial code
!     23Sep2008 Todling  Add stats summary and skip verbose
!     01Sep2009 Todling  Add obvar to handle attr other than xvec
!     13Nov2009 Todling  Accommulate global rms
!     13Apr2010 Gelaro   Accumulate numneg, increase NCMAX LVMAX for IASI 
!     15Aug2013 Todling/Daescu  Implement calc of imp from sensitvity
!     19Feb2014 Todling/Daescu  Implement calc of oma * xvec 
!     24Feb2014 Todling  Revisit Write_AccumStats (now to file and/or screen)
!     05Aug2014 Todling  - Replace oma * xvec w/ amo * xvec (more meaningful)
!                        - Ability to write out R-scaling factors based on amo*xvec 
!     11Dec2018 Todling Add DFS as in Lupu et al. 2011; eq. (10) - see also odsmatch
!     18Dec2020 Todling Add ncf opt to allow running from diag-bin files
!
!EOP
!BOC

      character(len=*), parameter :: myname = 'odsstats'

      integer, parameter :: NFILES_MAX = 1488! Max number of input files

      integer, parameter ::  NOMAX = 5       ! max no. of operations
      integer, parameter ::  NCMAX = 999     ! max no. of classes
      integer, parameter ::  NTMAX = 12*31*4 ! 4 syn times a day, for a year (roughly)
      integer, parameter ::  LVMAX = 800     ! Max number of levels/channels

      logical, parameter :: debug = .false.
      logical            :: verb  = .false.

!     Local variables
!     ---------------
      logical :: log_transf = .false.
      integer nfiles, ifile, lf, isyn, ksyn, nymd, nhms, nkt, nkx, nlv, ncfound
      integer i, j, ic, nc, nt, ierr, nobs, nsel, synhour, nop, nops
      integer nymdb, nhmsb
      logical aodsigofix
      logical jedi

      character(len=255) :: opers(nomax)          ! operations
      character(len=255) :: oclass(ncmax)         ! observations classes
      character(len=255) :: obvar                 ! get this attribute from ODS (usually, xvec)     
      character(len=255) :: RCfile
      integer            :: kttype(ktmax,ncmax)   !
      integer            :: kxtype(kxmax,ncmax)   !
      real               :: levlst(lvmax,ncmax)   ! levels/channels of selected obs
      real, pointer      :: ptr(:,:)              ! number of obs giving positive attribute
      real, target       :: obssum(ncmax,ntmax)   ! sum per class and syn hour
      real, target       :: negsum(ncmax,ntmax)   ! number of obs for which attr are negative
      real, target       :: neusum(ncmax,ntmax)   ! number of obs for which attr are neutral
      real, target       :: obsgms(ncmax,ntmax)   ! global mean square of obs for each attr
      real               :: obserr(ncmax,ntmax)   ! observations error
      integer            :: obsnum(ncmax,ntmax)   ! number of obs per class and syn hour
      integer            :: nymda(ntmax)          ! all dates
      integer            :: nhmsa(ntmax)          ! all times
      integer            :: trange(2)             ! time interval of selected obs
      real               :: latrange(2)           ! latitudes  of selected obs
      real               :: lonrange(2)           ! longitudes of selected obs
      real               :: accum_obssum(ncmax)   ! accumulated number of obs per class
      real               :: accum_obsgms(ncmax)   ! accumulated global mean square of obs per class
      integer            :: accum_obsnum(ncmax)   ! accumulated impact per class
      integer            :: accum_obsneg(ncmax)   ! accumulated number of obs with neg attribute
      integer            :: accum_obsneu(ncmax)   ! accumulated number of obs with neutral attribute
 
      integer,allocatable:: igood(:)

      character*255 infile (NFILES_MAX) ! input filenames
      character*255 fileout, outfile, outodsfn   ! output filename
      character*80  ftype
      logical       anotherclass, allkxs, allevs, lrms, ncf

!     storage for ODS:
!     ---------------
      type ( ods_vect ) ods
      type ( ods_vect ) odss

!     Option flags:
!     ------------
      synhour = -1       ! DEFAULT: process all synoptic times

      oclass = 'NONE'
      kxtype = 0
      kttype = 0
      levlst = -1.0
      lrms = .false.

!     Parse in command line
!     ---------------------
      call init ( infile, nfiles_max, nfiles, RCfile, outfile, trange, 
     .            latrange, lonrange, verb, jedi, aodsigofix, ncf )

!     Read in resource file
!     ---------------------
      call ObsOperSet_ ( RCfile, ierr )
      if(ierr>0) then
         print *, 'Trouble handling RC file tables, ierr = ', ierr
         print *, 'Aborting ...'
         stop(1)
      endif

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      nt   = 0
      accum_obsnum = 0
      accum_obssum = 0.0
      accum_obsgms = 0.0
      accum_obsneg = 0
      accum_obsneu = 0
      obssum = 0.0
      negsum = 0.0
      neusum = 0.0
      obsgms = 0.0
      nymdb = -1
      nhmsb = -1
      if(ncf) ksyn = 1

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!       Loop over all synoptic times on this file
!       -----------------------------------------
        do isyn = 1, ksyn

!         Read all data for this synoptic time
!         ------------------------------------

	  nymd = -1            ! get data for next synoptic time on file
	  nhms =  0
          call ODSNxTime ( trim(infile(ifile)), nymd, nhms )
          if ( nymd .eq. -1 ) then 
               if(verb) print *, 'End-Of-File'
               exit     ! skip to next file
          end if
	  
          if(verb) print *, 'calling ODS_Get'

          call ODS_Get ( trim(infile(ifile)), nymd, nhms, ftype, ods, ierr, ncf=ncf )

          if(verb) print *, 'completed ODS_Get'

          if ( ierr .gt. 0 ) then
               print *, 'ODS_Get error: ierr = ', ierr
               exit     ! skip to next file
          end if

!         Set number of observations found in the file
!         --------------------------------------------
          nobs   =  ods%data%nobs

          if ( nobs .eq. 0 ) then
               print *, 'No data for this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if ( synhour .ge. 0 .and. nhms .ne. synhour*10000 ) then
               if(verb) print *, 'Skipping this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          nt = nt + 1
          if ( nt>ntmax) then
              print *, 'Max number of syn hours reached'
              print *, 'Need to lower number if processed files'
              stop
          endif

          if ( nymdb<0.and.nhmsb<0 ) then
               nymdb=nymd; nhmsb=nhms
          endif
          if (verb) then
              print *, 'number of obs read = ', nobs
              print *, 'date = ', nymd
              print *, 'time = ', nhms
          endif

!         A little "QC" to avoid including outliers in stats
!         --------------------------------------------------
          if(aodsigofix .and. trim(obvar) == 'dfs' .or. trim(obvar) == 'imp0hr') call aodfix_sigo (ods,log_transf)

!         Loop over observation classes
!         -----------------------------
          do nc = 1, ncfound
             if(verb) print *, trim(oclass(nc))

!            Select observations need for operating with present class
!            ---------------------------------------------------------
             call icount_this_ ( kttype(:,nc), nkt )
             if ( allkxs .and. allevs ) then
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0, 
     .                                        kt_list=kttype(1:nkt,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             else if ( allevs .and. (.not. allkxs ) ) then
                  call icount_this_ ( kxtype(:,nc), nkx )
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0, kx_list=kxtype(1:nkx,nc),
     .                                        kt_list=kttype(1:nkt,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             else 
                  call icount_this_ ( kxtype(:,nc), nkx )
                  call rcount_this_ ( levlst(:,nc), nlv )
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0,  kx_list=kxtype(1:nkx,nc),
     .                                         kt_list=kttype(1:nkt,nc),
     .                                        lev_list=levlst(1:nlv,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             endif

!            A little "QC" to avoid including outliers in stats
!            --------------------------------------------------
             allocate(igood(nsel))
             do i = 1,nsel
                igood(i) = i
             enddo
!            igood = ((/i/),i=1,nsel)
             if (.not. ncf) then
                if(trim(obvar)/='xvec')  call no_outliers (odss,nsel,log_transf,igood) !  sigo must be available in this case
             endif

             obsnum(nc,nt) = nsel
             nymda(nt)     = nymd
             nhmsa(nt)     = nhms
             if (jedi) then
                 call sumattr ( verb, odss, nsel, igood, obvar, obserr(nc,nt) )
             endif
             do nop = 1, nops
                if(opers(nop)=='sum')    call sumattr ( verb, odss, nsel, igood, obvar, obssum(nc,nt) )
                if(opers(nop)=='numneg') call negattr ( verb, odss, nsel, igood, obvar, negsum(nc,nt) )
                if(opers(nop)=='numneu') call neuattr ( verb, odss, nsel, igood, obvar, neusum(nc,nt) )
                if(opers(nop)=='rms')  then
                   call gmsattr ( verb, odss, nsel, 1, igood, obvar, obsgms(nc,nt) )
                   lrms = .true.
                 endif
                if(opers(nop)=='stddev')  then
                   call gmsattr ( verb, odss, nsel, 2, igood, obvar, obsgms(nc,nt) )
                 endif
             enddo
             deallocate(igood)

!            Accumulate obs counts and impacts
!            ---------------------------------
             accum_obsnum(nc) = accum_obsnum(nc) + obsnum(nc,nt)
             accum_obssum(nc) = accum_obssum(nc) + obssum(nc,nt)
             accum_obsgms(nc) = accum_obsgms(nc) + obsgms(nc,nt)
             accum_obsneg(nc) = accum_obsneg(nc) + negsum(nc,nt)
             accum_obsneu(nc) = accum_obsneu(nc) + neusum(nc,nt)

!            When debugging, write out select ODS entries
!            ---------------------------------------------
             if ( debug ) then
                  write(outodsfn,'(2a,i8.8,a,i2.2,a)') trim(oclass(nc)), '.obs.', nymd, '_', nhms/10000, 'z.ods'
                  print *, 'nsel, nobs', nsel, odss%data%nobs
                  print *, 'Writing ods file for debug purposes: ', trim(outodsfn)
                  call ODS_Put ( trim(outodsfn), ftype, nymd, nhms, odss, ierr, append=.true. )
             endif

             call ODS_Clean ( odss, ierr ) 
          end do

!         Write to ASC file
!         -----------------
             if(verb) print *, 'calling Write_Stats'
          do nop = 1, nops
             if(opers(nop)=='sum') then
                ptr => obssum
             endif
             if(opers(nop)=='numneg') then
                ptr => negsum
             endif
             if(opers(nop)=='rms') then
                ptr => obsgms
             endif
             fileout = trim(outfile) // '_' // trim(opers(nop)) // '.txt'
             call Write_Stats ( fileout, ptr, obsnum, obserr, oclass, ncfound,
     .                          ncmax, nt, nymda, nhmsa, verb, jedi, ierr ) 
             if(opers(nop)=='sum') then
                fileout = trim(outfile) // '.' // 'rcov4gsi'
                call Write_gsiRfactor ( RCfile, fileout, obssum, oclass, ncfound, ncmax, nt, 
     .                                  nymda, nhmsa, verb, ierr ) 
             endif
          enddo
             if(verb) print *, 'completed Write_Stats'

          call ODS_Clean ( ods, ierr )
             if ( ierr .ne. 0 ) then
                  print *, 'ODS_Clean error: ierr = ', ierr
             endif

        end do  ! loop over synoptic times

      end do  ! loop over files

!     Write out accumulated results
!     ------------------------------
      if (.not.jedi) then
         fileout = trim(outfile) // '_' // 'all.txt'
         call Write_AccumStats ( fileout, accum_obssum, accum_obsgms, accum_obsnum, accum_obsneg, 
     .                           accum_obsneu, oclass, ncfound, ncmax, nymdb, nhmsb, lrms, verb, ierr ) 
      endif


      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ObsOpertSet --- General paramater setup for singular vector calc
! 
! !INTERFACE:
!
      subroutine ObsOperSet_ ( RCfile, stat )
 
! !USES: 
    
      Implicit None

! !INPUT PARAMETERS: 
!
      character(len=*), intent(in) :: RCfile  ! resource filename

! !OUTPUT PARAMETERS:

      integer,          intent(out) :: stat                  ! return error code


! !DESCRIPTION:  Initialize observation operations program.
!
! !REVISION HISTORY: 
!
!   23Jan2008  Todling    Initial code.
!   05Aug2008  Todling    Add table of levs/channels.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'ObsOperSet_'

      character(len=255) token, obsoperrc, tablename

      integer irow, j, n, kt, kx, jcnt, maxn, iret
      integer ii, lt, k1, k2, knext
      real    r1, r2, rnext

      stat = 0

!     Load resources file
!     -------------------
      obsoperrc = ' '
      call getenv('obs_opers.rc',OBSOPERRC)     ! Unix binding
      if(obsoperrc.eq.' ') obsoperrc=RCfile     ! default name
      call i90_loadf (trim(obsoperrc), iret)
      if( iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
          stat = 1
          return
      end if
      write(stdout,'(1x,  a  )') '---------------------------------------------------'
      write(stdout,'(1x, 2a  )') myname_, ': Reading resource file'
      write(stdout,'(1x,  a,/)') '---------------------------------------------------'

!     Define type of operations
!     -------------------------
      obvar = 'xvec'
      call I90_label('OBS*Variable:', iret)
      if (iret .eq. 0) then
          call I90_Gtoken ( token, iret )
          if ( iret==0 ) then
               obvar = lowercase(trim(token))
          endif
      end if
      write(stdout,'(2a)') ' Will extract information from this attribute: ', trim(obvar)

!     Define type of operations
!     -------------------------
      opers(1:nomax) = 'NONE'
      call I90_label('OBS*Operations:', iret)
      if (iret .eq. 0) then
          nops = 0
          do while ( iret==0 .and. jcnt<nomax )
             call I90_Gtoken ( token, iret )
             if ( iret==0 ) then
                  nops = nops + 1
                  opers(nops) = lowercase(trim(token))
                  write(stdout,'(2a)') ' Will perform the following operations: ', trim(opers(nops))
             endif
          enddo
      else
         nops      = 1
         opers(1)  = 'sum'   ! default: sum residuals
                             ! other possibilities are:
                             !  ave - average residuals
      end if

!     Read table with variable types
!     ------------------------------
      tablename = 'OBS*Class*Kts::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(2a,i5,2a)') myname_, ': I90_label error, iret=', iret,
     .                                       ': trying to read ', trim(tablename)
         stat = 2; return
      end if
      irow = 0
      ncfound = 0
      write(stdout,'(a)') ' Will process the following obs classes: '
      do while (iret==0)                   ! read table entries
         call I90_GLine ( iret )           ! iret=-1: end of file; +1: end of table
         if (iret==0.and.irow<ncmax) then  ! OK, we have next row of table
             irow = irow + 1

             call I90_GToken(token, iret ) ! obs class name
             if (iret/=0) then
                 write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
             end if
             oclass(irow) = trim(token)
             write(stdout,'(1x,a)') trim(oclass(irow))

             jcnt=0
             ierr=0
             do  j = 1, ktmax
               call I90_GToken(token, ierr )
               if(ierr/=0) exit
               ii = index(token,':') ! token is single entry or range of entries
               lt = len_trim(token)
               if (ii==0) then       ! no colon, therefore single entry
                   read(token,*) knext
                   jcnt = jcnt + 1
                   kttype(jcnt,irow) = knext
               else                  ! colon, therefore kx1:kx2
                   read(token(1:ii-1),*) k1
                   read(token(ii+1:lt),*) k2
                   do knext = k1, k2
                      if (jcnt==KTMAX) then    ! check space
                        write(stderr,'(2a,i5)') myname,': increase KTMAX'
                        stat = 4; return
                      else if (ierr==0) then
                        jcnt = jcnt + 1
                        kttype(jcnt,irow) = knext
                      end if
                   end do
               end if
             enddo

             ncfound = max(ncfound,irow)
         end if
      end do
      nkt = jcnt
      print *, 'Will process the following Kt''s:'
      do j = 1, ncfound
         call icount_this_ ( kttype(:,j), n )
         print *, kttype(1:n,j)
      enddo

!     Read table with instruments types
!     ---------------------------------
      allkxs = .false.
      tablename = 'OBS*Class*Kxs::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(4a)') myname_, ': table ', trim(tablename),
     .                        'not found in RC file ... will take all KXs'
         allkxs = .true.
      end if
      if ( .not. allkxs ) then
        irow = 0
        do while (iret==0)                     ! read table entries
           call I90_GLine ( iret )             ! iret=-1: end of file; +1: end of table
           if (iret==0.and.irow<ncmax) then    ! OK, we have next row of table
               irow = irow + 1
  
               call I90_GToken ( token, iret ) ! obs class name
               if (iret/=0) then
                   write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
               end if
               if(trim(token)/=trim(oclass(irow))) then     ! tables of kx and kt must be ordered 
                                                            ! in the same way, otherwise drop it
                  stat =4; return
               endif
  
               jcnt=0
               ierr=0
               do  j = 1, kxmax
                 call I90_GToken(token, ierr )
                 if(ierr/=0) exit
                 ii = index(token,':') ! token is single entry or range of entries
                 lt = len_trim(token)
                 if (ii==0) then       ! no colon, therefore single entry
                     read(token,*) knext
                     jcnt = jcnt + 1
                     kxtype(jcnt,irow) = knext
                 else                  ! colon, therefore kx1:kx2
                     read(token(1:ii-1),*) k1
                     read(token(ii+1:lt),*) k2
                     do knext = k1, k2
                        if (jcnt==KXMAX) then    ! check space
                          write(stderr,'(2a,i5)') myname,': increase KXMAX'
                          stat = 4; return
                        else if (ierr==0) then
                          jcnt = jcnt + 1
                          kxtype(jcnt,irow) = knext
                        end if
                     end do
                 end if
               enddo

           end if
        end do
        nkx = jcnt
        print *, 'Will process the following Kx''s:'
        do j = 1, ncfound
           call icount_this_ ( kxtype(:,j), n )
           print *, kxtype(1:n,j)
        enddo
      endif ! < all KXs >

!     Read table with levels/channels
!     -------------------------------
      allevs = .false.
      tablename = 'OBS*Class*Lev::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(4a)') myname_, ': table ', trim(tablename),
     .                        'not found in RC file ... will take all levels/channels'
         allevs = .true.
      end if
      if ( .not. allevs ) then
        irow = 0
        do while (iret==0)                     ! read table entries
           call I90_GLine ( iret )             ! iret=-1: end of file; +1: end of table
           if (iret==0.and.irow<ncmax) then    ! OK, we have next row of table
               irow = irow + 1
  
               call I90_GToken ( token, iret ) ! obs class name
               if (iret/=0) then
                   write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
               end if
               if(trim(token)/=trim(oclass(irow))) then     ! tables of kx and kt must be ordered 
                                                            ! in the same way, otherwise drop it
                  stat =5; return
               endif
  
               jcnt=0
               ierr=0
               do  j = 1, lvmax
                 call I90_GToken(token, ierr )
                 if(ierr/=0) exit
                 ii = index(token,':') ! token is single entry or range of entries
                 lt = len_trim(token)
                 if (ii==0) then       ! no colon, therefore single entry
                     read(token,*) rnext
                     jcnt = jcnt + 1
                     levlst(jcnt,irow) = rnext
                 else                  ! colon, therefore kx1:kx2
                     read(token(1:ii-1),*) r1
                     read(token(ii+1:lt),*) r2
                     k1=nint(r1); k2=nint(r2)    ! levs/chanels converted to integer
                     do knext = k1, k2
                        if (jcnt==LVMAX) then    ! check space
                          write(stderr,'(2a,i5)') myname,': increase LVMAX'
                          stat = 4; return
                        else if (ierr==0) then
                          jcnt = jcnt + 1
                          levlst(jcnt,irow) = knext
                        end if
                     end do
                 end if
               enddo

           end if
        end do
        nlv = jcnt
        print *, 'Will process the following lev/chanel''s:'
        do j = 1, ncfound
           call rcount_this_ ( levlst(:,j), n )
           print *, levlst(1:n,j)
        enddo
      endif ! < all Levs >

!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine ObsOperSet_

      subroutine icount_this_ ( indx, n )
      implicit none
      integer, intent(in)  :: indx(:)
      integer, intent(out) :: n
      integer  i,m
      m=size(indx,1)
      n=0
      do i = 1, m
         if(indx(i)/=0) n=n+1 
      enddo
      end subroutine icount_this_

      subroutine rcount_this_ ( array, n )
      implicit none
      real,    intent(in)  :: array(:)
      integer, intent(out) :: n
      integer  i,m
      m=size(array,1)
      n=0
      do i = 1, m
         if(array(i)>0.0) n=n+1 
      enddo
      end subroutine rcount_this_


!EOC

      end program odsstats ! program odsstats

      subroutine no_outliers ( ods, nobs, log_transf, igood )
       
      use m_odsmeta, only: X_TOO_HIGH
      use m_ods
      implicit none
      integer, intent(inout) :: nobs
      logical, intent(in) :: log_transf
      integer, intent(inout) :: igood(nobs)
      type(ods_vect) ods
      real, parameter :: sigtol    = 4.0
      real, parameter :: aodsigtol = 3.0
      real :: diag, mean, rms, stddev
      integer i,ii, ngood

      if(nobs==0) return

      mean = 0.0; rms = 0.0; ngood = 0
!     if ( trim(attr) == 'dfs' ) then   ! calculate DFS [h(xa)-h(xb)]*oma/R
         do i=1,nobs
            diag=((ods%data%omf(i)-ods%data%oma(i))*ods%data%oma(i))/(ods%data%xvec(i))**2
            mean=mean+diag
            rms = rms+diag*diag
         enddo
!     endif
!     do ii = 1,nobs    
!         mean = mean + ods%data%omf(ii)
!         rms  = rms  + ods%data%omf(ii)*ods%data%omf(ii)
!     enddo
      mean = mean / nobs
      rms  = rms  / nobs
      stddev = sqrt(nobs*(abs(rms - mean*mean))/(nobs-1))    
   
      ngood = 0
      igood = 0
      do ii=1,nobs
          diag=((ods%data%omf(ii)-ods%data%oma(ii))*ods%data%oma(ii))/(ods%data%xvec(ii))**2
         if (ods%data%kt(ii) == 43) then 
            if(    ods%data%qcexcl(ii)==0  .and.
!    .         (abs(ods%data%omf(ii))>aodsigtol*ods%data%xvec(ii) .or.
!    .          abs(ods%data%oma(ii))>aodsigtol*ods%data%xvec(ii)) ) then
!    .         (abs(ods%data%omf(ii))>aodsigtol*stddev .or.
!    .          abs(ods%data%oma(ii))>aodsigtol*stddev) ) then
     .          abs(diag)>aodsigtol ) then
               ods%data%qcexcl(ii) = X_TOO_HIGH
            else
               ngood=ngood+1
               igood(ngood) = ii
            endif
         else
            if(    ods%data%qcexcl(ii)==0  .and.
!    .         (abs(ods%data%omf(ii))>sigtol*ods%data%xvec(ii) .or.
!    .          abs(ods%data%oma(ii))>sigtol*ods%data%xvec(ii)) ) then
!    .         (abs(ods%data%omf(ii))>sigtol*stddev .or.
!    .          abs(ods%data%oma(ii))>sigtol*stddev) ) then
     .          abs(diag)>sigtol ) then
               ods%data%qcexcl(ii) = X_TOO_HIGH
             else
               ngood=ngood+1
               igood(ngood) = ii
             endif
          endif
      enddo
      if(ngood<nobs)
     .  print *, 'Number of obs and outliers found: ', nobs, nobs-ngood
      nobs=ngood

      end subroutine no_outliers

      subroutine aodfix_sigo ( ods, log_transf )
       
      use m_ods
      implicit none
      logical,intent(in) :: log_transf
      type(ods_vect) ods
      real :: aodeps = 0.01
      
      if (log_transf) then
         where (ods%data%kt == 43 .and. ods%data%qcexcl==0) 
             ods%data%xvec = exp(0.2)
         endwhere
         where (ods%data%kt == 43 .and. ods%data%qcexcl==0) 
             ods%data%oma = max(0.0,exp(ods%data%obs-ods%data%oma) - aodeps)
             ods%data%omf = max(0.0,exp(ods%data%obs-ods%data%omf) - aodeps)
             ods%data%obs = max(0.0,exp(ods%data%obs))
             ods%data%omf = ods%data%obs-ods%data%omf
             ods%data%oma = ods%data%obs-ods%data%oma
         endwhere
      else
         where (ods%data%kt == 43 .and. ods%data%qcexcl==0) 
             ods%data%xvec = 0.2
         endwhere
      endif

      end subroutine aodfix_sigo

      subroutine sumattr ( verb, ods, nobs, igood, attr, rsum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      integer, intent(in) :: igood(nobs)
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  rsum
      logical, intent(in) :: verb
      integer  i,ii

      rsum = 0.0
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'obs' ) then
         do i=1,nobs
            rsum=rsum+ods%data%obs(i)
         enddo
      endif
      if ( trim(attr) == 'omf' ) then
         do i=1,nobs
            rsum=rsum+ods%data%omf(i)
         enddo
      endif
      if ( trim(attr) == 'bbcomf' ) then
         do i=1,nobs
            rsum=rsum+ods%data%omf(i)+ods%data%xm(i)
         enddo
      endif
      if ( trim(attr) == 'oma' ) then
         do i=1,nobs
            rsum=rsum+ods%data%oma(i)
         enddo
      endif
      if ( trim(attr) == 'xm' ) then
         do i=1,nobs
            rsum=rsum+ods%data%xm(i)
         enddo
      endif
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
            rsum=rsum+ods%data%xvec(i)
         enddo
      endif
      if ( trim(attr) == 'imp0hr' ) then ! calculate 0hr impact as Jo(a)-Jo(b)
         do i=1,nobs
            rsum=rsum+((ods%data%oma(i))**2  -
     .                 (ods%data%omf(i))**2) /(ods%data%xvec(i))**2 ! calculate impact
         enddo
      endif
      if ( trim(attr) == 'dfs' ) then   ! calculate DFS [h(xa)-h(xb)]*oma/R
         do ii=1,nobs
            i=igood(ii)
            if (i==0) cycle
            rsum=rsum+((ods%data%omf(i)-ods%data%oma(i))*ods%data%oma(i))/(ods%data%xvec(i))**2 ! calculate dfs
         enddo
      endif
      if ( trim(attr) == 'imp0hrxm' ) then ! calculate 0hr impact as Jo(a)-Jo(b) (for when sigo stored in xm)
         do i=1,nobs
            rsum=rsum+((ods%data%oma(i))**2  -
     .                 (ods%data%omf(i))**2) /(ods%data%xm(i))**2 ! calculate impact
         enddo
      endif
      if ( trim(attr) == 'imp_from_sens' .or.
     .     trim(attr) == 'xvecxomf'    ) then  ! this handles the case when xvec holds sensitivities
                                               ! rather than impacts themselves
         do i=1,nobs
            rsum=rsum+ods%data%xvec(i)*ods%data%omf(i) ! calculate impact
         enddo
      endif
      if ( trim(attr) == 'xvecxamo' ) then  ! this handles the case when xvec holds sensitivities
                                            ! rather than impacts themselves
         do i=1,nobs
            rsum=rsum-ods%data%xvec(i)*ods%data%oma(i) ! calculate impact to sigO
         enddo
      endif
      if ( trim(attr) == 'omfxomf' ) then
         do i=1,nobs
            rsum=rsum+ods%data%omf(i)*ods%data%omf(i)
         enddo
      endif
      if ( trim(attr) == 'omaxoma' ) then
         do i=1,nobs
            rsum=rsum+ods%data%oma(i)*ods%data%oma(i)
         enddo
      endif
      if ( trim(attr) == 'xmxxm' ) then
         do i=1,nobs
            rsum=rsum+ods%data%xm(i)*ods%data%xm(i)
         enddo
      endif
      if ( trim(attr) == 'omf2byxvec2' .or. trim(attr) == 'job' ) then  ! Jo(b) when sigO in xvec slot
         do i=1,nobs
            rsum=rsum+0.5*ods%data%omf(i)*ods%data%omf(i)/(ods%data%xvec(i)*ods%data%xvec(i))
         enddo
      endif
      if ( trim(attr) == 'oma2byxvec2' .or. trim(attr) == 'joa' ) then  ! Jo(a) when sigO in xvec slot
         do i=1,nobs
            rsum=rsum+0.5*ods%data%oma(i)*ods%data%oma(i)/(ods%data%xvec(i)*ods%data%xvec(i))
         enddo
      endif
      if ( trim(attr) == 'omf2byxm2' ) then    ! Jo(b) when sigO in xm slot
         do i=1,nobs
            rsum=rsum+ods%data%omf(i)*ods%data%omf(i)/(ods%data%xm(i)*ods%data%xm(i))
         enddo
      endif
      if ( trim(attr) == 'oma2byxm2' ) then    ! Jo(a) when sigO in xm slot
         do i=1,nobs
            rsum=rsum+ods%data%oma(i)*ods%data%oma(i)/(ods%data%xm(i)*ods%data%xm(i))
         enddo
      endif

      if ( verb ) then
          print *, 'Obs impact(sum):      ', rsum
          print *, 'Obs impact(sum)/Nobs: ', rsum/nobs
      endif

      end subroutine sumattr

      subroutine negattr ( verb, ods, nobs, igood, attr, rsum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      integer, intent(in) :: igood(nobs)
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  rsum
      logical, intent(in) :: verb
      integer  i,ii
      real     imp
      real,parameter:: eps=0.0

      rsum = 0.0
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
            if(ods%data%xvec(i)<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp0hr' ) then
         do i=1,nobs
            imp = (ods%data%oma(i))**2 -
     .            (ods%data%omf(i))**2
            imp = imp/ods%data%xvec(i)**2
            if(imp<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'dfs' ) then ! note: this is here for completeness; typically DFS>0
         do ii=1,nobs
            i=igood(ii)
            if (i==0) cycle
            imp = (ods%data%omf(i)-ods%data%oma(i))*ods%data%oma(i)
            imp = imp/ods%data%xvec(i)**2
            if(imp>eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp0hrxm' ) then
         do i=1,nobs
            imp = (ods%data%oma(i))**2 -
     .            (ods%data%omf(i))**2
            imp = imp/ods%data%xm(i)**2
            if(imp<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp_from_sens' .or.
     .     trim(attr) == 'xvecxomf'    ) then  ! this handles the case when xvec holds sensitivities
         do i=1,nobs
            if(ods%data%xvec(i)*ods%data%omf(i)<eps) rsum=rsum+1.0
         enddo
      endif
      if ( verb ) then
          print *, '    Number of obs with positive impact:  ',  rsum
          print *, 'Percentage of obs with positive impact:  ', (rsum/nobs)*100.
      endif

      end subroutine negattr

      subroutine neuattr ( verb, ods, nobs, igood, attr, rsum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      integer, intent(in) :: igood(nobs)
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  rsum
      logical, intent(in) :: verb
      integer  i
      real     imp
      real,parameter:: eps=1.0e-10

      rsum = 0.0
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
            if(abs(ods%data%xvec(i))<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp0hr' ) then
         do i=1,nobs
            imp = (ods%data%oma(i))**2 -
     .            (ods%data%omf(i))**2
            imp = imp/ods%data%xvec(i)*2
            if(abs(imp)<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp0hrxm' ) then
         do i=1,nobs
            imp = (ods%data%oma(i))**2 -
     .            (ods%data%omf(i))**2
            imp = imp/ods%data%xm(i)**2
            if(abs(imp)<eps) rsum=rsum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp_from_sens' .or.
     .     trim(attr) == 'xvecxomf'    ) then  ! this handles the case when xvec holds sensitivities
         do i=1,nobs
            if(abs(ods%data%xvec(i)*ods%data%omf(i))<eps) rsum=rsum+1.0
         enddo
      endif
      if ( verb ) then
          print *, '    Number of obs with positive impact:  ',  rsum
          print *, 'Percentage of obs with positive impact:  ', (rsum/nobs)*100.
      endif

      end subroutine neuattr

      subroutine gmsattr ( verb, ods, nobs, iopt, igood, attr, rsum )
      use m_ods
      use m_die, only: die
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      integer, intent(in) :: iopt  ! 1=rms; 2=stdev
      integer, intent(in):: igood(nobs)
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  rsum
      logical, intent(in) :: verb

      character(len=*), parameter :: myname_ = "gmsattr"
      real dfs,mean
      integer  i,ii

      rsum = 0.0; mean = 0.0
      if (nobs==0) return ! nothing to do

      if (iopt==2) then
         select case (trim(attr))
          case ( 'obs' )
            do i=1,nobs
               mean=mean+ods%data%obs(i)
            enddo
          case ( 'omf' )
            do i=1,nobs
               mean=mean+ods%data%omf(i)
            enddo
          case ( 'bbcomf' )
            do i=1,nobs
               mean=mean+ods%data%omf(i)+ods%data%xm(i)
            enddo
          case ( 'oma' )
            do i=1,nobs
               mean=mean+ods%data%oma(i)
            enddo
          case ( 'xm' )
            do i=1,nobs
               mean=mean+ods%data%xm(i)
            enddo
          case ( 'dfs' )
            do ii=1,nobs
               i=igood(ii)
               if (i==0) cycle
               dfs = ((ods%data%omf(i)-ods%data%oma(i))*ods%data%oma(i))/(ods%data%xvec(i))**2
               mean=mean+dfs
            enddo
          case ( 'xvec' )
            do i=1,nobs
               mean=mean+ ods%data%xvec(i)
            enddo
          case default
            call die(myname_,'Cannot find attr ' // trim(attr) // ' aborting ...',99)
         end select
         mean = mean/nobs
      endif ! mean calculation done

      select case (trim(attr))
       case ( 'obs' )
         do i=1,nobs
            rsum=rsum+(ods%data%obs(i)-mean)**2
         enddo
       case ( 'omf' )
         do i=1,nobs
            rsum=rsum+(ods%data%omf(i)-mean)**2
         enddo
       case ( 'bbcomf' )
         do i=1,nobs
            rsum=rsum+(ods%data%omf(i)+ods%data%xm(i)-mean)**2
         enddo
       case ( 'oma' )
         do i=1,nobs
            rsum=rsum+(ods%data%oma(i)-mean)**2
         enddo
       case ( 'xm' )
         do i=1,nobs
            rsum=rsum+(ods%data%xm(i)-mean)**2
         enddo
       case ( 'dfs' )
         do ii=1,nobs
            i=igood(ii)
            if(i==0) cycle
            dfs = ((ods%data%omf(i)-ods%data%oma(i))*ods%data%oma(i))/(ods%data%xvec(i))**2
            rsum=rsum+(dfs-mean)*(dfs-mean)
         enddo
       case ( 'xvec' )
         do i=1,nobs
!           rsum=rsum+(ods%data%xvec(i))**2
            rsum=rsum+(max(1.e-15,abs(ods%data%xvec(i)-mean)))**2 ! hack: "-fp-model precise" flag
                                                                  ! of intel compiler version 11 gives
                                                                  ! underflow here in some cases - non-sense
         enddo
       case default
         call die(myname_,'Cannot find attr ' // trim(attr) // ' aborting ...',99)
      end select
      if ( verb ) then
          print *, 'MS of impacts (', trim(attr), '): ', rsum
      endif

      if (iopt==2) then
         if(nobs>1) rsum = sqrt(rsum/(nobs-1.0))
      endif
      end subroutine gmsattr

      subroutine Write_Stats ( outfile, obssum, obsnum, obserr, oclass, ncfound, ncmax, nt,
     .                         nymd, nhms, verb, jedi, stat ) 
      use m_ioutil, only : luavail
      implicit none
      integer, intent(in)          :: ncfound, ncmax, nt
      integer, intent(in)          :: nymd(nt), nhms(nt)
      character(len=*), intent(in) :: outfile
      character(len=*), intent(in) :: oclass(ncfound)
      real,    intent(in)          :: obssum(ncmax,nt)
      real,    intent(in)          :: obserr(ncmax,nt)
      integer, intent(in)          :: obsnum(ncmax,nt)
      logical, intent(in)          :: verb
      logical, intent(in)          :: jedi
      integer, intent(out)         :: stat

      integer i,j,ios,lu

      stat = 0
      lu=luavail()
      open(lu, file=outfile, form='formatted', iostat=ios)
      if ( ios .ne. 0 ) then
           print *, 'Cannot open file ', trim(outfile)
           stat = 1
           return
      else
           if(verb) print *, 'Writing ASC file ', trim(outfile)
      end if

      do j = 1, nt ! loop over times
         do i = 1, ncfound  ! loop ob classes
            if ( jedi ) then
              if ( obsnum(i,j) > 0 ) then
!                write(lu,'(i8.8,1x,i6.6,1x,3a,1x,1p,e11.4,a,i7,1x,a,f8.5,a,e11.4)') 
                 write(lu,'(i8.8,1x,i6.6,1x,3a,1x,1p,e11.4,a,i7,1x,a,e11.4,a,e11.4)') 
     .                 nymd(j), nhms(j), 
     .                 ', Jo(', trim(oclass(i)), ') =', obssum(i,j),
     .                 ', nobs = ', obsnum(i,j),
     .                 ', Jo/n = ', obssum(i,j)/obsnum(i,j),
     .                 ', err = ',  obserr(i,j)/obsnum(i,j)
              endif
            else            ! intentionally truncate oclass to 20 chars
              write(lu,'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4)') 
     .              nymd(j), nhms(j), oclass(i), obsnum(i,j), obssum(i,j)
            endif
         enddo
      enddo
      close(lu)

      end subroutine Write_Stats

      subroutine Write_gsiRfactor ( RCfile, outfile, obssum, oclass, ncfound, ncmax, nt, 
     .                              nymd, nhms, verb, stat ) 

      use m_ioutil, only: luavail
      use m_chars,  only: lowercase
      use m_stdio, only : stdout,stderr
      use m_inpak90
      use m_die, only: die
      implicit none
      integer, intent(in)          :: ncfound, ncmax, nt
      integer, intent(in)          :: nymd(nt), nhms(nt)
      character(len=*), intent(in) :: RCfile
      character(len=*), intent(in) :: outfile
      character(len=*), intent(in) :: oclass(ncfound)
      real,    intent(in)          :: obssum(ncmax,nt)
      logical, intent(in)          :: verb
      integer, intent(out)         :: stat

      character(len=*), parameter :: myname = 'Write_gsiRfactor'
      character(len=255) fname, token, obsoperrc
      integer i,j,ios,lu,nc,iret
      integer igsiRfacPrec  ! precision to writeout GSI "Rcov"
      integer,allocatable :: ichan(:)
      real    rn, alpha
      real    snrm2
      real,   allocatable :: osum_nrmzd(:)
      real(4),allocatable :: R4(:,:)
      real(8),allocatable :: R8(:,:)

      stat = 0

!     Load resources file
!     -------------------
      obsoperrc = ' '
      call getenv('obs_opers.rc',OBSOPERRC)     ! Unix binding
      if(obsoperrc.eq.' ') obsoperrc=RCfile     ! default name
      call i90_loadf (trim(obsoperrc), iret)
      if( iret .ne. 0) then
          write(stderr,'(2a,i5)') myname,': I90_loadf error, iret =',iret
          stat = 1
          return
      end if

!     Inquire about norm only
!     ---------------------------------------------
      rn = -1.0
      call I90_label('GSI*Rsens*Norm:', iret)
      if (iret .eq. 0) then
          call I90_Gtoken ( token, iret )
          if ( iret==0 ) then
               read(token,*) rn
          endif
      else
          call I90_release()
          return   ! nothing to do
      end if

!     In case doing more than simply recovering norm of sensitivity vector ...
!     ------------------------------------------------------------------------
      if ( rn > 0.0 ) then

!         Read scaling factor for sensitivities to sigO
!         ---------------------------------------------
          call I90_label('GSI*Rsens*Scale:', iret)
          if (iret .eq. 0) then
              call I90_Gtoken ( token, iret )
              if ( iret==0 ) then
                   read(token,*) alpha
              endif
          else
              call I90_release()
              return   ! nothing to do
          end if
          write(stdout,'(a,1p,e9.2)') ' Will scale R-sensitivities with this factor: ',alpha

!         Read option to write out GSI-ready R matrix re-scaling factor
!         -------------------------------------------------------------
          igsiRfacPrec = 8  ! 4 = real(4)
                            ! 8 = real(8) (used for testing)
          call I90_label('GSI*Rfactor*Precision:', iret)
          if (iret .eq. 0) then
              call I90_Gtoken ( token, iret )
              if ( iret==0 ) then
                   read(token,*) igsiRfacPrec
              endif
          end if
          write(stdout,'(a,i3)') ' Will write out GSI R-scaling factor based on sens, prec= ',igsiRfacPrec

      endif ! norm check
    
!     release resource file:
!     ---------------------
      call I90_release()

!     Handle R matrix
!     ---------------
      lu=luavail()
      do j = 1, nt ! loop over times

         nc=len_trim(oclass(1))

         ! Normalize sensitivity vector
         if (rn<0.0) then
             rn = dot_product(obssum(:,j),obssum(:,j))
             write(stdout,'(a,a10,a,1p,e10.3)') ' Obsclass: ', trim(oclass(1)(1:nc-3)) , 
     .                                          ' norm square of sensitivity vector: ', rn  
             return ! all done
         endif
         allocate(osum_nrmzd(ncmax))
         osum_nrmzd = obssum(:,j)/rn
    
         ! open file for this time ...
         write(fname,'(2a,i8.8,a,i2.2,a)') trim(outfile), '.', nymd(j), '_', nhms(j)/10000, 'z.bin'
         open(lu, file=trim(fname), form='unformatted', convert='little_endian', iostat=ios)
         if ( ios .ne. 0 ) then
              print *, 'Cannot open file ', trim(fname)
              stat = 1
              return
         else
              if(verb) print *, 'Writing BIN file ', trim(fname)
         end if

         ! extract last three digits in name class (supposedly channel indexes)
         allocate(ichan(ncfound))
         do i=1,ncfound
            read(oclass(i)(nc-2:nc),'(i3)') ichan(i)
         enddo

         ! create R matrix and copy factors to its diagonal
         write(lu) ncfound, igsiRfacPrec
         write(lu) ichan
         if (igsiRfacPrec==4) then
             allocate(R4(ncfound,ncfound))
             R4=0.0
             do i=1,ncfound
                R4(i,i)=1.0-2.0*alpha*osum_nrmzd(i)
                if ( R4(i,i)<0.0 ) then
                   R4(i,i)=0.5
                endif
                R4(i,i)=R4(i,i)*R4(i,i)
                if(verb) print *, R4(i,i)
             enddo
             write(lu) R4
             deallocate(R4)
         endif
         if (igsiRfacPrec==8) then
             allocate(R8(ncfound,ncfound))
             R8=0.d0
             do i=1,ncfound
                R8(i,i)=1.d0-2.d0*alpha*osum_nrmzd(i)
                if ( R8(i,i)<0.d0 ) then
                   R8(i,i)=0.5d0
                endif
                R8(i,i)=R8(i,i)*R8(i,i)
                if(verb) print *, R8(i,i)
             enddo
             write(lu) R8
             deallocate(R8)
         endif
         deallocate(ichan)

         ! wrap it up
         close(lu)
         deallocate(osum_nrmzd)
      enddo

      end subroutine Write_gsiRfactor

      subroutine Write_AccumStats ( fname, obssum, obsgms, obsnum, obsneg, obsneu, oclass, ncfound, ncmax, nymd, nhms, lrms, 
     .                              verb, stat )
      use m_ioutil, only : luavail
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in)          :: ncfound, ncmax
      character(len=*), intent(in) :: oclass(ncfound)
      real,    intent(in)          :: obssum(ncmax)
      real,    intent(in)          :: obsgms(ncmax)
      integer, intent(in)          :: obsnum(ncmax)
      integer, intent(in)          :: obsneg(ncmax)
      integer, intent(in)          :: obsneu(ncmax)
      integer, intent(in)          :: nymd, nhms
      logical, intent(in)          :: lrms
      logical, intent(in)          :: verb
      integer, intent(out)         :: stat

      integer,allocatable :: lu(:)
      integer i,j,ios,ll,n,nu
      real    xmean,rmean,tmp,xstdv

      stat = 0
      nu = 1
      if(verb) nu=2
      allocate(lu(nu))
      lu(1)=luavail()
      if(verb) lu(2)=6
      open(lu(1), file=fname, form='formatted', iostat=ios)

!     Write out accumulated results
!     -----------------------------
      if (lrms) then
         do ll=1,nu
            write(lu(ll),'(a,1p,e11.4,25x,5a)') '#Total: ', sum(obssum), ' obs ', '        sum ', 
     .                                               '     numneg ', '     numneu ', '      rms'
         enddo
         do i = 1, ncfound  ! loop ob classes
            n = obsnum(i)
                            ! intentionally truncate oclass to 20 chars
            do ll=1,nu
            write(lu(ll),'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4,1x,i11,1x,i11,1x,1p,e11.4)') nymd, nhms, 
     .                                       oclass(i), obsnum(i), obssum(i), obsneg(i), obsneu(i), obsgms(i)
            enddo
         enddo
      else
         do ll=1,nu
            write(lu(ll),'(a,1p,e11.4,25x,4a)') '#Total: ', sum(obssum), ' obs ', ' sum_impact ', 
     .                                          '     numneg ', '     numneu '
            do i = 1, ncfound  ! loop ob classes
                               ! intentionally truncate oclass to 20 chars
               write(lu(ll),'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4,1x,i11,1x,i11)') 
     .               nymd, nhms, oclass(i), obsnum(i), obssum(i), obsneg(i), obsneu(i)
         enddo
      enddo
      endif

      close(lu(1))
      deallocate(lu)

      end subroutine Write_AccumStats

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odstats
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, nfiles_max, nfiles, RCfile, outfile, 
     .                  trange, latrange, lonrange, verb, jedi, aodsigofix, ncf )

! !USES:

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max

! !OUTPUT PARAMETERS:

      character(len=*), intent(out) :: RCfile
      character(len=*), intent(out) :: infile(nfiles_max)
      integer,          intent(out) :: nfiles
      character(len=*), intent(out) :: outfile
      integer,          intent(out) :: trange(2)
      real,             intent(out) :: latrange(2)
      real,             intent(out) :: lonrange(2)
      logical,          intent(out) :: jedi
      logical,          intent(out) :: verb
      logical,          intent(out) :: aodsigofix
      logical,          intent(out) :: ncf
!
!
! !REVISION HISTORY:
!     22Jan2008 Todling - Initial code (stripped off odsselect)
!     19Feb2008 Todling - Add trange
!     16Apr2008 Todling - Add lat/lon ranges
!
!EOP
!BOC

      character*4, parameter :: myname_ = 'init'

      integer iret, i, ic, lt, lv, iarg, argc, iargc
      real swap
      character*255 argv
      character*10 SS

      RCfile  = 'odsstats.rc'
      outfile = 'odsstats'
      trange  = (/-9999,+9999/)   ! Default: includes obs in all time ranges
      latrange  = (/-90.,+90./)   ! Default: includes obs in all latitude ranges
      lonrange  = (/-180.,+180./) ! Default: includes obs in all longitude ranges 
      verb = .false.
      jedi = .false.
      aodsigofix = .false.
      ncf  = .false.

!     Parse command line
!     ------------------

      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
      nfiles = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, outfile )
         else if (index(argv,'-verbose' ) .gt. 0 ) then
            verb = .true.
         else if (index(argv,'-jediformat' ) .gt. 0 ) then
            jedi = .true.
         else if (index(argv,'-ncf' ) .gt. 0 ) then
            ncf = .true.
         else if (index(argv,'-aodsigofix' ) .gt. 0 ) then
            aodsigofix = .true.
         else if (index(argv,'-rc' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, RCfile )
         elseif (index(argv,'-time') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) trange(1)
                trange(2) = trange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) trange(1)
                read(SS(ic+1:lt),*) trange(2)
            end if
            if (trange(2)<trange(1)) then    ! let's be nice to the user..
                swap = trange(2)
                trange(2) = trange(1)
                trange(1) = swap
            end if
         elseif (index(argv,'-lat') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) latrange(1)
                latrange(2) = latrange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) latrange(1)
                read(SS(ic+1:lt),*) latrange(2)
            end if
            if (latrange(2)<latrange(1)) then    ! let's be nice to the user..
                swap = latrange(2)
                latrange(2) = latrange(1)
                latrange(1) = swap
            end if
         elseif (index(argv,'-lon') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) lonrange(1)
                lonrange(2) = lonrange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) lonrange(1)
                read(SS(ic+1:lt),*) lonrange(2)
            end if
            if (lonrange(2)<lonrange(1)) then    ! let's be nice to the user..
                swap = lonrange(2)
                lonrange(2) = lonrange(1)
                lonrange(1) = swap
            end if
         else
            nfiles = nfiles + 1
            if ( nfiles .gt. nfiles_max ) then
               print *, 'Maximum number of input files = ', nfiles_max
               stop
            end if
            infile(nfiles) = argv
         end if
      end do
 111  continue
      if ( nfiles .lt. 1 ) call usage()

      print *
      print *, 'Input files: ', nfiles
      print *
      do i = 1, nfiles
         lv = len_trim(infile(i))
         print *, ' o ', infile(i)(1:lv)
      end do
      print *
      print *, 'Output filename: ', trim(outfile)

      return

      end ! subroutine init

!EOC

!-------------------------------------------------------------------------

      subroutine usage()
      print *
      print *, 'Usage:'
      print *
      print *, 'obsstats [-o ID] -rc RCfile odsfile(s)'
      print *
      print *, 'where'
      print *
      print *,'-o  ID         use ID for naming output files'
      print *,'                (default: obsstat.txt)'
      print *,'-aodsigofix    wire sigO for AOD residuals (see Notes)'
      print *,'-jediformat    format output in similar as JEDI'
      print *,'-verbose       sets verbose on (default: off)'
      print *,'-rc RCfile      resource file'
      print *,'                (default: obsstat.rc)'
      print *,'-time ti:tf    specify window of time to select obs from'
      print *,'                (default: -999:+999, all)'
      print *,'-lat  latA:latB specify latitude range'
      print *,'                (default: -90:+90, all)'
      print *,'-lon  lonA:lonB specify longitude range'
      print *,'                (default: -180:+180, all)'
      print *
      print *
      print *,' odsfile(s)    ODS file(s)'
      print *
      print *, 'NOTES: '
      print *, '-----  '
      print *
      print *, ' The following are entries allowed in the rc file'
      print *
      print *, '  1. Known OBS*Operations: (default: sum)'
      print *, '     sum     - sum all variables '
      print *, '     rms     - root mean square '
      print *, '     stddev  - stdandard deviation from mean of variables '
      print *, '     numneg  - sum all non-negative variables '
      print *
      print *, '  2. Known OBS*Variable: (default: xvec)'
      print *, '     obs,omf,oma,xm,xvec,xvecxomf,xvecxamo '
      print *, '     2a. xvecxomf,xvecxamo are only meaningful when xvec holds '
      print *, '         the sensitivities instead of the impacts.'
      print * 
      print *, '  3. Earlier versions of the DAS did not fill in the xvec(sigO) slot'
      print *, '     from AOD ODS files. This knob allows for calculation of DFS and '
      print *, '     IMP0HR to take place by having the present code wired the typical'
      print *, '     (constant) value used in PSAS.'
      print *
      stop
      end
