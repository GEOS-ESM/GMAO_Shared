      program dyndiff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: dyndiff: find the dynamic vector difference between two hdf files
!
! !USAGE: see the routine usage() below
!
! !USES:
!
      use m_dyn
      use m_dyn_util, only: Dyn_Util_Tv2T
      use m_dyn_util, only: Dyn_Scale_by_TotEne
      use m_dyn_util, only: Dyn_Get_Energy

      implicit NONE

! !DESCRIPTION: Uses statistic subroutine dyn_stat to check the dynamic
!               vector dfference between two alike hdf files.
!
! !REVISION HISTORY:
!
!  2002.01.07  E. Yeh:    Initial code.
!  2005.03.28  Elena N.:  Added an option to save resulting dynamic vector in a file.
!  01May2007   Todling    Support to new GEOS-5 dyn-vector files
!  06May2007   Todling    Add opt to print ave/min/max when one file is given as input.
!  05Mar2009   Todling    Add land fractions
!  08Mar2012   Todling    Add acoeff [to allow calc of inc as -(b-a)]
!  12Sep2012   Todling    Quick fix for LM issue in dyn-vector
!  10Jun2020   Todling    Add ability to add rh diff and mean rh to diff file
!
!-------------------------------------------------------------------------
!EOP

      character(len=*), parameter :: myname = 'dyndiff'

!     File names
!     ----------
      integer, parameter :: MFILES = 2 ! max.   number of input files
      character(len=255) :: files(MFILES), binfiles(MFILES), etafile, st
      integer            :: nfiles       ! actual no. of input files

      character(len=255) :: dyn_dout      ! difference output file name
      logical, parameter :: fix_top_rh = .true.

!     Dynamics/simulator vectors
!     --------------------------
      type(dyn_vect) dyn(MFILES)   ! dynamics vector in eta

!     Locals
!     ------
      character(len=255) :: egress
      character(len=20)  :: anorm
      character(len=255) :: jnorm
      integer, parameter :: READ_ONLY = 1
      integer fid, nvars, ngatts
      integer ios, rc, iopt, ifile
      integer ntimes, n, freq, nymd, nhms, prec
      integer freq_d, nymd_d, nhms_d, prec_d    !Timetag for newly created diff in *.hdf format  
      integer im, jm, km, lm, system, dyntype, irh
      logical dominmax,verb,sbyene,tv2t
      integer addrh
      integer vnorm
      logical normlz
      character(len=3) ntype ! norm type (when applicable)
      real, allocatable :: ps  (:,:)
      real, allocatable :: delp(:,:,:)
      real, allocatable :: rh1(:,:,:)
      real, allocatable :: rh2(:,:,:)
      real    acoeff
      real    eps_eer
      real    projlat(2), projlon(2)
      integer projlev(2)
      
!  Initialize
!  ----------     
   call Init_ ( dyntype, mfiles, files, dominmax, verb, egress, eps_eer, anorm, jnorm,  &
                tv2t, projlon, projlat, projlev, normlz, ntype, addrh, vnorm )

!  Loop over input eta files
!  -------------------------
   nfiles = 1
   do ifile = 1, nfiles

      etafile = files(ifile)

!     Determine how many time levels on file
!     --------------------------------------
      call GFIO_Open ( etafile, READ_ONLY, fid, rc )
      if ( rc .ne. 0 ) then
         call die(myname,'cannot open GFIO file '//trim(etafile))
      end if
      call GFIO_DimInquire ( fid, im, jm, km, ntimes, nvars, ngatts, rc)
      if ( rc .ne. 0 ) then
         call die(myname,'problems getting dimensions' )
      end if
      call GFIO_Close ( fid, rc )
      
!     For each time on file...
!     ------------------------
      do n = 1, ntimes

!        Get ETA data for this time
!        --------------------------
         call dyn_get ( etafile, nymd, nhms, dyn(1), rc, timidx=n, freq=freq, vectype=dyntype )
         nymd_d = nymd     !Newly created diff file with nymd from first file dyn(1)
         nhms_d = nhms     !Newly created diff file with nhms from first file dyn(1)

         if ( rc .ne. 0 ) then
            call die(myname,'cannot read dynamics vector file')
         end if
         call dyn_get ( files(2), nymd, nhms, dyn(2), rc, timidx=n, freq=freq, vectype=dyntype )
         if ( rc .ne. 0 ) then
            call die(myname,'cannot read dynamics vector file')
         end if

         ! check dims
         if ( dyn(1)%grid%im/=dyn(2)%grid%im .or. &
              dyn(1)%grid%jm/=dyn(2)%grid%jm .or. &
              dyn(1)%grid%km/=dyn(2)%grid%km ) then
            write(6,'(a,3(i6,2x))') 'dyn-1:', dyn(1)%grid%im,dyn(1)%grid%jm,dyn(1)%grid%km
            write(6,'(a,3(i6,2x))') 'dyn-2:', dyn(2)%grid%im,dyn(2)%grid%jm,dyn(2)%grid%km
            call die(myname,'error, incompatible dims')
         endif
         print *, "> nymd, nhms: ", nymd, nhms, " (diff)"
         lm = min(dyn(1)%grid%lm,dyn(2)%grid%lm)
         if ( .not. dominmax ) then
           if (abs(addrh)>0) then
              allocate(rh1(im,jm,km))
              allocate(rh2(im,jm,km))
              call getrh_(rh1,dyn(1)%pt,dyn(1)%q(:,:,:,1),dyn(1)%ps,dyn(1)%grid%ak,dyn(1)%grid%bk)
              call getrh_(rh2,dyn(2)%pt,dyn(2)%q(:,:,:,1),dyn(2)%ps,dyn(2)%grid%ak,dyn(2)%grid%bk)
           endif

           print *, "scaling difference by: ", acoeff
           if (sbyene) then
              allocate(ps  (dyn(1)%grid%im,dyn(1)%grid%jm))
              allocate(delp(dyn(1)%grid%im,dyn(1)%grid%jm,dyn(1)%grid%km))
              ps  =dyn(1)%ps
              delp=dyn(1)%delp
           endif
           dyn(1)%ps     = acoeff*(dyn(1)%ps - dyn(2)%ps)
           dyn(1)%ts     = acoeff*(dyn(1)%ts - dyn(2)%ts)
           dyn(1)%phis   = acoeff*(dyn(1)%phis - dyn(2)%phis)
           dyn(1)%lwi    = acoeff*(dyn(1)%lwi - dyn(2)%lwi)
           dyn(1)%frland = acoeff*(dyn(1)%frland - dyn(2)%frland)
           dyn(1)%frlandice = acoeff*(dyn(1)%frlandice - dyn(2)%frlandice)
           dyn(1)%frlake =  acoeff*(dyn(1)%frlake - dyn(2)%frlake)
           dyn(1)%frocean=  acoeff*(dyn(1)%frocean - dyn(2)%frocean)
           dyn(1)%frseaice= acoeff*(dyn(1)%frseaice - dyn(2)%frseaice)
           dyn(1)%hs_stdv = acoeff*(dyn(1)%hs_stdv - dyn(2)%hs_stdv)
           dyn(1)%delp    = acoeff*(dyn(1)%delp - dyn(2)%delp)
           dyn(1)%u       = acoeff*(dyn(1)%u - dyn(2)%u)
           dyn(1)%v       = acoeff*(dyn(1)%v - dyn(2)%v)
           if (tv2t) then
              ! convert virtual temperature to temperature
              call Dyn_Util_Tv2T (dyn(1)%pt,dyn(1)%q(:,:,:,1))
              call Dyn_Util_Tv2T (dyn(2)%pt,dyn(2)%q(:,:,:,1))
           endif
           dyn(1)%pt      = acoeff*(dyn(1)%pt - dyn(2)%pt)
           dyn(1)%q(:,:,:,1:lm) = acoeff*(dyn(1)%q(:,:,:,1:lm) - dyn(2)%q(:,:,:,1:lm))
           if (sbyene) then
               call Dyn_Scale_by_TotEne(dyn(1),eps_eer,anorm,jnorm,projlon,projlat,projlev, &
                                        nymd,nhms,ntype=ntype,vnorm=vnorm,normlz=normlz,&
                                        ps=ps,delp=delp)
               deallocate(ps,delp)
               call Dyn_Get_Energy (dyn(1), nymd, nhms )
           endif
           if (abs(addrh)>0) then
               if (addrh<0) then
                  rh2 = (rh1 + rh2)
                  rh1 = acoeff*(2.*rh1 - rh2)
                  rh2 = 0.5*acoeff*rh2
               else
                  rh1 = acoeff*(rh1 - rh2)
               endif
           endif
         endif

!       If so, echo result to standard out
!       ----------------------------------
         if (verb) then
            call dyn_stat(6, dyn(1), rc)
            if ( rc .ne. 0 ) then
               call die(myname,'cannot process dyn_stat')
            endif
         endif

!       Do some cleaning
!       ----------------
        call dyn_clean ( dyn(2) )

!       If requested write *.hdf file with a header from dyn(1)
!       -------------------------------------------------------
        if ( trim(dyn_dout) .ne. 'NONE' ) then
             if ( abs(addrh)>0 ) then
                irh=1
                if(addrh<0) irh=2
                call dyn_init ( dyn(1), dyn(2), rc, copy=.true., vectype=dyntype, lm=dyn(1)%grid%lm+irh )
                dyn(2)%q(:,:,:,dyn(1)%grid%lm+1) = rh1
                dyn(2)%qm(dyn(1)%grid%lm+1)%name = 'rh';  dyn(2)%qm(dyn(1)%grid%lm+1)%long_name = 'Relative Humidity '
                dyn(2)%qm(dyn(1)%grid%lm+1)%units = '%'
                if(irh==2) then 
                  dyn(2)%q(:,:,:,dyn(1)%grid%lm+2) = rh2
                  dyn(2)%qm(dyn(1)%grid%lm+2)%name = 'mrh';  dyn(2)%qm(dyn(1)%grid%lm+2)%long_name = 'Mean Relative Humidity '
                  dyn(2)%qm(dyn(1)%grid%lm+2)%units = '%'
                endif
                call dyn_put ( trim(dyn_dout), nymd_d, nhms_d, 0, dyn(2), rc, freq=freq, vectype=dyntype, skip_setvec=.true. )
                deallocate(rh1,rh2)
                call dyn_clean ( dyn(2) )
             else
                dyn(1)%grid%lm = lm
                call dyn_put ( trim(dyn_dout), nymd_d, nhms_d, 0, dyn(1), rc, freq=freq, vectype=dyntype )
             endif
        endif 

!       Clean up mess
!       -------------
        call dyn_clean ( dyn(1) )

      end do

   end do ! loop over files

   st = "rm -f " // trim(binfiles(1)) // " " // trim(binfiles(2))  
   rc = system(st)      ! rc == 0 for success
   if (rc .ne. 0 ) then
     print *, "Unable to remove binary files."
   end if

!  All done
!  --------
!  All done
!  --------
   close(999)
   open (999,file=trim(egress),form='formatted')
   close(999)
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
      subroutine Init_ ( dyntype, mfiles, files, dominmax, verb, egress, &
                         eps_eer, anorm, jnorm, tv2t, &
                         projlon, projlat, projlev, normlz, ntype, addrh, &
                         vnorm )

      use m_inpak90
      use m_chars,   only: lowercase

      implicit NONE

      integer,       intent(out) :: dyntype ! 4=geos4, 5=geos5
      integer,       intent(in)  :: mfiles  ! max. number of eta files
                                            ! dynamics file names (eta)
      character(len=*), intent(out) :: files(mfiles) 
      character(len=*), intent(out) :: egress
      character(len=*), intent(out) :: anorm
      character(len=*), intent(out) :: jnorm
      integer, intent(out) :: vnorm
      real,    intent(out) :: eps_eer
      real,    intent(out) :: projlat(2), projlon(2)
      integer, intent(out) :: projlev(2)
      logical, intent(out) :: dominmax
      logical, intent(out) :: tv2t
      logical, intent(out) :: verb
      logical, intent(out) :: normlz
      integer, intent(out) :: addrh
      character(len=*), intent(out) :: ntype
      
!
! !REVISION HISTORY:
!       2002.01.07  E. Yeh   Initial code.
!       05Oct2012   Todling  Add verb
!       09Jun2020   Todling  Add rh option
!
!EOP
!BOC

      character*4, parameter :: myname_ = trim(myname)//'*init'

      character(len=256) :: rcfile
      character(len=255) :: etafile, argv
      integer iret, i, iarg, argc, iargc
      real    pnext
      logical dout
      logical invalid

      dout = .false.
      verb = .false.

      anorm    = 'twe'
      jnorm    = 'NULL'
      dyn_dout = 'NONE'
      dyntype  = 4        ! default is GEOS-4 files
      dominmax = .false.
      acoeff = 1.0
      egress = 'DYNDIFF_EGRESS'
      rcfile = 'NULL'
      addrh  =  0
      sbyene = .false.
      tv2t   = .false.
      eps_eer = 1.0
      normlz = .false.
      ntype  = 'ene'
      vnorm  = 0

!     Default LPO boundaries
!     ----------------------
      projlon(1) =    0.
      projlon(2) =  360.
      projlat(1) =  -90.
      projlat(2) =   90.
      projlev(1) =    1
      projlev(2) =   99

      print *
      print *, '     ---------------------------------------------------------'
      print *, '     dyndiff - dynamic vector difference between two hdf files'
      print *, '     ---------------------------------------------------------'
      print *

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()

      iarg = 0
      nfiles = 0

      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iarg, argv )
         
         select case (argv)
           case ("-g5")
             dyntype = 5
           case ("-verb")
             verb = .true.
           case ("-a")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, argv )
             read(argv,*) acoeff
           case ("-egress")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, egress )
           case ("-o")
             dout = .true.
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, dyn_dout )
           case ("-txt")
             dout = .true.
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, jnorm )
           case ("-addrh")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, argv )
             read(argv,*) addrh
           case ("-tv2t")
             tv2t = .true.
           case ("-h")
             if ( iarg+1 .gt. argc ) call usage()
           case ("-ene_scale")
             sbyene=.true.
           case ("-normlz")
             normlz=.true.
           case ("-ntype")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, ntype )
           case ("-rc")
             if ( iarg+1 .gt. argc ) call usage()
             iarg = iarg + 1
             call GetArg ( iArg, rcfile )
           case default
             nfiles = nfiles + 1
             if ( nfiles .gt. mfiles ) call usage()
             files(nfiles) = argv

         end select

      end do

      if ( nfiles .lt. 1 ) call usage()
      if ( nfiles .eq. 1 ) then
           nfiles = 2
           files(2) = files(1)
           dominmax = .true.
      endif
      if ( nfiles .gt. 2 ) call usage()

!     If rc-file specificied and scale-by-energy, read relevant options
!     -----------------------------------------------------------------
      if (trim(rcfile)/='NULL' .and. sbyene ) then

!        Load RC file
!        ------------
         call i90_loadf (trim(rcfile), iret)
         if( iret .ne. 0) then
             write(6,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
             call exit (1)
         endif

!        Read in norm type
!        -----------------
         call I90_label('pert_norm:', iret)
         if (iret .eq. 0) then
            call I90_Gtoken ( argv, iret )
            if ( iret==0 ) then
                anorm = lowercase(trim(argv))
                write(6,'(2a)') 'Norm for perturbation: ', trim(anorm)
             endif
         endif

!        Read Ehrendorfer, Errico, and Raeder's epsilon factor (apply to Q-component)
!        ----------------------------------------------------------------------------
         call I90_label('ehrendorfer_errico_raedder_eps:', iret)
         if (iret .ne. 0) then
           write(6,'(2a,i5)') myname_, ': I90_label error, iret =',iret
         else
           eps_eer = I90_GFloat(iret)
           if( iret .ne. 0) then
              write(6,'(3a,i5)') myname_,': I90_GFloat error, ', ' iret =',iret
              call exit (1)
           end if
         end if
         write(6,'(a,1p,e13.6)') 'Ehrendorfer, Errico, and Raeder eps: ',eps_eer

!        Read in LPO information
!        -----------------------
         call I90_label('local_svec_projection:', iret)
         if (iret .ne. 0) then
           if(verb) write(6,'(2a)') myname_, &
                        ': Using default local projection operator, i.e., none'
         else
             invalid = .false.
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlon(1) = pnext
             end if
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlon(2) = pnext
             end if
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlat(1) = pnext
             end if
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlat(2) = pnext
             end if
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlev(1) = pnext
             end if
             pnext = I90_GFloat(iret)
             if (iret/=0) invalid = .true.
             if (.not.invalid) then
                projlev(2) = pnext
             end if

!            Quick sanity check
!            ------------------
             if ( projlat(1) > projlat(2) ) invalid = .true.
             if ( projlev(1) > projlev(2) ) invalid = .true.
   
             if ( invalid ) then
                  projlon(1) =    0.
                  projlon(2) =  360.
                  projlat(1) =  -90.
                  projlat(2) =   90.
                  projlev(1) =    1
                  projlev(2) =   99
                  write(6,'(2a,/,a)') myname_, &
                     ': Something went wrong while setting projection box.', &
                     '  Taking default local projection (-180,180) (-90,90) (All Levels)'
             else
                 if(verb) write(6,'(a,/,2(a,f7.2,a,f7.2,/),2(a,i7),/)')    &
                    'User specified local projection: ',               &
                    '  From Lon ', projlon(1), ' to Lon ', projlon(2), &
                    '  From Lat ', projlat(1), ' to Lat ', projlat(2), &
                    '  From Lev ', projlev(1), ' to Lev ', projlev(2)
             endif
         endif
         print *
         print *, "Norm type: ", trim(ntype)

!        Decide between E- and V-norms
!        -----------------------------
         call I90_label('vnorm:', iret)
         if (iret .ne. 0) then
           write(6,'(2a)') myname_,    &
                        ': I90_label not found will use default '
         else
           vnorm = I90_GInt(iret)
           if (iret .ne. 0) then
               write(6,'(2a,i5)') myname_,    &
                            ': I90_Gtoken error, iret =',iret
               call exit(2)
           end if
         end if

         tv2t = .true.
      endif
      if (tv2t .and. verb) print *, 'Will convert virtual T to T'

!     Echo the parameters
!     -------------------
      print *
      print *, '------------------------------------------------------------------'
      print *, '  Eta     Dynamics state files: '
      do i = 1, nfiles
        print *, i, ': ', trim(files(i))
      end do
      if(dout) then
         select case (abs(addrh))
         case (1)
           print *, 'Adding GEOS-qsat-based RH to ouput file'
         case (2)
           print *, 'Adding BeCov-qsat-based RH to ouput file'
         case (3)
           print *, 'Adding GSI-qsat-based RH to ouput file'
         case default
           ! nothing to say
         end select
      endif

      end subroutine Init_

!.................................................................

      subroutine usage()
      print *
      print *,'Usage: '
      print *
      print *,'  dyndiff.x [-h] [-g5] [-verb] etafile_1 etafile_2 [-o diff_file]'
      print *
      print *, 'where'
      print *
      print *, '-h          Help (optional)'
      print *, '-verb       Echo different to standard out'
      print *, '              (default: FALSE) '
      print *, '-g5         Treats files as GEOS-5 files'
      print *, '-tv2t       Converts diff in Tv to diff in T'
      print *, '-egress     Name of EGRESS file for successful finalization'
      print *, '-a  coeff   Scale difference by this coefficient (see note)'
      print *, '-addrh N    Add rh diff to file: 1=use GEOS  qsat'
      print *, '                                 2=use BeCov qsat'
      print *, '                                 3=use GSI   qsat'
      print *
      print *
      print *, 'Where etafile_1 and etafile_2 are two (required)'
      print *, 'input dynamics vector files in hybrid (eta) coordinates'
      print *
      print *, '-o diff_file      Optional to save binary output'
      print *, '                  where diff_file - output dynamics vector'
      print *, '                  file in hybrid (eta) coordinates'
      print *, '                  with timestamp from the etafile_1'
      print *, '                  and 32-bit precision'
      print *, '                  output = file1 - file2'
      print *, '                  (default: difference is shown at the screen in ascii)'
      print *
      print *, 'Notes: '
      print *, '  1. User can save the ascii output by the following command:'
      print *, '     dyndiff.x etafile_1 etafile_2 > OutPutFileName' 
      print *, '  2. Ability to scale dff by a coefficient must be exercized'
      print *, '     with caution, since all entries will be scaled - which is '
      print *, '     in general meaningless. This ability is added to compensate'
      print *, '     for the fact that sometimes file1-file2 not possible, but'
      print *, '     file2-file1 is possible due to nc4-header issues'
      print *, '  3. If addrh<0, mean rh is added to the file (serves BeCov code)'
      call exit(1)
      end subroutine usage
      
!.................................................................

      subroutine die ( myname, msg )
      character(len=*) :: myname, msg
      write(*,'(a)') trim(myname) // ': ' // trim(msg)
      call exit(1)
      end subroutine die

      subroutine getrh_(rh,tv,qv,ps,ak,bk)
      use m_const, only: zvir
      use GEOS_UtilsMod, only: GEOS_Qsat
      use m_dyn_util, only: dyn_qsat
      implicit none
      real,intent(out):: rh(:,:,:)
      real,intent(in) :: tv(:,:,:), qv(:,:,:), ps(:,:)
      real,intent(in) :: ak(:), bk(:)
      real(4),allocatable :: tmp(:,:,:),pmk(:,:,:),qs(:,:,:)
      integer i,j,k,kb
      select case (abs(addrh))
      case (1)
         allocate(tmp(im,jm,1),pmk(im,jm,1),qs(im,jm,1))
         do k = 1, km
            pmk(:,:,1) = 0.5 * (  ak(k)+ak(k+1)  + &
                              ps*(bk(k)+bk(k+1)) )
            tmp(:,:,1) = tv(:,:,k)/(1.0+zvir*qv(:,:,k))
            qs (:,:,1) = GEOS_Qsat(tmp(:,:,1), pmk(:,:,1), PASCALS=.true.)
            rh(:,:,k) = qv(:,:,k)/qs(:,:,1)
         end do
      case (2)
         allocate(tmp(im,jm,1),pmk(im,jm,1),qs(im,jm,1))
         if (fix_top_rh) then
            k=1
            pmk(:,:,1) = 0.5 * (  ak(k)+ak(k+1)  + &
                              ps*(bk(k)+bk(k+1)) )
            tmp(:,:,1) = tv(:,:,k)/(1.0+zvir*qv(:,:,k))
            qs (:,:,1) = GEOS_Qsat(tmp(:,:,1), pmk(:,:,1), PASCALS=.true.)
            rh(:,:,k) = qv(:,:,k)/qs(:,:,1)
            kb = 2
         else
            kb = 1
         endif
         do k = kb, km
            pmk(:,:,1) = 0.5 * (  ak(k)+ak(k+1)  + &
                              ps*(bk(k)+bk(k+1)) )
            tmp(:,:,1) = tv(:,:,k)
            qs (:,:,1) = qv(:,:,k)
            call dyn_qsat(qs(:,:,1),tmp(:,:,1),pmk(:,:,1),im,jm,.true.)
            rh(:,:,k) = qv(:,:,k)/qs(:,:,1)
         end do
      case (3)
         allocate(tmp(im,jm,km),pmk(im,jm,km),qs(im,jm,km))
         do k = 1, km
             pmk(:,:,k) = 0.5 * (  ak(k)+ak(k+1)  + &
                               ps*(bk(k)+bk(k+1)) )
             tmp(:,:,k) = tv(:,:,k)/(1.0+zvir*qv(:,:,k))
         enddo
         call dyn_qsat(qs,tmp,pmk,im,jm,km,.true.,ntop=1)
         rh = qv/qs
         if (fix_top_rh) then
            k=1
            pmk(:,:,1) = 0.5 * (  ak(k)+ak(k+1)  + &
                              ps*(bk(k)+bk(k+1)) )
            tmp(:,:,1) = tv(:,:,k)/(1.0+zvir*qv(:,:,k))
            qs (:,:,1) = GEOS_Qsat(tmp(:,:,1), pmk(:,:,1), PASCALS=.true.)
            rh(:,:,1)=qv(:,:,1)/qs(:,:,1)
         endif
      case default
         return
      end select
      deallocate(tmp,pmk,qs)
      end subroutine getrh_
!.................................................................

  end program dyndiff
