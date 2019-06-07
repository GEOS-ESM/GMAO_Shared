!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ss2fv: convert spectral fields to eta
!
! !INTERFACE:
!
   program ss2fv

!
! !USAGE: ss2fv.x ncepfile etafile [options]
!
! !USES:
!
   use m_ss
   use m_dyn
   use m_mapz, only : set_eta
   use m_xform, only : xform
   use m_ioutil, only : luavail
   use make_surf, only : ncep2fv_ts
   use m_die  , only : die
   use util
   implicit NONE

! !DESCRIPTION: Convert NCEP's SSI analysis spectral fields to FvDAS
!               eta format.
!
! !REMARKS:
! 1) Implementaion of Tskin analysis assumes that same Tskin(bkg.sfc)=Ts(bkg.eta)
!    This analysis updates Tskin in ana.eta to be GSI-analyzed Tskin acconting for
!    the no-data/no-harm case.
!
! !REVISION HISTORY:
!
!  2002.10.30  C. Cruz:   Initial Code
!  04Aug2003   Todling    Removed attribute "new" for out eta/GFIO file
!  08Aug2003   Todling    Added logic to keep lwi, ts, and hs_stdv from bkg
!  06Nov2003   Todling    - Restated references refs to "restart" to eta-file
!                         - Added nstep in/out
!                         - Bug fix in write out of GCM rst file(q was bad)
!  26Mar2004   Todling    Added oldyn to avoid reading nstep.
!  26Sep2004   Todling    Rearranged dyn-vect initializations and clean-ups.
!  14Dec2004   Todling    Added passing of GSI Tskin analysis to ana.eta
!  12Jan2005   Todling    Adjusted set_eta interface to be consistent w/ m_mapz
!                         replaced ks=km in dyn_init to ks=ks
!  15Feb2005   Dee        Handle spectral-eta file produced by fv2ss
!  25Feb2005   Todling    Bug fix in combo eta opt and lack of ana.lwi/hs_stdv
!  23Jan2005   B Zhang    Updated for total cloud water fraction
!
!-------------------------------------------------------------------------
!EOP

   character(len=*), parameter :: myname = 'ss2fv'

! fv dynamics vectors

   type(dyn_vect) :: z_f    ! transformed eta first guess with error
   type(dyn_vect) :: w_ff   ! first guess
   type(dyn_vect) :: w_fs   ! w_ff remapped 
   type(dyn_vect) :: w_as   ! output dynamics vector

! spectral space vector

   type(ss_vect)  :: ss     ! analysis vector in spectral space

! various flags

   logical        :: lua    ! flag for blending NCEP's with
                            ! fvDAS' stratosphere
   logical        :: llwi   ! flag to use land-water-ice mask
   logical        :: lts    ! flag to use surface T
   logical        :: wrs    ! flag to  write a restart
   logical        :: fex    ! file exists flag 
   logical        :: pick   ! whether to read a specified
                            !  time from first guess file
   logical        :: oldyn     ! control eta file age (e.g., older ana files don't have nstep)
   logical        :: notracer  ! control tracer write out
   logical        :: bkgfrq ! when .t. uses input file freq to output file
   logical        :: docw =.false.  ! controls whether or not calculate cloud-water cond
   
   integer, parameter :: FV_PREC = 1  ! dyn_put default precision (64 bits)

!  Grid

   logical        :: eta  ! true if input is spectral-eta
   real, pointer  :: ak(:), bk(:) ! vertical grid coefficients
   real           :: ptop, pint, p1, p2
   integer        :: ks

!  Other locals

   integer :: lu, n,iargc,rc
   integer :: nymd, nhms
   integer :: nymd0, nhms0
   integer :: nymd1, nhms1
   integer :: im,jm,km,lm
   integer :: freq, ofreq, nstep
   integer :: prec

   real, allocatable :: cwmr(:,:,:)
!  File names

   character(len=255) :: outfile, ssifile, infile, insfcf, &
                         fv_rst, lwifile, tsfile, uafile, ndnhfn

! start

   nstep = 18000   ! default so that all files from this point on have it specified
   ofreq = 060000  ! default freq of output file
   prec  = FV_PREC ! default precision for output file

   call init ( outfile, ssifile, infile, insfcf, ndnhfn, pick, nymd, nhms, oldyn, notracer, &
               prec, p1, p2, lua, llwi, lts, lwifile, tsfile, uafile, wrs, bkgfrq, rc )

   call myout ( 6, myname, 'Initialize ss vector')
   call ss_init ( ssifile, ss, rc )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot initialize ss vector', rc )
   end if
   
   call myout ( 6, myname, 'Get ss vector')
   call ss_get ( ssifile, ss, rc )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot get ss vector', rc )
   end if

   ! read eta infile file into w_ff
   call myout ( 6, myname, 'read w_f' )
   if ( pick ) then
      if (oldyn) then
          call dyn_get ( infile, nymd, nhms, w_ff, rc, timidx=0, freq=freq )
       else
          call dyn_get ( infile, nymd, nhms, w_ff, rc, timidx=0, freq=freq, nstep=nstep )
       endif
   else
      if (oldyn) then
          call dyn_get ( infile, nymd, nhms, w_ff, rc, freq=freq )
      else
          call dyn_get ( infile, nymd, nhms, w_ff, rc, freq=freq, nstep=nstep )
      endif
   end if
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot read dynamics vector wff', rc )
   else
     print *, myname//': read dynamics state file '//trim(infile)
     print *, myname//': nymd, nhms = ', nymd, nhms
     print *
   end if

! inquire dimensions and initialize meta data for output w_as
! -----------------------------------------------------------
   call Dyn_GetDim ( w_ff, im, jm, km, lm )
   print *, myname // ': resolution : (',im,'x',jm,'x',km,'x',lm,')'

   call myout ( 6, myname, 'Compute grid')
   call myalloc(ak,km+1)
   call myalloc(bk,km+1)
   call set_eta ( km, ks, ptop, pint, ak, bk )

   call myout ( 6, myname, 'Initialize fv vector' )
   call dyn_init ( im, jm, km, lm, w_as, rc , &
                   ptop = ptop, &
                   ks = ks, &
                   ak = ak, &
                   bk = bk )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot initialize fv dynamics vector', rc )
   end if
   
   if ( w_ff%grid%lm >= 3 ) docw = .true.
   if ( docw ) then
        print *, myname//': Extracting cloud-water condensate mixing ratio from background'
        allocate ( cwmr(w_ff%grid%im,w_ff%grid%jm,w_ff%grid%km) )
   endif

! check vertical coordinate definition of input spectral file
! -----------------------------------------------------------   
   eta = ss%meta%eta
   if (eta) then
      call myout ( 6, myname, 'ss vector has eta coordinates')
      if (size(ss%meta%ak) .ne. km+1 .or. size(ss%meta%bk) .ne. km+1) then
         call myexit ( myname,'eta coordinates of ss and fv vectors are different',99 ) 
      end if
   end if

! no data no harm

   ! check if first ss2fv of FVSSI cycle was performed
   ! this would imply a file "ndnhfn" was created
   inquire ( FILE = trim(ndnhfn), EXIST = fex )

   if ( fex ) then

     call myout ( 6, myname, 'will ensure no data no harm' )

     ! read z_f created in first ss2fv step
     call myout ( 6, myname, 'read z_f' )
     call dyn_get ( trim(ndnhfn), nymd0, nhms0, z_f, rc ) 
     if ( rc .ne. 0 ) then
       call myexit ( myname,'cannot read dynamics vector z_f', rc )
     else
       print *, myname//': read dynamics state file '// trim(ndnhfn)
       print *, myname//': nymd, nhms = ', nymd0, nhms0
       print *
     end if

   end if

   call myout ( 6, myname, 'Transform ss vector to fv vector...please wait' )
   if ( lua .or. eta ) then 

     if ( fex ) then  ! no data no harm case

       call myout ( 6, myname, 'initialize w_fs' )
     ! initialize vector that will hold remapped w_f = w_fs
       call dyn_init ( im, jm, km, lm, w_fs, rc , &
                     ptop = ptop, &
                     ks = ks, &
                     ak = ak, &
                     bk = bk )
       if ( rc .ne. 0 ) then
         call myexit ( myname,'cannot initialize fv dynamics vector', rc )
       end if

       if ( eta ) then
          call xform ( nymd, nhms, pick, &
                    ss, w_as, ak, bk, w_ff, rc, &
                    w_fs = w_fs, z_f = z_f, &
                    eta=eta, cwmr=cwmr )
       else   ! blending upper-air fields 
          call xform ( nymd, nhms, pick, &
                    ss, w_as, ak, bk, w_ff, rc, &
                    w_fs = w_fs, z_f = z_f, &
                    uafile = uafile, pabove = p1, pbelow = p2, eta=eta, cwmr=cwmr )
       end if
         if ( rc .ne. 0 ) then
           call myexit ( myname,'could not transform ss vector to fv vector', rc )
         end if

!      Handle Tskin analysis in no-data no-harm mode (see remarks)
!      ---------------------------------------------
       w_as%ts = w_ff%ts
       if ( trim(insfcf)/='NONE') then
            call ncep2fv_ts ( w_as%grid%im, w_as%grid%jm, ss%meta%lonf, ss%meta%latf-2, &
                              w_as%ts, insfcf, rc, Tskin_Ge=z_f%ts )
              if(rc/=0) call myexit ( myname,'cannot perform Tskin analysis', rc )
       end if
 
       call dyn_clean ( w_fs )

     else             ! step 2 w/ blending

       print *, 'Doing blend-only case ...'; call flush(6)
       call xform ( nymd, nhms, pick, &
                    ss, w_as, ak, bk, w_ff, rc, &
                    uafile = uafile, pabove = p1, pbelow = p2, eta=eta, cwmr=cwmr )
       if ( rc .ne. 0 ) then
         call myexit ( myname,'could not transform ss vector to fv vector', rc )
       end if

!      Handle Tskin analysis in no-data no-harm mode (see remarks)
!      ---------------------------------------------
       w_as%ts = w_ff%ts
       if ( trim(insfcf)/='NONE') then
            call ncep2fv_ts ( w_as%grid%im, w_as%grid%jm, ss%meta%lonf, ss%meta%latf-2, &
                              w_as%ts, insfcf, rc )
              if(rc/=0) call myexit ( myname,'cannot perform Tskin analysis', rc )
       endif

     end if
   else               ! step 2 w/o blending
     print *, 'Doing no-blend case'
     if ( .not. eta ) print *, '(NO no-data-no-harm)'
     call flush(6)

     call xform ( nymd, nhms, pick, ss, w_as, ak, bk, w_ff, rc, eta=eta, cwmr=cwmr )
     if ( rc .ne. 0 ) then
       call die(myname,'could not transform ss vector to fv vector')
     end if
   end if

   if ( lts  ) then
     call myout ( 6, myname, 'Read NCEPs Ts' )
     call get_ncep_tsx( tsfile, w_as%ts, w_as%grid%im, w_as%grid%jm, nymd, nhms )
   end if
   w_as%lwi     = w_ff%lwi
   w_as%hs_stdv = w_ff%hs_stdv
   if ( llwi ) then
     call myout ( 6, myname, 'Read LWI' )
     lu = luavail()
     open ( lu,file=trim(lwifile),form='unformatted',access='sequential')
     read ( lu ) w_as%lwi
     close( lu )
   endif

! output

   if(notracer) w_as%grid%lm = 1  ! force write out of spec hum only

   call myout ( 6, myname, 'Write out eta (hdf)' )
   if ( bkgfrq ) ofreq = freq
   write(*,*) '  nymd = ',nymd  
   write(*,*) '  nhms = ',nhms 
   write(*,*) '  freq = ',ofreq 
   write(*,*) '  nstep= ',nstep 
   call dyn_put ( outfile, nymd, nhms, prec, w_as, rc,  &
                  freq=ofreq, nstep=nstep, verbose = .false.)
   if ( rc .ne. 0 ) then
      call myexit ( myname,'could not write out eta', rc )
   end if
   ! write a restart if so desired
   if ( wrs ) then
      call myout ( 6, myname, 'Write out restart (bin)' )
      call put_fvrst ( outfile, w_as, nymd, nhms, 6 )
   end if

! clean up

   call myout ( 6, myname, 'Clean up' )
   deallocate     ( ak, bk )
   if(docw) deallocate ( cwmr )
   call ss_clean  ( ss )
   call dyn_clean ( w_as )
   call dyn_clean ( w_ff )
   if ( fex ) then
     call dyn_clean ( z_f  )
   end if

! done
   
   call myout ( 6, myname, '-- ss2fv.x has successfully ended --' )
   call exit(0)

   CONTAINS

!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !IROUTINE: init --- Initialize ss2fv
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
   subroutine init ( outfile, ssifile, infile, insfcf, ndnhfn, pick, nymd, nhms, & 
                     oldyn, notracer,  prec, &
                     p1, p2, lua, llwi, lts, lwifile, tsfile, uafile, wrs, & 
                     bkgfrq, rc ) 

   implicit NONE

   integer, intent(out) :: rc
   character*255, intent(out) :: outfile  ! eta-hdf output filename
   character*255, intent(out) :: ndnhfn   ! dyn-vec filename for no-data no-harm
   character*255, intent(out) :: ssifile
   character*255, intent(out) :: infile
   character*255, intent(out) :: insfcf   ! sfc-ncep-like (gaussian) filename
   integer, intent(out) :: nymd ! date for first guess
   integer, intent(out) :: nhms ! time "   "     "
   integer, intent(inout) :: prec ! output precision (32/64 bits)
   real, intent(out) :: p1
   real, intent(out) :: p2
   character*255, intent(out) :: lwifile
   character*255, intent(out) :: tsfile
   character*255, intent(out) :: uafile
   logical, intent(out) :: llwi ! flag to use land-water-ice mask
   logical, intent(out) :: lts  ! flag to use surface T
   logical, intent(out) :: lua  ! flag for blending NCEP's with
                                ! fvDAS' ptop
   logical, intent(out) :: wrs  ! write restart flag
   logical, intent(out) :: pick ! if true,
                                ! (nymd,nhms) are input
                                ! parameters to dyn_get. 
   logical, intent(out) :: oldyn   ! controls older eta-type files attributes
   logical, intent(out) :: notracer ! controls whether output will have tracers or not
   logical, intent(out) :: bkgfrq ! when .t. uses input file freq to output file

!
! !REVISION HISTORY:
!
!  2002.10.30  C. Cruz  Initial code.
!  2003.11.21  Todling  Removed im,jm,km,lm initialization from here
!  19Mar2004   D Kokron Added test to generate a warning if p1 is greater than p2
!  26Mar2004   Todling  Added parameter oldyn/notracer
!  13Dec2004   Todling  Added insfcf to handle surface analysis
!  07Nov2005   Todling  Turned name of no-data no-harm filename into opt argument
!
!-----------------------------------------------------------------------
!

   character*4, parameter :: myname = 'init'
   integer :: n,iargc,nargs,iaux
   real    :: a2
   character*255 :: argv
   character*255, allocatable :: arg(:)

   nargs =  iargc()
   if( nargs.lt.2 ) call usage()
   allocate ( arg(nargs) )

! defaults

   insfcf  = 'NONE'
   outfile = 'ss2fv_eta_out.nc4'
   ndnhfn  = 'eta_out.nc4'
   llwi = .false.
   lts  = .false.
   lua  = .false.
   wrs  = .false.
   pick = .false.
   oldyn    = .false.
   notracer = .false.
   bkgfrq = .false.
   p1 = 10.0
   p2 = 30.0

! process options

   do n=1,nargs
     call getarg(n,arg(n))
   enddo
   do n=1,nargs
     if( trim(arg(n)).eq.'-pa'  ) read(arg(n+1), * ) p1
     if( trim(arg(n)).eq.'-pb'  ) read(arg(n+1), * ) p2
     if( trim(arg(n)).eq.'-o'   ) outfile = trim(arg(n+1))
     if( trim(arg(n)).eq.'-i'   ) ndnhfn = trim(arg(n+1))
     if( trim(arg(n)).eq.'-is'  ) insfcf = trim(arg(n+1))
     if( trim(arg(n)).eq.'-res' ) print *, myname//': -res is obsolete'
     if( trim(arg(n)).eq.'-oldyn'    ) oldyn    = .true.
     if( trim(arg(n)).eq.'-notracer' ) notracer = .true.
     if( trim(arg(n)).eq.'-bkgfrq' ) bkgfrq = .true.
     if( trim(arg(n)).eq.'-ua'  ) then
       uafile  = trim(arg(n+1))
       lua     = .true.
     endif
     if( trim(arg(n)).eq.'-rs'  ) wrs = .true.
     if( trim(arg(n)).eq.'-ts'  ) then
       tsfile  = trim(arg(n+1))
       lts     = .true.
     endif
     if( trim(arg(n)).eq.'-lwi' ) then
       lwifile  = trim(arg(n+1))
       llwi     = .true.
     endif
     if( trim(arg(n)).eq.'-pick' ) then
        read(arg(n+1),*) nymd
        read(arg(n+2),*) nhms
        pick = .true.
     end if
     if( trim(arg(n)).eq.'-prec' ) then
        read(arg(n+1),*) iaux
        if(iaux==32) prec = 0
        if(iaux==64) prec = 1
     endif
   enddo
   if (p1 .gt. p2) then
     print *, myname//': WARNING pabove must be a smaller number than pbelow'
     print *, myname//': in order for blending to work properly'
     print *, myname//': reversing input blending levels'
     a2 = p1
     p1 = p2
     p2 = a2
   endif

   ssifile=trim(arg(1))
   infile=trim(arg(2))
   !if (.not. lua) uafile  = infile

!  Echo the parameters
!  -------------------
   print *
   print *, '    ss2fv.x : convert spectral fields to eta'
   print *, '--------------------------------------------------'
   print *, ' analysis filename   : ',trim(ssifile)
   print *, ' eta_out  filename   : ',trim(outfile)
   print *, ' eta_in   filename   : ',trim(infile)
   print *, ' sfc_in   filename   : ',trim(insfcf)
   print *, ' no-data no-harm file: ',trim(ndnhfn)
   if (pick) then
   print *, ' NHMS                : ',nhms
   print *, ' NYMD                : ',nymd
   end if
   print *

   deallocate(arg)
   rc = 0

   end subroutine init

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
   subroutine usage
!-------------------------------------------------------------------------

   print *
   print *,'Usage:  '
   print *
   print *, "ss2fv.x ncep_fname eta_name"
   print *
   print *, "options :"
   print *
   print *, " [-o   fv_eta_fname] (default: ss2fv_eta_out.nc4)"
   print *, " [-i   fname       ] (default: eta_out.nc4)"
   print *, " [-time            ] (default: 1)"
   print *, " [-lwi fname       ] (default: none)"
   print *, " [-ts  fname       ] (default: none)"
   print *, " [-ua  ua_eta_fname] (default: none)"
   print *, " [-pa  pabove      ] (default: 10 mb for T62L28)"
   print *, " [-pb  pbelow      ] (default: 30 mb for T62L28)"
   print *, " [-rs              ] (default: none)"
   print *, " [-oldyn           ] (default: new dyn)"
   print *, " [-notracer        ] (default: all vars)"
   print *
   print *, "where:"
   print *
   print *, " ncep_fname      : Filename of NCEP sigma-level analysis in spectral space"
   print *, " eta_fname       : Filename of FvDAS dynamics (eta file)"
   print *, " -o fv_eta_fname : Filename of FvDAS dynamics in eta format"
   print *, " -i fname        : Filename of dyn-vect used for no-data no-harm procedure"
   print *, " -os fv_sfc_fname: Output filename with FvDAS surface fields"
   print *, " -is ncep_sfcfn  : Input  filename with FvDAS surface fields"
   print *, " -time time idx  : time index: 1=0z, 2=6z, 3=12z, 4=18z"
   print *, " -lwi fname      : Filename of fv phys restart" 
   print *, " -ts fname       : Filename of ncep surface datafile" 
   print *, " -ua ts_fname    : Filename of FV to get upper air fields"
   print *, " -pa pabove      : Pressure above which blending is not performed"
   print *, " -pb pbelow      : Pressure below which blending is not performed"
   print *, " -rs             : Write restart (bin) file"
   print *, " -oldyn          : Specify when using eta files w/o nstep"
   print *, " -notracer       : Specify when no tracers are needed on output"
   print *, " -prec     PREC  : (32/64) Precision of output; default 64 bits"
   print *
   print *, "NOTE: The file ncep_fname can be a spectral-eta file produced by"
   print *, "      fv2ss.x. In that case the -ua, -pa, -pb options are ignored" 
   call exit(1)
   end subroutine usage
 
!-------------------------------------------------------------------------

   subroutine put_fvrst ( etafile, w_f, nymd,nhms,nstep )
   use m_ioutil, only : luavail
   implicit none
   character(len=*), parameter :: myname = 'put_fvrst'
   character*255 :: binfile
   character*255 :: etafile
   type(dyn_vect) :: w_f 
   integer   :: im,jm,km,lm,nymd,nhms,nstep,leta
   integer   :: i ,j ,k ,l
   integer   :: lu
   
   im = w_f%grid%im
   jm = w_f%grid%jm
   km = w_f%grid%km
   lm = w_f%grid%lm
   lu=luavail()
   leta = len(trim(etafile))-3
   binfile = etafile(1:leta)//'bin'
   open ( unit=lu, file = binfile, form = 'unformatted' )
   write (lu) nstep,nymd,nhms
   write (lu) w_f%ps,w_f%delp,w_f%u,w_f%v,w_f%pt
   do l = 1, lm
      write(lu) (((w_f%q(i,j,k,l),i=1,im),j=1,jm),k=1,km)
   end do
   close (lu)

   end subroutine

   end program ss2fv
