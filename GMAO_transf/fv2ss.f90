!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP        
!
! !ROUTINE: fv2ss: convert eta fields to spectral
!
! !INTERFACE:
!
   program fv2ss

!
! !USAGE: fv2ss etafile [options]
!
! !USES:
!
   use m_ss
   use m_dyn
   use m_xform, only : xform
   use util
   use m_const, only : undef
   use m_set_eta
   implicit NONE
!
! !DESCRIPTION: This program performs the conversion of an eta restart
!               to a spectral restart used by the SSI.
! 
! !REVISION HISTORY:
!
!  2002.10.30  C. Cruz:   Initial Code
!  15Mar2004   Todling    Added cloud-water fraction capability
!  19Oct2004   WGu/RT     Bug fix: added NCEP sfc filename
!  15Feb2005   Dee        Added option to inherit fv eta vertical coordinate
!                         and produce a spectral-eta restart file.
!                         To activate, invoke fv2ss without the -nsig flag 
!  19Jan2006   B. Zhang   Updated cloud-water based on changed FV variable
!  2009.06.22  Todling    Add GEOS-5 (new vec) support
!  2009.10.07  Sienkiewicz allow hybrid vertical coordinate on non-fv 
!                             levels with remapping (via -hybrid flag)
!                             nsig is still the number of levels
!
!-------------------------------------------------------------------------
!EOP

   character(len=*), parameter :: myname = 'fv2ss'
   type(dyn_vect) :: w_f     ! fv dynamics vector in eta
   type(ss_vect)  :: ss      ! analysis vector in spectral space
   integer :: jcap
   integer :: nsig
   integer :: fhour

!  Locals
   
   integer :: n, rc, ks
   integer :: nymd, nhms
   integer :: im,jm,km,lm
   integer :: vectype
   logical :: pick           ! whether to read a specified
                             !  time from first guess file
   logical :: ncep_phis      ! use ncep_phis when .t.
   logical :: use_sigio      ! write output with sigio_module if .t.
   logical :: eta = .true.   ! inherit fv eta coordinate
   logical :: docw =.false.  ! controls whether or not pass cloud-water cond to GSI
   logical :: hybrid         ! controls whether output is sigma or hybrid sig-p

   real, allocatable :: phis(:,:)
   real, allocatable :: cwmr(:,:,:)
   real, allocatable :: ako(:), bko(:)
   real  ptop, pint

!  File names

   character(len=255) :: etafile, sfcfvfn, ssiupaf, ssisfcf
   character(len=255) evalue

! start

   call init( etafile, sfcfvfn, ssiupaf, ssisfcf, jcap, nsig, fhour, & 
              pick, nymd, nhms, ncep_phis, use_sigio, vectype, hybrid, rc )
   if ( nsig>0 ) eta = .false.
   if ( hybrid .and. w_f%grid%km .eq. nsig ) then
      print *,'Using input vertical coord without remapping (eta=.true.)'
      eta = .true.
   end if

   if ( eta ) ncep_phis = .false.  ! just making sure

   call myout ( 6, myname, 'Get eta fields' )
   if ( pick ) then
      call dyn_get ( etafile, nymd, nhms, w_f, rc, timidx=0, vectype=vectype  )
   else
      call dyn_get ( etafile, nymd, nhms, w_f, rc, vectype=vectype  )
   end if
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot read dynamics vector file', rc )
   else
     print *, myname//': read dynamics state file '//trim(etafile)
     print *, myname//': nymd, nhms = ', nymd, nhms
     print *
   end if
   
   if ( w_f%grid%lm >= 3 ) docw = .true.
   if ( docw ) then
        print *, myname//': Extracting cloud-water condensate mixing ratio from background' 
        allocate ( cwmr(w_f%grid%im,w_f%grid%jm,w_f%grid%km) )
        where ( abs(w_f%q(:,:,:,3)-undef)>0.001 ) 
                cwmr(:,:,:) = reshape(w_f%q(:,:,:,3),(/w_f%grid%im,w_f%grid%jm,w_f%grid%km/))
        end where
   endif

   call myout ( 6, myname, 'Initialize ss vector' )
   if ( eta ) then    ! even with eta coordinates, nsig is the number of layers
        nsig = size(w_f%grid%ak)-1
        call ss_init( ssiupaf, ss, rc, &
                 jcap, nsig, nymd, nhms, fhour=fhour, ak=w_f%grid%ak, bk=w_f%grid%bk )
   else if ( hybrid ) then
        allocate(ako(nsig+1),bko(nsig+1))
        print *, myname//': setting ak/bk values for ',nsig,' levels'
        call set_eta( nsig, ks, ptop, pint, ako, bko )
        call ss_init( ssiupaf, ss, rc, &
                 jcap, nsig, nymd, nhms, fhour=fhour, ak=ako, bk=bko )       
   else
        call ss_init( ssiupaf, ss, rc, &
                 jcap, nsig, nymd, nhms, fhour=fhour )
   end if
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot initialize ss vector', rc)
   end if

   call myout ( 6, myname, 'Transform fv vector to ss vector...please wait' )
   if ( ncep_phis ) then
        if ( docw ) then
           call xform( w_f, ss, nymd, nhms, sfcfvfn, ssisfcf, rc, cwmr=cwmr, vectype=vectype, nhybrid=hybrid )
        else 
           call xform( w_f, ss, nymd, nhms, sfcfvfn, ssisfcf, rc, vectype=vectype, nhybrid=hybrid )
        endif
   else
        allocate ( phis(w_f%grid%im,w_f%grid%jm) )
        phis = w_f%phis    
        if ( docw ) then
          call xform( w_f, ss, nymd, nhms, sfcfvfn, ssisfcf, rc, fvphis=phis, cwmr=cwmr, eta=eta, vectype=vectype, nhybrid=hybrid )
        else
          call xform( w_f, ss, nymd, nhms, sfcfvfn, ssisfcf, rc, fvphis=phis, eta=eta, vectype=vectype, nhybrid=hybrid )
        endif
        deallocate ( phis )
   end if
   if ( rc .ne. 0 ) then
     call myexit ( myname, 'could not transform fv vector to ss vector', rc )
   end if

   if (use_sigio) then
      call myout ( 6, myname, 'Write out restart using sigio_module' )
      call ss_put_sigio( ssiupaf, ss, rc )
   else
      call myout ( 6, myname, 'Write out restart ' )
      call ss_put( ssiupaf, ss, rc )
   endif
   if ( rc .ne. 0 ) then
     call myexit ( myname, 'could not write out ss vector to a file', rc )
   end if
 
! clean up

   call myout ( 6, myname, 'Clean up ' )
   call ss_clean(ss)
   if(docw) deallocate ( cwmr )
   call dyn_clean(w_f)

! done
   
   call myout ( 6, myname, ' -- fv2ss.x has successfully ended -- ' )
   call exit(0)

   CONTAINS



!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !IROUTINE: init --- Initialize fv2ss
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
   subroutine init ( etafile, sfcfvfn, ssiupaf, ssisfcf, jcap, nsig, fhour, pick, & 
                     nymd, nhms, ncep_phis, use_sigio, vectype, hybrid, rc )
!
! !USAGE:
!
! !USES:
!
   implicit NONE
!
   character*255, intent(out) :: etafile
   character*255, intent(out) :: sfcfvfn
   character*255, intent(out) :: ssiupaf
   character*255, intent(out) :: ssisfcf
   integer, intent(out) :: jcap
   integer, intent(out) :: nsig
   integer, intent(out) :: fhour
   integer, intent(out) :: rc
   logical, intent(out) :: pick       ! if true,
                                      ! (nymd,nhms) are input
                                      ! parameters to dyn_get.
   integer, intent(out) :: nymd, nhms ! date/time for first guess
                                      ! if pick = .true.
   logical, intent(out)   :: ncep_phis
   logical, intent(out)   :: use_sigio
   integer, intent(out) :: vectype    ! 5 for GEOS-5; 4 for GEOS-4
   logical, intent(out) :: hybrid

!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  2002.10.30  C. Cruz:   Initial Code
!  2003.11.21  Todling    Added option to skip using ncep phis'
!  2003.12.08  Todling    Added option to set forecast hour (fhour)
!  2009.06.22  Todling    Add GEOS-5 (new vec) support
!  2009.09.25  Sienkiewicz Flag to use sigio_module for output
!  2009.10.07  Sienkiewicz 'hybrid' flag changes output from sigma
!
!-------------------------------------------------------------------------
!

   character*4, parameter :: myname = 'init'
   integer :: n,iargc,nargs
   character*255 :: argv
   character*255, allocatable :: arg(:)

   nargs =  iargc()
   if( nargs.eq.0 ) call usage()
   allocate ( arg(nargs) )

! defaults

   ssiupaf = 'fv2ss_sig_out.dat'
   ssisfcf = 'sfcf06'
   jcap = 62
   nsig = 0
   pick = .false.
   ncep_phis = .false.
   use_sigio = .false.
   hybrid = .false.
   fhour = 6
   vectype = 4

   do n=1,nargs
     call getarg(n,arg(n))
   enddo
   do n=1,nargs
     if( trim(arg(n)).eq.'-o'    ) ssiupaf = trim(arg(n+1))
     if( trim(arg(n)).eq.'-os'   ) ssisfcf = trim(arg(n+1))
     if( trim(arg(n)).eq.'-jcap' ) read(arg(n+1), * ) jcap
     if( trim(arg(n)).eq.'-nsig' ) read(arg(n+1), * ) nsig
     if( trim(arg(n)).eq.'-fhour') read(arg(n+1), * ) fhour
     if( trim(arg(n)).eq.'-ncep_phis' ) ncep_phis = .true.
     if( trim(arg(n)).eq.'-g5'   ) vectype = 5
     if( trim(arg(n)).eq.'-sigio') use_sigio = .true.
     if( trim(arg(n)).eq.'-pick' ) then
        read(arg(n+1),*) nymd
        read(arg(n+2),*) nhms
        pick = .true.
     end if
     if( trim(arg(n)).eq.'-hybrid') hybrid = .true.
   enddo
   etafile = trim(arg(1))
   sfcfvfn = trim(arg(2))

   if (nsig==0) ncep_phis = .false.

!  Echo the parameters
!  -------------------
   print *
   print *, '    fv2ss.x : convert eta fields to spectral'
   print *, '--------------------------------------------------'
   print *, ' fv_eta   filename     : ',trim(etafile)
   print *, ' fv_sfc   filename     : ',trim(sfcfvfn)
   print *, ' ana upa  filename     : ',trim(ssiupaf)
   print *, ' ana sfc  filename     : ',trim(ssisfcf)
   print *, ' triangular truncation : ',jcap
   if (nsig>0) then
      if (hybrid) then
         print *, ' No. hybrid layers     : ',nsig
      else
         print *, ' No. sigma layers      : ',nsig
      end if
   else 
      print *, ' ***** Inherit fv layers *****'
   end if
   if (pick) then
   print *, ' NHMS                  : ',nhms
   print *, ' NYMD                  : ',nymd
   end if
   if (ncep_phis) then
   print *, ' Using NCEP Topography '
   else
   print *, ' Not using NCEP Topography (fv-instead)'
   endif
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
   print *, " fv2ss.x fv_eta_fname fv_sfc_fname [...options...]"
   print *
   print *, "options :"
   print *
   print *, " [-o  ssiupafile  ] (default: fv2ss_sig_out.dat)"
   print *, " [-os ssisfcfile  ] (default: sfcf06)"
   print *, " [-jcap truncation] (default: 170)"
   print *, " [-nsig # sigma   ] (default: 0; see below)"
   print *, " [-fhour          ] (default: 6)"
   print *, " [-pick nymd nhms ] (default: last)"
   print *, " [-ncep_phis      ] "
   print *, " [-g5             ] "
   print *, " [-sigio          ] "
   print *, " [-hybrid         ] "
   print *
   print *, "where:"
   print *
   print *, " fv_eta_fname     : Filename of FvDAS dynamics in eta format"
   print *, " fv_sfc_fname     : Filename of FvDAS surface fields"
   print *
   print *, " -o  ssiupafile   : Filename of NCEP sigma-level analysis in spectral space"
   print *, " -os ssisfcfile   : Filename of NCEP surface     analysis field (guassian)"
   print *, " -jcap truncation : model triangular truncation"
   print *, " -nsig  #sigma    : number of sigma layers"
   print *, "             NOTE : The default is to inherit fv eta coordinates,"
   print *, "                    using fv topography."
   print *, " -fhour           : forecast hours (usually 6hr fcst for ana)" 
   print *, " -pick nymd nhms  : date and time to be read from first guess file"
   print *, " -ncep_phis       : to use ncep phis instead of fvgcm - ignored when nsig=0"
   print *, " -g5              : using g5 eta format (default g4 format)"
   print *, " -sigio           : use sigio_module for ss output"
   print *, " -hybrid          : use hybrid coord. (nsig lvls) instead of sigma for output"
   print *
   call exit(1)

   end subroutine usage
 
   end program fv2ss

