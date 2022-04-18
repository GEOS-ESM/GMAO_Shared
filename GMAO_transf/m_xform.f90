!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_xform --- FVSSI transformation module 

!
! !INTERFACE:
!
   module m_xform

! !USES:

   use m_ss
   use m_gg
   use m_dyn
   use m_fv2fv, only : fv2fv
   use m_fv2ll, only : fv2ll
   use m_ll2fv, only : ll2fv
   use m_ioutil, only : luavail, opnieee, clsieee
   use util
   use m_die, only : die
   use make_surf, only : fv2ncep_surf

   use m_const, only : grav
   use m_const, only : kappa
   use m_const, only : rvap
   use m_const, only : rgas
   use m_const, only : undef
   use m_const, only : o3_ppmv2gpg
   use m_const, only : o3_gpg2ppmv

   use m_dynp, only : dynp_add
   implicit none

   private

!
! !PUBLIC TYPES:
!

! !PUBLIC MEMBER FUNCTIONS:
!
   public xform      ! transforms between FV eta and SSI vectors
   
! !DESCRIPTION: This module contains the methods necessary for 
!               tranforming FV dynamics vectors in eta format to
!               NCEP's SSI vectors in spectral space
!
! !REVISION HISTORY:
!
!  30Sep2002  Cruz     Initial code.
!  15Mar2004  Todling  Added cloud-water capability for the direction fv->ss;
!                      need to see how to best go about adding the other direction.
!  17Dec2004  Kokron   Initialize pointers to null to get myalloc to work on Altix
!  12Jan2005  Todling  Replaced getcon by m_const
!  15Feb2005  Dee      Added eta option (to support spectral-eta vectors)
!  10May2005  Todling  Renamed var cwf to cwmr for consistency w/ GSI.
!
!EOP
!-------------------------------------------------------------------------

!  Interfaces
!  ----------

   interface xform
     module procedure xform_f2s, xform_s2f
   end interface
  
   interface
     subroutine interp_h ( q_cmp,im,jm,lm, &
                            dlam,dphi,rotation,tilt,precession, &
                            q_geo,irun,lon_geo,lat_geo, &
                            msgn,norder,check,undef )
     integer im,jm,lm,irun,norder,msgn
     logical check
     real      q_geo(irun,lm)
     real    lon_geo(irun)
     real    lat_geo(irun)
     real    q_cmp(im,jm,lm)
     real     dlam(im)
     real     dphi(jm)
     end subroutine
   end interface
  
   interface
     subroutine get_ncep_phis ( filename,phis,im,jm,err )
     integer     im,jm,err
     real   phis(im,jm)
     character(len=*) filename
     end subroutine
   end interface

   interface
     subroutine get_ncep_tsx( filename,ts,im,jm,nymd,nhms )
     integer  im,jm,nymd,nhms
     real  ts(im,jm)
     character*80 filename
     end subroutine
   end interface

   logical, parameter :: ncep_phis_out = .false. ! in case need to get NCEP phis

!---------------------------------------------------------------------------

   contains

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xform_f2s: transform FV to SSI vector

!
! !INTERFACE:
!
   subroutine xform_f2s ( fv, ss, nymd, nhms, sfcfvfn, sfcncfn, rc, &
                          cwmr, fvphis, eta, vectype, nhybrid  )   ! optionals
!
! !INPUT PARAMETERS:
!
   type(dyn_vect), intent(in)          :: fv      ! FV vector
   integer, intent(in)                 :: nymd    ! year-month-day
   integer, intent(in)                 :: nhms    ! hour-min-sec
   character(len=*), intent(in)        :: sfcfvfn ! bkg.sfc filename
   character(len=*), intent(in)        :: sfcncfn ! sfcfHH (NCEP-like) sfc filename

   real, intent(in), dimension(:,:),   optional :: fvphis ! fv-provided phis
   real, intent(in), dimension(:,:,:), optional :: cwmr   ! fv-cloud-water condensate mixing ratio

   logical, intent(in), optional       :: eta      ! if present and .true., ss is spectral-eta
   integer, intent(in), optional       :: vectype  ! differentiate GEOS-4 and GEOS-5 vectors
   logical, intent(in), optional       :: nhybrid  ! if present and .true. 
!
! !OUTPUT PARAMETERS:
!
   type(ss_vect),intent(inout)         :: ss   ! SSI vector
   integer, intent(out)                :: rc   ! return code
!
! !DESCRIPTION: This routine transforms an FV vector into an SSI vector
!
! !REVISION HISTORY:
!
!  05Nov2002  Cruz     Initial code.
!  21Nov2003  Todling  Modified handling of phis
!  18Dec2003  Cruz     Call to dtoa(v) had bad last argument in call
!  15Mar2004  Todling  Added cloud-water capability
!  15Mar2004  Kokron   Bug fix: dtoa was passing fvgg instead of fv.
!  16Mar2004  Kokron/RT Bug fix: initialization of fvgg w/ wrong #levs.
!  21Oct2004  Todling  Asyn fix: added output filename for surface fields.
!  15Feb2005  Dee      Added eta option (bypass vertical remap)
!  22Jun2009  Todling  Properly treat GEOS-5 vectors
!  30Sep2009  Sienkiewicz  1st grid longitude determines rotation in ll2gg
!   7Oct2009  Sienkiewicz  seta->fvseta (our levels,no remapping); add hybrid
!                           grid option using different (NCEP 64) levels with
!                           remap.
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'xform_f2s'
   type(gg_vect)               :: gg   ! gaussian grid vector
   type(gg_vect)               :: fvgg ! fv vector on a gaussian grid 
   real, pointer               :: phis(:,:) => null()
   real, pointer               :: ncep_phis(:,:) => null()
   real, pointer               :: akbk8(:,:) => null()
   integer                     :: id,err,k,L
   integer                     :: rlon, rlat, mlev, glon, glat, nsig
   logical 		       :: docwmr
   logical                     :: g5vec=.false.
   logical                     :: havetv=.false.
   logical 		       :: fvseta = .false.
   logical                     :: hybrid = .false.
   real                        :: eps
   real, allocatable           :: pe(:,:,:), pke(:,:,:), pk(:,:)
   real                        :: rot
! start

   rc=0
   rlon = fv%grid%im
   rlat = fv%grid%jm
   mlev = fv%grid%km
   glon = ss%meta%lonf
   glat = ss%meta%latf
   nsig = ss%meta%nsig
   docwmr= present(cwmr)
   eps  = rvap/rgas-1.0
   
   if (present(eta)) fvseta = eta
   if (present(nhybrid)) hybrid = nhybrid
    
!  Check dim's of optional arguments
!  ---------------------------------
   if (present(fvphis)) then
       if (size(fvphis,1)/=rlon .and. size(fvphis,2)/=rlat) then
           rc = 1
           return
       endif
   endif
   if (present(cwmr)) then
       if (size(cwmr,1)/=rlon .and. size(cwmr,2)/=rlat .and. size(cwmr,3)/=mlev) then
           rc = 2
           return
       endif
   endif
   if (present(vectype)) then
      if(vectype==5) then
         g5vec=.true.
         havetv=.true.
      endif
   endif
      
   ! put winds on 'A' grid before remapping
   
   if(g5vec)then
      print *, "Winds assumed to already be on A-grid"
   else  
      call dtoa ( fv%u,fv%u,rlon,rlat,mlev,2 )
      call dtoa ( fv%v,fv%v,rlon,rlat,mlev,1 )
   endif

   call gg_init ( nsig, rlon, rlat, fvgg, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems initializing gg vector' ) 
     return
   end if

   if (fvseta) then  ! ss will be fv-spectral-eta: no vertical remap
   
!     must have FV topography

      if (.not. present(fvphis)) then
         rc = 3
         return
      endif
      
      fvgg%hs = fvphis
      fvgg%ps = fv%ps
      fvgg%dp = fv%delp     
      fvgg%u  = fv%u
      fvgg%v  = fv%v
      fvgg%q  = fv%q(:,:,:,1)
      fvgg%o  = fv%q(:,:,:,2)  
      if(present(cwmr)) fvgg%cwmr  = cwmr
      
!     compute dry temperature (NOTE: converted to virtual temperature in gg2ss)
      
      if (g5vec) then

          do L=1,mlev ! convert virtual temperature to dry temperature
             fvgg%t(:,:,L) = fv%pt(:,:,L)/( 1.0 + eps*fvgg%q(:,:,L) )
          end do

      else ! < g5vec >

          allocate(pe(rlon,rlat,mlev+1),pke(rlon,rlat,mlev+1),pk(rlon,rlat),stat=err)
          if (err .ne. 0) then
             rc = 4
             call mywarn ( myname, 'problem allocating pe, pke, pL' ) 
             return
          end if
         
          pe(:,:,mlev+1) = fvgg%ps
          pke(:,:,mlev+1) = pe(:,:,mlev+1)**kappa
          do L=mlev,1,-1
             pe(:,:,L) = pe(:,:,L+1) - fvgg%dp(:,:,L)
             pke(:,:,L) = pe(:,:,L)**kappa
          end do
          do L=1,mlev
             pk = ( pke(:,:,L+1)-pke(:,:,L) )/( kappa*log(pe(:,:,L+1)/pe(:,:,L)) )
             fvgg%t(:,:,L) = fv%pt(:,:,L)*pk/( 1.0 + eps*fvgg%q(:,:,L) )
          end do
      
          deallocate(pe,pke,pk)

      endif ! < g5vec >

   else            ! ss will be spectral-sigma: vertical remap
   
!     we need NCEP topography for remapping

      call myalloc(ncep_phis, rlon, rlat)
      call myalloc(     phis, rlon, rlat)
      if (present(fvphis)) then
         fvgg%hs = fvphis
      else
         call get_ncep_phis ( 'ncep_out.dat', ncep_phis, rlon, rlat, err )
         if (err .ne. 0) then
            rc = err
            call mywarn ( myname, 'problem reading NCEP phis')
            return
         endif
         fvgg%hs = ncep_phis
      endif
   
!     we also need NCEP vertical grid

      call myalloc(akbk8,nsig+1,2)


      if (.not. hybrid) then

!     for sigma, ordering of levels is NCEPs, so swap vertical order temporarily to do remapping
         call myout ( 6, myname,'swap order of sig interface values' )
         akbk8(1:nsig+1,2) = ss%meta%si(nsig+1:1:-1)
         akbk8(1:nsig+1,1) = 0.0
      else

!     for hybrid, we're using the GMAO ordering (matching the fv-72eta in old ss format)
         akbk8(1:nsig+1,1) = ss%meta%ak(1:nsig+1)
         akbk8(1:nsig+1,2) = ss%meta%bk(1:nsig+1)
      endif

!     remap based on new vertical grid (ncep_si) and topography (phis) 

      call myout ( 6, myname,'remap' )
      if (present(cwmr)) then
       call fv2ll ( fvgg%ps, fvgg%dp, fvgg%u, fvgg%v, fvgg%t, fvgg%q, fvgg%o, &
                    fvgg%hs, nsig, akbk8, &
                    fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                    fv%q(:,:,:,1), &
                    fv%q(:,:,:,2), &
                    fv%phis, mlev, rlon, rlat, havetv, hybrid, &
                    w1=fvgg%cwmr, w2=cwmr )
      else
       call fv2ll ( fvgg%ps, fvgg%dp, fvgg%u, fvgg%v, fvgg%t, fvgg%q, fvgg%o, &
                    fvgg%hs, nsig, akbk8, &
                    fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                    fv%q(:,:,:,1), &
                    fv%q(:,:,:,2), &
                    fv%phis, mlev, rlon, rlat, havetv, hybrid )
      endif
   
!     call myout ( 6, myname,'swap order of sig interface values' )
!     ss%meta%si(nsig+1:1:-1) = si8(1:nsig+1)
      
      deallocate ( ncep_phis, phis, akbk8 ) 
      
   end if
   
! regrid

   call gg_init ( nsig, glon, glat, gg, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems initializing gg vector' ) 
     return
   end if

   call myout ( 6, myname,'regrid' )
   rot = -fv%grid%lon(1)             ! rotation needed for fv-grid

   call ll2gg ( fvgg, gg, rlon, rlat, nsig, glon, glat, docwmr, rot, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems regridding' ) 
     return
   end if
   
! now that we have the gg vector make the SSI surface restart

   call myout ( 6, myname,'Create SSI surface restart' )
   call fv2ncep_surf (gg, ss, nymd, nhms, sfcfvfn, sfcncfn, rc)
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems in fv2ncep_surf' ) 
     return
   end if

! now transform gaussian to spectral via FFT and Legendre transforms

   call myout ( 6, myname,'transform' )
   call gg2ss ( gg, ss, docwmr, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems in gg2ss' ) 
     return
   end if
 
   call gg_clean ( gg )  
   call gg_clean ( fvgg ) 
   
   end subroutine xform_f2s

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xform_s2f: transform from SSI to FV

!
! !INTERFACE:
!
   subroutine xform_s2f ( nymd, nhms, pick,        &
                          ss, fv, ak, bk, w_f, rc, &
                          w_fs, z_f, uafile,       &
                          pabove, pbelow, cwmr, eta, vectype )

!
! !INPUT PARAMETERS:
!
   logical, intent(in) :: pick
   real, dimension(:), intent(in) :: ak   ! vertical grid coefficients
   real, dimension(:), intent(in) :: bk   ! vertical grid coefficients

!    pbelow ... pressure below which analysis is used completely
!    pabove ... pressure above which model    is used completely
!               Note: a blend is used in-between pbelow and pabove
!               If pbelow=pabove, blending code is disabled
   real, intent(in), optional  :: pabove, pbelow 

   real, intent(in), dimension(:,:,:), optional :: cwmr   ! fv-cloud-water condensate mixing ratio
   logical, intent(in), optional       :: eta  ! if present and .true., ss is spectral-eta
   integer, intent(in), optional       :: vectype  ! differentiate GEOS-4 and GEOS-5 vectors
!
! !INPUT/OUTPUT PARAMETERS:
!

   integer, intent(inout) :: nymd, nhms
   
   type(dyn_vect), intent(inout)  :: fv   ! on output, contains final eta analysis
   type(dyn_vect), intent(inout)  :: w_f  ! eta background state vector
   type(ss_vect) , intent(inout)  :: ss   ! sigma spectral (analysis) state vector

   type(dyn_vect), intent(inout), optional :: z_f  ! background transformed over and back
   type(dyn_vect), intent(inout), optional :: w_fs ! remapped w_f

   character(len=*), intent(in), optional :: uafile ! filename of state vector to blend with

!
! !OUTPUT PARAMETERS:
!
   integer, intent(out)          :: rc ! return code
!
! !DESCRIPTION: This routine transforms a SSI vector to an FV vector 
!               Note this routine is configured the GMAO local format
!               of spectral (eta) files and is incompatible with the
!               NCEP hybrid grids written by 'sigio_module'
!
! !REVISION HISTORY:
!
!  05Nov2002  Cruz     Initial code.
!  15Jun2004  Todling  Turned pabove/pblow intent(in) only.
!  01Sep2004  Todling  Added shaving from hermes.
!  09Sep2004  Todling  Setting fv%phis to w_f%phis; passed zero-diff test.
!  25Oct2004  Takacs   Bug fix: added dtoa to handle fvua (blending field)
!  15Feb2005  Dee      Added eta option (bypass vertical remap)
!  11May2005  Todling  For now, cwmr does not feed back to model, therefore bkg 
!                      cloud condensate info simply copied to analysis field
!  24Jan2006  B Zhang  Updated to put cwmr in ss back to fv
!  22Jun2009  Todling  Properly treat GEOS-5 vectors
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter         :: myname = 'xform_s2f'
   
   type(dyn_vect) :: fvua    ! state vector to blend with
   type(dyn_vect) :: z_as    ! remapped analysis state on analysis grid
   type(dyn_vect) :: z_fs    ! remapped z_f (hold original analyzed field)
   type(gg_vect)  :: gg      ! guassian-grid vector
   type(gg_vect)  :: ll      ! holds gg vector on a lat-lon grid 

!  define arrays to hold 'safe' quantities

   real, allocatable ::  psa(:,:)
   real, allocatable :: dpsa(:,:,:)

   real    :: pa, pb, eps
   integer :: k, err, id, lu, ier, L, li
   integer :: rlon, rlat, mlev, glon, glat, nsig
   integer, parameter :: IGRID = 0 ! integer grid identifier that describes
                                   ! spectral to grid conversion:
                                   ! =4 for gaussian grid
                                   ! =0 for equally-spaced grid including poles,
                                   ! =256 for equally-spaced grid excluding poles
   logical :: seta  = .false.
   logical :: g5vec = .false.
   logical                     :: docwmr
   real, allocatable           :: pe(:,:,:), pke(:,:,:), pk(:,:)
   
! start

   eps  = rvap/rgas-1.0
   if (present(eta)) seta = eta
   if (present(vectype)) then
      if(vectype==5) g5vec=.true.
   endif
   rc=0
   rlon = fv%grid%im
   rlat = fv%grid%jm
   mlev = fv%grid%km
   ! if IGRID ==0 then ss2gg converts the spectral state to a state
   ! on a latlon grid, else to a gaussian grid
   if(IGRID==0) then
     glon = rlon
     glat = rlat
   else
     glon = ss%meta%lonf
     glat = ss%meta%latf
   end if
   docwmr= present(cwmr)
   nsig = ss%meta%nsig

   if (.not. seta) then
      call myout ( 6, myname,'swap order of sig interface values' )
      ss%meta%si(1:nsig+1) = ss%meta%si(nsig+1:1:-1)
   end if

   call gg_init ( nsig, glon, glat, gg, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems initializing gg vector' ) 
     return 
   end if

! convert spectral vector to a gaussian or latlon vector
 
   call myout ( 6, myname, 'transform' )
   call ss2gg ( ss, gg, glon, glat, docwmr, rc, eta=seta )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems in ss2gg' ) 
     return
   end if

   if (.not. seta) then
      call myout ( 6, myname,'swap order of sig interface values' )
      ss%meta%si(1:nsig+1) = ss%meta%si(nsig+1:1:-1)
   end if

   call gg_init ( nsig, rlon, rlat, ll, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems initializing ll vector' ) 
     return
   end if

! regrid gaussian or latlon vector to lat-lon
   
   call myout ( 6, myname, 'regrid' )
   call gg2ll ( ss, gg, rlon, rlat, IGRID, ll, docwmr, rc )
   if ( rc .ne. 0 ) then
     call mywarn ( myname, 'problems in gg2ll' ) 
     return
   end if

!  Structure not needed beyond this point, clear memory ...

   call gg_null ( gg )

   fv%phis = w_f%phis
   if ( ncep_phis_out ) then
        lu = luavail()
        open(lu,file='ncep_phis_out.dat',form='unformatted')
        write(lu) ll%hs
        close(lu)
   endif

!  Remap based on fv vertical grid

   if ( present (uafile) .and. .not. seta ) then
     
   ! if so desired read an eta field to be used for upper air fields  
   ! this eta file will typically be the same as w_ff otherwise it should
   ! the same number of time levels as w_ff
     
     call myout ( 6, myname, 'get eta fields for upper air blending' )
     if ( pick ) then
       call dyn_get ( uafile, nymd, nhms, fvua, err, timidx=0, vectype=vectype )
     else
       call dyn_get ( uafile, nymd, nhms, fvua, err, vectype=vectype )
     end if
     if ( err .ne. 0 ) then
       call myexit ( myname,'cannot read dynamics vector file', rc )
     else
       print *, myname//': read dynamics state file '//trim(uafile)
       print *, myname//': nymd, nhms = ', nymd, nhms
       print *
     end if
      
! put winds on 'A' grid before blending

     if (g5vec) then
        print *, "Winds assumed to already be on A-grid"
     else
        call dtoa ( fvua%u,fvua%u,rlon,rlat,mlev,2 )
        call dtoa ( fvua%v,fvua%v,rlon,rlat,mlev,1 )
     endif

     pa = pabove*100.0 ! convert to hPa
     pb = pbelow*100.0 ! convert to hPa 
     call myout ( 6, myname, 'remap with UA blending' )
     if (present(cwmr)) then
     call ll2fv ( fvua%ps, fvua%delp, fvua%u, fvua%v, fvua%pt, & 
                  fvua%q(:,:,:,1), &
                  fvua%q(:,:,:,2), &
                  fv%phis, mlev, &
                  ll%ps, ll%dp, ll%u, ll%v, ll%t, ll%q, ll%o, ll%hs, &
                  nsig, rlon, rlat, pb, pa, &
                  w1=fvua%q(:,:,:,3), w2=ll%cwmr )
     else
     call ll2fv ( fvua%ps, fvua%delp, fvua%u, fvua%v, fvua%pt, & 
                  fvua%q(:,:,:,1), &
                  fvua%q(:,:,:,2), &
                  fv%phis, mlev, &
                  ll%ps, ll%dp, ll%u, ll%v, ll%t, ll%q, ll%o, ll%hs, &
                  nsig, rlon, rlat, pb, pa )
     endif
     fv%ps   = fvua%ps
     fv%delp = fvua%delp
     fv%pt   = fvua%pt
     fv%u    = fvua%u
     fv%v    = fvua%v
     do li = 1, fv%grid%lm
        fv%q(:,:,:,li) = fvua%q(:,:,:,li)
     enddo
                         ! take care of variables not "analyzed"
     fv%ts      = fvua%ts
     fv%lwi     = fvua%lwi
     fv%hs_stdv = fvua%hs_stdv
     
     call dyn_null ( fvua )
     
   else if ( .not. seta ) then    ! this is not desirable
     
     pa = 0.00
     pb = pa
     call myout ( 6, myname, 'remap' )

     if (present(cwmr)) then
     call ll2fv ( fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                  fv%q(:,:,:,1), &
                  fv%q(:,:,:,2), &
                  fv%phis, mlev, &
                  ll%ps, ll%dp,   ll%u, ll%v, ll%t,  ll%q, ll%o, ll%hs, &
                  nsig, rlon, rlat, pb, pa, & 
                  w1=fvua%q(:,:,:,3), w2=ll%cwmr )
     else
     call ll2fv ( fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                  fv%q(:,:,:,1), &
                  fv%q(:,:,:,2), &
                  fv%phis, mlev, &
                  ll%ps, ll%dp,   ll%u, ll%v, ll%t,  ll%q, ll%o, ll%hs, &
                  nsig, rlon, rlat, pb, pa )
     endif
   else  ! spectral-eta case
     
      fv%phis       = ll%hs
      fv%ps         = ll%ps
      fv%delp       = ll%dp   
      fv%u          = ll%u
      fv%v          = ll%v
      fv%q(:,:,:,1) = ll%q
      fv%q(:,:,:,2) = ll%o 
      if (present(cwmr)) then
      fv%q(:,:,:,3) = ll%cwmr
      endif
      
!     compute virtual potential temperature
      
      allocate(pe(rlon,rlat,mlev+1),pke(rlon,rlat,mlev+1),pk(rlon,rlat),stat=err)
      if (err .ne. 0) then
         rc = 4
         call mywarn ( myname, 'problem allocating pe, pke, pk' ) 
         return
      end if
         
      pe(:,:,mlev+1) = fv%ps
      pke(:,:,mlev+1) = pe(:,:,mlev+1)**kappa
      do L=mlev,1,-1
         pe(:,:,L) = pe(:,:,L+1) - fv%delp(:,:,L)
         pke(:,:,L) = pe(:,:,L)**kappa
      end do
      do L=1,mlev
         pk = ( pke(:,:,L+1)-pke(:,:,L) )/( kappa*log(pe(:,:,L+1)/pe(:,:,L)) )
         fv%pt(:,:,L) = ll%t(:,:,L)*( 1.0 + eps*fv%q(:,:,L,1) )/pk
      end do
      
      deallocate(pe,pke,pk)
     
   end if

!  Structure not needed beyond this point, clear memory ...

   call gg_null ( ll )

!  Put u,v winds on 'D' grid

   if (g5vec) then
      print *, "Keeping winds on A-grid"
   else
      call atod ( fv%u,fv%u,rlon,rlat,mlev,2 )
      call atod ( fv%v,fv%v,rlon,rlat,mlev,1 )
   endif

   if ( present (w_fs) .and. present(z_f) .and. .not. seta ) then
   
! Calculate 'safe' pressure and delta pressure
! --------------------------------------------
     allocate( psa(rlon,rlat),     stat=ier); if(ier/=0) call die(myname,'error in alloc(psa )')
     allocate(dpsa(rlon,rlat,mlev),stat=ier); if(ier/=0) call die(myname,'error in alloc(dpsa)')
     psa  = w_f%ps   + ( fv%ps   - z_f%ps   )
     dpsa = w_f%delp + ( fv%delp - z_f%delp )
   
!  Remap background w_f to safe grid
!  ---------------------------------
     call myout ( 6, myname, 'remap w_f state to analysis grid w/ psafe' )
     if (present(cwmr)) then
     call fv2fv ( psa, dpsa, w_fs%u, w_fs%v, w_fs%pt, &
                  w_fs%q(:,:,:,1), &
                  w_fs%q(:,:,:,2), &
                  fv%phis, mlev, &
                  w_f%ps, w_f%delp, w_f%u, w_f%v, w_f%pt, &
                  w_f%q(:,:,:,1), &
                  w_f%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat, &
                  w1=w_fs%q(:,:,:,3),w2=w_f%q(:,:,:,3)) 
     else
     call fv2fv ( psa, dpsa, w_fs%u, w_fs%v, w_fs%pt, &
                  w_fs%q(:,:,:,1), &
                  w_fs%q(:,:,:,2), &
                  fv%phis, mlev, &
                  w_f%ps, w_f%delp, w_f%u, w_f%v, w_f%pt, &
                  w_f%q(:,:,:,1), &
                  w_f%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat )
     endif
     w_fs%ps = psa 
     w_fs%delp = dpsa 
     w_fs%phis = fv%phis
       
     call myout ( 6, myname, 'initialize z_fs to hold original analyzed field' )
     call dyn_init ( rlon, rlat, mlev, fv%grid%lm, z_fs, rc , &
                     ptop = fv%grid%ptop, &
                     ks   = fv%grid%ks,   &
                     ak   = fv%grid%ak,   &
                     bk   = fv%grid%bk    )

!  Remap transformed over and back background z_f to safe grid
!  -----------------------------------------------------------
     call myout ( 6, myname, 'remap z_f state to analysis grid w/ psafe' )
     if (present(cwmr)) then
     call fv2fv ( psa, dpsa, z_fs%u, z_fs%v, z_fs%pt, &
                  z_fs%q(:,:,:,1), &
                  z_fs%q(:,:,:,2), &
                  fv%phis, mlev, &
                  z_f%ps, z_f%delp, z_f%u, z_f%v, z_f%pt, &
                  z_f%q(:,:,:,1), &
                  z_f%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat, &
                  w1=z_fs%q(:,:,:,3),w2=z_f%q(:,:,:,3)) 
     else
     call fv2fv ( psa, dpsa, z_fs%u, z_fs%v, z_fs%pt, &
                  z_fs%q(:,:,:,1), &
                  z_fs%q(:,:,:,2), &
                  fv%phis, mlev, &
                  z_f%ps, z_f%delp, z_f%u, z_f%v, z_f%pt, &
                  z_f%q(:,:,:,1), &
                  z_f%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat )
     endif
     z_fs%ps   = psa 
     z_fs%delp = dpsa 
     z_fs%phis = fv%phis


     call myout ( 6, myname, 'initialize z_as ...' )
     call dyn_init ( rlon, rlat, mlev, fv%grid%lm, z_as, rc , &
                     ptop = fv%grid%ptop, &
                     ks   = fv%grid%ks,   &
                     ak   = fv%grid%ak,   &
                     bk   = fv%grid%bk    )

 
!  Remap analysis z_as to safe grid
!  --------------------------------
     call myout ( 6, myname, 'remap z_as state to analysis grid w/ psafe' )
     if (present(cwmr)) then
     call fv2fv ( psa, dpsa, z_as%u, z_as%v, z_as%pt, & 
                  z_as%q(:,:,:,1), &
                  z_as%q(:,:,:,2), &
                  fv%phis, mlev, &
                  fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                  fv%q(:,:,:,1), &
                  fv%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat, &
                  w1=z_as%q(:,:,:,3),w2=fv%q(:,:,:,3))
     else
     call fv2fv ( psa, dpsa, z_as%u, z_as%v, z_as%pt, & 
                  z_as%q(:,:,:,1), &
                  z_as%q(:,:,:,2), &
                  fv%phis, mlev, &
                  fv%ps, fv%delp, fv%u, fv%v, fv%pt, &
                  fv%q(:,:,:,1), &
                  fv%q(:,:,:,2), &
                  fv%phis, mlev, rlon, rlat )
     endif
     z_as%ps = psa 
     z_as%delp = dpsa 
     z_as%phis = fv%phis

!  Construct analyzed state vector by adding incremental diff between remapped
!  analysis and remapped transformed over and back background state
!  ---------------------------------------------------------------------------
     fv%ps   = psa 
     fv%delp = dpsa
     fv%phis = fv%phis 
     fv%u    = w_fs%u
     fv%v    = w_fs%v
     fv%pt   = w_fs%pt
     fv%q    = w_fs%q
                                         ! from here on, use z_as as "dw"
     z_as%phis = z_as%phis - z_fs%phis
     z_as%ts   = z_as%ts   - z_fs%ts

     z_as%u    = z_as%u    - z_fs%u
     z_as%v    = z_as%v    - z_fs%v
     z_as%pt   = z_as%pt   - z_fs%pt
     z_as%q    = z_as%q    - z_fs%q
    
     call dyn_null ( z_fs )             ! clean up

     call dynp_add ( fv, 1.0, z_as, rc, verbose=.true. ) 

! clean up

     call dyn_clean ( z_as ) 
     deallocate (dpsa)
     deallocate ( psa)

   elseif ( present (w_fs) .and. present(z_f) .and. seta ) then
   
     fv%phis  = fv%phis  + (w_f%phis - z_f%phis )
     fv%ps    = fv%ps    + (w_f%ps   - z_f%ps   )
     fv%ts    = fv%ts    + (w_f%ts   - z_f%ts   )
     fv%u     = fv%u     + (w_f%u    - z_f%u    )
     fv%v     = fv%v     + (w_f%v    - z_f%v    )
     fv%delp  = fv%delp  + (w_f%delp - z_f%delp )
     fv%pt    = fv%pt    + (w_f%pt   - z_f%pt   )
     fv%q     = fv%q     + (w_f%q    - z_f%q    )  
     do li = 3, fv%grid%lm
        fv%q(:,:,:,li) = w_f%q(:,:,:,li)   ! equal to background for now
     enddo

   end if

   end subroutine xform_s2f

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ss2gg - transform SS vector to GG vector

!
! !INTERFACE:
!
   subroutine ss2gg ( ss, gg, glon, glat, docwmr, rc, &
                                          eta ) ! optional
!
! !USES:
!
   implicit none
!
! !INPUT PARAMETERS:
!
   type(ss_vect),intent(in)  :: ss          ! SSI vector
   integer, intent(in)       :: glon, glat  ! gaussian grid dimensions 
   
   logical, intent(in)         :: docwmr ! controls cloud-water cond mix ratio
   logical, intent(in), optional    :: eta  ! if present and .true., ss is spectral-eta
!
! !OUTPUT PARAMETERS:
!
   type(gg_vect),intent(inout) :: gg  ! gaussian or latlon grid vector
   integer, intent(out)      :: rc  ! return code
!
! !DESCRIPTION: this routine converts a spectral space vector to a gaussian
!               or latlon grid vector. Right now it is set to latlon (see
!               IGRID below).

! !BUGS:

! !REVISION HISTORY:
!
!  01Aug2002 Cruz     Adapted from NCEP to FVSSI
!  26Mar2004 Todling  Zeroing out of O3 when not available in ss
!  27Sep2004 Todling  A bit more memory savvy
!  15Feb2005 Dee      Added eta option
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'ss2gg'
   integer, parameter :: IROMB = 0 ! integer spectral domain shape
                                   ! (0 for triangular, 1 for rhomboidal)
   integer, parameter :: IGRID = 0 ! integer grid identifier
                                   ! spectral to grid conversion:
                                   ! =4 for gaussian grid
                                   ! =0 for equally-spaced grid including poles,
                                   ! =256 for equally-spaced grid excluding poles
   integer, parameter :: IDIR =  1 ! integer transform flag
                                   ! idir>0 for wave to grid, idir<0 for grid to wave

   integer                         :: j,k,jcap,nsig,ntrac,err
   real                            :: con_fvirt
   logical                         :: seta = .false.
   real, pointer :: si8(:) => null()
   real, pointer :: ak(:) => null()
   real, pointer :: bk(:) => null()
   real, pointer :: dumg(:,:,:) => null()
   real, pointer :: dumgv(:,:,:) => null()
  
! Start

   if (present(eta)) seta = eta

! swap north-south pole here and swap vertical levels in ss_get

   rc=0
   jcap=ss%meta%jcap
   nsig=ss%meta%nsig
   ntrac=ss%meta%ntrac
   con_fvirt = 4.6150e+2/2.8705e+2 - 1
   call myalloc(dumg,glon,glat,nsig)
   call myalloc(dumgv,glon,glat,nsig)

! topography

   call myout ( 6, myname, 'topography' )
   call sptez(IROMB,jcap,IGRID,glon,glat,ss%hs,dumg(:,:,1),IDIR)
   gg%hs(:,glat:1:-1) = dumg(:,1:glat,1)

! surface pressure

   call myout ( 6, myname, 'surface pressure' )
   call sptez(IROMB,jcap,IGRID,glon,glat,ss%ps,dumg(:,:,1),IDIR)
   dumg(:,:,1)=exp(dumg(:,:,1))*1.e3
   gg%ps(:,glat:1:-1) = dumg(:,1:glat,1)

! delta pressure

   call myout ( 6, myname, 'delta pressure' )
   
   if (seta) then
      call myalloc(ak,nsig+1)
      call myalloc(bk,nsig+1)
      ak = ss%meta%ak
      bk = ss%meta%bk
      do k=1,nsig
         gg%dp(:,:,k) = (ak(k+1)-ak(k)) + (bk(k+1)-bk(k))*gg%ps
      end do
      deallocate (ak,bk)
   else
      call myalloc(si8,nsig+1)
      si8=ss%meta%si
      do k=1,nsig
         gg%dp(:,:,k) = gg%ps*(si8(k+1)-si8(k))
      end do
      deallocate (si8)
   end if
 
! temperature

   call myout ( 6, myname, 'temperature' )
   call sptezm(IROMB,jcap,IGRID,glon,glat,nsig,ss%t,dumg,IDIR)
   gg%t(:,glat:1:-1,:) = dumg(:,1:glat,:)

! water-vapor mixing ratio 

   call myout ( 6, myname, 'water-vapor mixing ratio' )
   call sptezm(IROMB,jcap,IGRID,glon,glat,nsig,ss%q(:,:,1),dumg,IDIR)
   gg%q(:,glat:1:-1,:) = dumg(:,1:glat,:)
      
! convert virtual T to dry bulb T

   gg%t=gg%t/(1+con_fvirt*gg%q)

! cloud-water fraction
                                                                                                                                                             
   if ( docwmr ) then
   call myout ( 6, myname, 'cloud-water fraction' )
   call sptezm(IROMB,jcap,IGRID,glon,glat,nsig,ss%q(:,:,3),dumg,IDIR)
   gg%cwmr(:,glat:1:-1,:) = dumg(:,1:glat,:)
   endif

! ozone mixing ratio 

   if ( ntrac>=2 ) then
        call myout ( 6, myname, 'ozone' )
        call sptezm(IROMB,jcap,IGRID,glon,glat,nsig,ss%q(:,:,2),dumg,IDIR)
        gg%o(:,glat:1:-1,:) = dumg(:,1:glat,:)
   else
        call myout ( 6, myname, 'ozone zeroed out' )
        gg%o(:,glat:1:-1,:) = 0.0
   endif
      
! winds

   call myout ( 6, myname, 'winds' )
   call sptezmv(IROMB,jcap,IGRID,glon,glat,nsig,ss%d,ss%z,dumg,dumgv,IDIR)
   gg%u(:,glat:1:-1,:) = dumg(:,1:glat,:)
   gg%v(:,glat:1:-1,:) = dumgv(:,1:glat,:)

   deallocate(dumgv)
   deallocate(dumg )

   end subroutine ss2gg

!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gg2ll: regrid gaussian or latlon to latlon

!
! !INTERFACE:
!
   subroutine  gg2ll ( ss, gg, im, jm, igrid, ll, docwmr, rc )
!
! !INPUT PARAMETERS:
!
   type(ss_vect), intent(in)   :: ss    ! spectral space vector
   type(gg_vect), intent(in)   :: gg    ! gaussian grid vector vector
   integer, intent(in)         :: im    ! FV lon-dimension
   integer, intent(in)         :: jm    ! FV lat-dimension
   integer, intent(in)         :: igrid ! flag for regular or gaussian grid 
   logical, intent(in)        :: docwmr     ! controls cloud-water cond mix ratio
!
! !OUTPUT PARAMETERS:
!
   type(gg_vect),intent(inout) :: ll  ! gaussian grid vector on lat-lon grid
   integer, intent(out)        :: rc  ! return code
!
! !DESCRIPTION: routine that regrids NCEP's sigma state vectors on
!               gaussian or latlon grids to the user specified latlon resolution
!
! !BUGS:
!
! !REVISION HISTORY:
!
!  01Sep2002 Cruz     Created
!  30Jul2004 Todling  Updated interface call to comput_gaus.
!  28Sep2004 Todling  First-in, Last-out for memory alloc()
!  13Dec2004 Todling  Bug fix; glat not was set when igrid/=0; 
!                              glon was set incorrectly.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'gg2ll'
   integer :: k,glon,glat,nsig,err
   real, pointer :: dlam(:) => null()
   real, pointer :: dphi(:) => null()
   real, pointer :: lons(:) => null()
   real, pointer :: lats(:) => null()
   real, pointer :: dum(:,:,:) => null()
   real :: tnp,tsp
   integer, parameter :: ONE=1  
   integer, parameter :: SFL=1      ! flag for scalar 
   integer, parameter :: VFL=-1     ! flag for vector 
   integer, parameter :: ORD=3      ! order of interpolation
   real, parameter    :: ROT = 0.0  ! Rotation parameter for lam
   real, parameter    :: TLT = 90.0 ! Rotation parameter for phi 
   real, parameter    :: PRE = 0.0  ! Rotation parameter lam_0
   logical, parameter :: FLG = .true.

! Start

   rc=0
   if(igrid==0) then
     glon=im
     glat=jm
   else
     glon=ss%meta%lonf
     glat=ss%meta%latf
   end if
   write (*,*) 'will regrid from ',glon,'x',glat,' to ',im,'x',jm
   nsig=ss%meta%nsig
  
   call myalloc(dum,glon,glat,nsig) 
   call myalloc(dlam,glon)
   call myalloc(dphi,glat)
   call myalloc(lons,im*jm)
   call myalloc(lats,im*jm)

   if(igrid==0) then
     call setup_grid(glon,glat,im,jm,dphi,dlam,lons,lats)
   else
     call comput_gaus(glon,glat,im,jm,dphi,dlam,lons,lats,.false.)
   end if

! regrid from NCEP's glon*glat gaussian or latlon grid to a
! regular lat-lon grid used by FVDAS
! note that the SSI fields are pole-reversed as well as vertical-reversed
! i.e. the SSI's north pole is FV's south pole and vice versa
! likewise SSI's bottom layer is FV's top.
 
! surface geopotential

   call myout ( 6, myname, 'surface geopotential' )
   dum(:,:,1) = gg%hs*grav
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum(:,:,1),glon,glat,ONE,dlam,dphi,ROT,TLT,PRE,ll%hs, &
                     im*jm,lons,lats,SFL,ORD,.false.,undef)
   else
     ll%hs = dum(:,:,1)
   end if

! surface pressure

   call myout ( 6, myname, 'surface pressure' )
   dum(:,:,1) = gg%ps
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum(:,:,1),glon,glat,ONE,dlam,dphi,ROT,TLT,PRE,ll%ps, &
                     im*jm,lons,lats,SFL,ORD,.false.,undef)
   else
     ll%ps = dum(:,:,1)
   end if

! Pressure Thickness

   call myout ( 6, myname, 'pressure thickness' )
   dum = gg%dp
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%dp, &
                     im*jm,lons,lats,SFL,ORD,FLG,undef)
   else
     ll%dp = dum
   end if

! dry bulb Temperature

   call myout ( 6, myname, 'dry bulb temperature' )
   dum = gg%t
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%t, &
                     im*jm,lons,lats,SFL,ORD,FLG,undef)
   else
     ll%t = dum
   end if

! Specific Humidity

   call myout ( 6, myname, 'specific humidity' )
   dum = gg%q
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%q, &
                     im*jm,lons,lats,SFL,ORD,FLG,undef)
   else
     ll%q = dum
   end if

! Ozone

   call myout ( 6, myname, 'ozone' )
   dum = gg%o
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%o, &
                     im*jm,lons,lats,SFL,ORD,FLG,undef)
   else
     ll%o = dum
   end if
   ! at this point ll%o is in kg/kg, we need to convert to ppmv
   ll%o = o3_gpg2ppmv * ll%o

! Cloud-water fraction

   if (docwmr) then
   call myout ( 6, myname, 'cloud-water fraction' )
   dum = gg%cwmr
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%cwmr, &
                     im*jm,lons,lats,SFL,ORD,FLG,undef)
   else
     ll%cwmr = dum
   end if
   end if

! U-Wind

   call myout ( 6, myname, 'u-wind' )
   dum = gg%u 
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%u, &
                     im*jm,lons,lats,VFL,ORD,FLG,undef)
   else
     ll%u = dum
   end if

! V-Wind

   call myout ( 6, myname, 'v-wind' )
   dum = gg%v
   if (glon.ne.im .or. glat.ne.jm) then
     call interp_h ( dum,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,ll%v, &
                     im*jm,lons,lats,VFL,ORD,FLG,undef)
   else
     ll%v = dum
   end if

   deallocate(lats)
   deallocate(lons)
   deallocate(dphi)
   deallocate(dlam)
   deallocate(dum)
    
   end subroutine gg2ll

!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ll2gg: regrid lat-lon to gaussian
!
! !INTERFACE:
!
   subroutine  ll2gg ( fvgg, gg, rlon, rlat, nsig, glon, glat, docwmr, rot, rc )
!
! !INPUT PARAMETERS:
!
   type(gg_vect), intent(in)  :: fvgg       ! gaussian grid vector
                                            ! with sigma vertical levels
   integer, intent(in)        :: rlon, rlat ! regular FV grid dimensions
   integer, intent(in)        :: nsig       ! number of sigma levels
   integer, intent(in)        :: glon, glat ! gaussian grid dimensions
   logical, intent(in)        :: docwmr     ! controls cloud-water cond mix ratio
   real,    intent(in)        :: rot        ! rotation for fv-grid (from lon(1))
                                            !  (Rotation parameter for lam)
!
! !OUTPUT PARAMETERS:
!
   type(gg_vect), intent(inout)  :: gg  ! gaussian grid vector
                                        ! with FV vertical structure
   integer, intent(out)          :: rc  ! return code
!
! !DESCRIPTION: this routine regrids fields on a regular lat-lon
!               grid to fields on a gaussian grid
!               interpolation is done on each level - no vertical
!               interpolation is performd
!
! !BUGS:
!
! !REVISION HISTORY:
!
!  01Sep2002 Cruz     Created
!  15Mar2004 Todling  Added cloud-water fraction capability
!  30Jul2004 Todling  Updated interface call to comput_gaus
!  28Sep2004 Todling  First-in, Last-out for memory alloc()
!
! !TO DO:
!      Remove dum and dumg , use fv (change intent) and fvgg instead
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'll2gg'
   integer :: L, gxg
   real, pointer :: dlam(:) => null()
   real, pointer :: dphi(:) => null()
   real, pointer :: lons(:,:) => null()
   real, pointer :: lats(:,:) => null()
   real, pointer :: dum(:,:,:) => null()
   real, pointer :: dumg(:,:,:) => null()
   integer, parameter :: ONE=1
   integer, parameter :: SFL=1      ! flag for scalar
   integer, parameter :: VFL=-1     ! flag for vector
   integer, parameter :: ORD=3      ! order of interpolation
   real, parameter    :: TLT = 90.0 ! Rotation parameter for phi
   real, parameter    :: PRE = 0.0  ! Rotation parameter lam_0
   logical, parameter :: FLG = .true.

! Start

   write (*,*) 'will regrid from ',rlon,'x',rlat,' to ',glon,'x',glat
   if (rot .ne. 0.0) write(*,*) 'will rotate grid ',rot,' degrees'
   rc=0
   gxg = glon*glat
 
   call myalloc(dum,rlon,rlat,nsig) 
   call myalloc(dumg,glon,glat,nsig) 
   call myalloc(dlam,rlon)
   call myalloc(dphi,rlat)
   call myalloc(lons,glon,glat)
   call myalloc(lats,glon,glat)

   call comput_gaus(rlon,rlat,glon,glat,dphi,dlam,lons,lats,.false.)

! regrid from latlon grid to a gaussian grid

! surface geopotential

   call myout ( 6, myname,'surface geopotential' )
   dum(:,:,1) = fvgg%hs
   call interp_h ( dum,rlon,rlat,ONE,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,.false.,undef)
   gg%hs = dumg(:,:,1)/grav

! surface pressure

   call myout ( 6, myname,'surface pressure' )
   dum(:,:,1) = fvgg%ps
   call interp_h ( dum,rlon,rlat,ONE,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,.false.,undef)
   gg%ps = dumg(:,:,1)

! delta P

   call myout ( 6, myname,'delta pressure' )
   dum = fvgg%dp
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gg%dp = dumg  

! Temperature

   call myout ( 6, myname,'temperature' )
   dum = fvgg%t
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gg%t = dumg 

! Specific Humidity

   call myout ( 6, myname,'specific humidity' )
   dum = fvgg%q
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gg%q = dumg

! Ozone

   call myout ( 6, myname,'ozone' )
   dum = fvgg%o
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gg%o = dumg
   ! ozone from eta file is in ppmv. we need to convert to kg/kg
   gg%o = o3_ppmv2gpg * gg%o

! Cloud-water fraction

   if (docwmr) then
   call myout ( 6, myname,'cloud-water fraction' )
   dum = fvgg%cwmr
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gg%cwmr = dumg
   endif

! U-Wind

   call myout ( 6, myname,'u-wind' )
   dum = fvgg%u
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,VFL,ORD,FLG,undef)
   gg%u = dumg

! V-Wind

   call myout ( 6, myname,'v-wind' )
   dum = fvgg%v
   call interp_h ( dum,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,VFL,ORD,FLG,undef)
   gg%v = dumg


   deallocate(lats)
   deallocate(lons)
   deallocate(dphi)
   deallocate(dlam)
   deallocate(dumg)
   deallocate(dum)
   
   end subroutine ll2gg

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gg2ss - transform GG to SS vectors
!
!
! !INTERFACE:
!
   subroutine gg2ss ( gg, ss, docwmr, rc )
   implicit none
!
! !INPUT PARAMETERS:
!
   type(gg_vect),intent(in)    :: gg     ! gaussian grid vector
   logical, intent(in)         :: docwmr ! controls cloud-water cond mix ratio
!
! !OUTPUT PARAMETERS:
!
   type(ss_vect),intent(inout) :: ss ! specdtral space vector
   integer, intent(out)        :: rc ! return code
!
! !DESCRIPTION: this routine converts a gaussian grid vector to a 
!               spectral space vector
!
! !BUGS:
!
! !REVISION HISTORY:
!
!  01Aug2002 Cruz     Created
!  02Aug2003 Todling  Initialized rc
!  18Aug2003 Todling  Turned ss into intent(inout) for Halem
!  15Mar2004 Todling  Added reference to cloud-water fraction
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'gg2ss'
   integer, parameter :: IROMB = 0 ! integer spectral domain shape
                                   ! (0 for triangular, 1 for rhomboidal)
   integer, parameter :: IGRID = 4 ! integer grid identifier
                                   ! spectral to grid conversion:
                                   ! =4 for gaussian grid
                                   ! =0 for equally-spaced grid including poles,
                                   ! =256 for equally-spaced grid excluding poles
   integer, parameter :: IDIR = -1 ! integer transform flag
                                   ! idir>0 for wave to grid, idir<0 for grid to wave

   integer :: k,levs,glon,glat,jcap
   real, pointer :: du(:,:,:) => null()
   real, pointer :: dv(:,:,:) => null()
   real, pointer :: virt_t(:,:,:) => null()
   real, pointer :: dumg(:,:,:) => null()
   real, pointer :: logp(:,:) => null()
   real :: con_fvirt

! start

   rc = 0
   jcap=ss%meta%jcap
   levs=ss%meta%nsig
   glon=ss%meta%lonf
   glat=ss%meta%latf
   call myalloc(du,glon,glat,levs)
   call myalloc(dv,glon,glat,levs)
   call myalloc(dumg,glon,glat,levs)
   call myalloc(virt_t,glon,glat,levs)
   call myalloc(logp,glon,glat)
   con_fvirt = 4.6150e+2/2.8705e+2 - 1

! swap north-south pole here and swap vertical levels in ss_put

! topography

   call myout ( 6, myname, 'topography' )
   dumg(:,1:glat,1) = gg%hs(:,glat:1:-1)
   call sptez(IROMB,jcap,IGRID,glon,glat,ss%hs,dumg(:,:,1),IDIR)

! pressure

   call myout ( 6, myname, 'pressure' )
   logp=log(gg%ps/1.e3)
   dumg(:,1:glat,1) = logp(:,glat:1:-1) 
   call sptez(IROMB,jcap,IGRID,glon,glat,ss%ps,dumg(:,:,1),IDIR)

! convert dry bulb T to virtual T

   virt_t=gg%t*(1+con_fvirt*gg%q)
   call myout ( 6, myname, 'temperature' )
   dumg(:,1:glat,:) = virt_t(:,glat:1:-1,:)
   call sptezm(IROMB,jcap,IGRID,glon,glat,levs,ss%t,dumg,IDIR)

! mixing ratio

   call myout ( 6, myname, 'mixing ratio' )
   dumg(:,1:glat,:) = gg%q(:,glat:1:-1,:)
   call sptezm(IROMB,jcap,IGRID,glon,glat,levs,ss%q(:,:,1),dumg,IDIR)

! ozone

   call myout ( 6, myname, 'ozone' )
   dumg(:,1:glat,:) = gg%o(:,glat:1:-1,:)
   call sptezm(IROMB,jcap,IGRID,glon,glat,levs,ss%q(:,:,2),dumg,IDIR)

! cloud-water fraction

   if ( docwmr ) then
   call myout ( 6, myname, 'cloud-water fraction' )
   dumg(:,1:glat,:) = gg%cwmr(:,glat:1:-1,:)
   call sptezm(IROMB,jcap,IGRID,glon,glat,levs,ss%q(:,:,3),dumg,IDIR)
   endif

! u,v winds

   call myout ( 6, myname, 'u,v winds' )
   du(:,1:glat,:) = gg%u(:,glat:1:-1,:)
   dv(:,1:glat,:) = gg%v(:,glat:1:-1,:)
   call sptezmv(IROMB,ss%meta%jcap,IGRID,glon,glat,levs,ss%d,ss%z,du,dv,IDIR)
   
! clean up

   deallocate(du,dv,dumg)
   deallocate(virt_t,logp)

   end subroutine gg2ss

   end module m_xform
   
