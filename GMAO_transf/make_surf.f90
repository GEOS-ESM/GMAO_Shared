!-------------------------------------------------------------------------
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, DAS  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: make_surf --- Module to handle surface fields of SSI/GSI
!
! !INTERFACE:
!
   module make_surf
!
!USES:
!
   use m_gg, only : gg_vect
   use m_ss, only : ss_vect
   use util, only : myexit
   use util, only : myout

   use sfcio_module, only : sfcio_head
   use sfcio_module, only : sfcio_data
   use sfcio_module, only : sfcio_srohdc
   use sfcio_module, only : sfcio_swohdc
   use sfcio_module, only : sfcio_axdata
   use sfcio_module, only : sfcio_realfill

   use m_stdio, only: stdout

   implicit none

   PRIVATE

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC fv2ncep_surf
   PUBLIC ncep2fv_ts
!  PUBLIC ncep_rwsurf

   interface fv2ncep_surf
      module procedure fv2ncep_surf_
   end interface
   interface ncep2fv_ts
      module procedure ncep2fv_surf_
   end interface
   interface ncep_rwsurf
      module procedure ncep_rwsurf_
   end interface

!
! !DESCRIPTION: This module handles surface fields form the NCEP surface files.
!
! !REVISION HISTORY:
!
!   13Dec2004 Todling  Transformed into f90-module; added read-surf 
!                      capability
!   18Feb2005 Todling  Revampped I/O to NCEP SFC file by using sfcio lib
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'make_surf'

   CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fv2ncep_surf_ - create NCEP surface background file
!
!
! !INTERFACE:
!
   subroutine fv2ncep_surf_ ( gg, ss, nymd, nhms, fvsurf, ncsurf, rc )

   implicit none
!
! !INPUT PARAMETERS:
!
   type(gg_vect), intent(in) :: gg    ! gaussian state vector
   type(ss_vect), intent(in) :: ss    ! spectral state vector
   integer, intent(in)       :: nymd  ! year-month-day, e.g., 19990701
   integer, intent(in)       :: nhms  ! hour-min-sec,   e.g., 120000

   character(len = *) :: fvsurf       !   FV-surface filename
   character(len = *) :: ncsurf       ! NCEP-surface filename
!
! !OUTPUT PARAMETERS:
!
   integer, intent(out)        :: rc  ! return code
!
! !DESCRIPTION: this routine reads/regrids fv surface data to generate a
!               surface file for the SSI
!
! !REMARKS:  ncsurf: NCEP's surface file name is just sfcanl, expected by old SSI
!                   (note 6/18/03: latest version expects sfcf06)
!
! !REVISION HISTORY:
!
!  31Dec2002 Cruz     Created
!  14Nov2003 Todling  Changed compilation message.
!  18Dec2003 Cruz     Call to dtoa(v) had bad last argument in call
!  30Jul2004 Todling  Flipping lats from comput_gaus.
!  22Oct2004 Todling  Added NCEP surf filename (asyn fix).
!  21Jul2005 Todling  Added auto-recognition of GEOS-5 files
!  09Aug2005 Todling  Apply dtoa only to fvgcm fields (are they really on d-grid?)
!  13Oct2009 Sienkiewicz few changes for new GSI - set snow to 0 where undef
!                        add z0m, redefine fhour based on ss-file time,
!                        set orography=undef-forces GSI to use uprlev orog.
!
!EOP
!-------------------------------------------------------------------------

! local variables

   character(len=*), parameter         :: myname_ = myname//'fv2ncep_surf_'
   integer :: im,jm,km,lm,nvars,ngatt,timinc,idel
   integer :: glon, glat, glat2, glath
   integer :: iid, err
   integer :: i,j,k,l,n, nsig
   integer :: lxl, nhmss, nymds
   integer, allocatable :: yymmdd(:)
   integer, allocatable :: hhmmss(:)
   integer, allocatable ::  kmvar(:)
   integer :: idate(4)
   integer, parameter :: nrec=20
   character(len=255) :: cvt(1)
   logical fliplon, lz0m

   real :: undef, undef4, undef_dummy, rnd, hourg
   real, allocatable, dimension(:,:) :: u10m, u10mg, &
                                        v10m, v10mg, &
                                        tskin, tsoil, &
                                        sli, &
                                        sno, &
                                        soil_wet, &
                                        f10mg, alb, z0m
   real, allocatable, dimension(:,:,:) :: albedo
   real, allocatable, dimension(:,:) :: vrange, prange, lons, lats
   real, allocatable, dimension(:)   :: dlam, dphi, lon, lat, lev
   real(4) :: fhour4
   real(4), allocatable, dimension (:,:)   :: fld4, ifld4
   real(4), allocatable, dimension (:,:,:) :: fld4_2, ifld4_2
   real(4), allocatable, dimension (:,:,:) :: fld4_4, ifld4_4
   real(8), allocatable, dimension (:,:)   :: fld8

   type(sfcio_head) :: sfchead
   type(sfcio_data) :: sfcdata

   character(len = 255) :: msg

   character(len = 256) :: title
   character(len = 256) :: source
   character(len = 256) :: contact
   character(len = 256) :: levunits
   character(len = 256), allocatable :: vname(:)
   character(len = 256), allocatable :: vtitle(:)
   character(len = 256), allocatable :: vunits(:)

! start

   glon = ss%meta%lonf
   glat = ss%meta%latf
   nsig = ss%meta%nsig
   glat2 = glat-2   ! exclude poles
   glath = glat/2   
   undef = 1.e25    ! fvGCM undef instead of GEOS-3 m_const/getcon ('UNDEF')
   undef4= -9.99e33 ! SSI undef
   rc = 0

! Open fv surface file

   call myout ( 6, myname_, 'Open FV surface background' )
   call Gfio_Open ( trim(fvsurf), 1, iid, err )
   if ( err .LT. 0 ) then
     rc = 1
     call myexit ( myname_, 'could not open surface file', rc )
   endif

   call GFIO_DimInquire ( iid,im,jm,km,lm,nvars,ngatt,rc)

   !write(*,*)myname_,' : Note:'
   !write(*,*)myname_,' : This routine will create a surface background for the SSI'
   !write(*,*)myname_,' : by interpolation the GMAO surface fields from'
   !write(*,*)myname_,' : ',im,'x',jm,' to ',glon,'x',glat2
   !write(*,*)myname_,' : withe the exception of veg frac and veg and soil type'
   !write(*,*)myname_,' : which are copied from an NCEP surface file'

   allocate (dlam(im))
   allocate (dphi(jm))
   allocate (u10m(im,jm))
   allocate (v10m(im,jm))
   allocate (tskin(im,jm))
   allocate (tsoil(im,jm))
   allocate (sno(im,jm))
   allocate (sli(im,jm))
   allocate (soil_wet(im,jm))
   allocate (alb(im,jm))
   allocate (z0m(im,jm))
   allocate (albedo(im,jm,4))

   allocate (lons(glon,glat2))
   allocate (lats(glon,glat2))
   allocate (u10mg(glon,glat2))
   allocate (v10mg(glon,glat2))
   allocate (f10mg(glon,glat2))

   allocate ( lon(im) )
   allocate ( lat(jm) )
   allocate ( lev(km) )
   allocate ( yymmdd(lm) )
   allocate ( hhmmss(lm) )
   allocate ( vname(nvars) )
   allocate ( vtitle(nvars) )
   allocate ( vunits(nvars) )
   allocate (  kmvar(nvars) )
   allocate ( vrange(2,nvars) )
   allocate ( prange(2,nvars) )

   call gfio_inquire ( iid,im,jm,km,lm,nvars, &
                       title,source,contact,undef_dummy, &
                       lon,lat,lev,levunits, &
                       yymmdd,hhmmss,timinc, &
                       vname,vtitle,vunits,kmvar, &
                       vrange,prange,rc )

   if ( abs(lon(1)+180.0)<1.e-5 ) then
        write(6,'(a)') 'Fliping Longitudes of Surface File'
        fliplon = .true.            
        call hflip2_ ( lon, im, 1 )
   else if ( abs(lon(1))<1.e-5 ) then
        fliplon = .false.            
   else
        call Gfio_Close (iid, err)
        rc = 3
        call myexit ( myname_, 'longitudes not compatible w/ GEOS-5/horiz flip', rc )
   endif

   !write(*,*)myname_,' : dates available :'
   !do k=1,lm
   !   write(*,*)myname_,' : ',yymmdd(k),hhmmss(k)
   !end do
   !write(*,*)myname_,' : date requested :',nymd,nhms

! Define variable names.
!
! "U10M"    ! 10 meter u wind (m/s)
! "V10M"    ! 10 meter V wind (m/s)
! "TSKIN"   ! Surface skin temp (K)
! "SNOWDP"  ! Snow depth (2D)
! "ORO"     ! Surface type flag (flag)
! "GWETTOP" ! Top Soil Layer Wetness (fraction)
! "TSOIL1"  ! Top layer soil temperature (K)
! "ALDIF"   ! Albedo: longwave, diffuse
! "ALDIR"   ! Albedo: longwave, direct
! "ASDIF"   ! Albedo: shortwave, diffuse
! "ASDIR"   ! Albedo: shortwave, direct
! "Z0M"     ! surface roughness
! ** Read veg fraction, soil type and veg type from FIX file

   nymds = nymd ; nhmss = nhms 

!  call myout ( 6, myname_, '---get ALDIF' )
!  call Gfio_GetVar ( iid, "ALDIF", nymds, nhmss, im, jm, 1, &
!                   1, alb, rc )
!  if (rc .LT. 0) then
!    call myexit ( myname_, 'error reading ALDIF', rc )
!  else
!     call Undef_fv2ssi ( .true., im*jm, alb, undef, undef4, 'ALDIF' )
!     albedo(1:im,1:jm,1:1) = reshape(alb(1:im,1:jm),(/im,jm,1/))
!  endif
!
!  call myout ( 6, myname_, '---get ALDIR' )
!  call Gfio_GetVar ( iid, "ALDIR", nymds, nhmss, im, jm, 1, &
!                   1, alb, rc )
!  if (rc .LT. 0) then
!    call myexit ( myname_, 'error reading ALDIR', rc )
!  else
!     call Undef_fv2ssi ( .true., im*jm, alb, undef, undef4, 'ALDIR' )
!     albedo(1:im,1:jm,2:2) = reshape(alb(1:im,1:jm),(/im,jm,1/))
!  endif
!
!  call myout ( 6, myname_, '---get ASDIF' )
!  call Gfio_GetVar ( iid, "ASDIF", nymds, nhmss, im, jm, 1, &
!                   1, alb, rc )
!  if (rc .LT. 0) then
!    call myexit ( myname_, 'error reading ASDIF', rc )
!  else
!     call Undef_fv2ssi ( .true., im*jm, alb, undef, undef4, 'ASDIF' )
!     albedo(1:im,1:jm,3:3) = reshape(alb(1:im,1:jm),(/im,jm,1/))
!  endif
!
!  call myout ( 6, myname_, '---get ASDIR' )
!  call Gfio_GetVar ( iid, "ASDIR", nymds, nhmss, im, jm, 1, &
!                   1, alb, rc )
!  if (rc .LT. 0) then
!    call myexit ( myname_, 'error reading ASDIR', rc )
!  else
!     call Undef_fv2ssi ( .true., im*jm, alb, undef, undef4, 'ASDIR' )
!     albedo(1:im,1:jm,4:4) = reshape(alb(1:im,1:jm),(/im,jm,1/))
!  endif

   call myout ( 6, myname_, '---get U10M' )
   call Gfio_GetVar ( iid, "U10M", nymds, nhmss, im, jm, 1, &
                    1, u10m, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading U10M', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, u10m, undef, undef4, 'U10M' )
      if ( fliplon ) call hflip2_ ( u10m,im,jm )
   endif
   
   call myout ( 6, myname_, '---get V10M' )
   call Gfio_GetVar ( iid, "V10M", nymds, nhmss, im, jm, 1, &
                     1, v10m, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading V10M', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, v10m, undef, undef4, 'V10M' )
      if ( fliplon ) call hflip2_ ( v10m,im,jm )
   endif

   call myout ( 6, myname_, '---get TSKIN' )
   call Gfio_GetVar ( iid, "TSKIN", nymds, nhmss, im, jm, 1, &
                     1, tskin, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading TSKIN', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, tskin, undef, undef4, 'TSKIN' )
      if ( fliplon ) call hflip2_ ( tskin,im,jm )
   endif
   
   call myout ( 6, myname_, '---get ORO' )
   call Gfio_GetVar ( iid, "ORO", nymds, nhmss, im, jm, 1,&
                     1, sli, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading ORO', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, sli, undef, undef4, 'ORO' )
      if ( fliplon ) call hflip2_ ( sli,im,jm )
   endif

   call myout ( 6, myname_, '---get SNOWDP' )
   call Gfio_GetVar ( iid, "SNOWDP", nymds, nhmss, im, jm, 1,&
                     1, sno, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading SNOWDP', rc )
   else
      where(sno==undef_dummy) 
         sno=0.
      elsewhere(sno < 0.0)
         sno=0.
      endwhere
!      call Undef_fv2ssi ( .true., im*jm, sno, undef, undef4, 'SNOWDP' )
      if ( fliplon ) call hflip2_ ( sno,im,jm )
   endif

   call myout ( 6, myname_, '---get GWETTOP' )
   call Gfio_GetVar ( iid, "GWETTOP", nymds, nhmss, im, jm, 1,&
                     1, soil_wet, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading GWETTOP', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, soil_wet, undef, undef4, 'GWETTOP' )
      if ( fliplon ) call hflip2_ ( soil_wet,im,jm )
   endif

   call myout ( 6, myname_, '---get TSOIL1' )
   call Gfio_GetVar ( iid, "TSOIL1", nymds, nhmss, im, jm, 1,&
                     1, tsoil, rc )
   if (rc .LT. 0) then
     call myexit ( myname_, 'error reading TSOIL1', rc )
   else
      call Undef_fv2ssi ( .true., im*jm, tsoil, undef, undef4, 'TSOIL1' )
      if ( fliplon ) call hflip2_ ( tsoil,im,jm )
   endif

   call myout ( 6, myname_, '---get Z0M' )
   call Gfio_GetVar ( iid, "Z0M", nymds, nhmss, im, jm, 1,&
                     1, z0m, rc )
   if (rc .LT. 0) then
     call myout (6, myname_, 'error reading Z0M - take field from input sfc')
     lz0m = .false.
   else
      call Undef_fv2ssi ( .true., im*jm, tsoil, undef, undef4, 'Z0M' )
      if ( fliplon ) call hflip2_ ( z0m,im,jm )
      lz0m = .true.
   endif

   call Gfio_Close (iid, rc)

   allocate(fld8(glon,glat2))

!  calculate gaussian grid 

   call comput_gaus(im,jm,glon,glat2,dphi,dlam,lons,lats,.true.)
   lxl=glon*glat2

!  calculate 10m wind factors and regrid 

   call myout (6, myname_, 'Interpolate fields, perform corrections and create sfcf06') 
   call myout ( 6, myname_, '---regrid u,v and calculate 10m wind factors' )

! u-wind

   ! put u10m winds on 'A' grid before calculating f10m ( gg is on A grid )
   if ( .not. fliplon ) call dtoa ( u10m,u10m,im,jm,1,2 )
   call interp_h (u10m,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
   lxl,lons,lats,1,3,.true.,undef4)
   u10mg = fld8
   
! v-wind

   ! put v10m winds on 'A' grid before calculating f10m ( gg is on A grid )
   if ( .not. fliplon ) call dtoa ( v10m,v10m,im,jm,1,1 )
   call interp_h (v10m,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
   lxl,lons,lats,1,3,.true.,undef4)
   v10mg = fld8

! compute 10m factors 

   f10mg = sqrt ( u10mg**2 + v10mg**2 ) / &
           sqrt ( gg%u(:,2:glat-1,nsig)**2 + gg%v(:,2:glat-1,nsig)**2)
   where (f10mg>1.0)
      f10mg = 1.0
   end where

   !call myout ( 6, myname_, 'regrid other fields and write out surface file' )

! write out fields to SSI surface file in correct order

   ! compute date - now use ss%meta%fhour to match ss file
   idel = -3600 * ss%meta%fhour
   call tick (nymds, nhmss, idel)
   fhour4 = ss%meta%fhour
   idate(1) = nhmss/10000
   idate(3) = nymds-100  *(nymds/100)
   idate(2) = nymds-10000*(nymds/10000)-idate(3)
   idate(2) = idate(2)/100
   idate(4) = nymds/10000

! Read in NCEP SFC file (file number wired-in still)
! ---------------------
   call sfcio_srohdc(34,'ncepsfc',sfchead,sfcdata,rc)
     if(rc/=0)then
       call myexit ( myname_, 'error reading NCEP SFC file', rc )
     endif

   print *, 'fv2ncep_surf: original NCEP surf file with lon= ', sfchead%lonb, ' and lat= ', sfchead%latb
   print *, 'fv2ncep_surf: writing  GMAO surf file with lon= ', glon, ' and lat2= ', glat2

   !call myout ( 6, myname_, 'open files')

!  Start redefining surface fields w/ a few GMAO provided fields
!  -------------------------------------------------------------
   sfchead%fhour = fhour4
   sfchead%idate = idate

!    surface fields ordering:   SSI         GSI        
! 3  surface T                 1layer      1layer       same
! 4  soil moisture             1layer      2layers      same
! 5  snow depth                1layer      1layer       same
! 6  soil T                    1layer      2layers      same
! 7  not used                  1layer      1layer       tg3
! 8  roughness (cm)            1layer      1layer       zor
! 9  not used                  1layer      1layer       cv
! 10 not used                  1layer      1layer       cvb
! 11 not used                  1layer      1layer       cvt
! 12 not used                  1layer      4layers      albedo
! 13 SLI flag                  1layer      1layer       same
! 14 veg fraction (from NCEP)  1layer      1layer       same
! 15 not used                  1layer      1layer       plantr
! 16 10 m wind fraction        1layer      1layer       same
! 17 not in SSI?! (from NCEP)   ---        1layer       canopy water cont.
! 18 veg type (from NCEP)      1layer      1layer       same
! 19 soil type(from NCEP)      1layer      1layer       same
! 20    not in SSI?!            ---        2layers      zenith angle veg. frac
!
! NOTES: 
!   from NCEP means read from a fixed ncepsfc file for now
!   New format ivs>=200501 includes orog (orography) field used by 2009 GSI
!   (GSI will interpolate terrain from atm file if not avail in sfc file)
!   2009 GSI reads: tsea sheleg zorl slmsk f10m orog  & for radassim: smc stc vfrac vtype stype  
!

   allocate(fld4(im,jm),fld4_2(im,jm,2),fld4_4(im,jm,4))   ! _RT: This is a waste of memory 
   allocate(ifld4(glon,glat2),ifld4_2(glon,glat2,2),ifld4_4(glon,glat2,4))
   if(im/=glon) &
      print *, ' Case when im=', im, ' != glon= ', glon, &
               ' interpolation will take place'

! surface T

   call myout ( 6, myname_, '---surface T')
   fld4 = tskin 
   if(im/=glon) then 
      call interp_h (tskin,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      ifld4 = fld8
   else
      ifld4 = fld4
   end if
   sfcdata%tsea = ifld4

! soil moisture

   call myout ( 6, myname_, '---soil moisture')
   fld4 = soil_wet 
   if(im/=glon) then 
      call interp_h (soil_wet,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      ifld4 = fld8
      where (ifld4<=0.0)
        ifld4 = 0.05
      end where
      where (ifld4>1.0)
        ifld4 = 1.0
      end where
      ifld4_2(1:glon,1:glat2,1:1) = reshape(ifld4,(/glon,glat2,1/))
      ifld4_2(1:glon,1:glat2,2:2) = reshape(ifld4,(/glon,glat2,1/))
   else
      ifld4_2(1:im,1:jm,1:1) = reshape(fld4,(/im,jm,1/))
      ifld4_2(1:im,1:jm,2:2) = reshape(fld4,(/im,jm,1/))
   end if
   sfcdata%smc = ifld4_2   ! need to be careful w/ dim here

! snow depth  <---- need to check units

   call myout ( 6, myname_, '---snow depth')
   fld4 = sno 
   if(im/=glon) then 
      call interp_h (sno,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      ifld4 = fld8
      where (ifld4<0.0)
        ifld4 = 0.0
      end where
   else
      ifld4 = fld4
   end if
   sfcdata%sheleg = ifld4

! soil T 

   call myout ( 6, myname_, '---soil T')
   fld4 = tsoil
   if(im/=glon) then 
      call interp_h (tsoil,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      ifld4 = fld8
      ifld4_2(1:glon,1:glat2,1:1) = reshape(ifld4,(/glon,glat2,1/))
      ifld4_2(1:glon,1:glat2,2:2) = reshape(ifld4,(/glon,glat2,1/))
   else
      ifld4_2(1:im,1:jm,1:1) = reshape(fld4,(/im,jm,1/))
      ifld4_2(1:im,1:jm,2:2) = reshape(fld4,(/im,jm,1/))
   end if
   sfcdata%stc = ifld4_2   ! need to be careful w/ dim

   if (lz0m) then
! surface roughness - convert from m  

   call myout ( 6, myname_, '---z0m')
   fld4 = z0m
   if(im/=glon) then 
      call interp_h (z0m,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      ifld4 = fld8
   else
      ifld4 = fld4
   end if
   sfcdata%zorl = ifld4*100.
   end if

!   SLI flag <---- NEED TO CHECK INTERPOLATION

   ! special case to interpolate SLI flag 

   !  S = 0, L = 1, I = 2
   call myout ( 6, myname_, '---SLI flag')
   fld4 = sli 
   if(im/=glon) then 
      call interp_h (sli,im,jm,1,dlam,dphi,0.0,90.0,0.0,fld8, &
      lxl,lons,lats,1,3,.true.,undef4)
      !call interpint (sli,im,jm,dlam,dphi,fld8,lxl,lons,lats)
      ifld4 = fld8
      where (ifld4<0.5)
        ifld4 = 0.0
      end where
      where (ifld4>=0.5 .and. ifld4<1.5)
        ifld4 = 1.0
      end where
      where (ifld4>=1.5)
        ifld4 = 2.0
      end where
   else
      ifld4 = fld4
   end if
   sfcdata%slmsk = ifld4

   call myout ( 6, myname_, '---10m wind frac')
   ifld4 = f10mg
   sfcdata%f10m = ifld4
   
!  fill in orography with sigio missing value 
!  ------------------------------------------

   sfcdata%orog = sfcio_realfill

!  Write out surface file with GMAO entries
!  ----------------------------------------
   call sfcio_swohdc(33,trim(ncsurf),sfchead,sfcdata,rc)
     if(rc/=0)then
       call myexit ( myname_, 'error writing NCEP SFC file', rc )
     endif
   call sfcio_axdata(sfcdata,rc)  ! release memory
     if (rc/=0) then
          print*, trim(myname_), ': cannot release mem for NCEP surface field, ier= ', rc
          return
     endif

! clean up

   deallocate ( fld4, fld8 )
   deallocate ( dlam, dphi, lons, lats, lon, lat, lev )
   deallocate ( u10m, v10m, tskin, tsoil, sno, sli, soil_wet )
   deallocate ( f10mg, u10mg, v10mg )
   deallocate ( alb, albedo )
   deallocate ( yymmdd, hhmmss )
   deallocate ( vname, vtitle, vunits, kmvar, vrange, prange )


   msg = 'Finished creating file: ' // trim(ncsurf)
   call myout ( 6, myname_, trim(msg))

   end subroutine fv2ncep_surf_

!-------------------------------------------------------------------------
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, DAS  !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ncep2fv_surf_ - update FV surface background field(s)
!
! !INTERFACE:

   subroutine ncep2fv_surf_ ( im, jm, glon, glat2, Tskin_Ga, ncsurf, rc, &
                              Tskin_Ge )

   implicit none

! !INPUT PARAMETERS:

   integer, intent(in) :: im, jm
   integer, intent(in) :: glon, glat2
   character(len=*), intent(in) :: ncsurf

! !INPUT/OUTPUT PARAMETERS:

   real, intent(inout)        :: Tskin_Ga(im,jm) ! GMAO analysis/background

! !OUTPUT PARAMETERS:

   real, intent(in), optional :: Tskin_Ge(im,jm) ! GMAO back/forth converted background
   integer, intent(out) :: rc

! !DESCRIPTION: Updates surface fields with SSI/GSI surface analysis.
!  
! !REMARKS: For now, operates only on Tskin.
!
! !REVISION HISTORY:
!
!   13Dec2004 Todling  Initial code.
!   12Jan2005 Todling  Undef defined per fv-gcm.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname_ = myname//'ncep2fv_surf_'

   character(len=*), parameter :: ndnhfn = 'eta_out.nc4' ! wired-in non-data no-harm filename

   real, allocatable :: Tskin_a (:,:) ! NCEP analysis on gaussian-grid
   real, allocatable :: Tskin_Na(:,:) ! NCEP analysis on fv-grid
   real, allocatable :: dlam(:)
   real, allocatable :: dphi(:)
   real, allocatable :: lons(:)
   real, allocatable :: lats(:)

   type(sfcio_head) :: sfchead
   type(sfcio_data) :: sfcdata

   logical :: verbose = .false.
   logical :: noharm, iexist
   integer :: j, ncid, ier
   real    :: undef

   noharm = .false. ! default 
   rc     = 0
   undef = 1.e25    ! fvGCM undef instead of GEOS-3 m_const/getcon ('UNDEF')

!  Allocate memory for local arrays
!  --------------------------------
   allocate( Tskin_a(glon,glat2), Tskin_Na(im,jm), stat=ier )
     if ( ier/=0) then
          print*, trim(myname_), ': Error in Alloc(Tskin)'
          rc = 1
          return
     endif

   allocate( dlam(glon), dphi(glat2), lons(im*jm), lats(im*jm), stat=ier )
     if ( ier/=0 ) then
          print*, trim(myname_), ': Error in Alloc(aux)'
          rc = 1
          return
     endif

!  Read in Tksin for NCEP surface analysis file (third record, jrec=3)
!  ------------------------------------------------------------------
   call sfcio_srohdc(34,trim(ncsurf),sfchead,sfcdata,ier)
   Tskin_a = sfcdata%tsea
!  call ncep_rwsurf_ ( verbose, trim(ncsurf), glat2, glon, ier, jrec=3, fld=Tskin_a )
     if ( ier/=0 ) then
          print*, trim(myname_), ': cannot read NCEP surface analysis file, ier= ', ier
          rc = 2
          return
     endif
   call sfcio_axdata(sfcdata,ier)  ! release memory
     if ( ier/=0 ) then
          print*, trim(myname_), ': cannot release mem for NCEP surface field, ier= ', ier
          rc = 2
          return
     endif

!  Interpolate Tskin to proper fv-grid
!  -----------------------------------
   call comput_gaus(glon,glat2,im,jm,dphi,dlam,lons,lats,.true.)
   if(im/=glon) then
     call interp_h ( Tskin_a,glon,glat2,1,dlam,dphi,0.0,90.0,0.0, Tskin_Na, &
                     im*jm,lons,lats,1,3,.false.,undef)
   else
      do j = 1, jm
         Tskin_Na(1:im,jm-j+1) = Tskin_a(1:im,j)  ! swap latitudes
      end do
   endif

   deallocate( dlam, dphi, lons, lats, stat=ier )
     if ( ier/=0) then
          print*, trim(myname_), ': Error in Dealloc(aux)'
          rc = 99
          return
     endif

!  Check on possibly existing guassian converted file
!  --------------------------------------------------
   if (present(Tskin_Ge) ) then
       print *, trim(myname_), ': handling no-data no-harm Tskin analysis case ...'

!      Update Tskin
!      ------------
       Tskin_Ga = Tskin_Ga + ( Tskin_Na - Tskin_Ge )

   else
       print *, trim(myname_), ': updating Tskin field with converted gaussian grid Tskin ...'
       Tskin_Ga = Tskin_Na
   endif

   deallocate( Tskin_a, Tskin_Na, stat=ier )
     if ( ier/=0) then
          print*, trim(myname_), ': Error in Dealloc()'
          rc = 99
          return
     endif

   end subroutine ncep2fv_surf_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
   subroutine interpint (fin,im,jm,dlam,dphi,fout,lxl,lons,lats) 

! HISTORY:
!
!  30Jan2003  Cruz     Initial code.
!
!-------------------------------------------------------------------------
!BOC

   integer lxl
   integer im,jm,i,j
   integer, allocatable :: ip1(:), ip0(:), im1(:), im2(:)
   integer, allocatable :: jp1(:), jp0(:), jm1(:), jm2(:)
   real      fout(lxl)
   real    lons(lxl)
   real    lats(lxl)
   real    fin(im,jm)
   real     dlam(im)
   real     dphi(jm)
   real    lon_tmp(im)
   real    lat_tmp(jm)
   real    f_tmp(lxl)
   real sinnp, cosnp, lam_tmp,phi_tmp, p1,p2,p3,lam_0,lam_np
   integer im1_tmp,jm1_tmp ,itmp,jtmp

   real pi

! start

   allocate (ip1(lxl) ,ip0(lxl) ,im1(lxl) ,im2(lxl) )
   allocate (jp1(lxl) ,jp0(lxl) ,jm1(lxl) ,jm2(lxl) )

   cosnp = 0.0 ; sinnp = 1.0
   pi = 4.*atan(1.)
  
   lon_tmp(1) = -pi
   do i=2,im
     lon_tmp(i) = lon_tmp(i-1) + dlam(i-1)
   enddo
   lat_tmp(1) = -pi*0.5
   do j=2,jm-1
     lat_tmp(j) = lat_tmp(j-1) + dphi(j-1)
   enddo
   lat_tmp(jm) =  pi*0.5
 
! determine indexing based on input grid
 
   lam_0 = 0.0 
   lam_np = 0.0 
   do i=1,lxl
     p1 = cosnp*cos(lats(i))*cos(lons(i)+lam_0-pi) &
        + sin(lats(i))*sinnp
     p1 = min(p1, 1.0)
     p1 = max(p1,-1.0)
     phi_tmp = asin( p1 )

     p2 = sinnp*cos(lons(i)+lam_0-pi)
     p2 = acos( p2 )

     p3 = cos(lats(i))*sin(lons(i)+lam_0-pi)
     if( p3.lt.0.0 ) p2 = -p2
     p2 = p2 + lam_np - pi
     lam_tmp = mod( p2+3.0*pi,2.0*pi ) - pi

! Determine Indexing Based on Computational Grid
! ----------------------------------------------
     im1_tmp = 1
     do itmp = 2,im
       if( lon_tmp(itmp).lt.lam_tmp ) im1_tmp = itmp
     enddo
     jm1_tmp = 1
     do jtmp = 2,jm
       if( lat_tmp(jtmp).lt.phi_tmp ) jm1_tmp = jtmp
     enddo

     im1(i) = im1_tmp
     ip0(i) = im1(i) + 1
     ip1(i) = ip0(i) + 1
     im2(i) = im1(i) - 1

     jm1(i) = jm1_tmp
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
   end do 

   do i=1,lxl
     f_tmp(i) = ( fin( im1(i),jm1(i) ) + fin( ip0(i),jm1(i) ) + &
                fin( im1(i),jp0(i) ) + fin( ip0(i),jp0(i) ) ) / 4.0
   end do

   fout = f_tmp 

   deallocate (ip1,ip0,im1,im2)
   deallocate (jp1,jp0,jm1,jm2)

   end subroutine interpint

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Undef_fv2ssi - Replace undef FV values on input by SSI undef
! 
! !INTERFACE:
!
      subroutine Undef_fv2ssi ( verbose, ndim, var, undef, undef4, vname )
!
! !USES:


      implicit none


! !INPUT PARAMETERS: 
 
      logical  verbose
      integer  ndim
      real     undef
      real     undef4
      character*(*) vname

! !INPUT/OUTPUT PARAMETERS: 
 
      real     var(ndim)

! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  13Apr00   Todling    - Initial code.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'Undef_fv2ssi'

      integer  i
      integer  nonmiss
      real     tol
      real     val_max_in, val_max_out

      tol   = 0.01 * undef

!     Make sure there is no precision problem with undef's
!     ----------------------------------------------------
      nonmiss = 0
      val_max_in  = maxval(var)
      do i = 1, ndim
         if ( abs(var(i)-undef) .lt. tol ) then
              nonmiss = nonmiss + 1
              var(i) = undef4
         end if
      end do
      val_max_out = maxval(abs(var))
      if (nonmiss.ne.0) then
         if(verbose) &
         write(stdout,'(2a,i10,3a)') myname_, ': Fixed ', nonmiss, &
                                   ' values from input with strange undef', &
                                   ' for variable ', trim(vname)
         if ( val_max_out .ne. abs(undef4) ) then
            write(stdout,*) ' Largest  value on  input: ',  val_max_in
            write(stdout,*) ' Largest  value on output: ',  val_max_out
            write(stdout,*) ' Undef        value spec.: ',  undef
            write(stdout,*) ' Correction not done. Aborting ... '
            call exit (7)
         end if
      end if

      return
      end subroutine Undef_fv2ssi

!-------------------------------------------------------------------------
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, DAS  !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ncep_rwsurf_ - create SSI surface background file
!
! !INTERFACE:

   subroutine ncep_rwsurf_ ( verbose, infile, glat2, glon, stat, & 
                             outfile, wout, jrec, fld )

! !USES:

   implicit none

! INPUT PARAMETERS:

   logical :: verbose
   integer :: glat2, glon
   character(len=*), intent(in)           :: infile
   character(len=*), intent(in), optional :: outfile
   logical,          intent(in), optional :: wout
   integer,                      optional :: jrec     ! set to record of fld to be returned

! OUTPUT PARAMETERS:

   real, intent(out), optional :: fld(glon,glat2)  ! requested output field from file
   integer, intent(out) :: stat                    ! return error code

! !DESCRIPTION: This module handles surface fields form the NCEP surface files.
!
! !REVISION HISTORY:
!   
!  13Dec2004 Todling  Initial code (already existed somewhere else).
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname_ = myname // 'ncep_rwsurf_'

   real(4) :: sdum4
   real(8), allocatable :: junk(:,:)
   real(4), allocatable :: vfrac(:,:)
   real(4), allocatable :: vtype(:,:)
   real(4), allocatable :: stype(:,:)
   real(4), allocatable :: fdum4(:,:)
   real(4), allocatable :: fdum42(:,:,:)
   real(4), allocatable :: albedo4(:,:,:)
   real(4), allocatable :: zenith2(:,:,:)

   logical :: writeout
   integer :: icount, irec
   character(8),dimension(4):: labfix
   real(4) yhour
   integer latd,lonl,version
   integer,dimension(4):: igdate,iadate
   integer,allocatable,dimension(:):: lonsperlat

   stat     = 0
   irec     = 0
   writeout = .false.
   if (present(jrec)) then
       if(.not.present(fld))then
          print *, trim(myname_), ' Error, missing fld array'
          stat = 1
          return
       endif
       irec = jrec
   endif

   ! unit 34 is a temporary file that contains surface fields
   ! not available from GMAO and so the corresponding NCEP fields are used

   open(34,file=infile,form='unformatted')
   if(present(wout)) then
      if(wout)then
         writeout = .true.
         if (present(outfile)) then
             open(35,file=outfile,form='unformatted')
         else
             print *, trim(myname_), ' Error, no filename specified to write sfc flds'
             stat = 2
             return
         endif
      endif
   endif
   allocate ( junk   (glon,glat2)   )
   allocate ( vtype  (glon,glat2)   )
   allocate ( stype  (glon,glat2)   )
   allocate ( vfrac  (glon,glat2)   )
   allocate ( fdum4  (glon,glat2)   )
   allocate ( fdum42 (glon,glat2,2) )
   allocate ( albedo4(glon,glat2,4) )
   allocate ( zenith2(glon,glat2,2) )

   read(34) 
   read(34) yhour,IgDATE,LONL,LATD,version
   allocate ( lonsperlat(latd/2) )

   rewind(34)
   read(34) labfix
   read(34) yhour,IgDATE,LONL,LATD,version,lonsperlat
   icount = 2
   if(writeout) write(35) labfix
   if(writeout) write(35) yhour,IgDATE,LONL,LATD,version,lonsperlat
   read(34) fdum4                      ! record 1  tsf
       icount = icount + 1
       if(verbose) print*,icount, ' ',maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum42                     ! record 2  soilm
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum42)
       if(writeout) write(35) fdum42
   read(34) fdum4                      ! record 3  snow
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum42                     ! record 4  soilt
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum42)
       if(writeout) write(35) fdum42
   read(34) fdum4                      ! record 5  tg3
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 6  zor
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 7  cv
       icount = icount + 1
       if(verbose) print*,icount, ' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 8  cvb
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 9  cvt
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) albedo4                    ! record 10 albedo
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(albedo4)
       if(writeout) write(35) albedo4
   read(34) fdum4                      ! record 11 slimsk
       icount = icount + 1
       if(verbose) print*,icount,' ',maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) vfrac                      ! record 12 vegetation cover
       icount = icount + 1
       junk = vfrac
       vfrac= junk
       if(verbose) print*,icount,' ', maxval(vfrac)
       if(writeout) write(35) vfrac
   read(34) fdum4                      ! record 13 canopy water
       icount = icount + 1
       junk = fdum4
       fdum4 = junk
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 14 f10m
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) vtype                      ! record 15 vegetation type
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(vtype)
       if(writeout) write(35) vtype
   read(34) stype                      ! record 16 soil type
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(stype)
       if(writeout) write(35) stype
   read(34) zenith2                    ! record 17 zenith
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(zenith2)
       if(writeout) write(35) zenith2
   read(34) fdum4                      ! record 18 ustar
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 19 ffmm
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4
   read(34) fdum4                      ! record 20 ffhh
       icount = icount + 1
       if(verbose) print*,icount,' ', maxval(fdum4)
       if(writeout) write(35) fdum4
       if(irec==icount) fld = fdum4

   close(34)
   if(writeout)close(35)

   deallocate ( lonsperlat )
   deallocate ( junk  )
   deallocate ( vtype  )
   deallocate ( stype  )
   deallocate ( vfrac  )
   deallocate ( fdum4   )
   deallocate ( fdum42   )
   deallocate ( albedo4 )
   deallocate ( zenith2 )
   print * , trim(myname_), ': done'
   end subroutine ncep_rwsurf_

  subroutine hflip2_ ( q,im,jm )
      implicit none
      integer  im,jm,i,j
      real, intent(inout) :: q(im,jm)
      real, allocatable   :: dum(:)
      allocate ( dum(im) )
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j)
         dum(i+im/2) = q(i,j)
      enddo
         q(:,j) = dum(:)
      enddo
      deallocate ( dum )
  end subroutine hflip2_
   end module make_surf
