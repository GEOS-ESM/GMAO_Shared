      subroutine interp_h ( q_cmp,im,jm,lm,
     .                      dlam,dphi,rotation,tilt,precession,
     .                      q_geo,irun,lon_geo,lat_geo,
     .                      msgn,norder,check,undef )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Performs a horizontal interpolation from a field on a computational grid
C    to arbitrary locations.
C
C  INPUT:
C  ======
C    q_cmp ...... Field q_cmp(im,jm,lm) on the computational grid
C    im ......... Longitudinal dimension of q_cmp
C    jm ......... Latitudinal  dimension of q_cmp
C    lm ......... Vertical     dimension of q_cmp
C    dlam ....... Computational Grid Delta Lambda
C    dphi ....... Computational Grid Delta Phi
C    rotation ... Rotation parameter lam_np (Degrees)
C    tilt ....... Rotation parameter phi_np (Degrees)
C    precession . Rotation parameter lam_0  (Degrees)
C    irun ....... Number of Output Locations
C    lon_geo .... Longitude Location of Output
C    lat_geo .... Latitude  Location of Output
C    msgn ....... Flag for scalar field  ( msgn =  1 )
C                    or vector component ( msgn = -1 )
C    norder ..... Order of Interpolation:  Bi-Linear => abs(norder) = 1
C                                          Bi-Cubic  => abs(norder) = 3
C                 Note: If norder < 0, then check for positive definite
C    check ...... Logical Flag to check for Undefined values
C
C  OUTPUT:
C  =======
C    q_geo ...... Field q_geo(irun,lm) at arbitrary locations
C
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      implicit none

c Input Variables
c ---------------
      integer,intent(in) :: im,jm,lm,irun,norder,msgn
      logical,intent(in) :: check

      real,intent(out)::      q_geo(irun,lm)
      real,intent(in) ::    lon_geo(irun)
      real,intent(in) ::    lat_geo(irun)

      real,intent(in) ::    q_cmp(im,jm,lm)
      real,intent(in) ::     dlam(im)
      real,intent(in) ::     dphi(jm)

c Local Variables
c ---------------
      integer  i,j,l,m,n
      integer, allocatable       :: ip1(:), ip0(:), im1(:), im2(:)
      integer, allocatable       :: jp1(:), jp0(:), jm1(:), jm2(:)

      integer ip1_for_jp1, ip0_for_jp1, im1_for_jp1, im2_for_jp1
      integer ip1_for_jm2, ip0_for_jm2, im1_for_jm2, im2_for_jm2
      integer jm2_for_jm2, jp1_for_jp1


c Bi-Linear Weights
c -----------------
      real, allocatable       ::    wl_ip0jp0 (:)
      real, allocatable       ::    wl_im1jp0 (:)
      real, allocatable       ::    wl_ip0jm1 (:)
      real, allocatable       ::    wl_im1jm1 (:)

c Bi-Cubic Weights
c ----------------
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

      real, allocatable       ::    old_lon (:)
      real, allocatable       ::    old_lat (:)
      real, allocatable       ::    old_dlam(:)
      real, allocatable       ::    old_dphi(:)

      real    ux, ap1, ap0, am1, am2
      real    uy, bp1, bp0, bm1, bm2

      real    lon_cmp(im)
      real    lat_cmp(jm)
      real    q_tmp(irun)

      real    pi,cosnp,sinnp,p1,p2,p3,eps,d
      real    lam,lam_ip1,lam_ip0,lam_im1,lam_im2
      real    phi,phi_jp1,phi_jp0,phi_jm1,phi_jm2
      real    dl,dp,lam_np,phi_np,lam_0,eps_np
      real,intent(in) :: rotation , tilt , precession
      real    lam_geo, lam_cmp
      real    phi_geo, phi_cmp
      real,intent(in) :: undef
      integer im1_cmp,icmp
      integer jm1_cmp,jcmp

      logical compute_weights
      real    old_rotation
      real    old_tilt
      real    old_precession
      data    old_rotation   /-999.9/
      data    old_tilt       /-999.9/
      data    old_precession /-999.9/

      parameter ( eps = 1.e-10 )

c Initialization
c --------------
      pi = 4.*atan(1.)
      dl = 2*pi/ im     ! Uniform Grid Delta Lambda
      dp =   pi/(jm-1)  ! Uniform Grid Delta Phi

c Allocate Memory for Weights and Index Locations
c -----------------------------------------------
      if(.not.allocated(old_lon)) then

      allocate (   old_dlam(im)  ,   old_dphi(jm)  )
      allocate (   old_lon(irun) ,   old_lat(irun) )
      allocate ( wl_ip0jp0(irun) , wl_im1jp0(irun) )
      allocate ( wl_ip0jm1(irun) , wl_im1jm1(irun) )
      allocate ( wc_ip1jp1(irun) , wc_ip0jp1(irun) ,
     &		 wc_im1jp1(irun) , wc_im2jp1(irun) )
      allocate ( wc_ip1jp0(irun) , wc_ip0jp0(irun) ,
     &		 wc_im1jp0(irun) , wc_im2jp0(irun) )
      allocate ( wc_ip1jm1(irun) , wc_ip0jm1(irun) ,
     &		 wc_im1jm1(irun) , wc_im2jm1(irun) )
      allocate ( wc_ip1jm2(irun) , wc_ip0jm2(irun) ,
     &		 wc_im1jm2(irun) , wc_im2jm2(irun) )
      allocate (       ip1(irun) ,       ip0(irun) ,
     &		       im1(irun) ,       im2(irun) )
      allocate (       jp1(irun) ,       jp0(irun) ,
     &		       jm1(irun) ,       jm2(irun) )
      do i=1,irun
      old_lon(i) = -999.9
      old_lat(i) = -999.9
      enddo
      do i=1,im
      old_dlam(i) = 0.0
      enddo
      do j=1,jm
      old_dphi(j) = 0.0
      enddo

      else
             i =  size (old_dlam)
             j =  size (old_dphi)
             m =  size (old_lon)
          if(i.ne.im .or. j.ne.jm .or. m.ne.irun) then
          deallocate (  old_dlam ,  old_dphi )
          deallocate (   old_lon ,   old_lat )
          deallocate ( wl_ip0jp0 , wl_im1jp0 )
          deallocate ( wl_ip0jm1 , wl_im1jm1 )
          deallocate ( wc_ip1jp1 , wc_ip0jp1 , wc_im1jp1 , wc_im2jp1 )
          deallocate ( wc_ip1jp0 , wc_ip0jp0 , wc_im1jp0 , wc_im2jp0 )
          deallocate ( wc_ip1jm1 , wc_ip0jm1 , wc_im1jm1 , wc_im2jm1 )
          deallocate ( wc_ip1jm2 , wc_ip0jm2 , wc_im1jm2 , wc_im2jm2 )
          deallocate (       ip1 ,       ip0 ,       im1 ,       im2 )
          deallocate (       jp1 ,       jp0 ,       jm1 ,       jm2 )
          allocate (   old_dlam(im)  ,   old_dphi(jm)  )
          allocate (   old_lon(irun) ,   old_lat(irun) )
          allocate ( wl_ip0jp0(irun) , wl_im1jp0(irun) )
          allocate ( wl_ip0jm1(irun) , wl_im1jm1(irun) )
          allocate ( wc_ip1jp1(irun) , wc_ip0jp1(irun) ,
     &               wc_im1jp1(irun) , wc_im2jp1(irun) )
          allocate ( wc_ip1jp0(irun) , wc_ip0jp0(irun) ,
     &               wc_im1jp0(irun) , wc_im2jp0(irun) )
          allocate ( wc_ip1jm1(irun) , wc_ip0jm1(irun) ,
     &               wc_im1jm1(irun) , wc_im2jm1(irun) )
          allocate ( wc_ip1jm2(irun) , wc_ip0jm2(irun) ,
     &               wc_im1jm2(irun) , wc_im2jm2(irun) )
          allocate (       ip1(irun) ,       ip0(irun) ,
     &                     im1(irun) ,       im2(irun) )
          allocate (       jp1(irun) ,       jp0(irun) ,
     &                     jm1(irun) ,       jm2(irun) )
          do i=1,irun
          old_lon(i) = -999.9
          old_lat(i) = -999.9
          enddo
          do i=1,im
          old_dlam(i) = 0.0
          enddo
          do j=1,jm
          old_dphi(j) = 0.0
          enddo
          endif
      endif

c Compute Input Computational-Grid Latitude and Longitude Locations
c -----------------------------------------------------------------
      lon_cmp(1) = -pi
      do i=2,im
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_cmp(1) = -pi*0.5
      do j=2,jm-1
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_cmp(jm) =  pi*0.5

c Check for Co-incident Grid-Point Latitude and Pole Locations
c ------------------------------------------------------------
                                          eps_np = 0.0
      do j=1,jm
      phi_cmp = lat_cmp(j)*180./pi
      if( abs( phi_cmp-tilt ).lt.1.e-3 )  eps_np =  1.e-3
      if( tilt+eps_np .gt. 90. )          eps_np = -1.e-3
      enddo

      lam_np = pi/180.*rotation
      phi_np = pi/180.*(tilt+eps_np)
      lam_0  = pi/180.*precession

      if( tilt.eq.90. ) then 
      cosnp = 0.0
      sinnp = 1.0
      else if(tilt.eq.-90.0) then
      cosnp = 0.0
      sinnp =-1.0
      else
      cosnp = cos(phi_np)
      sinnp = sin(phi_np)
      endif

c Determine if Weights Need to be Updated
c ---------------------------------------
      compute_weights =  rotation.ne.old_rotation .or.
     .                       tilt.ne.old_tilt     .or.
     .                 precession.ne.old_precession

      m = 1
      do while ( .not.compute_weights .and. m.le.irun )
      compute_weights = (lon_geo(m).ne.old_lon(m)) .or.
     .                  (lat_geo(m).ne.old_lat(m))
      m = m+1
      enddo
      i = 1
      do while ( .not.compute_weights .and. i.le.im )
      compute_weights = dlam(i).ne.old_dlam(i)
      i = i+1
      enddo
      j = 1
      do while ( .not.compute_weights .and. j.le.jm-1 )
      compute_weights = dphi(j).ne.old_dphi(j)
      j = j+1
      enddo

c Compute Weights for Computational to Geophysical Grid Interpolation
c -------------------------------------------------------------------
      if( compute_weights ) then

      old_rotation   = rotation
      old_tilt       = tilt
      old_precession = precession

#if   (openmp)
!$omp  parallel do 
!$omp& default (shared)
!$omp& private (i,lam_geo,phi_geo,lam_cmp,phi_cmp,lam,phi)
!$omp& private (p1,p2,p3,d,icmp,jcmp,im1_cmp,jm1_cmp)
!$omp& private (lam_im2, lam_im1, lam_ip0, lam_ip1)
!$omp& private (phi_jm2, phi_jm1, phi_jp0, phi_jp1)
!$omp& private (ap1, ap0, am1, am2)
!$omp& private (bp1, bp0, bm1, bm2)
#endif
      do i=1,irun
      old_lon(i) = lon_geo(i)
      old_lat(i) = lat_geo(i)
      lam_geo    = lon_geo(i)
      phi_geo    = lat_geo(i)

      p1 = cosnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   + sin(phi_geo)*sinnp
      p1 = min(p1, 1.0)
      p1 = max(p1,-1.0)
      phi_cmp = asin( p1 )

      if( tilt.eq.90.0 .or. tilt.eq.-90.0 ) then
      p2 = sinnp*cos(lam_geo+lam_0-pi)
      else
      p2 = sinnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   - sin(phi_geo)*cosnp
      p2 = p2 / max( cos(phi_cmp),eps )
      p2 = min(p2, 1.0)
      p2 = max(p2,-1.0)
      endif
      p2 = acos( p2 )

      p3 = cos(phi_geo)*sin(lam_geo+lam_0-pi)
      if( p3.lt.0.0 ) p2 = -p2
      p2 = p2 + lam_np - pi
      lam_cmp = mod( p2+3.0*pi,2.0*pi ) - pi

c Determine Indexing Based on Computational Grid
c ----------------------------------------------
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

c Fix Longitude Index Boundaries
c ------------------------------
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


c Compute Immediate Surrounding Coordinates
c -----------------------------------------
      lam     =  lam_cmp
      phi     =  phi_cmp

c Compute and Adjust Longitude Weights
c ------------------------------------
      lam_im2 =  lon_cmp(im2(i))
      lam_im1 =  lon_cmp(im1(i))
      lam_ip0 =  lon_cmp(ip0(i))
      lam_ip1 =  lon_cmp(ip1(i))

      if( lam_im2.gt.lam_im1 ) lam_im2 = lam_im2 - 2*pi
      if( lam_im1.gt.lam_ip0 ) lam_ip0 = lam_ip0 + 2*pi
      if( lam_im1.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi
      if( lam_ip0.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi


c Compute and Adjust Latitude Weights   
c Note:  Latitude Index Boundaries are Adjusted during Interpolation
c ------------------------------------------------------------------
      phi_jm2 =  lat_cmp(jm2(i))
      phi_jm1 =  lat_cmp(jm1(i))
      phi_jp0 =  lat_cmp(jp0(i))
      phi_jp1 =  lat_cmp(jp1(i))

      if( jm2(i).eq.0    ) phi_jm2 = phi_jm1 - dphi(1)
      if( jm1(i).eq.jm   ) then
                           phi_jp0 = phi_jm1 + dphi(jm-1)
                           phi_jp1 = phi_jp0 + dphi(jm-2)
      endif
      if( jp1(i).eq.jm+1 ) phi_jp1 = phi_jp0 + dphi(jm-1)


c Bi-Linear Weights
c -----------------
              d    = (lam_ip0-lam_im1)*(phi_jp0-phi_jm1)
      wl_im1jm1(i) = (lam_ip0-lam    )*(phi_jp0-phi    )/d
      wl_ip0jm1(i) = (lam    -lam_im1)*(phi_jp0-phi    )/d
      wl_im1jp0(i) = (lam_ip0-lam    )*(phi    -phi_jm1)/d
      wl_ip0jp0(i) = (lam    -lam_im1)*(phi    -phi_jm1)/d

c Bi-Cubic Weights
c ----------------
      ap1 = ( (lam    -lam_ip0)*(lam    -lam_im1)*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_ip0)*(lam_ip1-lam_im1)*(lam_ip1-lam_im2) )
      ap0 = ( (lam_ip1-lam    )*(lam    -lam_im1)*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_ip0)*(lam_ip0-lam_im1)*(lam_ip0-lam_im2) )
      am1 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_im1)*(lam_ip0-lam_im1)*(lam_im1-lam_im2) )
      am2 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam_im1-lam    ) )
     .    / ( (lam_ip1-lam_im2)*(lam_ip0-lam_im2)*(lam_im1-lam_im2) )

      bp1 = ( (phi    -phi_jp0)*(phi    -phi_jm1)*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jp0)*(phi_jp1-phi_jm1)*(phi_jp1-phi_jm2) )
      bp0 = ( (phi_jp1-phi    )*(phi    -phi_jm1)*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jp0)*(phi_jp0-phi_jm1)*(phi_jp0-phi_jm2) )
      bm1 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jm1)*(phi_jp0-phi_jm1)*(phi_jm1-phi_jm2) )
      bm2 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi_jm1-phi    ) )
     .    / ( (phi_jp1-phi_jm2)*(phi_jp0-phi_jm2)*(phi_jm1-phi_jm2) )

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
      endif

c Interpolate Computational-Grid Quantities to Geophysical Grid Using Bi-Linear
c -----------------------------------------------------------------------------
      if( abs(norder).eq.1 ) then

      if( check  ) then
#if   (openmp)
!$omp  parallel do 
!$omp& default (shared)
!$omp& private (L,i,q_tmp)
#endif
      do L=1,lm
      do i=1,irun

      if( q_cmp( im1(i),jm1(i),L ).ne.undef  .and.
     .    q_cmp( ip0(i),jm1(i),L ).ne.undef  .and.
     .    q_cmp( im1(i),jp0(i),L ).ne.undef  .and.
     .    q_cmp( ip0(i),jp0(i),L ).ne.undef ) then

      q_tmp(i) = wl_im1jm1(i) * q_cmp( im1(i),jm1(i),L )
     .         + wl_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )
     .         + wl_im1jp0(i) * q_cmp( im1(i),jp0(i),L )
     .         + wl_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )

      else
      q_tmp(i) = undef
      endif
      enddo

c Load Temp array into Output array
c ---------------------------------
      do i=1,irun
      q_geo(i,L) = q_tmp(i)
      enddo
      enddo
      endif

      if( .not.check  ) then
#if   (openmp)
!$omp  parallel do 
!$omp& default (shared)
!$omp& private (L,i,q_tmp)
#endif
      do L=1,lm
      do i=1,irun
      q_tmp(i) = wl_im1jm1(i) * q_cmp( im1(i),jm1(i),L )
     .         + wl_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )
     .         + wl_im1jp0(i) * q_cmp( im1(i),jp0(i),L )
     .         + wl_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )
      enddo

c Load Temp array into Output array
c ---------------------------------
      do i=1,irun
      q_geo(i,L) = q_tmp(i)
      enddo
      enddo
      endif

      endif      ! End Check for Bi-Linear Interpolation


c Interpolate Computational-Grid Quantities to Geophysical Grid Using Bi-Cubic
c ----------------------------------------------------------------------------
      if( abs(norder).eq.3 ) then

      if( check  ) then
#if   (openmp)
!$omp  parallel do 
!$omp& default (shared)
!$omp& private (L,i,m,n,q_tmp)
!$omp& private (ip1_for_jp1, ip0_for_jp1, im1_for_jp1, im2_for_jp1)
!$omp& private (ip1_for_jm2, ip0_for_jm2, im1_for_jm2, im2_for_jm2)
!$omp& private (jp1_for_jp1, jm2_for_jm2)
#endif
      do L=1,lm
      do i=1,irun

      ip1_for_jp1 = ip1(i)
      ip0_for_jp1 = ip0(i)
      im1_for_jp1 = im1(i)
      im2_for_jp1 = im2(i)
      jp1_for_jp1 = jp1(i)
                m = 1

      if( jp0(i).eq.jm ) then
      ip1_for_jp1 = 1 + mod( ip1_for_jp1 + im/2 -1, im )
      ip0_for_jp1 = 1 + mod( ip0_for_jp1 + im/2 -1, im )
      im1_for_jp1 = 1 + mod( im1_for_jp1 + im/2 -1, im )
      im2_for_jp1 = 1 + mod( im2_for_jp1 + im/2 -1, im )
      jp1_for_jp1 = jm-1
      if(msgn.eq.-1) m=-1
      endif

      ip1_for_jm2 = ip1(i)
      ip0_for_jm2 = ip0(i)
      im1_for_jm2 = im1(i)
      im2_for_jm2 = im2(i)
      jm2_for_jm2 = jm2(i)
                n = 1

      if( jm1(i).eq.1 ) then
      ip1_for_jm2 = 1 + mod( ip1_for_jm2 + im/2 -1, im )
      ip0_for_jm2 = 1 + mod( ip0_for_jm2 + im/2 -1, im )
      im1_for_jm2 = 1 + mod( im1_for_jm2 + im/2 -1, im )
      im2_for_jm2 = 1 + mod( im2_for_jm2 + im/2 -1, im )
      jm2_for_jm2 = 2
      if(msgn.eq.-1) n=-1
      endif


      if( q_cmp( ip1(i),jp0(i),L ).ne.undef .and.
     .    q_cmp( ip0(i),jp0(i),L ).ne.undef .and.
     .    q_cmp( im1(i),jp0(i),L ).ne.undef .and.
     .    q_cmp( im2(i),jp0(i),L ).ne.undef .and.

     .    q_cmp( ip1(i),jm1(i),L ).ne.undef .and.
     .    q_cmp( ip0(i),jm1(i),L ).ne.undef .and.
     .    q_cmp( im1(i),jm1(i),L ).ne.undef .and.
     .    q_cmp( im2(i),jm1(i),L ).ne.undef .and.

     .    q_cmp( ip1_for_jm2,jm2_for_jm2,L ).ne.undef  .and.
     .    q_cmp( ip0_for_jm2,jm2_for_jm2,L ).ne.undef  .and.
     .    q_cmp( im1_for_jm2,jm2_for_jm2,L ).ne.undef  .and.
     .    q_cmp( im2_for_jm2,jm2_for_jm2,L ).ne.undef  .and.

     .    q_cmp( ip1_for_jp1,jp1_for_jp1,L ).ne.undef  .and.
     .    q_cmp( ip0_for_jp1,jp1_for_jp1,L ).ne.undef  .and.
     .    q_cmp( im1_for_jp1,jp1_for_jp1,L ).ne.undef  .and.
     .    q_cmp( im2_for_jp1,jp1_for_jp1,L ).ne.undef ) then

      q_tmp(i) = wc_ip1jp1(i) * q_cmp( ip1_for_jp1,jp1_for_jp1,L )*m
     .         + wc_ip0jp1(i) * q_cmp( ip0_for_jp1,jp1_for_jp1,L )*m
     .         + wc_im1jp1(i) * q_cmp( im1_for_jp1,jp1_for_jp1,L )*m
     .         + wc_im2jp1(i) * q_cmp( im2_for_jp1,jp1_for_jp1,L )*m

     .         + wc_ip1jp0(i) * q_cmp( ip1(i),jp0(i),L )
     .         + wc_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )
     .         + wc_im1jp0(i) * q_cmp( im1(i),jp0(i),L )
     .         + wc_im2jp0(i) * q_cmp( im2(i),jp0(i),L )

     .         + wc_ip1jm1(i) * q_cmp( ip1(i),jm1(i),L )
     .         + wc_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )
     .         + wc_im1jm1(i) * q_cmp( im1(i),jm1(i),L )
     .         + wc_im2jm1(i) * q_cmp( im2(i),jm1(i),L )

     .         + wc_ip1jm2(i) * q_cmp( ip1_for_jm2,jm2_for_jm2,L )*n
     .         + wc_ip0jm2(i) * q_cmp( ip0_for_jm2,jm2_for_jm2,L )*n
     .         + wc_im1jm2(i) * q_cmp( im1_for_jm2,jm2_for_jm2,L )*n
     .         + wc_im2jm2(i) * q_cmp( im2_for_jm2,jm2_for_jm2,L )*n


      elseif( q_cmp( im1(i),jm1(i),L ).ne.undef  .and.
     .        q_cmp( ip0(i),jm1(i),L ).ne.undef  .and.
     .        q_cmp( im1(i),jp0(i),L ).ne.undef  .and.
     .        q_cmp( ip0(i),jp0(i),L ).ne.undef ) then

      q_tmp(i) = wl_im1jm1(i) * q_cmp( im1(i),jm1(i),L )
     .         + wl_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )
     .         + wl_im1jp0(i) * q_cmp( im1(i),jp0(i),L )
     .         + wl_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )

      else
      q_tmp(i) = undef
      endif
      enddo

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
      do i=1,irun
      if( q_tmp(i).ne.undef .and.
     .    q_tmp(i).lt.0.0 )  then
          q_tmp(i) = 0.0
      endif
      enddo
      endif

c Load Temp array into Output array
c ---------------------------------
      do i=1,irun
      q_geo(i,L) = q_tmp(i)
      enddo
      enddo
      endif

      if( .not.check  ) then
#if   (openmp)
!$omp  parallel do 
!$omp& default (shared)
!$omp& private (L,i,m,n,q_tmp)
!$omp& private (ip1_for_jp1, ip0_for_jp1, im1_for_jp1, im2_for_jp1)
!$omp& private (ip1_for_jm2, ip0_for_jm2, im1_for_jm2, im2_for_jm2)
!$omp& private (jp1_for_jp1, jm2_for_jm2)
#endif
      do L=1,lm
      do i=1,irun

      ip1_for_jp1 = ip1(i)
      ip0_for_jp1 = ip0(i)
      im1_for_jp1 = im1(i)
      im2_for_jp1 = im2(i)
      jp1_for_jp1 = jp1(i)
                m = 1

      if( jp0(i).eq.jm ) then
      ip1_for_jp1 = 1 + mod( ip1_for_jp1 + im/2 -1, im )
      ip0_for_jp1 = 1 + mod( ip0_for_jp1 + im/2 -1, im )
      im1_for_jp1 = 1 + mod( im1_for_jp1 + im/2 -1, im )
      im2_for_jp1 = 1 + mod( im2_for_jp1 + im/2 -1, im )
      jp1_for_jp1 = jm-1
      if(msgn.eq.-1) m=-1
      endif

      ip1_for_jm2 = ip1(i)
      ip0_for_jm2 = ip0(i)
      im1_for_jm2 = im1(i)
      im2_for_jm2 = im2(i)
      jm2_for_jm2 = jm2(i)
                n = 1

      if( jm1(i).eq.1 ) then
      ip1_for_jm2 = 1 + mod( ip1_for_jm2 + im/2 -1, im )
      ip0_for_jm2 = 1 + mod( ip0_for_jm2 + im/2 -1, im )
      im1_for_jm2 = 1 + mod( im1_for_jm2 + im/2 -1, im )
      im2_for_jm2 = 1 + mod( im2_for_jm2 + im/2 -1, im )
      jm2_for_jm2 = 2
      if(msgn.eq.-1) n=-1
      endif

      q_tmp(i) = wc_ip1jp1(i) * q_cmp( ip1_for_jp1,jp1_for_jp1,L )*m
     .         + wc_ip0jp1(i) * q_cmp( ip0_for_jp1,jp1_for_jp1,L )*m
     .         + wc_im1jp1(i) * q_cmp( im1_for_jp1,jp1_for_jp1,L )*m
     .         + wc_im2jp1(i) * q_cmp( im2_for_jp1,jp1_for_jp1,L )*m

     .         + wc_ip1jp0(i) * q_cmp( ip1(i),jp0(i),L )
     .         + wc_ip0jp0(i) * q_cmp( ip0(i),jp0(i),L )
     .         + wc_im1jp0(i) * q_cmp( im1(i),jp0(i),L )
     .         + wc_im2jp0(i) * q_cmp( im2(i),jp0(i),L )

     .         + wc_ip1jm1(i) * q_cmp( ip1(i),jm1(i),L )
     .         + wc_ip0jm1(i) * q_cmp( ip0(i),jm1(i),L )
     .         + wc_im1jm1(i) * q_cmp( im1(i),jm1(i),L )
     .         + wc_im2jm1(i) * q_cmp( im2(i),jm1(i),L )

     .         + wc_ip1jm2(i) * q_cmp( ip1_for_jm2,jm2_for_jm2,L )*n
     .         + wc_ip0jm2(i) * q_cmp( ip0_for_jm2,jm2_for_jm2,L )*n
     .         + wc_im1jm2(i) * q_cmp( im1_for_jm2,jm2_for_jm2,L )*n
     .         + wc_im2jm2(i) * q_cmp( im2_for_jm2,jm2_for_jm2,L )*n

      enddo

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
      do i=1,irun
      if( q_tmp(i).ne.undef .and.
     .    q_tmp(i).lt.0.0 )  then
          q_tmp(i) = 0.0
      endif
      enddo
      endif

c Load Temp array into Output array
c ---------------------------------
      do i=1,irun
      q_geo(i,L) = q_tmp(i)
      enddo
      enddo
      endif

      endif      ! End Check for Bi-Cubic Interpolation

          deallocate (  old_dlam ,  old_dphi )
          deallocate (   old_lon ,   old_lat )
          deallocate ( wl_ip0jp0 , wl_im1jp0 )
          deallocate ( wl_ip0jm1 , wl_im1jm1 )
          deallocate ( wc_ip1jp1 , wc_ip0jp1 , wc_im1jp1 , wc_im2jp1 )
          deallocate ( wc_ip1jp0 , wc_ip0jp0 , wc_im1jp0 , wc_im2jp0 )
          deallocate ( wc_ip1jm1 , wc_ip0jm1 , wc_im1jm1 , wc_im2jm1 )
          deallocate ( wc_ip1jm2 , wc_ip0jm2 , wc_im1jm2 , wc_im2jm2 )
          deallocate (       ip1 ,       ip0 ,       im1 ,       im2 )
          deallocate (       jp1 ,       jp0 ,       jm1 ,       jm2 )

      return
      end
