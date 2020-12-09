      PROGRAM  main
      implicit none
      character*120  title
      character*120  begdate

      real,   allocatable ::     uz(:,:)
      real,   allocatable ::     vz(:,:)
      real,   allocatable ::     tz(:,:)
      real,   allocatable ::     wz(:,:)
      real,   allocatable ::    thz(:,:)
      real,   allocatable ::  upvpz(:,:)
      real,   allocatable ::  upwpz(:,:)
      real,   allocatable ::  vptpz(:,:)
      real,   allocatable :: vpthpz(:,:)
      real,   allocatable ::     pl(:,:)
      real,   allocatable ::     pk(:,:)
      real,   allocatable ::   strm(:,:)
      real,   allocatable ::    res(:,:)
      real,   allocatable ::  vstar(:,:)
      real,   allocatable ::  wstar(:,:)
      real,   allocatable ::  wmean(:,:)
      real,   allocatable ::  weddy(:,:)
      real,   allocatable ::   psi1(:,:)  ! Residual Mass StreamFunction (Method 1)
      real,   allocatable ::   psi2(:,:)  ! Residual Mass StreamFunction (Method 2)
      real,   allocatable ::   psim(:,:)  ! Mean     Mass StreamFunction
      real,   allocatable ::   epfy(:,:)  ! Eliassen-Palm Flux in Northward Direction
      real,   allocatable ::   epfz(:,:)  ! Eliassen-Palm Flux in    Upward Direction
      real,   allocatable :: epfdiv(:,:)  ! Eliassen-Palm Flux Divergence
      real,   allocatable ::   vstr(:,:)
      real*4, allocatable ::    dum(:)

      real,allocatable :: upvp  (:,:)
      real,allocatable :: upwp  (:,:)
      real,allocatable :: dudp  (:,:)
      real,allocatable :: dudphi(:,:)
      real,allocatable :: psie  (:,:)
      real,allocatable :: dfdphi(:,:)
      real,allocatable :: dfdp  (:,:)
      real,allocatable :: plz   (:,:)
      real,allocatable :: delp  (:,:)

      character*120, allocatable :: arg(:)
      character*120  tag, output
      real         undef, lat0
      integer      im,jm,lm,tm
      integer      j,L,n,nt,lrec
      integer      rc,nargs,iargc

      undef = 1e15

      nargs = iargc()
      if( nargs.ne.0 ) then
          allocate( arg(nargs) )
          do n=1,nargs
          call getarg(n,arg(n))
          enddo
          do n=1,nargs
             if( trim(arg(n)).eq.'-tag' ) tag = '.' // arg(n+1)
          enddo
       else
          tag = ""
       endif

c Read TXT Files to determine resolution
c --------------------------------------
      open (10,file='LAT0' // trim(tag) // '.txt',form='formatted')
      read (10,*) lat0
      print *, 'LAT0 = ',lat0
      close(10)
      open (10,file='XDIM' // trim(tag) // '.txt',form='formatted')
      read (10,*) im
      print *, 'IM = ',im
      close(10)
      open (10,file='YDIM' // trim(tag) // '.txt',form='formatted')
      read (10,*) jm
      print *, 'JM = ',jm
      close(10)
      open (10,file='ZDIM' // trim(tag) // '.txt',form='formatted')
      read (10,*) lm
      print *, 'LM = ',lm
      close(10)
      open (10,file='TDIM' // trim(tag) // '.txt',form='formatted')
      read (10,*) tm
      print *, 'TDIM = ',tm
      close(10)
      open (10,file='BEGDATE' // trim(tag) // '.txt',form='formatted')
      read (10,*) begdate
      close(10)

      allocate(     uz(jm,lm) )
      allocate(     vz(jm,lm) )
      allocate(     tz(jm,lm) )
      allocate(     wz(jm,lm) )
      allocate(    thz(jm,lm) )
      allocate(  upvpz(jm,lm) )
      allocate(  upwpz(jm,lm) )
      allocate(  vptpz(jm,lm) )
      allocate( vpthpz(jm,lm) )
      allocate(     pl(jm,lm) )
      allocate(     pk(jm,lm) )
      allocate(   strm(jm,lm) )
      allocate(    res(jm,lm) )
      allocate(  vstar(jm,lm) )
      allocate(  wstar(jm,lm) )
      allocate(  wmean(jm,lm) )
      allocate(  weddy(jm,lm) )
      allocate(   psi1(jm,lm) )
      allocate(   psi2(jm,lm) )
      allocate(   psim(jm,lm) )
      allocate(   epfy(jm,lm) )
      allocate(   epfz(jm,lm) )
      allocate( epfdiv(jm,lm) )
      allocate(   vstr(jm,lm) )
      allocate(    dum(jm)    )
                                                                                                                                
      allocate     (upvp  (jm,LM)  )
      allocate     (upwp  (jm,LM)  )
      allocate     (dudp  (jm,LM)  )
      allocate     (dudphi(jm,LM)  )
      allocate     (psie  (jm,LM)  )
      allocate     (dfdphi(jm,LM)  )
      allocate     (dfdp  (jm,LM)  )
      allocate     (plz   (jm,LM)  )
      allocate     (delp  (jm,LM)  )

c GRADS Datasets
c --------------
      open (10,file=   'grads' // trim(tag) // '.fwrite',form='unformatted',access='direct',recl=jm*4)
      open (20,file='residual' // trim(tag) // '.data'  ,form='unformatted',access='sequential')
      open (30,file='residual' // trim(tag) // '.ctl'   ,form='formatted')

c Read data from grads.fwrite
c ---------------------------
        nt = 0
        rc = 0
      lrec = 1
      do n=1,tm
         if(n.eq.1) print *, 'Reading VZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                vz(:,L) = dum(:)
            else
                vz(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading TZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                tz(:,L) = dum(:)
            else
                tz(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading VPTPZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                vptpz(:,L) = dum(:)
            else
                vptpz(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading PL'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                pl(:,L) = dum(:)
            else
                pl(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading PK'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                pk(:,L) = dum(:)
            else
                pk(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading UZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                uz(:,L) = dum(:)
            else
                uz(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading UPVPZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                upvpz(:,L) = dum(:)
            else
                upvpz(:,L) = undef
            endif
            lrec = lrec+1
         enddo
         if(n.eq.1) print *, 'Reading UPWPZ'
         do L=1,lm
            read(10,rec=lrec,iostat=rc)  dum
            if( rc.eq.0 ) then
                upwpz(:,L) = dum(:)
            else
                upwpz(:,L) = undef
            endif
            lrec = lrec+1
         enddo

      nt = nt+1
      do L=1,lm
      do j=1,jm
         if( abs(tz(j,L)-undef).gt.0.1 ) then
                thz(j,L) = tz(j,L)/pk(j,L)          ! K/mb**kappa
         else
                thz(j,L) = undef
         endif
         if( abs(vptpz(j,L)-undef).gt.0.1 ) then
                 vpthpz(j,L) = vptpz(j,L)/pk(j,L)   ! m/sec K/mb**kappa
         else
                 vpthpz(j,L) = undef
         endif
      enddo
      enddo

c Compute Meridional Streamfunction
c ---------------------------------
      call stream ( vz,pl,jm,lm,strm,undef )
 
      call make_psi ( uz,vz,thz,upvpz,upwpz,vpthpz,pl,jm,lm,psi1,psi2,psim,epfy,epfz,epfdiv,undef,
     .                     upvp,upwp,dudp,dudphi,psie,dfdphi,dfdp,plz,delp)

c Compute Mean Vertical Velocity from Continuity
c ----------------------------------------------
      call make_w ( vz,pl,jm,lm,wz,undef )

c Compute Residual Circulation
c ----------------------------
      call residual ( vz,vpthpz,thz,wz,pl,jm,lm,res,vstar,wstar,wmean,weddy,undef )
 
      call GLAWRT  (strm  ,jm,LM,20)
      call GLAWRT  (res   ,jm,LM,20)
      call GLAWRT  (vstar ,jm,LM,20)
      call GLAWRT  (wstar ,jm,LM,20)
      call GLAWRT  (wmean ,jm,LM,20)
      call GLAWRT  (weddy ,jm,LM,20)
      call GLAWRT  (psi1  ,jm,LM,20)
      call GLAWRT  (psi2  ,jm,LM,20)
      call GLAWRT  (psim  ,jm,LM,20)
      call GLAWRT  (epfy  ,jm,LM,20)
      call GLAWRT  (epfz  ,jm,LM,20)
      call GLAWRT  (epfdiv,jm,LM,20)

      call GLAWRT  (upvp  ,jm,LM,20)
      call GLAWRT  (upwp  ,jm,LM,20)
      call GLAWRT  (dudp  ,jm,LM,20)
      call GLAWRT  (dudphi,jm,LM,20)
      call GLAWRT  (psie  ,jm,LM,20)
      call GLAWRT  (dfdphi,jm,LM,20)
      call GLAWRT  (dfdp  ,jm,LM,20)
      call GLAWRT  (plz   ,jm,LM,20)
      call GLAWRT  (delp  ,jm,LM,20)

      enddo
      close(10)

c Write Grads Control File
c ------------------------
      output = '^residual' // trim(tag) // '.data'
      title  = 'Streamfunction and Residual Circulation'

      write(30,101) trim(output),trim(title),undef,jm,lat0,2*abs(lat0)/(jm-1),lm
      do L=1,lm
      print *, 'Pressure = ',pl(1,L)
      write(30,102) pl(1,L)
      enddo
      print *, 'Finished   , nt = ',nt
      write(30,103) nt,trim(begdate),lm,lm,lm,lm,lm,lm,lm,lm,lm,lm,lm,lm,
     .                               lm,lm,lm,lm,lm,lm,lm,lm,lm
  101 format('dset  ',a,/,
     .       'title ',a,/,
     .       'options sequential ',/,
     .       'undef ',e15.6,/,
     .       'xdef  1 linear -180 1',/,
     .       'ydef ',i4,' linear ',f8.3,2x,f8.3,/,
     .       'zdef ',i3,' levels ')
  102 format(10x,f8.3)
  103 format('tdef ',i3,' linear ',a,' 1mo',/,
     .       'vars 21',/,
     .       'str   ',i3,' 0 Streamfunction',/,
     .       'res   ',i3,' 0 Residual Circulation',/,
     .       'vstar ',i3,' 0 Vstar',/,
     .       'wstar ',i3,' 0 wstar',/,
     .       'wmean  ',i3,' 0 wmean ',/,
     .       'weddy  ',i3,' 0 weddy ',/,
     .       'psi1   ',i3,' 0 Res1 streamfunction ',/,
     .       'psi2   ',i3,' 0 Res2 streamfunction ',/,
     .       'psim   ',i3,' 0 Mass streamfunction ',/,
     .       'epfy   ',i3,' 0 Eliassen-Palm flux y',/,
     .       'epfz   ',i3,' 0 Eliassen-Palm flux z',/,
     .       'epfdiv ',i3,' 0 Eliassen-Palm flux Divergence',/,
     .       'upvp   ',i3,' 0 Uprime Vprim',/,
     .       'upwp   ',i3,' 0 Uprime Omegaprime',/,
     .       'dudp   ',i3,' 0 DuDp',/,
     .       'dudphi ',i3,' 0 DuDphi',/,
     .       'psie   ',i3,' 0 Eddy Streamfunction',/,
     .       'dfdphi ',i3,' 0 DfDphi',/,
     .       'dfdp   ',i3,' 0 DfDp',/,
     .       'plz    ',i3,' 0 Pressure',/,
     .       'delp   ',i3,' 0 Pressure Thickness',/,
     .       'endvars')

      stop
      end

      SUBROUTINE ZONAL ( A,AZ,IM,JNP,undef )
      DIMENSION  A(IM,JNP), AZ(JNP)

      DO J=1,JNP
      AZ(J) = 0.0
      IC    = 0
      DO   I = 1, IM
      IF( abs(A(I,J)-UNDEF).gt.0.1 )  THEN
      AZ(J) = AZ(J) + A(I,J)
      IC    = IC + 1
      ENDIF
      ENDDO
      IF(IC.NE.0) AZ(J) = AZ(J) / IC
      IF(IC.eq.0) AZ(J) = UNDEF
      ENDDO

      RETURN
      END
      SUBROUTINE  GLAWRT  (A, IM,LM, KTP)
      real        A   (IM,LM)
      real*4      TEM (IM)
      DO  L=1,LM
      DO  I=1,IM
      if( abs(a(i,L)).gt.1.e-20 ) then
      TEM (I) = A(I,L)
      else
      TEM (I) = 0.
      endif
      ENDDO
      WRITE(KTP)  TEM
      ENDDO
      RETURN
      END

      subroutine make_psi( u0,v0,th0,upvp0,upwp0,vpthp0,p0,jm,lm,psi1,psi2,psim,epfy,epfz,epfdiv,undef,
     .                     upvp,upwp,dudp,dudphi,psie,dfdphi,dfdp,p,delp)
      use MAPL_ConstantsMod
      implicit none
      integer j,k,L,jm,lm
      real undef,dphi,a,g,pi,phi,pk0
      logical defined
       
      real    th0(jm,lm),    th(jm,lm)
      real  upvp0(jm,lm),  upvp(jm,lm)
      real  upwp0(jm,lm),  upwp(jm,lm)
      real vpthp0(jm,lm), vpthp(jm,lm)
      real     u0(jm,lm),     u(jm,lm)
      real     v0(jm,lm),     v(jm,lm)
      real     p0(jm,lm),     p(jm,lm)

      real   psi1(jm,lm)
      real   psi2(jm,lm)
      real   psim(jm,lm)
      real   epfy(jm,lm)
      real   epfz(jm,lm)
      real epfdiv(jm,lm)

      real   dudp(jm,lm)
      real   dfdp(jm,lm)
      real  dthdp(jm,lm)
      real dudphi(jm,lm)
      real dfdphi(jm,lm)
      real   psie(jm,lm)
      real   delp(jm,lm)
      real  veddy(jm,lm)
      real  vstar(jm,lm)
      real  stuff(jm,lm)
      real    the(jm,0:lm) ! theta_edge
      real    ple(jm,0:lm) !     p_edge
      real     ue(jm,0:lm) !     u_edge
      real  epfze(jm,0:lm) !  epfz_edge
      real      f(jm)
      real    dum(jm)
      integer method

c Define Constants
c ----------------
        pi = 4.*atan(1.)
      dphi = pi/(jm-1)
        a  = MAPL_RADIUS
        g  = MAPL_GRAV

      Method = 0

c Invert level index (in order to be top=>bottom)
c -----------------------------------------------
      do L=1,lm
          u(:,L) =     u0(:,lm-L+1)  ! m/sec
          v(:,L) =     v0(:,lm-L+1)  ! m/sec
          p(:,L) =     p0(:,lm-L+1)  ! mb
         th(:,L) =    th0(:,lm-L+1)  ! K/mb**kappa
       upvp(:,L) =  upvp0(:,lm-L+1)  ! m/sec m/sec
       upwp(:,L) =  upwp0(:,lm-L+1)  ! m/sec Pa/sec
      vpthp(:,L) = vpthp0(:,lm-L+1)  ! m/sec K/mb**kappa
      enddo

      pk0 = (1000.0)**(2.0/7.0)      ! mb**kappa

      where( abs(p    -undef).gt.0.1 ) ; p     =     p*100  ; endwhere
      where( abs(th   -undef).gt.0.1 ) ; th    =    th*pk0  ; endwhere
      where( abs(vpthp-undef).gt.0.1 ) ; vpthp = vpthp*pk0  ; endwhere

c Compute PLE Edge Values
c -----------------------
        ple(:,0) = max( 0.0, p(:,1) - 0.5*( p(:,2)-p(:,1) ) )
      do L=1,lm-1
      do j=1,jm
        ple(j,L) = (  p(j,L+1)+ p(j,L) )*0.5
      enddo
      enddo
      ple(:,lm) =  p(:,lm) + 0.5*( p(:,lm)-p(:,lm-1) )

      do L=1,lm
        delp(:,L) = ple(:,L)-ple(:,L-1)
      enddo

c Compute Mass Streamfunction
c ---------------------------
        pi = 4.*atan(1.)
      dphi = pi/(jm-1)
        a  = MAPL_RADIUS
        g  = MAPL_GRAV

      do L=1,LM
         dum(:) = 0.0
         do k=1,L
            where( abs(v(:,k)-undef).gt.0.1 )
                     dum(:) = dum(:) + v(:,k)*delp(:,k)
            endwhere
         enddo
         do j=1,jm
            phi       = -pi/2 + (j-1)*dphi
            psim(j,L) = 2*pi*a*cos(phi)/g * dum(j)
         enddo
      enddo

c Define Eddy Streamfunction = vpthp/dthdp
c ----------------------------------------

   !  call compute_edge( th,p,ple,jm,lm,undef,the )
      call map1_cubic( lm,p,th, lm+1,ple,the, jm, Method, undef)
      call compute_dqdp( the,delp,jm,lm,undef,dthdp )

      do L=1,lm
      do j=1,jm
         if( defined(dthdp(j,L),undef)  .and. 
     .       defined(vpthp(j,L),undef) ) then
              dthdp(j,L) = min( -0.003*pk0/100, dthdp(j,L) )
               psie(j,L) = vpthp(j,L) / dthdp(j,L)
         else
               psie(j,L) = undef
         endif
      enddo
      enddo

c Compute Veddy = D/Dp[ psie ]
c ----------------------------
      do L=2,lm-1
      do j=1,jm
      if( defined(psie(j,L+1),undef)  .and.
     .    defined(psie(j,L-1),undef) ) then
               veddy(j,L) = ( psie(j,L+1)-psie(j,L-1) )/ ( 2*(ple(j,L)-ple(j,L-1)) )
      else
               veddy(j,L) = undef
      endif
      enddo
      enddo
      do j=1,jm
      veddy(j,1)  = veddy(j,2)
      veddy(j,lm) = veddy(j,lm-1)
      enddo

c Compute Vstar = v - veddy
c -------------------------
      do L=1,lm
      do j=1,jm
         if( defined( veddy(j,L),undef)  .and.
     .       defined(     v(j,L),undef) ) then
             vstar(j,L) = v(j,L) - veddy(j,L)
         else
             vstar(j,L) = undef
         endif
      enddo
      enddo


c Construct Residual Streamfunction from Vstar (Method 1)
c -------------------------------------------------------
      do L=1,LM
         dum(:) = 0.0
         do k=1,L
            where( abs(vstar(:,k)-undef).gt.0.1 )
                     dum(:) = dum(:) + vstar(:,k)*delp(:,k)
            endwhere
         enddo
         do j=1,jm
            phi       = -pi/2 + (j-1)*dphi
            psi1(j,L) = 2*pi*a*cos(phi)/g * dum(j)
         enddo
      enddo


c Compute Residual Streamfunction (Method 2)
c ------------------------------------------
      do L=1,lm
      do j=1,jm
             phi = -pi/2 + (j-1)*dphi
         if( defined(psie(j,L),undef) ) then
             psi2(j,L) = 2*pi*a*cos(phi)/g * psie(j,L)
         else
             psi2(j,L) = undef
         endif
      enddo
      enddo

      do L=1,lm
         where( abs(psim(:,L)-undef).gt.0.1 .and. 
     .          abs(psi2(:,L)-undef).gt.0.1 )
                    psi2(:,L) = psim(:,L) - psi2(:,L)
         elsewhere
                    psi2(:,L) = undef
         endwhere
      enddo


c Compute Eliassen-Palm Flux
c --------------------------
      do j=1,jm
         phi  = -pi/2 + (j-1)*dphi
         f(j) = 2*MAPL_OMEGA*sin(phi)
      enddo

      !------------------------- Compute du/dp --------------------------------

   !  call compute_edge( u,p,ple,jm,lm,undef,ue )
      call map1_cubic( lm,p,u, lm+1,ple,ue, jm, Method, undef)
      call compute_dqdp( ue,delp,jm,lm,undef,dudp )

      !--------------------- Compute d(u*cos)/(a*cos*dphi) ---------------------

      do L=1,lm
      do j=1,jm
                  phi = -pi/2 + (j-1)*dphi
          if( defined(u(j,L),undef) ) then
                  stuff(j,L) = u(j,L)*cos(phi)
          else
                  stuff(j,L) = undef
          endif
      enddo
      enddo

      do L=1,lm
         dudphi(1 ,L) = undef
         dudphi(jm,L) = undef
      do j=2,jm-1
          phi = -pi/2 + (j-1)*dphi
          if( defined(stuff(j+1,L),undef)  .and. 
     .        defined(stuff(j-1,L),undef) ) then 
              dudphi(j,L) = ( stuff(j+1,L)-stuff(j-1,L) )/(a*cos(phi)*2*dphi)
          else
              dudphi(j,L) = undef
          endif
      enddo
      enddo

      !----------------------- Compute epfy & epfz ----------------------------

      do L=1,lm
      do j=1,jm
               phi = -pi/2 + (j-1)*dphi
          if( defined( dudp(j,L),undef)  .and. 
     .        defined( psie(j,L),undef)  .and. 
     .        defined( upvp(j,L),undef) ) then 
                       epfy(j,L) = a*cos(phi)*( dudp(j,L)*psie(j,L) - upvp(j,L) )
          else
                       epfy(j,L) = undef
          endif
          if( defined( dudphi(j,L),undef)  .and. 
     .        defined(   psie(j,L),undef)  .and. 
     .        defined(   upwp(j,L),undef) ) then 
                         epfz(j,L) = a*cos(phi)*( (f(j)-dudphi(j,L))*psie(j,L) - upwp(j,L) )
          else
                         epfz(j,L) = undef
          endif
      enddo
      enddo

      !----------------------- Compute d(epfy*cos)/(a*cos*dphi) -----------------------

      do L=1,lm
      do j=1,jm
                  phi = -pi/2 + (j-1)*dphi
          if( defined(epfy(j,L),undef) ) then
                  stuff(j,L) = epfy(j,L)*cos(phi)
          else
                  stuff(j,L) = undef
          endif
      enddo
      enddo

      do L=1,lm
         dfdphi(1 ,L) = undef
         dfdphi(jm,L) = undef
      do j=2,jm-1
         phi = -pi/2 + (j-1)*dphi
         if( defined(stuff(j+1,L),undef)  .and. 
     .       defined(stuff(j-1,L),undef) ) then 
             dfdphi(j,L) = ( stuff(j+1,L)-stuff(j-1,L) )/(a*cos(phi)*2*dphi)
         else
             dfdphi(j,L) = undef
         endif
      enddo
      enddo

      !------------------------- Compute d(epfz)/dp ---------------------------

   !  call compute_edge( epfz,p,ple,jm,lm,undef,epfze )
      call map1_cubic( lm,p,epfz, lm+1,ple,epfze, jm, Method, undef)
      call compute_dqdp( epfze,delp,jm,lm,undef,dfdp )

      !----------------------- Compute EPFlux Divergence ----------------------

      do L=1,lm
      do j=1,jm
          if( defined(   dfdp(j,L),undef)  .and. 
     .        defined( dfdphi(j,L),undef) ) then 
                       epfdiv(j,L) = dfdphi(j,L) + dfdp(j,L)
          else
                       epfdiv(j,L) = undef
          endif
      enddo
      enddo

c Invert Streamfunction for grads output (in order to be bottom=>top)
c -------------------------------------------------------------------

      call flipz( psim  ,jm,lm,1.0e-10       ,undef )
      call flipz( psi1  ,jm,lm,1.0e-10*2.4892,undef )
      call flipz( psi2  ,jm,lm,1.0e-10*2.4892,undef )
      call flipz( epfy  ,jm,lm,1.0          ,undef )
      call flipz( epfz  ,jm,lm,1.0          ,undef )
      call flipz( epfdiv,jm,lm,1.0          ,undef )

      call flipz( upvp  ,jm,lm,1.0            ,undef )
      call flipz( upwp  ,jm,lm,1.0            ,undef )
      call flipz( dudp  ,jm,lm,1.0            ,undef )
      call flipz( dudphi,jm,lm,1.0            ,undef )
      call flipz( psie  ,jm,lm,1.0            ,undef )
      call flipz( dfdphi,jm,lm,1.0            ,undef )
      call flipz( dfdp  ,jm,lm,1.0            ,undef )
      call flipz( p     ,jm,lm,1.0            ,undef )
      call flipz( delp  ,jm,lm,1.0            ,undef )

      return
      end

      subroutine flipz( q,jm,lm,scale,undef )
      implicit none
      integer j,L,jm,lm
      real undef,scale
      logical defined
      real q(jm,lm)
      real z(jm,lm)
      do L=1,lm
         where( abs(q(:,LM-L+1)-undef).gt.0.1 )
                z(:,L) = q(:,LM-L+1)*scale
         elsewhere
                z(:,L) = undef
         endwhere
      enddo
      do L=1,lm
         q(:,L) = z(:,L)
      enddo

      return
      end

      subroutine stream ( v0,p0,jm,lm,s,undef )
      use MAPL_ConstantsMod
      implicit none
      integer j,k,L,jm,lm
      real pi,dp,a,g,const,phi,undef
       
      real  v(jm,lm), v0(jm,lm)
      real  s(jm,lm)
      real p0(jm,lm),  p(jm,lm)
      real dum(jm)

      real  ple(jm,0:lm)
      real delp(jm,  lm)

c Invert VWND and P level index (in order to be top=>bottom)
c ----------------------------------------------------------
      do L=1,lm
      p(:,L) = p0(:,lm-L+1)
      v(:,L) = v0(:,lm-L+1)
      enddo

c Compute Edge Pressures and Thickness
c ------------------------------------
        ple(:,0) = max( 0.0, p(:,1) - 0.5*( p(:,2)-p(:,1) ) )
      do L=1,lm-1
        ple(:,L) = (  p(:,L)+ p(:,L+1) )*0.5
      enddo
        ple(:,lm) =  p(:,lm) + 0.5*( p(:,lm)-p(:,lm-1) )
      do L=1,lm
        delp(:,L) = ple(:,L)-ple(:,L-1)
      enddo

      pi = 4.*atan(1.)
      dp = pi/(jm-1)
      a  = MAPL_RADIUS
      g  = MAPL_GRAV

      const = 2*pi*a/g * 1.0e-8

      do k=1,lm
         dum(:) = 0.0
      do L=1,k
      do j=1,jm
      phi = -pi/2+(j-1)*dp
      if( abs(v(j,L)-undef).gt.0.1 ) then
                  dum(j) = dum(j) + v(j,L)*cos(phi)*delp(j,L)
      endif
      enddo
      enddo
         s(:,k) = dum(:)*const
      enddo

c Invert Streamfunction for grads output (in order to be bottom=>top)
c -------------------------------------------------------------------
      do k=1,lm
      do j=1,jm
      v(j,k) = s(j,lm-k+1)
      enddo
      enddo

      do k=1,lm
      do j=1,jm
      s(j,k) = v(j,k)
      enddo
      enddo

      return
      end
      subroutine residual ( v0,vpthp0,th0,w0,p0,jm,lm,res,vstar,wstar,wmean,weddy,undef )
      use MAPL_ConstantsMod
      implicit none
      integer j,k,L,jm,lm
      real pi,dp,a,g,H,ps,ts,rhos,z,phi,undef
      real airmw,runiv,cpd,rgas,akap
       
      real     v0(jm,lm),   v(jm,lm)
      real     w0(jm,lm),   w(jm,lm)
      real vpthp0(jm,lm), th0(jm,lm)
      real vpthp (jm,lm), th (jm,lm), dthdp(jm,lm)
      real stuff (jm,lm)
      real  res  (jm,lm)
      real  vtlda(jm,lm)
      real  vstar(jm,lm)
      real  wstar(jm,lm), wmean(jm,lm), weddy(jm,lm)
      real    s  (jm,lm)
      real p0(jm,lm), p(jm,lm), rho0(jm,lm)
      real   cosp(jm), dum(jm,lm)
      real ddcosp(jm,lm)
      real    the(jm,0:lm)
      real    ple(jm,0:lm)
      real   delp(jm,  lm)
      logical defined

      PARAMETER ( AIRMW  = MAPL_AIRMW  )
      PARAMETER ( RUNIV  = MAPL_RUNIV  )
      PARAMETER ( CPD    = MAPL_CP     )
      PARAMETER ( RGAS   = RUNIV/AIRMW )
      PARAMETER ( AKAP   = MAPL_KAPPA  )

c Invert v,th,vpthp, and P level index (in order to be top=>bottom)
c -----------------------------------------------------------------
      do L=1,lm
          w(:,L) =     w0(:,lm-L+1)
          v(:,L) =     v0(:,lm-L+1)
          p(:,L) =     p0(:,lm-L+1)
         th(:,L) =    th0(:,lm-L+1)
      vpthp(:,L) = vpthp0(:,lm-L+1)
      enddo

      pi = 4.*atan(1.)
      dp = pi/(jm-1)
      a  = MAPL_RADIUS
      g  = MAPL_GRAV
      H  = 7000.0
      ps = 1000.0
      ts =  240.0
      rhos = ps/(rgas*ts)

      print *, '  rhos = ',    rhos
      print *, '1/rhos = ',1.0/rhos

c Compute Mean Air Density
c ------------------------
      do L=1,lm
      do j=1,jm
             z  = -H*log(p(j,L)/ps)
      rho0(j,L) = rhos*exp(-z/H)
      enddo
      enddo
      
      do j=1,jm
      phi = -pi/2 + (j-1)*dp
      cosp(j) = cos(phi)
      enddo

c Compute Edge Pressures and Thickness
c ------------------------------------
        the(:,0) = th(:,1)
        ple(:,0) = max( 0.0, p(:,1) - 0.5*( p(:,2)-p(:,1) ) )
      do L=1,lm-1
      do j=1,jm
        ple(j,L) = (  p(j,L)+ p(j,L+1) )*0.5
                                                                    the(j,L) =   undef
        if( defined(th(j,L  ),undef)                              ) the(j,L) =   th(j,L)
        if( defined(th(j,L+1),undef)                              ) the(j,L) =   th(j,L+1)
        if( defined(th(j,L+1),undef) .and. defined(th(j,L),undef) ) the(j,L) = ( th(j,L+1)+th(j,L) )*0.5
      enddo
      enddo
      ple(:,lm) =  p(:,lm) + 0.5*( p(:,lm)-p(:,lm-1) )

                 the(:,lm) = undef
      where(  abs(th(:,lm)-undef).gt.0.1 .and. abs(the(:,lm-1)-undef).gt.0.1 )
                 the(:,lm) = the(:,lm-1) + ( th(:,lm)-the(:,lm-1) ) * ( ple(:,lm)-ple(:,lm-1) )/( p(:,lm)-ple(:,lm-1) )
      endwhere

      do L=1,lm
        delp(:,L) = ple(:,L)-ple(:,L-1)
      enddo


c Compute D(Theta)/DZ (with a forced minimum to prevent dthdz => 0)
c -----------------------------------------------------------------
      do L=1,lm
      do j=1,jm
         if( defined(the(j,L  ),undef )  .and.
     .       defined(the(j,L-1),undef ) ) then
             dthdp(j,L) = min( -0.003, ( the(j,L)-the(j,L-1) )/ delp(j,L) )
         else
             dthdp(j,L) = undef
         endif
      enddo
      enddo

c Compute Vtlda based on D(rho*vpthp/dthdz)/DZ
c --------------------------------------------
      do L=1,lm
      do j=1,jm
      if( defined(dthdp(j,L),undef) .and.
     .    defined(vpthp(j,L),undef) .and.
     .    dthdp(j,L).ne.0.0       ) then
          stuff(j,L) = rho0(j,L)*vpthp(j,L)/(p(j,L)*dthdp(j,L))
      else
          stuff(j,L) = undef
      endif
      enddo
      enddo

      do L=2,lm-1
      do j=1,jm
      if( defined(stuff(j,L+1),undef)  .and. 
     .    defined(stuff(j,L-1),undef) ) then 
               vtlda(j,L) = p(j,L)/rho0(j,L) * ( stuff(j,L+1)-stuff(j,L-1) )/ ( 2*(ple(j,L)-ple(j,L-1)) )
      else
               vtlda(j,L) = undef
      endif
      enddo
      enddo
      do j=1,jm
      vtlda(j,1)  = vtlda(j,2)
      vtlda(j,lm) = vtlda(j,lm-1)
      enddo

c Compute Vstar
c -------------
      do L=1,lm
      do j=1,jm
      if( defined( vtlda(j,L),undef)  .and.
     .    defined(     v(j,L),undef) ) then
          vstar(j,L) = v(j,L) - vtlda(j,L)
      else
          vstar(j,L) = undef
      endif
      enddo
      enddo

c Construct Residual Streamfunction from Vstar
c --------------------------------------------
      do k=1,lm
         dum(:,1) = 0.0
         do L=1,k
         do j=1,jm
         if( defined(vstar(j,L),undef) ) then
             dum(j,1) = dum(j,1) + vstar(j,L)*delp(j,L)*cosp(j)*rho0(j,L)*H/p(j,L)
         endif
         enddo
         enddo
         res(:,k) = dum(:,1)
      enddo

c Invert Streamfunction and Vstar for grads output (in order to be bottom=>top)
c -----------------------------------------------------------------------------
      do L=1,lm
      do j=1,jm
      dum(j,L) = res(j,LM-L+1)
      enddo
      enddo
      do L=1,lm
      do j=1,jm
      res(j,L) = dum(j,L)
      enddo
      enddo

      do L=1,lm
      do j=1,jm
      dum(j,L) = vstar(j,LM-L+1)
      enddo
      enddo
      do L=1,lm
      do j=1,jm
      vstar(j,L) = dum(j,L)
      enddo
      enddo

c Compute D(cos*vpthp/dthdz)/Dphi
c -------------------------------
      do L=1,lm
      do j=1,jm
      if( defined(vpthp(j,L),undef)  .and. 
     .    defined(dthdp(j,L),undef)  .and.
     .            dthdp(j,L).ne.0.0 ) then
                  stuff(j,L) = -H*cosp(j)*vpthp(j,L)/(p(j,L)*dthdp(j,L))
      else
                  stuff(j,L) = undef
      endif
      enddo
      enddo

      do L=1,lm
      do j=1,jm
      if( j.gt.1 .and. j.lt.jm ) then
          if( defined(stuff(j+1,L),undef)  .and. 
     .        defined(stuff(j-1,L),undef) ) then 
              ddcosp(j,L) = ( stuff(j+1,L)-stuff(j-1,L) )/(2*dp)
          else
              ddcosp(j,L) = undef
          endif
      else
              ddcosp(j,L) = undef
      endif
      enddo
      enddo

c Compute Wstar
c -------------
      do L=1,lm
      do j=1,jm
      if( defined(ddcosp(j,L),undef) ) then
            wmean(j,Lm-L+1) = w(j,L)
            weddy(j,Lm-L+1) = ddcosp(j,L)/(a*cosp(j))
            wstar(j,Lm-L+1) = w(j,L) + ddcosp(j,L)/(a*cosp(j))
      else
            wstar(j,Lm-L+1) = undef
            wmean(j,Lm-L+1) = undef
            weddy(j,Lm-L+1) = undef
      endif
      enddo
      enddo


      return
      end

      subroutine make_w ( v0,p0,jm,lm,w,undef )
      use MAPL_ConstantsMod
      implicit none
      integer j,k,L,jm,lm
       
      real     v0(jm,lm),    v(jm,lm)
      real     p0(jm,lm),    p(jm,lm)
      real      w(jm,lm), rho0(jm,lm)
      real      s(jm,lm), cosp(jm)
      real    dum(jm)
      real  dvcos_dphi(jm,lm)
      real         ple(jm,0:lm)
      real        delp(jm,  lm)
      logical defined
      real airmw,runiv,cpd,rgas,akap
      real pi,dp,a,g,H,ps,ts,rhos,phi,z,undef

      PARAMETER ( AIRMW  = MAPL_AIRMW  )
      PARAMETER ( RUNIV  = MAPL_RUNIV  )
      PARAMETER ( CPD    = MAPL_CP     )
      PARAMETER ( RGAS   = RUNIV/AIRMW )
      PARAMETER ( AKAP   = MAPL_KAPPA  )

c Invert v and P level index (in order to be top=>bottom)
c -------------------------------------------------------
      do L=1,lm
         v(:,L) = v0(:,lm-L+1)
         p(:,L) = p0(:,lm-L+1)
      enddo

      pi = 4.*atan(1.)
      dp = pi/(jm-1)
      a  = MAPL_RADIUS
      g  = MAPL_GRAV
      H  = 7000.0
      ps = 1000.0
      ts =  240.0
      rhos = ps/(rgas*ts)

      do L=1,lm
      do j=1,jm
      phi = -pi/2 + (j-1)*dp
      cosp(j)   =  cos(phi)
             z  = -H*log(p(j,L)/ps)
      rho0(j,L) = rhos*exp(-z/H)
      enddo
      enddo

      do L=1,lm
      do j=1,jm
      if( j.gt.1 .and. j.lt.jm ) then
          if( defined(v(j+1,L),undef)  .and.
     .        defined(v(j-1,L),undef) ) then
              dvcos_dphi(j,L) = ( v(j+1,L)*cosp(j+1)-v(j-1,L)*cosp(j-1) )/(2*dp)
          else
              dvcos_dphi(j,L) = undef
          endif
      else
              dvcos_dphi(j,L) = undef
      endif
      enddo
      enddo

c Compute Edge Pressures and Thickness
c ------------------------------------
        ple(:,0) = max( 0.0, p(:,1) - 0.5*( p(:,2)-p(:,1) ) )
      do L=1,lm-1
        ple(:,L) = (  p(:,L)+ p(:,L+1) )*0.5
      enddo
        ple(:,lm) =  p(:,lm) + 0.5*( p(:,lm)-p(:,lm-1) )

      do L=1,lm
        delp(:,L) = ple(:,L)-ple(:,L-1)
      enddo

c Construct W from Continuity
c ---------------------------
      do k=1,lm
         dum(:) = 0.0
         do L=1,k
         do j=1,jm
         phi = -pi/2+(j-1)*dp
         if( dvcos_dphi(j,L).ne.undef ) then
             dum(j) = dum(j) + dvcos_dphi(j,L)*delp(j,L)*rho0(j,L)*H/(p(j,L)*a*cosp(j))
         endif
         enddo
         enddo
         s(:,k) = dum(:)/rho0(:,k)
      enddo

c Invert Streamfunction for grads output (in order to be bottom=>top)
c -------------------------------------------------------------------
      do k=1,lm
      do j=1,jm
      w(j,k) = s(j,lm-k+1)
      enddo
      enddo

      return
      end

    ! ************************************************************************************************************

      function defined ( q,undef )
      implicit none
      logical  defined
      real     q,undef
      defined = abs(q-undef).gt.0.1*undef
      return
      end function defined

    ! ************************************************************************************************************

      subroutine compute_edge( q,p,pe,jm,lm,undef,qe )
      implicit none
      integer j,L,jm,lm
      real undef
      logical defined
      real  q(jm,  lm),  p(jm,  lm)
      real qe(jm,0:lm), pe(jm,0:lm)

      qe(:,0) = q(:,1)
      do L=1,lm-1
      do j=1,jm
                                                                     qe(j,L) =   undef
        if( defined( q(j,L  ),undef)                              )  qe(j,L) =   q(j,L)
        if( defined( q(j,L+1),undef)                              )  qe(j,L) =   q(j,L+1)
      ! if( defined( q(j,L+1),undef) .and. defined( q(j,L),undef) )  qe(j,L) = ( q(j,L+1)+ q(j,L) )*0.5

      ! Linear Interpolation to Pressure Edge
      ! -------------------------------------
        if( defined( q(j,L+1),undef) .and. defined( q(j,L),undef) ) then
               qe(j,L) = q(j,L)   + ( q(j,L+1)-q(j,L) )*( log(pe(j,L)/p(j,L)) )/( log(p(j,L+1)/p(j,L)) )
          !    qe(j,L) = q(j,L)   + ( q(j,L+1)-q(j,L) )*( pe(j,L)  - p(j,L) )/( p(j,L+1)-p(j,L) )
          ! or qe(j,L) = q(j,L+1) - ( q(j,L+1)-q(j,L) )*(  p(j,L+1)-pe(j,L) )/( p(j,L+1)-p(j,L) )
        endif

      enddo
      enddo
                  qe(:,lm) = undef
      where(  abs( q(:,lm)-undef).gt.0.1 .and. abs( qe(:,lm-1)-undef).gt.0.1 )
                  qe(:,lm) =  qe(:,lm-1) + (  q(:,lm)- qe(:,lm-1) ) * ( log(pe(:,lm)/pe(:,lm-1)) )/( log(p(:,lm)/pe(:,lm-1)) )
               !  qe(:,lm) =  qe(:,lm-1) + (  q(:,lm)- qe(:,lm-1) ) * ( pe(:,lm)-pe(:,lm-1) )/( p(:,lm)-pe(:,lm-1) )
      endwhere

      return
      end

    ! ************************************************************************************************************

      subroutine compute_dqdp( qe,dp,jm,lm,undef,dqdp )
      implicit none
      integer j,L,jm,lm
      real undef
      logical defined
      real dp(jm,  lm)
      real qe(jm,0:lm)
      real dqdp(jm,lm)

      do L=1,lm
      do j=1,jm
         if( defined(qe(j,L-1),undef) .and. defined(qe(j,L),undef) ) then
             dqdp(j,L) = ( qe(j,L)-qe(j,L-1) )/ dp(j,L)
         else
             dqdp(j,L) = undef
         endif
      enddo
      enddo

      return
      end

    ! ************************************************************************************************************

      subroutine map1_cubic( km, pe1, q1, kn, pe2, q2, jm, Method, undef)
      use MAPL_ConstantsMod
      implicit none

      real,    intent(in) :: undef
      integer, intent(in) :: Method            ! 0: Linear in P
                                               ! 1: Linear in Log(P)
                                               ! 2: Linear in P**kappa
      integer, intent(in) :: jm                ! Latitude dimension
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(jm,km)          ! pressure at mid-layers
                                               ! in the original vertical coordinate

      real, intent(in) ::  pe2(jm,kn)          ! pressure at mid-layers
                                               ! in the new vertical coordinate

      real, intent(in)   ::  q1(jm,km)         ! Field input
      real, intent(inout)::  q2(jm,kn)         ! Field output

! !DESCRIPTION:
!
!     Perform Cubic Interpolation in the vertical
!     -------------------------------------------
!      pe1: pressure associated with q1
!      pe2: pressure associated with q2
!
!-----------------------------------------------------------------------
!
      real       qx(jm,km)
      real   logpl1(jm,km)
      real   logpl2(jm,kn)
      real   dlogp1(jm,km)
      real   am2,am1,ap0,ap1,P,PLP1,PLP0,PLM1,PLM2,DLP0,DLM1,DLM2

      integer j, k, LM2,LM1,LP0,LP1
      logical defined

      real airmw,runiv,cpd,rgas,akap
      PARAMETER ( AIRMW  = MAPL_AIRMW  )
      PARAMETER ( RUNIV  = MAPL_RUNIV  )
      PARAMETER ( CPD    = MAPL_CP     )
      PARAMETER ( RGAS   = RUNIV/AIRMW )
      PARAMETER ( AKAP   = MAPL_KAPPA  )

! Initialization
! --------------

      select case (Method)

    ! Linear in P
    ! -----------
      case(0)
        do k=1,km
            qx(:,k) =  q1(:,k)
        logpl1(:,k) = pe1(:,k)
        enddo
        do k=1,kn
        logpl2(:,k) = pe2(:,k)
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

    ! Linear in Log(P)
    ! ----------------
      case(1)
        do k=1,km
            qx(:,k) = q1(:,k)
        logpl1(:,k) = log( pe1(:,k) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( pe2(:,k) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

    ! Linear in P**kappa
    ! ------------------
      case(2)
        do k=1,km
            qx(:,k) = q1(:,k)
        logpl1(:,k) = exp( akap*log( pe1(:,k) ))
        enddo
        do k=1,kn
        logpl2(:,k) = exp( akap*log( pe2(:,k) ))
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      end select

! Interpolate Q1 onto target Pressures
! ------------------------------------
      do j=1,jm
      do k=1,kn
         LM1 = 1
         LP0 = 1
         do while( LP0.le.km )
            if (logpl1(j,LP0).lt.logpl2(j,k)) then
               LP0 = LP0+1
            else
               exit
            endif
         enddo
         LM1 = max(LP0-1,1)
         LP0 = min(LP0, km)

! Extrapolate Linearly above first model level
! --------------------------------------------
         if( LM1.eq.1 .and. LP0.eq.1 ) then
                                        q2(j,k) = qx(j,1)
           if( defined(qx(j,2),undef) ) q2(j,k) = qx(j,2)
           if( defined(qx(j,1),undef) .and. defined(qx(j,2),undef) ) then
               q2(j,k) = qx(j,1) + ( qx(j,2)-qx(j,1) )*( logpl2(j,k)-logpl1(j,1) )
     .                                                /( logpl1(j,2)-logpl1(j,1) )
           endif

! Extrapolate Linearly below last model level
! -------------------------------------------
         else if( LM1.eq.km .and. LP0.eq.km ) then
                                           q2(j,k) = qx(j,km)
           if( defined(qx(j,km-1),undef) ) q2(j,k) = qx(j,km-1)
           if( defined(qx(j,km-1),undef) .and. defined(qx(j,km),undef) ) then
               q2(j,k) = qx(j,km) + ( qx(j,km)-qx(j,km-1) )*( logpl2(j,k )-logpl1(j,km  ) )
     .                                                     /( logpl1(j,km)-logpl1(j,km-1) )
           endif

! Interpolate Linearly between levels 1 => 2 and km-1 => km
! ---------------------------------------------------------
         else if( LM1.eq.1 .or. LP0.eq.km ) then
                                          q2(j,k) = qx(j,LP0)
           if( defined(qx(j,LM1),undef) ) q2(j,k) = qx(j,LM1)
           if( defined(qx(j,LP0),undef) .and. defined(qx(j,LM1),undef) ) then
             q2(j,k) = qx(j,LP0) + ( qx(j,LM1)-qx(j,LP0) )*( logpl2(j,k  )-logpl1(j,LP0) )
     .                                                    /( logpl1(j,LM1)-logpl1(j,LP0) )
           endif

! Interpolate Cubicly between other model levels
! ----------------------------------------------
         else
              LP1 = LP0+1
              LM2 = LM1-1
             P    = logpl2(j,k)
             PLP1 = logpl1(j,LP1)
             PLP0 = logpl1(j,LP0)
             PLM1 = logpl1(j,LM1)
             PLM2 = logpl1(j,LM2)
             DLP0 = dlogp1(j,LP0)
             DLM1 = dlogp1(j,LM1)
             DLM2 = dlogp1(j,LM2)

              ap1 = (P-PLP0)*(P-PLM1)*(P-PLM2)/( DLP0*(DLP0+DLM1)*(DLP0+DLM1+DLM2) )
              ap0 = (PLP1-P)*(P-PLM1)*(P-PLM2)/( DLP0*      DLM1 *(     DLM1+DLM2) )
              am1 = (PLP1-P)*(PLP0-P)*(P-PLM2)/( DLM1*      DLM2 *(DLP0+DLM1     ) )
              am2 = (PLP1-P)*(PLP0-P)*(PLM1-P)/( DLM2*(DLM1+DLM2)*(DLP0+DLM1+DLM2) )

           if( defined(qx(j,LP1),undef)  .and. 
     .         defined(qx(j,LP0),undef)  .and. 
     .         defined(qx(j,LM1),undef)  .and. 
     .         defined(qx(j,LM2),undef) ) then
             q2(j,k) = ap1*qx(j,LP1) + ap0*qx(j,LP0) + am1*qx(j,LM1) + am2*qx(j,LM2)

           else if( defined(qx(j,LP0),undef) .and. defined(qx(j,LM1),undef) ) then
             q2(j,k) = qx(j,LP0) + ( qx(j,LM1)-qx(j,LP0) )*( logpl2(j,k  )-logpl1(j,LP0) )
     .                                                    /( logpl1(j,LM1)-logpl1(j,LP0) )

           else
             q2(j,k) = undef
           endif

         endif

      enddo
      enddo

      return
      end subroutine map1_cubic
