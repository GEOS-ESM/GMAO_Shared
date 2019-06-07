      module m_fv2ll
      use m_const, only : kappa
      use m_const, only : rgas
      use m_const, only : rvap
      use m_const, only : grav
      use m_const, only : cpm
      implicit none
      PRIVATE
      PUBLIC fv2ll
      interface fv2ll
      module procedure fv2ll_
      end interface
      CONTAINS
      subroutine fv2ll_ ( ps1,dp1,u1,v1,t1  ,q1,o1,phis1,lm1,akbk,
     .                    ps2,dp2,u2,v2,thv2,q2,o2,phis2,lm2,im,jm,
     .                    havetv, hybrid,
     .                    w1, w2 )  ! optional arguments

C***********************************************************************
C  
C  Purpose
C     Driver for remapping fv_eta to NCEP sigma levels
C
C  Argument Description
C     ps1 ...... ncep surface  pressure
C     dp1 ...... ncep pressure thickness
C     u1 ....... ncep zonal      wind
C     v1 ....... ncep meridional wind
C     t1 ....... ncep dry-bulb   temperature
C     q1 ....... ncep specific   humidity
C     o1 ....... ncep ozone
C     phis1 .... ncep surface geopotential
C     lm1 ...... ncep vertical   dimension
C     akbk ..... ncep hybrid or sigma edge levels
C     w1 ....... ncep cloud condensate mixing ratio
C
C     ps2 ...... fv surface  pressure
C     dp2 ...... fv pressure thickness
C     u2 ....... fv zonal      wind
C     v2 ....... fv meridional wind
C     thv2 ..... fv virtual potential  temperature
C     q2 ....... fv specific   humidity
C     o2 ....... fv ozone
C     phis2 .... fv surface geopotential
C     lm2 ...... fv vertical   dimension
C     w2 ....... fv cloud water fraction (optional)
C
C     im ....... zonal      dimension
C     jm ....... meridional dimension
C
C     havetv ... if .true. g5 format input with tv in place of theta
C     hybrid ... if .true. output hybrid vcoord instead of sigma
C
!     29Sep2004  Guo - bug fix: gmap were called w/ w1 and w2 reversed
!     08Oct2009  Sienkiewicz - hybrid option for other eta-type grids
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
      use util, only : myout
      use m_gmap, only : gmap
      implicit none
      character(len=*), parameter :: myname = 'fv2ll'
      integer  im,jm,lm1,lm2
      logical  havetv, hybrid

c NCEP variables
c --------------
      real     dp1(im,jm,lm1)
      real      u1(im,jm,lm1)
      real      v1(im,jm,lm1)
      real      t1(im,jm,lm1)
      real    thv1(im,jm,lm1)
      real      q1(im,jm,lm1)
      real      o1(im,jm,lm1)
      real,optional:: w1(:,:,:)
      real     ps1(im,jm)

      real   phis1(im,jm)
      real    akbk(lm1+1,2)

c fv variables
c ------------
      real     dp2(im,jm,lm2)
      real      u2(im,jm,lm2)
      real      v2(im,jm,lm2)
      real    thv2(im,jm,lm2)
      real      q2(im,jm,lm2)
      real      o2(im,jm,lm2)
      real,optional:: w2(:,:,:)
      real     ps2(im,jm)
      real   phis2(im,jm)

c Local variables
c ---------------
      real   pe1(im,jm,lm1+1)
      real  pke1(im,jm,lm1+1)
      real   pe2(im,jm,lm2+1)
      real  pke2(im,jm,lm2+1)
      real  phi2(im,jm,lm2+1)
      real  ppk
      real, allocatable:: tv2(:,:,:)

      real    cp,ptop,akap,dqsdt,prs,pk,p,qs,t0,q0,dt,dq
      real    pref,tref,pkref,tstar,eps
      integer i,j,L,n
      logical rhbad, debug
      data    debug /.false./

c Set constants
c -------------
      cp    = cpm ! used to be getcon('CP') (ie, dry, but not consistent w/ fv)
      eps   = rvap/rgas-1.0
      akap  = 1.0/kappa
      ptop  = 1.0  ! 0.01 mb

c Size check
c ----------
      if (present(w1) .and. present(w2) ) then
          if (size(w1,1)/=im .and. size(w1,2)/=jm .and. size(w1,3)/=lm1) then
              print *, 'fv2ll: bad dimension for w1 '
              print *, 'fv2ll: program aborting ... '
              stop
          endif
          if (size(w2,1)/=im .and. size(w2,2)/=jm .and. size(w2,3)/=lm2) then
              print *, 'fv2ll: bad dimension for w2 '
              print *, 'fv2ll: program aborting ... '
              stop
          endif
      endif

c Iniitialize local vars in parallel 
c -------------------------------
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=lm1+1,1,-1
      do j=1,jm
      do i=1,im
      pe1(i,j,L) = 0.0
      pke1(i,j,L) = 0.0
      enddo
      enddo
      enddo

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=lm2+1,1,-1
      do j=1,jm
      do i=1,im
      pe2(i,j,L) = 0.0
      pke2(i,j,L) = 0.0
      phi2(i,j,L) = 0.0 
      enddo
      enddo
      enddo

c Construct fv pressure variables
c -------------------------------

      call myout ( 6, myname, 'construct fv pressure and heights' )
      do j=1,jm
      do i=1,im
      pe2(i,j,lm2+1) = ps2(i,j)
      enddo
      enddo

      do L=lm2,1,-1
      do j=1,jm
      do i=1,im
      pe2(i,j,L) = pe2(i,j,L+1) - dp2(i,j,L)
      enddo
      enddo
      enddo

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm2+1
      do j=1,jm
      do i=1,im
      pke2(i,j,L) = pe2(i,j,L)**kappa
      enddo
      enddo
      enddo


! Construct virtual potential temperature from virtual temperature
! ----------------------------------------------------------------
      if ( havetv ) then
        allocate(tv2(im,jm,lm2))
        do L=1,lm2
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif
        do j=1,jm
        do i=1,im
        ppk         = ( pke2(i,j,L+1)-pke2(i,j,L) )/( kappa*log(pe2(i,j,L+1)/pe2(i,j,L)) )
        tv2(i,j,L)  = thv2(i,j,L)      ! store input tv2
        thv2(i,j,L) =  tv2(i,j,L)/ppk
        enddo
        enddo
        enddo
      endif


c Construct fv heights
c --------------------
      
      phi2(:,:,lm2+1) = phis2(:,:)
      do L=lm2,1,-1
      phi2(:,:,L) = phi2(:,:,L+1) + cp*thv2(:,:,L)*( pke2(:,:,L+1)-pke2(:,:,L) )
      enddo
      
c Compute NCEP surface pressure consistent with NCEP topography
c -------------------------------------------------------------
      call myout ( 6, myname, 'compute NCEP surface pressure consistent with NCEP topography' )
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do j=1,jm
      do i=1,im
           L = lm2
           do while ( phi2(i,j,L).lt.phis1(i,j) )
           L = L-1
           enddo
           ps1(i,j) = pe2(i,j,L+1)*( 1+(phi2(i,j,L+1)-phis1(i,j))/(cp*thv2(i,j,L)*pke2(i,j,L+1)) )**(1.0/kappa)
      enddo
      enddo

c Construct NCEP pressure variables using NCEP surface pressure
c -------------------------------------------------------------
      call myout ( 6, myname, 'construct NCEP pressure variables using NCEP surface pressure' )
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      if ( .not. hybrid ) then
         do L=1,lm1+1
            do j=1,jm
               do i=1,im
                  pe1(i,j,L) = akbk(L,2)*( ps1(i,j)-pe2(i,j,1) ) + pe2(i,j,1)
                  pke1(i,j,L) = pe1(i,j,L)**kappa
               enddo
            enddo
         enddo
      else
         do L = 1,lm1+1
            do j = 1,jm
               do i = 1,im
                  pe1(i,j,L) = akbk(L,1) + akbk(L,2)*ps1(i,j)
                  pke1(i,j,L) = pe1(i,j,L)**kappa
               end do
            end do
         end do
!
!  force fill in original field's pressure for top level if new top level pressure is less than original field
!  this is to prevent strange extrapolations of top level values
         do j = 1,jm
            do i = 1,im
               if (pe2(i,j,1) > pe1(i,j,1) ) then
                  pe1(i,j,1) = pe2(i,j,1)
                  pke1(i,j,1) = pe1(i,j,1)**kappa
               end if
            end do
         end do

      end if
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm1
      do j=1,jm
      do i=1,im
       dp1(i,j,L) = pe1(i,j,L+1)-pe1(i,j,L)
      enddo
      enddo
      enddo

c Map fv_eta onto NCEP pressure levels
c ------------------------------------
      call myout ( 6, myname, 'map fv_eta onto NCEP pressure levels' )
      if (present(w1) .and. present(w2) ) then
         call gmap ( im,jm,1, kappa,
     .               lm2, pke2, pe2, u2, v2, thv2, q2,      ! <= Input  fv_eta
     .               lm1, pke1, pe1, u1, v1, thv1, q1,      ! <= Output NCEP sigma
     .               o_m=o2, o_n=o1,                        ! <= Optional ozone
     .               w_m=w2, w_n=w1 )                       ! <= Optional I/O cloud water frac
      else
         call gmap ( im,jm,1, kappa,
     .               lm2, pke2, pe2, u2, v2, thv2, q2,      ! <= Input  fv_eta
     .               lm1, pke1, pe1, u1, v1, thv1, q1,      ! <= Output NCEP sigma
     .               o_m=o2, o_n=o1 )                       ! <= Optional ozone
      endif

c Construct NCEP moisture and dry-bulb temperature
c ------------------------------------------------
      call myout ( 6, myname, 'construct NCEP moisture and dry-bulb temperature' )
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L,pk,p,prs,q0,t0,qs,dqsdt,rhbad,n,dq,dt)
#endif
      do L=1,lm1
      do j=1,jm
      do i=1,im
              pk = ( pke1(i,j,L+1)-pke1(i,j,L) )/( kappa*log(pe1(i,j,L+1)/pe1(i,j,L)) )
               p = pk**akap
             prs = p/100
              q0 = max( q1(i,j,L),0.0 )
              t0 = thv1(i,j,L)*pk/( 1.0+eps*q0 )
             call qsat (t0,prs,qs,dqsdt,.true.)

c if supersaturated, find dt & dq which equals saturation while preserving thv
c ----------------------------------------------------------------------------
            if( q0.gt.qs ) then
                   rhbad = q0/qs.gt.1.5 .and. debug
                if(rhbad) then
                   n = 0
                   print *
                   print *, 'Extremely large RH for (i,j,L): ',i,j,L,'  prs: ',prs
                   print *
                   write(6,1000) n,t0,q0,q0/qs*100,thv1(i,j,L)
                endif
                do n = 1,4 
                 dq = (qs-q0)/( 1+dqsdt*eps*t0/(1+eps*q0) )
                 dt = -eps*t0*dq/(1+eps*q0)
                 t0 = t0 + dt
                 q0 = q0 + dq
                 call qsat (t0,prs,qs,dqsdt,.true.)
                 if(rhbad) write(6,1000) n,t0,q0,q0/qs*100,t0*(1+eps*q0)/pk
                enddo
                 dq = (qs-q0)/( 1+dqsdt*eps*t0/(1+eps*q0) )
                 dt = -eps*t0*dq/(1+eps*q0)
                 t0 = t0 + dt
                 q0 = q0 + dq
                 if(rhbad) then
                    call qsat (t0,prs,qs,dqsdt,.true.)
                    write(6,1000) n,t0,q0,q0/qs*100,t0*(1+eps*q0)/pk
                 endif
            endif

       q1(i,j,L) = q0
       t1(i,j,L) = t0
      enddo
      enddo
      enddo

! Place virtual temperature back into input vector
! ------------------------------------------------
      if ( havetv ) then
        thv2=tv2
        deallocate(tv2)
      endif

 1000 format(1x,'n: ',i2,2x,'tmpu: ',f10.6,2x,'sphu: ',e12.6,2x,'rh: ',f10.6,2x,'theta_v: ',f10.6)
      return
      end subroutine fv2ll_
      end module m_fv2ll
