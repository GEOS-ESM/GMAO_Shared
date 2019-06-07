      module m_ll2fv
      use m_mapz, only : set_eta
      use m_const, only : kappa
      use m_const, only : rgas
      use m_const, only : rvap
      use m_const, only : grav
      use m_const, only : cpm
      implicit none
      PRIVATE
      PUBLIC ll2fv
      interface ll2fv
      module procedure ll2fv_
      end interface
      CONTAINS
      subroutine ll2fv_ ( ps1,dp1,u1,v1,thv1,q1,o1,phis1,lm1,
     .                    ps2,dp2,u2,v2,t2  ,q2,o2,phis2,lm2,
     .                    im,jm,pbelow,pabove,
     .                    w1, w2 )  ! optional arguments

C***********************************************************************
C  
C  Purpose
C     Driver for remapping of target analysis to fv model levels
C
C  Argument Description
C     ps1 ...... model surface  pressure
C     dp1 ...... model pressure thickness
C     u1 ....... model zonal      wind
C     v1 ....... model meridional wind
C     thv1 ..... model virtual potential  temperature
C     q1 ....... model specific   humidity
C     o1 ....... model ozone
C     phis1 .... model surface geopotential
C     lm1 ...... model vertical   dimension
C     w1 ....... fv cloud water fraction (optional)
C
C     ps2 ...... analysis surface  pressure
C     dp2 ...... analysis pressure thickness
C     u2 ....... analysis zonal      wind
C     v2 ....... analysis meridional wind
C     t2 . ..... analysis dry-bulb temperature
C     q2 ....... analysis specific   humidity
C     o2 ....... analysis ozone
C     phis2 .... analysis surface geopotential
C     lm2 ...... analysis vertical   dimension
C     w2 ....... ncep cloud condensate mixing ratio
C
C     im ....... zonal      dimension
C     jm ....... meridional dimension
C     pbelow ... pressure below which analysis is used completely
C     pabove ... pressure above which model    is used completely
C                Note: a blend is used in-between pbelow and pabove
C                If pbelow=pabove, blending code is disabled
C
! ! REVISION HISTORY:
!   12Jan2005  Todling  Updated interface to set_eta (now from m_mapz) 
!   23Jan2006  B Zhang  Updated for total cloud water fraction
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
      use m_gmap, only : gmap
      use util, only : myout
      implicit none
      character(len=*), parameter :: myname = 'll2fv'
      integer  im,jm,lm1,lm2

c fv-DAS variables
c ----------------
      real      dp1(im,jm,lm1),  dp0(im,jm,lm1)
      real       u1(im,jm,lm1),   u0(im,jm,lm1)
      real       v1(im,jm,lm1),   v0(im,jm,lm1)
      real     thv1(im,jm,lm1), thv0(im,jm,lm1)
      real       q1(im,jm,lm1),   q0(im,jm,lm1)
      real       o1(im,jm,lm1),   o0(im,jm,lm1)
      real,optional:: w1(im,jm,lm1)
      real      ps1(im,jm),      ps0(im,jm)

      real   phis1(im,jm)
      real      ak(lm1+1)
      real      bk(lm1+1)

c Target analysis variables
c -------------------------
      real     dp2(im,jm,lm2)
      real      u2(im,jm,lm2)
      real      v2(im,jm,lm2)
      real      t2(im,jm,lm2)
      real    thv2(im,jm,lm2)
      real      q2(im,jm,lm2)
      real      o2(im,jm,lm2)
      real,optional:: w2(im,jm,lm2)
      real     ps2(im,jm)
      real   phis2(im,jm)

c Local variables
c ---------------
      real   pe0(im,jm,lm1+1)
      real   pe1(im,jm,lm1+1)
      real   pe2(im,jm,lm2+1)
      real   pk (im,jm,lm2  )
      real  pke0(im,jm,lm1+1)
      real  pke1(im,jm,lm1+1)
      real  pke2(im,jm,lm2+1)
      real  phi2(im,jm,lm2+1)

      real    cp,ptop,pbelow,pabove,pl,alf
      real    pref,tref,pkref,tstar,eps
      real    pint
      integer i,j,L
      integer ks

      real    w0(im,jm,lm1)

c Size check
c ----------
      if (present(w1) .and. present(w2) ) then
          print *, 'll2fv: w1,w2= ',present(w1),present(w2)
          if (size(w1,1)/=im .and. size(w1,2)/=jm .and. size(w1,3)/=lm1) then
              print *, 'll2fv: bad dimension for w1 '
              print *, 'll2fv: program aborting ... '
              stop
          endif
          if (size(w2,1)/=im .and. size(w2,2)/=jm .and. size(w2,3)/=lm2) then
              print *, 'll2fv: bad dimension for w2 '
              print *, 'll2fv: program aborting ... '
              stop
          endif
      endif

      cp    = cpm ! used to be getcon('CP') (ie, dry, but not consistent w/ fv)
      eps   = rvap/rgas-1.0

      call set_eta ( lm1,ks,ptop,pint,ak,bk )
      
c Compute edge-level pressures
c ----------------------------
      pe1(:,:,lm1+1) = ps1(:,:)                
      do L=lm1,1,-1
      pe1(:,:,L) = pe1(:,:,L+1)-dp1(:,:,L)    
      enddo                                  
      
c  this is needed if we have an original fv state...

c Copy input fv state into local variables
c ----------------------------------------
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif
      do j=1,jm
      do i=1,im
        ps0(i,j) = ps1(i,j)
      enddo
      enddo 

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm1
      do j=1,jm
      do i=1,im
        dp0(i,j,L) = dp1(i,j,L) 
        u0(i,j,L) = u1(i,j,L) 
        v0(i,j,L) = v1(i,j,L) 
        thv0(i,j,L) = thv1(i,j,L) 
        q0(i,j,L) = q1(i,j,L) 
        o0 (i,j,L) = o1(i,j,L) 
        if (present(w1)) then
        w0 (i,j,L) = w1(i,j,L) 
        endif

      enddo
      enddo
      enddo

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm1+1
      do j=1,jm
      do i=1,im
        pe0(i,j,L) = pe1(i,j,L)
        pke0(i,j,L) = pe0(i,j,L)**kappa
      enddo
      enddo
      enddo

c Construct target analysis pressure variables
c --------------------------------------------
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif     
      do j=1,jm
      do i=1,im
      pe2(i,j,lm2+1) = ps2(i,j)
      enddo
      enddo

      do L=lm2,1,-1
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif
      do j=1,jm
      do i=1,im
      pe2(i,j,L) = pe2(i,j,L+1) - dp2(i,j,L)
      enddo
      enddo
      enddo

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif
      do j=1,jm
      do i=1,im
      pe2(i,j,1) = 1.0  ! Set ptop = 0.01 mb (rather than 0.0 mb from NCEP)
! DSK Isn't our model top .015 mb?
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
      
c Construct target virtual potential temperature
c ----------------------------------------------
      
      do L=1,lm2
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j)
#endif
      do j=1,jm
      do i=1,im
      pk  (i,j,L) = ( pke2(i,j,L+1)-pke2(i,j,L) )/( kappa*log(pe2(i,j,L+1)/pe2(i,j,L)) )
      thv2(i,j,L) =     t2(i,j,L)*( 1.0+eps*max(0.0,q2(i,j,L)) )/pk(i,j,L)
      enddo
      enddo
      enddo
      
c Construct target analysis heights
c ---------------------------------
      phi2(:,:,lm2+1) = phis2(:,:)
      do L=lm2,1,-1
      phi2(:,:,L) = phi2(:,:,L+1) + cp*thv2(:,:,L)*( pke2(:,:,L+1)-pke2(:,:,L) )
      enddo
      
c Compute new surface pressure (ps1) consistent with fv topography
c ----------------------------------------------------------------
      
      call myout ( 6, myname, 'compute new surface pressure consistent with fv topography' )
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
           if (L.lt.1) then
              write(6,*) 'll2fv: Level variable is less than 1.'
              write(6,*) 'll2fv: i,j,l have values ',i,j,l
              write(6,*) 'll2fv: ABORT'
              call flush(6)
              call abort()
           endif
           enddo
           ps1(i,j) = pe2(i,j,L+1)*( 1+(phi2(i,j,L+1)-phis1(i,j))/(cp*thv2(i,j,L)*pke2(i,j,L+1)) )**(1.0/kappa)
      enddo
      enddo
      
c Construct fv pressure variables (pe1,pke1,dp1) using new surface pressure
c -------------------------------------------------------------------------
      
      call myout ( 6, myname, 'construct fv pressure variables using new surface pressure' )
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm1+1
      do j=1,jm
      do i=1,im
       pe1(i,j,L) = ak(L) + bk(L)*ps1(i,j)
      pke1(i,j,L) = pe1(i,j,L)**kappa
      enddo
      enddo
      enddo

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

c Map original fv state onto new eta grid
c ---------------------------------------
      call myout ( 6, myname, 'map original fv state onto new eta grid' )
! DSK What does original mean?
      if (present(w1) .and. present(w2) ) then
      call gmap ( im,jm,1, kappa,
     .            lm1, pke0, pe0, u1, v1, thv1, q1,
     .            lm1, pke1, pe1, u0, v0, thv0, q0,
     .            o_m=o1, o_n=o0,
     .            w_m=w1, w_n=w0 )                       ! <= Optional I/O cloud water frac
      else
      call gmap ( im,jm,1, kappa,
     .            lm1, pke0, pe0, u1, v1, thv1, q1,
     .            lm1, pke1, pe1, u0, v0, thv0, q0,
     .            o_m=o1, o_n=o0 )
      endif

c Map target analysis onto fv grid
c --------------------------------
      call myout ( 6, myname, 'map target analysis onto fv grid' )
      if (present(w1) .and. present(w2) ) then
      call gmap ( im,jm,1, kappa,
     .            lm2, pke2, pe2, u2, v2, thv2, q2,
     .            lm1, pke1, pe1, u1, v1, thv1, q1, 
     .            o_m=o1, o_n=o0,
     .            w_m=w2, w_n=w1 )                       ! <= Optional I/O cloud water frac
      else
      call gmap ( im,jm,1, kappa,
     .            lm2, pke2, pe2, u2, v2, thv2, q2,
     .            lm1, pke1, pe1, u1, v1, thv1, q1, 
     .            o_m=o2, o_n=o1 )
      endif

c Blend result with original fv state
c -----------------------------------
      if( pbelow.ne.pabove ) then
      call myout ( 6, myname, 'blend result with original fv state' )
      print *,' -> blend between ',pabove/100.0,' and ',pbelow/100.0,' mb'
      if (pabove .gt. pbelow) then
          print *, myname//': WARNING pabove must be a smaller number than pbelow'
          print *, myname//': in order for blending to work properly'
      endif

#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L,pl,alf)
#endif
      do L=1,lm1
      do j=1,jm
      do i=1,im
        pl=0.5*(pe1(i,j,L+1)+pe1(i,j,L))
       alf=(pl-pabove)/(pbelow-pabove)
      if( pl.lt.pabove ) then
         u1(i,j,L) =   u0(i,j,L)
         v1(i,j,L) =   v0(i,j,L)
       thv1(i,j,L) = thv0(i,j,L)
         q1(i,j,L) =   q0(i,j,L)
         o1(i,j,L) =   o0(i,j,L)
       if (present(w1) .and. present(w2) ) then
         w1(i,j,L) =   w0(i,j,L)
       endif
      else if( pl.lt.pbelow ) then
         u1(i,j,L) =   u1(i,j,L)*alf +   u0(i,j,L)*(1-alf)
         v1(i,j,L) =   v1(i,j,L)*alf +   v0(i,j,L)*(1-alf)
       thv1(i,j,L) = thv1(i,j,L)*alf + thv0(i,j,L)*(1-alf)
         q1(i,j,L) =   q1(i,j,L)*alf +   q0(i,j,L)*(1-alf)
         o1(i,j,L) =   o1(i,j,L)*alf +   o0(i,j,L)*(1-alf)
       if (present(w1) .and. present(w2) ) then
         w1(i,j,L) =   w1(i,j,L)*alf +   w0(i,j,L)*(1-alf)
       endif
      endif
      enddo
      enddo
      enddo
      endif
! DSK lm1 has a value of 55 and represents the surface.  Thus as L goes from 
! 1 to lm1, this loop travels from pl = 1.5 hPa at the model top to the 
! surface at ~100,000 hPa.  pabove is a larger number than pbelow.  Thus,
! the loop should never enter the 'else' part of the 'if' test. I verified
! this with  print statement and a run using pabove = 30 and pbelow = 10 
      return
      end subroutine ll2fv_
      end module m_ll2fv
