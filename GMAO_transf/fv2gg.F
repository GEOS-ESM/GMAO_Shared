      program  main
      implicit none

c **********************************************************************
c **********************************************************************
c ****                                                              ****
c ****   Program to create gridded NCEP files from fv_eta files     ****
c ****                                                              ****
! RTodling, 12Jan2005 Replaced getcon with m_const
c **********************************************************************
c **********************************************************************

      integer  im,jm,lm
      integer  nymd,nhms,nstep
      integer  trunc,im_out,jm_out

      character*1  blank(80), string(80)
      character*80 ncep_data, fv_data, tag, gg_rst
      data         blank /80*' '/
      equivalence( blank (01),tag )

      character*2  cnhms
      character*8  cnymd
      data         gg_rst   /'fv2gg.yyyymmdd_hhz.data'/

      equivalence ( blank (01),tag   )
      equivalence ( string(01),gg_rst)
      equivalence ( string(07),cnymd )
      equivalence ( string(16),cnhms )

c fv variables
c ------------
      real, allocatable ::   dp(:,:,:)
      real, allocatable ::    u(:,:,:)
      real, allocatable ::    v(:,:,:)
      real, allocatable ::  thv(:,:,:)
      real, allocatable ::    q(:,:,:)
      real, allocatable :: phis(:,:)
      real, allocatable ::  std(:,:)
      real, allocatable ::  lwi(:,:)
      real, allocatable ::   ts(:,:)
      real, allocatable ::   ps(:,:)

c NCEP variables
c --------------
      real, allocatable :: phisx(:,:)
      real, allocatable ::   psx(:,:)
      real, allocatable ::    ux(:,:,:)
      real, allocatable ::    vx(:,:,:)
      real, allocatable ::    tx(:,:,:)
      real, allocatable ::    qx(:,:,:)
      real, allocatable ::   dpx(:,:,:)
      real, allocatable ::  sige(:)

      integer ID,mlev,rc,grid
      integer ntime,nvars,ngatts

      character*120, allocatable :: arg(:)

      integer n,nargs,iargc

C **********************************************************************
C ****                    Initialization                            ****
C **********************************************************************

          trunc = -999
          nargs = iargc()
      if( nargs.eq.0 ) then
          call usage()
      else
          allocate ( arg(nargs) )
          do n=1,nargs
          call getarg(n,arg(n))
          enddo
          do n=1,nargs
             if( trim(arg(n)).eq.'-trunc') read(arg(n+1),*)   trunc
             if( trim(arg(n)).eq.'-date' ) read(arg(n+1),100) nymd
             if( trim(arg(n)).eq.'-time' ) read(arg(n+1),200) nhms
             if( trim(arg(n)).eq.'-mlev' ) read(arg(n+1), * ) mlev
             if( trim(arg(n)).eq.'-grid' ) read(arg(n+1), * ) grid
             if( trim(arg(n)).eq.'-eta'  )   fv_data = trim(arg(n+1))
             if( trim(arg(n)).eq.'-ncep' ) ncep_data = trim(arg(n+1))
             if( trim(arg(n)).eq.'-tag'  )       tag = trim(arg(n+1))
          enddo
          print *
          print *, 'fv_eta filename: ',trim(   fv_data)
          print *, '  ncep filename: ',trim( ncep_data)
          print *, '           nymd: ',nymd
          print *, '           nhms: ',nhms
      endif
          if( trim(tag).ne.'' ) tag = trim(tag) // '.'

      write( cnymd,100 ) nymd
      write( cnhms,300 ) nhms/10000

  100 format(i8.8)
  200 format(i6.6)
  300 format(i2.2)

C **********************************************************************
C ****                     Read input files                         ****
C **********************************************************************

      call gfio_open ( trim(fv_data),1,ID,rc )
      call gfio_diminquire ( ID,im,jm,lm,ntime,nvars,ngatts,rc )
      print *, '           rslv: ',im,jm,lm
      print *, '    output  tag: ',trim( tag )
      if( grid.eq.0 ) then
                         im_out = im
                         jm_out = jm
                      print *, '    output grid: uniform (',im_out,',',jm_out,')'
      endif
      if( grid.eq.1 ) then
                      if(trunc.eq.-999) then
                         call nglat (jm,jm_out)
                         im_out = 2*jm_out
                      else
                         jm_out = (trunc+1)*3/2
                         im_out = 2*jm_out
                      endif
                      print *, '    output grid: guassian (',im_out,',',jm_out,')'
      endif
      print *
      print *, 'Reading fv_eta file   for Date: ',nymd,' Time: ',nhms
      allocate (   dp(im,jm,lm) )
      allocate (    u(im,jm,lm) )
      allocate (    v(im,jm,lm) )
      allocate (  thv(im,jm,lm) )
      allocate (    q(im,jm,lm) )
      allocate ( phis(im,jm)    )
      allocate (  std(im,jm)    )
      allocate (  lwi(im,jm)    )
      allocate (   ts(im,jm)    )
      allocate (   ps(im,jm)    )
      call get_fveta  ( ID,ps,dp,u,v,thv,q,phis,std,lwi,ts,im,jm,lm,nymd,nhms,nstep )
      call gfio_close ( ID,rc )

      if( mlev.ne.42  .and.
     .    mlev.ne.64 ) then
      print *
      print *, 'Valid NCEP mlev values are 42,64'
      stop
      endif

      allocate (  dpx(im,jm,mlev) )
      allocate (   ux(im,jm,mlev) )
      allocate (   vx(im,jm,mlev) )
      allocate (   tx(im,jm,mlev) )
      allocate (   qx(im,jm,mlev) )
      allocate (  psx(im,jm) )
      allocate (phisx(im,jm) )
      allocate ( sige(mlev+1))

      if( mlev.eq.42 ) then
      sige(43) = 1.00000000
      sige(42) = 0.99197000
      sige(41) = 0.98273998
      sige(40) = 0.97215998
      sige(39) = 0.96006000
      sige(38) = 0.94625998
      sige(37) = 0.93061000
      sige(36) = 0.91293001
      sige(35) = 0.89306003
      sige(34) = 0.87085998
      sige(33) = 0.84619999
      sige(32) = 0.81902999
      sige(31) = 0.78930998
      sige(30) = 0.75708002
      sige(29) = 0.72245997
      sige(28) = 0.68564999
      sige(27) = 0.64691001
      sige(26) = 0.60661000
      sige(25) = 0.56515998
      sige(24) = 0.52305001
      sige(23) = 0.48076999
      sige(22) = 0.43886000
      sige(21) = 0.39780000
      sige(20) = 0.35804999
      sige(19) = 0.32001001
      sige(18) = 0.28400999
      sige(17) = 0.25029001
      sige(16) = 0.21901000
      sige(15) = 0.19025999
      sige(14) = 0.16406000
      sige(13) = 0.14036000
      sige(12) = 0.11906000
      sige(11) = 0.10005000
      sige(10) = 0.08316000
      sige(09) = 0.06824000
      sige(08) = 0.05512000
      sige(07) = 0.04362000
      sige(06) = 0.03357000
      sige(05) = 0.02482000
      sige(04) = 0.01722000
      sige(03) = 0.01063000
      sige(02) = 0.00492000
      sige(01) = 0.00000000
      endif

      if( mlev.eq.64 ) then
      sige(mlev+2- 1) = 1.00000000
      sige(mlev+2- 2) = 0.99467099
      sige(mlev+2- 3) = 0.98863202
      sige(mlev+2- 4) = 0.98180002
      sige(mlev+2- 5) = 0.97408301
      sige(mlev+2- 6) = 0.96538502
      sige(mlev+2- 7) = 0.95560300
      sige(mlev+2- 8) = 0.94463098
      sige(mlev+2- 9) = 0.93235999
      sige(mlev+2-10) = 0.91867799
      sige(mlev+2-11) = 0.90347999
      sige(mlev+2-12) = 0.88666302
      sige(mlev+2-13) = 0.86813903
      sige(mlev+2-14) = 0.84783000
      sige(mlev+2-15) = 0.82568502
      sige(mlev+2-16) = 0.80167699
      sige(mlev+2-17) = 0.77581102
      sige(mlev+2-18) = 0.74813300
      sige(mlev+2-19) = 0.71872902
      sige(mlev+2-20) = 0.68773103
      sige(mlev+2-21) = 0.65531600
      sige(mlev+2-22) = 0.62170500
      sige(mlev+2-23) = 0.58715999
      sige(mlev+2-24) = 0.55197400
      sige(mlev+2-25) = 0.51646298
      sige(mlev+2-26) = 0.48095500
      sige(mlev+2-27) = 0.44577801
      sige(mlev+2-28) = 0.41124901
      sige(mlev+2-29) = 0.37765899
      sige(mlev+2-30) = 0.34526899
      sige(mlev+2-31) = 0.31430000
      sige(mlev+2-32) = 0.28492799
      sige(mlev+2-33) = 0.25728399
      sige(mlev+2-34) = 0.23145400
      sige(mlev+2-35) = 0.20748200
      sige(mlev+2-36) = 0.18537199
      sige(mlev+2-37) = 0.16509899
      sige(mlev+2-38) = 0.14660800
      sige(mlev+2-39) = 0.12982300
      sige(mlev+2-40) = 0.11465500
      sige(mlev+2-41) = 0.10100200
      sige(mlev+2-42) = 0.08875600
      sige(mlev+2-43) = 0.07780800
      sige(mlev+2-44) = 0.06804900
      sige(mlev+2-45) = 0.05937000
      sige(mlev+2-46) = 0.05167100
      sige(mlev+2-47) = 0.04485500
      sige(mlev+2-48) = 0.03883100
      sige(mlev+2-49) = 0.03351500
      sige(mlev+2-50) = 0.02883000
      sige(mlev+2-51) = 0.02470800
      sige(mlev+2-52) = 0.02108400
      sige(mlev+2-53) = 0.01790100
      sige(mlev+2-54) = 0.01510700
      sige(mlev+2-55) = 0.01265800
      sige(mlev+2-56) = 0.01051100
      sige(mlev+2-57) = 0.00863100
      sige(mlev+2-58) = 0.00698500
      sige(mlev+2-59) = 0.00554400
      sige(mlev+2-60) = 0.00428400
      sige(mlev+2-61) = 0.00318300
      sige(mlev+2-62) = 0.00222000
      sige(mlev+2-63) = 0.00137800
      sige(mlev+2-64) = 0.00064200
      sige(mlev+2-65) = 0.00000000
      endif

      print *, 'Reading NCEP topography'
      call get_ncep_phis ( ncep_data,phisx,im,jm )

C **********************************************************************
C ****                       Remap fv Data                          ****
C **********************************************************************

      call remap ( psx,dpx,ux,vx,tx ,qx,phisx,mlev,sige,
     .             ps ,dp ,u ,v ,thv,q ,phis ,lm,im,jm )

C **********************************************************************
C ****           Write Remapped fv Data in NCEP gridded format      ****
C **********************************************************************

      gg_rst = trim(tag)//trim(gg_rst)
      print *
      print *, 'Creating NCEP data: ',trim(gg_rst)
      print *, '                    ',trim(gg_rst)//'.ctl'

     
      open (50, file=trim(gg_rst),form='unformatted',access='sequential')
      open (81, file=trim(gg_rst)//'.ctl',form='formatted',access='sequential')
      write(81,5001) gg_rst
 5001 format('dset ^',a80,/,
     .       'options yrev sequential',/,
     .       'undef  -9.99E+33',/,
     .       'title ff2gg')

      call put_ncep ( psx,dpx,ux,vx,tx,qx,phisx,sige,im,jm,mlev,nymd,nhms,im_out,jm_out,50 )

      close (50)

      stop
      end

      subroutine get_fveta ( ID,ps,dp,u,v,thv,q,phis,std,lwi,ts,im,jm,lm,nymd,nhms,nstep )
      implicit none
      integer   im,jm,lm,nymd,nhms,ID,rc,nstep
      real   dp(im,jm,lm)
      real    u(im,jm,lm)
      real    v(im,jm,lm)
      real  thv(im,jm,lm)
      real    q(im,jm,lm)
      real phis(im,jm)
      real  std(im,jm)
      real  lwi(im,jm)
      real   ts(im,jm)
      real   ps(im,jm)
      call gfio_getvar ( id,'phis'   ,nymd,nhms,im,jm,0,1,phis,rc )
      call gfio_getvar ( id,'hs_stdv',nymd,nhms,im,jm,0,1,std ,rc )
      call gfio_getvar ( id,'lwi'    ,nymd,nhms,im,jm,0,1,lwi ,rc )
      call gfio_getvar ( id,'ts'     ,nymd,nhms,im,jm,0,1,ts  ,rc )
      call gfio_getvar ( id,'ps'     ,nymd,nhms,im,jm,0,1,ps  ,rc )
      call gfio_getvar ( id,'delp'   ,nymd,nhms,im,jm,1,lm,dp ,rc )
      call gfio_getvar ( id,'uwnd'   ,nymd,nhms,im,jm,1,lm,u  ,rc )
      call gfio_getvar ( id,'vwnd'   ,nymd,nhms,im,jm,1,lm,v  ,rc )
      call gfio_getvar ( id,'theta'  ,nymd,nhms,im,jm,1,lm,thv,rc )
      call gfio_getvar ( id,'sphu'   ,nymd,nhms,im,jm,1,lm,q  ,rc )
      call dtoa ( u,u,im,jm,lm,2 )
      call dtoa ( v,v,im,jm,lm,1 )
      nstep = 100
      return
      end

      subroutine get_ncep_phis ( filename,phis,im,jm )
      use m_const, only : grav
      implicit none
      integer     im,jm,i,j
      real   phis(im,jm)
      real*4  dum(im,jm)

      character*80 filename

c surface geopotential
c --------------------
      open (30,file=trim(filename),form='unformatted',access='sequential')
      read (30) ((dum(i,j),i=1,im),j=jm,1,-1)
      close(30)

      phis(:,:) = dum(:,:)*grav

      return
      end

      subroutine put_ncep ( ps,dp,u,v,t,q,phis,sige,im,jm,mlev,nymd,nhms,im_out,jm_out,ku )
      use m_const, only : kappa
      use m_const, only : grav
      implicit none
      integer  im,jm,mlev,nymd,nhms,im_out,jm_out,ku

      real     ps(im,jm)
      real     dp(im,jm,mlev)
      real      u(im,jm,mlev)
      real      v(im,jm,mlev)
      real      t(im,jm,mlev)
      real      q(im,jm,mlev)
      real      p(im,jm,mlev)
      real     pk(im,jm,mlev)
      real     rh(im,jm,mlev)
      real    vor(im,jm,mlev)
      real    div(im,jm,mlev)
      real    tmp(im,jm,mlev)
      real     pe(im,jm,mlev+1)
      real    pke(im,jm,mlev+1)
      real   phis(im,jm)
      real   topo(im,jm)
      real   sige(mlev+1)
      real   sig (mlev)

      integer hour,day,month,year
      character*3  cmon(12)
      data cmon  /'jan','feb','mar','apr','may','jun',
     .            'jul','aug','sep','oct','nov','dec'/

      real,   parameter :: undef = -9.99E+33
      real    dum,prs,qs,kncep
      integer L,i,j,m,n

      real, allocatable :: glat(:)

      kncep = 287.05/1004.5

      topo(:,:)   = phis(:,:)/grav
       tmp(:,:,:) = undef

c compute edge-level pressures
c ----------------------------
      pe (:,:,mlev+1) = ps(:,:)
      pke(:,:,mlev+1) = ps(:,:)**kappa
      do L=mlev,1,-1
      pe (:,:,L) = pe(:,:,L+1)-dp(:,:,L)
      pke(:,:,L) = pe(:,:,L  )**kappa
      enddo

c compute mid-level pressures and sigma
c -------------------------------------
      do L=1,mlev
      pk(:,:,L) = ( pke(:,:,L+1)-pke(:,:,L) )/( kappa*log(pe(:,:,L+1)/pe(:,:,L)) )
       p(:,:,L) =    pk(:,:,L)**(1.0/kappa)
         sig(L) = ( ( sige(L+1)**(kncep+1)-sige(L)**(kncep+1) )/((kncep+1)*(sige(L+1)-sige(L))) )**(1.0/kncep)
      enddo

c compute relative humidity
c -------------------------
      do L=1,mlev
      do j=1,jm
      do i=1,im
      prs = p(i,j,L)/100
      call qsat (t(i,j,L),prs,qs,dum,.false.)
      rh(i,j,L) = 100*q(i,j,L)/qs
      enddo
      enddo
      enddo

c compute vorticity and divergence
c --------------------------------
      do L=1,mlev
      call vordiv ( u(1,1,L),v(1,1,L),vor(1,1,L),div(1,1,L),im,jm )
      enddo

c write surface fields
c --------------------
      call writit ( topo,im,jm,1,1,im_out,jm_out,ku,undef,.false. )
      call writit (   ps,im,jm,1,1,im_out,jm_out,ku,undef,.false. )

c write upper-air fields
c ----------------------
      call writit (   p,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit (  dp,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit (   t,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit (   q,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit (  rh,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit (   u,im,jm,mlev,-1,im_out,jm_out,ku,undef,.false. )
      call writit (   v,im,jm,mlev,-1,im_out,jm_out,ku,undef,.false. )
      call writit ( div,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit ( vor,im,jm,mlev, 1,im_out,jm_out,ku,undef,.false. )
      call writit ( tmp,im,jm,mlev, 1,im_out,jm_out,ku,undef,.true.  )
      call writit ( tmp,im,jm,mlev, 1,im_out,jm_out,ku,undef,.true.  )

C **********************************************************************
C ****                  Write Grads Control File                    ****
C **********************************************************************

      hour  = nhms/10000
      day   = nymd-100  *(nymd/100)
      month = nymd-10000*(nymd/10000)-day
      month = month/100
      year  = nymd/10000

      write(81,5002) im_out,360.0/im_out

      if( jm.eq.jm_out) write(81,5003) jm,180.0/(jm-1)
      if( jm.ne.jm_out) then
          allocate ( glat(jm_out)    )
          call gauss_lat_nmc (glat,jm_out)
          write(81,6003) jm_out
          do j=1,jm_out,8
             m=j
             n=j+7
            if(n.gt.jm_out) n=jm_out
            write(81,6004) (glat(i),i=m,n)
          enddo
          deallocate ( glat )
      endif

      write(81,5004) mlev,(sig(L),L=mlev,1,-1)
      write(81,5005) hour,day,cmon(month),year,
     .               mlev,mlev,mlev,mlev,mlev,
     .               mlev,mlev,mlev,mlev,mlev,mlev

      close (81)

 5002 format('xdef  ',i3,' linear     0.0 ',f9.6)
 5003 format('ydef  ',i3,' linear   -90.0 ',f9.6)
 6003 format('ydef  ',i3,' levels  ')
 6004 format(8f10.5)
 5004 format('zdef  ',i3,' levels  ',5(f8.6,3x),/,
     .    50(18x,                    5(f8.6,3x),/) )
 5005 format('tdef  1 linear ',i2.2,'Z',i2.2,a3,i4,'  6hr',/,
     .       'vars  13',/,
     .       'HS    1 99 surface orography (m)',/,
     .       'PS    1 99 surface pressure (Pa)',/,
     .       'P    ',i2,' 99 pressure (Pa)',/,
     .       'DP   ',i2,' 99 delta pressure (Pa)',/,
     .       'T    ',i2,' 99 temperature (K)',/,
     .       'Q    ',i2,' 99 specific humidity (kg/kg)',/,
     .       'RH   ',i2,' 99 relative humidity (%)',/,
     .       'U    ',i2,' 99 zonal wind (m/s)',/,
     .       'V    ',i2,' 99 meridional wind (m/s)',/,
     .       'DIV  ',i2,' 99 divergence (m/s**2)',/,
     .       'VOR  ',i2,' 99 vorticity (m/s**2)',/,
     .       'Q2   ',i2,' 99 tracer 2 (kg/kg)',/,
     .       'Q3   ',i2,' 99 tracer 3 (kg/kg)',/,
     .       'endvars')

      return
      end

      subroutine qsat (tt,p,q,dqdt,ldqdt)
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Compute Saturation Specific Humidity
C
C  INPUT:
C  ======
C    TT ......... Temperature (Kelvin)
C    P .......... Pressure (mb)
C    LDQDT ...... Logical Flag to compute QSAT Derivative
C
C  OUTPUT:
C  =======
C    Q .......... Saturation Specific Humidity
C    DQDT ....... Saturation Specific Humidity Derivative wrt Temperature
C
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      IMPLICIT NONE
      REAL TT, P, Q, DQDT
      LOGICAL LDQDT
      REAL AIRMW, H2OMW
      
      PARAMETER ( AIRMW  = 28.97      )                                         
      PARAMETER ( H2OMW  = 18.01      )                                         

      REAL ESFAC, ERFAC
      PARAMETER ( ESFAC = H2OMW/AIRMW       )
      PARAMETER ( ERFAC = (1.0-ESFAC)/ESFAC )

      real aw0, aw1, aw2, aw3, aw4, aw5, aw6
      real bw0, bw1, bw2, bw3, bw4, bw5, bw6
      real ai0, ai1, ai2, ai3, ai4, ai5, ai6
      real bi0, bi1, bi2, bi3, bi4, bi5, bi6

      real d0, d1, d2, d3, d4, d5, d6
      real e0, e1, e2, e3, e4, e5, e6
      real f0, f1, f2, f3, f4, f5, f6
      real g0, g1, g2, g3, g4, g5, g6

c ********************************************************
c ***  Polynomial Coefficients WRT Water (Lowe, 1977) ****
c ***              (Valid +50 C to -50 C)             ****
c ********************************************************

      parameter ( aw0 =  6.107799961e+00 * esfac )
      parameter ( aw1 =  4.436518521e-01 * esfac )
      parameter ( aw2 =  1.428945805e-02 * esfac )
      parameter ( aw3 =  2.650648471e-04 * esfac )
      parameter ( aw4 =  3.031240396e-06 * esfac )
      parameter ( aw5 =  2.034080948e-08 * esfac )
      parameter ( aw6 =  6.136820929e-11 * esfac )

      parameter ( bw0 = +4.438099984e-01 * esfac )
      parameter ( bw1 = +2.857002636e-02 * esfac )
      parameter ( bw2 = +7.938054040e-04 * esfac )
      parameter ( bw3 = +1.215215065e-05 * esfac )
      parameter ( bw4 = +1.036561403e-07 * esfac )
      parameter ( bw5 = +3.532421810e-10 * esfac )
      parameter ( bw6 = -7.090244804e-13 * esfac )


c ********************************************************
c ***   Polynomial Coefficients WRT Ice  (Lowe, 1977) ****
c ***              (Valid  +0 C to -50 C)             ****
c ********************************************************

      parameter ( ai0 = +6.109177956e+00 * esfac )
      parameter ( ai1 = +5.034698970e-01 * esfac )
      parameter ( ai2 = +1.886013408e-02 * esfac )
      parameter ( ai3 = +4.176223716e-04 * esfac )
      parameter ( ai4 = +5.824720280e-06 * esfac )
      parameter ( ai5 = +4.838803174e-08 * esfac )
      parameter ( ai6 = +1.838826904e-10 * esfac )

      parameter ( bi0 = +5.030305237e-01 * esfac )
      parameter ( bi1 = +3.773255020e-02 * esfac )
      parameter ( bi2 = +1.267995369e-03 * esfac )
      parameter ( bi3 = +2.477563108e-05 * esfac )
      parameter ( bi4 = +3.005693132e-07 * esfac )
      parameter ( bi5 = +2.158542548e-09 * esfac )
      parameter ( bi6 = +7.131097725e-12 * esfac )


c ********************************************************
c ***         Polynomial Coefficients WRT Ice         ****
c ***   Starr and Cox (1985) (Valid -40 C to -70 C)   ****
c ********************************************************


      parameter ( d0 = 0.535098336e+01 * esfac )
      parameter ( d1 = 0.401390832e+00 * esfac )
      parameter ( d2 = 0.129690326e-01 * esfac )
      parameter ( d3 = 0.230325039e-03 * esfac )
      parameter ( d4 = 0.236279781e-05 * esfac )
      parameter ( d5 = 0.132243858e-07 * esfac )
      parameter ( d6 = 0.314296723e-10 * esfac )

      parameter ( e0 = 0.469290530e+00 * esfac )
      parameter ( e1 = 0.333092511e-01 * esfac )
      parameter ( e2 = 0.102164528e-02 * esfac )
      parameter ( e3 = 0.172979242e-04 * esfac )
      parameter ( e4 = 0.170017544e-06 * esfac )
      parameter ( e5 = 0.916466531e-09 * esfac )
      parameter ( e6 = 0.210844486e-11 * esfac )


c ********************************************************
c ***         Polynomial Coefficients WRT Ice         ****
c ***   Starr and Cox (1985) (Valid -65 C to -95 C)   ****
c ********************************************************

      parameter ( f0 = 0.298152339e+01 * esfac )
      parameter ( f1 = 0.191372282e+00 * esfac )
      parameter ( f2 = 0.517609116e-02 * esfac )
      parameter ( f3 = 0.754129933e-04 * esfac )
      parameter ( f4 = 0.623439266e-06 * esfac )
      parameter ( f5 = 0.276961083e-08 * esfac )
      parameter ( f6 = 0.516000335e-11 * esfac )

      parameter ( g0 = 0.312654072e+00 * esfac )
      parameter ( g1 = 0.195789002e-01 * esfac )
      parameter ( g2 = 0.517837908e-03 * esfac )
      parameter ( g3 = 0.739410547e-05 * esfac )
      parameter ( g4 = 0.600331350e-07 * esfac )
      parameter ( g5 = 0.262430726e-09 * esfac )
      parameter ( g6 = 0.481960676e-12 * esfac )

      REAL        TMAX, TICE
      PARAMETER ( TMAX=323.15, TICE=273.16)
      
      REAL T, D, W, QX, DQX
      T = MIN(TT,TMAX) - TICE
      DQX = 0.
      QX  = 0.

c Fitting for temperatures above 0 degrees centigrade
c ---------------------------------------------------
      if(t.gt.0.) then
       qx = aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6)))))
      if (ldqdt)  then
      dqx = bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6)))))
      endif
      endif

c Fitting for temperatures between 0 and -40
c ------------------------------------------
      if( t.le.0. .and. t.gt.-40.0 ) then
        w = (40.0 + t)/40.0
       qx =     w *(aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6))))))
     .    + (1.-w)*(ai0+T*(ai1+T*(ai2+T*(ai3+T*(ai4+T*(ai5+T*ai6))))))
      if (ldqdt)  then
      dqx =     w *(bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6))))))
     .    + (1.-w)*(bi0+T*(bi1+T*(bi2+T*(bi3+T*(bi4+T*(bi5+T*bi6))))))
      endif
      endif

c Fitting for temperatures between -40 and -70
c --------------------------------------------
      if( t.le.-40.0 .and. t.ge.-70.0 ) then
       qx = d0+T*(d1+T*(d2+T*(d3+T*(d4+T*(d5+T*d6)))))
      if (ldqdt) then
      dqx = e0+T*(e1+T*(e2+T*(e3+T*(e4+T*(e5+T*e6)))))
      endif
      endif

c Fitting for temperatures less than -70
c --------------------------------------
      if(t.lt.-70.0) then
       qx = f0+t*(f1+t*(f2+t*(f3+t*(f4+t*(f5+t*f6)))))
      if (ldqdt) then
      dqx = g0+t*(g1+t*(g2+t*(g3+t*(g4+t*(g5+t*g6)))))
      endif
      endif

c Compute Saturation Specific Humidity
c ------------------------------------
      D = (P-ERFAC*QX)
      IF(D.LT.0.) THEN
       Q = 1.0
       IF (LDQDT)  DQDT = 0.
      ELSE
       D = 1.0 / D
       Q = MIN(QX * D,1.0)
       IF (LDQDT)  DQDT = (1.0 + ERFAC*Q) * D * DQX
      ENDIF
      RETURN
      END
#ifdef _USE_STANDALONE_
      subroutine atod ( qa,qd,im,jm,lm,itype )
 
C ******************************************************************
C ****                                                          ****
C ****  This program converts 'A' gridded data                  ****
C ****                     to 'D' gridded data.                 ****
C ****                                                          ****
C ****  The D-Grid Triplet is defined as:                       ****
C ****                                                          ****
C ****              u(i,j+1)                                    ****
C ****                |                                         ****
C ****     v(i,j)---delp(i,j)---v(i+1,j)                        ****
C ****                |                                         ****
C ****              u(i,j)                                      ****
C ****                                                          ****
C ****  Thus, v is shifted left (westward),                     ****
C ****        u is shifted down (southward)                     ****
C ****                                                          ****
C ****  An FFT shift transformation is made in x for itype = 1  ****
C ****  An FFT shift transformation is made in y for itype = 2  ****
C ****                                                          ****
C ******************************************************************

      real qa   (im,jm,lm)
      real qd   (im,jm,lm)

      real qax  (   im+2 ,lm)
      real  cx  (2*(im+2),lm)
      real qay  (   2*jm ,lm)
      real  cy  (2*(2*jm),lm)

      real    cosx (im/2), sinx(im/2)
      real    cosy (jm)  , siny(jm)
      real      trigx(3*(im+1))
      real      trigy(3*(2*jm))

      integer IFX (100)
      integer IFY (100)

      jmm1  = jm-1
      jp    = 2*jmm1

      imh   = im/2
      pi    = 4.0*atan(1.0)
      dx    = 2*pi/im
      dy    = pi/jmm1

C *********************************************************
C ****             shift left (-dx/2)                  ****
C *********************************************************

      if (itype.eq.1) then

         call fftfax (im,ifx,trigx)

         do k=1,imh 
         thx = k*dx*0.5
         cosx(k) = cos(thx)
         sinx(k) = sin(thx)
         enddo

      do j=1,jm
         do L=1,lm
         do i=1,im
         qax(i,L) = qa(i,j,L)
         enddo
         enddo
         call rfftmlt (qax,cx,trigx,ifx,1 ,im+2,im,lm,-1)

         do L=1,lm
         do k=1,imh 
         kr = 2*k+1
         ki = 2*k+2
         crprime = qax(kr,L)*cosx(k) + qax(ki,L)*sinx(k)
         ciprime = qax(ki,L)*cosx(k) - qax(kr,L)*sinx(k)
         qax(kr,L) = crprime
         qax(ki,L) = ciprime
         enddo
         enddo

         call rfftmlt (qax,cx,trigx,ifx,1 ,im+2,im,lm, 1)
         do L=1,lm
         do i=1,im
         qd(i,j,L) = qax(i,L)
         enddo
         enddo
      enddo

      endif

C *********************************************************
C ****             shift down (-dy/2)                  ****
C *********************************************************

      if (itype.eq.2) then

         call fftfax (jp,ify,trigy)

         do L=1,jmm1
         thy = L*dy*0.5
         cosy(L) = cos(thy)
         siny(L) = sin(thy)
         enddo

      do i=1,imh 
         do L=1,lm
         do j=1,jmm1
         qay(j,L)      =  qa(i,j+1,L)
         qay(j+jmm1,L) = -qa(i+imh,jm-j,L)
         enddo
         enddo

         call rfftmlt (qay,cy,trigy,ify,1 ,jp+2,jp,lm,-1)

         do L=1,lm
         do k=1,jmm1
         kr = 2*k+1
         ki = 2*k+2
         crprime = qay(kr,L)*cosy(k) + qay(ki,L)*siny(k)
         ciprime = qay(ki,L)*cosy(k) - qay(kr,L)*siny(k)
         qay(kr,L) = crprime
         qay(ki,L) = ciprime
         enddo
         enddo

         call rfftmlt (qay,cy,trigy,ify,1 ,jp+2,jp,lm, 1)

         do L=1,lm
         do j=1,jmm1
         qd(i,j+1,L)        =  qay(j,L)
         qd(i+imh,jm-j+1,L) = -qay(j+jmm1,L)
         enddo
         enddo
      enddo

      endif

      return
      end

      subroutine dtoa ( qd,qa,im,jm,lm,itype )
 
C ******************************************************************
C ****                                                          ****
C ****  This program converts 'D' gridded data                  ****
C ****                     to 'A' gridded data.                 ****
C ****                                                          ****
C ****  The D-Grid Triplet is defined as:                       ****
C ****                                                          ****
C ****              u(i,j+1)                                    ****
C ****                |                                         ****
C ****     v(i,j)---delp(i,j)---v(i+1,j)                        ****
C ****                |                                         ****
C ****              u(i,j)                                      ****
C ****                                                          ****
C ****  Thus, v is shifted right (eastward),                    ****
C ****        u is shifted up   (northward)                     ****
C ****                                                          ****
C ****  An FFT shift transformation is made in x for itype = 1  ****
C ****  An FFT shift transformation is made in y for itype = 2  ****
C ****                                                          ****
C ******************************************************************

      real qa   (im,jm,lm)
      real qd   (im,jm,lm)

      real qax  (   im+2 ,lm)
      real  cx  (2*(im+2),lm)
      real qay  (   2*jm ,lm)
      real  cy  (2*(2*jm),lm)

      real    cosx (im/2), sinx(im/2)
      real    cosy (jm)  , siny(jm)
      real      trigx(3*(im+1))
      real      trigy(3*(2*jm))

      integer IFX (100)
      integer IFY (100)

      jmm1  = jm-1
      jp    = 2*jmm1

      imh   = im/2
      pi    = 4.0*atan(1.0)
      dx    = 2*pi/im
      dy    = pi/jmm1

C *********************************************************
C ****             shift right (dx/2)                  ****
C *********************************************************

      if (itype.eq.1) then

         call fftfax (im,ifx,trigx)

         do k=1,imh 
         thx = k*dx*0.5
         cosx(k) = cos(thx)
         sinx(k) = sin(thx)
         enddo

      do j=1,jm
         do L=1,lm
         do i=1,im
         qax(i,L) = qd(i,j,L)
         enddo
         enddo
         call rfftmlt (qax,cx,trigx,ifx,1 ,im+2,im,lm,-1)

         do L=1,lm
         do k=1,imh 
         kr = 2*k+1
         ki = 2*k+2
         crprime = qax(kr,L)*cosx(k) - qax(ki,L)*sinx(k)
         ciprime = qax(ki,L)*cosx(k) + qax(kr,L)*sinx(k)
         qax(kr,L) = crprime
         qax(ki,L) = ciprime
         enddo
         enddo

         call rfftmlt (qax,cx,trigx,ifx,1 ,im+2,im,lm, 1)
         do L=1,lm
         do i=1,im
         qa(i,j,L) = qax(i,L)
         enddo
         enddo
      enddo

      endif

C *********************************************************
C ****               shift up (dy/2)                   ****
C *********************************************************

      if (itype.eq.2) then

         call fftfax (jp,ify,trigy)

         do L=1,jmm1
         thy = L*dy*0.5
         cosy(L) = cos(thy)
         siny(L) = sin(thy)
         enddo

      do i=1,imh 
         do L=1,lm
         do j=1,jmm1
         qay(j,L)      =  qd(i,j+1,L)
         qay(j+jmm1,L) = -qd(i+imh,jm-j+1,L)
         enddo
         enddo

         call rfftmlt (qay,cy,trigy,ify,1 ,jp+2,jp,lm,-1)

         do L=1,lm
         do k=1,jmm1
         kr = 2*k+1
         ki = 2*k+2
         crprime = qay(kr,L)*cosy(k) - qay(ki,L)*siny(k)
         ciprime = qay(ki,L)*cosy(k) + qay(kr,L)*siny(k)
         qay(kr,L) = crprime
         qay(ki,L) = ciprime
         enddo
         enddo

         call rfftmlt (qay,cy,trigy,ify,1 ,jp+2,jp,lm, 1)

         do L=1,lm
         do j=1,jmm1
         qa(i,j+1,L)      =  qay(j,L)
         qa(i+imh,jm-j,L) = -qay(j+jmm1,L)
         enddo
         enddo

      enddo

         do L=1,lm
         do i=1,imh
         qa(i+imh,jm,L) = -qa(i,jm,L)
         qa(i,1,L)      = -qa(i+imh,1,L)
         enddo
         enddo
      endif

      return
      end

      subroutine fftfax (n,rfax,trigs)
      real trigs(1)
      integer*4 nn
      nn=n
      call dzfftm1dui (nn,trigs)
      return
      end

      subroutine rfftmlt (x,work,trigs,rfax,inc1x,inc2x,n,lot,isign)
      real x(inc2x,lot),trigs(1)

      integer*4 iisign,in,iinc1x,iinc2x,ilot

c-----convert to i4
      iisign=isign
      iinc1x=inc1x
      iinc2x=inc2x
      in=n
      ilot=lot

c-----forward
      if(iisign.lt.0) then
          scale=1.0/float(n)
          call dzfftm1du (iisign,in,ilot,x,iinc1x,iinc2x,trigs)
          do j=1,lot
          do i=1,inc2x
          x(i,j)=scale*x(i,j)
          enddo
          enddo
          endif
c-----backward
      if(iisign.gt.0) then
          scale=float(n)
          do j=1,lot
          do i=1,inc2x
          x(i,j)=scale*x(i,j)
          enddo
          enddo
          call zdfftm1du (iisign,in,ilot,x,iinc1x,iinc2x,trigs)
          scale=1.0/float(n)
          do j=1,lot
          do i=1,inc2x
          x(i,j)=scale*x(i,j)
          enddo
          enddo
          endif

      return
      end
#endif

      SUBROUTINE VORDIV ( UZ,VZ,VOR,DIV,IMAX,JMAX )                              
C ********************************************************************          
C ****                                                            ****          
C ****  THIS PROGRAM CALCULATES VORTICITY AND DIVERGENCE          ****          
C ****  AT EACH LEVEL FOR A NON-STAGGERED A-GRID                  ****          
C ****                                                            ****          
C ****  INPUT:                                                    ****          
C ****    UZ ...... ZONAL WIND                                    ****          
C ****    VZ ...... MERIDIONAL WIND                               ****          
C ****    IMAX .... NUMBER OF LONGITUDE POINTS                    ****          
C ****    JMAX .... NUMBER OF LATITUDE  POINTS                    ****          
C ****                                                            ****          
C ****  OUTPUT:                                                   ****          
C ****    VOR  (IMAX,JMAX) .... VORTICITY                         ****          
C ****    DIV  (IMAX,JMAX) .... DIVERGENCE                        ****          
C ****                                                            ****          
C ********************************************************************          
                                                                                
      real     UZ(IMAX,JMAX)                                               
      real     VZ(IMAX,JMAX)                                               
      real    DIV(IMAX,JMAX)                                               
      real    VOR(IMAX,JMAX)                                               
                                                                                
      real   P1X (IMAX,JMAX)                                               
      real   P2X (IMAX,JMAX)                                               
      real   P1Y (IMAX,JMAX)                                               
      real   P2Y (IMAX,JMAX)                                               
      real   TMP1(IMAX,JMAX)                                               
      real   TMP2(IMAX,JMAX)                                               
      real  cosij(IMAX,JMAX)                                               
                                                                                
      PARAMETER ( ZERO    =  0.0 )                                              
      PARAMETER ( ONE     =  1.0 )                                              
      PARAMETER ( TWO     =  2.0 )                                              
      PARAMETER ( THREE   =  3.0 )                                              
      PARAMETER ( FOUR    =  4.0 )                                              
      PARAMETER ( TWELVE  = 12.0 )                                              
                                                                                
      DIMENSION MSGN(2)                                                         
                                                                                
      DATA MSGN /-1,1/                                                          
                                                                                
C *********************************************************                     
C ****          INITIALIZATION FOR DIVERGENCE          ****                     
C *********************************************************                     
                                                                                
      im      = imax
      imm1    = im-1

      IMR     = IMAX
      JMR     = JMAX-1
      IMJMR   = IMR*(JMR+1)
      IMJMR   = IMR*(JMR+1)
      IMJMM1R = IMR* JMR
      IMJMM2R = IMR*(JMR-1)

      A  = 6.372e6
      pi = 4.*atan(1.)
      dl = 2*pi/imr
      dp =   pi/jmr
                                                                                
      C11 =  FOUR / (THREE *A*IM*(ONE-COS(DP))    )                           
      C12 =  ONE  / (THREE *A*IM*(ONE-COS(TWO*DP)))                           
      CX1 =  TWO  / (THREE *A*DL)                                             
      CX2 =  ONE  / (TWELVE*A*DL)                                             
      CY1 =  TWO  / (THREE *A*DP)                                             
      CY2 =  ONE  / (TWELVE*A*DP)                                             
                                                                                
      do j=1,jmax
      phi = -pi/2.+(j-1)*dp
      cosphi = cos(phi)
      do i=1,imax
      cosij(i,j) = cosphi
      enddo
      enddo

C ********************************************************                     
C ****          CALCULATE AVERAGE QUANTITIES           ****                     
C *********************************************************                     
                                                                                
      DO   j=1,jmax                                                           
      DO   I=1,imax                                                           
      P1X (I,j) = ZERO                                                          
      P2X (I,j) = ZERO                                                          
      P1Y (I,j) = ZERO                                                          
      P2Y (I,j) = ZERO                                                          
      TMP1(I,j) = UZ(I,j)                                             
      TMP2(I,j) = VZ(I,j)*COSIJ(I,j)                                  
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      P1X(I  ,j) = TMP1(I+1,j) + TMP1(I,j)                                     
      P2X(I+1,j) = TMP1(I+2,j) + TMP1(I,j)                                     
      ENDDO                                                                     
      ENDDO                                                                     

      DO   J=2,JMR                                                             
      P1X(IM,J) = TMP1(1,J) + TMP1(IM  ,J)                                     
      P2X(IM,J) = TMP1(1,J) + TMP1(IMM1,J)                                     
      P2X( 1,J) = TMP1(2,J) + TMP1(IM  ,J)                                     
      ENDDO                                                                     
                                                                                
      DO   j=1,jmax-1                                                           
      DO   I=1,imax                                                           
      P1Y(I, j) = TMP2(I,J+1) + TMP2(I,j)                                        
      ENDDO                                                                     
      ENDDO                                                                     
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      P2Y(I, j) = TMP2(I,j+1) + TMP2(I,j-1)                                        
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
      DO   I=1,IMR                                                             
      P2Y(I,  1) = ZERO                                                        
      P2Y(I,jmax) = ZERO                                                        
      ENDDO                                                                     
                                                                                
C *********************************************************                     
C ****        CALCULATE HORIZONTAL DIVERGENCE          ****                     
C *********************************************************                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      TMP1(I+1,j) = ( P1X(I+1,j)-P1X(I,j) )*CX1                               
     .            - ( P2X(I+2,j)-P2X(I,j) )*CX2                               
      ENDDO                                                                     
      ENDDO                                                                     
      DO    J=2,JMR                                                             
      TMP1( 1,J) = ( P1X( 1,J)-P1X(IM  ,J) )*CX1                              
     .           - ( P2X( 2,J)-P2X(IM  ,J) )*CX2                              
      TMP1(IM,J) = ( P1X(IM,J)-P1X(IMM1,J) )*CX1                              
     .           - ( P2X( 1,J)-P2X(IMM1,J) )*CX2                              
      ENDDO                                                                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      TMP2(I,j) = ( P1Y(I,j)  -P1Y(I,j-1) )*CY1                                  
     .          - ( P2Y(I,j+1)-P2Y(I,j-1) )*CY2                                  
      DIV (I,j) = ( TMP1(I,j)+TMP2(I,j) )/(cosij(i,j))                        
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
C *********************************************************                     
C ****     CALCULATE HORIZONTAL DIVERGENCE AT POLES    ****                     
C *********************************************************                     
                                                                                
      DO 100 M=1,2                                                              
      JPOLE  = 1 + (M-1)*(jmax-1)                                                
      JPH    = 1 + (M-1)*(jmax-2)                                               
      JSTAR1 = 2 + (M-1)*(jmax-3)                                               
      JSTAR2 = 3 + (M-1)*(jmax-5)
                                                                                
      SUM11  = ZERO                                                             
      SUM12  = ZERO                                                             
      DO    I=1,IMR                                                            
      SUM11  = SUM11 + P1Y(I,JPH   )                                           
      SUM12  = SUM12 + P2Y(I,JSTAR1)                                           
      ENDDO                                                                     
                                                                                
      DO  I=1,IMR                                                            
      DIV(I,JPOLE) = - MSGN(M) * ( C11*SUM11 - C12*SUM12 )                    
      ENDDO                                                                     
                                                                                
  100 CONTINUE                                                                  
                                                                                
C ********************************************************                     
C ****          CALCULATE AVERAGE QUANTITIES           ****                     
C *********************************************************                     
                                                                                
      DO   j=1,jmax                                                           
      DO   I=1,imax                                                           
      P1X (I,j) = ZERO                                                          
      P2X (I,j) = ZERO                                                          
      P1Y (I,j) = ZERO                                                          
      P2Y (I,j) = ZERO                                                          
      TMP1(I,j) = VZ(I,j)                                             
      TMP2(I,j) = UZ(I,j)*COSIJ(I,j)                                  
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      P1X(I  ,j) = TMP1(I+1,j) + TMP1(I,j)                                     
      P2X(I+1,j) = TMP1(I+2,j) + TMP1(I,j)                                     
      ENDDO                                                                     
      ENDDO                                                                     

      DO   J=2,JMR                                                             
      P1X(IM,J) = TMP1(1,J) + TMP1(IM  ,J)                                     
      P2X(IM,J) = TMP1(1,J) + TMP1(IMM1,J)                                     
      P2X( 1,J) = TMP1(2,J) + TMP1(IM  ,J)                                     
      ENDDO                                                                     
                                                                                
      DO   j=1,jmax-1                                                           
      DO   I=1,imax                                                           
      P1Y(I, j) = TMP2(I,J+1) + TMP2(I,j)                                        
      ENDDO                                                                     
      ENDDO                                                                     
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      P2Y(I, j) = TMP2(I,j+1) + TMP2(I,j-1)                                        
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
      DO   I=1,IMR                                                             
      P2Y(I,  1) = ZERO                                                        
      P2Y(I,jmax) = ZERO                                                        
      ENDDO                                                                     
                                                                                
C *********************************************************                     
C ****        CALCULATE HORIZONTAL VORTICITY           ****                     
C *********************************************************                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      TMP1(I+1,j) = ( P1X(I+1,j)-P1X(I,j) )*CX1                               
     .            - ( P2X(I+2,j)-P2X(I,j) )*CX2                               
      ENDDO                                                                     
      ENDDO                                                                     
      DO    J=2,JMR                                                             
      TMP1( 1,J) = ( P1X( 1,J)-P1X(IM  ,J) )*CX1                              
     .           - ( P2X( 2,J)-P2X(IM  ,J) )*CX2                              
      TMP1(IM,J) = ( P1X(IM,J)-P1X(IMM1,J) )*CX1                              
     .           - ( P2X( 1,J)-P2X(IMM1,J) )*CX2                              
      ENDDO                                                                     
                                                                                
      DO   j=2,jmax-1                                                           
      DO   I=1,imax                                                           
      TMP2(I,j) = ( P1Y(I,j)  -P1Y(I,j-1) )*CY1                                  
     .          - ( P2Y(I,j+1)-P2Y(I,j-1) )*CY2                                  
      VOR (I,j) = ( TMP1(I,j)-TMP2(I,j) )/(cosij(i,j))                        
      ENDDO                                                                     
      ENDDO                                                                     
                                                                                
C *********************************************************                     
C ****     CALCULATE HORIZONTAL VORTICITY  AT POLES    ****                     
C *********************************************************                     
                                                                                
      DO 200 M=1,2                                                              
      JPOLE  = 1 + (M-1)*(jmax-1)                                                
      JPH    = 1 + (M-1)*(jmax-2)                                               
      JSTAR1 = 2 + (M-1)*(jmax-3)                                               
      JSTAR2 = 3 + (M-1)*(jmax-5)
                                                                                
      SUM11  = ZERO                                                             
      SUM12  = ZERO                                                             
      DO    I=1,IMR                                                            
      SUM11  = SUM11 + P1Y(I,JPH   )                                           
      SUM12  = SUM12 + P2Y(I,JSTAR1)                                           
      ENDDO                                                                     
                                                                                
      DO  I=1,IMR                                                            
      VOR(I,JPOLE) =   MSGN(M) * ( C11*SUM11 - C12*SUM12 )                    
      ENDDO                                                                     
                                                                                
  200 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

      subroutine remap  ( ps1,dp1,u1,v1,t1  ,q1,phis1,lm1,sige,
     .                    ps2,dp2,u2,v2,thv2,q2,phis2,lm2,im,jm )

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
C     phis1 .... ncep surface geopotential
C     lm1 ...... ncep vertical   dimension
C     sige ..... ncep sigma edge levels
C
C     ps2 ...... fv surface  pressure
C     dp2 ...... fv pressure thickness
C     u2 ....... fv zonal      wind
C     v2 ....... fv meridional wind
C     thv2 ..... fv virtual potential  temperature
C     q2 ....... fv specific   humidity
C     phis2 .... fv surface geopotential
C     lm2 ...... fv vertical   dimension
C
C     im ....... zonal      dimension
C     jm ....... meridional dimension
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use m_const, only: kappa
      use m_const, only: rgas
      use m_const, only: rvap
      use m_const, only: grav
      use m_const, only: cpm
      implicit none
      integer  im,jm,lm1,lm2

c NCEP variables
c --------------
      real     dp1(im,jm,lm1)
      real      u1(im,jm,lm1)
      real      v1(im,jm,lm1)
      real      t1(im,jm,lm1)
      real    thv1(im,jm,lm1)
      real      q1(im,jm,lm1)
      real     ps1(im,jm)

      real   phis1(im,jm)
      real    sige(lm1+1)

c fv variables
c ------------
      real     dp2(im,jm,lm2)
      real      u2(im,jm,lm2)
      real      v2(im,jm,lm2)
      real    thv2(im,jm,lm2)
      real      q2(im,jm,lm2)
      real     ps2(im,jm)
      real   phis2(im,jm)

c Local variables
c ---------------
      real   pe1(im,jm,lm1+1)
      real  pke1(im,jm,lm1+1)
      real   pe2(im,jm,lm2+1)
      real  pke2(im,jm,lm2+1)
      real  phi2(im,jm,lm2+1)

      real    cp,akap,dqsdt,prs,pk,p,qs,t0,q0,dt,dq
      real    rgas,pref,tref,pkref,tstar,eps
      integer i,j,L,n
      logical rhbad, debug
      data    debug /.false./

c Set constants
c -------------
      cp    = cpm ! _RT getcon('CP') dry cp not correct for fv
      eps   = rvap/rgas-1.0
      akap  = 1.0/kappa

c Construct fv pressure variables
c -------------------------------
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

c Construct fv heights
c --------------------
      phi2(:,:,lm2+1) = phis2(:,:)
      do L=lm2,1,-1
      phi2(:,:,L) = phi2(:,:,L+1) + cp*thv2(:,:,L)*( pke2(:,:,L+1)-pke2(:,:,L) )
      enddo
      
c Compute NCEP surface pressure consistent with NCEP topography
c -------------------------------------------------------------
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

c Construct NCEP pressure variables using NCEP surface pressure and fv Model Top pressure
c ---------------------------------------------------------------------------------------
#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L)
#endif
      do L=1,lm1+1
      do j=1,jm
      do i=1,im
       pe1(i,j,L) = sige(L)*( ps1(i,j)-pe2(i,j,1) ) + pe2(i,j,1)
      pke1(i,j,L) = pe1(i,j,L)**kappa
      enddo
      enddo
      enddo

      do L=1,lm1
      do j=1,jm
      do i=1,im
       dp1(i,j,L) = pe1(i,j,L+1)-pe1(i,j,L)
      enddo
      enddo
      enddo

c Map fv_eta onto NCEP sigma
c --------------------------
      call gmap ( im,jm,1, kappa,
     .            lm2,  pke2,  pe2, u2,  v2,  thv2,  q2,   ! <= Input  fv_eta
     .            lm1,  pke1,  pe1, u1,  v1,  thv1,  q1 )  ! <= Output NCEP sigma

c Construct NCEP moisture and dry-bulb temperature
c ------------------------------------------------
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

 1000 format(1x,'n: ',i2,2x,'tmpu: ',f10.6,2x,'sphu: ',e12.6,2x,'rh: ',f10.6,2x,'theta_v: ',f10.6)
      return
      end

      subroutine usage()
      print *, "Usage:  "
      print *
      print *, " fv2gg.x [-eta  fv_eta_filename]"
      print *, "         [-ncep   ncep_filename]"
      print *, "         [-mlev   ncep vertical resolution]"
      print *, "         [-tag    output tag]"
      print *, "         [-grid   0 or 1]"
      print *, "         [-trunc  spectral truncation]"
      print *, "         [-date   nymd]"
      print *, "         [-time   nhms]"
      print *
      print *, "where:"
      print *
      print *, "  -eta fv_eta_filename:  Filename of fv background in analysis eta format"
      print *, "  -ncep  ncep_filename:  Filename of sample (eg.previous) NCEP data containing PHIS"
      print *, "  -mlev           mlev:  Number of NCEP sigma levels (42 or 64)"
      print *, "  -tag      output_tag:  Tag identifier for output data"
      print *, "  -grid              0:  For Uniform   grid output"
      print *, "                     1:  For Gaussian  grid output"
      print *, "  -trunc           nnn:  Spectral Truncation for Guassian Grid (Default is nearest 2**n > jm)"
      print *
      print *, "  -date        nymd:  Desired date in yyyymmdd format"
      print *, "  -time        nhms:  Desired time in   hhmmss format"
      print *
      call exit(7)
      end

      subroutine nglat (jm,jm2)
      implicit none
      integer jm,jm2
      integer m,n
      m = jm
      n = 1
      do while ( m/2.gt.0 )
      n = n+1
      m = m/2
      enddo
      jm2 = 2**n
      return
      end

      subroutine writit (q,im,jm,lm,msgn,im_out,jm_out,ku,undef,flag)
      integer   im,jm,lm,msgn,im_out,jm_out
      integer   i,j,L,m,n
      real   q (im,jm,lm)
      real  pi,dl,dp,lon,undef
      logical flag

      real*4,allocatable ::  dum(:,:)
      real,  allocatable ::   q2(:,:,:)
      real,  allocatable :: glat(:)
      real,  allocatable :: dlam(:)
      real,  allocatable :: dphi(:)
      real,  allocatable :: lats(:)
      real,  allocatable :: lons(:)

c uniform grid
c ------------
      if( im.eq.im_out .and. jm.eq.jm_out) then
      allocate ( dum(im,jm) )
      do L=lm,1,-1
      dum(:,:) = q(:,:,L)
      write(ku) ( (dum(i,j),i=1,im),j=jm,1,-1 )
      enddo
      deallocate ( dum )
      return
      endif

c gaussian grid
c -------------
      if( im.ne.im_out .or. jm.ne.jm_out) then

      pi = 4.0*atan(1.0)
      dl = 2.0*pi/im
      dp = pi/(jm-1)

      allocate ( dphi(jm) )
      allocate ( dlam(im) )

      dlam(:) = dl  ! Uniform grid input
      dphi(:) = dp  ! Uniform grid input

      allocate ( glat(jm_out)    )
      allocate ( lats(im_out*jm_out) )
      allocate ( lons(im_out*jm_out) )
      allocate (  dum(im_out,jm_out) )
      allocate (   q2(im_out,jm_out,lm) )

      call gauss_lat_nmc (glat,jm_out)
      dl = 2.0*pi/im_out

      loc = 0
      do j=1,jm_out
      do i=1,im_out
      loc = loc + 1
      lons (loc) = -pi + (i-1)*dl
      enddo
      enddo

      loc = 0
      do j=1,jm_out
      do i=1,im_out
      loc = loc + 1
      lats (loc) = glat(j)*pi/180.0
      enddo
      enddo

      call interp_h ( q,im,jm,lm,
     .                dlam,dphi,0.0,90.0,0.0,
     .                q2,im_out*jm_out,lons,lats,
     .                msgn,3,flag,undef )

      do L=lm,1,-1
      dum(:,:) = q2(:,:,L)
      write(ku) ( (dum(i,j),i=1,im_out),j=jm_out,1,-1 )
      enddo

      deallocate ( dlam,dphi,glat,lats,lons,q2,dum )
      return
      endif

      end

      subroutine gauss_lat_nmc(gaul,k)
      implicit double precision (a-h,o-z)
      dimension a(500)
      real gaul(1)
      save
      esp=1.d-14
      c=(1.d0-(2.d0/3.14159265358979d0)**2)*0.25d0
      fk=k
      kk=k/2
      call bsslz1(a,kk)
      do 30 is=1,kk
      xz=cos(a(is)/sqrt((fk+0.5d0)**2+c))
      iter=0
   10 pkm2=1.d0
      pkm1=xz
      iter=iter+1
      if(iter.gt.10) go to 70
      do 20 n=2,k
      fn=n
      pk=((2.d0*fn-1.d0)*xz*pkm1-(fn-1.d0)*pkm2)/fn
      pkm2=pkm1
   20 pkm1=pk
      pkm1=pkm2
      pkmrk=(fk*(pkm1-xz*pk))/(1.d0-xz**2)
      sp=pk/pkmrk
      xz=xz-sp
      avsp=abs(sp)
      if(avsp.gt.esp) go to 10
      a(is)=xz
   30 continue
      if(k.eq.kk*2) go to 50
      a(kk+1)=0.d0
      pk=2.d0/fk**2
      do 40 n=2,k,2
      fn=n
   40 pk=pk*fn**2/(fn-1.d0)**2
   50 continue
      do 60 n=1,kk
      l=k+1-n
      a(l)=-a(n)
   60 continue
      radi=180./(4.*atan(1.))
      do 211 n=1,k
      gaul(n)=acos(a(n))*radi-90.0
  211 continue
      return
   70 write(6,6000)
 6000 format(//5x,14herror in gauaw//)
      stop
      end

      subroutine bsslz1(bes,n)
      implicit double precision (a-h,o-z)
      dimension bes(n)
      dimension bz(50)
      data pi/3.14159265358979d0/
      data bz  / 2.4048255577d0, 5.5200781103d0,
     $  8.6537279129d0,11.7915344391d0,14.9309177086d0,18.0710639679d0,
     $ 21.2116366299d0,24.3524715308d0,27.4934791320d0,30.6346064684d0,
     $ 33.7758202136d0,36.9170983537d0,40.0584257646d0,43.1997917132d0,
     $ 46.3411883717d0,49.4826098974d0,52.6240518411d0,55.7655107550d0,
     $ 58.9069839261d0,62.0484691902d0,65.1899648002d0,68.3314693299d0,
     $ 71.4729816036d0,74.6145006437d0,77.7560256304d0,80.8975558711d0,
     $ 84.0390907769d0,87.1806298436d0,90.3221726372d0,93.4637187819d0,
     $ 96.6052679510d0,99.7468198587d0,102.888374254d0,106.029930916d0,
     $ 109.171489649d0,112.313050280d0,115.454612653d0,118.596176630d0,
     $ 121.737742088d0,124.879308913d0,128.020877005d0,131.162446275d0,
     $ 134.304016638d0,137.445588020d0,140.587160352d0,143.728733573d0,
     $ 146.870307625d0,150.011882457d0,153.153458019d0,156.295034268d0/
      nn=n
      if(n.le.50) go to 12
      bes(50)=bz(50)
      do 5 j=51,n
    5 bes(j)=bes(j-1)+pi
      nn=49
   12 do 15 j=1,nn
   15 bes(j)=bz(j)
      return
      end

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
      integer im,jm,lm,irun,norder,msgn
      logical check

      real      q_geo(irun,lm)
      real    lon_geo(irun)
      real    lat_geo(irun)

      real    q_cmp(im,jm,lm)
      real     dlam(im)
      real     dphi(jm)

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
      real    rotation , tilt , precession
      real    lam_geo, lam_cmp
      real    phi_geo, phi_cmp
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
      allocate ( wc_ip1jp1(irun) , wc_ip0jp1(irun) , wc_im1jp1(irun) , wc_im2jp1(irun) )
      allocate ( wc_ip1jp0(irun) , wc_ip0jp0(irun) , wc_im1jp0(irun) , wc_im2jp0(irun) )
      allocate ( wc_ip1jm1(irun) , wc_ip0jm1(irun) , wc_im1jm1(irun) , wc_im2jm1(irun) )
      allocate ( wc_ip1jm2(irun) , wc_ip0jm2(irun) , wc_im1jm2(irun) , wc_im2jm2(irun) )
      allocate (       ip1(irun) ,       ip0(irun) ,       im1(irun) ,       im2(irun) )
      allocate (       jp1(irun) ,       jp0(irun) ,       jm1(irun) ,       jm2(irun) )
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
          allocate ( wc_ip1jp1(irun) , wc_ip0jp1(irun) , wc_im1jp1(irun) , wc_im2jp1(irun) )
          allocate ( wc_ip1jp0(irun) , wc_ip0jp0(irun) , wc_im1jp0(irun) , wc_im2jp0(irun) )
          allocate ( wc_ip1jm1(irun) , wc_ip0jm1(irun) , wc_im1jm1(irun) , wc_im2jm1(irun) )
          allocate ( wc_ip1jm2(irun) , wc_ip0jm2(irun) , wc_im1jm2(irun) , wc_im2jm2(irun) )
          allocate (       ip1(irun) ,       ip0(irun) ,       im1(irun) ,       im2(irun) )
          allocate (       jp1(irun) ,       jp0(irun) ,       jm1(irun) ,       jm2(irun) )
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
