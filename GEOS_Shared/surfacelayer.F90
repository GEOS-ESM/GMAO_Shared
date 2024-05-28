! $Id$
! $Name$
#include "unused_dummy.H"
module sfclayer

! !USES:

use MAPL
use DragCoefficientsMod

implicit none
private
PUBLIC louissurface,z0sea,helfsurface,psi,phi,linadj,zcsub

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: louissurface
! !INTERFACE:
  subroutine louissurface(ISTYPE,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI, &
                               Z0,DZ,CM,CN,RI,ZT,ZQ,CH,CQ,UUU,UCN,RE,DCH,DCQ)
      integer,           intent(IN ) :: N
      integer,           intent(IN ) :: ISTYPE
      real,    intent(IN)    :: UU (:)
      real,    intent(IN)    :: WW (:,:)
      real,    intent(IN)    :: PS (:)
      real,    intent(IN)    :: TA (:)
      real,    intent(IN)    :: TS (:,:)
      real,    intent(IN)    :: QA (:)
      real,    intent(IN)    :: QS (:,:)
      real,    intent(IN)    :: PCU(:)
      real,    intent(IN)    :: LAI(:)
      real,    intent(INOUT) :: Z0 (:,:)
      real,    intent(IN)    :: DZ (:)
      real,    intent(INOUT) :: CM (:,:)
      real,    intent(OUT)   :: CN (:)
      real,    intent(OUT)   :: RI (:)
      real,    intent(OUT)   :: ZT (:)
      real,    intent(OUT)   :: ZQ (:)
      real,    intent(OUT)   :: CH (:,:)
      real,    intent(OUT)   :: CQ (:,:)
      real,    intent(OUT)   :: UUU (:)
      real,    intent(OUT)   :: UCN (:)
      real,    intent(OUT)   :: RE  (:)
      real,    optional, intent(OUT)   :: DCH (:,:)
      real,    optional, intent(OUT)   :: DCQ (:,:)


      real, parameter :: OCEANICEZ0 = 1.0e-3
      real, parameter :: LAKEZ0     = 1.0E-5
      real, parameter :: LAKEICEZ0  = 0.01
      real, parameter :: LANDICEZ0  = 0.002
      integer, parameter :: ICE   = 1
      integer, parameter :: WATER = 2
      integer,parameter  :: FSNW=4  !  Snowcover subtile
      integer  :: irun

      real, allocatable  :: UST(:)
      real, allocatable  :: TVS(:)
      real, allocatable  :: URA(:)
      real, allocatable  :: TVA(:)

      real :: LAICRT

      irun = size(UU,1)
      LAICRT = 0.3

      allocate( UST(1:irun)  )
      allocate( TVS(1:irun)  )
      allocate( URA(1:irun)  )
      allocate( TVA(1:irun)  )

      if(ISTYPE==1 .or. ISTYPE==3) then
        UCN = PCU*8640. !  cm/day
        UCN = ((19.8*UCN*UCN)/(1.5 + UCN + UCN*UCN))**0.4 !  REDELSPERGER et al. (2000)
        UCN = UCN*UCN
!       UCN = 0.0  ! Remove Redelsperger
      endif
!
      TVA = TA*( 1.0 + MAPL_VIREPS*QA)
      if(ISTYPE==3)TVA = TA*(1+MAPL_VIREPS*QA)
!
      if(ISTYPE==1 .or. ISTYPE==3) then
        UUU = max( sqrt(UU*UU + WW(:,N) + UCN), MAPL_USMIN )
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        UUU = max(UU,MAPL_USMIN)
      endif

      if(ISTYPE==1) then
      URA = UUU*( PS/(MAPL_RGAS*TVA) )
      elseif (ISTYPE==2 .or. ISTYPE==3 .or. ISTYPE==4) then
      URA = UUU*PS/(MAPL_RGAS*TVA)
      endif
      TVS = TS(:,N)*(1+MAPL_VIREPS*QS(:,N))

!  Estimate of friction velocity from previous step's drag coefficient
!---------------------------------------------------------------------

      UST = UUU*sqrt(CM(:,N)/URA)
      UST = min( max(UST,1.e-6) , 5.0 )

!  Iterate for aerodynamic roughness, based on previous step's stability
!-----------------------------------------------------------------------

      if(ISTYPE==1) then
        if (N==WATER) then
          call Z0SEA (Z0(:,N),UST,DZ)
        else
          Z0(:,N) = OCEANICEZ0
        end if
      elseif (ISTYPE==2) then
        if(N==ICE) then
          Z0(:,N) = LAKEICEZ0
        else
          Z0(:,N) = LAKEZ0
        end if
      elseif (ISTYPE==4) then
        Z0(:,N) = LANDICEZ0
      end if

! Get a new drag coefficient for the updated roughness
!-----------------------------------------------------

      call GETCDM(TVA,UUU,DZ,TVS,Z0(:,N), CM(:,N),CN,RI)
      if(ISTYPE==1) UST = UUU*sqrt(CM(:,N))

!  Reynolds number at the top of the viscous sublayer
!----------------------------------------------------

      if(ISTYPE==1) then
        RE  = max(min(Z0(:,N)*UST/MAPL_NUAIR, 50.),0.1)
      elseif (ISTYPE==3) then
        RE = max(min(Z0(:,N)*(sqrt(CM(:,N))*UUU) / MAPL_NUAIR, 1000.),0.4)
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        RE = min(Z0(:,N) * (sqrt(CM(:,N))*UUU) / MAPL_NUAIR, 1000.)
      end if


!  Calculate scalar roughnesses
!------------------------------

      if(ISTYPE==1) then
        if(N==WATER) then
! This is based on Zilitinkevich et al.
!         ZT = Z0(:,N) * exp(MAPL_KARMAN*(3.2 - 4.0*sqrt(max(RE,0.1))))
! This is based on Garrett 1992
!        ZT = Z0(:,N) * exp(2.0 - 2.48*RE**.25)
! This is as used in ecmwf model, based on:
!     Beljaars, A. C. M., 1994: 
!    The parametrization of surface fluxes in large-scale models under
!    free convection. Q. J. R. Meteorol. Soc., 121, 255-270.
!
          ZT = Z0(:,N) * (0.4/RE)
        else
! This is used over snow and bare soil in CSM
          ZT = Z0(:,N) * exp(-0.13*sqrt(RE))
        end if
      elseif (ISTYPE==2) then
        if(N==ICE) then 
          ZT = Z0(:,N) * exp(-0.13*sqrt(RE))
        else
          ZT = Z0(:,N) * (0.4/RE)
        end if
      elseif (ISTYPE==3) then
! Viscous sub-layer effects on Z0 vanish as the surface
!  becomes very vegetated (LAI -> infinity), assuming they are handled
!  in the vegetation resistance to heat and moisture fluxes.

!!!  Note: Comment out special code for SNOW condition to alleviate excessive night-time cooling
!!!  -------------------------------------------------------------------------------------------
!!!     if (N==FSNW) then
!!!       ZT = Z0(:,N)*exp( -0.13 * sqrt(RE) )
!!!       ZQ = ZT
!!!     else
          ZT = 1.0 - (1.0-exp(-0.13*sqrt(RE)))*exp(-max(LAI,LAICRT))
          where(LAI<LAICRT)
!!  Brutsaert, W., 1982: Evaporation into the Atmosphere. D. Reidel, 299pp.
            ZT = ZT + (exp( -2.5*sqrt(sqrt(RE)) + log(7.4) ) - ZT)*(LAICRT-LAI)/LAICRT
          end where

          ZT = Z0(:,N)*ZT
          ZQ = ZT
!!!     end if
      elseif (ISTYPE==4) then
        ZT = Z0(:,N)*exp( -0.13 * sqrt(RE) )
      endif

      if (ISTYPE==1) then
        ZQ = ZT * 1.5
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        ZQ = ZT
      endif

!  Compute the the heat and moisture surface coefficients
!--------------------------------------------------------

      if (present(DCH) .and. present(DCQ)) then

         call GETCDH(TVA,UUU,DZ,TVS,ZT,ZQ,CN,RI,  CH(:,N),CQ(:,N),DCH(:,N),DCQ(:,N))

      else

         call GETCDH(TVA,UUU,DZ,TVS,ZT,ZQ,CN,RI,  CH(:,N),CQ(:,N))

      end if

      CH(:,N) = CH(:,N)*URA
      CM(:,N) = CM(:,N)*URA
      CQ(:,N) = CQ(:,N)*URA

      deallocate( UST )
      deallocate( TVS )
      deallocate( URA )
      deallocate( TVA )

   end subroutine louissurface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Z0SEA - Computes $z_o$ as a function of $u^*$ over water surfaces
! !INTERFACE:

   subroutine Z0SEA (Z0, USTAR, DZ)

! !ARGUMENTS:

     real,    intent(INOUT) :: Z0   (:)
     real,    intent(INOUT) :: USTAR(:)
     real,    intent(IN   ) :: DZ   (:)

! !DESCRIPTION:
!        Compute roughness length for ocean points
!          based on functions of Large and Pond
!          and of Kondo

     real               :: UF(size(Z0))
     integer, parameter :: NumIter = 3
     integer            :: K


! Begin

! UF is Karman*|U|*sqrt[f_m(Ri)] = |U|*sqrt[C_D]*log(DZ/Z0+1). 
! The stability factor is from the previous time step.
!------------------------------------------------------------------------------------

     UF    = USTAR * ALOG(DZ/Z0 + 1.0)

! Iterate roughness and U*
!-------------------------

     ITERATION: do K=1,NumIter

! Evaluate z0 using Beljaars' approx.
!------------------------------------

       Z0 = (0.11*MAPL_NUAIR)/USTAR + (0.018/MAPL_GRAV)*USTAR**2

! Update U* including stability effects
!--------------------------------------

       USTAR = min(max(UF/ALOG( DZ/Z0 + 1.0 ) , 1.e-6),5.0)

     end do ITERATION

   end subroutine Z0SEA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: helfsurface
! !INTERFACE:
   SUBROUTINE helfsurface(VUS,VVS,VT1,VT2,VSH1,VSH2,VP,VPE, &
        VZ0,LAI,IVWATER,VHS,N,IRUN, &
        VRHO,VKH,VKM,VUSTAR,VXX,VYY,VCU,VCT,VRIB,VZETA,VWS, &
        t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0,WMCHARNOCK)
!**********************************************************************
!  SUBROUTINE helfsurface - COMPUTES SURFACE TRANSFER COEFFICIENTS
!
!   ARGUMENTS ::
!
!     INPUT:
!     ------
!    US            -         U - COMPONENT OF SURFACE WIND
!    VS            -         V - COMPONENT OF SURFACE WIND
!    THV1          -         VIRTUAL POTENTIAL TEMPERATURE AT NLAY
!    THV2          -         VIRTUAL POTENTIAL TEMPERATURE AT GROUND
!    TH1           -         POTENTIAL TEMPERATURE AT NLAY
!    TH2           -         POTENTIAL TEMPERATURE AT GROUND
!    SH1           -         SPECIFIC HUMIDITY AT NLAY
!    SH2           -         SPECIFIC HUMIDITY AT GROUND
!    PK            -         EVEN LEVEL PRESSURE ** KAPPA AT LEVEL NLAY
!    PKE           -         EDGE LEVEL PRESSURE ** KAPPA AT GROUND
!    PE            -         SURFACE PRESSURE
!    Z0            -         SURFACE ROUGHNESS
!    WATER         -         ARRAY WITH '1' OVER OCEANS
!    HS            -         DEPTH OF SURFACE LAYER
!    N             -         NUMBER OF helfsurface ITERATIONS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!                                          4 - wave model Charnock coefficient
!     OUTPUT:
!     -------
!    RHO           -         DENSITY AT SURFACE
!    KH            -         HEAT TRANSFER COEFFICIENT (CT*USTAR)
!    KM            -         MOMENTUM TRANSFER COEFFICIENT (CU*USTAR)
!    USTAR         -         FRICTION VELOCITY
!    XX            -         PHIM(ZETA) - DIMENSIONLESS WIND SHEAR
!    YY            -         PHIH(ZETA) - DIMENSIONLESS TEMP GRADIENT
!    CU            -         MOMENTUM TRANSPORT COEFFICIENT
!    CT            -         HEAT TRANSPORT COEFFICIENT
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer, intent(in)   :: n,irun,CHOOSEZ0
      real,    intent(in)   :: VUS(:),VVS(:),VT1(:),VT2(:),VSH1(:),VSH2(:)
      real,    intent(in)   :: VPE(:),VP(:),LAI(:),VHS(:)
      real,    intent(inout):: VZ0(:)
      integer, intent(in)   :: IVWATER(:)
      real,    intent(out)  :: VRHO(:)
      real,    intent(out)  :: VKM(:),VKH(:),VUSTAR(:),VXX(:)
      real,    intent(out)  :: VYY(:),VCU(:),VCT(:),VRIB(:)
      real,    intent(out)  :: VZETA(:),VWS(:)
      real,    intent(OUT)  :: t2m(:),q2m(:),u2m(:),v2m(:)
      real,    intent(OUT)  :: t10m(:),q10m(:),u10m(:),v10m(:)
      real,    intent(OUT)  :: u50m(:),v50m(:)
      real, optional, intent(in) :: WMCHARNOCK(:)

! Local Variables
      integer IVBITRIB(irun)
      LOGICAL LWATER
      real VHZ(irun),VPSIM(irun),VAPSIM(irun),VPSIG(irun),VPSIHG(irun)
      real VTEMP(irun),VDZETA(irun),VDZ0(irun),VDPSIM(irun)
      real VDPSIH(irun),VZH(irun),VXX0(irun),VYY0(irun)
      real VAPSIHG(irun),VRIB1(irun)
      real VPSIH(irun),VPSIH2(irun),VH0(irun)
      real VX0PSIM(irun),VG(irun),VG0(irun),VR1MG0(irun)
      real VZ2(irun),VDZSEA(irun),VAZ0(irun),VXNUM1(irun)
      real VPSIGB2(irun),VDX(irun),VDXPSIM(irun),VDY(irun)
      real VXNUM2(irun),VDEN(irun),VAWS1(irun),VXNUM3(irun)
      real VXNUM(irun),VDZETA1(irun),VDZETA2(irun)
      real VZCOEF2(irun),VZCOEF1(irun),VTEMPLIN(irun)
      real VDPSIMC(irun),VDPSIHC(irun),VAHS(irun)
      real VTHV1(IRUN),VTHV2(IRUN),VTH1(IRUN),VTH2(IRUN),VPKE(IRUN),VPK(IRUN)

      real vz0h(irun),vh0h(irun),dummy1(irun),dummy2(irun),dummy3(irun),dummy4(irun),dummy5(irun), minVZ0(irun)

! Local Variables
      real USTMX3,USTZ0S,Z0MIN,H0BYZ0,USTH0S,H0VEG,Z0VEGM,PRFAC,Z0MAX
      real XPFAC,DIFSQT
      PARAMETER ( USTMX3 =   0.0632456)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( Z0MIN  =  USTZ0S/USTMX3)
      PARAMETER ( Z0MAX  =  USTZ0S/USTMX3)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )
      PARAMETER ( Z0VEGM =   0.005    )
      PARAMETER ( H0VEG  =  H0BYZ0*Z0VEGM )  !! This prevents discontinuity
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)

      real psihdiag(irun),psimdiag(irun)
      real rvk,vk2,bmdl(irun)
      integer itype
      integer iter
      logical call_psi

      real VCH(irun)

!
      if (present(WMCHARNOCK)) then
         VCH = WMCHARNOCK
      else
         VCH = 0.018
      end if

!
      _UNUSED_DUMMY(LAI)
      rvk = 1./MAPL_KARMAN
      vk2 = MAPL_KARMAN*MAPL_KARMAN
      where (ivwater == 3 ) 
         BMDL = 0.
!scale BMDL(i)    = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT) * exp(-lai(i)*2.)
      elsewhere
         BMDL = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT)
      endwhere

!     INITIALIZATION 

      VAHS = 1. / VHS
      VPKE = VPE ** MAPL_KAPPA
      VPK  = VP ** MAPL_KAPPA
      VTH1 = VT1/VPK
      VTH2 = VT2/VPKE
      VTHV1= VTH1*( 1.0 + MAPL_VIREPS*VSH1)
      VTHV2= VTH2*( 1.0 + MAPL_VIREPS*VSH2)

!     DETERMINE SURFACE WIND MAGNITUDE AND BULK RICHARDSON NUMBER
!
      VWS = max(VUS*VUS + VVS*VVS, 1.e-4)
      VRIB= MAPL_CP*(VPKE-VPK)*(VTHV1-VTHV2) / VWS
      VWS = SQRT( VWS )

!  INITIAL GUESS FOR ROUGHNESS LENGTH Z0 OVER WATER
!
      LWATER = any(IVWATER == 1)
!
      IF(LWATER)THEN
        where (IVWATER == 1) VZ0 = 0.0003
      ENDIF

      where(vz0 >= z0vegm) 
        vh0 = h0veg
      elsewhere
        vh0 = h0byz0 * vz0
      endwhere
      VZ0H = 0.001

!     CU AND PSIHG FOR NEUTRALLY STRATIFIED FLOW
!
      VHZ    = VHS / VZ0 + 1.
      VPSIM  = LOG( VHZ )
      VAPSIM = 1. / VPSIM
      VCU    = MAPL_KARMAN * VAPSIM
      VUSTAR = VCU * VWS
!
      VPSIG  = BMDL*sqrt(max(VH0*VUSTAR-USTH0S,0.))
      VPSIHG = VPSIM + VPSIG

!
!     LINEAR CORRECTION FOR ERROR IN ROUGHNESS LENGTH Z0
!
      IF(LWATER)THEN
        VTEMP = 0.
        CALL LINADJ(VRIB,VRIB,VWS,VWS,VZ0,VUSTAR,IVWATER,VAPSIM, &
                    VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,1,.TRUE.,IRUN,VDZETA, &
                    VDZ0,VDPSIM,VDPSIH,IVBITRIB, &
                    VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
                    VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
                    VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0,VCH)
        where ( IVWATER == 1) 
          VCU   = VCU * (1. - VDPSIM*VAPSIM)
          VZ0   = VZ0 + VDZ0
          VZ0   = max(VZ0, Z0MIN) 
          vh0   = h0byz0 * vz0
          VPSIG = VH0 * VCU * VWS - USTH0S
          VPSIG = max(VPSIG, 0.)
          VPSIG = SQRT( VPSIG  )
          VPSIG = BMDL * VPSIG
          VPSIHG= VPSIM + VDPSIH + VPSIG
        endwhere  
      ENDIF
!
!  INITIAL GUESS FOR STABILITY PARAMETER ZETA
!
      VZETA = VK2 * VRIB / (VCU * VCU * VPSIHG)
!
!  RECOMPUTE CU, ESTIMATE PSIHG AND UPDATE ZETA AND Z0
!
!      VZH(I) = VZ0(I) * VAHS(I)
      VZH = VZ0 / (VHS + VZ0)
      CALL PSI (VZETA,VZH,VPSIM,VTEMP,IRUN,VXX,VXX0,VYY,VYY0,2)
      VCU    = MAPL_KARMAN / VPSIM
      VPSIG  = VH0 * VCU * VWS - USTH0S
      VPSIG  = max(VPSIG, 0.)
      VPSIG  = SQRT(VPSIG)
      VPSIG  = BMDL * VPSIG
      VPSIHG = VPSIM + VPSIG
      VZETA  = VK2 * VRIB / (VCU * VCU * VPSIHG)
!
      IF(LWATER)THEN
        where (IVWATER.EQ.1) VUSTAR = VCU * VWS
        CALL ZCSUB ( VUSTAR,VCH,VHZ,IVWATER,.FALSE.,IRUN,VTEMP,CHOOSEZ0)
        CALL ZCSUB ( VUSTAR,VCH,VHZ,IVWATER,.FALSE.,IRUN,vz0h,2)
        where (IVWATER.EQ.1 )
          VZ0 = VTEMP
          VZ0   = max(VZ0,  Z0MIN) 
          VZ0H  = max(VZ0H, Z0MIN) 
          vh0  = h0byz0 * vz0
          vh0h = h0byz0 * vz0h
         endwhere
      ENDIF
!
!  ITERATIVE LOOP - N ITERATIONS
!     COMPUTE CU AND CT
!
      call_psi = (choosez0.eq.3 .AND. Lwater)
      ITYPE = 3
      DO ITER = 1,N

!       VZH(I) = VZ0(I) * VAHS(I)
        VZH = VZ0 / (VHS + VZ0)
        CALL PSI (VZETA,VZH,VPSIM,VPSIH,IRUN,VXX,VXX0,VYY,VYY0,1)
!       VZH(I) = VZ0H(I) * VAHS(I)
        VZH = VZ0H / (VHS + VZ0H)
        if ( call_psi ) CALL PSI (VZETA,VZH,dummy1,VPSIH,IRUN,dummy2,dummy3,dummy4,dummy5,3)

        VCU = MAPL_KARMAN / VPSIM
        VUSTAR = VCU * VWS
!
        VPSIG  = VH0 * VUSTAR - USTH0S
        VPSIG  = max(VPSIG, 0.)
        VPSIG  = SQRT(VPSIG)
        VPSIG  = BMDL * VPSIG
        VPSIHG = VPSIH + VPSIG
!
!  LINEAR CORRECTIONS FOR CU, CT, ZETA, AND Z0
!
        VAPSIM = VCU * RVK
        VAPSIHG= 1. / VPSIHG
        VRIB1  = VAPSIM * VAPSIM * VPSIHG * VZETA
!
        IF(ITER.EQ.N) ITYPE = 5
!
        CALL LINADJ(VRIB1,VRIB,VWS, &
                    VWS,VZ0,VUSTAR,IVWATER, &
                    VAPSIM,VAPSIHG,VPSIH, &
                    VPSIG,VXX,VXX0, &
                    VYY,VYY0,ITYPE,LWATER,IRUN,VDZETA, &
                    VDZ0,VDPSIM,VDPSIH, &
                    IVBITRIB, &
                    VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
                    VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
                    VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0,VCH)
!
!  UPDATES OF ZETA, Z0, CU AND CT
!
        VZETA = VZETA * ( 1. + VDZETA )
        where (IVBITRIB.EQ.1 ) VZETA = VPSIM * VPSIM * VRIB * VAPSIHG
!
        IF ( LWATER ) THEN
          where (IVWATER.EQ.1 ) 
            VZ0  = VZ0 * ( 1. + VDZ0)
            VZ0H = VZ0H * ( 1. + VDZ0 )
            VZ0   = max(VZ0,  Z0MIN) 
            VZ0H  = max(VZ0H, Z0MIN) 
            vh0  = h0byz0 * vz0
            vh0h = h0byz0 * vz0h
          endwhere
        ENDIF
!
        IF ( ITER .EQ. N ) THEN
          VPSIM  = VPSIM + VDPSIM
          VCU    = MAPL_KARMAN / VPSIM
          VUSTAR = VCU * VWS
!
          VPSIG  = VH0 * VUSTAR - USTH0S
          VPSIG  = max(VPSIG, 0.)
          VPSIG  = SQRT(VPSIG)
          VPSIG  = BMDL * VPSIG
          VPSIHG = VPSIH + VDPSIH + VPSIG
          VCT    = MAPL_KARMAN / VPSIHG
        ENDIF

!
!  SAVE VALUES OF RIB AND WS
!
        VRIB1 = VRIB
!
      enddo ! ITER 
!
!  CALCULATE RHO-SURFACE ( KG / M**3 )
!
      VTEMP =  10. * VAHS * VZETA
!       VZH(I) = VZ0(I) * 0.1
      VZH   = VZ0 / (10. + VZ0)

      CALL PSI (VTEMP,VZH,VHZ,VPSIH2,IRUN,VHZ,VHZ,VHZ,VHZ,3)

      VTEMP = min(( VPSIH2 + VPSIG ) / VPSIHG,1.)
      VRHO  = VPKE*( VTH2 + VTEMP * (VTH1-VTH2) )
      VRHO  = VPE *100. / ( MAPL_RGAS * VRHO )
!
! interpolate uvtq to 2, 10 and 50 meters for diagnostic output
!  use psih and psim which represent non-dim change from ground
!                 to specified level
! and multiply theta by surface p**kappa to get temperatures
!
      vtemp = 2. * vahs * vzeta
!        vzh(i) = min(vz0(i),2.) * 0.5
      minVZ0= min(VZ0,2.)
      VZH   = minvz0 / (2. + minVZ0)
      call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)

      vtemp = min(( psihdiag + vpsig ) / vpsihg, 1.)
      t2m   = (vth2 + vtemp* (vth1-vth2)) * vpke
      q2m   = vsh2 + vtemp* (vsh1-vsh2)
      vtemp = psimdiag/vpsim
      u2m   = vtemp * vus
      v2m   = vtemp * vvs

      vtemp = 10. * vahs * vzeta
!       vzh(i) = vz0(i) * 0.1
      VZH   = VZ0 / (10. + VZ0)
      call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)

      vtemp = min(( psihdiag + vpsig ) / vpsihg, 1.)
      t10m  = (vth2 + vtemp* (vth1-vth2))  * vpke
      q10m  = vsh2 + vtemp* (vsh1-vsh2)
      vtemp = psimdiag/vpsim
      u10m  = vtemp * vus
      v10m  = vtemp * vvs

      vtemp = 50. * vahs * vzeta
!       vzh(i) = vz0(i) * 0.02
      VZH   = VZ0 / (50. + VZ0)
      call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
      
      vtemp = psimdiag/vpsim
      u50m  = vtemp * vus
      v50m  = vtemp * vvs
!
!  EVALUATE TURBULENT TRANSFER COEFFICIENTS
!

!!     VKH(I) = VUSTAR(I) * VCT(I)
!!     VKM(I) = VUSTAR(I) * VCU(I)
      VKH = VUSTAR * VCT * VRHO
      VKM = VUSTAR * VCU * VRHO

      VRIB = MAPL_CP*(VPKE-VPK)*(VTHV1-VTHV2) / max(VUS*VUS + VVS*VVS,1.e-1)

   end subroutine helfsurface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: phi
! !INTERFACE:
   SUBROUTINE PHI(Z,PHIM,PHIH,IFLAG,N)
!**********************************************************************
!
!  FUNCTION PHI - SOLVES KEYPS EQUATIONS
!               - CALLED FROM PSI
!
!  DESCRIPTION OF PARAMETERS
!     Z     -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!               TIMES APPROPRIATE CONSTANT
!     PHIM  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR MOMENTUM
!     PHIH  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR SCALARS
!     IFLAG -  FLAG TO DETERMINE IF X IS NEEDED (IFLAG=2), Y IS NEEDED
!                  (IFLAG=3), OR BOTH (IFLAG=1)
!     N     -  LENGTH OF VECTOR TO BE SOLVED
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer, intent(in) :: n,iflag
      real,    intent(in) :: Z(:)
      real,    intent(out):: PHIM(:),PHIH(:)

! Local Variables
      integer I1(N),I2(N)
      real ZSTAR(N),E1(N),E2(N),TEMP1(N)
!
      real PHIM0(385),ZLINM1(75),ZLINM2(75),ZLINM3(36)
      real ZLOGM1(74),ZLOGM2(75),ZLOGM3(50)
      real PHIH0(385),ZLINH1(75),ZLINH2(75),ZLINH3(36)
      real ZLOGH1(74),ZLOGH2(75),ZLOGH3(50)
      EQUIVALENCE (PHIM0(1),ZLINM1(1)),(PHIM0(76),ZLINM2(1))
      EQUIVALENCE (PHIM0(151),ZLINM3(1))
      EQUIVALENCE (PHIM0(187),ZLOGM1(1)),(PHIM0(261),ZLOGM2(1))
      EQUIVALENCE (PHIM0(336),ZLOGM3(1))
      EQUIVALENCE (PHIH0(1),ZLINH1(1)),(PHIH0(76),ZLINH2(1))
      EQUIVALENCE (PHIH0(151),ZLINH3(1))
      EQUIVALENCE (PHIH0(187),ZLOGH1(1)),(PHIH0(261),ZLOGH2(1))
      EQUIVALENCE (PHIH0(336),ZLOGH3(1))
!
       DATA ZLOGM1/ &
                   0.697894,0.678839,0.659598,0.640260, &
        0.620910,0.601628,0.582486,0.563550,0.544877, &
        0.526519,0.508516,0.490903,0.473708,0.456951, &
        0.440649,0.424812,0.409446,0.394553,0.380133, &
        0.366182,0.352695,0.339664,0.327082,0.314938, &
        0.303222,0.291923,0.281029,0.270528,0.260409, &
        0.250659,0.241267,0.232221,0.223509,0.215119, &
        0.207041,0.199264,0.191776,0.184568,0.177628, &
        0.170949,0.164519,0.158331,0.152374,0.146641, &
        0.141123,0.135813,0.130702,0.125783,0.121048, &
        0.116492,0.112107,0.107887,0.103826,0.0999177, &
        0.0961563,0.0925364,0.0890528,0.0857003,0.0824739, &
        0.0793690,0.0763810,0.0735054,0.0707380,0.0680749, &
        0.0655120,0.0630455,0.0606720,0.0583877,0.0561895, &
        0.0540740,0.0520382,0.0500790,0.0481936,0.0463791/
       DATA ZLOGM2/ &
        0.0446330,0.0429526,0.0413355,0.0397792,0.0382816, &
        0.0368403,0.0354533,0.0341185,0.0328340,0.0315978, &
        0.0304081,0.0292633,0.0281616,0.0271013,0.0260809, &
        0.0250990,0.0241540,0.0232447,0.0223695,0.0215273, &
        0.0207168,0.0199369,0.0191862,0.0184639,0.0177687, &
        0.0170998,0.0164560,0.0158364,0.0152402,0.0146664, &
        0.0141142,0.0135828,0.0130714,0.0125793,0.0121057, &
        0.0116499,0.0112113,0.0107892,0.0103830,0.999210E-2, &
        0.961590E-2,0.925387E-2,0.890547E-2,0.857018E-2,0.824752E-2, &
        0.793701E-2,0.763818E-2,0.735061E-2,0.707386E-2,0.680754E-2, &
        0.655124E-2,0.630459E-2,0.606722E-2,0.583880E-2,0.561897E-2, &
        0.540742E-2,0.520383E-2,0.500791E-2,0.481937E-2,0.463792E-2, &
        0.446331E-2,0.429527E-2,0.413355E-2,0.397793E-2,0.382816E-2, &
        0.368403E-2,0.354533E-2,0.341185E-2,0.328340E-2,0.315978E-2, &
        0.304082E-2,0.292633E-2,0.281616E-2,0.271013E-2,0.260809E-2/
       DATA ZLOGM3/ &
        0.250990E-2,0.241541E-2,0.232447E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199369E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170998E-2,0.164560E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999210E-3, &
        0.961590E-3,0.925387E-3,0.890547E-3,0.857018E-3,0.824752E-3, &
        0.793701E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
       DATA ZLOGH1/ &
                   0.640529,0.623728,0.606937,0.590199, &
        0.573552,0.557032,0.540672,0.524504,0.508553, &
        0.492843,0.477397,0.462232,0.447365,0.432809, &
        0.418574,0.404670,0.391103,0.377878,0.364999, &
        0.352468,0.340284,0.328447,0.316954,0.305804, &
        0.294992,0.284514,0.274364,0.264538,0.255028, &
        0.245829,0.236933,0.228335,0.220026,0.211999, &
        0.204247,0.196762,0.189537,0.182564,0.175837, &
        0.169347,0.163088,0.157051,0.151231,0.145620, &
        0.140211,0.134998,0.129974,0.125133,0.120469, &
        0.115975,0.111645,0.107475,0.103458,0.995895E-1, &
        0.958635E-1,0.922753E-1,0.888199E-1,0.854925E-1,0.822886E-1, &
        0.792037E-1,0.762336E-1,0.733739E-1,0.706208E-1,0.679704E-1, &
        0.654188E-1,0.629625E-1,0.605979E-1,0.583217E-1,0.561306E-1, &
        0.540215E-1,0.519914E-1,0.500373E-1,0.481564E-1,0.463460E-1/
       DATA ZLOGH2/ &
        0.446034E-1,0.429263E-1,0.413120E-1,0.397583E-1,0.382629E-1, &
        0.368237E-1,0.354385E-1,0.341053E-1,0.328222E-1,0.315873E-1, &
        0.303988E-1,0.292550E-1,0.281541E-1,0.270947E-1,0.260750E-1, &
        0.250937E-1,0.241494E-1,0.232405E-1,0.223658E-1,0.215240E-1, &
        0.207139E-1,0.199342E-1,0.191839E-1,0.184618E-1,0.177669E-1, &
        0.170981E-1,0.164545E-1,0.158351E-1,0.152390E-1,0.146653E-1, &
        0.141133E-1,0.135820E-1,0.130707E-1,0.125786E-1,0.121051E-1, &
        0.116494E-1,0.112108E-1,0.107888E-1,0.103826E-1,0.999177E-2, &
        0.961561E-2,0.925360E-2,0.890523E-2,0.856997E-2,0.824733E-2, &
        0.793684E-2,0.763803E-2,0.735048E-2,0.707375E-2,0.680743E-2, &
        0.655114E-2,0.630450E-2,0.606715E-2,0.583873E-2,0.561891E-2, &
        0.540737E-2,0.520379E-2,0.500787E-2,0.481933E-2,0.463789E-2, &
        0.446328E-2,0.429524E-2,0.413353E-2,0.397790E-2,0.382814E-2, &
        0.368401E-2,0.354532E-2,0.341184E-2,0.328338E-2,0.315977E-2, &
        0.304081E-2,0.292632E-2,0.281615E-2,0.271012E-2,0.260809E-2/
       DATA ZLOGH3/ &
        0.250990E-2,0.241540E-2,0.232446E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199368E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170997E-2,0.164559E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999209E-3, &
        0.961590E-3,0.925387E-3,0.890546E-3,0.857018E-3,0.824752E-3, &
        0.793700E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
 
       DATA ZLINM1/ &
        0.964508,0.962277,0.960062,0.957863,0.955680, &
        0.953512,0.951359,0.949222,0.947100,0.944992, &
        0.942899,0.940821,0.938758,0.936709,0.934673, &
        0.932652,0.930645,0.928652,0.926672,0.924706, &
        0.922753,0.920813,0.918886,0.916973,0.915072, &
        0.913184,0.911308,0.909445,0.907594,0.905756, &
        0.903930,0.902115,0.900313,0.898522,0.896743, &
        0.894975,0.893219,0.891475,0.889741,0.888019, &
        0.886307,0.884607,0.882917,0.881238,0.879569, &
        0.877911,0.876264,0.874626,0.872999,0.871382, &
        0.869775,0.868178,0.866591,0.865013,0.863445, &
        0.861887,0.860338,0.858798,0.857268,0.855747, &
        0.854235,0.852732,0.851238,0.849753,0.848277, &
        0.846809,0.845350,0.843900,0.842458,0.841025, &
        0.839599,0.838182,0.836774,0.835373,0.833980/
       DATA ZLINM2/ &
        0.832596,0.831219,0.829850,0.828489,0.827136, &
        0.825790,0.824451,0.823121,0.821797,0.820481, &
        0.819173,0.817871,0.816577,0.815289,0.814009, &
        0.812736,0.811470,0.810210,0.808958,0.807712, &
        0.806473,0.805240,0.804015,0.802795,0.801582, &
        0.800376,0.799176,0.797982,0.796794,0.795613, &
        0.794438,0.793269,0.792106,0.790949,0.789798, &
        0.788652,0.787513,0.786380,0.785252,0.784130, &
        0.783014,0.781903,0.780798,0.779698,0.778604, &
        0.777516,0.776432,0.775354,0.774282,0.773215, &
        0.772153,0.771096,0.770044,0.768998,0.767956, &
        0.766920,0.765888,0.764862,0.763840,0.762824, &
        0.761812,0.760805,0.759803,0.758805,0.757813, &
        0.756824,0.755841,0.754862,0.753888,0.752918, &
        0.751953,0.750992,0.750035,0.749083,0.748136/
       DATA ZLINM3/ &
        0.747192,0.746253,0.745318,0.744388,0.743462, &
        0.742539,0.741621,0.740707,0.739798,0.738892, &
        0.737990,0.737092,0.736198,0.735308,0.734423, &
        0.733540,0.732662,0.731788,0.730917,0.730050, &
        0.729187,0.728328,0.727472,0.726620,0.725772, &
        0.724927,0.724086,0.723248,0.722414,0.721584, &
        0.720757,0.719933,0.719113,0.718296,0.717483, &
        0.716673/
       DATA ZLINH1/ &
        0.936397,0.932809,0.929287,0.925827,0.922429, &
        0.919089,0.915806,0.912579,0.909405,0.906284, &
        0.903212,0.900189,0.897214,0.894284,0.891399, &
        0.888558,0.885759,0.883001,0.880283,0.877603, &
        0.874962,0.872357,0.869788,0.867255,0.864755, &
        0.862288,0.859854,0.857452,0.855081,0.852739, &
        0.850427,0.848144,0.845889,0.843662,0.841461, &
        0.839287,0.837138,0.835014,0.832915,0.830841, &
        0.828789,0.826761,0.824755,0.822772,0.820810, &
        0.818869,0.816949,0.815050,0.813170,0.811310, &
        0.809470,0.807648,0.805845,0.804060,0.802293, &
        0.800543,0.798811,0.797095,0.795396,0.793714, &
        0.792047,0.790396,0.788761,0.787141,0.785535, &
        0.783945,0.782369,0.780807,0.779259,0.777724, &
        0.776204,0.774696,0.773202,0.771720,0.770251/
       DATA ZLINH2/ &
        0.768795,0.767351,0.765919,0.764499,0.763091, &
        0.761694,0.760309,0.758935,0.757571,0.756219, &
        0.754878,0.753547,0.752226,0.750916,0.749616, &
        0.748326,0.747045,0.745775,0.744514,0.743262, &
        0.742020,0.740787,0.739563,0.738348,0.737141, &
        0.735944,0.734755,0.733574,0.732402,0.731238, &
        0.730083,0.728935,0.727795,0.726664,0.725539, &
        0.724423,0.723314,0.722213,0.721119,0.720032, &
        0.718952,0.717880,0.716815,0.715756,0.714704, &
        0.713660,0.712621,0.711590,0.710565,0.709547, &
        0.708534,0.707529,0.706529,0.705536,0.704549, &
        0.703567,0.702592,0.701623,0.700660,0.699702, &
        0.698750,0.697804,0.696863,0.695928,0.694998, &
        0.694074,0.693155,0.692241,0.691333,0.690430, &
        0.689532,0.688639,0.687751,0.686868,0.685990/
       DATA ZLINH3/ &
        0.685117,0.684249,0.683386,0.682527,0.681673, &
        0.680824,0.679979,0.679139,0.678303,0.677472, &
        0.676645,0.675823,0.675005,0.674191,0.673381, &
        0.672576,0.671775,0.670978,0.670185,0.669396, &
        0.668611,0.667830,0.667054,0.666281,0.665512, &
        0.664746,0.663985,0.663227,0.662473,0.661723, &
        0.660977,0.660234,0.659495,0.658759,0.658027, &
        0.657298/

!
!
        where ( Z .GT.1.78e10 ) 
           ZSTAR = 384.9999
        elsewhere( Z .GT. 2. )
           TEMP1  = Z*0.5
           TEMP1 = LOG10(TEMP1)
           TEMP1 = (TEMP1 + 9.3) * 20.
           ZSTAR = TEMP1
        elsewhere
           ZSTAR    = 100. * Z - 14.
        endwhere
!
        I1 = ZSTAR
        I2 = I1 + 1
        TEMP1 = ZSTAR - I1
        ZSTAR = -Z
!
!

        IF ( IFLAG  <= 2 ) then
           where (z .ge. 0.15)
              E1 = PHIM0( I1 )
              E2 = PHIM0( I2 )
              PHIM  = TEMP1 * ( E2-E1 )
              PHIM  = PHIM +   E1
           elsewhere
              PHIM = 1. + ZSTAR &
                  *(0.25+ZSTAR*(0.09375+ZSTAR* &
                  (0.03125+0.00732422 * ZSTAR)))
           endwhere
        endif

        IF ( IFLAG  /=2 ) then
           where (z .ge. 0.15)
              E1 = PHIH0( I1 )
              E2 = PHIH0( I2 )
              PHIH  = TEMP1 * ( E2-E1 )
              PHIH  = PHIH +   E1
           elsewhere
              PHIH = 1.+ Z * (0.5+ZSTAR*(0.375+ZSTAR* &
                  (0.5+ZSTAR*(0.8203125+ZSTAR* &
                  (1.5+2.93262*ZSTAR)))))
              PHIH = 1. / PHIH
           endwhere
        endif
       
   end subroutine phi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: psi
! !INTERFACE:
   SUBROUTINE PSI(VZZ,VZH,VPSIM,VPSIH,IRUN,VX,VXS,VY,VYS,IFLAG)
!**********************************************************************
!
!  SUBROUTINE PSI - DETERMINES DIMENSIONLESS WIND AND
!                    SCALAR PROFILES IN SURFACE LAYER
!                 - CALLED FROM helfsurface
!
!  DESCRIPTION OF PARAMETERS
!     ZZ   -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!     ZH   -  INPUTED VALUE OF Z0 DIVIDED BY SFC LAYER HEIGHT
!     PSIM -  OUTPUTED VALUE OF DIMENSIONLESS WIND
!     PSIH -  OUTPUTED VALUE OF DIMENSIONLESS SCALAR
!     X    -  OUTPUTED VALUE OF PHIM(ZETA)
!     XS   -  OUTPUTED VALUE OF PHIM(ZETA0)
!     Y    -  OUTPUTED VALUE OF PHIH(ZETA)
!     YS   -  OUTPUTED VALUE OF PHIH(ZETA0)
!     IFLAG-  FLAG TO DETERMINE IF CU IS NEEDED (IFLAG=2),
!                  IF CT IS NEEDED (IFLAG=3), OR BOTH (IFLAG=1)
!  SUBPROGRAMS NEEDED
!     PHI  -  COMPUTES SIMILARITY FUNCTION FOR MOMENTUM AND SCALARS
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer, intent(in) :: irun,iflag
      real,    intent(in) :: VZZ(:),VZH(:)
      real,    intent(out) :: VPSIM(:),VPSIH(:),VX(:),VXS(:),VY(:),VYS(:)
 
! Local Variables
      real ZWM,RZWM,Z0M,ZCM,RZCM,CM1,CM2,CM6,CM7,CM8ARG,YCM
      PARAMETER ( ZWM     =    1.    )
      PARAMETER ( RZWM    =  1./ZWM  )
      PARAMETER ( Z0M     =    0.2    )
      PARAMETER ( ZCM     =    42.    )
      PARAMETER ( RZCM    =  1./ZCM  )
      PARAMETER ( CM1     =  1./126. )
      PARAMETER ( CM2     =  1./(6.*CM1)  )
      PARAMETER ( CM6     =  6. / ( 1. + 6.*CM1 )  )
      PARAMETER ( CM7     =  CM2 + ZWM  )
      PARAMETER ( CM8ARG  =  CM7*ZCM*RZWM / (CM2+ZCM)  )
      PARAMETER ( YCM     =  6. / ( 1. + 6.*CM1*ZCM )  )

      integer INTSTB(irun),INTZ0(irun)
      real, allocatable :: ZZ0(:),Z(:),Z2(:),Z1(:),Z0(:)
      real, allocatable :: X0(:),X1(:),Y0(:),Y1(:)
      real, allocatable :: PSI2(:),TEMP(:)
      real, allocatable :: HZ(:),ARG0(:),ARG1(:),DX(:)
      real, allocatable :: X0NUM(:),X1NUM(:),X0DEN(:)
      real, allocatable :: X1DEN(:),Y1DEN(:),Z2ZWM(:)
      real cm3,cm4,cm5,cm8
      integer, allocatable :: indxs(:)
      integer ibit
      integer k
!
      CM3 =   sqrt( 0.2/CM1-0.01 )
      CM4 =   1./CM3
      CM5 =  (10.-CM1) / (10.*CM1*CM3)
      CM8 =   6. * LOG(CM8ARG)
!
      VPSIM = 0.
      VPSIH = 0.
      VX    = 0.
      VXS   = 0.
      VY    = 0.
      VYS   = 0.
      ZZ0   = VZH*VZZ

      indxs = pack([(k,k=1,irun)], VZZ .LE. -1.e-7)
      IBIT  = size(indxs) 
      if (IBIT == 0) return
!
! ****************************************
! *****    UNSTABLE SURFACE LAYER    *****
! ****************************************
!
      Z  = VZZ(indxs)
      Z0 = ZZ0(indxs)
      Z  = -18. * Z
      Z0 = -18. * Z0
      allocate(X1(IBIT), Y1(IBIT), X0(IBIT), Y0(IBIT)) 

      CALL PHI( Z,X1,Y1,IFLAG,IBIT )
      CALL PHI( Z0,X0,Y0,IFLAG,IBIT )
 
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF(IFLAG <3 ) then
!
         ARG1 = 1. - X1
         where( Z .LT. 0.013 ) 
           ARG1 = Z * ( 0.25 -  0.09375 * Z )
         endwhere
!         
         ARG0  = 1. - X0
         where( Z0 .LT. 0.013 )
            ARG0 = Z0 * ( 0.25 -  0.09375 * Z0 )
         endwhere
! 
         ARG1 = ARG1 * ( 1.+X0 )
         ARG0 = ARG0 * ( 1.+X1 )
         DX   = X1 - X0
         ARG1 = ARG1 / ARG0
         ARG0 = -DX / ( 1. + X1*X0 )
         ARG0 = ATAN( ARG0 )
         ARG1 = LOG( ARG1 )
         PSI2 = 2. * ARG0 + ARG1
         PSI2 = PSI2 + DX
!
         VPSIM(indxs) = PSI2
         VX(indxs)    = X1
         VXS(indxs)   = X0
      endif ! iflag<3
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
!

      IF (IFLAG /=2) then

        ARG1 = 1. - Y1
        where( Z .LT. 0.0065 )
           ARG1 = Z * ( 0.5 -  0.625 * Z )
        endwhere
!
        ARG0  = 1. - Y0
        where( Z0 .LT. 0.0065 )
           ARG0 = Z0 * ( 0.5 -  0.625 * Z0 )
        endwhere
! 
        ARG1 = ARG1 * ( 1. + Y0 )
        ARG0 = ARG0 * ( 1. + Y1 )
        ARG1 = ARG1 / ARG0
        PSI2 = LOG( ARG1 )
        PSI2 = PSI2 - Y1 + Y0
!
        VPSIH(indxs) = PSI2
        VY(indxs)    = Y1
        VYS(indxs)   = Y0
      ENDIF
!
! **************************************
! *****    STABLE SURFACE LAYER    *****
! **************************************
!
      Z   = VZZ(indxs)
      Z0  = ZZ0(indxs)
      ARG1= VZH(indxs)

      HZ = 1. / ARG1
      Z1 = Z
      Z2 = Z ! just for allocatation purpose
      Z2 = ZWM
!
      where( Z .GT. ZWM ) 
        Z1 = ZWM
        Z2 = Z
      endwhere
!
      Z0 = min(Z0, Z0M ) 
!
      X1NUM = 1. + 5. * Z1
      X0NUM = 1. + 5. * Z0
      X1DEN = 1. / (1. + CM1 * (X1NUM * Z1) )
      X0DEN = 1. + CM1 * (X0NUM * Z0)
!
      where( Z0 .GT. Z0M .OR. Z.GT.ZWM ) HZ = Z1 / Z0
      ARG1 = HZ*HZ*X0DEN*X1DEN
      ARG1 = LOG( ARG1 )
      ARG1 = 0.5 * ARG1
      ARG0 = (Z1 + 0.1) * (Z0 + 0.1)
      ARG0 = CM3 + ARG0 * CM4
      ARG0 = ( Z1 - Z0 ) / ARG0
      ARG0 = ATAN( ARG0 )
      TEMP = ARG1 + CM5 * ARG0
!
      where( Z0 .GT. Z0M) 
        X0 = 0.
      elsewhere
        X0 = X0NUM / X0DEN
      endwhere
  
      Z2ZWM = Z2 * RZWM
!
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF ( IFLAG < 3 ) then
        X1   = X1NUM * X1DEN
        ARG1 = LOG( Z2ZWM )
        PSI2 = TEMP + CM6 * ARG1
        VPSIM(indxs) = PSI2
        VX(indxs)    = X1
        VXS(indxs)   = X0
      endif !i flag<3)
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
!
      if (IFLAG /=2) then
        Y1DEN = 1. + CM1 * ( X1NUM * Z )
        Y1    = X1NUM / Y1DEN
        ARG1  = CM7 * Z2ZWM / ( CM2 + Z2 )
        ARG0  = 6.
        where ( Z2 .GT. ZCM ) 
           Y1   = YCM
           ARG1 = Z2 * RZCM
           ARG0 = YCM
           TEMP = TEMP + CM8
        ENDwhere
        ARG1 = LOG( ARG1 )
        PSI2 = TEMP + ARG0 * ARG1

        VPSIH(indxs) = PSI2
        VY(indxs)    = Y1
        VYS(indxs)   = X0
      endif ! iflag /=2
!
   end subroutine psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: linadj
! !INTERFACE:
   SUBROUTINE LINADJ ( VRIB1,VRIB2,VWS1,VWS2,VZ1,VUSTAR,IWATER, &
       VAPSIM, VAPSIHG,VPSIH,VPSIG,VX,VX0,VY,VY0,ITYPE,LWATER,IRUN, &
       VDZETA,VDZ0,VDPSIM,VDPSIH,INTRIB, &
       VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
       VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
       VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,vk,bmdl,CHOOSEZ0,VCHARNOCK)
!
!**********************************************************************
!
!  ARGUMENTS ::
!
!     INPUT:
!     ------
!    RIB1          -         BULK RICHARDSON NUMBER OF INPUT STATE
!    RIB2          -         DESIRED BULK RICH NUMBER OF OUTPUT STATE
!    WS1           -         SURFACE WIND SPEED OF INPUT STATE
!    WS2           -         DESIRED SURFACE WIND SPEED OF OUTPUT STATE
!    Z1            -         INPUT VALUE OF ROUGHNESS HEIGHT
!    USTAR         -         INPUT VALUE OF CU * WS
!    WATER         -         BIT ARRAY - '1' WHERE OCEAN
!    APSIM         -         (1/PSIM)
!    APSIHG        -         ( 1 / (PSIH+PSIG) )
!    PSIH          -         NON-DIM TEMP GRADIENT
!    PSIG          -         PSIH FOR THE MOLECULAR LAYER
!    X             -         PHIM(ZETA) - DERIVATIVE OF PSIM
!    X0            -         PHIM(ZETA0)
!    Y             -         PHIH(ZETA) - DERIVATIVE OF PSIH
!    Y0            -         PHIH(ZETA0)
!    ITYPE         -         INTEGER FLAG :
!                               1    = NEUTRAL ADJUSTMENT
!                               3, 5 = ADJUSTMENT INSIDE LOOP
!                               5    = ADJUST CU AND CT
!    LWATER        -         LOGICAL - .TRUE. IF THERE ARE WATER POINTS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!                                          4 - wave model Charnock coefficient
!
!     OUTPUT:
!     -------
!    DZETA         -         D LOG ZETA
!    DZ0           -         D Z0 (ITYPE 1) OR D LOG Z0 (ITYPE 2-5)
!    DPSIM         -         D PSIM
!    DPSIH         -         D PSIH
!    BITRIB        -         BIT ARRAY - '1' WHERE RIB1 = 0
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer, intent(in)   :: irun,itype,CHOOSEZ0
      real,    intent(in)   :: VRIB1(:),VRIB2(:)
      real,    intent(in)   :: VWS1(:),VWS2(:),VZ1(:)
      real,    intent(inout):: VUSTAR(:)
      integer, intent(in)   :: IWATER(:)
      real,    intent(in)   :: VAPSIM(:),VAPSIHG(:)
      real,    intent(in)   :: VPSIH(:),VPSIG(:),VX(:)
      real,    intent(in)   :: VX0(:),VY(:),VY0(:)
      LOGICAL, intent(in)   :: LWATER
      real,    intent(out)  :: VDZETA(:),VDZ0(:),VDPSIM(:)
      real,    intent(out)  :: VDPSIH(:)
      integer, intent(out)  :: INTRIB(:)
      real,    intent(out)  :: VX0PSIM(:),VG(:),VG0(:),VR1MG0(:)
      real,    intent(out)  :: VZ2(:),VDZSEA(:),VAZ0(:)
      real,    intent(out)  :: VPSIGB2(:),VDX(:),VDXPSIM(:),VDY(:)
      real,    intent(out)  :: VXNUM1(:), VXNUM2(:),VDEN(:),VAWS1(:),VXNUM3(:)
      real,    intent(out)  :: VXNUM(:),VDZETA1(:),VDZETA2(:)
      real,    intent(out)  :: VZCOEF2(:),VZCOEF1(:),VTEMPLIN(:)
      real,    intent(out)  :: VDPSIMC(:),VDPSIHC(:)
      real,    intent(in)   :: bmdl(:)
      real,    intent(in)   :: VCHARNOCK(:)

! Local Variables
      real xx0max,prfac,xpfac,difsqt,ustz0s,h0byz0,usth0s
      PARAMETER ( XX0MAX  =   1.49821 )
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )

      integer VINT1(irun),VINT2(irun)
      real vk,b2uhs(irun)
!
      _UNUSED_DUMMY(VWS2)
      _UNUSED_DUMMY(vk)

      B2UHS   = BMDL * BMDL * USTH0S

!   COMPUTE X0/PSIM, 1/Z0, G, G0, 1/(1-G0),
!     DEL LOG Z0, D LOG ZO / D USTAR
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
        where (IWATER.EQ.1) VX0PSIM = VAPSIM
      ENDIF
      IF ( ITYPE .GE. 3 ) THEN
        VX0PSIM = VX0 * VAPSIM
      ENDIF
!
      VDZ0   = 0.
      VG     = 0.
      VG0    = 0.
      VR1MG0 = 1.
!
      IF ( LWATER ) THEN
        CALL ZCSUB ( VUSTAR,VCHARNOCK,VDZSEA,IWATER,.TRUE.,IRUN,VZ2,CHOOSEZ0)

        VDZSEA = min( VDZSEA, 0.2*VZ1/VAPSIM ) ! To prevent Divide by Zero as VG0 => 1.0
!
        where ( IWATER.EQ.1) 
          VAZ0   = 1. / VZ1
          VG     = VDZSEA  * VAZ0
          VG0    = VX0PSIM * VG
          VR1MG0 = 1. / ( 1. - VG0 )
          VDZ0   = ( VZ2 - VZ1 ) * VR1MG0
        ENDwhere
        IF (ITYPE.GE.3) where ( IWATER.EQ.1) VDZ0 = VDZ0 * VAZ0 
      ENDIF
!
!
!   COMPUTE NUM1,NUM2,NUM3, DEN
!
      IF (ITYPE.GE.3) THEN
        where(VRIB1.EQ.0.) 
          INTRIB = 1
          VXNUM1 = 0.
        ELSEwhere
          INTRIB = 0
          VXNUM1 = 1. / VRIB1
        ENDwhere

        VPSIGB2 = 0.
        where(vpsig.gt.0.) VPSIGB2 = &
              0.5 * ( vpsig*vpsig + b2uhs ) / vpsig
        VDX      = VX - VX0
        VDXPSIM  = VDX * VAPSIM
        VDY      = VY - VY0
        VXNUM3   = - VPSIGB2
!
        IF ( LWATER ) THEN
          where (IWATER.EQ.1) 
             VDXPSIM = VDXPSIM * VR1MG0
             VXNUM3  = VXNUM3 + VG * ( VY0 - VPSIGB2)
             VXNUM2  = VY0 - VPSIGB2 - VX0PSIM * VPSIGB2
             VXNUM2  = (VXNUM2 * VAPSIHG) - 2. * VX0PSIM
             VXNUM2  = VXNUM2 * VDZ0
          endwhere
        ENDIF
!
        VDEN = VDY + VDXPSIM * VXNUM3
        VDEN = ( 1. + VDEN * VAPSIHG ) - 2. * VDXPSIM
      ENDIF
!
      IF (ITYPE.EQ.5) THEN
        VAWS1  = VR1MG0 / VWS1
        VXNUM3 = VXNUM3 * VAPSIHG
!
        IF ( LWATER ) THEN
          where(IWATER.EQ.1) 
             VXNUM3 = VXNUM3 - 2. * VG0
             VXNUM3 = VAWS1 * VXNUM3
          endwhere
        ENDIF
      ENDIF
!
!   COMPUTE D LOG ZETA
!
      IF (ITYPE.GE.3) THEN
        VXNUM = VRIB2 - VRIB1
        where( (VX0.GT.XX0MAX).AND.(VXNUM.GE.0.) )VXNUM = 0.
        VXNUM = VXNUM1 * VXNUM
!
        VDZETA1 = VXNUM
        IF(LWATER) then
          where(IWATER.EQ.1) VXNUM = VXNUM + VXNUM2
        endif

        VDEN = max(VDEN, 0.1)
!
        VDZETA = VXNUM / VDEN
        where((VRIB2.EQ.0.).OR.(VDZETA.LE.-1.)) VDZETA = VDZETA1
      ENDIF
!
!   COMPUTE D LOG Z0
!
      IF ( LWATER .AND. (ITYPE.GE.3) )THEN
        where ( IWATER.EQ.1 ) 
         VZCOEF2 = VG * VDXPSIM
         VDZ0    = VDZ0 - VZCOEF2 * VDZETA
        ENDwhere
      ENDIF
!
      IF ( LWATER .AND. (ITYPE.EQ.5) ) THEN
        where(IWATER.EQ.1) VZCOEF1 = VG * VAWS1
      ENDIF
!
!   CALCULATE D PSIM AND D PSIH
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
        where (IWATER.EQ.1) 
         VDPSIM = - VDZ0 * VAZ0
         VDPSIH = VDPSIM
        ENDwhere
      ENDIF
!
      IF (ITYPE.GE.3) THEN
        VDPSIM = VDX * VDZETA
        VDPSIH = VDY * VDZETA
        IF ( LWATER ) THEN
          where (IWATER.EQ.1 ) 
             VDPSIM = VDPSIM - VX0 * VDZ0
             VDPSIH = VDPSIH - VY0 * VDZ0
          ENDwhere
        ENDIF
      ENDIF
!
!   PREVENT OVERCORRECTION OF PSIM OR PSIH FOR UNSTABLE CASE
!
      IF (ITYPE.GE.4) THEN
        VDPSIMC = -0.9 - VDPSIM * VAPSIM
        VDPSIHC = -0.9 *  VPSIH - VDPSIH
        where( VDPSIMC.GT.0.  ) 
          VINT1 = 1
        ELSEwhere
          VINT1 = 0
        ENDwhere

        where ( VDPSIHC.GT.0.  )
          VINT2 = 1
        ELSEwhere
          VINT2 = 0
        ENDwhere

        VDZETA1 = 0.
        where(VINT1.EQ.1) VDZETA1 = VDPSIMC / VDXPSIM
        where((VINT1.EQ.1).OR.(VINT2.EQ.1)) VTEMPLIN = &
              VDY + VY0 * VG * VDXPSIM
!AMM    IF (VINT2(I).EQ.1 .and. VTEMPLIN(I).GT.tiny(1.0)) then
        where (VINT2.EQ.1) VDZETA2 =  VDPSIHC / VTEMPLIN
        where (VDZETA2.LT.VDZETA1 ) VDZETA1 = VDZETA2
        where ((VINT1.EQ.1).OR.(VINT2.EQ.1)) 
          VDZETA = VDZETA1 + VDZETA
          VDPSIM = VDPSIM + VDX * VR1MG0 * VDZETA1
          VDPSIH = VDPSIH + VTEMPLIN * VDZETA1
          where ( IWATER.EQ.1 ) &
            VDZ0 = VDZ0 - VG * VDXPSIM * VDZETA1
        ENDwhere
      ENDIF
!
   end subroutine linadj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: zcsub
! !INTERFACE:
   SUBROUTINE ZCSUB (VUSTAR,VCHARNOCK,VDZSEA,IWATER,LDZSEA,IRUN,VZSEA,CHOOSEZ0)
!**********************************************************************
!  FUNCTION ZSEA
!  PURPOSE
!     COMPUTES Z0 AS A FUNCTION OF USTAR OVER WATER SURFACES
!  USAGE
!     CALLED BY helfsurface
!  DESCRIPTION OF PARAMETERS
!     USTAR    -  INPUTED VALUE OF SURFACE-STRESS VELOCITY
!     DZSEA    -  OUTPUTED VALUE OF DERIVATIVE  D(ZSEA)/D(USTAR)
!     WATER    -  INPUTED BIT VECTOR TO DETERMINE WATER POINTS
!     LDZSEA   -  LOGICAL FLAG TO DETERMINE IF DZSEA SHOULD BE COMPUTED
!     ZSEA     -  OUTPUTED VALUE OF ROUGHNESS LENGTH
!     CHOOSEZ0 -  INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                               1 - Edson Z0 for mom. and heat, high wind limit
!                               2 - L&P Z0, high wind limit
!                               3 - Edson Z0 for mom. only, high wind limit
!                               4 - wave-model Charnock coefficient
!  SUBPROGRAMS NEEDED
!     NONE
!  RECORD OF MODIFICATIONS
!   Molod 6/8/2011 - Implement new choozez0 options (expand from 0,1 choice)
!  REMARKS:
!        COMPUTE ROUGHNESS LENGTH FOR OCEAN POINTS
!          BASED ON FUNCTIONS OF LARGE AND POND
!          AND OF KONDO --- DESIGNED FOR K = .4
! *********************************************************************
      implicit none 

! Argument List Delcarations
      integer, intent(in)    :: irun, CHOOSEZ0
      real,    intent(in)    :: VCHARNOCK(:)
      real,    intent(inout) :: VUSTAR(:)
      real,    intent(out)   :: VZSEA(:),VDZSEA(:)
      integer, intent(in)    :: IWATER(:)
      LOGICAL, intent(in)    :: LDZSEA

! Local Variables
      real USTMX1_OLD,USTMX2_OLD
      real USTMX1_NEW,USTMX2_NEW
      real USTMX1,USTMX2,USTMX3

      PARAMETER ( USTMX1_NEW =   0.80 )
      PARAMETER ( USTMX2_NEW =   0.80 )
      PARAMETER ( USTMX1_OLD =   1.1  )
      PARAMETER ( USTMX2_OLD =   0.381844 )
      PARAMETER ( USTMX3     =   0.0632456)

      real AA(IRUN,5),TEMP(IRUN)
      integer INT2(IRUN),INT3(IRUN),INT4(IRUN)
      integer k
      real ustloc(irun)

      real AA1(5),AA2(5),AA3(5),AA4(5)
      real AA2_NEW(5),AA3_NEW(5),AA4_NEW(5)
      real AA2_OLD(5),AA3_OLD(5),AA4_OLD(5)

      DATA AA1/.2030325E-5,0.0,0.0,0.0,0.0/

      DATA AA2_NEW/-1.102451E-08,0.1593E-04,0.1E-03,2.918E-03, &
               0.695649E-04/
      DATA AA3_NEW/-1.102451E-08,0.12E-04,0.1E-03,2.918E-03, &
               1.5649E-04/
      DATA AA4_NEW/0.085E-03,1.5E-03,-0.210E-03,0.215E-02, &
               -0.0/

      DATA AA2_OLD/-0.402451E-08,0.239597E-04,0.117484E-03,0.191918E-03, &
               0.395649E-04/
      DATA AA3_OLD/-0.237910E-04,0.228221E-03,-0.860810E-03,0.176543E-02, &
               0.784260E-04/
      DATA AA4_OLD/-0.343228E-04,0.552305E-03,-0.167541E-02,0.250208E-02, &
               -0.153259E-03/

      ustloc = max(1.e-6, vustar)
      CHARNOCK: if ( CHOOSEZ0 == 4 ) then
          VZSEA = (0.11*MAPL_NUAIR)/ustloc + (VCHARNOCK/MAPL_GRAV)*ustloc**2
          
          DERIVATIVE: if ( LDZSEA ) then
              VDZSEA = -(0.11*MAPL_NUAIR)/ustloc**2 + (VCHARNOCK/MAPL_GRAV)*2*ustloc
          end if DERIVATIVE

          return
      end if CHARNOCK

      if( CHOOSEZ0.eq.0 .OR. CHOOSEZ0.eq.2) then
          USTMX1 = USTMX1_OLD
          USTMX2 = USTMX2_OLD
             AA2 =    AA2_OLD
             AA3 =    AA3_OLD
             AA4 =    AA4_OLD
      else
          USTMX1 = USTMX1_NEW
          USTMX2 = USTMX2_NEW
             AA2 =    AA2_NEW
             AA3 =    AA3_NEW
             AA4 =    AA4_NEW
      endif
!
!**********************************************************************
!*****              LOWER CUTOFF CONDITION FOR USTAR                ***
!**********************************************************************
!
      vustar = ustloc
!
!***********************************
!*****  LOAD THE ARRAY A(I,K)  *****
!***********************************
!

      where( (ustloc .GT. USTMX1) .AND. (IWATER .EQ.1) )
        AA(:,1) = AA4(1)
        AA(:,2) = AA4(2)
        AA(:,3) = AA4(3)
        AA(:,4) = AA4(4)
        AA(:,5) = AA4(5)
      elsewhere(ustloc .GT. USTMX2)
        AA(:,1) = AA3(1)
        AA(:,2) = AA3(2)
        AA(:,3) = AA3(3)
        AA(:,4) = AA3(4)
        AA(:,5) = AA3(5)
      elsewhere(ustloc .GE. USTMX3)
        AA(:,1) = AA2(1)
        AA(:,2) = AA2(2)
        AA(:,3) = AA2(3)
        AA(:,4) = AA2(4)
        AA(:,5) = AA2(5)
      elsewhere
        AA(:,1) = AA1(1)
        AA(:,2) = AA1(2)
        AA(:,3) = AA1(3)
        AA(:,4) = AA1(4)
        AA(:,5) = AA1(5)
      endwhere
      if( CHOOSEZ0.gt.0 ) where( (ustloc .GT. USTMX1) .AND. (IWATER .EQ.1) ) ustloc = ustmx1

!
!********************************************************
!*****  EVALUATE THE ENHANCED POLYNOMIAL FOR ZSEA  *****
!********************************************************
!
      VDZSEA  =  ( AA(:,4) + AA(:,5) * ustloc ) * ustloc
      VZSEA   =  AA(:,2) + ( AA(:,3) + VDZSEA ) * ustloc
      TEMP    =  AA(:,1) / ustloc
      VZSEA   =  VZSEA + TEMP
!
!**********************************************************************
!*****        EVALUATE THE DERIVATIVE DZSEA IF LDZSEA IS TRUE       ***
!**********************************************************************
!
      IF( LDZSEA ) THEN
        VDZSEA  =  3. * VDZSEA -(AA(:,4)*ustloc - AA(:,3))
        VDZSEA  =  VDZSEA * ustloc - TEMP
      ENDIF

!
   end subroutine zcsub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sfclayer
