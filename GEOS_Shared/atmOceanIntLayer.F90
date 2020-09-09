module atmOcnIntlayer

! !USES:

use MAPL

implicit none
private SIMPLE_SW_ABS, AOIL_SST
public  ALBSEA, COOL_SKIN

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_SST - Computes skin SST using AOIL (single column)

!  !DESCRIPTION:
!        Based on Akella and Suarez, 2018
!        "The Atmosphere-Ocean Interface Layer of the NASA
!         Goddard Earth Observing System Model and Data Assimilation System."
!         GMAO Tech Memo, Vol 51.

! !INTERFACE:
  subroutine AOIL_SST  ( DO_DATASEA, DT, epsilon_d, F_PHI, STOKES_SPEED, &
             AOIL_depth, MUSKIN, SWN, PEN, PEN_ocean, LHF, SHF, LWDNSRF, &
             ALW, BLW, USTARW, TDROP, TS_FOUNDi, SWWARM, QWARM, ZETA_W,  &
             PHIW, LANGM, TAUTW, DELTC, TBAR, TDEL, DTS, TS, TWMTF)

! !ARGUMENTS:

    integer, intent(IN) :: DO_DATASEA  ! =0: coupled with ocean model (CGCM)
                                       ! =1: uncoupled                (AGCM)

    real, intent(IN)    :: DT          ! model time step
    real, intent(IN)    :: epsilon_d   ! ratio: (depth of AOIL)/(ocean model top level)

    real, intent(IN)    :: F_PHI       ! a scaler (tunable parameter) used to compute stability function
    real, intent(IN)    :: STOKES_SPEED! should be input from wave model, dummy for now

    real, intent(IN)    :: AOIL_depth  ! depth of the AOIL
    real, intent(IN)    :: MUSKIN      ! exponent in the prescribed temperature profile within AOIL

    real, intent(IN)    :: SWN         ! net shortwave radiation at surface
    real, intent(IN)    :: PEN         ! net shortwave radiation below the AOIL
    real, intent(IN)    :: PEN_ocean   ! net shortwave radiation below the AOIL with ocean model
    real, intent(IN)    :: LHF         ! latent   heat flux
    real, intent(IN)    :: SHF         ! sensible heat flux
    real, intent(IN)    :: LWDNSRF     ! downward longwave radiation
    real, intent(IN)    :: ALW         ! upward   longwave = ALW + BLW * TS
    real, intent(IN)    :: BLW         ! upward   longwave = ALW + BLW * TS

    real, intent(IN)    :: USTARW      ! friction velocity over water
    real, intent(IN)    :: TDROP       ! temperature drop due to cool skin layer
    real, intent(IN)    :: TS_FOUNDi   ! temperature at base of the AOIL

    real, intent(OUT)   :: SWWARM      ! net shortwave radiation absorbed within the AOIL
    real, intent(OUT)   :: QWARM       ! net heat flux within the AOIL
    real, intent(OUT)   :: ZETA_W      ! similarity parameter = AOIL_depth/(MO length scale)
    real, intent(OUT)   :: PHIW        ! stability function
    real, intent(OUT)   :: LANGM       ! Langmuir `number', dummy for now
    real, intent(OUT)   :: TAUTW       ! time scale at which diurnal warming -> 0, or TWMTF -> 0.
    real, intent(OUT)   :: DELTC       ! temperature drop due to cool skin layer = (above) TDROP
    real, intent(OUT)   :: TBAR        ! depth averaged mean AOIL temperature
    real, intent(OUT)   :: TDEL        ! temperature at top of warm layer (within the AOIL)
    real, intent(OUT)   :: DTS         ! temperature change: next time step - previous

    real, intent(INOUT) :: TS          ! skin SST
    real, intent(INOUT) :: TWMTF       ! AOIL state variable

!  !LOCAL VARIABLES

    real         :: ALPH


         ALPH   = (0.6 + 0.0935*(TS-MAPL_TICE))*1.E-4

         SWWARM = SWN - PEN
         QWARM  = SWWARM - (LHF + SHF - (LWDNSRF - ALW - BLW*TS))

         ZETA_W = (AOIL_depth*MAPL_KARMAN/USTARW**3.)*MAPL_GRAV*ALPH*QWARM/(MAPL_RHO_SEAWATER*MAPL_CAPWTR)

         PHIW  = 1.+SQRT(1.+4.*MAPL_KARMAN**2.*(1.+MUSKIN)*F_PHI*AOIL_depth*MAPL_GRAV*ALPH*TWMTF/USTARW**2.)
         PHIW  = 0.5*PHIW

         LANGM = SQRT( USTARW/STOKES_SPEED)

         TAUTW = (AOIL_depth*PHIW)/(MAPL_KARMAN*USTARW*(1.+MUSKIN))

         if (DO_DATASEA == 0) then ! with an OGCM, as in coupled GCM
           QWARM = QWARM - (epsilon_d/(1.-epsilon_d))* (PEN-PEN_ocean)
           TAUTW = (1.- epsilon_d) * TAUTW   ! compare \tau_{\sigma} in Eq.(22) and that in section 2.2 of AS2018
         endif

         TWMTF = (TWMTF+DT*(QWARM/(AOIL_depth*MAPL_RHO_SEAWATER*MAPL_CAPWTR)))/(1.+DT/TAUTW)
         TWMTF = max( TWMTF, 0.)

         DELTC = TDROP

         TBAR  = TWMTF + TS_FOUNDi

         if (DO_DATASEA == 1) then ! Atmospheric GCM, ocean surface is from "data"
          TDEL = TS_FOUNDi + ((1.+MUSKIN)/MUSKIN) * TWMTF
         else                      ! Coupled GCM
          TDEL = TS_FOUNDi + (1./MUSKIN + (1.-epsilon_d)) * TWMTF
         endif

         DTS = (TDEL- TDROP) - TS
         TS  = TDEL - TDROP       ! updated skin temperature

  end subroutine AOIL_SST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: COOL_SKIN - Computes variables related to the Cool Skin Layer

!  !DESCRIPTION:
!        Based on Fairall et al, 1996
!        And described in section 2.1 of
!        https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.2988

! !INTERFACE:
  subroutine COOL_SKIN (NT,CM,UUA,VVA,UW,VW,SWN,LHF,SHF,LWDNSRF,    &
                        ALW,BLW,TXW,TYW,USTARW_,                    &
                        DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_, &
                        TS,WATER,FR,n_iter_cool,fr_ice_thresh)

! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
    real,    intent(IN)    :: FR     (:,:)   ! fraction of surface (water/ice)
    integer, intent(IN)    :: WATER          ! subtile  number assigned to surface type: "WATER" 
    real,    intent(IN)    :: CM     (:,:)   ! transfer coefficient for wind
    real,    intent(IN)    :: UUA    (:)     ! zonal       wind
    real,    intent(IN)    :: VVA    (:)     ! meridional  wind
    real,    intent(IN)    :: UW     (:)     ! u-current
    real,    intent(IN)    :: VW     (:)     ! v-current
    real,    intent(IN)    :: SWN    (:)     ! net shortwave radiation incident at surface
    real,    intent(IN)    :: LHF    (:)     ! latent   heat flux
    real,    intent(IN)    :: SHF    (:)     ! sensible heat flux
    real,    intent(IN)    :: LWDNSRF(:)     ! downward longwave at surface
    real,    intent(IN)    :: ALW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: BLW    (:)     ! for linearized \sigma T^4
    integer, intent(IN)    :: n_iter_cool    ! number of iterations to compute cool-skin layer 
    real,    intent(IN)    :: fr_ice_thresh  ! threshold on ice fraction, sort of defines Marginal Ice Zone
    real,    intent(IN)    :: TS     (:,:)   ! skin temperature

    real,    intent(OUT)   :: USTARW_(:)     ! u_{*,w} 
    real,    intent(OUT)   :: DCOOL_ (:)     ! depth of cool-skin layer
    real,    intent(OUT)   :: TDROP_ (:)     ! temperature drop across cool-skin
    real,    intent(OUT)   :: SWCOOL_(:)     ! shortwave radiation absorbed in cool-skin 
    real,    intent(OUT)   :: QCOOL_ (:)     ! net heat flux in cool layer
    real,    intent(OUT)   :: BCOOL_ (:)     ! bouyancy in cool layer
    real,    intent(OUT)   :: LCOOL_ (:)     ! Saunder's parameter in cool layer

    real,    intent(INOUT) :: TXW    (:)     ! zonal      stress
    real,    intent(INOUT) :: TYW    (:)     ! meridional stress

!  !LOCAL VARIABLES

    integer         :: N, iter_cool
    real            :: ALPH, Qb, fC

    real, parameter :: NU_WATER        = 1.0E-6  ! kinematic viscosity of water  [m^2/s]
    real, parameter :: TherCond_WATER  = 0.563   ! Thermal conductivity of water [W/m/ K]
    real, parameter :: bigC            = &
          (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2


    do N = 1, NT  ! N is now looping over all tiles (NOT sub-tiles).

!      Stress over "open" water (or Marginal Ice Zone) depends on ocean currents
!      --------------------------------------------------------------------------
       TXW(N) = CM(N,WATER)*(UUA(N) - UW(N))
       TYW(N) = CM(N,WATER)*(VVA(N) - VW(N))

       if( FR(N,WATER) > fr_ice_thresh ) then 

!        Ustar in water has a floor of 2 \mu m/s
!        ----------------------------------------
         USTARW_(N) = max( 2.e-6, sqrt(sqrt(TXW(N)*TXW(N)+TYW(N)*TYW(N))/MAPL_RHO_SEAWATER) )

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! Cool skin layer- heat loss and temperature drop  @ top of interface layer !
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          DCOOL_(N)  = 1.e-3           ! initial guess for cool-skin layer thickness
          TDROP_(N)  = 0.2             ! guess for cool-skin tdrop. FINAL TDROP IS SENSITIVE TO INITIAL CHOICE. 3 ITER ENOUGH?

          cool_iter: do iter_cool = 1, n_iter_cool

!         Short wave absorbed in the cool layer. This is a modified version of Zeng and Beljaars, 2005
!         ----------------------------------------------------------------------------------------------

             fC  = 0.0685 + 11.0*DCOOL_(N) - (3.3e-5/DCOOL_(N))*(1.0-exp(-DCOOL_(N)/8.0E-4))
             fC  = max( fC, 0.01)        ! absorb at least 1% of shortwave in cool layer
             SWCOOL_(N) = SWN(N)*fC

!            Heat loss at top of skin (cool) layer
!            --------------------------------------

             QCOOL_(N)  = &
                          LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*( TS(N,WATER)-TDROP_(N)))) - &
                          SWCOOL_(N)

!            Bouyancy production in cool layer depends on surface cooling
!            and evap-salinity effect from surface. It does not depend on solar
!            heating, which is assumed to be uniform in cool layer. This last assumption
!            could be improved by including some NIR. For this calculation, we include
!            temperature dependence of the thermal expansion coefficient.
!            -------------------------------------------------------------------------------

             ALPH   = (0.6 + 0.0935*(TS(N,WATER)-MAPL_TICE))*1.E-4
             Qb     = QCOOL_(N) + ( (0.026*MAPL_CAPWTR)/(ALPH*MAPL_ALHL) )*LHF(N)
             BCOOL_(N) = (ALPH*MAPL_GRAV*Qb) / (MAPL_RHO_SEAWATER*MAPL_CAPWTR)

!            Saunders parameter
!            BigC = (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2  
!            -------------------------------------------------------------------------------

             if ( BCOOL_(N) > 0.0) then  ! Eqn(14) of F96
                LCOOL_(N)  = 6.0/( 1.0 + ( BCOOL_(N)*bigC / USTARW_(N)**4 )**0.75 )**(1./3.)
                DCOOL_(N)  = LCOOL_(N)*NU_WATER/USTARW_(N)
             else 
                LCOOL_(N)  = 6.0
                DCOOL_(N)  = min( LCOOL_(N)*NU_WATER/USTARW_(N), 1.e-2)    ! Prevent very thick cool layer depth
             end if

             TDROP_(N)    = max( 0.0, DCOOL_(N)*QCOOL_(N)/TherCond_WATER ) ! Eqn(4) & (13) of F96

          end do cool_iter

       else            ! FR(N, WATER) <= fr_ice_thresh
          USTARW_(N)     = MAPL_UNDEF
          DCOOL_ (N)     = MAPL_UNDEF
          TDROP_ (N)     = 0.0
          SWCOOL_(N)     = MAPL_UNDEF
          QCOOL_ (N)     = MAPL_UNDEF
          BCOOL_ (N)     = MAPL_UNDEF
          LCOOL_ (N)     = MAPL_UNDEF
       end if
    end do

  end subroutine COOL_SKIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: ALBSEA - Computes albedos as a function of $cos(\zeta)$ over ocean surfaces

!  !DESCRIPTION:
!        Compute albedo for ocean points
!          based on ceres

!  CERES ocean albedo at zth=.5 is 0.052. Our formulation gives .077
!    thus the scaling. The diffuse albedo is given by computing
!    the zth weighted average of the albedo over the hemisphere and
!    then applying the same scaling to match CERES.

!LLT: CERESFAC = 1           reduces to old formulation 1-5-05
!     CERESFAC = 0.052/0.077 is the Original CERES Factor
!     CERESFAC = 0.068/0.077 is the EROS Tuned Value

! !INTERFACE:
  subroutine ALBSEA (ALBVR,ALBVF,ALBNR,ALBNF,ZTH)

! !ARGUMENTS:

    real,    intent(IN)  :: ZTH  (:)
    real,    intent(OUT) :: ALBVR(:) ! visible beam    albedo
    real,    intent(OUT) :: ALBVF(:) ! visible diffuse albedo
    real,    intent(OUT) :: ALBNR(:) ! nearIr  beam    albedo
    real,    intent(OUT) :: ALBNF(:) ! nearIr  diffuse albedo

    real, parameter :: CERESFAC   = 0.068/0.077
!   real, parameter :: CERESFAC   = 1.0

    real, parameter :: OCNALBVF   = .08*CERESFAC
    real, parameter :: OCNALBNF   = .08*CERESFAC

    real, parameter :: A0         = 0.40670980*CERESFAC
    real, parameter :: A1         =-1.23236340*CERESFAC
    real, parameter :: A2         = 1.42240510*CERESFAC
    real, parameter :: A3         =-0.55573341*CERESFAC

! Beam albedos
!-------------

    ALBVR = A0+(A1+(A2+A3*ZTH)*ZTH)*ZTH
    ALBNR = ALBVR

! Diffuse albedos
!----------------

    ALBVF = OCNALBVF
    ALBNF = OCNALBNF

  end subroutine ALBSEA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: SIMPLE_SW_ABS - 
!  Implements two simple ways for the absorption of shortwave radiation,
!  as an alternative (for "TESTING" purposes) to using KPAR & KUVR based on PEN,
!  into an interface layer of typical depth = 2m, that is also called "depth of AOIL"

! !INTERFACE:
  subroutine SIMPLE_SW_ABS(USE_KPAR, depth, ZTH, SWN, PEN)

! !ARGUMENTS:

    integer, intent(IN)    :: USE_KPAR  ! absorption profile option
    real,    intent(IN)    :: ZTH       ! cosine of solar zenith angle
    real,    intent(IN)    :: depth     ! depth up to which shortwave needs to be absorbed
    real,    intent(IN)    :: SWN       ! net shortwave at surface of ocean, or at top of air/sea interface
    real,    intent(OUT)   :: PEN       ! shortwave penetrated below the depth    

!  local variables
    real  :: fW

    fW  = 0.0
    PEN = 0.0                ! initialize to zero

    if (USE_KPAR == -1) then
       ! Soloviev, 1982 shortwave absorption profile
       ! --------------------------------------------
       fW = 0.28*exp(-71.5*depth) + 0.27*exp(-2.8*depth) + 0.45*exp(-0.07*depth)

    else if (USE_KPAR == -2) then
       ! Paulson & Simpson, 1981- Taken from Gentemann et al, 2009
       ! ----------------------------------------------------------
       fW    = 0.237*exp(-(depth*ZTH)/34.84)  +  0.36*exp(-(depth*ZTH)/2.266)   + &
               0.179*exp(-(depth*ZTH)/0.0315) + 0.087*exp(-(depth*ZTH)/0.0055)  + &
                0.08*exp(-(depth*ZTH)/8.32e-4)+ 0.025*exp(-(depth*ZTH)/1.26e-4) + &
               0.025*exp(-(depth*ZTH)/3.13e-4)+ 0.007*exp(-(depth*ZTH)/7.82e-4) + &
              0.0004*exp(-(depth*ZTH)/1.44e-5)
    else
       if(MAPL_AM_I_ROOT()) print *, 'ERROR! Unknown use_kpar option: ', USE_KPAR
    end if

    PEN   = SWN * fW

  end subroutine SIMPLE_SW_ABS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module atmOcnIntlayer
