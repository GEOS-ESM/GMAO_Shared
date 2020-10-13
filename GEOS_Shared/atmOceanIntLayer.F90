module atmOcnIntlayer

! !USES:

use MAPL
use GEOS_UtilsMod, only: GEOS_QSAT, GEOS_DQSAT

implicit none
private

public  ALBSEA
public  AOIL_sfcLayer_T
public  water_RHO
public  AOIL_Shortwave_abs
public  AOIL_v0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_v0_HW - update mass in AOIL

!  !DESCRIPTION:
!        Computes an update to mass in AOIL (version 0)

! !INTERFACE:
  subroutine AOIL_v0_HW ( NT, DT, DO_DATASEA,           &
                          MaxWaterDepth, MinWaterDepth, &
                          FRWATER, SNO, EVP, RAIN, HW)

! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
    real,    intent(IN)    :: DT             ! time-step
    integer, intent(IN)    :: DO_DATASEA     ! 1:uncoupled (AGCM);   0:coupled (AOGCM)
    real,    intent(IN)    :: MaxWaterDepth  ! maximum depth of AOIL
    real,    intent(IN)    :: MinWaterDepth  ! minimum depth of AOIL

    real,    intent(IN)    :: FRWATER(:)     ! fr of water
    real,    intent(IN)    :: SNO(:)         ! snow fall   rate
    real,    intent(IN)    :: EVP(:)         ! evaporation rate
    real,    intent(IN)    :: RAIN(:)        ! rain        rate= liquid_water_convective_precipitation + liquid_water_large_scale_precipitation

    real,    intent(INOUT) :: HW(:)          ! mass of AOIL

!  !LOCAL VARIABLES
    !real         :: FRESH                    ! this should include freshwater flux from ??
    real         :: FRESHATM(NT)

    !FRESH = 0.

    ! Layer thickness; liquid precip goes right thru ice.
    ! FRESHATM is useful for mass flux balance.
    ! freshwater flux from atmosphere needs to be added to HW here since it carries zero enthalpy 
    !---------------------------------------------------------------------------------------------
    FRESHATM  = FRWATER*(SNO - EVP) + RAIN 

    !HW   = HW + DT*(FRESHATM + FRESH)
    HW   = HW + DT*(FRESHATM)

    if (DO_DATASEA == 0) then               ! coupled   mode
      HW = max( min(HW, (MaxWaterDepth*water_RHO('salt_water'))),  (MinWaterDepth*water_RHO('salt_water')))
    else                                    ! uncoupled mode
      HW = max( min(HW, (MaxWaterDepth*water_RHO('fresh_water'))), (MinWaterDepth*water_RHO('fresh_water')))
    endif

  end subroutine AOIL_v0_HW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_v0_S - update salinity in AOIL

!  !DESCRIPTION:
!        Computes an update to salinity (version 0)

! !INTERFACE:
  subroutine AOIL_v0_S (DT, HW, S_old, S_new)

! !ARGUMENTS:

    real,    intent(IN)    :: DT             ! time-step
    real,    intent(IN)    :: HW(:)          ! mass of AOIL
    real,    intent(IN)    :: S_old(:)       ! salinity * mass of AOIL
    real,    intent(OUT)   :: S_new(:)       ! salinity * mass of AOIL

!  !LOCAL VARIABLES
    real         :: FSALT                    ! this should come from CICE as import: salt flux due to sea ice melt

    FSALT = 0.
    S_new = (S_old+DT*1.e3*FSALT)/HW         ! multiply by 1000 to account for g->kg conversion

  end subroutine AOIL_v0_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_Shortwave_abs - shortwave radiation absorption in AOIL

!  !DESCRIPTION:
!        Computes the shortwave radiation that penetrates
!        through the AOIL and also the ocean top layer

! !INTERFACE:
  subroutine AOIL_Shortwave_abs (NT, DO_SKIN_LAYER, DO_DATASEA, &
                                 AOIL_depth, OGCM_top_thickness, HW, KUVR, KPAR, &
                                 ALBVRO, ALBVFO, DRUVR, DFUVR, DRPAR, DFPAR, &
                                 PEN, PEN_ocean)
! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
    integer, intent(IN)    :: DO_SKIN_LAYER  ! 0:No interface layer, 1:active, and accounts for change in SST
    integer, intent(IN)    :: DO_DATASEA     ! 1:uncoupled (AGCM);   0:coupled (AOGCM)
    real,    intent(IN)    :: AOIL_depth     ! depth of the AOIL
    real,    intent(IN)    :: OGCM_top_thickness ! thickness of OGCM top layer

    real,    intent(IN)    :: KUVR           ! UV radiation extinction_coefficient
    real,    intent(IN)    :: KPAR(:)        ! PAR extinction_coefficient
    real,    intent(IN)    :: HW(:)          ! mass of AOIL

    real,    intent(IN)    :: ALBVRO(:)      ! visible beam albedo
    real,    intent(IN)    :: ALBVFO(:)      ! visible diffuse albedo
    real,    intent(IN)    :: DRUVR(:)       ! surface_downwelling_uvr_beam_flux
    real,    intent(IN)    :: DFUVR(:)       ! surface_downwelling_uvr_diffuse_flux
    real,    intent(IN)    :: DRPAR(:)       ! surface_downwelling_par_beam_flux
    real,    intent(IN)    :: DFPAR(:)       ! surface_downwelling_par_diffuse_flux

    real,    intent(OUT)   :: PEN(:)         ! shortwave flux penetrated through AOIL_depth
    real,    intent(OUT)   :: PEN_ocean(:)   ! shortwave flux penetrated through OGCM_top_thickness

!  !LOCAL VARIABLES
    real         :: PUR(NT), PUF(NT), PPR(NT), PPF(NT)

!   init local variables
    PUR   = 0.0; PUF   = 0.0; PPR   = 0.0; PPF   = 0.0

    if (DO_SKIN_LAYER==0) then
      PEN       = 0.0
      PEN_ocean = 0.0
      RETURN
    endif

!   UV penetration
    if (DO_DATASEA /= 0) then                ! UNcoupled mode
        PEN = exp(-(KUVR/water_RHO('fresh_water'))*HW)
    else                                     ! coupled mode
        PEN = exp(-KUVR*AOIL_depth)
    endif
    PUR = (1.-ALBVRO)*DRUVR*PEN
    PUF = (1.-ALBVFO)*DFUVR*PEN

!   near-IR ("blue light" 490nm?)
    if (DO_DATASEA /= 0) then                ! UNcoupled mode
        PEN = exp(-(KPAR/water_RHO('fresh_water'))*HW)
    else                                     ! coupled mode
        PEN = exp(-KPAR*AOIL_depth)                 
    endif
    PPR = (1.-ALBVRO)*DRPAR*PEN
    PPF = (1.-ALBVFO)*DFPAR*PEN

    PEN = PUR + PUF + PPR + PPF              ! penetrated flux through AOIL_depth

    if (DO_DATASEA /= 0) then                ! UNcoupled mode
      PEN_ocean = 0.0
    else                                     ! coupled mode
      PEN_ocean =  exp(-KUVR*OGCM_top_thickness)
      PUR       =  (1.-ALBVRO)*DRUVR*PEN_ocean
      PUF       =  (1.-ALBVFO)*DFUVR*PEN_ocean
      PEN_ocean =  exp(-KPAR*OGCM_top_thickness)
      PPR       =  (1.-ALBVRO)*DRPAR*PEN_ocean
      PPF       =  (1.-ALBVFO)*DFPAR*PEN_ocean
      PEN_ocean =  PUR + PUF + PPR + PPF    ! penetrated flux through OGCM_top_thickness
    endif

  end subroutine AOIL_Shortwave_abs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: water_RHO - returns density of water: fresh/sea water from MAPL_Constants

!  !DESCRIPTION:

! !INTERFACE:
  function water_RHO (WATER_TYPE)

! !ARGUMENTS:

    character(len=*), intent(IN) :: WATER_TYPE     ! 'fresh_water' or 'sea_water'
    real                         :: water_RHO

    select case(trim(WATER_TYPE))
    case ('fresh_water')
      water_RHO = MAPL_RHOWTR
    case ('salt_water')
      water_RHO = MAPL_RHO_SEAWATER
    case default
      water_RHO = MAPL_UNDEF
      print *, ' Unknown option in water_RHO.'
    end select
    
  end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_sfcLayer_T - Connection to the surface layer: which "temperature" is to be passed from AOIL to sfclayer

!  !DESCRIPTION:
!        It is not exact "clear" which temperature is to be input to the sfclayer, that computes the bulk flux formulae 
!        Hence it may be considered that the temperature that is input to sfclayer is a "tunable parameter"

!        Based on:
!        (i) Akella and Suarez, 2018 "The Atmosphere-Ocean Interface Layer of the NASA
!         Goddard Earth Observing System Model and Data Assimilation System." GMAO Tech Memo, Vol 51.
!
!        Options for the temperatures are:
!        Num.     Var name       Definition                                      Source
!        --------------------------------------------------------------------------------
!        1.       TF             Temperature at the bottom of AOIL               from ocean
!        2.       TS             Temperature at the top of AOIL                  calculated from temperature profile in AOIL
!        3.       TW             Depth averaged mean temperature in the AOIL     from internal state
!        4.       TW             Depth averaged mean temperature in the AOIL     calculated from temperature profile in AOIL
!        5.       T(z)           Temperature at specified depth (z) within AOIL  calculated from temperature profile in AOIL and input depth

! !INTERFACE:
  subroutine AOIL_sfcLayer_T (WHICH_OPTION, depth, DO_DATASEA, MUSKIN, epsilon_d,  &
                              AOIL_depth, TW, TS_FOUND, TWMTF, DELTC,              &
                              T_OUT)

! !ARGUMENTS:

    character(len=*), intent(IN) :: WHICH_OPTION   ! See above description of options for the temperature
    integer,          intent(IN) :: DO_DATASEA     ! =1:uncoupled (AGCM); =0:coupled (AOGCM)
    real,             intent(IN) :: MUSKIN         ! exponent in T(z) profile in warm layer, based on Zeng & Beljaars, 2005, typically <= 1.
    real,             intent(IN) :: epsilon_d      ! (thickness of AOIL)/(thickness of OGCM top level), typically < 1.
    real,             intent(IN) :: depth          ! specified depth (z) used in above option 5
    real,             intent(IN) :: AOIL_depth     ! depth of AOIL

    real,             intent(IN) :: TW (:)         ! depth averaged mean temperature in the AOIL
    real,             intent(IN) :: TS_FOUND(:)    ! temperature at the bottom of AOIL
    real,             intent(IN) :: TWMTF(:)       ! difference: TW - TS_FOUND
    real,             intent(IN) :: DELTC(:)       ! temperature drop due to cool skin layer
    real,             intent(OUT):: T_OUT (:)      ! temperature that is to be input to the sfclayer

    select case(trim(WHICH_OPTION))
    case ('TF')
      T_OUT = TS_FOUND
    case ('TS')
      if (DO_DATASEA == 1) then
        T_OUT = TS_FOUND + ((1.+MUSKIN)/MUSKIN) * TWMTF            ! Eqn.(14) of Akella and Suarez, 2018
      else
        T_OUT = TS_FOUND + (1./MUSKIN + (1.-epsilon_d)) * TWMTF    ! RHS is from Eqn.(15) of Akella and Suarez, 2018
      endif
      T_OUT = T_OUT - DELTC                                        ! Eqn.(16) of Akella and Suarez, 2018
    case ('TW_from_internal')
      T_OUT = TW
    case ('TW_from_Tprof')
      T_OUT = TWMTF + TS_FOUND                                     ! Eqn.(5) of Akella and Suarez, 2018
    case ('T_at_depth')
!     ! T_{\delta} : top of warm layer or base of cool layer
      if (DO_DATASEA == 1) then
        T_OUT = TS_FOUND + ((1.+MUSKIN)/MUSKIN) * TWMTF            ! Eqn.(14) of Akella and Suarez, 2018
      else
        T_OUT = TS_FOUND + (1./MUSKIN + (1.-epsilon_d)) * TWMTF    ! RHS is from Eqn.(15) of Akella and Suarez, 2018
      endif

      where ( depth <= DELTC)    ! within cool layer
        T_OUT = T_OUT - ( 1.-depth/(DELTC + 1.e-6)) * DELTC        ! Eqn.(16) of Akella and Suarez, 2018; avoiding divide by zero.
      elsewhere                  !     in warm layer
        T_OUT = T_OUT - ( (depth/AOIL_depth)**MUSKIN) * ((1.+MUSKIN)/MUSKIN) * TWMTF ! neglect \delta/d since it is \ll 1.
      endwhere
    case default
      T_OUT = MAPL_UNDEF
      print *, ' Unknown option in AOIL_sfcLayer_T.'
    end select

  end subroutine AOIL_sfcLayer_T
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: SKIN_SST - Computes changes to SST in interface layer due to Cool Skin & Diurnal Warming 

!  !DESCRIPTION:
!        Described in section 2 of
!        https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.2988
!        Zeng and Beljaars, 2005

! !INTERFACE:
  subroutine SKIN_SST (DO_SKIN_LAYER, DO_DATASEA, NT,CM,UUA,VVA,UW,VW,HW,SWN,LHF,SHF,LWDNSRF,      &
                       ALW,BLW,PEN, PEN_ocean, STOKES_SPEED,DT,MUSKIN,TS_FOUNDi,DWARM_,TBAR_,      &
                       TXW,TYW,USTARW_,  &
                       DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_,TDEL_,SWWARM_,QWARM_,ZETA_W_,    &
                       PHIW_,LANGM_,TAUTW_,uStokes_,TS,TWMTS,TWMTF,DELTC,TW,FR,n_iter_cool,fr_ice_thresh, &
                       epsilon_d, do_grad_decay)

! !ARGUMENTS:

    integer, intent(IN)    :: DO_SKIN_LAYER  ! 0: No interface layer,     1: active, and accounts for change in SST
    integer, intent(IN)    :: DO_DATASEA     ! =1:uncoupled (AGCM); =0:coupled (AOGCM)
    integer, intent(IN)    :: NT             ! number of tiles
    real,    intent(IN)    :: FR     (:)     ! fraction of surface water
    real,    intent(IN)    :: CM     (:)     ! transfer coefficient for wind
    real,    intent(IN)    :: UUA    (:)     ! zonal       wind
    real,    intent(IN)    :: VVA    (:)     ! meridional  wind
    real,    intent(IN)    :: UW     (:)     ! u-current
    real,    intent(IN)    :: VW     (:)     ! v-current
    real,    intent(IN)    :: HW     (:)     ! mass of AOIL
    real,    intent(IN)    :: SWN    (:)     ! net shortwave radiation incident at surface
    real,    intent(IN)    :: LHF    (:)     ! latent   heat flux
    real,    intent(IN)    :: SHF    (:)     ! sensible heat flux
    real,    intent(IN)    :: LWDNSRF(:)     ! longwave at surface
    real,    intent(IN)    :: ALW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: BLW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: PEN    (:)     ! shortwave radiation that penetrates below interface layer
    real,    intent(IN)    :: PEN_ocean(:)   ! shortwave radiation that penetrates below top ocean model layer
    real,    intent(IN)    :: STOKES_SPEED   ! scalar value set for Stokes speed- place holder for output from Wave model
    real,    intent(IN)    :: DT             ! time-step
    real,    intent(IN)    :: MUSKIN         ! exponent of temperature: T(z) profile in warm layer
    real,    intent(IN)    :: TS_FOUNDi(:)   ! bulk SST (temperature at base of warm layer)
    integer, intent(IN)    :: n_iter_cool    ! number of iterations to compute cool-skin layer 
    real,    intent(IN)    :: fr_ice_thresh  ! threshold on ice fraction, sort of defines Marginal Ice Zone
    real,    intent(IN)    :: epsilon_d      ! (thickness of AOIL)/(thickness of OGCM top level), typically < 1.
    character(len=*), intent(IN) :: do_grad_decay   ! simulate a gradual decay of diurnal warming? yes or no. Follows Zeng and Beljaars, 2005.

    real,    intent(OUT)   :: DWARM_ (:)     ! depth of AOIL
    real,    intent(OUT)   :: TBAR_  (:)     ! copy of TW (also internal state) to export out
    real,    intent(OUT)   :: USTARW_(:)     ! u_{*,w} 
    real,    intent(OUT)   :: DCOOL_ (:)     ! depth of cool-skin layer
    real,    intent(OUT)   :: TDROP_ (:)     ! temperature drop across cool-skin
    real,    intent(OUT)   :: SWCOOL_(:)     ! shortwave radiation absorbed in cool-skin 
    real,    intent(OUT)   :: QCOOL_ (:)     ! net heat flux in cool layer
    real,    intent(OUT)   :: BCOOL_ (:)     ! bouyancy in cool layer
    real,    intent(OUT)   :: LCOOL_ (:)     ! Saunder's parameter in cool layer

    real,    intent(OUT)   :: TDEL_  (:)     ! temperature at top of warm layer
    real,    intent(OUT)   :: SWWARM_(:)     ! shortwave radiation absorbed in warm layer
    real,    intent(OUT)   :: QWARM_ (:)     ! net heat flux in warm layer
    real,    intent(OUT)   :: ZETA_W_(:)     ! stability parameter = dwarm/(Obukhov length)
    real,    intent(OUT)   :: PHIW_  (:)     ! similarity function
    real,    intent(OUT)   :: LANGM_ (:)     ! Langmuir number
    real,    intent(OUT)   :: TAUTW_ (:)     ! time-scale of relaxation to bulk SST (i.e., TS_FOUND)
    real,    intent(OUT)   :: uStokes_(:)    ! Stokes speed
    real,    intent(OUT)   :: TXW    (:)     ! zonal      stress
    real,    intent(OUT)   :: TYW    (:)     ! meridional stress

    real,    intent(OUT)   :: DELTC  (:)     ! "internal state" variable that has: TDROP_
    real,    intent(INOUT) :: TWMTS  (:)     ! "internal state" variable that has: TW - TS
    real,    intent(INOUT) :: TWMTF  (:)     ! "internal state" variable that has: TW - bulk SST
    real,    intent(INOUT) :: TW     (:)     ! "internal state" variable that has: TW
    real,    intent(INOUT) :: TS     (:)     ! skin temperature

!  !LOCAL VARIABLES

    integer         :: N, iter_cool
    real            :: ALPH, Qb, fC, fLA, X1, X2, dTw

    real, parameter :: RHO_SEAWATER    = 1022.0  ! sea water density             [kg/m^3]    ! Replace Usage of RHO_SEAWATER with MAPL_RHO_SEAWATER
    real, parameter :: NU_WATER        = 1.0E-6  ! kinematic viscosity of water  [m^2/s]
    real, parameter :: TherCond_WATER  = 0.563   ! Thermal conductivity of water [W/m/ K]
    real, parameter :: bigC            = &
          (16.0 * (MAPL_CAPWTR*MAPL_RHOWTR)**2 * NU_WATER**3) / TherCond_WATER**2

    do N = 1, NT  ! N is now looping over all tiles (NOT sub-tiles).

! Stress over "open" water (or Marginal Ice Zone) depends on ocean currents
!--------------------------------------------------------------------------

       TXW(N) = CM(N)*(UUA(N) - UW(N))
       TYW(N) = CM(N)*(VVA(N) - VW(N))

       if( FR(N) > fr_ice_thresh) then 

! Depth and mean temperature of interface layer
!----------------------------------------------

          DWARM_(N) = HW(N)/MAPL_RHOWTR                                                   ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER
          TBAR_(N)  = TS(N) + TWMTS(N)

! Ustar in water has a floor of 2 \mu m/s
!----------------------------------------

          USTARW_(N) = max( 2.e-6, sqrt(sqrt(TXW(N)*TXW(N)+TYW(N)*TYW(N))/MAPL_RHOWTR) )  ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cool skin layer- heat loss and temperature drop  @ top of interface layer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          DCOOL_(N)  = 1.e-3           ! initial guess for cool-skin layer thickness
          TDROP_(N)  = 0.2             ! guess for cool-skin tdrop. FINAL TDROP IS SENSITIVE TO INITIAL CHOICE. 3 ITER ENOUGH?

          COOL_SKIN: do iter_cool = 1, n_iter_cool

! Short wave absorbed in the cool layer. This is a modified version of Zhang and Beljaars, 2005
!----------------------------------------------------------------------------------------------

             fC  = 0.0685 + 11.0*DCOOL_(N) - (3.3e-5/DCOOL_(N))*(1.0-exp(-DCOOL_(N)/8.0E-4))
             fC  = max( fC, 0.01)        ! absorb at least 1% of shortwave in cool layer
             SWCOOL_(N) = SWN(N)*fC

! Heat loss at top of skin (cool) layer
!--------------------------------------

             X1 = LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*( TS(N)-TDROP_(N))))
             QCOOL_(N)  = X1 - SWCOOL_(N)

! Bouyancy production in cool layer depends on surface cooling
! and evap-salinity effect from surface. It does not depend on solar
! heating, which is assumed to be uniform in cool layer. This last assumption
! could be improved by including some NIR. For this calculation, we include
! temperature dependence of the thermal expansion coefficient.
!-------------------------------------------------------------------------------

             ALPH   = (0.6 + 0.0935*(TBAR_(N)-MAPL_TICE))*1.E-4
             Qb     = QCOOL_(N) + ( (0.026*MAPL_CAPWTR)/(ALPH*MAPL_ALHL) )*LHF(N)
             BCOOL_(N) = (ALPH*MAPL_GRAV*Qb) / (RHO_SEAWATER*MAPL_CAPWTR)                 ! replace RHO_SEAWATER with MAPL_RHO_SEAWATER

! Saunders parameter
! BigC = (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2  
!-------------------------------------------------------------------------------

             if ( BCOOL_(N) > 0.0) then  ! Eqn(14) of F96
                LCOOL_(N)  = 6.0/( 1.0 + ( BCOOL_(N)*bigC / USTARW_(N)**4 )**0.75 )**(1./3.)
                DCOOL_(N)  = LCOOL_(N)*NU_WATER/USTARW_(N)
             else 
                LCOOL_(N)  = 6.0
                DCOOL_(N)  = min( LCOOL_(N)*NU_WATER/USTARW_(N), 1.e-2)  ! Prevent very thick cool layer depth
             end if

             TDROP_(N)    = max( 0.0, DCOOL_(N)*QCOOL_(N)/TherCond_WATER ) ! Eqn(4) & (13) of F96

          end do COOL_SKIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Done with Cool skin layer.  Now turbluent heat flux at base of interface layer  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WARM_LAYER: if(DO_SKIN_LAYER==0) then   ! Warm layer temperature increase calculated based on definition of mean interface temperature.

             TDEL_(N)    = TS_FOUNDi(N)
             TS(N)       = TS_FOUNDi(N)
             TW(N)       = TS_FOUNDi(N)
             TBAR_(N)    = TS_FOUNDi(N)
             TWMTS(N)    = 0.
             TWMTF(N)    = 0.
             DELTC(N)    = 0.

             SWWARM_(N)  = MAPL_UNDEF
             QWARM_ (N)  = MAPL_UNDEF
             ZETA_W_(N)  = MAPL_UNDEF
             PHIW_(N)    = MAPL_UNDEF
             LANGM_(N)   = MAPL_UNDEF
             TAUTW_(N)   = MAPL_UNDEF

          else  ! use Takaya et al 2012

! Compute warm layer temperature increase based on Takaya et al, 2012
!--------------------------------------------------------------------

             ALPH        = (0.6 + 0.0935*(TBAR_(N)-MAPL_TICE))*1.E-4

! Short wave absorbed in the warm layer.
!--------------------------------------

             X1         = LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*TS(N)))

             if (DO_DATASEA == 0) then ! coupled   mode
              SWWARM_(N) = SWN(N) - PEN(N) - (epsilon_d/(1.-epsilon_d))* (PEN(N)-PEN_ocean(N))
             else                      ! uncoupled mode
              SWWARM_(N) = SWN(N) - PEN(N)
             endif
             QWARM_(N)  = SWWARM_(N) - X1

! Stability parameter & Similarity function
!------------------------------------------

             if ( trim(do_grad_decay) == 'yes') then
               dTw = ((1.+MUSKIN)/MUSKIN) * (TBAR_(N)-TS_FOUNDi(N))
               if ( (dTw > 0.) .and. (QWARM_(N)<0.)) then
                 ZETA_W_(N)  = (sqrt(DWARM_(N)*dTw*MUSKIN*MAPL_GRAV*ALPH)*MAPL_KARMAN)/ & ! F_d formulation
                               (sqrt(5.) * USTARW_(N))
               else
                 ZETA_W_(N)  = (DWARM_(N)*MAPL_KARMAN*MAPL_GRAV*ALPH*QWARM_(N)) / & ! zeta_w = dwarm/obukhov length
                               (RHO_SEAWATER*MAPL_CAPWTR*USTARW_(N)**3)             ! replace RHO_SEAWATER with MAPL_RHO_SEAWATER
               endif
             else
               ZETA_W_(N)  = (DWARM_(N)*MAPL_KARMAN*MAPL_GRAV*ALPH*QWARM_(N)) / & ! zeta_w = dwarm/obukhov length
                             (RHO_SEAWATER*MAPL_CAPWTR*USTARW_(N)**3)             ! replace RHO_SEAWATER with MAPL_RHO_SEAWATER
             endif 

             if ( ZETA_W_(N) >= 0.0) then   ! Takaya: Eqn(5)
                PHIW_(N) = 1. + (5*ZETA_W_(N) + 4.*ZETA_W_(N)**2)/(1+3.*ZETA_W_(N)+0.25*ZETA_W_(N)**2)
             else
                PHIW_(N) = 1.0/sqrt(1.-16.*ZETA_W_(N))
             end if

! Langmuir number- need imports from Wave Model
!----------------------------------------------

             uStokes_(N) = STOKES_SPEED
             LANGM_(N)   = sqrt(USTARW_(N)/uStokes_(N))
             fLA         = LANGM_(N)**(-0.66667)           ! Takaya: Eqn(6)

             IF (fLA       <= 1.0) fLA = 1.0               ! Limit range of fLa to be >=1
             IF (ZETA_W_(N)<= 0.0) fLA = 1.0               ! Apply fLa to stable conditions only

             if (DO_DATASEA == 0) then ! coupled   mode
               TAUTW_(N)   = &
                 ((1.-epsilon_d)*DWARM_(N)*PHIW_(N))/(MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))
               X2          = DT * &
                 ( (MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))/(DWARM_(N)*PHIW_(N) * (1.-epsilon_d)) )
             else                      ! uncoupled mode
               TAUTW_(N)   = &
                 (DWARM_(N)*PHIW_(N))/(MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))
               X2          = DT * &
                 ( (MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))/(DWARM_(N)*PHIW_(N)) )
             endif

! We DO NOT include cool-skin tdrop in TW, therefore, we now save TW

             TW(N) = TS_FOUNDi(N) + ( 1.0/(1.+X2))  *    (TBAR_(N) - TS_FOUNDi(N))
             TS(N) = TS(N)  + ((1.0+MUSKIN)/MUSKIN) *    (TW(N)    - TBAR_(N))

             TDEL_(N)    = TS_FOUNDi(N) + ((1.0+MUSKIN)/MUSKIN) * MAX(TW(N)    - TS_FOUNDi(N), 0.0)
             TBAR_(N)    = TW(N)

             TS(N)    = TDEL_(N) - TDROP_(N)
             TWMTS(N) = TW(N)    - TS(N)
             TWMTF(N) = 0.0
!            TWMTF(N) = TW(N)    - TS_FOUNDi(N)  ! This will cause non-zero diff in internal/checkpoint, but ZERO DIFF in OUTPUT.
          end if WARM_LAYER

       else            ! FR(N) <= fr_ice_thresh
          DCOOL_ (N)     = MAPL_UNDEF
          LCOOL_ (N)     = MAPL_UNDEF
          DWARM_ (N)     = MAPL_UNDEF
          TBAR_  (N)     = MAPL_UNDEF
          TDROP_ (N)     = MAPL_UNDEF
          QCOOL_ (N)     = MAPL_UNDEF
          USTARW_(N)     = MAPL_UNDEF
          SWCOOL_(N)     = MAPL_UNDEF
          BCOOL_ (N)     = MAPL_UNDEF
          TDEL_  (N)     = MAPL_UNDEF
          TWMTS  (N)     = 0.0
          TWMTF  (N)     = 0.0
          QWARM_ (N)     = MAPL_UNDEF
          SWWARM_(N)     = MAPL_UNDEF
          PHIW_  (N)     = MAPL_UNDEF
          LANGM_ (N)     = MAPL_UNDEF
          TAUTW_ (N)     = MAPL_UNDEF
          ZETA_W_(N)     = MAPL_UNDEF
       end if
    end do

    DELTC = 0.0
!   DELTC = TDROP_   ! This will cause non-zero diff in internal/checkpoint, but ZERO DIFF in OUTPUT.

  end subroutine SKIN_SST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: SURFACE_FLUX_UPDATE - update using surface (atmospheric) fluxes only

!  !DESCRIPTION:
!        Computes an update using surface fluxes from atmospheric model

! !INTERFACE:
  subroutine surf_hflux_update (NT,DO_SKIN_LAYER,DO_DATASEA,MUSKIN,DT,epsilon_d,  &
             CFT,CFQ,SH,EVAP,DSH,DEV,THATM,QHATM,PS,FRWATER,HW,SNO,               &
             LWDNSRF,SWN,ALW,BLW,SHF,LHF,EVP,DTS,DQS,TS,QS,TWMTS,TWMTF)

! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
    integer, intent(IN)    :: DO_SKIN_LAYER  ! 0:No interface layer, 1:active, and accounts for change in SST
    integer, intent(IN)    :: DO_DATASEA     ! 1:uncoupled (AGCM);   0:coupled (AOGCM)
    real,    intent(IN)    :: MUSKIN         ! exponent of temperature: T(z) profile in warm layer
    real,    intent(IN)    :: DT             ! time-step
    real,    intent(IN)    :: epsilon_d      ! (thickness of AOIL)/(thickness of OGCM top level), typically < 1.

    real,    intent(IN)    :: CFT(:)         ! sensible    heat transfer coefficient
    real,    intent(IN)    :: CFQ(:)         ! evaporation heat transfer coefficient
    real,    intent(IN)    :: SH(:)          ! upward_sensible_heat_flux 
    real,    intent(IN)    :: EVAP(:)        ! evaporation
    real,    intent(IN)    :: DSH(:)         ! derivative_of_upward_sensible_heat_flux
    real,    intent(IN)    :: DEV(:)         ! derivative_of_evaporation
    real,    intent(IN)    :: THATM(:)       ! effective_surface_skin_temperature
    real,    intent(IN)    :: QHATM(:)       ! effective_surface_specific_humidity
    real,    intent(IN)    :: PS(:)          ! surface pressure
    real,    intent(IN)    :: FRWATER(:)     ! 1. - fr of sea ice
    real,    intent(IN)    :: HW(:)          ! mass of AOIL
    real,    intent(IN)    :: SNO(:)         ! snow fall 
    real,    intent(IN)    :: LWDNSRF(:)     ! surface_downwelling_longwave_flux
    real,    intent(IN)    :: SWN(:)         ! shortwave radiation absorbed in AOIL
    real,    intent(IN)    :: ALW(:)         ! linearization_of_surface_upwelling_longwave_flux
    real,    intent(IN)    :: BLW(:)         ! linearization_of_surface_upwelling_longwave_flux

    real,    intent(OUT)   :: SHF (:)        ! sensible heat flux
    real,    intent(OUT)   :: LHF (:)        ! latent   heat flux
    real,    intent(OUT)   :: EVP (:)        ! evaporation `heat' flux
    real,    intent(OUT)   :: DTS (:)        ! change in skin temperature
    real,    intent(OUT)   :: DQS (:)        ! change in specific humidity

    real,    intent(INOUT) :: TS  (:)        ! skin temperature
    real,    intent(INOUT) :: QS  (:)        ! specific humidity
    real,    intent(INOUT) :: TWMTS  (:)     ! "internal state" variable that has: TW - TS
    real,    intent(INOUT) :: TWMTF  (:)     ! "internal state" variable that has: TW - TF

!  !LOCAL VARIABLES
    real         :: SHD(NT), EVD(NT), QFLX(NT), DTX(NT), DTY

    EVP = CFQ* (EVAP + DEV * (QS-QHATM))   ! evaporation "flux"
    SHF = CFT* (SH   + DSH * (TS-THATM))   ! sensible heat flux

    SHD = CFT*DSH                                                 ! d (sensible heat flux)/d Ts
    EVD = CFQ*DEV*GEOS_DQSAT(TS, PS, RAMP=0.0, PASCALS=.TRUE.)    ! d (evap)/ d Ts

    QFLX = LWDNSRF - (ALW + BLW*TS) - SHF                         ! net longwave - sensible heat flux

    ! FR accounts for fraction of water/ice
    DTX = DT*FRWATER / (MAPL_CAPWTR*HW)

    if (DO_SKIN_LAYER == 0) then
      DTX = 0.
    else
      DTX = DTX*((MUSKIN+1.-MUSKIN*epsilon_d)/MUSKIN) ! note: epsilon_d is = 0. in uncoupled mode (DO_DATASEA == 1)
    endif

    ! DTY accounts for ice on top of water. Part of Shortwave is absorbed by ice and rest goes to warm water.
    ! skin layer only absorbs the portion of SW radiation passing thru the bottom of ice MINUS the portion passing thru the skin layer    
    ! Penetrated shortwave from sea ice bottom + associated ocean/ice heat flux
    DTY = 0. ! Revisit above with CICE6 [Nov, 2019] and compute DTY = DT / (SALTWATERCAP*HW) * (PENICE * FI + FHOCN)

    DTS = DTX * ( QFLX + SWN - EVP*MAPL_ALHL - MAPL_ALHF*SNO ) + DTY ! add net SW, minus latent heat flux to QFLX. Keep DTY (=0) for ZERO DIFF in OUTPUT.
    DTS = DTS / (1.0 + DTX*(BLW+SHD+EVD*MAPL_ALHL))     ! implicit solution for an update in temperature: DTS. Tnew = Told + DTS

    EVP = EVP + EVD * DTS                               ! update evaporation
    SHF = SHF + SHD * DTS                               ! update sensible heat flux
    LHF = EVP * MAPL_ALHL                               ! update latent   heat flux

    ! update temperature, moisture and mass
    TS   = TS + DTS
    DQS  = GEOS_QSAT(TS, PS, RAMP=0.0, PASCALS=.TRUE.) - QS
    QS   = QS + DQS

    ! updated TS implies an update in TWMTS(=TW-TS) and TWMTF(=TW-TF)
    if(DO_SKIN_LAYER == 0) then
      TWMTS = 0.
      TWMTF = 0.
    else
      TWMTS = TWMTS - (1.0/(MUSKIN+1.0))*DTS
      if (DO_DATASEA == 0) then            ! coupled   mode
         TWMTF = TWMTF + (MUSKIN/(1.+MUSKIN-MUSKIN*epsilon_d))*DTS
      else                                 ! uncoupled mode
         TWMTF = 0.
!        TWMTF = TWMTF + (MUSKIN/(1.+MUSKIN))*DTS ! This will cause non-zero diff in internal/checkpoint, but ZERO DIFF in OUTPUT.
      end if
    end if         

  end subroutine surf_hflux_update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: AOIL_V1 - version number: 1 of the Atmosphere Ocean Interface Layer (AOIL)

!  !DESCRIPTION:
!        Wraps the sequence of steps involved in updating the internal (state) variables

! !INTERFACE:
  subroutine AOIL_v0 (NT, DO_SKIN_LAYER, DO_DATASEA, n_iter_cool, fr_ice_thresh, do_grad_decay,                 &
                      DT, MUSKIN, epsilon_d, MaxWaterDepth, MinWaterDepth, MaxSalinity, MinSalinity,            &
                      STOKES_SPEED, CM, CFT, CFQ, SH, EVAP, DSH, DEV, THATM, QHATM, PS, SNO, RAIN, UUA, VVA,    &
                      UW, VW, FRWATER, SWN, SWN_surf, PEN, PEN_ocean, LWDNSRF, ALW, BLW,                        &
                      HH_in, TS_in, SS_in, QS_in, TS_FOUND,                                                     &
                      DWARM, TBAR, USTARW, DCOOL, TDROP, SWCOOL, QCOOL, BCOOL, LCOOL,                           &
                      TDEL, SWWARM, QWARM, ZETA_W, PHIW, LANGM, TAUTW, uStokes, TXW, TYW,                       &
                      SHF, LHF, EVP, DTS, DQS, DELTC, HW, TW, SW, TWMTS, TWMTF,                                 &
                      do_update_fluxes_AOIL_second_step)

! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
    integer, intent(IN)    :: DO_SKIN_LAYER  ! 0: No interface layer, 1: active, and accounts for change in SST
    integer, intent(IN)    :: DO_DATASEA     ! =1:uncoupled (AGCM); =0:coupled (AOGCM)
    integer, intent(IN)    :: n_iter_cool    ! number of iterations to compute cool-skin layer
    real,    intent(IN)    :: fr_ice_thresh  ! threshold on ice fraction, sort of defines Marginal Ice Zone

    character(len=*), intent(IN) :: do_grad_decay   ! simulate a gradual decay of diurnal warming? yes or no. Follows Zeng and Beljaars, 2005.
    character(len=*), intent(IN) :: do_update_fluxes_AOIL_second_step   ! update DTS, DQS, EVP, SHF, LHF after 2nd implicit update of temperature

    real,    intent(IN)    :: DT             ! time-step
    real,    intent(IN)    :: MUSKIN         ! exponent of temperature: T(z) profile in warm layer
    real,    intent(IN)    :: epsilon_d      ! (thickness of AOIL)/(thickness of OGCM top level), typically < 1.
    real,    intent(IN)    :: MaxWaterDepth  ! maximum depth of AOIL
    real,    intent(IN)    :: MinWaterDepth  ! minimum depth of AOIL
    real,    intent(IN)    :: MaxSalinity    ! maximum salinity in AOIL
    real,    intent(IN)    :: MinSalinity    ! minimum salinity in AOIL
    real,    intent(IN)    :: STOKES_SPEED   ! scalar value set for Stokes speed- place holder for output from Wave model

    real,    intent(IN)    :: CM(:)          ! transfer coefficient for wind
    real,    intent(IN)    :: CFT(:)         ! sensible    heat transfer coefficient
    real,    intent(IN)    :: CFQ(:)         ! evaporation heat transfer coefficient
    real,    intent(IN)    :: SH(:)          ! upward_sensible_heat_flux
    real,    intent(IN)    :: EVAP(:)        ! evaporation    
    real,    intent(IN)    :: DSH(:)         ! derivative_of_upward_sensible_heat_flux
    real,    intent(IN)    :: DEV(:)         ! derivative_of_evaporation 
    real,    intent(IN)    :: THATM(:)       ! effective_surface_skin_temperature
    real,    intent(IN)    :: QHATM(:)       ! effective_surface_specific_humidity
    real,    intent(IN)    :: PS(:)          ! surface pressure 
    real,    intent(IN)    :: SNO(:)         ! snow fall
    real,    intent(IN)    :: RAIN(:)        ! rain = PCU+PLS = liquid_water_convective_precipitation + liquid_water_large_scale_precipitation
    real,    intent(IN)    :: UUA    (:)     ! zonal       wind
    real,    intent(IN)    :: VVA    (:)     ! meridional  wind
    real,    intent(IN)    :: UW     (:)     ! u-current
    real,    intent(IN)    :: VW     (:)     ! v-current

    real,    intent(IN)    :: FRWATER(:)     ! 1. - fr of sea ice
    real,    intent(IN)    :: SWN(:)         ! shortwave radiation absorbed in AOIL
    real,    intent(IN)    :: SWN_surf(:)    ! net shortwave radiation incident at surface
    real,    intent(IN)    :: PEN    (:)     ! shortwave radiation that penetrates below interface layer
    real,    intent(IN)    :: PEN_ocean(:)   ! shortwave radiation that penetrates below top ocean model layer 
    real,    intent(IN)    :: LWDNSRF(:)     ! surface_downwelling_longwave_flux
    real,    intent(IN)    :: ALW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: BLW    (:)     ! for linearized \sigma T^4 

    real,    intent(IN)    :: TS_FOUND(:)    ! bulk SST (temperature at base of warm layer)

    real,    intent(OUT)   :: DWARM (:)      ! depth of AOIL
    real,    intent(OUT)   :: TBAR  (:)      ! copy of TW (also internal state) to export out
    real,    intent(OUT)   :: USTARW(:)      ! u_{*,w} 
    real,    intent(OUT)   :: DCOOL (:)      ! depth of cool-skin layer
    real,    intent(OUT)   :: TDROP (:)      ! temperature drop across cool-skin
    real,    intent(OUT)   :: SWCOOL(:)      ! shortwave radiation absorbed in cool-skin 
    real,    intent(OUT)   :: QCOOL (:)      ! net heat flux in cool layer
    real,    intent(OUT)   :: BCOOL (:)      ! bouyancy in cool layer
    real,    intent(OUT)   :: LCOOL (:)      ! Saunder's parameter in cool layer

    real,    intent(OUT)   :: TDEL  (:)      ! temperature at top of warm layer
    real,    intent(OUT)   :: SWWARM(:)      ! shortwave radiation absorbed in warm layer
    real,    intent(OUT)   :: QWARM (:)      ! net heat flux in warm layer
    real,    intent(OUT)   :: ZETA_W(:)      ! stability parameter = dwarm/(Obukhov length)
    real,    intent(OUT)   :: PHIW  (:)      ! similarity function
    real,    intent(OUT)   :: LANGM (:)      ! Langmuir number
    real,    intent(OUT)   :: TAUTW (:)      ! time-scale of relaxation to bulk SST (i.e., TS_FOUND)
    real,    intent(OUT)   :: uStokes(:)     ! Stokes speed
    real,    intent(OUT)   :: TXW    (:)     ! zonal      stress
    real,    intent(OUT)   :: TYW    (:)     ! meridional stress

    real,    intent(OUT)   :: SHF (:)        ! sensible heat flux
    real,    intent(OUT)   :: LHF (:)        ! latent   heat flux
    real,    intent(OUT)   :: EVP (:)        ! evaporation `heat' flux
    real,    intent(OUT)   :: DTS (:)        ! change in skin temperature
    real,    intent(OUT)   :: DQS (:)        ! change in specific humidity

    real,    intent(OUT)   :: DELTC  (:)     ! "internal state" variable that has: TDROP_

    real,    intent(INOUT) :: HH_in(:)       ! initial mass of AOIL, will be updated
    real,    intent(INOUT) :: TS_in(:)       ! initial skin temperature, will be updated
    real,    intent(INOUT) :: SS_in(:)       ! initial salinity * mass of AOIL, will be updated
    real,    intent(INOUT) :: QS_in(:)       ! initial specific humidity, will be updated

    real,    intent(INOUT) :: HW     (:)     ! "internal state" variable that has: HW  (analog of HH)
    real,    intent(INOUT) :: TW     (:)     ! "internal state" variable that has: TW  (analog of TS)
    real,    intent(INOUT) :: SW     (:)     ! "internal state" variable that has: SW  (analog of SS)
    real,    intent(INOUT) :: TWMTS  (:)     ! "internal state" variable that has: TW - TS
    real,    intent(INOUT) :: TWMTF  (:)     ! "internal state" variable that has: TW - TF

!  !LOCAL VARIABLES
    real         :: SHD(NT), EVD(NT), TS0(NT), QS0(NT)!, HH0(NT), SS0(NT)

    TS0 = TS_in
    QS0 = QS_in
!   HH0 = HH_in
!   SS0 = SS_in

!   Update TS, QS, DTS, DQS; fluxes (EVP, SHF, LHF), internals (TWMTS, TWMTF) due to surface fluxes from atmosphere
!   ---------------------------------------------------------------------------------------------------------------
    call surf_hflux_update (NT,DO_SKIN_LAYER,DO_DATASEA,MUSKIN,DT,epsilon_d,  &
             CFT,CFQ,SH,EVAP,DSH,DEV,THATM,QHATM,PS,FRWATER,HH_in,SNO,           &
             LWDNSRF,SWN,ALW,BLW,SHF,LHF,EVP,DTS,DQS,TS_in,QS_in,TWMTS,TWMTF)

!   Update mass of AOIL using fresh water flux from atmosphere
!   ----------------------------------------------------------
    call AOIL_v0_HW (NT,DT,DO_DATASEA,MaxWaterDepth,MinWaterDepth, &
                      FRWATER,SNO,EVP,RAIN,HH_in)

!   Copy back to internal variables
!   -------------------------------
    TW = TS_in + TWMTS
    HW = HH_in

!   Update salinity
!   ---------------
    call AOIL_v0_S (DT,HW,SS_in,SW)
    where (.not. (abs(UW) > 0.0 .or. abs(VW) > 0.0))
       SW = max(min(SW,MAXSALINITY),MINSALINITY)
    endwhere

!   Cool-skin and diurnal warm layer, latter uses ocean variables (TS_FOUND)
!   it updates temperatures: TS, TW, DELTC and differences: TWMTF, TWMTS but
!   not the associated turbulent heat fluxes
!   ------------------------------------------------------------------------
    call SKIN_SST (DO_SKIN_LAYER,DO_DATASEA,NT,CM,UUA,VVA,UW,VW,HW,SWN_surf,LHF,SHF,LWDNSRF,          &
                   ALW,BLW,PEN,PEN_OCEAN,STOKES_SPEED,DT,MUSKIN,TS_FOUND,DWARM,TBAR,TXW,TYW,USTARW,   &
                   DCOOL,TDROP,SWCOOL,QCOOL,BCOOL,LCOOL,TDEL,SWWARM,QWARM,ZETA_W,                     &
                   PHIW,LANGM,TAUTW,uStokes,TS_in,TWMTS,TWMTF,DELTC,TW,FRWATER,n_iter_cool,           &
                   fr_ice_thresh,epsilon_d,do_grad_decay)


!   Following updates variables/fluxes using the final temperature from cool-skin and warm-layer effects
!   In coupled mode, this is really important
!   ----------------------------------------------------------------------------------------------------
    if (trim(do_update_fluxes_AOIL_second_step) == 'yes') then 
      DTS = TS_in - TS0                                            ! final change in TS
      DQS  = GEOS_QSAT(TS_in, PS, RAMP=0.0, PASCALS=.TRUE.) - QS0  ! final change in QS
      QS_in= QS0 + DQS
!     DHH  = HW - HH0
!     DSS  = SW - SS0

!     Updating variables in a sequence of two (implicit) updates is not desirable, as one can see, 
!     DSH and DEV are _out of sync_. The next version of AOIL will implement a single implicit update
!     -----------------------------------------------------------------------------------------------
      SHD = CFT*DSH                                                 ! d (sensible heat flux)/d Ts
      EVD = CFQ*DEV*GEOS_DQSAT(TS_in, PS, RAMP=0.0, PASCALS=.TRUE.) ! d (evap)/ d Ts
      EVP = EVP + EVD * DTS                                         ! update evaporation
      SHF = SHF + SHD * DTS                                         ! update sensible heat flux
      LHF = EVP * MAPL_ALHL                                         ! update latent   heat flux
    endif

  end subroutine AOIL_v0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module atmOcnIntlayer
