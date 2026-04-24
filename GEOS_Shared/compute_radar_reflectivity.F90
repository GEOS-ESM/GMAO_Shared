module GEOS_RadarMod

   implicit none

   private
   public :: compute_radar_reflectivity

contains

FUNCTION compute_radar_reflectivity(PRS, TMK, QVP, QRAIN, QSNOW, &
                                    QGRAUPEL, QHAIL, disable_variable_intercept_params, &
                                    liqskin_snow, liqskin_graupel, liqskin_hail) RESULT(DBZ)

  IMPLICIT NONE

  ! State and Required Hydrometeor Arguments (1D column slices)
  REAL, INTENT(IN)    :: PRS(:)
  REAL, INTENT(IN)    :: TMK(:)
  REAL, INTENT(INOUT) :: QVP(:)
  REAL, INTENT(INOUT) :: QRAIN(:)
  REAL, INTENT(INOUT) :: QSNOW(:)
  
  ! Optional Hydrometeor Arguments (1D column slices)
  REAL, INTENT(INOUT), OPTIONAL :: QGRAUPEL(:)
  REAL, INTENT(INOUT), OPTIONAL :: QHAIL(:)

  ! Optional Configuration Arguments
  LOGICAL, INTENT(IN), OPTIONAL :: disable_variable_intercept_params
  LOGICAL, INTENT(IN), OPTIONAL :: liqskin_snow
  LOGICAL, INTENT(IN), OPTIONAL :: liqskin_graupel
  LOGICAL, INTENT(IN), OPTIONAL :: liqskin_hail

  ! Return Array (Sized automatically based on the vertical column size)
  REAL :: DBZ(SIZE(PRS, 1))

  ! Local Configuration Variables
  LOGICAL :: use_var_intercepts
  LOGICAL :: has_graupel, has_hail
  LOGICAL :: do_liq_snow, do_liq_graupel, do_liq_hail

  ! Local Variables
  INTEGER :: K, LDIM
  REAL :: tmk_val, prs_val, qv_val, qr_val, qs_val, qg_val, qh_val
  REAL :: TEMP_C, VIRTUAL_T, RHOAIR, Z_E
  REAL :: GONV, HONV, RONV, SONV
  REAL :: FACTOR_G, FACTOR_H, FACTOR_R, FACTOR_S
  REAL :: FACTOR_G_LIQ, FACTOR_H_LIQ, FACTOR_S_LIQ
  REAL :: FACTORB_G, FACTORB_H, FACTORB_S
  REAL :: ZE_R, ZE_S, ZE_G, ZE_H
  REAL :: PI_RHO_G, PI_RHO_H, INV_RON_DELQR0

  ! Constants used to calculate variable intercepts
  REAL :: R1, RON, RON2, SON
  REAL :: GON_CONST, HON_CONST
  REAL :: RON_MIN, RON_QR0, RON_DELQR0
  REAL :: RON_CONST1R, RON_CONST2R
  
  ! Constant intercepts
  REAL :: RN0_R, RN0_S, RN0_G, RN0_H
  
  ! Other constants
  REAL :: RHO_R, RHO_S, RHO_G, RHO_H
  REAL :: GAMMA_SEVEN, ALPHA, RHOWAT, CELKEL, PI, RD

  ! Extract vertical dimension dynamically
  LDIM = SIZE(PRS, 1)

  ! Evaluate OPTIONAL arguments exactly once
  has_graupel = PRESENT(QGRAUPEL)
  has_hail    = PRESENT(QHAIL)
  
  ! Variable intercepts default to ON
  use_var_intercepts = .TRUE.
  IF (PRESENT(disable_variable_intercept_params)) THEN
      IF (disable_variable_intercept_params) use_var_intercepts = .FALSE.
  END IF

  ! Liquid skin defaults to OFF
  do_liq_snow    = .FALSE.
  do_liq_graupel = .FALSE.
  do_liq_hail    = .FALSE.
  IF (PRESENT(liqskin_snow))    do_liq_snow    = liqskin_snow
  IF (PRESENT(liqskin_graupel)) do_liq_graupel = liqskin_graupel
  IF (PRESENT(liqskin_hail))    do_liq_hail    = liqskin_hail

  ! Constants
  R1 = 1.D-15
  RON = 8.D6
  RON2 = 1.D10
  SON = 2.D7
  RON_MIN = 8.D6
  RON_QR0 = 0.00010D0
  RON_DELQR0 = 0.25D0 * RON_QR0
  RON_CONST1R = (RON2 - RON_MIN) * 0.5D0
  RON_CONST2R = (RON2 + RON_MIN) * 0.5D0
  
  ! Precompute Division Constant
  INV_RON_DELQR0 = 1.D0 / RON_DELQR0

  GAMMA_SEVEN = 720.D0
  RHOWAT = 1000.D0
  RHO_R = RHOWAT
  RHO_S = 100.D0
  ALPHA = 0.224D0
  CELKEL = 273.15D0
  PI = 3.141592653589793D0
  RD = 287.04D0

  ! Graupel Constants
  RHO_G     = 400.D0
  GON_CONST = 5.D7
  RN0_G     = 4.D6
  PI_RHO_G  = PI * RHO_G

  ! Hail Constants
  RHO_H     = 917.D0
  HON_CONST = 4.D4
  RN0_H     = 4.D4
  PI_RHO_H  = PI * RHO_H

  ! Precalculate Base Factors (hoisted outside loops)
  FACTOR_R = GAMMA_SEVEN * 1.D18 * (1.D0 / (PI * RHO_R))**1.75D0
  FACTOR_S = GAMMA_SEVEN * 1.D18 * (1.D0 / (PI * RHO_S))**1.75D0 * &
             (RHO_S / RHOWAT)**2 * ALPHA
  FACTOR_G = GAMMA_SEVEN * 1.D18 * (1.D0 / (PI_RHO_G))**1.75D0 * &
             (RHO_G / RHOWAT)**2 * ALPHA
  FACTOR_H = GAMMA_SEVEN * 1.D18 * (1.D0 / (PI_RHO_H))**1.75D0 * &
             (RHO_H / RHOWAT)**2 * ALPHA

  ! Precalculate Liquid Skin Division Factors
  FACTOR_S_LIQ = FACTOR_S / ALPHA
  FACTOR_G_LIQ = FACTOR_G / ALPHA
  FACTOR_H_LIQ = FACTOR_H / ALPHA

  ! Main Compute Loop - Pure Sequential over 1D Column
  DO K = 1, LDIM

      ! 1. Apply negative bounds checks AND extract to local scalars
      IF (QVP(K) < 0.0)   QVP(K) = 0.0
      IF (QRAIN(K) < 0.0) QRAIN(K) = 0.0
      IF (QSNOW(K) < 0.0) QSNOW(K) = 0.0
      
      qv_val  = QVP(K)
      qr_val  = QRAIN(K)
      qs_val  = QSNOW(K)
      tmk_val = TMK(K)
      prs_val = PRS(K)

      ! Initialize optional species
      ZE_G = 0.D0
      ZE_H = 0.D0

      ! Thermodynamics
      VIRTUAL_T = tmk_val * (0.622D0 + qv_val) / (0.622D0 * (1.D0 + qv_val))
      RHOAIR = prs_val / (RD * VIRTUAL_T)

      ! 2. Brightband Factor Assignment
      FACTORB_S = FACTOR_S
      FACTORB_G = FACTOR_G
      FACTORB_H = FACTOR_H

      IF (tmk_val > CELKEL) THEN
          IF (do_liq_snow)    FACTORB_S = FACTOR_S_LIQ
          IF (do_liq_graupel) FACTORB_G = FACTOR_G_LIQ
          IF (do_liq_hail)    FACTORB_H = FACTOR_H_LIQ
      END IF

      ! 3. Variable Intercept Parameters
      IF (use_var_intercepts) THEN
          TEMP_C = DMIN1(-0.001D0, tmk_val - CELKEL)
          SONV = DMIN1(2.0D8, 2.0D6 * EXP(-0.12D0 * TEMP_C))

          ! Simulate large snow aggregates by halving the intercept parameter
          IF (qs_val > R1) THEN
              SONV = SONV * 0.5D0
          END IF

          RONV = RON2
          IF (qr_val > R1) THEN
              RONV = RON_CONST1R * TANH((RON_QR0 - qr_val) * &
                     INV_RON_DELQR0) + RON_CONST2R
          END IF
      ELSE
          RONV = RN0_R
          SONV = RN0_S
      END IF

      ! Calculate Required equivalent reflectivity factors
      ZE_R = (RHOAIR * qr_val)**1.75D0
      ZE_R = FACTOR_R * ZE_R / RONV**.75D0

      ZE_S = (RHOAIR * qs_val)**1.75D0
      ZE_S = FACTORB_S * ZE_S / SONV**.75D0

      ! 4. Graupel Processing
      IF (has_graupel) THEN
          IF (QGRAUPEL(K) < 0.0) QGRAUPEL(K) = 0.0
          qg_val = QGRAUPEL(K)
          
          IF (use_var_intercepts) THEN
              GONV = GON_CONST
              IF (qg_val > R1) THEN
                  GONV = 2.38D0 * (PI_RHO_G / (RHOAIR * qg_val))**0.92D0
                  GONV = MAX(1.D4, MIN(GONV, GON_CONST))
              END IF
          ELSE
              GONV = RN0_G
          END IF
          
          ZE_G = (RHOAIR * qg_val)**1.75D0
          ZE_G = FACTORB_G * ZE_G / GONV**.75D0
      END IF

      ! 5. Hail Processing
      IF (has_hail) THEN
          IF (QHAIL(K) < 0.0) QHAIL(K) = 0.0
          qh_val = QHAIL(K)

          IF (use_var_intercepts) THEN
              HONV = HON_CONST
              IF (qh_val > R1) THEN
                  HONV = 2.38D0 * (PI_RHO_H / (RHOAIR * qh_val))**0.92D0
                  HONV = MAX(1.D4, MIN(HONV, HON_CONST))
              END IF
          ELSE
              HONV = RN0_H
          END IF
          
          ZE_H = (RHOAIR * qh_val)**1.75D0
          ZE_H = FACTORB_H * ZE_H / HONV**.75D0

          ! --- PSEUDO-MIE PLATEAU FOR WET HAIL ---
          IF (do_liq_hail .AND. tmk_val > CELKEL) THEN
              ZE_H = MIN(ZE_H, 1.0D7)
          END IF
      END IF

      ! Total Z_e
      Z_E = MAX(ZE_R + ZE_S + ZE_G + ZE_H, .001D0)

      ! Convert to dBZ
      DBZ(K) = 10.D0 * LOG10(Z_E)

  END DO

END FUNCTION compute_radar_reflectivity
end module GEOS_RadarMod
