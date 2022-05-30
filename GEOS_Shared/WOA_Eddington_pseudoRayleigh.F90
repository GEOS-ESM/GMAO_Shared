#include "MAPL_Generic.h"

module WOA_Eddington_pseudoRayleigh_Mod

  ! =============================================================
  ! Whole-of-atmosphere Eddington (Ed) plus pseudo-Rayleigh (pR)
  ! layer module, designed for quick changes to extern forcings:
  !   FdirDnTOA, mu0, adir, adif.
  ! Construct a Eddington_pseudoRayleigh_Optics then construct a
  ! WOA_Eddington_pseudoRayleigh with it.
  ! =============================================================
  ! Note: this module can be used to represent a DELTA-Eddington
  ! plus pseudo-Rayleigh model as well. You should do the delta-
  ! scaling outside and just construct the Eddington_pseudoRayleigh
  ! _Optics with the effective (delta-scaled) taue and ssae.
  ! =============================================================

  use ESMF
  use MAPL_Mod

  implicit none
  private

? ! constants
  double precision, parameter :: zero   =  0.0d0
  double precision, parameter :: one    =  1.0d0
  double precision, parameter :: two    =  2.0d0
  double precision, parameter :: three  =  3.0d0
  double precision, parameter :: four   =  4.0d0
  double precision, parameter :: five   =  5.0d0
  double precision, parameter :: seven  =  7.0d0
  double precision, parameter :: eight  =  8.0d0
  double precision, parameter :: nine   =  9.0d0
  double precision, parameter :: ten    = 10.0d0

? double precision, parameter :: half = one / two
  double precision, parameter :: third = one / three
  double precision, parameter :: quarter = one / four
  double precision, parameter :: eighth = one / eight
  double precision, parameter :: tenth = one / ten

? double precision, parameter :: three_sevenths = three / seven
  double precision, parameter :: five_eighths = five / eight
  double precision, parameter :: seven_ninths = seven / nine
  double precision, parameter :: fourteen_ninths = two * seven_ninths

  double precision, parameter :: pi = MAPL_PI_R8

  ! Constant for external use to avoid very low
  !   surface incoming direct flux and biased analyses
  double precision, parameter, public :: mu0_min = 0.1d0

  ! Constant for external use to avoid clash of
  !   -1/(any eigenvalue) too close to mu0
  double precision, parameter, public :: mu0_davoid = 1.0d-6

  ! constants to avoid extreme cases, which shouldnt occur
  ! in real cases and which can cause analyze_optics problems
  double precision, parameter, public :: ssae_min = 0.000001d0
  double precision, parameter, public :: ssae_max = 0.999999d0
  double precision, parameter, public :: g_Es_min = -third + 0.000001d0
  double precision, parameter, public :: g_Es_max = +third - 0.000001d0
  ! ranges
  double precision, parameter :: ssae_range = ssae_max - ssae_min
 r ouble precision, parameter :: g_Es_range = g_Es_max - g_Es_min

  ! fixed matricies
  double precision, parameter, dimension(2,2) :: matA = &
    reshape( [zero, one, one, zero], shape(matA), order=[2,1]) ! row-filled

  type, public :: Eddington_pseudoRayleigh_Optics

    ! Ed + pR optical props for a single layer

    private

    ! fundamental optical properties ...
    ! ==================================

    ! scaled optical thickness of WHOLE LAYER (0,)
    double precision, public :: taue

    ! scaled whole-media single scattering albedo (0,1)
    double precision, public :: ssae

    ! Eddington-star asymmetry parameter (-1/3,+1/3)
    double precision, public :: g_Es

    ! derived optical properties (see constructor) ...
    ! ================================================

    double precision :: ssaec
    double precision, dimension(2) :: vecla, vecT, vecnrla
    double precision, dimension(2,2) :: matK, matC

  end type Eddington_pseudoRayleigh_Optics

  type, public :: WOA_Eddington_pseudoRayleigh

    ! Whole-of-atmosphere Ed + pR layer
    ! (optical properties + external forcing)

    private

    ! optical state
    class(Eddington_pseudoRayleigh_Optics), allocatable, public :: optics

    ! solar forcing
    double precision, public :: FdirDnTOA  ! TOA downward solar flux (W/m2)
    double precision, public :: mu0        ! cos(solar_zenith_ang) in (0,1]
    ! (FdirDnTOA = mu0 * Fsol, where Fsol = TOA solar irrad normal to beam)
    ! (This module assumes and enforces zero downward diffuse flux at TOA)

    ! actual mu0 used since can be adjusted slightly to avoid an eigenvalue
    double precision, public :: mu0_used

    ! surface albedos 
    double precision, public :: adir  ! for direct solar beam flux
    double precision, public :: adif  ! for diffuse flux
    ! (in both cases, reflected flux is always diffuse)

    ! derived variables (see constructor)
    double precision :: Tsol, nrmu0_used, vecP(2), vecC(2)
 
  contains

    procedure :: compute_fluxes
    procedure :: analyze_optics

  end type WOA_Eddington_pseudoRayleigh

  interface Eddington_pseudoRayleigh_Optics
    module procedure Constructor_Eddington_pseudoRayleigh_Optics
  end interface

  interface WOA_Eddington_pseudoRayleigh
    module procedure Constructor_WOA_Eddington_pseudoRayleigh
  end interface

contains

  ! ============
  ! CONSTRUCTORS
  ! ============

  function Constructor_Eddington_pseudoRayleigh_Optics( &
    taue_, ssae_, g_Es_, rc) result (this)

    __Iam__('Constructor_Eddington_pseudoRayleigh_Optics')

    type(Eddington_pseudoRayleigh_Optics) :: this
    integer, intent(out), optional :: rc

    ! Eddington_pseudoRayleigh optical properties
    ! (see type for details)
    double precision, intent(in) :: taue_
    double precision, intent(in) :: ssae_
    double precision, intent(in) :: g_Es_

    ! locals
    double precision :: k2fac, ksquared

    ! rectify extreme optics
    this%taue = taue_
    this%ssae = min(max(ssae_,ssae_min),ssae_max)
    this%g_Es = min(max(g_Es_,g_Es_min),g_Es_max)

    ! warn on rectification
!   if (this%ssae /= ssae_) print *, trim(Iam) // ': warning: ssae rectified'
!   if (this%g_Es /= g_Es_) print *, trim(Iam) // ': warning: g_Es rectified'

    ! backstop: crash for unphysical or extreme optics
    _ASSERT(zero < this%taue,                       'require taue > 0')
    _ASSERT(zero < this%ssae .and. this%ssae < one, 'require ssae in (0,1)')

    ! Normally require an asymmetry factor g in (-1,+1).
    ! But a stricter (-1/3,+1/3) ensures the Eddington
    ! phase function is strictly positive at all angles.
    _ASSERT(-third < this%g_Es .and. this%g_Es < third, 'require g_Es in (-1/3,+1/3)')

    ! Solve for the homogeneous solution eigenvalues/vectors ...
    ! ==========================================================

    ! intermediates
    this%ssaec = one - this%ssae
    k2fac = three * (one - this%g_Es * this%ssae)
    ksquared = this%ssaec * k2fac
    _ASSERT(ksquared > zero, 'require ksquared > 0')

    ! eigenvalues (in increasing order)
    this%vecla(1) = -sqrt(ksquared)  ! -ve
    this%vecla(2) = -this%vecla(1)   ! +ve

    ! in case sqrt() underflows to zero
    _ASSERT(this%vecla(1) < zero, 'eigenvalues not distinct')

    ! mu0 equivalent of the eigenvalues
    this%vecnrla = -one / this%vecla

    ! eigen-transmissions
    this%vecT = exp(this%vecla * this%taue)

    ! eigenvectors
    this%matK(:,1) = eigenvec(this%vecla(1))
    this%matK(:,2) = eigenvec(this%vecla(2))
    ! (a 2x2 matrix with each COLUMN being an eigenvector)

    ! diagonal matrix needed to solve for particular solution later ...
    this%matC = zero
    this%matC(1,1) = this%ssaec
    this%matC(2,2) = k2fac

    RETURN_(ESMF_SUCCESS)

  contains

    function eigenvec(la) result(vecK)
      double precision, intent(in) :: la
      double precision :: vecK(2)
      vecK(1) = one
      vecK(2) = -this%ssaec / la
    end function eigenvec

  end function Constructor_Eddington_pseudoRayleigh_Optics

  function Constructor_WOA_Eddington_pseudoRayleigh( &
    optics, FdirDnTOA, mu0, adir, adif, rc) result (this)

    __Iam__('Constructor_WOA_Eddington_pseudoRayleigh')

    type(WOA_Eddington_pseudoRayleigh) :: this
    integer, intent(out), optional :: rc

    class(Eddington_pseudoRayleigh_Optics), intent(in) :: optics

    ! solar forcing (require mu0 in (0,1], i.e., strictly downward)
    double precision, intent(in) :: FdirDnTOA  ! TOA downward solar flux (W/m2)
    double precision, intent(in) :: mu0        ! cosine of solar zenith angle

    ! surface albedos 
    double precision, intent(in) :: adir  ! for direct solar beam flux
    double precision, intent(in) :: adif  ! for diffuse flux
    ! (in both cases, reflected flux is always diffuse)

    ! locals
    double precision :: Qsol, vec(2), mat(2,2), matQ(2,2)
    double precision :: opadif, omadif, fac, delta
    integer :: ipiv(2), info, n
    logical :: mu0_rectified

    ! check inputs
    _ASSERT(zero < FdirDnTOA,               'require FdirDnTOA > 0')
    _ASSERT(zero < mu0 .and. mu0 <= one,    'require mu0 in (0,1]')
    _ASSERT(zero <= adir .and. adir <= one, 'require adir in [0,1]')
    _ASSERT(zero <= adif .and. adif <= one, 'require adif in [0,1]')

    ! copy optical properties
    this%optics = optics

    ! copy forcings
    this%FdirDnTOA = FdirDnTOA
    this%mu0       = mu0
    this%adir      = adir
    this%adif      = adif

    ! avoid mu0 too close to -1/(any eigenvalue)
    mu0_rectified = .false.
    do n = 1,2
      delta = this%mu0 - optics%vecnrla(n)
      if (abs(delta) < mu0_davoid) then
        this%mu0_used = optics%vecnrla(n) + 1.01d0 * sign(mu0_davoid, delta)
        this%mu0_used = min(max(this%mu0_used,zero),one)
        mu0_rectified = .true.
        exit ! assume eigenvalues well separated
      end if
    end do
    if (mu0_rectified) then
      ! backstop
      _ASSERT(minval(abs(this%mu0_used - optics%vecnrla)) >= mu0_davoid, 'eigenvalue resonance')
    else
      ! default
      this%mu0_used = this%mu0
    end if

    ! Note: we do not adjust FdirDnTOA for mu0_used,
    ! but keep the downward energy flux at TOA unchanged.

    ! Solve for the particular solution ...
    ! =====================================

    ! direct beam transmission of whole layer (i.e., to SFC)
    this%nrmu0_used = -one / this%mu0_used
    this%Tsol = exp(optics%taue * this%nrmu0_used)

    ! intermediate variables
    Qsol = this%FdirDnTOA / pi  ! Qsol = "mu0 * Ssol"
    vec = half * Qsol * optics%ssae * &
      [one, three * optics%g_Es * this%mu0_used]

    ! solve for vecP: (mu0_used * matC - matA) * vecP = vec
    ! using LAPACK routine dgesv
    mat = this%mu0_used * optics%matC - matA
    call dgesv(2, 1, mat, 2, ipiv, vec, 2, info)
    _ASSERT(info == 0, 'dgesv fail in solve for vecP')
    this%vecP = vec

    ! Apply boundary conditions ...
    ! =============================

    ! solve for vecC = [C_1, C_2]
    omadif  = one - this%adif
    opadif  = one + this%adif
    matQ = reshape([half, one, half*omadif, -opadif], [2,2], order=[2,1]) ! row-filled
    mat = matmul(matQ,optics%matK)
    mat(2,:) = mat(2,:) * optics%vecT
    vec = -matmul(matQ,this%vecP)
    fac = this%adir * Qsol
    vec(2) = (vec(2) + fac) * this%Tsol
    ! solve mat * vecC = vec
    call dgesv(2, 1, mat, 2, ipiv, vec, 2, info)
    _ASSERT(info == 0, 'dgesv fail in solve for vecC')
    this%vecC = vec

    RETURN_(ESMF_SUCCESS)

  end function Constructor_WOA_Eddington_pseudoRayleigh

  ! =====================
  ! TYPE-BOUND PROCEDURES
  ! =====================

  subroutine compute_fluxes(this, &
    Rsol, FdirDn, FdifDn, FdifUp, rc)

    ! compute fluxes in a WOA Eddington + Rayleigh layer

    __Iam__('compute_fluxes')

    class(WOA_Eddington_pseudoRayleigh) :: this
    integer, intent(out), optional :: rc
    
    ! Rsol is the input within-layer-position parameter [0,1].
    ! This is the RATIO of the ATTENUATION of the direct solar beam
    ! at the position in question to that of the whole layer (to SFC).
    ! Thus Rsol in [0,1], with 0=>TOA, 1=>SFC, and, e.g., Rsol=0.25,
    ! is the position at which only one quarter of the attenuation of
    ! the whole layer occurs.
    double precision, intent(in ) :: Rsol

    ! vertical fluxes (W/m2)
    double precision, intent(out) :: FdirDn  ! downward solar beam flux
    double precision, intent(out) :: FdifDn  ! downward diffuse flux
    double precision, intent(out) :: FdifUp  ! upward diffuse flux

    ! locals
    double precision :: Tsol, taue, ID(2)

    ! must select position somewhere in WOA layer
    _ASSERT(zero <= Rsol .and. Rsol <= one, 'require Rsol in [0,1]')

    ! formally Rsol = (one - Tsol) / (one - this%Tsol)
    Tsol = one - Rsol * (one - this%Tsol)
    _ASSERT(zero <= Tsol .and. Tsol <= one, 'require Tsol in [0,1]')
    FdirDn = this%FdirDnTOA * Tsol
    taue = -log(Tsol) * this%mu0_used
    call Idif(taue, ID)
    call Fdif(ID, FdifDn, FdifUp)

    RETURN_(ESMF_SUCCESS)

  contains

    ! azimuthally averaged diffuse radiance components
    subroutine Idif(taue, ID)

      double precision, intent(in ) :: taue
      double precision, intent(out) :: ID(2)

      ID = exp(this%optics%vecla * taue) * this%vecC
      ID = matmul(this%optics%matK, ID) + this%vecP * exp(this%nrmu0_used * taue)

    end subroutine Idif

    ! diffuse vertical fluxes
    subroutine Fdif(ID, FdifDn, FdifUp)

      double precision, intent(in ) :: ID(2)
      double precision, intent(out) :: FdifDn, FdifUp

      double precision :: fac

      fac = half * ID(1)
      FdifDn = pi * (fac + ID(2))
      FdifUp = pi * (fac - ID(2))

    end subroutine Fdif
  
  end subroutine compute_fluxes

  function analyze_optics(bkg,    &
    nlay, FdifUp, FtotDn, FdirDn, &
    lmdif_info, frmse,            &
    lmdif_epsfcn,                 &
    lmdif_nprint,                 &
    rc) result (ana)

    ! Analyzes a Eddington_pseudoRayleigh_Optics consistent with the provided
    ! vertical fluxes. The passed-object dummy argument is called bkg because
    ! it serves as the background state (initial condition) for the analysis.
    ! It also provides the forcings, FdirDnTOA, mu0, adir and adif, which are
    ! held fixed in the analysis. The returned ana is a WOA_Eddington_pseudo-
    ! Rayleigh with these same forcings, but with the analyzed Eddington_pseudo-
    ! Rayleigh_Optics, and should yield computed fluxes close to those provided.
    !   This version uses non-linear least squares (minpack's lmdif) to constrain
    ! ssae and g_Es. taue is constrained directly from the atmospheric atten-
    ! uation of the direct beam.
    !   For data we use FdifUp(1:nlay) and FtotDn(2:nlay+1), with the corres-
    ! ponding FdirDn used to place the levels on an appropriate Rsol scale.
    ! FtotDn(1) is assumed to be FdirDn(1), consistent with our model's assum-
    ! ption of zero incoming diffuse. Likewise, FdifUp(nlay+1) is assumed by
    ! our model to be a diagnostic based on the reflected surface fluxes and
    ! given albedos. This is why FtotDn(1) and FdifUp(nlay+1) are not used as
    ! constraints.

    __Iam__('analyze_optics')

    class(WOA_Eddington_pseudoRayleigh) :: bkg
    integer, intent(out), optional :: rc

    type(WOA_Eddington_pseudoRayleigh) :: ana

    ! number of model layers
    integer, intent(in) :: nlay

    ! vertical fluxes (W/m2) to constrain against [nlay+1]
    double precision, intent(in) :: FdifUp(nlay+1)  ! upwelling diffuse flux
    double precision, intent(in) :: FtotDn(nlay+1)  ! downwelling total flux
    double precision, intent(in) :: FdirDn(nlay+1)  ! solar beam flux

    ! report how lmdif_driver did
    integer, intent(out) :: lmdif_info      ! info from lmdif solver
    double precision, intent(out) :: frmse  ! final RMS error in fluxes

    ! used in determining step length for frwd-dif approx
    double precision, intent(in), optional :: lmdif_epsfcn

    ! level of reporting
    integer, intent(in), optional :: lmdif_nprint

    ! locals
    type(Eddington_pseudoRayleigh_Optics) :: optics
    double precision :: Rsol(nlay+1), taue, ssae, g_Es
    double precision :: FdifUpSFC_expected
    character(len=4) :: istr

    ! lmdif_driver locals
    ! there are 2*nlay constraints, as per subroutine header
    ! there are 2 parameters to constrain: ssae and g_Es.
    !   (taue is constrained directly)
    double precision :: x(2), fvec(2*nlay), epsfcn_
    integer :: nprint_, k

    ! the default output in case of error before lmdif call
    lmdif_info = -9

    ! number of constraints should equal or exceed number of parameters to constrain
    _ASSERT(size(fvec) >= size(x), 'not enough constraints')

    ! check on assumptions:
    ! (1) no downward diffuse at TOA
    _ASSERT(FtotDn(1) == FdirDn(1), 'FtotDn @ TOA must be all direct')
    ! (2) direct at TOA should be in background (forcing) state
    _ASSERT(FdirDn(1) == bkg%FdirDnTOA, 'FDirDn @ TOA not consistent with bkg')
    ! (3) check using correct albedos
    FdifUpSFC_expected = &
      bkg%adir * FdirDn(nlay+1) + bkg%adif * (FtotDn(nlay+1)-FdirDn(nlay+1))
    _ASSERT(abs(FdifUp(nlay+1) - FdifUpSFC_expected) < 1.d-6, 'FdifUp @ SFC not consistent')

    ! further checks:
    ! (1) total downwelling fluxes should be positive at all levels
    _ASSERT(all( FtotDn > zero ), 'require all FtotDn > 0')
    ! (2) exclude also exponential underflow to zero for beam attenuation
    ! (partial redundancy: have already ensured FdirDn(1) == FtotDn(1) > 0)
    _ASSERT(all( FdirDn > zero ), 'require all FdirDn > 0')
    ! (3) check that beam never increases downward
    _ASSERT(all( FdirDn(1:nlay) >= FdirDn(2:nlay+1) ), 'FdirDn should only attenuate')
    ! (4) upward fluxes should be positive at all levels,
    ! but allow zero up at surface since albedos can be zero
    _ASSERT(     FdifUp(nlay+1) >= zero  , 'require FdifUp @ SFC >= 0')
    _ASSERT(all( FdifUp(1:nlay) >  zero ), 'require FdifUp above SFC > 0')
    ! (5) outgoing shouldn't exceed incoming
    _ASSERT( FdifUp(1) <= FtotDn(1), 'outgoing should not exceed incoming @ TOA' )

    ! set taue directly from solar transmission
    taue = -bkg%mu0 * log(FdirDn(nlay+1) / bkg%FdirDnTOA)
    _ASSERT(taue > zero, 'require taue > 0')

    ! set Rsol of data points
    ! remember that bkg FdirDnTOA held fixed
    Rsol(1)      = zero
    Rsol(2:nlay) = (bkg%FdirDnTOA - FdirDn(2:nlay)) / (bkg%FdirDnTOA - FdirDn(nlay+1))
    Rsol(nlay+1) = one

    ! normalize to [0,1] on rectified range
    ssae = (bkg%optics%ssae - ssae_min) / ssae_range
    g_Es = (bkg%optics%g_Es - g_Es_min) / g_Es_range

    ! transformation requires strictly within rectified limits
    _ASSERT(zero < ssae .and. ssae < one, 'require normalized ssae in (0,1)')
    _ASSERT(zero < g_Es .and. g_Es < one, 'require normalized g_Es in (0,1)')

    ! project initial conditions to the nD real space of solver
    ! Note: (0,1) -> (-Inf,+Inf) using inverse logistic
    x(1) = -log(one / ssae - one)
    x(2) = -log(one / g_Es - one)

    ! non-linear least-squared solution using MINPACK's lmdif
    !   (which uses an approximate non-analytic Jacobian).
    ! The direct call is via simplified lmdif_driver, which
    ! is very similar to MINPACK's simplified lmdif1 driver,
    ! but with modified parameters.
    if (present(lmdif_epsfcn)) then
      epsfcn_ = lmdif_epsfcn
    else
      epsfcn_ = zero	! what lmdif1 uses
!h    epsfcn_ = 1.0d-3	! seem to be fast and accurate
    endif
    if (present(lmdif_nprint)) then
      nprint_ = lmdif_nprint
    else
      nprint_ = 0	! no reporting
    endif
    call lmdif_driver(lmdif_residuals,size(fvec),size(x),x,fvec,lmdif_info,epsfcn_,nprint_)
    if (lmdif_info < 1 .or. lmdif_info > 3) then
      write(istr,'(i0)') lmdif_info
      _ASSERT(.FALSE., 'non-linear solver failing with info = ' // istr)
    end if
    frmse = sqrt(sum(fvec**2)/size(fvec))

    ! extract solution from nD real space
    ssae = one / (one + exp(-x(1)))
    g_Es = one / (one + exp(-x(2)))
    ssae = ssae * ssae_range + ssae_min
    g_Es = g_Es * g_Es_range + g_Es_min

    ! reset Eddington_pseudoRayleigh optics and WOA layer state
    optics = Eddington_pseudoRayleigh_Optics(taue, ssae, g_Es, __RC__)
    ana = WOA_Eddington_pseudoRayleigh(optics, bkg%FdirDnTOA, bkg%mu0, bkg%adir, bkg%adif,__RC__)

    RETURN_(ESMF_SUCCESS)

  contains

    subroutine lmdif_driver(fcn,m,n,x,fvec,info,epsfcn,nprint)

      integer,          intent(in   ) :: m, n
      double precision, intent(inout) :: x(n)
      double precision, intent(out  ) :: fvec(m)
      integer,          intent(out  ) :: info
      double precision, intent(in   ) :: epsfcn
      integer,          intent(in   ) :: nprint
      external fcn

      ! epsfcn: used in determining a suitable step length
      !   for the forward-difference approximation

      ! nprint: if /= 0, report progress on first, last,
      !   and every nprint iterations

      integer :: iwa(n)
      double precision :: xtol
      double precision :: wa(m*n+5*n+m)

      ! dimension of solution space
      _ASSERT(n > 0, 'invalid solution space dimension')

      ! precision in x space
      xtol = sqrt(epsilon(one))
      _ASSERT(xtol >= zero, 'lmdif1 requires non-negative tolerance')

      ! for now try lmdif1
      call lmdif1(fcn,m,n,x,fvec,xtol,info,iwa,wa,size(wa))

    end subroutine lmdif_driver

    subroutine lmdif_residuals(m, n, x, fvec, iflag)

      ! Flux residuals @ x = [ssae, g_Es], 
      !   but transformed to 2D real space.
      ! Should set iflag negative on any error.

      integer :: m, n, k, iflag, status
      double precision :: x(n), fvec(m)

      double precision :: x_ssae, x_g_Es
      type(Eddington_pseudoRayleigh_Optics) :: x_opt
      type(WOA_Eddington_pseudoRayleigh) :: x_WOA
      double precision :: x_FdirDn, x_FdifDn, x_FdifUp

      ! convert to physical variables
      x_ssae = one / (one + exp(-x(1)))
      x_g_Es = one / (one + exp(-x(2)))
      x_ssae = x_ssae * ssae_range + ssae_min
      x_g_Es = x_g_Es * g_Es_range + g_Es_min

      if (iflag == 0) then
        print *, 'x:', x
        print *, 'x_:', x_ssae, x_g_Es
      end if

      ! make trial optics
      x_opt = Eddington_pseudoRayleigh_Optics( &
        taue, x_ssae, x_g_Es, rc=status)
      if (status == ESMF_FAILURE) then
        iflag = -1
        return
      end if

      ! solve with provided fixed (bkg) forcings
      x_WOA = WOA_Eddington_pseudoRayleigh(x_opt, &
        bkg%FdirDnTOA, bkg%mu0, bkg%adir, bkg%adif, rc=status)
      if (status == ESMF_FAILURE) then
        iflag = -2
        return
      end if

      ! solve for flux residuals
! can rework to do a single set of evaluations for greater efficiency TODO PMN
      ! (1) upwelling fluxes
      do k = 1, nlay
        call x_WOA%compute_fluxes(Rsol(k), x_FdirDn, x_FdifDn, x_FdifUp, rc=status)
        if (status == ESMF_FAILURE) then
          iflag = -3
          return
        end if
        fvec(k) = FdifUp(k) - x_FdifUp
      end do
      ! (2) downwelling fluxes
      do k = 1, nlay
        call x_WOA%compute_fluxes(Rsol(k+1), x_FdirDn, x_FdifDn, x_FdifUp, rc=status)
        if (status == ESMF_FAILURE) then
          iflag = -4
          return
        end if
        fvec(nlay+k) = FtotDn(k+1) - (x_FdirDn + x_FdifDn)
      end do

      if (iflag == 0) then
        print *, 'fvec:', fvec
        print *
      end if

    end subroutine lmdif_residuals

  end function analyze_optics

end module WOA_Eddington_pseudoRayleigh_Mod
