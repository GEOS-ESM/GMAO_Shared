module svecs_cf

#include "MAPL_Generic.h"

use ESMF
use MAPL_mod

implicit none
public


! !DESCRIPTION: Configuration and parameters for singular vector calculation.
!
! !REVISION HISTORY:
!
!  28Apr2017  Holdaway   Ported from g5pert environment.

!Constants
real(8), parameter :: zero = 0.0_8
real(8), parameter :: half = 0.5_8
real(8), parameter :: one  = 1.0_8
real(8), parameter :: d2r = MAPL_PI / 180.0_8
real(8), parameter :: r2d = 180.0_8 / MAPL_PI

!Grid info
character(len=ESMF_MAXSTR)          :: gridtype
integer                             :: npx, npy, npz
integer                             :: im, jm, lm, nq
integer                             :: ifirst, ilast, jfirst, jlast
real(8), allocatable, dimension(:)  :: lambda, theta
real(8), allocatable, dimension(:)  :: lats, lons
integer, allocatable, dimension(:)  :: levs
real(8)                             :: ptop

!Screen for norm
real(8), allocatable, dimension(:,:,:) :: coefftpI         !Projection coefficients start time
real(8), allocatable, dimension(:,:,:) :: coefftpF         !Projection coefficients final time

!Area/level weighting
real(8), allocatable, dimension(:,:,:) :: gridweightI
real(8), allocatable, dimension(:,:,:) :: gridweightF

!Pointers to the trajectory bundle
type trajPointer
 real, pointer, dimension(:,:,:)   :: u,v,t,delp,sphu
end type trajPointer

type pertState
 real(8), allocatable, dimension(:,:,:)   :: u,v,t,delp,sphu
 real(8), allocatable, dimension(:,:,:)   :: qitot,qltot,ozone
end type pertState

type pertStateR4
 real(4), allocatable, dimension(:,:,:)   :: u,v,t,delp,sphu
 real(4), allocatable, dimension(:,:,:)   :: qitot,qltot,ozone
end type pertStateR4

integer, save                    ::   isolve        ! selects eignen solver: 1=ARPACK, 2=NAG
logical, save                    ::   rvec          ! calculate eigenvectors
logical, save                    ::   eigchk        ! check on accuracy of eigen-decomposition
logical, save                    ::   propsvec      ! determines whether or not to evolve svecs
character(len=ESMF_MAXSTR), save :: svecnormI       ! initial state norm
character(len=ESMF_MAXSTR), save :: svecnormF       ! final state norm for eigenvectors
                                                    !   L2  - for L2-norm
                                                    !   KE  - for Kinetic Energy norm
                                                    !   TE  - for Total   Energy norm
                                                    !   Pa  - for analysis error cov-based norm
                                                    !   Da  - for analysis error variances norm
                                                    ! local projection operator box limits
logical, save                    :: lclproj         !   .t. when local proj is applied
real(8), save                    :: projlonI(2)     !   longitudes
real(8), save                    :: projlatI(2)     !   latitudes
integer, save                    :: projlevI(2)     !   vertical levels
real(8), save                    :: projlonF(2)     !   longitudes
real(8), save                    :: projlatF(2)     !   latitudes
integer, save                    :: projlevF(2)     !   vertical levels
real(8), save                    :: projlongridI(2) !   longitudes but on grid
real(8), save                    :: projlatgridI(2) !   latitudes but on grid
real(8), save                    :: projlongridF(2) !   longitudes but on grid
real(8), save                    :: projlatgridF(2) !   latitudes but on grid
integer, save                    :: test_norm       !   levs of norm tests: 0, 1, or 2
integer, save                    :: perc_var        ! percentage variance captured in innov cov
integer, save                    :: maxitr          ! maximum no. of Lanczos iterations per job
integer, save                    :: maxd = 1        ! default maximum no. Lanczos iterations per job
integer, save                    :: ncalls          ! counter for lanczos iterations
integer, save                    :: lanunit         ! unit number for lanczos restart data
character(len=ESMF_MAXSTR), save :: pflabel         ! perturbation file tag name
real(8), save                    :: eps_eer         ! eps from Ehrendorfer, Errico and Raeder (1999)
                                                    ! this controls the extent to which q
                                                    ! influences the "wet-energy" norm

logical, save                    :: cp_int_fv = .false. ! Keep internal checkpoints for FV3 adjoint in memory

! 1. ARPACK-specific parameters
! -----------------------------
integer, save                    ::   nevd = 1      ! default number of eigenvalues/vectors
integer, save                    ::   ncvd = 2      ! default number of Lanczos basis vectors
character*1, save                ::   bmat	    ! define eigen-problem type in ARPACK sense
integer, save                    ::   iparam(11)    ! array of specific options for ARPACK
real(8), save                    ::   tol           ! relative accuracy of eigenvalues

! 2. NAG-specific parameters
! --------------------------
integer, save                    :: lanmaxd = 1     ! default total number of Lanczos iterations
real(8), save                    ::   kappa         ! relative accuracy of eigenvalues

! 3. CNOP-specific parameters
! ---------------------------
integer, save                    :: spg_mfls  = 1000! number of function evaluations used in SPG line search
integer, save                    :: spg_maxfc = 2000! max number of function evaluations
integer, save                    :: spg_miterd= 100 ! default max number of iterations
real(8), save                    :: spg_einf  = 0.0 ! SPG inf norm criterium
logical, save                    :: spg_verb  = .true.
integer, save                    :: spg_miter       ! max number of iterations
real(8), save                    :: cnop_sigma      ! magnitude of norm of growing perturbation;
                                                    ! typical LSV total energy
real(8), save                    :: cnop_tol

interface allocate_xpert
  module procedure allocate_xpert_r4
  module procedure allocate_xpert_r8
end interface allocate_xpert

contains

subroutine SetSvecs(SVFILE)

 implicit none

 character(len=ESMF_MAXSTR) :: SVFILE
 type (ESMF_Config)         :: cf
 integer                    :: rc

 character(len=ESMF_MAXSTR)   :: Iam="SetSvecs"

 cf = ESMF_ConfigCreate(rc=rc)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = rc )

! isolve
 call ESMF_ConfigGetAttribute( cf, isolve, label='eigensolver_package:', default=2, rc = rc)
 if (MAPL_AM_I_ROOT()) then
    if (isolve == 1) then
       print*, 'Eigensolver: ',isolve, ' ARPACK'
    elseif (isolve == 2) then
       print*, 'Eigensolver: ',isolve, ' NAG'
    elseif (isolve == 3) then
       print*, 'Eigensolver: ',isolve, ' Will solve CNOP problem'
    endif
 endif

! svecnormF
 call ESMF_ConfigGetAttribute( cf, svecnormF, label='final_svec_norm:',default= 'te', rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Final time svec norm: ', trim(svecnormF)

!svecnormI
 call ESMF_ConfigGetAttribute( cf, svecnormI, label='start_svec_norm:', default=svecnormF, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Initial time svec norm: ', trim(svecnormI)

!test_norm
 call ESMF_ConfigGetAttribute( cf, test_norm, label='do_norm_test:', default=0, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Test norm value: ', test_norm

!rvec
 call ESMF_ConfigGetAttribute( cf, rvec, label='calculate_eigenvectors:', default=.true., rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Calculate eigenvectors: ', rvec

!propsvec
 call ESMF_ConfigGetAttribute( cf, propsvec, label='evolve_svec:', default=.false., rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Evolve svecs: ', propsvec 

!maxitr
 call ESMF_ConfigGetAttribute( cf, maxitr, label='maximum_iterations_per_job:', default=1, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Max iterations for this job: ',maxitr

!eps_eer
 if (trim(svecnormF)=='we' .or. trim(svecnormI)=='we') then
    call ESMF_ConfigGetAttribute( cf, eps_eer, label='ehrendorfer_errico_raedder_eps:', default=1.0_8, rc = rc)
    if (MAPL_AM_I_ROOT()) print*, 'Ehrendorfer, Errico, and Raeder eps: ',eps_eer 
 endif

!Hold internal FV checkpoints in memory 
 call ESMF_ConfigGetAttribute( cf, cp_int_fv, label='hold_fv_checkpoints:', default=.false., rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Hold FV checkpoints in memory: ', cp_int_fv

 VERIFY_(rc)

endsubroutine SetSvecs


subroutine SetSvecs_ARPACK(SVFILE,nev,ncv,which)

 implicit none

 character(len=ESMF_MAXSTR) :: SVFILE
 type (ESMF_Config)         :: cf
 integer                    :: rc

 character(len=ESMF_MAXSTR)   :: Iam="SetSvecs_ARPACK"

 integer,     intent(out) :: nev     ! number of eigenvalues/eigenvectors to calculate
 integer,     intent(out) :: ncv     ! number of Lanczos basis vectors
 character*2, intent(out) :: which

 integer :: j, maxupd

 cf = ESMF_ConfigCreate(rc=rc)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = rc )

!nev
 call ESMF_ConfigGetAttribute( cf, nev, label='number_eigenvectors:', default=nevd, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Number of singular values(nev): ',nev

!ncv
 call ESMF_ConfigGetAttribute( cf, ncv, label='number_lanczos_basis_vectors:', default=ncvd, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Number of Lanczos basis vectors (ncv): ', ncv

!which
 call ESMF_ConfigGetAttribute( cf, which, label='which_eigenvalues:', default='LM', rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Which eigenvalues: ', which

!tol
 call ESMF_ConfigGetAttribute( cf, tol, label='eigenvalue_relative_accuracy:', default=0.05_8, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Tolerance: ', tol

!maxupd
 call ESMF_ConfigGetAttribute( cf, maxupd, label='maximum_arnoldi_iterations:', default=10, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Maximum allowed number of Arnoldi update iterations: ',maxupd

!eigchk
 call ESMF_ConfigGetAttribute( cf, eigchk, label='eigen_decomposition_accuracy:', default=.false., rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Check on accuracy of eigen-decomposition: ', eigchk

!now do some setting up:
 do j = 1, 11
    iparam(j) = 0
 end do
 iparam(1) = 1		!  shifting option; set to exact shifts
 iparam(3) = maxupd	!  maximum number of Arnoldi update iterations
 if ( trim(svecnormI) == 'ke' .or. trim(svecnormI) == 'te' .or. trim(svecnormI) == 'we') then
    iparam(7) =  1	!  type of eigen-problem: standard eigen-problem
    bmat      = 'I'
 else
    if ( MAPL_AM_I_ROOT() ) print*, ': No such scheme implemented. Aborting ...'
    stop
 end if

 VERIFY_(rc)

endsubroutine SetSvecs_ARPACK


subroutine SetSvecs_NAG(SVFILE,lanmax)

 implicit none

 character(len=ESMF_MAXSTR) :: SVFILE
 type (ESMF_Config)         :: cf
 integer                    :: rc
 integer, intent(out)       :: lanmax

 character(len=ESMF_MAXSTR)   :: Iam="SetSvecs_NAG"

 cf = ESMF_ConfigCreate(rc=rc)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = rc )

!lanmax
 call ESMF_ConfigGetAttribute( cf, lanmax, label='number_lanczos_iterations:', default=lanmaxd, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Number of Lanczos steps (lanmax): ',lanmax

!kappa
 call ESMF_ConfigGetAttribute( cf, kappa, label='eigenvector_accuracy:', default=0.05_8, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Eigenvector accuracy: ',kappa

 VERIFY_(rc)

endsubroutine SetSvecs_NAG


subroutine SetSvecs_CNOP(SVFILE,what,ncv)

 implicit none

 character(len=ESMF_MAXSTR) :: SVFILE
 type (ESMF_Config)         :: cf
 integer                    :: rc

 character(len=ESMF_MAXSTR)   :: Iam="SetSvecs_CNOP"

 integer,          intent(out) :: ncv  ! number of CNOPs to get
 character(len=4), intent(out) :: what ! dummy to trick interface

 cf = ESMF_ConfigCreate(rc=rc)
 call ESMF_ConfigLoadFile( cf, SVFILE, rc = rc )

!what
 what = 'dummy'

!ncv
 call ESMF_ConfigGetAttribute( cf, ncv, label='ncv:', default=ncvd, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Number of CNOPs to get: ',ncv

!spg_miter
 call ESMF_ConfigGetAttribute( cf, spg_miter, label='number_spg_iterations:', default=spg_miterd, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'Number of SPG max iterations (spg_miter): ',spg_miter

!cnop_tol
 call ESMF_ConfigGetAttribute( cf, cnop_tol, label='cnop_tolerance:', default=0.32_8, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'CNOP norm tolerance: ',cnop_tol

!cnop_sigma
 call ESMF_ConfigGetAttribute( cf, cnop_sigma, label='cnop_norm_value:', default=28._8, rc = rc)
 if (MAPL_AM_I_ROOT()) print*, 'CNOP target norm value: ',cnop_sigma

 VERIFY_(rc)

endsubroutine SetSvecs_CNOP

subroutine allocate_config()

 implicit none

 if (.not.allocated(lambda)) allocate(lambda(im)) 
 if (.not.allocated(theta)) allocate(theta(jm))
 if (.not.allocated(lats)) allocate(lats(jm))
 if (.not.allocated(lons)) allocate(lons(im))
 if (.not.allocated(levs)) allocate(levs(lm))
 if (.not.allocated(coefftpI)) allocate(coefftpI(im,jm,lm))
 if (.not.allocated(coefftpF)) allocate(coefftpF(im,jm,lm))
 if (.not.allocated(gridweightI)) allocate(gridweightI(im,jm,lm))
 if (.not.allocated(gridweightF)) allocate(gridweightF(im,jm,lm))

end subroutine allocate_config

subroutine deallocate_config()

 implicit none

 if (allocated(lambda)) deallocate(lambda) 
 if (allocated(theta)) deallocate(theta)
 if (allocated(lats)) deallocate(lats)
 if (allocated(lons)) deallocate(lons)
 if (allocated(levs)) deallocate(levs)
 if (allocated(coefftpI)) deallocate(coefftpI)
 if (allocated(coefftpF)) deallocate(coefftpF)
 if (allocated(gridweightI)) deallocate(gridweightI)
 if (allocated(gridweightF)) deallocate(gridweightF)

end subroutine deallocate_config

subroutine allocate_xpert_r8(xpert)

 implicit none

 type(pertState), intent(inout) :: xpert

 allocate(xpert%u(im,jm,lm)) 
 allocate(xpert%v(im,jm,lm)) 
 allocate(xpert%t(im,jm,lm)) 
 allocate(xpert%delp(im,jm,lm)) 
 allocate(xpert%sphu(im,jm,lm)) 
 allocate(xpert%qitot(im,jm,lm)) 
 allocate(xpert%qltot(im,jm,lm)) 
 allocate(xpert%ozone(im,jm,lm)) 

 xpert%u = 0.0_8
 xpert%v = 0.0_8
 xpert%t = 0.0_8
 xpert%delp = 0.0_8
 xpert%sphu = 0.0_8
 xpert%qitot = 0.0_8
 xpert%qltot = 0.0_8
 xpert%ozone = 0.0_8

end subroutine allocate_xpert_r8

subroutine allocate_xpert_r4(xpert)

 implicit none

 type(pertStateR4), intent(inout) :: xpert

 allocate(xpert%u(im,jm,lm)) 
 allocate(xpert%v(im,jm,lm)) 
 allocate(xpert%t(im,jm,lm)) 
 allocate(xpert%delp(im,jm,lm)) 
 allocate(xpert%sphu(im,jm,lm)) 
 allocate(xpert%qitot(im,jm,lm)) 
 allocate(xpert%qltot(im,jm,lm)) 
 allocate(xpert%ozone(im,jm,lm)) 

 xpert%u = 0.0_4
 xpert%v = 0.0_4
 xpert%t = 0.0_4
 xpert%delp = 0.0_4
 xpert%sphu = 0.0_4
 xpert%qitot = 0.0_4
 xpert%qltot = 0.0_4
 xpert%ozone = 0.0_4

end subroutine allocate_xpert_r4

subroutine deallocate_xpert(xpert)

 implicit none

 type(pertState), intent(inout) :: xpert

 deallocate(xpert%u) 
 deallocate(xpert%v) 
 deallocate(xpert%t) 
 deallocate(xpert%delp) 
 deallocate(xpert%sphu) 

end subroutine deallocate_xpert

end module svecs_cf
