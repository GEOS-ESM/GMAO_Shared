!  $Id$

#include "MAPL_Generic.h"

#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData
#define PTOP 1

!=============================================================================
!BOP

! !MODULE: 

! GEOS\_PertSharedMod -- A container module for constants and types used
!   Pert model.

! !INTERFACE:

module GEOS_PertSharedMod

! !USES:

  use ESMF
  use MAPL_Mod

  implicit none
! private

  public:: phase		! thisPhase = phase(adjoint=.false.,import=.true.,export=.true.)
  public:: phase_opAdjoint	! isAdjoint = phase_opAdjoint(phase)
  public:: phase_addImport	! addImport = phase_addImport(phase)
  public:: phase_getExport	! setExport = phase_getExport(phase)
  public:: DefVertGrid
  public:: IceFraction
  public:: GetShapiroCoeff

! public:: T_1DVar
! public:: T_2DVar
! public:: T_3DVar
!EOP

  integer, save :: DynTestCase 

  integer, save :: ObserverMode

  !Component switches
  integer, save :: DO_DYNAMICS
  integer, save :: DO_BL_PHYS
  integer, save :: DO_MOIST_PHYS
  integer, save :: DO_RAD_PHYS
  integer, save :: DO_GWD_PHYS
  integer, save :: DO_PHYS_D2A

  integer, save :: TrT_ver
  real   , save :: TrT_h0
  real   , save :: TrT_R
  real   , save :: TrT_latc
  real   , save :: TrT_lonc

  integer, parameter :: TLMPhase    = 1
  integer, parameter :: ADJPhase    = 2

  logical, save :: BothPhases = .false.

  integer, parameter :: TLMPhase_yIMPORT_yEXPORT = 1	! 000 +1
  integer, parameter :: ADJPhase_yIMPORT_yEXPORT = 2	! 001 +1
  integer, parameter :: TLMPhase_nIMPORT_yEXPORT = 3	! 010 +1
  integer, parameter :: ADJPhase_nIMPORT_yEXPORT = 4	! 011 +1
  integer, parameter :: TLMPhase_yIMPORT_nEXPORT = 5	! 100 +1
  integer, parameter :: ADJPhase_yIMPORT_nEXPORT = 6	! 101 +1
  integer, parameter :: TLMPhase_nIMPORT_nEXPORT = 7	! 110 +1
  integer, parameter :: ADJPhase_nIMPORT_nEXPORT = 8	! 111 +1

! The type for the private internal state. All fields here
!   are on the cubed sphere grid, are real*4, and nothing
!   is ghosted. 

  integer, save :: NUM_GCM3Dvars
  integer, save :: NUM_GCM2Dvars
  integer, parameter :: NUM_GCM3DvarsRead = 7
  integer, parameter :: NUM_GCM3DvarsReadMoist = 8
  integer, parameter :: NUM_GSIvars   = 8
  integer, parameter :: NameSize    = 8


  !PERTURBATION TRAJECTORY, RANK 3 VARIABLES
  character(len=NameSize), parameter :: GSIvars(NUM_GSIvars) = &
       (/'U ','V ','TV','DP','QV','QI','QL','O3'/)

  !GCM TRAJECTORY, RANK 3 VARIABLES
  character(len=NameSize), parameter :: GCM3Dvars(11) = &
       [ character(NameSize) :: 'U ','V ','PT','DP','QV','QI','QL','O3','QLS','QCN','CFCN' ]

  !Variables that will be read when dry
  character(len=NameSize), parameter :: GCM3DvarsRead(NUM_GCM3DvarsRead) = &
       [ character(NameSize) :: 'U ','V ','PT',     'QV','QI','QL','O3' ]
  !Varaibles that will be read when doing moist physics.
  character(len=NameSize), parameter :: GCM3DvarsReadMoist(NUM_GCM3DvarsReadMoist) = &
       [ character(NameSize) :: 'U ','V ','PT',     'QV',          'O3','QLS','QCN','CFCN' ]

  !GCM TRAJECTORY, RANK 2 VARIABLES
  !KEEP PHIS AND PS IN 1ST AND 2ND POSITION.
  character(len=NameSize), parameter :: GCM2Dvars(24) =                  &
       [ character(NameSize) ::                                          &
         'PHIS','PS','FRLAND','HS_STDV','FROCEAN',                       &  !Dynamics
         'VARFLT','USTAR','BSTAR','ZPBL','CM','CT','CQ',                 &  !Turbulence/BL
         'KCBL','TS','KHL','KHU',                                        &  !Moist Physics
         'DELTIRD','SLR','EMIS','COSZ','RGBUV','RGFUV','RGBIR','RGFIR'  ]   !Radiation
  
  type T_1DVar
     character(len=NameSize) :: Name
     real,   pointer :: X (:)=>null()
     real*8, pointer :: X8(:)=>null()
  end type T_1DVar

  type T_2DVar
     character(len=NameSize) :: Name
     real,   pointer :: X (:,:)=>null()
     real*8, pointer :: X8(:,:)=>null()
  end type T_2DVar

  type T_3DVar
     character(len=NameSize) :: Name
     real,   pointer :: X (:,:,:)=>null()
     real*8, pointer :: X8(:,:,:)=>null()
  end type T_3DVar

  real,pointer :: pert_ak(:)
  real,pointer :: pert_bk(:)
  real,pointer :: pert_shapiro_coeff(:)

  interface apert_dot_product
        module procedure dot_product1_
        module procedure dot_product2_
        module procedure dot_product3_
        module procedure dot_product4_
  end interface

  interface apert_state2state
        module procedure state2state_
  end interface
  interface apert_scal
        module procedure scal_
  end interface
  interface DefVertGrid
        module procedure DefVertGrid_
  end interface
  interface IceFraction
        module procedure IceFraction_r4
        module procedure IceFraction_r8
  end interface

  interface phase
  	module procedure phase_; end interface
  interface phase_opAdjoint
  	module procedure adjoint_; end interface
  interface phase_addImport
  	module procedure import_; end interface
  interface phase_getExport
  	module procedure export_; end interface

  interface GetShapiroCoeff
        module procedure GetShapiroCoeff_
  end interface

contains
  
!!!===================================================================


  subroutine DP2PK( DPT, PKET, PPET, BX, CX, DPP, PKP ,TRANSPOSE, PLP)

    real,    intent(IN   )  :: DPT(:,:,:), PPET(:,:,0:), PKET(:,:,0:)
    real,    intent(IN   )  :: BX(:,:,:), CX(:,:,:)
    real,    intent(INOUT)  :: PKP(:,:,:), DPP(:,:,:)
    logical, intent(IN   )  :: TRANSPOSE
    real, optional, intent(OUT) :: PLP(:,:,0:)

    real, dimension(size(DPT,1),size(DPT,2),0:size(DPT,3)) :: PEP, PLPx, PKE
    integer :: L, LM

    LM = size(DPT,3)

    if(.not.TRANSPOSE) then

       ! Compute linearized edge p (PE)
       !-------------------------------

       PEP(:,:,0) = 0.0

       do L=0,LM-1
          PEP(:,:,L+1) = PEP(:,:,L) + DPP(:,:,L+1)
       end do

       ! Compute linearized ln(p) (PL)
       !------------------------------

       PLPX = PEP/PPET

       ! Compute linearized p**kappa (PK)
       !---------------------------------

       PKE = MAPL_KAPPA*(PKET/PPET)*PEP   ! perturbation edge pressure to kappa

       PKP = BX*(PKE (:,:,1:LM) - PKE (:,:,0:LM-1)) &
           - CX*(PLPX(:,:,1:LM) - PLPX(:,:,0:LM-1))
    else

       ! Compute influence of PK on PKE and PLP
       !---------------------------------------

       PKE             = 0.0
       PKE(:,:,1:LM  ) = PKE(:,:,1:LM  ) + PKP*BX
       PKE(:,:,0:LM-1) = PKE(:,:,0:LM-1) - PKP*BX

       PLPX             = 0.0
       PLPX(:,:,1:LM  ) = PLPX(:,:,1:LM  ) - PKP*CX
       PLPX(:,:,0:LM-1) = PLPX(:,:,0:LM-1) + PKP*CX

       ! Include influence on PE from PK AND PL
       !---------------------------------------

       PEP = MAPL_KAPPA*(PKET/PPET)*PKE + PLPX/PPET

       ! Finally, compute influence on DP
       !---------------------------------

       do L=LM-1,0,-1
          PEP(:,:,L  ) = PEP(:,:,L  ) + PEP(:,:,L+1)
          DPP(:,:,L+1) = DPP(:,:,L+1) + PEP(:,:,L+1)
       enddo

    end if

    if(present(PLP)) then
       PLP = PLPX
    endif

    return
  end subroutine DP2PK


  subroutine GetPressVarsTraj(DP, PE, PKE, PKT, BX, CX, PLT)
    real,  intent(IN ) :: DP(:,:,:)
    real,  intent(OUT) :: PE(:,:,0:),  PKE(:,:,0:), PKT(:,:,:)
    real,  intent(OUT) :: BX(:,:,:),  CX (:,:,:)
    real, optional, intent(OUT) :: PLT(:,:,0:)
    integer :: L, LM

    LM = size(DP,3)

    PE (:,:,0) = PTOP
    do L=0,LM-1
       PE(:,:,L+1) = PE(:,:,L) + DP(:,:,L+1)
    end do

    PKE = log(PE)
    BX  = PKE(:,:,1:LM) - PKE(:,:,0:LM-1)

    if(present(PLT)) then
       PLT = PKE
    endif

    PKE = PE**MAPL_KAPPA
    CX  = PKE(:,:,1:LM) - PKE(:,:,0:LM-1)

    BX  = 1.0/(MAPL_KAPPA*BX)
    PKT = BX*CX
    CX  = BX*PKT*MAPL_KAPPA

    return
  end subroutine GetPressVarsTraj

  subroutine GetPressVarsTrajR8(DP, PE, PKE, PKT, BX, CX, PLT)
    real(8),  intent(IN ) :: DP(:,:,:)
    real(8),  intent(OUT) :: PE(:,:,0:),  PKE(:,:,0:), PKT(:,:,:)
    real(8),  intent(OUT) :: BX(:,:,:),  CX (:,:,:)
    real(8), optional, intent(OUT) :: PLT(:,:,0:)
    integer :: L, LM
    real(8) :: MAPL8_KAPPA

    MAPL8_KAPPA = dble(MAPL_KAPPA)

    LM = size(DP,3)

    PE (:,:,0) = PTOP
    do L=0,LM-1
       PE(:,:,L+1) = PE(:,:,L) + DP(:,:,L+1)
    end do

    PKE = log(PE)
    BX  = PKE(:,:,1:LM) - PKE(:,:,0:LM-1)

    if(present(PLT)) then
       PLT = PKE
    endif

    PKE = PE**MAPL8_KAPPA
    CX  = PKE(:,:,1:LM) - PKE(:,:,0:LM-1)

    BX  = 1.0/(MAPL8_KAPPA*BX)
    PKT = BX*CX
    CX  = BX*PKT*MAPL8_KAPPA

    return
  end subroutine GetPressVarsTrajR8

  subroutine PERT_RefVarsFromBundle(Bundle,Vars,RC)
    type (ESMF_FieldBundle),  intent(INout) :: BUNDLE
    type (T_3DVar),           intent(INOUT) :: Vars(:)
    integer,                  intent(  OUT) :: RC

    integer :: STATUS, i, NQ
    character(len=ESMF_MAXSTR) :: IAm="PERT_RefVarsFromBundle"
    character(len=ESMF_MAXSTR), allocatable :: Names(:)

    call ESMF_FieldBundleGet (Bundle, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

    allocate(Names(NQ))

    call ESMF_FieldBundleGet ( Bundle, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)

    do i=1,size(Vars)
       if(any(Names==Vars(i)%name)) then
          call MAPL_FieldBundleGetPointer(BUNDLE, Vars(i)%name, Vars(i)%X, rc=STATUS)
          VERIFY_(STATUS)
       endif
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_RefVarsFromBundle

  subroutine PERT_RefVarsFromBundle_2D(Bundle,Vars,RC)
    type (ESMF_FieldBundle),  intent(INout) :: BUNDLE
    type (T_2DVar),           intent(INOUT) :: Vars(:)
    integer,                  intent(  OUT) :: RC

    integer :: STATUS, i, NQ
    character(len=ESMF_MAXSTR) :: IAm="PERT_RefVarsFromBundle_2D"
    character(len=ESMF_MAXSTR), allocatable :: Names(:)

    call ESMF_FieldBundleGet (Bundle, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)
    
    allocate(Names(NQ))

    call ESMF_FieldBundleGet ( Bundle, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)

    do i=1,size(Vars)
       if(any(Names==Vars(i)%name)) then
          call MAPL_FieldBundleGetPointer(BUNDLE, Vars(i)%name, Vars(i)%X, rc=STATUS)
          VERIFY_(STATUS)
       endif
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_RefVarsFromBundle_2D

  subroutine PERT_RefVarsFromState(State,Vars,RC)
    type (ESMF_State),  intent(INout) :: State
    type (T_3Dvar),     intent(INOUT) :: Vars(:)
    integer,            intent(  OUT) :: RC

    integer :: Status, i
    character(len=ESMF_MAXSTR) :: IAm="PERT_GetPointers"

    do i=1,size(Vars)
       call MAPL_GetPointer(STATE, Vars(i)%X, Vars(i)%name, rc=STATUS)
       VERIFY_(STATUS)
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_RefVarsFromState

  integer function PERT_GetIndex(Vars,name) Result(N)
    type (T_3Dvar),           intent(IN) :: Vars(:)
    character(len=*),         intent(IN) :: Name

    do n=1,size(Vars)
       if(trim(Vars(n)%name)==trim(Name)) exit
    end do

  end function PERT_GetIndex

  integer function PERT_GetIndex2D(Vars,name) Result(N)
    type (T_2Dvar),           intent(IN) :: Vars(:)
    character(len=*),         intent(IN) :: Name

    do n=1,size(Vars)
       if(trim(Vars(n)%name)==trim(Name)) exit
    end do

  end function PERT_GetIndex2D

  subroutine PERT_CopyVarsFromBundle(Bundle,Vars,RC)
    type (ESMF_FieldBundle),  intent(INout) :: BUNDLE
    type (T_3Dvar),           intent(INOUT) :: Vars(:)
    integer,                  intent(  OUT) :: RC

    integer           :: i, STATUS
    real, pointer     :: X(:,:,:)

    character(len=ESMF_MAXSTR) :: IAm="PERT_CopyVarsFromBundle"

    do i=1,size(Vars)
       ASSERT_(associated(Vars(i)%X))

       call MAPL_FieldBundleGetPointer(BUNDLE, Vars(i)%name, X, rc=STATUS)
       VERIFY_(STATUS)

       Vars(i)%X  = X
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_CopyVarsFromBundle
    
  subroutine PERT_CopyVarsFromBundle_2D(Bundle,Vars,RC)
    type (ESMF_FieldBundle),  intent(INout) :: BUNDLE
    type (T_2Dvar),           intent(INOUT) :: Vars(:)
    integer,                  intent(  OUT) :: RC

    integer           :: i, STATUS
    real, pointer     :: X(:,:)

    character(len=ESMF_MAXSTR) :: IAm="PERT_CopyVarsFromBundle_2D"

    do i=1,size(Vars)
       ASSERT_(associated(Vars(i)%X))

       call MAPL_FieldBundleGetPointer(BUNDLE, Vars(i)%name, X, rc=STATUS)
       VERIFY_(STATUS)

       Vars(i)%X  = X
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_CopyVarsFromBundle_2D

  subroutine PERT_Null(Vars,RC)
    type (T_3Dvar),           intent(INOUT) :: Vars(:)
    integer,                  intent(  OUT) :: RC

    integer           :: i, STATUS
    real, pointer     :: X(:,:,:)

    character(len=ESMF_MAXSTR) :: IAm="PERT_Null"

    do i=1,size(Vars)
       ASSERT_(associated(Vars(i)%X))
       nullify(Vars(i)%X)
    end do
    RC=0

    RETURN_(ESMF_SUCCESS)
  end subroutine PERT_Null

  function vec2list (vec) result(var)
  implicit none
  character(len=ESMF_MAXSTR)   :: var
  character(len=*),intent(in)  :: vec(:)
!
  character(len=ESMF_MAXSTR),allocatable :: work(:)
  character(len=ESMF_MAXSTR) list,desc
  integer is,ie,i,i0, nt,istatus
     istatus=-1
     nt=size(vec)
     if(nt>0) then
        allocate(work(size(vec)))
        work=vec
        list=trim(work(1))
        do i=2,nt
           i0=len_trim(list)
           is=i0+1
           ie=is+len_trim(work(i))+1
           list(is:ie)=','//work(i)
        enddo
        if(nt>1.and.list(1:1)==',') list=list(2:ie)
        var = trim(list)
        if(var/='') istatus=0
        deallocate(work)
     endif
  end function vec2list

   function dot_product1_ (XBundle,YBundle,rc) result(dprd)
   implicit none
   type(ESMF_FieldBundle) :: XBundle
   type(ESMF_FieldBundle) :: YBundle
   integer, optional,  intent(OUT)   :: rc
   real*8 dprd
   character(len=ESMF_MAXSTR) :: IAm="dot_product1_"
   real, pointer    :: xptr(:,:,:), yptr(:,:,:)
   type(ESMF_VM) :: vm
   real*8 this_dot(1)
   integer iv,i,j,k,imax,jmax,kmax
   integer ier,status
   this_dot=0.d0
   do iv = 1, NUM_GSIvars
      status=0
      call MAPL_FieldBundleGetPointer(XBundle, trim(GSIvars(iv)), xptr, rc=ier )
      status=ier+status
      call MAPL_FieldBundleGetPointer(YBundle, trim(GSIvars(iv)), yptr, rc=ier )
      status=ier+status
      VERIFY_(status)
      imax=size(xptr,1)
      jmax=size(xptr,2)
      kmax=size(xptr,3)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         this_dot(1) = this_dot(1) + xptr(i,j,k)*yptr(i,j,k)
      enddo
      enddo
      enddo
   enddo
   call ESMF_VmGetCurrent(vm)
   call MAPL_CommsAllReduceSum(vm, sendbuf=sum(this_dot), &
                                   recvbuf=dprd, cnt= 1, rc=status)
   end function dot_product1_

   function dot_product2_ (XState,YState,rc) result(dprd)
   implicit none
   type(ESMF_State) :: XState
   type(ESMF_State) :: YState
   integer, optional,  intent(OUT)   :: rc
   real*8 dprd
   character(len=ESMF_MAXSTR) :: IAm="dot_product2_"
   type(ESMF_VM) :: vm
   real,   pointer    :: xptr(:,:,:), yptr(:,:,:)
   integer iv,i,j,k,imax,jmax,kmax
   integer ier,status
   real*8 this_dot(1)
   this_dot=0.d0
   do iv = 1, NUM_GSIvars
      status=0
      call MAPL_GetPointer(XState, xptr, trim(GSIvars(iv)), rc=ier)
      status=ier+status
      call MAPL_GetPointer(YState, yptr, trim(GSIvars(iv)), rc=ier)
      status=ier+status
      VERIFY_(status)
      imax=size(xptr,1)
      jmax=size(xptr,2)
      kmax=size(xptr,3)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         this_dot(1) = this_dot(1) + xptr(i,j,k)*yptr(i,j,k)
      enddo
      enddo
      enddo
   enddo
   call ESMF_VmGetCurrent(vm)
   call MAPL_CommsAllReduceSum(vm, sendbuf=sum(this_dot), &
                                   recvbuf=dprd, cnt= 1, rc=status)
   end function dot_product2_

   function dot_product3_ (XState,YBundle,rc) result(dprd)
   implicit none
   type(ESMF_State)       :: XState
   type(ESMF_FieldBundle) :: YBundle
   integer, optional,  intent(OUT)   :: rc
   real*8 dprd
   character(len=ESMF_MAXSTR) :: IAm="dot_product3_"
   type(ESMF_VM) :: vm
   real,   pointer    :: xptr(:,:,:), yptr(:,:,:)
   integer iv,i,j,k,imax,jmax,kmax
   integer ier,status
   real*8 this_dot(1)
   this_dot=0.d0
   do iv = 1, NUM_GSIvars
      status=0
      call MAPL_GetPointer(XState, xptr, trim(GSIvars(iv)), rc=ier)
      status=ier+status
      call MAPL_FieldBundleGetPointer(YBundle, trim(GSIvars(iv)), yptr, rc=ier )
      status=ier+status
      VERIFY_(status)
      imax=size(xptr,1)
      jmax=size(xptr,2)
      kmax=size(xptr,3)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         this_dot(1) = this_dot(1) + xptr(i,j,k)*yptr(i,j,k)
      enddo
      enddo
      enddo
   enddo
   call ESMF_VmGetCurrent(vm)
   call MAPL_CommsAllReduceSum(vm, sendbuf=sum(this_dot), &
                                   recvbuf=dprd, cnt= 1, rc=status)
   end function dot_product3_

   function dot_product4_ (xvars,yvars,rc) result(dprd)
   implicit none
   type(T_3DVar) :: xvars(:)
   type(T_3DVar) :: yvars(:)
   integer, optional,  intent(OUT)   :: rc
   real*8 dprd
   character(len=ESMF_MAXSTR) :: IAm="dot_product4_"
   type(ESMF_VM) :: vm
   integer iv,i,j,k,imax,jmax,kmax
   integer ier,status
   real*8 this_dot(1)
   this_dot=0.d0
   do iv = 1, size(xvars)
      imax=size(xvars(iv)%X,1)
      jmax=size(xvars(iv)%X,2)
      kmax=size(xvars(iv)%X,3)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         this_dot(1) = this_dot(1) + xvars(iv)%X(i,j,k)*yvars(iv)%X(i,j,k)
      enddo
      enddo
      enddo
   enddo
   rc=0
   call ESMF_VmGetCurrent(vm)
   call MAPL_CommsAllReduceSum(vm, sendbuf=sum(this_dot), &
                                   recvbuf=dprd, cnt= 1, rc=status)
   end function dot_product4_

   subroutine state2state_ (XState,YState,rc)
   implicit none
   type(ESMF_State) :: XState
   type(ESMF_State) :: YState
   integer, optional,  intent(OUT)   :: rc
   character(len=ESMF_MAXSTR) :: IAm="state2state_"
!  real*8, pointer    :: xptr(:,:,:), yptr(:,:,:)   !_RT this is what they should be
   real,   pointer    :: xptr(:,:,:), yptr(:,:,:)
   integer iv
   integer ier,status
   do iv = 1, NUM_GSIvars
      status=0
      call MAPL_GetPointer(XState, xptr, trim(GSIvars(iv)), rc=ier)
      status=ier+status
      call MAPL_GetPointer(YState, yptr, trim(GSIvars(iv)), rc=ier)
      status=ier+status
      VERIFY_(status)
      yptr = xptr
   enddo
   end subroutine state2state_

   subroutine scal_ (a,XState,rc)
   implicit none
   real,intent(in)  :: a
   type(ESMF_State) :: XState
   integer, optional,  intent(OUT)   :: rc
   character(len=ESMF_MAXSTR) :: IAm="state2state_"
!  real*8, pointer    :: xptr(:,:,:), yptr(:,:,:)   !_RT this is what they should be
   real,   pointer    :: xptr(:,:,:)
   integer iv
   integer ier,status
   do iv = 1, NUM_GSIvars
      call MAPL_GetPointer(XState, xptr, trim(GSIvars(iv)), rc=status)
      VERIFY_(status)
      xptr = a
   enddo
   end subroutine scal_

!_RT: TAF routines here a for debug purpose only and will be removed when all sorted all
subroutine getpressvarsGCM2GSI_tl( pt, pt_tl, qv, qv_tl, dp, dp_tl, pe, pe_tl, pke, pke_tl, pkt, pkt_tl, tv_tl, plt, plt_tl )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 2.1.10  **
!******************************************************************
!******************************************************************
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

real, parameter :: EPS = MAPL_VIREPS

!==============================================
! declare arguments
!==============================================
real, intent(in) :: dp(:,:,:)
real, intent(in) :: dp_tl(:,:,:)
real, intent(out) :: pe(:,:,0:)
real, intent(out) :: pe_tl(:,:,0:)
real, intent(out) :: pke(:,:,0:)
real, intent(out) :: pke_tl(:,:,0:)
real, intent(out) :: pkt(:,:,:)
real, intent(out) :: pkt_tl(:,:,:)
real,optional, intent(out) :: plt(:,:,0:)
real,optional, intent(out) :: plt_tl(:,:,0:)
real, intent(in) :: pt(:,:,:)
real, intent(in) :: pt_tl(:,:,:)
real, intent(in) :: qv(:,:,:)
real, intent(in) :: qv_tl(:,:,:)
!real, intent(out) :: tv(:,:,:)
real, intent(out) :: tv_tl(:,:,:)

!==============================================
! declare local variables
!==============================================
real, allocatable :: bx(:,:,:)
real, allocatable :: bx_tl(:,:,:)
real, allocatable :: cx(:,:,:)
real, allocatable :: cx_tl(:,:,:)
integer :: l
integer :: lm

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
lm = size(dp,3)
allocate( bx_tl(size(dp,1),size(dp,2),lm) )
allocate( bx(size(dp,1),size(dp,2),lm) )
allocate( cx_tl(size(dp,1),size(dp,2),lm) )
allocate( cx(size(dp,1),size(dp,2),lm) )
pe_tl(:,:,0) = 0.
pe(:,:,0) = PTOP
do l = 0, lm-1
  pe_tl(:,:,l+1) = dp_tl(:,:,l+1)+pe_tl(:,:,l)
  pe(:,:,l+1) = pe(:,:,l)+dp(:,:,l+1)
end do
pke_tl = pe_tl*(1./pe)
pke = log(pe)
bx_tl = pke_tl(:,:,1:lm)-pke_tl(:,:,0:lm-1)
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
if (present(plt)) then
  plt_tl = pke_tl
  plt = pke
endif
pke_tl = pe_tl*mapl_kappa*pe**(mapl_kappa-1)
pke = pe**mapl_kappa
cx_tl = pke_tl(:,:,1:lm)-pke_tl(:,:,0:lm-1)
cx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx_tl = -(bx_tl*(1.*mapl_kappa/(mapl_kappa*bx*mapl_kappa*bx)))
bx = 1./(mapl_kappa*bx)
pkt_tl = bx_tl*cx+cx_tl*bx
pkt = bx*cx
deallocate( cx_tl )
deallocate( cx )
deallocate( bx_tl )
deallocate( bx )
!tv_tl = pkt_tl*(1+eps*qv)*pt+pt_tl*(1+eps*qv)*pkt+qv_tl*eps*pt*pkt
!tv = (1+eps*qv)*pt*pkt
TV_TL = PT*PKT*EPS*QV_TL + (1.0 + EPS*QV)*(PKT_TL*PT + PKT*PT_TL)

end subroutine getpressvarsGCM2GSI_tl


subroutine getpressvarsGCM2GSI_ad( pt, pt_ad, qv, qv_ad, dp, dp_ad, pe, pe_ad, pke, pke_ad, pkt, pkt_ad, tv, tv_ad, plt, plt_ad )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 2.1.10  **
!******************************************************************
!******************************************************************
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

real, parameter :: EPS = MAPL_VIREPS

!==============================================
! declare arguments
!==============================================
real, intent(in) :: dp(:,:,:)
real, intent(inout) :: dp_ad(:,:,:)
real, intent(out) :: pe(:,:,0:)
real, intent(inout) :: pe_ad(:,:,0:)
real, intent(out) :: pke(:,:,0:)
real, intent(inout) :: pke_ad(:,:,0:)
real, intent(out) :: pkt(:,:,:)
real, intent(inout) :: pkt_ad(:,:,:)
real,optional, intent(out) :: plt(:,:,0:)
real,optional, intent(inout) :: plt_ad(:,:,0:)
real, intent(in) :: pt(:,:,:)
real, intent(inout) :: pt_ad(:,:,:)
real, intent(in) :: qv(:,:,:)
real, intent(inout) :: qv_ad(:,:,:)
real, intent(out) :: tv(:,:,:)
real, intent(inout) :: tv_ad(:,:,:)

!==============================================
! declare local variables
!==============================================
real, allocatable :: bx(:,:,:)
real, allocatable :: bx_ad(:,:,:)
real, allocatable :: cx(:,:,:)
real, allocatable :: cx_ad(:,:,:)
integer :: l
integer :: lm
real, allocatable :: pei_ad(:,:)
real :: pej(lbound(pe,1):ubound(pe,1),lbound(pe,2):ubound(pe,2),lbound(pe,3):ubound(pe,3))
real :: pkth(lbound(pkt,1):ubound(pkt,1),lbound(pkt,2):ubound(pkt,2),lbound(pkt,3):ubound(pkt,3))

!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
!----------------------------------------------
! FUNCTION AND TAPE COMPUTATIONS
!----------------------------------------------
lm = size(dp,3)
allocate( bx(size(dp,1),size(dp,2),lm) )
allocate( cx(size(dp,1),size(dp,2),lm) )
allocate( pei_ad(size(pe,1),size(pe,2)) )
pe(:,:,0) = PTOP
do l = 0, lm-1
  pe(:,:,l+1) = pe(:,:,l)+dp(:,:,l+1)
end do
pke = log(pe)
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
if (present(plt)) then
  plt = pke
endif
pke = pe**mapl_kappa
cx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx = 1./(mapl_kappa*bx)
pkt = bx*cx
tv = (1+eps*qv)*pt*pkt

!----------------------------------------------
! SAVE DEPENDEND VARIABLES
!----------------------------------------------
pkth(:,:,:) = pkt(:,:,:)
pej(:,:,:) = pe(:,:,:)

!----------------------------------------------
! ADJOINT COMPUTATIONS
!----------------------------------------------
pkt_ad = pkt_ad+tv_ad*(1+eps*qv)*pt
pt_ad = pt_ad+tv_ad*(1+eps*qv)*pkt
qv_ad = qv_ad+tv_ad*eps*pt*pkt
tv_ad = 0.
allocate( bx_ad(size(dp,1),size(dp,2),lm) )
bx_ad = 0.
allocate( cx_ad(size(dp,1),size(dp,2),lm) )
cx_ad = 0.
bx_ad = bx_ad+pkt_ad*cx
cx_ad = cx_ad+pkt_ad*bx
pkt_ad = 0.
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx_ad = -(bx_ad*(1.*mapl_kappa/(mapl_kappa*bx*mapl_kappa*bx)))
pke_ad(:,:,1:lm) = pke_ad(:,:,1:lm)+cx_ad
pke_ad(:,:,0:lm-1) = pke_ad(:,:,0:lm-1)-cx_ad
cx_ad = 0.
pe_ad = pe_ad+pke_ad*mapl_kappa*pe**(mapl_kappa-1)
pke_ad = 0.
if (present(plt)) then
  pke_ad = pke_ad+plt_ad
  plt_ad = 0.
endif
pke_ad(:,:,1:lm) = pke_ad(:,:,1:lm)+bx_ad
pke_ad(:,:,0:lm-1) = pke_ad(:,:,0:lm-1)-bx_ad
bx_ad = 0.
pe_ad = pe_ad+pke_ad*(1./pe)
pke_ad = 0.
do l = lm-1, 0, -1
  pei_ad = pe_ad(:,:,l+1)
  pe_ad(:,:,l+1) = 0.
  dp_ad(:,:,l+1) = dp_ad(:,:,l+1)+pei_ad
  pe_ad(:,:,l) = pe_ad(:,:,l)+pei_ad
end do
deallocate( cx_ad )
deallocate( bx_ad )
!----------------------------------------------
! GET DEPENDEND VARIABLES
!----------------------------------------------
pe(:,:,:) = pej(:,:,:)
pkt(:,:,:) = pkth(:,:,:)


!----------------------------------------------
! DEALLOCATE STATEMENTS
!----------------------------------------------
if (allocated(cx)) then
  deallocate( cx )
endif
if (allocated(bx)) then
  deallocate( bx )
endif
if (allocated(pei_ad)) then
  deallocate( pei_ad )
endif
end subroutine getpressvarsGCM2GSI_ad


subroutine getpressvarsGSI2GCM_tl( pt, pt_tl, qv, qv_tl, dp, dp_tl, pe, pe_tl, pke, pke_tl, pkt, pkt_tl, tv, tv_tl, plt, plt_tl )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 2.1.10  **
!******************************************************************
!******************************************************************
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

real, parameter :: EPS = MAPL_VIREPS

!==============================================
! declare arguments
!==============================================
real, intent(in) :: dp(:,:,:)
real, intent(in) :: dp_tl(:,:,:)
real, intent(out) :: pe(:,:,0:)
real, intent(out) :: pe_tl(:,:,0:)
real, intent(out) :: pke(:,:,0:)
real, intent(out) :: pke_tl(:,:,0:)
real, intent(out) :: pkt(:,:,:)
real, intent(out) :: pkt_tl(:,:,:)
real,optional, intent(out) :: plt(:,:,0:)
real,optional, intent(out) :: plt_tl(:,:,0:)
!real, intent(out) :: pt(:,:,:)
real, intent(in) :: pt(:,:,:)
real, intent(out) :: pt_tl(:,:,:)
real, intent(in) :: qv(:,:,:)
real, intent(in) :: qv_tl(:,:,:)
real, intent(in) :: tv(:,:,:)
real, intent(in) :: tv_tl(:,:,:)

!==============================================
! declare local variables
!==============================================
real, allocatable :: bx(:,:,:)
real, allocatable :: bx_tl(:,:,:)
real, allocatable :: cx(:,:,:)
real, allocatable :: cx_tl(:,:,:)
real, allocatable, dimension(:,:,:)::AY,BY,CY
integer :: l
integer :: lm

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
lm = size(dp,3)
allocate( bx_tl(size(dp,1),size(dp,2),lm) )
allocate( bx(size(dp,1),size(dp,2),lm) )
allocate( cx_tl(size(dp,1),size(dp,2),lm) )
allocate( cx(size(dp,1),size(dp,2),lm) )
pe_tl(:,:,0) = 0.
pe(:,:,0) = PTOP
do l = 0, lm-1
  pe_tl(:,:,l+1) = dp_tl(:,:,l+1)+pe_tl(:,:,l)
  pe(:,:,l+1) = pe(:,:,l)+dp(:,:,l+1)
end do
pke_tl = pe_tl*(1./pe)
pke = log(pe)
bx_tl = pke_tl(:,:,1:lm)-pke_tl(:,:,0:lm-1)
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
if (present(plt)) then
  plt_tl = pke_tl
  plt = pke
endif
pke_tl = pe_tl*mapl_kappa*pe**(mapl_kappa-1)
pke = pe**mapl_kappa
cx_tl = pke_tl(:,:,1:lm)-pke_tl(:,:,0:lm-1)
cx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx_tl = -(bx_tl*(1.*mapl_kappa/(mapl_kappa*bx*mapl_kappa*bx)))
bx = 1./(mapl_kappa*bx)
pkt_tl = bx_tl*cx+cx_tl*bx
pkt = bx*cx
deallocate( cx_tl )
deallocate( cx )
deallocate( bx_tl )
deallocate( bx )
!pt_tl = pkt_tl*(tv/(1+eps*qv))-qv_tl*tv*eps/((1+eps*qv)*(1+eps*qv))*pkt+tv_tl/(1+eps*qv)*pkt
!pt = tv/(1+eps*qv)*pkt
    allocate(AY (size(DP,1),size(DP,2),1:LM))
    allocate(BY (size(DP,1),size(DP,2),1:LM))
    allocate(CY (size(DP,1),size(DP,2),1:LM))
AY  = 1.0/(1.0 + EPS*QV)
CY  = AY*PT/(PKT*PKT)
AY  = AY/PKT
BY  = AY*PT*EPS
PT_TL = AY*TV_TL - BY*QV_TL - CY*PKT_TL
    deallocate(AY,BY,CY)

end subroutine getpressvarsGSI2GCM_tl


subroutine getpressvarsGSI2GCM_ad( pt, pt_ad, qv, qv_ad, dp, dp_ad, pe, pe_ad, pke, pke_ad, pkt, pkt_ad, tv, tv_ad, plt, plt_ad )
!******************************************************************
!******************************************************************
!** This routine was generated by Automatic differentiation.     **
!** FastOpt: Transformation of Algorithm in Fortran, TAF 2.1.10  **
!******************************************************************
!******************************************************************
!==============================================
! all entries are defined explicitly
!==============================================
implicit none

real, parameter :: EPS = MAPL_VIREPS

!==============================================
! declare arguments
!==============================================
real, intent(in) :: dp(:,:,:)
real, intent(inout) :: dp_ad(:,:,:)
real, intent(out) :: pe(:,:,0:)
real, intent(inout) :: pe_ad(:,:,0:)
real, intent(out) :: pke(:,:,0:)
real, intent(inout) :: pke_ad(:,:,0:)
real, intent(out) :: pkt(:,:,:)
real, intent(inout) :: pkt_ad(:,:,:)
real,optional, intent(out) :: plt(:,:,0:)
real,optional, intent(inout) :: plt_ad(:,:,0:)
real, intent(in) :: pt(:,:,:)
real, intent(inout) :: pt_ad(:,:,:)
real, intent(in) :: qv(:,:,:)
real, intent(inout) :: qv_ad(:,:,:)
!real, intent(in) :: tv(:,:,:)
real, intent(out) :: tv(:,:,:) ! RT
real, intent(inout) :: tv_ad(:,:,:)

!==============================================
! declare local variables
!==============================================
real, allocatable :: bx(:,:,:)
real, allocatable :: bx_ad(:,:,:)
real, allocatable :: cx(:,:,:)
real, allocatable :: cx_ad(:,:,:)
integer :: l
integer :: lm
real, allocatable :: pei_ad(:,:)
real :: pej(lbound(pe,1):ubound(pe,1),lbound(pe,2):ubound(pe,2),lbound(pe,3):ubound(pe,3))
real :: pkth(lbound(pkt,1):ubound(pkt,1),lbound(pkt,2):ubound(pkt,2),lbound(pkt,3):ubound(pkt,3))

!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
!----------------------------------------------
! FUNCTION AND TAPE COMPUTATIONS
!----------------------------------------------
lm = size(dp,3)
allocate( bx(size(dp,1),size(dp,2),lm) )
allocate( cx(size(dp,1),size(dp,2),lm) )
allocate( pei_ad(size(pe,1),size(pe,2)) )
pe(:,:,0) = PTOP
do l = 0, lm-1
  pe(:,:,l+1) = pe(:,:,l)+dp(:,:,l+1)
end do
pke = log(pe)
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
if (present(plt)) then
  plt = pke
endif
pke = pe**mapl_kappa
cx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx = 1./(mapl_kappa*bx)
pkt = bx*cx
!pt = tv/(1+eps*qv)*pkt
tv = pt*(1+eps*qv)*pkt

!----------------------------------------------
! SAVE DEPENDEND VARIABLES
!----------------------------------------------
pkth(:,:,:) = pkt(:,:,:)
pej(:,:,:) = pe(:,:,:)

!----------------------------------------------
! ADJOINT COMPUTATIONS
!----------------------------------------------
pkt_ad = pkt_ad+pt_ad*(tv/(1+eps*qv))
qv_ad = qv_ad-pt_ad*tv*eps/((1+eps*qv)*(1+eps*qv))*pkt
tv_ad = tv_ad+pt_ad/(1+eps*qv)*pkt
!? pt_ad = 0.
allocate( bx_ad(size(dp,1),size(dp,2),lm) )
bx_ad = 0.
allocate( cx_ad(size(dp,1),size(dp,2),lm) )
cx_ad = 0.
bx_ad = bx_ad+pkt_ad*cx
cx_ad = cx_ad+pkt_ad*bx
pkt_ad = 0.
bx = pke(:,:,1:lm)-pke(:,:,0:lm-1)
bx_ad = -(bx_ad*(1.*mapl_kappa/(mapl_kappa*bx*mapl_kappa*bx)))
pke_ad(:,:,1:lm) = pke_ad(:,:,1:lm)+cx_ad
pke_ad(:,:,0:lm-1) = pke_ad(:,:,0:lm-1)-cx_ad
cx_ad = 0.
pe_ad = pe_ad+pke_ad*mapl_kappa*pe**(mapl_kappa-1)
pke_ad = 0.
if (present(plt)) then
  pke_ad = pke_ad+plt_ad
  plt_ad = 0.
endif
pke_ad(:,:,1:lm) = pke_ad(:,:,1:lm)+bx_ad
pke_ad(:,:,0:lm-1) = pke_ad(:,:,0:lm-1)-bx_ad
bx_ad = 0.
pe_ad = pe_ad+pke_ad*(1./pe)
pke_ad = 0.
do l = lm-1, 0, -1
  pei_ad = pe_ad(:,:,l+1)
  pe_ad(:,:,l+1) = 0.
  dp_ad(:,:,l+1) = dp_ad(:,:,l+1)+pei_ad
  pe_ad(:,:,l) = pe_ad(:,:,l)+pei_ad
end do
deallocate( cx_ad )
deallocate( bx_ad )
!----------------------------------------------
! GET DEPENDEND VARIABLES
!----------------------------------------------
pe(:,:,:) = pej(:,:,:)
pkt(:,:,:) = pkth(:,:,:)


!----------------------------------------------
! DEALLOCATE STATEMENTS
!----------------------------------------------
if (allocated(cx)) then
  deallocate( cx )
endif
if (allocated(bx)) then
  deallocate( bx )
endif
if (allocated(pei_ad)) then
  deallocate( pei_ad )
endif

end subroutine getpressvarsGSI2GCM_ad

 
function phase_(adjoint,import,export)
  implicit none
  logical,optional,intent(in):: adjoint
  logical,optional,intent(in):: import
  logical,optional,intent(in):: export
  integer:: phase_

  logical:: adjoint_,import_,export_

  adjoint_=.false.; if(present(adjoint)) adjoint_=adjoint
  import_ =.true. ; if(present(import )) import_ =import
  export_ =.true. ; if(present(export )) export_ =export

  phase_=0
  if(    adjoint_) phase_=ior(phase_,1)
  if(.not.import_) phase_=ior(phase_,2)
!  if(.not.export_) phase_=ior(phase_,4)
  phase_=phase_+1
end function phase_

function adjoint_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: adjoint_
  adjoint_ = iand(phase-1,1)==1
end function adjoint_
function import_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: import_
  import_  = iand(phase-1,2)==0
end function import_
function export_(phase)
  implicit none
  integer,intent(in):: phase
  logical:: export_
  export_  = iand(phase-1,4)==0
end function export_

   subroutine DefVertGrid_(CF,ak,bk,nsig,verbose,rc)
   use MAPL_Mod, only: MAPL_ROOT
   use m_mpif90,only : MP_REAL8
   use m_set_eta,only : set_eta
!-------------------------------------------------------------------------
!
! !REVISION HISTORY:
!
!-------------------------------------------------------------------------
   type(ESMF_Config)      :: CF
   integer, intent(in   ) :: nsig
   real,    intent(inout) :: ak(0:nsig), bk(0:nsig)
   logical, intent(in)    :: verbose
   integer, intent(out)   :: rc

!  local variables
   character(len=*), parameter       :: IAm='DefVertGrid_'
   character(len=20)                 :: vgridlabl
   character(len=3)                  :: cnsig
   integer                           :: comm
   integer                           :: i,k,ks
   real*8                            :: ptop,pint
   real*8 , allocatable              :: ak8(:), bk8(:)
   type(ESMF_VM)                     :: vm       ! ESMF Virtual Machine

! start

   if (verbose) then
      if(MAPL_AM_I_ROOT()) print *,trim(Iam),': Get GCM vertical grid '
   endif

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    call ESMF_VMGet(vm,mpiCommunicator=comm,rc=rc)

! Create the label to be searched for in the RC file based on nsig

   allocate(ak8(nsig+1),bk8(nsig+1))
   if(MAPL_AM_I_ROOT()) then
     call set_eta ( nsig,ks,ptop,pint,ak8,bk8 )
   endif
   call mpi_bcast(ak8, nsig+1,MP_REAL8,MAPL_root,comm,rc)
   call mpi_bcast(bk8, nsig+1,MP_REAL8,MAPL_root,comm,rc)
   do i=0,nsig
      ak(i)=ak8(i+1) 
      bk(i)=bk8(i+1) 
   end do
   deallocate(ak8,bk8)

   if (verbose) then
      if(MAPL_AM_I_ROOT()) then
         print *,trim(IAm),' - lev, ak, bk - '
         do i=0,nsig
            write(*,'(1x,i3,2f16.6)') i,ak(i),bk(i)
         end do
      end if
   end if
   rc=0

   end subroutine DefVertGrid_

   subroutine IceFraction_r4(TEMP, ICEFRCT)

   IMPLICIT NONE

   !Inputs
   real(4), intent(in) :: TEMP

   !Outputs
   real(4), intent(out) :: ICEFRCT

   !Locals
   real(4), parameter :: T_ICE_ALL = 233.16, T_ICE_MAX = 273.16
   integer, parameter :: ICEFRPWR = 4


    ICEFRCT  = 0.00
    if ( TEMP <= T_ICE_ALL ) then
       ICEFRCT = 1.000
    else if ( (TEMP > T_ICE_ALL) .AND. (TEMP <= T_ICE_MAX) ) then
       ICEFRCT = 1.00 -  ( TEMP - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) 
    end if
 
    ICEFRCT = MIN(ICEFRCT,1.00)
    ICEFRCT = MAX(ICEFRCT,0.00)

    ICEFRCT = ICEFRCT**ICEFRPWR

   end subroutine IceFraction_r4

   subroutine IceFraction_r8(TEMP, ICEFRCT)

   IMPLICIT NONE

   !Inputs
   real(8), intent(in) :: TEMP

   !Outputs
   real(8), intent(out) :: ICEFRCT

   !Locals
   real(8), parameter :: T_ICE_ALL = 233.16, T_ICE_MAX = 273.16
   integer, parameter :: ICEFRPWR = 4

    ICEFRCT  = 0.00
    if ( TEMP <= T_ICE_ALL ) then
       ICEFRCT = 1.000
    else if ( (TEMP > T_ICE_ALL) .AND. (TEMP <= T_ICE_MAX) ) then
       ICEFRCT = 1.00 -  ( TEMP - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) 
    end if
 
    ICEFRCT = MIN(ICEFRCT,1.00)
    ICEFRCT = MAX(ICEFRCT,0.00)

    ICEFRCT = ICEFRCT**ICEFRPWR

   end subroutine IceFraction_r8

   subroutine GetShapiroCoeff_ (CF,shapiro_coeff_,nsig,verbose,rc)
!-------------------------------------------------------------------------
!
! !REVISION HISTORY:
!    03Dec2015 Todling - label for Shapiro Coeffs now holds number of levels
!
!-------------------------------------------------------------------------
   type(ESMF_Config)      :: CF
   real,    intent(inout) :: shapiro_coeff_(1:nsig)
   integer, intent(in   ) :: nsig
   logical, intent(in)    :: verbose
   integer, intent(out)   :: rc

!  local variables
   character(len=*), parameter       :: IAm='GetShapiroCoeff_'
   character(len=40)                 :: shaplabel
   integer                           :: i,k,status
   real*8                            :: shapiro_coeffr4

! start

   if (verbose) then
      if(MAPL_AM_I_ROOT()) print *,trim(Iam),': Get PERT Shapiro filter coefficients '
   endif

! Create the label to be searched for in the RC file based on nsig

   write(shaplabel,'(a,i3.3,a)') 'AGCM_SHAPIRO_COEFF',nsig,':'
   CALL ESMF_ConfigFindLabel(CF, label = trim(shaplabel), rc=status )
   VERIFY_(STATUS)

   DO i = 1, nsig
      CALL ESMF_ConfigNextLine    (CF, rc=status)
      VERIFY_(STATUS)
      CALL ESMF_ConfigGetAttribute(CF, k, rc=status)
      VERIFY_(STATUS)
      CALL ESMF_ConfigGetAttribute(CF, shapiro_coeffr4, rc=status)
      VERIFY_(STATUS)
      shapiro_coeff_(i)=shapiro_coeffr4
   END DO

   if (verbose) then
      if (MAPL_AM_I_ROOT()) then
         print *,trim(IAm),' - lev, filter coeff - '
         do i=1, nsig
            write(*,'(1x,i3,f16.6)') i,shapiro_coeff_(i)
         end do
      end if
   end if

   end subroutine GetShapiroCoeff_

end module GEOS_PertSharedMod

