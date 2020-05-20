
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE:  Aerosol_Callbacks --- Call back methods for aerosol optics.
!                             
!
! !INTERFACE:
!
module  Chem_AeroGeneric

! !USES:
   use ESMF
   use MAPL
   USE Chem_MieMod2G

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:
   public  add_aero
   public  run_aerosol_optics1  ! Method that runs the aerosol_optics method for each child of GOCART2G.
   public append_to_bundle
   public determine_data_driven
!
! !DESCRIPTION:
!
!  These modules compute aerosol optical properties for GOCART2G.
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco - Initial code.
!  11Jul2005 da Silva   Standardization.
!  30Dec2019 Sherman, da Silva, Darmenov, Clune - 2nd Gen. Made ESMF compliant.
!                                                 No longer relies on Chem_Reg
!
!EOP
!-------------------------------------------------------------------------
contains


!====================================================================================
  subroutine add_aero (state, label, label2, grid, typekind, ptr, rc)

!   Description: Adds fields to aero state for aerosol optics calcualtions. 

    implicit none

    type (ESMF_State),                          intent(inout)     :: state
    character (len=*),                          intent(in   )     :: label
    character (len=*),                          intent(in   )     :: label2
    type (ESMF_Grid),                           intent(inout)     :: grid
    integer,                                    intent(in   )     :: typekind
    real, pointer, dimension(:,:,:), optional,  intent(in   )     :: ptr
    integer,                                    intent(  out)     :: rc

    ! locals
    type (ESMF_Field)                                             :: field
    character (len=ESMF_MAXSTR)                                   :: field_name

    __Iam__('add_aero')

!----------------------------------------------------------------------------------
!   Begin...

    call ESMF_AttributeSet (state, name=trim(label), value=trim(label2),  __RC__)

    call ESMF_AttributeGet (state, name=trim(label), value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=typekind, hw=0, __RC__)
        call MAPL_StateAdd (state, field, __RC__)
    end if

!   if (field_name /= '') then
!       field = ptr
!       call MAPL_StateAdd (state, field, __RC__)
!   end if

    RETURN_(ESMF_SUCCESS)

  end subroutine add_aero

!=====================================================================================

  subroutine determine_data_driven(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,               intent(  out)   :: RC          ! Error code:
    character (len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                         intent(  out)   :: data_driven

    !Local
    integer                                          :: i

!   Description: Determines whether gridded component is data driven or not.

     __Iam__('determine_data_driven')

!   Begin... 

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine determine_data_driven

!=====================================================================================

  subroutine append_to_bundle(varName, providerState, prefix, bundle, rc)

    implicit none

!   !ARGUMENTS:
    character (len=*),           intent(in   )   :: varName, prefix
    type (ESMF_State),           intent(in   )   :: providerState
    type (ESMF_FieldBundle),     intent(inout)   :: bundle
    integer,                     intent(  out)   :: rc  ! return code

!   !Local
    type (ESMF_Field)                              :: field

!   Description: Adds deposition variables to deposition bundle

     __Iam__('append_to_bundle')

!   Dry deposition
!   ---------------
    call ESMF_StateGet (providerState, trim(prefix)//trim(varName), field, __RC__)
    call MAPL_AllocateCoupling (field, __RC__)
    call MAPL_FieldBundleAdd (bundle, field, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine append_to_bundle

!===================================================================================
  subroutine run_aerosol_optics1 (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,                intent(out)              :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: PLE
    real, dimension(:,:,:), pointer                  :: RH
    real, dimension(:,:,:), pointer                  :: var

    character (len=ESMF_MAXSTR)                      :: fld_name

    real, dimension(:,:,:),pointer                   :: ext_, ssa_, asy_      ! (lon:,lat:,lev:,band:)
    real, dimension(:,:,:), allocatable              :: ext,  ssa,  asy       ! (lon:,lat:,lev:,band:)

    integer                                          :: i, n, b, j
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band
    integer, parameter                               :: n_bands = 1

    character (len=ESMF_MAXSTR), allocatable         :: itemList(:), AEROlist(:)
    type (ESMF_State)                                :: child_state
    real, pointer,     dimension(:,:,:)              :: AS_PTR_3D

    type (ESMF_StateItem_Flag), allocatable          :: itemTypes(:)

!   integer :: status
!   character (len=ESMF_MAXSTR)                 :: IAm
!    IAm = 'run_aerosol_optics1'

    __Iam__('run_aerosol_optics1')


!   Begin...

!   Radiation band
!   --------------
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, RH, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, PLE, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)


    allocate(ext(i1:i2,j1:j2,km),  &
             ssa(i1:i2,j1:j2,km),  &
             asy(i1:i2,j1:j2,km), __STAT__)


!   Get list of child states within state
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes,__RC__)

    b=0
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
        end if
    end do

    allocate (AEROlist(b), __STAT__)

    j = 1
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            AEROlist(j) = trim(itemList(i))
            j = j + 1
        end if
    end do

    ext = 0.0d0
    ssa = 0.0d0
    asy = 0.0d0


!   do i = 1, size(AEROlist)
   do i = 1, 1
        call ESMF_StateGet(state, trim(AEROlist(i)), child_state, __RC__)

!       ! set RH for aerosol optics
        call ESMF_AttributeGet(child_state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, AS_PTR_3D, trim(fld_name), __RC__)
            AS_PTR_3D = RH
        end if

!       ! set PLE for aerosol optics
        call ESMF_AttributeGet(child_state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, AS_PTR_3D, trim(fld_name), __RC__)
            AS_PTR_3D = PLE
        end if

        call ESMF_AttributeSet(child_state, name='band_for_aerosol_optics', value=band, __RC__)

!       ! execute the aero provider's optics method
        call ESMF_MethodExecute(child_state, label="aerosol_optics", __RC__)

!       ! EXT from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ext_, trim(fld_name), __RC__)
        end if

!       ! SSA from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ssa_, trim(fld_name), __RC__)
        end if

!       ! ASY from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, asy_, trim(fld_name), __RC__)
        end if

        ext = ext + ext_
        ssa = ssa + ssa_
        asy = asy + asy_
    end do



    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy(:,:,:)
    end if

    deallocate(ext, ssa, asy, __STAT__)



   RETURN_(ESMF_SUCCESS)

  end subroutine run_aerosol_optics1

!-----------------------------------------------------------------------------------------------

end module  Chem_AeroGeneric


