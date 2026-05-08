! $Id$
#include "MAPL_ErrLog.h"

!BOP

! !MODULE: GEOS_TopoGetMod -- The topography computation

! !INTERFACE:

module GEOS_TopoGetMod

! !USES:

  use ESMF
  use MAPL
  use MAPL2, only: MAPL_VarRead  ! TODO: move NCIOMod to mapl_base3g so this comes from MAPL
  
  implicit none
  private

! !PUBLIC ROUTINES:

  public GEOS_TopoGet

! !DESCRIPTION:

!  This module computes the Earth's topography associated with
!  an input ESMF grid.  The mean height topography was averaged
!  from the 30"x30" GTOPO30 dataset obtained from:
!
!          http://edcdaac.usgs.gov/gtopo30/gtopo30.html 
!
!  using a box stencil with gaussian weights, retaining scales >= 100 km.
!  In addition to mean heights, the routine can also return isotropic and
!  directional variances associated with Gravity-Wave-Drag scales (10-100 km) 
!  and turbulence scales (0-10 km).

!EOP

contains

!==========================================================================

  ! Helper: get a string value from HConfig with a default
  subroutine hconfig_get_string(cf, label, value, default, rc)
    type(ESMF_HConfig), intent(in)  :: cf
    character(len=*),   intent(in)  :: label
    character(len=*),   intent(out) :: value
    character(len=*),   intent(in)  :: default
    integer, optional,  intent(out) :: rc
    integer :: status
    logical :: defined
    character(len=:), allocatable :: tmp
    value = default
    defined = ESMF_HConfigIsDefined(cf, keyString=label, _RC)
    if (defined) then
       tmp = ESMF_HConfigAsString(cf, keyString=label, _RC)
       value = tmp
    end if
    _RETURN(ESMF_SUCCESS)
  end subroutine hconfig_get_string

  ! Helper: get a real(4) value from HConfig with a default
  subroutine hconfig_get_r4(cf, label, value, default, rc)
    type(ESMF_HConfig), intent(in)  :: cf
    character(len=*),   intent(in)  :: label
    real,               intent(out) :: value
    real,               intent(in)  :: default
    integer, optional,  intent(out) :: rc
    integer :: status
    logical :: defined
    value = default
    defined = ESMF_HConfigIsDefined(cf, keyString=label, _RC)
    if (defined) value = ESMF_HConfigAsR4(cf, keyString=label, _RC)
    _RETURN(ESMF_SUCCESS)
  end subroutine hconfig_get_r4

!==========================================================================

!BOP

! !IROUTINE GEOS_TopoGet -- Gets Topographic Variables

! !INTERFACE

  subroutine GEOS_TopoGet ( cf,                          &
                            MEAN,                        &
                            GWDVAR,   GWDVARX,  GWDVARY, &
                            GWDVARXY, GWDVARYX, TRBVAR, RC)

! !ARGUMENTS

    type(ESMF_HConfig),                          intent(in   ) :: cf
    type(ESMF_Field),    optional, intent(INOUT) :: MEAN
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVAR
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARX
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARY
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARXY
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARYX
    type(ESMF_Field),    optional, intent(INOUT) :: TRBVAR
    integer,             optional, intent(OUT)   :: RC

!EOP

! Locals

    character(len=ESMF_MAXSTR), parameter :: IAm = "TopoGet"
    character(len=ESMF_MAXSTR)            :: filename(7)

    real, pointer :: ptr(:,:)
    integer       :: STATUS
    real          :: GWDFAC
    real          :: TRBFAC
    type(Netcdf4_Fileformatter)           :: formatter

! Unsupported optional variables
! -------------------------------
    _ASSERT(.not. present(GWDVARX),  "GWDVARX not yet supported in NetCDF")
    _ASSERT(.not. present(GWDVARY),  "GWDVARY not yet supported in NetCDF")
    _ASSERT(.not. present(GWDVARXY), "GWDVARXY not yet supported in NetCDF")
    _ASSERT(.not. present(GWDVARYX), "GWDVARYX not yet supported in NetCDF")

! Get filenames for Get_Topo utility
! ----------------------------------
    call hconfig_get_string(cf, 'TOPO_MEAN_FILE',   filename(1), default='hmean.2.5x2.5min.nc',     _RC)
    call hconfig_get_string(cf, 'TOPO_GWDVAR_FILE', filename(2), default='hgrav_var.2.5x2.5min.nc', _RC)
    call hconfig_get_string(cf, 'TOPO_TRBVAR_FILE', filename(3), default='hturb_var.2.5x2.5min.nc', _RC)

  if( present(MEAN) ) then
! -------------------------
       call formatter%open(trim(filename(1)), pFIO_READ, _RC)
       call ESMF_FieldGet(MEAN, 0, PTR, rc=status)
       _VERIFY(STATUS)
       call MAPL_VarRead(formatter, 'z', ptr, _RC)
       call formatter%close(_RC)
       ptr = ptr*MAPL_GRAV
  endif

  if( present(GWDVAR) ) then
! --------------------------
       call hconfig_get_r4(cf, 'GWDVAR_FACTOR', GWDFAC, default=1.0, _RC)
       call formatter%open(trim(filename(2)), pFIO_READ, _RC)
       call ESMF_FieldGet(GWDVAR, 0, PTR, rc=status)
       _VERIFY(STATUS)
       call MAPL_VarRead(formatter, 'gwd', ptr, _RC)
       call formatter%close(_RC)
       ptr = sqrt( max(gwdfac*ptr, 0.0) )
  endif

  if( present(TRBVAR) ) then
! --------------------------
       call hconfig_get_r4(cf, 'TRBVAR_FACTOR', TRBFAC, default=1.0, _RC)
       call formatter%open(trim(filename(3)), pFIO_READ, _RC)
       call ESMF_FieldGet(TRBVAR, 0, PTR, rc=status)
       _VERIFY(STATUS)
       call MAPL_VarRead(formatter, 'trb', ptr, _RC)
       call formatter%close(_RC)
       ptr = max(trbfac*ptr, 0.0)
  endif

    _RETURN(ESMF_SUCCESS)
  end subroutine GEOS_TopoGet

end module GEOS_TopoGetMod
