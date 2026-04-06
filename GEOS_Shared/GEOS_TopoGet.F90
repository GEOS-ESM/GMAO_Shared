! $Id$
! VERIFY_ and RETURN_ macros for error handling

#include "MAPL.h"

!BOP

! !MODULE: GEOS_TopoGetMod -- The topography computation

! !INTERFACE:

module GEOS_TopoGetMod

! !USES:

  use ESMF
  use MAPL_ConstantsMod, only: MAPL_GRAV
  use MAPL_ErrorHandlingMod
  use mapl3g_generic, only: MAPL_GridCompGetResource
  
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

!BOP

! !IROUTINE GEOS_TopoGet -- Gets Topographic Variables

! !INTERFACE

  subroutine GEOS_TopoGet ( gc,                          &
                            MEAN,                        &
                            GWDVAR,   GWDVARX,  GWDVARY, &
                            GWDVARXY, GWDVARYX, TRBVAR, RC)

! !ARGUMENTS

    type(ESMF_GridComp)                          :: gc
    type(ESMF_Field),    optional, intent(INOUT) :: MEAN
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVAR
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARX
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARY
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARXY
    type(ESMF_Field),    optional, intent(INOUT) :: GWDVARYX
    type(ESMF_Field),    optional, intent(INOUT) :: TRBVAR
    integer,             optional, intent(OUT)   :: RC

! !DESCRIPTION

!  This subroutine creates topographic data associated with an input
!  ESMF grid.  The available topographic data types are:  MEAN, GWDVAR,
!  GWDVARX, GWDVARY, GWDVARXY, GWDVARYX, and TRBVAR.  The raw data
!  for each of these types has been pre-processed and stored at 2.5'x2.5'
!  resolution.  The resulting gridded data will be binned-averaged on the
!  input ESMF grid from the 2.5'x2.5' data.
!  The arguments are:
!
! \begin{description}
!   \item[GRID]
!                   The ESMF GRID which contains information about 
!                   horizontal grid structure.
!   \item[MEAN]
!                   The mean values of topography with scales >= 100 km.
!   \item[GWDVAR]
!                   The isotropic variance of the GWD topography data,
!                   (scales 10-100 km).
!   \item[GWDVARX]
!                   The variance of the GWD topography data in the
!                   East - West direction.
!   \item[GWDVARY]
!                   The variance of the GWD topography data in the
!                   North - South direction.
!   \item[GWDVARXY]
!                   The variance of the GWD topography data in the
!                   South_West - North_East direction.
!   \item[GWDVARYX]
!                   The variance of the GWD topography data in the
!                   North_West - South_East direction.
!   \item[TRBVAR]
!                   The isotropic variance of the Turbulence topography data,
!                   (scales 1-10 km).
!   \item[RC]
!                   Return code
! \end{description}
!

!EOP

! Locals

    character(len=ESMF_MAXSTR), parameter :: IAm = "TopoGet"
    ! Define a type to hold variable-length strings
    type :: string_type
        character(len=:), allocatable :: str
    end type string_type
    type(string_type)           :: filename(7)

    real, pointer :: ptr(:,:)
    integer       :: STATUS
    integer       :: unit
    real          :: GWDFAC
    real          :: GWDFACX
    real          :: GWDFACY
    real          :: GWDFACXY
    real          :: GWDFACYX
    real          :: TRBFAC

! Get filenames for Get_Topo utility
! ----------------------------------
    call MAPL_GridCompGetResource ( gc, "TOPO_MEAN_FILE", filename(1)%str, &
                                   default="hmean.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_GWDVAR_FILE", filename(2)%str, &
                                   default="hgrav_var.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_GWDVARX_FILE", filename(3)%str, &
                                   default="hgrav_varx.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_GWDVARY_FILE", filename(4)%str,  &
                                   default="hgrav_vary.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_GWDVARXY_FILE", filename(5)%str, &
                                   default="hgrav_varxy.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_GWDVARYX_FILE", filename(6)%str, &
                                   default="hgrav_varyx.2.5x2.5min.data", _RC)

    call MAPL_GridCompGetResource ( gc, "TOPO_TRBVAR_FILE", filename(7)%str, &
                                   default="hturb_var.2.5x2.5min.data", _RC)

  if( present(MEAN)  ) then
! -------------------------
       open(newunit=unit, file=filename(1)%str, form="unformatted")
       call MAPL_VarRead (UNIT,MEAN)
       close    (UNIT)
       call ESMF_FieldGet(MEAN, 0, PTR, _RC)
       ptr = ptr*MAPL_GRAV
  endif

  if( present(GWDVAR) ) then
! --------------------------
       call MAPL_GridCompGetResource( gc, "GWDVAR_FACTOR", GWDFAC, default = 1.0, _RC)
       open(newunit=unit, file=filename(2)%str, form="unformatted")
       call MAPL_VarRead (UNIT,GWDVAR)
       close    (UNIT)
       call ESMF_FieldGet (GWDVAR, 0, PTR, _RC)
       ptr = sqrt( max(gwdfac*ptr,0.0) )
  endif

  if( present(GWDVARX) ) then
! ---------------------------
       call MAPL_GridCompGetResource( gc, "GWDVARX_FACTOR", GWDFACX, default = 1.0, _RC)
       open(newunit=unit, file=filename(3)%str, form="unformatted")
       call MAPL_VarRead (UNIT,GWDVARX)
       close    (UNIT)
       call ESMF_FieldGet (GWDVARX, 0, PTR, _RC)
       ptr = sqrt( max(gwdfacx*ptr,0.0) )
  endif

  if( present(GWDVARY) ) then
! ---------------------------
       call MAPL_GridCompGetResource( gc, "GWDVARY_FACTOR", GWDFACY, default = 1.0, _RC)
       open(newunit=unit, file=filename(4)%str, form="unformatted")
       call MAPL_VarRead (UNIT,GWDVARY)
       close    (UNIT)
       call ESMF_FieldGet (GWDVARY, 0, PTR, _RC)
       ptr = sqrt( max(gwdfacy*ptr,0.0) )
  endif

  if( present(GWDVARXY) ) then
! ----------------------------
       call MAPL_GridCompGetResource( gc, "GWDVARXY_FACTOR", GWDFACXY, default = 1.0, _RC)
       open(newunit=unit, file=filename(5)%str, form="unformatted")
       call MAPL_VarRead (UNIT,GWDVARXY)
       close    (UNIT)
       call ESMF_FieldGet (GWDVARXY, 0, PTR, _RC)
       ptr = sqrt( max(gwdfacxy*ptr,0.0) )
  endif

  if( present(GWDVARYX) ) then
! ----------------------------
       call MAPL_GridCompGetResource( gc, "GWDVARYX_FACTOR", GWDFACYX, default = 1.0, _RC)
       open(newunit=unit, file=filename(6)%str, form="unformatted")
       call MAPL_VarRead (UNIT,GWDVARYX)
       close    (UNIT)
       call ESMF_FieldGet (GWDVARYX, 0, PTR, _RC)
       ptr = sqrt( max(gwdfacyx*ptr,0.0) )
  endif

  if( present(TRBVAR) ) then
! --------------------------
       call MAPL_GridCompGetResource( gc, "TRBVAR_FACTOR", TRBFAC, default = 1.0, _RC)
       open(newunit=unit, file=filename(7)%str, form="unformatted")
       call MAPL_VarRead (UNIT,TRBVAR)
       close    (UNIT)
       call ESMF_FieldGet (TRBVAR, 0, PTR, _RC)
       ptr = max(trbfac*ptr,0.0)
  endif

    return
  end subroutine GEOS_TopoGet

end module GEOS_TopoGetMod
