module atmOcnIntlayer

! !USES:

use MAPL

implicit none
private
public SIMPLE_SW_ABS

contains

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
