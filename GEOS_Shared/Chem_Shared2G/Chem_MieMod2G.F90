
#include "MAPL_Exceptions.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE:  Chem_MieMod2G --- Load and manipulate Mie tables. 2nd Generation.
!                             ESMF Compliant.
!
! !INTERFACE:
!
   module  Chem_MieMod2G

! !USES:
   use ESMF
   use MAPL_Mod
   use Chem_MieTableMod1G

   use m_chars, only : uppercase
   implicit none

! !PUBLIC TYPES:
   private
   public  Chem_Mie                ! Holds Mie Lookup Tables

!
! !PUBLIC MEMBER FUNCTIONS:
   public  Chem_MieCreateng        ! Constructor for GOCART2G. Does not use Chem_Registry
   public  Chem_MieQuery           ! Query the Mie table to return parameters (qname interface)
   public  Chem_MieQueryIdx        ! Query the index of the mie table given the qname


!
! !DESCRIPTION:
!
!  This module read the mie aerosol tables.
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

! Mie LUT table
! Will be reduced from input files to the desired channels
! --------
  type Chem_Mie
     integer :: nch                               ! number of channels
     integer :: nMom=0                            ! number of moments (phase function)
     integer :: nPol=0                            ! number of moments (phase function)
     real, pointer    :: channels(:)              ! wavelengths

     character(len=255) :: rcfile
     character(len=255) :: optics_file

                                                  ! mie tables -- dim(nch,nrh,nbin)
     type(Chem_MieTable), pointer :: mie_aerosol  => null()

     integer :: nq                                ! number of tracers
     character(len=255), pointer  :: vname(:)  => null()      ! possibly remove lines 67-71 vname,vindex,vtable
     integer, pointer             :: vindex(:) => null()
     type(Chem_MieTable), pointer :: vtable(:) => null()
                                                  ! mapping of vtable for given idx
     type(Chem_MieTable), pointer :: vtableUse => null()
  end type Chem_Mie


  interface Chem_MieQuery
     module procedure Chem_MieQueryByInt
     module procedure Chem_MieQueryByChar  !can remove
     module procedure Chem_MieQueryByIntWithpmom  !possibly remove
  end interface


contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieCreateng --- Construct Mie LUTs from CF object for GOCARTng.  !CHANGE ng to 2G
!
! !INTERFACE:
!
  function Chem_MieCreateng ( cf, NUM_BANDS, rc ) result(this)

! !INPUT PARAMETERS:
   type (ESMF_Config)             :: cf          ! Mie table file name
   integer                        :: NUM_BANDS   ! Number of wavelengths/radiation bands 

! !OUTPUT PARAMETERS:
   type (Chem_Mie) this
   integer, intent(out) ::  rc

! !DESCRIPTION:
!
!  This routine creates a LUT object from an ESMF configuration
!  attribute CF. This routine is usually called from GEOS-5.
!
! !REVISION HISTORY:
!
!  09Mar2005 da Silva  API, prologues.
!  25Nov2019 Sherman - No longer uses Chem_Reg.
!
!EOP
!-------------------------------------------------------------------------

   type (ESMF_Config)    :: cfg
   integer               :: iq, i, nCols
   __Iam__('Chem_MieCreateng')


!  Set up the hash table to map the variable names to the
!  corresponding Mie Table
   call ESMF_ConfigGetDim( CF, this%nq, nCols, label=('variable_table::'), __RC__ )
   allocate(this%vname(this%nq), this%vindex(this%nq), __STAT__ )
   allocate(this%vtable(this%nq), __STAT__ )

   call ESMF_ConfigFindLabel( CF, 'variable_table::', __RC__ )
   do iq = 1, this%nq
      this%vindex(iq) = iq
      call ESMF_ConfigNextLine( CF, __RC__ )
      call ESMF_ConfigGetAttribute( CF, this%vname(iq), __RC__ )
   enddo


!  Get file names for the optical tables
!  -------------------------------------
   call ESMF_ConfigGetAttribute( CF, this%optics_file, Label="radiation_optics:", __RC__ )

!  Set the number of bands and channels
!  -------------------------------------
  this%nch = NUM_BANDS

!  Make chanel = number of bands
!  --------------------------------------------------------------------
   allocate ( this%channels(this%nch), __STAT__ )

   do i = 1, this%nch
       this%channels(i) = i
   end do

   allocate(this%mie_aerosol, __STAT__)
   this%mie_aerosol = Chem_MieTableCreate( this%optics_file, __RC__ )
   call Chem_MieTableRead( this%mie_aerosol, this%nch, this%channels, __RC__)


!  Now map the mie tables to the hash table    !remove this%vtable(iq)....
!  -----------------------------------------
   do iq = 1, this%nq
       this%vtable(iq) = this%mie_aerosol
   end do


!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

 end function Chem_MieCreateng
!-----------------------------------------------------------------------------------

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQueryIdx --- Return the index of the mie table given
!                                  a qname requested
!
!
! !INTERFACE:
!
   Function Chem_MieQueryIdx ( this, qname, rc ) result(idx)

   implicit none

! !INPUT PARAMETERS:

   type(Chem_Mie), intent(inout) :: this   ! Input mie table structure
   character(len=*), intent(in)  :: qname  ! Variable name to find in table, e.g., du001

! !OUTPUT PARAMETERS:

   integer, optional, intent(out) ::  rc ! Error return code:
                                         !  0 - all is well
                                         !  1 - 
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!   24Apr2006, PRC
!
!EOP
!-------------------------------------------------------------------------
      character(len=255) :: NAME
      integer            :: idx         ! Index number in Mie table of qname
      integer            :: iq, i

!     Find the right table for this aerosol from its name

      NAME = trim(qname)
!     Remove qualifier from variable name: GOCART::du001 --> du001
!     ------------------------------------------------------------
      i = index(NAME,'::')
      if ( i > 0 ) then
         NAME = NAME(i+2:)
      end if

      idx = -1
      do iq = 1, this%nq
       if((trim((NAME))) .eq. (trim(this%vname(iq)))) then
        idx = this%vindex(iq)
        this%vtableUse => this%vtable(iq)
        exit
       endif
      enddo

      if(present(rc)) then
         if(idx .eq. -1) then
            rc = 1
         else
            rc = 0
         end if
      end if

      return

  end Function Chem_MieQueryIdx
!----------------------------------------------------------------------------------

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQueryByInt --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!
   impure elemental subroutine Chem_MieQueryByInt ( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   reff, p11, p22, gf, rhop, rhod, &
                                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

   type(Chem_Mie), target, intent(in ) :: this
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

! !OUTPUT PARAMETERS:

   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   real,    optional,      intent(out) :: gf    ! Growth factor (ratio of wet to dry radius)
   real,    optional,      intent(out) :: rhop  ! Wet particle density [kg m-3]
   real,    optional,      intent(out) :: rhod  ! Dry particle density [kg m-3]
   real,    optional,      intent(out) :: vol   ! Wet particle volume [m3 kg-1]
   real,    optional,      intent(out) :: area  ! Wet particle cross section [m2 kg-1]
   real,    optional,      intent(out) :: refr  ! Wet particle real part of ref. index
   real,    optional,      intent(out) :: refi  ! Wet particle imag. part of ref. index
   integer, optional,      intent(out) :: rc    ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function 
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------


      integer                      :: ICHANNEL, TYPE
      integer                      :: irh, irhp1, isnap
      real                         :: rhUse, arh
      real                         :: bextIn, bscaIn, bbckIn, gasymIn, p11In, p22In, &
                                      gfIn, rhopIn, rhodIn, volIn, areaIn, &
                                      refrIn, refiIn
      type(Chem_MieTable), pointer :: TABLE

      character(len=*), parameter  :: Iam = 'Chem_MieQuery'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TABLE => this%vtableUse
      TYPE = idx

!      ASSERT_(TYPE>0)
!      ASSERT_(ICHANNEL>=LBOUND(TABLE%bext,1))
!      ASSERT_(ICHANNEL<=UBOUND(TABLE%bext,1))

!     Now map the input RH to the high resolution hash table for RH
      rhUse = max(rh,0.)
      rhUse = min(rh,0.99)
      isnap = int((rhUse+0.001)*1000.)
      if(isnap .lt. 1) isnap = 1
      arh   = TABLE%rha( isnap )
      irh   = TABLE%rhi( isnap )
      irhp1 = irh+1
      if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.

      if(present(bext) .or. present(tau) .or. present(ssa) ) then
         bextIn =   TABLE%bext(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bext(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bsca) .or. present(ssa) ) then
         bscaIn =   TABLE%bsca(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bsca(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bbck)) then
         bbckIn =   TABLE%bbck(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bbck(irhp1,ichannel,TYPE) * arh
      endif

      if(present(gasym)) then
         gasymIn =  TABLE%g(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%g(irhp1,ichannel,TYPE) * arh
      endif

      if(present(rEff) ) then
         rEff =     TABLE%rEff(irh  ,TYPE) * (1.-arh) &
                  + TABLE%rEff(irhp1,TYPE) * arh
         rEff = 1.E6 * rEff ! convert to microns
      endif

!      if(present(pmom)) then
!         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
!                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh
!      endif


!      if(present(pmom)) then
!         call Chem_MieQueryByIntWithpmom(this, idx, channel, q_mass, rh, pmom)
!      endif


      if(present(p11) ) then
         p11In =   TABLE%pback(irh  ,ichannel,TYPE,1) * (1.-arh) &
                 + TABLE%pback(irhp1,ichannel,TYPE,1) * arh
      endif

      if(present(p22) ) then
         p22In =   TABLE%pback(irh  ,ichannel,TYPE,5) * (1.-arh) &
                 + TABLE%pback(irhp1,ichannel,TYPE,5) * arh
      endif

      if(present(gf) ) then
         gfIn =     TABLE%gf(irh  ,TYPE) * (1.-arh) &
                  + TABLE%gf(irhp1,TYPE) * arh
      endif

      if(present(rhod) ) then
         rhodIn =   TABLE%rhod(1  ,TYPE)
      endif

      if(present(vol) ) then
         volIn  =   TABLE%vol(irh  ,TYPE) * (1.-arh) &
                  + TABLE%vol(irhp1,TYPE) * arh
      endif

      if(present(area) ) then
         areaIn  =   TABLE%area(irh  ,TYPE) * (1.-arh) &
                  + TABLE%area(irhp1,TYPE) * arh
      endif

      if(present(refr) .or. present(tau) .or. present(ssa) ) then
         refrIn =   TABLE%refr(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%refr(irhp1,ichannel,TYPE) * arh
      endif

      if(present(refi) .or. present(tau) .or. present(ssa) ) then
         refiIn =   TABLE%refi(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%refi(irhp1,ichannel,TYPE) * arh
      endif

!     Fill the requested outputs
      if(present(tau  )) tau   = bextIn * q_mass
      if(present(ssa  )) ssa   = bscaIn/bextIn
      if(present(bext )) bext  = bextIn
      if(present(bsca )) bsca  = bscaIn
      if(present(bbck )) bbck  = bbckIn
      if(present(gasym)) gasym = gasymIn
      if(present(p11  )) p11   = p11In
      if(present(p22  )) p22   = p22In
      if(present(gf   )) gf    = gfIn
      if(present(rhop )) rhop  = rhopIn
      if(present(rhod )) rhod  = rhodIn
      if(present(vol ))  vol   = volIn
      if(present(area )) area  = areaIn
      if(present(refr )) refr  = refrIn
      if(present(refi )) refi  = refiIn

!  All Done
!----------
  end subroutine Chem_MieQueryByInt



  subroutine Chem_MieQueryByChar( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   rEff, pmom, p11, p22, rc )

!  ! INPUT parameters
   type(Chem_Mie), target, intent(in ) :: this     
   character(*),           intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

!  ! OUTPUT Parameters
   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,    optional,      intent(out) :: pmom(:,:)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   integer, optional,      intent(out) :: rc    ! error code

   integer :: iq, i

   character(len=*), parameter  :: Iam = 'Chem_MieQueryByChar'
   character(len=255) :: NAME

   if ( present(rc) ) rc = 0

!  Remove qualifier from variable name: GOCART::du001 --> du001
!   ------------------------------------------------------------
   NAME = trim(idx)
   i = index(NAME,'::')
   if ( i > 0 ) then
      NAME = NAME(i+2:)
   end if

   do iq = 1, this%nq
      if( uppercase(trim(NAME)) == uppercase(trim(this%vname(iq)))) then
         call  Chem_MieQueryByInt( this, iq, channel, q_mass, rh,     &
                             tau, ssa, gasym, bext, bsca, bbck, &
                             rEff, p11, p22, rc=rc )
         if ( rc /= 0 ) return
      endif
   enddo

 end subroutine Chem_MieQueryByChar


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQueryByIntWithpmom --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!
   subroutine Chem_MieQueryByIntWithpmom ( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   reff, pmom, p11, p22, gf, rhop, rhod, &
                                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

   type(Chem_Mie), target, intent(in ) :: this
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

! !OUTPUT PARAMETERS:

   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,                   intent(out) :: pmom(:,:)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   real,    optional,      intent(out) :: gf    ! Growth factor (ratio of wet to dry radius)
   real,    optional,      intent(out) :: rhop  ! Wet particle density [kg m-3]
   real,    optional,      intent(out) :: rhod  ! Dry particle density [kg m-3]
   real,    optional,      intent(out) :: vol   ! Wet particle volume [m3 kg-1]
   real,    optional,      intent(out) :: area  ! Wet particle cross section [m2 kg-1]
   real,    optional,      intent(out) :: refr  ! Wet particle real part of ref. index
   real,    optional,      intent(out) :: refi  ! Wet particle imag. part of ref. index
   integer, optional,      intent(out) :: rc    ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function 
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------


      integer                      :: ICHANNEL, TYPE
      integer                      :: irh, irhp1, isnap
      real                         :: rhUse, arh
      type(Chem_MieTable), pointer :: TABLE

      character(len=*), parameter  :: Iam = 'Chem_MieQueryByIntWithpmom'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TABLE => this%vtableUse
      TYPE = idx

!     Now map the input RH to the high resolution hash table for RH
      rhUse = max(rh,0.)
      rhUse = min(rh,0.99)
      isnap = int((rhUse+0.001)*1000.)
      if(isnap .lt. 1) isnap = 1
      arh   = TABLE%rha( isnap )
      irh   = TABLE%rhi( isnap )
      irhp1 = irh+1
      if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.
    call Chem_MieQuery ( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   reff, p11, p22, gf, rhop, rhod, &
                                   vol, area, refr, refi, rc )

         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh



!  All Done
!----------
  end subroutine Chem_MieQueryByIntWithpmom







 end module Chem_MieMod2G


