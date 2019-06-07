!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_geosapi - GEOS interfaces external to this component
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_geosapi
      implicit none
      private	! except

      public :: getcon		! get GEOS3 constants

    interface getcon
      function getcon(name)
	implicit none
        character(len=*),intent(in) :: name
	real :: getcon
      end function getcon
    end interface

! !REVISION HISTORY:
! 	12Nov04	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_geosapi'
end module m_geosapi
