#include "unused_dummy.H"
!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: m_topo_remap --- Wrapper to remap dyn vector due to topo change
!
! !INTERFACE:
module m_topo_remap

! !USES:

use m_dyn
use shared_topo_remap, only : dyn_topo_remap, gmap
implicit none
private

! !PUBLIC MEMBER FUNCTIONS:
public dyn_topo_remap  ! nothing else should be made public
public dyn_real_eta    ! nothing else should be made public

!
interface dyn_topo_remap
  module procedure dyn_topo_remap_
end interface
interface dyn_real_eta
  module procedure dyn_real_eta_
end interface

character(len=*), parameter :: myname = 'm_topo_remap'

! !DESCRIPTION: This module is a wrapper based on Larry Takacs wrapper
!               of S. J. Lin routine to remap dyn-vector give a change
!               in the topography.
! !REMARKS:
!   1. Nothing beyond dyn_topo_remap should ever be made public to 
!      to avoid conflict with similar routines elsewhere.
!   2. In all of these, wind-handling needs to be revisited since
!      SJ original routines are for d-grid winds (far as I know).
!
! !SEE ALSO: m_mapz, m_maph
!
! !REVISION HISTORY:
!
!  23Feb2013  Todling    Initial (wrapper) code.
!  20Nov2013  Todling    Add interface for converting to real-eta
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
      subroutine dyn_topo_remap_( w_f,phis_new, dyntype, info )
      implicit none
      integer,intent(in) :: dyntype
      type(dyn_vect) w_f
      real,intent(in) :: phis_new(w_f%grid%im,w_f%grid%jm)
      integer,optional::info

      character(len=*), parameter :: myname_ = myname//'*dyn_topo_remap'
      real,allocatable ::  pk(:,:,:)
      real,allocatable :: pke(:,:,:)
      real,allocatable :: ple(:,:,:)
      real,allocatable :: thv(:,:,:)
      real kappa
      integer k,im,jm,km

      kappa = 2.0/7.0
      im=w_f%grid%im
      jm=w_f%grid%jm
      km=w_f%grid%km

!     check to see if remap is needed
!     only remap when max diff in topo larger than 1% of max value of topo
      if(maxval(abs(phis_new-w_f%phis))<0.01*maxval(w_f%phis)) then
        print*, myname_, ': no remap necessary, returning ...'
        if (present(info)) then
            info = 1
        endif
        return
      endif
      if (present(info)) then
         info = 0
      endif
      allocate(thv(im,jm,km))
      if(dyntype==5) then
         allocate ( pk(im,jm,km  ))
         allocate (ple(im,jm,km+1))
         allocate (pke(im,jm,km+1))
         ple(:,:,1) = w_f%grid%ak(1)
         do k=1,km
            ple(:,:,k+1) = ple(:,:,k) + w_f%delp(:,:,k)
         enddo
       ! convert virtual temperature to virtual potential temperature
         pke(:,:,:) = ple(:,:,:)**kappa
         do k=1,km
            pk(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                       / ( kappa*log(ple(:,:,k+1)/ple(:,:,k)) )
         enddo
         thv = w_f%pt/pk
         deallocate (  pk )
         deallocate ( pke )
         deallocate ( ple )
      else 
          print *, ' not coded for old dyn'
          call exit(1)
          stop
      endif

      call dyn_topo_remap( w_f%ps,w_f%delp,w_f%u,w_f%v,thv,w_f%q,w_f%phis,phis_new, &
                         w_f%grid%ak,w_f%grid%bk,&
                         im,jm,km,4)

      if(dyntype==5) then
         allocate (ple(im,jm,km+1))
         allocate (pke(im,jm,km+1))
         allocate ( pk(im,jm,km  ))
         ple(:,:,1) = w_f%grid%ak(1)
         do k=1,km
            ple(:,:,k+1) = ple(:,:,k) + w_f%delp(:,:,k)
         enddo
         pke(:,:,:) = ple(:,:,:)**kappa
         do k=1,km
            pk(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                       / ( kappa*log(ple(:,:,k+1)/ple(:,:,k)) )
         enddo
       ! convert virtual potential temperature to virtual temperature
         w_f%pt = thv*pk
         deallocate (  pk )
         deallocate ( pke )
         deallocate ( ple )
      endif
      w_f%phis = phis_new ! for consistency

      deallocate(thv)
      end subroutine dyn_topo_remap_

      subroutine dyn_real_eta_( w_f, dyntype, info )
      implicit none
      integer,intent(in) :: dyntype
      type(dyn_vect) w_f
      integer,optional::info

      character(len=*), parameter :: myname_ = myname//'*dyn_topo_remap'
      real,allocatable ::  pk(:,:,:)
      real,allocatable :: pke(:,:,:)
      real,allocatable :: ple_cur(:,:,:)
      real,allocatable :: ple_new(:,:,:)
      real,allocatable :: thv(:,:,:)
      real kappa
      integer k,im,jm,km

      kappa = 2.0/7.0
      im=w_f%grid%im
      jm=w_f%grid%jm
      km=w_f%grid%km

!     check to see if remap is needed
!     only remap when max diff in ps larger than 1% of max value of ps
!     if(maxval(abs(ps_new-w_f%ps))<0.01*maxval(w_f%ps)) then
!       print*, myname_, ': no remap necessary, returning ...'
!       if (present(info)) then
!           info = 1
!       endif
!       return
!     endif
      if (present(info)) then
         info = 0
      endif
      allocate(ple_cur(im,jm,km+1))
      allocate(ple_new(im,jm,km+1))
      allocate(    thv(im,jm,km))
      if(dyntype==5) then
         allocate ( pk(im,jm,km  ))
         allocate (pke(im,jm,km+1))
         ple_cur(:,:,1) = w_f%grid%ak(1)
         ple_new(:,:,1) = w_f%grid%ak(1)
         do k=1,km
            ple_cur(:,:,k+1) = ple_cur(:,:,k) + w_f%delp(:,:,k)
            ! overwrite delp in w_f and force it to be eta-consistent
            w_f%delp(:,:,k) = w_f%grid%ak(k+1)-w_f%grid%ak(k) + &
                             (w_f%grid%bk(k+1)-w_f%grid%bk(k))*w_f%ps
            ple_new(:,:,k+1) = ple_new(:,:,k) + w_f%delp(:,:,k)
         enddo
       ! convert virtual temperature to virtual potential temperature
         pke(:,:,:) = ple_cur(:,:,:)**kappa
         do k=1,km
            pk(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                       / ( kappa*log(ple_cur(:,:,k+1)/ple_cur(:,:,k)) )
         enddo
         thv = w_f%pt/pk
         deallocate (  pk )
         deallocate ( pke )
      else 
          print *, ' not coded for old dyn'
          call exit(1)
          stop
      endif

      call ps_remap0_( ple_cur,ple_new,w_f%u,w_f%v,thv,w_f%q, &
                       w_f%grid%ak,w_f%grid%bk,&
                       im,jm,km,4)

      if(dyntype==5) then
         allocate (pke(im,jm,km+1))
         allocate ( pk(im,jm,km  ))
       ! convert virtual potential temperature to virtual temperature
         pke(:,:,:) = ple_new(:,:,:)**kappa
         do k=1,km
            pk(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                       / ( kappa*log(ple_new(:,:,k+1)/ple_new(:,:,k)) )
         enddo
         w_f%pt = thv*pk
         deallocate (  pk )
         deallocate ( pke )
      endif

      deallocate(thv)
      deallocate(ple_new)
      deallocate(ple_cur)
      end subroutine dyn_real_eta_

      subroutine ps_remap0_( ple,ple_out,u,v,thv,q,ak,bk,im,jm,lm,nq )

!***********************************************************************
!
!  Purpose
!     Driver for remapping fields lcv fields to consistent eta-coordinate
!
!  Argument Description
!     ple ...... model edge pressure
!     u  ....... model zonal      wind
!     v  ....... model meridional wind
!     thv  ..... model virtual potential  temperature
!     q  ....... model specific   humidity; ozone; others
!     ak ....... model vertical   dimension
!     bk ....... model vertical   dimension
!     ple_out .. target pressure levels
!
!     im ....... zonal      dimension
!     jm ....... meridional dimension
!     lm ....... meridional dimension
!     nq ....... number of tracers including spec. hum.
!
! 20Oct2013  Todling  Addapted for dyn-vect needs
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************



      implicit none
      integer  im,jm,lm,nq

! Input variables
! ---------------
      real      ple(im,jm,lm+1)
      real  ple_out(im,jm,lm+1)
      real        u(im,jm,lm)
      real        v(im,jm,lm)
      real      thv(im,jm,lm)
      real        q(im,jm,lm,nq)

      real       ak(lm+1)
      real       bk(lm+1)

! Local variables
! ---------------
      real, allocatable ::  pke    (:,:,:)
      real, allocatable ::  pke_out(:,:,:)

      real, allocatable ::    u_out(:,:,:)
      real, allocatable ::    v_out(:,:,:)
      real, allocatable ::  thv_out(:,:,:)
      real, allocatable ::    q_in (:,:,:,:)
      real, allocatable ::    q_out(:,:,:,:)

      real    kappa,cp,rgas,eps,rvap

      _UNUSED_DUMMY(ak)
      _UNUSED_DUMMY(bk)

      kappa = 2.0/7.0
      rgas  = 8314.3/28.97
      rvap  = 8314.3/18.01
      eps   = rvap/rgas-1.0
      cp    = rgas/kappa

      allocate(  pke    (im,jm,lm+1) )
      allocate(  pke_out(im,jm,lm+1) )

      allocate(    u_out(im,jm,lm)   )
      allocate(    v_out(im,jm,lm)   )
      allocate(  thv_out(im,jm,lm)   )
      allocate(    q_in (im,jm,lm,nq))
      allocate(    q_out(im,jm,lm,nq))

! Construct fv pressure variables using new surface pressure
! ----------------------------------------------------------
      pke(:,:,:)     = ple    (:,:,:)**kappa
      pke_out(:,:,:) = ple_out(:,:,:)**kappa

! Map original fv state onto new eta grid
! ---------------------------------------
      q_in(:,:,:,:) =  q(:,:,:,:)

      call gmap( im,jm,nq, kappa, &
                  lm, pke    ,ple    ,u    ,v    ,thv    ,q_in , &
                  lm, pke_out,ple_out,u_out,v_out,thv_out,q_out)

        u(:,:,:)   =   u_out(:,:,:)
        v(:,:,:)   =   v_out(:,:,:)
      thv(:,:,:)   = thv_out(:,:,:)
        q(:,:,:,:) =   q_out(:,:,:,:)

      deallocate(  pke_out )
      deallocate(  pke     )

      deallocate(    q_out )
      deallocate(    q_in  )
      deallocate(  thv_out )
      deallocate(    v_out )
      deallocate(    u_out )

      return
      end subroutine ps_remap0_

end module m_topo_remap
