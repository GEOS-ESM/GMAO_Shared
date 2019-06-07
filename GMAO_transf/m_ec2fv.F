        module m_ec2fv
        use m_const, only : kappa
        implicit none
        PRIVATE
        PUBLIC ec2fv
        interface ec2fv
        module procedure ec2fv_
        end interface
        CONTAINS
        subroutine ec2fv_(im_e,jm_e,km_e,nq, pke1, pe1,   
     .                    u, v, pt, q,o,qitot,qltot,kn_e,pke2,pe2,  
     .                    pke0,pe0,u0, v0, thv0, q0,o0,qitot0,qltot0,
     .                    dgrid,pabove,pbelow,
     .                    wsh,blend)  ! optional argument


! **********************************************************************
!
!  Purpose
!     Driver for remapping of target analysis to fv model levels
!
! ! REVISION HISTORY:
!   12Jan2005    Todling  Updated interface to set_eta (now from m_mapz)
!   23Jan2006    B Zhang  Updated for total cloud water fraction
!   22July2009   Ravi     Added optional dgrid for GEOS-5
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************
        use m_gmap, only : gmap
        implicit none
        character(len=*), parameter :: myname = 'ec2fv'
        integer  im_e,jm_e,lm_e,km_e,kn_e,nq
        integer  i,j,L
        real pke1(im_e,jm_e,km_e+1),pe1(im_e,jm_e,km_e+1)
        real u(im_e,jm_e,km_e),v(im_e,jm_e,km_e),pt(im_e,jm_e,km_e)
        real q(im_e,jm_e,km_e),o(im_e,jm_e,km_e)
        real qitot(im_e,jm_e,km_e),qltot(im_e,jm_e,km_e)

        real pke0(im_e,jm_e,kn_e+1),pe0(im_e,jm_e,kn_e+1)
        real pke2(im_e,jm_e,kn_e+1),pe2(im_e,jm_e,kn_e+1)
        real u0(im_e,jm_e,kn_e),v0(im_e,jm_e,kn_e),thv0(im_e,jm_e,kn_e)
        real q0(im_e,jm_e,kn_e),o0(im_e,jm_e,kn_e)
        real qitot0(im_e,jm_e,kn_e),qltot0(im_e,jm_e,kn_e)

        real u1(im_e,jm_e,kn_e),v1(im_e,jm_e,kn_e),thv1(im_e,jm_e,kn_e)
        real q1(im_e,jm_e,kn_e),o1(im_e,jm_e,kn_e)
        real qitot1(im_e,jm_e,kn_e),qltot1(im_e,jm_e,kn_e)
        real, intent(in), optional  :: pabove, pbelow
        real alf,pl
        logical,optional:: wsh,blend
 

        logical :: dgrid               ! default: in GEOS-4 dyn-vector winds are on D-grid


! put CERES winds on 'A' grid before remapping
                                                                                                            

        if(dgrid) then
         call dtoa ( u0,u0,im_e,jm_e,kn_e,2 )
         call dtoa ( v0,v0,im_e,jm_e,kn_e,1 )
        endif

        u1   = u0
        v1   = v0
        thv1 = thv0
        q1   = q0
        o1   = o0
        qitot1   = qitot0
        qltot1   = qltot0
c
! Map original fv state onto new eta grid
!        print *,' im_e,jm_e,km_e,kn_e ',im_e,jm_e,km_e,kn_e

         call gmap ( im_e,jm_e,nq, kappa,     
     .               kn_e, pke2, pe2, u0, v0, thv0, q0,
     .               kn_e, pke0, pe0, u1, v1, thv1, q1,
     .               o_m=o1, w_m=qitot1,wl_m=qltot1,o_n=o0 ,
     .               w_n=qitot0, wl_n=qltot0)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of U-wind  GMAO"
!        print *, '  '
!        call minmax_uv(u1,im_e,jm_e,kn_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of V-wind  GMAO"
!        print *, '  '
!        call minmax_uv(v1,im_e,jm_e,kn_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of thv1  GMAO"
!        print *, '  '
!        call minmax_uv(thv1,im_e,jm_e,kn_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of q1  GMAO"
!        print *, '  '
!        call minmax_uv(q1,im_e,jm_e,kn_e)


!     --------------------------------------------------------
!       Call a routine to remap ECMWF 60 or 91 level variables to
!       fv 55 or 72 level variable.
!     --------------------------------------------------------
!

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of U-wind  ECMWF 91"
!        print *, '  '
!        call minmax_uv(u,im_e,jm_e,km_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of V-wind  ECMWF"
!        print *, '  '
!        call minmax_uv(v,im_e,jm_e,km_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of pt  ECMWF"
!        print *, '  '
!        call minmax_uv(pt,im_e,jm_e,km_e)

!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of qltot  ECMWF ",km_e
!        print *, '  '
!        call minmax_uv(qltot,im_e,jm_e,km_e)

!        print *,'*** ec2fv: calling gmap ***'
!        print *,' im_e,jm_e,km_e,kn_e ',im_e,jm_e,km_e,kn_e
         call gmap ( im_e,jm_e,nq, kappa,     
     .               km_e, pke1, pe1, u, v, pt, q,
     .               kn_e, pke0, pe0, u0, v0, thv0, q0,
     .               o_m=o, w_m=qitot,wl_m=qltot,o_n=o0 ,
     .               w_n=qitot0, wl_n=qltot0)

!        print *,'*** ec2fv: out from gmap ***'
!        print *, '  '
!        print *, '  '
!        print *," MIN/MAX of qltot0 	ECMWF ",kn_e
!        print *, '  '
!        call minmax_uv(qltot0,im_e,jm_e,kn_e)



#if   (openmp)
!$omp  parallel do
!$omp& default (shared)
!$omp& private (i,j,L,pl,alf)
#endif

       if (present(blend)) then
        print *,' *** Blending is performed *** '
! -----------------------------------------------------------
!        New routine to blend the u,v and theta using windshaer.
!        Ravi 8/1/2006.
! -----------------------------------------------------------
!
          if (present(wsh)) then
           print *,' Calling blend_uvsh *** ',present(wsh)
           call blend_uvsh(u1,v1, u0,v0,thv1,thv0,pe0,pbelow,pabove,im_e,jm_e,kn_e)
          else
           print *,' *** blend_uvsh is not being used *** ', present(wsh)
          endif
!          ----------------------------------------
!              Blend the top .1 and above.
!          ----------------------------------------

            do  L = 1,kn_e
             do  j = 1,jm_e
              do  i = 1,im_e
                pl=0.5*(pe0(i,j,L+1)+pe0(i,j,L))
                alf=(pl-pabove)/(pbelow-pabove)
!               if ( i == 1 .and. j == 1) then
!                write(6,1000) i,j,L,pl,pe0(i,j,L+1),pe0(i,j,L),pbelow,pabove,alf
!               endif
 1000           format(3(2x,i3),5(2x,f10.3),2x,f12.4)
                if( pl.lt.pabove ) then
                 if (present(wsh)) then
                      q0(i,j,L) =   q1(i,j,L)
                      o0(i,j,L) =   o1(i,j,L)
                  qitot0(i,j,L) = qitot1(i,j,L)
                  qltot0(i,j,L) = qltot1(i,j,L)
                 else
                      u0(i,j,L) =   u1(i,j,L)
                      v0(i,j,L) =   v1(i,j,L)
                    thv0(i,j,L) = thv1(i,j,L)
                      q0(i,j,L) =   q1(i,j,L)
                      o0(i,j,L) =   o1(i,j,L)
                  qitot0(i,j,L) = qitot1(i,j,L)
                  qltot0(i,j,L) = qltot1(i,j,L)
                 endif
                else if( pl.lt.pbelow ) then
!                print *,' *** ec2fv:We are here **** '
                 if (present(wsh)) then
                    q0(i,j,L) =   q0(i,j,L)*alf +   q1(i,j,L)*(1-alf)
                    o0(i,j,L) =   o0(i,j,L)*alf +   o1(i,j,L)*(1-alf)
                  qitot0(i,j,L) = qitot0(i,j,L)*alf + qitot1(i,j,L)*(1-alf)
                  qltot0(i,j,L) = qltot0(i,j,L)*alf + qltot1(i,j,L)*(1-alf)
                 else
                    u0(i,j,L) =   u0(i,j,L)*alf +   u1(i,j,L)*(1-alf)
                    v0(i,j,L) =   v0(i,j,L)*alf +   v1(i,j,L)*(1-alf)
                  thv0(i,j,L) = thv0(i,j,L)*alf + thv1(i,j,L)*(1-alf)
                    q0(i,j,L) =   q0(i,j,L)*alf +   q1(i,j,L)*(1-alf)
                    o0(i,j,L) =   o0(i,j,L)*alf +   o1(i,j,L)*(1-alf)
                  qitot0(i,j,L) = qitot0(i,j,L)*alf + qitot1(i,j,L)*(1-alf)
                  qltot0(i,j,L) = qltot0(i,j,L)*alf + qltot1(i,j,L)*(1-alf)
                 endif
                endif
              end do
             end do
            end do

!           print *, '  '
!           print *, '  '
!           print *," MIN/MAX of thv0-wind  	BLENDED"
!           print *, '  '
!           call minmax_uv(thv0,im_e,jm_e,kn_e)
        else
           print *,' **** NO BLENDING IS DONE ****'
        endif

                                                                                                            
!        -------------------------------------------------
!           Convert u and v to 'D' grid.
!        -------------------------------------------------

       if(dgrid) then
        print *,'*** ec2fv: calling atod ***'
        call atod ( u0,u0,im_e,jm_e,kn_e,2 )
        call atod ( v0,v0,im_e,jm_e,kn_e,1 )
       endif

       return
       end subroutine ec2fv_
      end module m_ec2fv
