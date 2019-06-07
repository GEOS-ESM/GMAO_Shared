program  main

  use m_dyn, only: dyn_vect, dyn_get, dyn_put, dyn_init
  implicit none
  integer,parameter :: IPREC=0
  integer :: nstep, rc, ier
  integer :: nymd,nhms,freq
  integer :: n,iargc,nargs
  character*255, allocatable :: arg(:)
  character(len=255) :: infile,outfile
  type(dyn_vect) :: fv ! fv dynamics vector in eta
  type(dyn_vect) :: fv5 ! fv dynamics vector in eta
!
  nargs =  iargc()
  if( nargs.eq.0 ) then
     print *,'usage: g4tog5.x inetafile outetafile'
     stop
  end if
  allocate ( arg(nargs) )
  do n=1,nargs
     call getarg(n,arg(n))
  enddo
  infile = trim(arg(1))
  outfile = trim(arg(2))
  nymd=0;nhms=0
  print *, 'dyn_get ',trim(infile)
  call dyn_get(infile, nymd, nhms, fv, freq=freq, nstep=nstep, vectype=4, rc=ier )
  print *,'transform fields g4 to g5'
  call dyn_g4tog5 ( fv )
  call dyn_init ( fv%grid%im, fv%grid%jm, fv%grid%km, fv%grid%lm, fv5, rc,  &
                  fv%grid%ptop, fv%grid%ks, fv%grid%ak, fv%grid%bk, vectype=5 )
  fv5%phis    = fv%phis  
  fv5%hs_stdv = fv%hs_stdv  
  fv5%ts      = fv%ts  
  fv5%lwi     = fv%lwi  
  fv5%ps      = fv%ps  
  fv5%delp = fv%delp  
  fv5%u    = fv%u  
  fv5%v    = fv%v  
  fv5%pt   = fv%pt 
  fv5%q    = fv%q 

! hack fractions
  fv5%frland=fv%lwi
  
  print *,'dyn_put ',trim(outfile)
  call dyn_put(outfile,nymd,nhms,IPREC,fv5,ier,freq=freq,nstep=nstep,vectype=5,forceflip=.true.)

end program  main

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  dyn_g4tog5 --- creates a particular geos5 export state
!
! !INTERFACE:
!
  subroutine  dyn_g4tog5 ( w_f )
!
! !USES:
!
  use m_dyn
  use m_set_eta
  use m_remap, only: d2a => remap_d2a

  implicit NONE

!
! !INPUT/OUTPUT PARAMETERS:
!
  type(dyn_vect), intent(inout)        :: w_f    ! dynamics state vector

!
! !DESCRIPTION: Creates the GSI dynamics state vector {\tt w\_f} 
!               as produced by GEOS5. 
!
! !REVISION HISTORY:
!
!  07Aug2006 Cruz  Initial code.
!  17Apr2007 Todling Removed vertical flip and scaling.
!  01Jul2010 Ravi/daSilva Fix to pole-winds: use d2a instead of dtoa
!
!EOP
!-------------------------------------------------------------------------

   integer :: im, jm, km, lm,l, ks, km1
!  real, parameter :: kappa = 287.04/1004.64  ! r/cp
!  real, parameter :: zvir  =  0.61           ! rh2o/rair - 1.
   real, parameter :: grav  =  9.806          ! accel. of gravity
   real, parameter :: ppmv2du =  1.657e-6     ! ppmv to Dobson units
   real, parameter :: kPa_per_Pa=.001         ! converts Pa to kPa
   real, allocatable :: tmpu(:,:,:)
   real, allocatable :: ak(:), bk(:)
   real, allocatable :: swap(:,:,:)
   real, allocatable :: ua(:,:), va(:,:), coslon(:), sinlon(:)
   real              :: ptop, pint,dl
   integer           :: i,k

! start

!  Short hand for dimensions
!  -------------------------
   im = w_f%grid%im; jm = w_f%grid%jm
   km = w_f%grid%km; lm = w_f%grid%lm

   allocate(ak(km+1))
   allocate(bk(km+1))
!
!    Allocated working arrays for d2a
!
   allocate ( ua(im,jm) )
   allocate ( va(im,jm) )
   allocate ( coslon(im) )
   allocate ( sinlon(im) )

   ak=w_f%grid%ak
   bk=w_f%grid%bk
   call set_eta ( km, ks, ptop, pint, ak, bk )
   print *,'km1, ks, ptop, pint: ', km, ks, ptop, pint

!  print *,' dtoa(u)...'
!  call dtoa ( w_f%u,w_f%u,im,jm,km,2 )
!  print *,' dtoa(v)...'
!  call dtoa ( w_f%v,w_f%v,im,jm,km,1 )

! -------------------------------------------------------
!   Changes made to correct the pole.
!   Used the original d2a routine.
!
   dl = 8.*atan(1.0) / float(im)
   do  i = 1,im/2
     coslon(i)       = -cos((i-1)*dl)
     coslon(i+im/2) = -coslon(i)
     sinlon(i)       = -sin((i-1)*dl)
     sinlon(i+im/2) = -sinlon(i)
   enddo
   
   do k=1,km
    call d2a(w_f%u(:,:,k), w_f%v(:,:,k), ua,va, &
             im, jm, 1, jm, coslon, sinlon)
    w_f%u(:,:,k) = ua
    w_f%v(:,:,k) = va
   end do
! -------------------------------------------------------
 
   print *,' convert PT to Tv...'
   allocate(tmpu(im,jm,km))
! using SJ Lin's algorithm (see diag2dyn.F)
   call tmpu2pt ( tmpu, w_f%q(:,:,:,1), w_f%delp, im, jm, km, ks, &
        w_f%grid%ptop, w_f%pt )
   w_f%pt = tmpu

! calcutale total cloud fraction and place in 3rd slot of q
!  if(lm == 4) then
!     print *, 'replace 3rd slot of q w/ total cloud fraction'
!     w_f%q(:,:,:,3) = w_f%q(:,:,:,3) + w_f%q(:,:,:,4)
!     w_f%grid%lm = 3 ! force to avoid writing more than needed
!  endif
  
! print *,' swap vertical levels'
! allocate ( swap(im,jm,km) )
! swap = w_f%u
! w_f%u(:,:,km:1:-1) = swap(:,:,1:km:+1)
! swap = w_f%v
! w_f%v(:,:,km:1:-1) = swap(:,:,1:km:+1)
! swap = w_f%delp
! w_f%delp(:,:,km:1:-1) = swap(:,:,1:km:+1)
! swap = w_f%pt
! w_f%pt(:,:,km:1:-1) = swap(:,:,1:km:+1)
! do L=1,lm
!    swap = w_f%q(:,:,:,L)
!    w_f%q(:,:,km:1:-1,L) = swap(:,:,1:km:+1)
! enddo
! deallocate(swap)
! w_f%grid%ak(1:km+1:+1)=w_f%grid%ak(km+1:1:-1)
! w_f%grid%bk(1:km+1:+1)=w_f%grid%bk(km+1:1:-1)

! print *,' other adjustments...'
! w_f%grid%ak = w_f%grid%ak*kPa_per_Pa        ! from Pa to kPa

! w_f%ps = (w_f%ps/1000.0)                    ! ps to cb
! w_f%phis = w_f%phis/grav                    ! geopot to m
! w_f%q(:,:,:,2) = ppmv2du * w_f%q(:,:,:,2)   ! ozone in D.U.

   deallocate(ua,va,coslon,sinlon)
!  All done
!  --------
   return

  CONTAINS


!..................................................................................

      subroutine tmpu2pt ( tmpu, sphu, delp, im, jm, km, ks, ptop, pt )
!
!     Computes Scaled Potential Temperature (p**cappa) * T based on an
!     algorithm by SJ Lin (from m_insitu.F).
!
      USE m_const, only: cp    => cpm
      use m_const, only: rair  => rgas
      use m_const, only: cappa => kappa
      use m_const, only: zvir

      implicit none

      integer, intent(in) :: im, jm, km, ks
      real, intent(in) :: sphu(im,jm,km)
      real, intent(in) :: delp(im,jm,km)
      real, intent(in) :: pt(im,jm,km)
      real, intent(in) :: ptop

      real, intent(out) :: tmpu(im,jm,km)

!                      -------------

      integer :: i, j, k
      real   pkz(im,jm,km)
      real   pk(im,jm,km+1)            ! pressure at edges (hPa)
      real   pe(im,km+1,jm)            ! log(pe)
      real peln(im,km+1,jm)            ! log(pe)


!     Compute pk, pe
!     --------------
      call geopm ( ptop, pe, pk, delp, im, jm, km, 1, jm,     &
                  cp, cappa )


!     Compute pkz for conversion between pt and temperature
!     -----------------------------------------------------
      call pkez ( im, jm, km, 1, jm, ptop, &
                 pe, pk, cappa, ks, peln, pkz, .false.)

      do k=1,km
         do j=1,jm
            do i=1,im
!               pt(i,j,k) = tmpu(i,j,k) /pkz(i,j,k)/(1.+zvir*sphu(i,j,k))
                tmpu(i,j,k) = pt(i,j,k)*pkz(i,j,k)
            end do
         end do
      end do

      return

      end subroutine tmpu2pt

!..........................................................................

      subroutine geopm(ptop,pe,pk,delp,im,jm, km,jfirst, jlast, cp,akap )

! From SJ, slightly modified

      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, jfirst, jlast
      real akap, cp, ptop
      real delp(im,jm,km)

! !OUTPUT PARAMETERS
      real pe(im,km+1,jm)                ! only if id .ne. 0
      real pk(im,jm,km+1)                !

! Local:
      real pk2(im,km+1)
      integer i, j, k
      real p2d(im,km+1)
      real ptk


      do 1000 j=jfirst,jlast

        ptk  = ptop ** akap

        do i=1,im
          p2d(i,1) = ptop
          pk2(i,1) = ptk
        enddo

!c Top down
        do k=2,km+1
          do i=1,im
            p2d(i,k)  = p2d(i,k-1) + delp(i,j,k-1)
            pk2(i,k) = p2d(i,k) ** akap
          enddo
        enddo

!c Bottom up

          do k=1,km+1
            do i=1,im
              pe(i,k,j) = p2d(i,k)
              pk(i,j,k) = pk2(i,k)
            enddo
          enddo

1000  continue

      return
      end subroutine geopm

!..........................................................................

      subroutine pkez(im, jm, km, jfirst, jlast, ptop,&
                    pe, pk, akap,  ks, peln, pkz, eta)
!C
!C eta: true on eta coordinate; pk needs to be updated

!C true:
!C Input:  pe
!C Output: pk, pkz, peln

!C false:
!C Input:  pk, pe
!C Output: pkz

! WS 99.05.19 : Added im, jm, km as arguments
! WS 99.07.27 : Limited region to jfirst:jlast

      implicit none

! WS 99.05.19 : Removed fvcore.h

      integer im, jm, km, jfirst, jlast
      real  pe(im, km+1, jm)
      real  pk(im, jm, km+1)
      real  pkz(im, jm, km)
      real peln(im, km+1, jm)
      real ptop

      integer ks
      logical eta

      real akap

!C Local
      real  pk2(im, km+1)
      real pek
      real lnp

      integer i, j, k, j1, jmm0

      j1   = max(1,jfirst)
      jmm0 = min(jm,jlast)


! WS 99.07.27 : Limited region to jfirst:jlast

!!!   do 1000 j=1, jm
      do 1000 j=j1, jmm0

      if ( eta ) then

!C <<<<<<<<<<< Eta cordinate Coordinate  >>>>>>>>>>>>>>>>>>>

      pek =   ptop ** akap
      lnp = log(pe(1,1,j))

      do i=1,im
          pk2(i,1)   = pek
         peln(i,1,j) = lnp
      enddo

      if(ks .ne. 0) then
      do k=2, ks+1
             pek = pe(1,k,j)**akap
             lnp = log(pe(1,k,j))
         do i=1,im
             pk2(i,k)   = pek
            peln(i,k,j) =  lnp
         enddo
      enddo

      do k=1, ks
           pek = (pk2(1,k+1)- pk2(1,k))/(akap*(peln(1,k+1,j) - peln(1,k,j)) )
           do i=1,im
              pkz(i,j,k) = pek
           enddo
      enddo

      endif

      do k=ks+2,km
         do i=1,im
            pk2(i,k) = pe(i,k,j)**akap
         enddo
      enddo

      do i=1,im
         pk2(i,km+1) = pk(i,j,km+1)
      enddo

      do k=ks+2,km+1
         do i=1,im
            peln(i,k,j) =  log(pe(i,k,j))
         enddo
      enddo

      do k=ks+1,km
         do i=1,im
            pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k))/(akap*(peln(i,k+1,j) - peln(i,k,j)) )
         enddo
      enddo

      do k=2,km
         do i=1,im
            pk(i,j,k) = pk2(i,k)
         enddo
      enddo

      else

!C <<<<<<<<<<< General Coordinate  >>>>>>>>>>>>>>>>>>>

      pek =   ptop ** akap
      lnp = log(pe(1,1,j))

      do i=1,im
          pk2(i,1) = pek
         peln(i,1,j) = lnp
      enddo

      do k=2,km+1
         do i=1,im
            peln(i,k,j) =  log(pe(i,k,j))
             pk2(i,k) =  pk(i,j,k)
         enddo
      enddo

        do k=1,km
           do i=1,im
              pkz(i,j,k) = (       pk2(i,k+1) - pk2(i,k) )  / &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
           enddo
        enddo

      endif
1000  continue

      return
      end subroutine pkez
! -------------------------------------------------

  end Subroutine dyn_g4tog5
