!
      subroutine blend_uvsh(u0,v0, u1,v1,theta0,theta1, pe1,  &
                           pbelow,pabove,im,jm,lm)
      implicit none

      integer    i,j,L,im,jm,lm
      character(len=*), parameter :: myname = 'uvsh_incr'

      real       u1(im,jm,lm),   v1(im,jm,lm)
      real       u0(im,jm,lm),   v0(im,jm,lm)
      real       uz(im,jm,lm),   vz(im,jm,lm)
      real       theta0(im,jm,lm),   theta1(im,jm,lm)
      real       pe1(im,jm,lm+1)
      real       pbelow,pabove,pl

      real,allocatable ::  du(:,:,:)
      real,allocatable ::  dv(:,:,:)
      real,allocatable ::  duz(:,:,:)
      real,allocatable ::  dvz(:,:,:)
      real,allocatable ::  dpt(:,:,:)
      real,allocatable ::  gduz(:,:,:)
      real,allocatable ::  gdvz(:,:,:)
      real,allocatable ::  nduz(:,:,:)
      real,allocatable ::  ndvz(:,:,:)
      real*4,allocatable ::  dummyu(:,:)

      print *,'  ',trim(myname), ':lm ',lm

      
      allocate (dummyu(im,jm))
      allocate (duz(im,jm,lm+1))
      allocate (dvz(im,jm,lm+1))
      allocate (gduz(im,jm,lm+1))
      allocate (gdvz(im,jm,lm+1))
      allocate (nduz(im,jm,lm+1))
      allocate (ndvz(im,jm,lm+1))
      allocate (du(im,jm,lm+1))
      allocate (dv(im,jm,lm+1))
      allocate (dpt(im,jm,lm+1))
!     open(88,file='uwnd.sheer.data2',form='unformatted', &
!          access='sequential',status='unknown')

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of U-wind  GMAO"
!      print *, '  '
!      call minmax_uv(u0,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of V-wind  GMAO"
!      print *, '  '
!      call minmax_uv(v0,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of Theta  GMAO"
!      print *, '  '
!      call minmax_uv(theta0,im,jm,lm)

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of U-wind  NCEP"
!      print *, '  '
!      call minmax_uv(u1,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of V-wind  NCEP"
!      print *, '  '
!      call minmax_uv(v1,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of Theta  NCEP"
!      print *, '  '
!      call minmax_uv(theta1,im,jm,lm)

      do L=2,lm
       do j=1,jm
        do i=1,im

 
!           --------------------------------------------------
!            Step-1
!            Compute the wind sheer for GMAO and the NCEP u,v.
!           --------------------------------------------------

           gduz(i,j,L) =   u0(i,j,L) - u0(i,j,L-1)
           gdvz(i,j,L) =   v0(i,j,L) - v0(i,j,L-1)

           nduz(i,j,L) =   u1(i,j,L) - u1(i,j,L-1)
           ndvz(i,j,L) =   v1(i,j,L) - v1(i,j,L-1)

         end do
        end do
       end do
       gduz(:,:,1) = u0(:,:,1)
       nduz(:,:,1) = u1(:,:,1)

      do L=1,lm
       do j=1,jm
        do i=1,im
           dpt(i,j,L) =   theta1(i,j,L) - theta0(i,j,L)
         end do
        end do
       end do

!      do  L = lm,1,-1
!       dummyu(:,:) = u0(:,:,L) 
!       write(88) dummyu
!      end do

!      do  L = lm,1,-1
!       dummyu(:,:) = u1(:,:,L) 
!       write(88) dummyu
!      end do

!      do  L = lm,1,-1
!       dummyu(:,:) = gduz(:,:,L) 
!       write(88) dummyu
!      end do

!      do  L = lm,1,-1
!       dummyu(:,:) = nduz(:,:,L) 
!       write(88) dummyu
!      end do
        

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of U-wind increments between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(du,im,jm,lm)
!      print *," MIN/MAX of V-wind increments between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(dv,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of Theta increments between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(dpt,im,jm,lm)
 
      do L=2,lm
       do j=1,jm
        do i=1,im

!           --------------------------------------
!            Step-3
!            Compute u,v sheere increment.
!           --------------------------------------

          duz(i,j,L) =   nduz(i,j,L) - gduz(i,j,L)
          dvz(i,j,L) =   ndvz(i,j,L) - gdvz(i,j,L)

          du(i,j,L) =   nduz(i,j,L) - gduz(i,j,L)
          dv(i,j,L) =   ndvz(i,j,L) - gdvz(i,j,L)

        end do
       end do
      end do

!      do  L = lm,1,-1
!       dummyu(:,:) = duz(:,:,L) 
!       write(88) dummyu
!      end do

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX U-wind increments shear between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(duz,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of V-wind increments shear between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(dvz,im,jm,lm)
 
      do L=2,lm
       do j=1,jm
        do i=1,im
          pl=0.5*(pe1(i,j,L+1)+pe1(i,j,L))

!           --------------------------------------
!            Step-4
!            Tamperring the increments duz,dvz, and dpt.

!            reproduce du,dv and dpt from the the shear.
!           --------------------------------------

           if(pl <= pabove) then
            dpt(i,j,L) = 0.0
           elseif(pabove < pl .and. pl < pbelow) then
            dpt(i,j,L) = dpt(i,j,L) * (pl-pabove)/(pbelow-pabove)
           endif

           if(pe1(i,j,L) <= pabove) then
            duz(i,j,L) = 0.0
            dvz(i,j,L) = 0.0
            du(i,j,L) = 0.0
            dv(i,j,L) = 0.0
           elseif(pe1(i,j,L) > pabove .and. pe1(i,j,L) < pbelow) then
!           print *,'pabove, pe1(',i,j,L,'), ','pbelow ',pabove,pe1(i,j,L),pbelow
!           duz(i,j,L) =   duz(i,j,L) * (pl-pabove)/(pbelow-pabove)
!           dvz(i,j,L) =   dvz(i,j,L) * (pl-pabove)/(pbelow-pabove)

            duz(i,j,L) =   duz(i,j,L) * (pe1(i,j,L)-pabove)/(pbelow-pabove)
            dvz(i,j,L) =   dvz(i,j,L) * (pe1(i,j,L)-pabove)/(pbelow-pabove)
           endif
        end do
       end do
      end do

!    -------------------------------------------------------
!      Step-4a
!          Reconstruct the wind increment sarting at pbelow.
!    -------------------------------------------------------

      do L=lm,1,-1
       do j=1,jm
        do i=1,im
!          print *,' L,du,dz ',L,du(i,j,L+1),duz(i,j,L+1)
           du(i,j,L+1) = du(i,j,L+1) + duz(i,j,L+1)
           dv(i,j,L+1) = dv(i,j,L+1) + dvz(i,j,L+1)
        end do
       end do
      end do

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced U-wind increments from shear "
!      print *, '  '
!      call minmax_uv(duz,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced V-wind increments from shear "
!      print *, '  '
!      call minmax_uv(dvz,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced Theta increments between NCEP and GMAO"
!      print *, '  '
!      call minmax_uv(dpt,im,jm,lm)
!

!    -------------------------------------------------------
!      Step-5
!          Reconstruct the blended field.
!    -------------------------------------------------------

      do L=1,lm
       do j=1,jm
        do i=1,im
           pl=0.5*(pe1(i,j,L+1)+pe1(i,j,L))
           if(pl < pbelow) then
            u1(i,j,L) = u0(i,j,L) + du(i,j,L)
            v1(i,j,L) = v0(i,j,L) + dv(i,j,L)
            theta1(i,j,L) = theta0(i,j,L) + dpt(i,j,L)
           endif
        end do
       end do
      end do

!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced U-wind from shear "
!      print *, '  '
!      call minmax_uv(u1,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced V-wind  from shear "
!      print *, '  '
!      call minmax_uv(v1,im,jm,lm)
!      print *, '  '
!      print *, '  '
!      print *," MIN/MAX of reproduced Theta "
!      print *, '  '
!      call minmax_uv(theta1,im,jm,lm)
!      print *, '  '

      deallocate(dummyu)
      deallocate(du,dv,dpt,duz,dvz)
      deallocate(gduz,gdvz,nduz,ndvz)

      return
    end subroutine blend_uvsh

      subroutine minmax_uv(x,im,jm,lm)
      integer im,jm,km,imj,l
      real x(im,jm,lm),xmin,xmax

!         -------------------------------------------------
!            Compute min max for a given field.
!         -------------------------------------------------

      print *,'  '
      print *,'  '
      print *,'        Layer       MIN                    MAX '
      print *,'  '
      print *,'  '
      do  l = lm,1,-1
       xmin = 10.e+20
       xmax = -10.e+20
       do j = 1,jm
        do i = 1,im
         if(x(i,j,l) < xmin) xmin = x(i,j,l)
         if(x(i,j,l) > xmax) xmax = x(i,j,l)
        end do
       end do
       print *,'   ',l,xmin,xmax
      end do
     end subroutine minmax_uv
