!******************************************************************************
  program surf_bkg
!******************************************************************************
!
! surf_bkg
!
! description: reads/regrids fv surface data to generate a
!              surface file for the SSI
!
! input :
!   needs a FV-generated file that contains surface data
!   needs the FV dynamics vector to retrieve u,v winds
!
! output :
!   surface data file used by the SSI
!
! history:
!    12/2002  C. Cruz created
!
!-------------------------------------------------------------------------

   use m_dyn 
   use m_ioutil, only : luavail, opnieee, clsieee

   implicit none

! part I

   character (len = 132) :: myname   ! routine name
   real undef

   integer iid, & ! file id for input surface file
           oid, & ! file id for output surface file
           rc, err
   integer i,j,k,l

! part II

   real, allocatable, dimension(:)   :: dlam,dphi
   real, allocatable, dimension(:,:) :: lons,lats
   integer :: idate(4), glat2
   real :: hourg
   !integer, parameter :: glon = 192, glat = 94
   !integer, parameter :: im = 192, jm = 94
   !integer, parameter :: glon = 384, glat = 190
   !integer, parameter :: im = 384, jm = 190
   integer, parameter :: glon = 512, glat = 254
   integer, parameter :: im = 512, jm = 254
   real*4, allocatable, dimension (:,:) :: fld4, ifld4
   real(8), allocatable, dimension (:,:) :: fld8, ifld8
   integer, allocatable, dimension (:,:) :: ifld, iifld
   real(4) fhour4

! start

   myname='ssi_surf'

   write (*,*)
   write (*,*) "###############################################"
   write (*,*) "          SSI SURFACE Application "
   write (*,*) "###############################################"
   write(*,'(/)')

! try to read an NCEP surface file in the same way that the SSI does it
! (see read_guess.F90)
! The SSI reads the surface fields with size glonxglat2
! but expects them to be of size glonxglat !!
! Therefore the fields are written to sfcanl (see below) with size glonxglat

   !open(30,file='/scratch1/ccruz/nov/sfcanl',form='unformatted')
   !open(31,file='sfcanl',form='unformatted')

   !open(32,file='ncepsfc.62',form='unformatted')

   !allocate(ifld(im,jm))
   allocate(fld4(im,jm))
   allocate(fld8(im,jm))
   allocate (dlam(im))
   allocate (dphi(jm))
   glat2 = glat-2
   !allocate(iifld(glon,glat2))
   allocate(ifld8(glon,glat2))
   allocate(ifld4(glon,glat2))
   allocate (lons(glon,glat2))
   allocate (lats(glon,glat2))

   !j=0
   !read(30)  !; read(31)
   !j=j+1
   !read(30)fhour4,idate
   !write(*,*) 'read : ',fhour4, idate
   !read(31)fhour4,idate
   !write(*,*) 'read : ',fhour4, idate
   !j=j+1
   !j=0
   !do k=1,100
   !read(30, end=333)fld4
   !write(*,*)'NCEP : ',k,minval(fld4),maxval(fld4)
   !read(31, end=333)fld4
   !write(*,*)'FVSSI : ',k,minval(fld4),maxval(fld4)
   !if(k==3.or.k==12.or.k==15.or.k==16) then
   ! write(32)fld4
   !  write(*,*)' -- ',minval(fld4),maxval(fld4)
   !end if
   !j=j+1
   !end do

!  call comput_gaus(im,jm,glon,glat2,dphi,dlam,lons,lats,.false.)  ! orig code
   call comput_gaus(im,jm,glon,glat2,dphi,dlam,lons,lats,.true.)   ! what I think should be RT

   open(30,file='sfcanl',form='unformatted')
   open(31,file='ssisfc.dat',form='unformatted')
   open(34,file='ncepsfc',form='unformatted')

   read(30)
   read(30)fhour4,idate
   write(*,*) 'read : ',fhour4, idate
   write(32)
   write(32) fhour4, idate

   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)

   ifld4=ifld8 ; write(32) ifld4 !tskin
  
   read(30)fld4 ; write(31) fld4 !

   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !soil moisture
  
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !snow cover
  
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !tground
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?

   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?

   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !sli mask
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !veg frac
   write(34)fld4
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !?
 
   read(30)fld4 ; write(31) fld4 !
   fld8=fld4
   if(im/=glon) call interp_h (fld8,im,jm,1,dlam,dphi,0.0,90.0,0.0,ifld8, &
   glon*glat2,lons,lats,1,3,.true.,undef)
   ifld4=ifld8 ; write(32) ifld4 !10mwfrac
 
   read(30)fld4 ; write(31) fld4 !
   write(32) fld4 !veg type
   write(34)fld4
 
   read(30)fld4 ; write(31) fld4 !
   write(32) fld4 !soil type
   write(34)fld4

222 continue 

   !j=0
   !do k=1,100
   !read(30, end=333)fld4
   !j=j+1
   !write(32)fld4
   !end do

!333 write(*,*) 'read n= ',j

   close(34)

   open(34,file='ncepsfc',form='unformatted')
   read(34)fld4; print *,minval(fld4),maxval(fld4)
   read(34)fld4; print *,minval(fld4),maxval(fld4)
   read(34)fld4; print *,minval(fld4),maxval(fld4)
   close(34)

   end program 
