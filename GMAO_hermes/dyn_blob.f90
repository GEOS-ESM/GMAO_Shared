program dyn_blob

use m_dyn,  only: dyn_get
use m_dyn,  only: dyn_put
use m_dyn,  only: dyn_vect
use m_dyn,  only: dyn_clean
use m_const, only: radius_earth

implicit none

character(len=*), parameter :: fname = 'bkg.eta.nc4'
character(len=*), parameter :: pname = 'fsens.eta.nc4'
character(len=*), parameter :: bname = 'blob.eta.nc4'
integer, parameter :: dyntype=5
integer nymd, nhms, ier, im, jm, km, freq, nstep
integer ii, nymdw, nhmsw
integer, parameter :: npnts=4
real blocs(2,npnts)
logical :: advect=.false.

real,parameter:: adrate     = 0.75
real,parameter:: lenfcst    = 24
real,parameter:: corrlen    = 800.                      ! to be read from file
real,parameter:: corrlength = corrlen*1000/radius_earth

integer, parameter :: zlevout = 1 ! -1 won''t do anything

type(dyn_vect) :: wind
type(dyn_vect) :: pert

integer :: ntimes = 12
integer :: lu=10
real,allocatable :: covloc(:,:,:,:)
real,allocatable :: adcovloc(:,:,:,:)
real,allocatable :: ploc(:,:,:)

! read in perturbation fields
call dyn_get ( trim(pname), nymd, nhms, pert, ier, timidx=1, freq=freq, nstep=nstep, vectype=dyntype )
im=pert%grid%im
jm=pert%grid%jm
km=pert%grid%km

if (zlevout>0) then
   open (lu,file='tracer.grd',form='unformatted',access='sequential',convert='little_endian')
   call wrtout_ (lu,pert%q(:,:,zlevout,1))
   close(lu)
   call readin_(lu,pert%q(:,:,zlevout,1))
endif

pert%u = 0.0
pert%v = 0.0
pert%ps = 0.0
pert%delp = 0.0
pert%q = 0.0
pert%pt = 0.0

call set_blobs_
call make_blobs_
call dyn_put ( trim(bname), nymd, nhms, 0, pert, ier, freq=freq, nstep=nstep, vectype=dyntype )

! read in wind fields
 call dyn_get ( trim(fname), nymdw, nhmsw, wind, ier, timidx=1, freq=freq, nstep=nstep, vectype=dyntype )
if (advect) then
   call advect_blobs_
endif
if (zlevout>0) then
   open (lu,file='wind.grd',form='unformatted',access='sequential',convert='little_endian')
   call wrtout_ (lu,wind%u(:,:,zlevout))
   call wrtout_ (lu,wind%v(:,:,zlevout))
   close(lu)
endif

call clean_blobs_

if(zlevout<=0) then
  call dyn_put ( trim(pname), nymd, nhms, 0, pert, ier, freq=freq, nstep=nstep, vectype=dyntype )
endif

contains

subroutine wrtout_(lu,fld)
  integer, intent(in) :: lu
  real, intent(in) :: fld(:,:)
  real(4),allocatable:: fld4(:,:)
  integer myim,myjm,ndim
  myim=size(fld,1)
  myjm=size(fld,2)
  ndim = myim*myjm
  allocate(fld4(myim,myjm))
  fld4=fld
  call lon_shift(fld4,myim,myjm)
  print *, 'wrtout: sum field ', myim, myjm, sum(fld)/ndim, sum(fld4)/ndim
  write(lu) fld4
  deallocate(fld4) 
end subroutine wrtout_

subroutine readin_(lu,fld)
  integer, intent(in) :: lu
  real, intent(in) :: fld(:,:)
  real(4),allocatable:: fld4(:,:)
  integer myim,myjm,ndim
  myim=size(fld,1)
  myjm=size(fld,2)
  ndim = myim*myjm
  allocate(fld4(myim,myjm))
  open (lu,file='tracer.grd',form='unformatted',access='sequential',convert='little_endian')
  read(lu) fld4
  close(lu)
  call lon_shift(fld4,myim,myjm)
  print *, 'readin: sum field ', myim, myjm, sum(fld4)/ndim
  deallocate(fld4) 
end subroutine readin_

subroutine set_blobs_
  blocs(1,1)=10  ! :,1=lats, :,2=lons
  blocs(2,1)= 0

! blocs(1,2)=90
  blocs(1,2)=30
  blocs(2,2)=180

  blocs(1,3)=45
  blocs(2,3)=50

! blocs(1,4)=-90
  blocs(1,4)=-30
  blocs(2,4)=-60

  allocate(ploc(1,1,3))
  allocate(covloc(im,jm,3,npnts))
  if(advect) allocate(adcovloc(im,jm,3,npnts))
end subroutine set_blobs_
subroutine clean_blobs_
  if(advect) deallocate(adcovloc)
  deallocate(covloc)
  deallocate(ploc)
end subroutine clean_blobs_

subroutine make_blobs_
  integer nn,mm,ii,jj,iii,jjj
  real pi,cs,sn,dist

  mm=km
  do nn=1,npnts
     call globeloc ( ploc, blocs(1,nn:nn), blocs(2,nn:nn) )
     call globeloc ( covloc(:,:,:,nn), pert%grid%lat, pert%grid%lon )
     do jj=1,jm
        do ii=1,im
           dist = sqrt( (covloc(ii,jj,1,nn)-ploc(1,1,1))**2 + &
                        (covloc(ii,jj,2,nn)-ploc(1,1,2))**2 + &
                        (covloc(ii,jj,3,nn)-ploc(1,1,3))**2   )/corrlength
           if (dist<10*corrlength) pert%pt(ii,jj,:mm) = gc(dist)
        enddo
     enddo
     mm=mm-1
  enddo

end subroutine make_blobs_

subroutine advect_blobs_
  integer nn,ii,jj,kk,iii,jjj
  real fcstlen,dt,pi,cs,sn,dist

end subroutine advect_blobs_


real function gc(r)
 ! Gaspari-Cohn taper function.
 ! r should be positive, and normalized so taper = 0 at r = 1
 ! very close to exp(-(r/c)**2), where c = 0.388
 implicit none

 real a1,a2,a3,a4,a5,a6,a7,a8,a9
 parameter(a1 = -8.0)
 parameter(a2 = 8.0)
 parameter(a3 = 5.0)
 parameter(a4 = 20.0/3.0)
 parameter(a5 = 1.0)
 parameter(a6 = 8.0/3.0)
 parameter(a7 = -10.0)
 parameter(a8 =  4.0)
 parameter(a9 = -1.0/3.0)

 real, intent(in) :: r
 if(r < a5)then
   if(r > 0.5)then
      gc = ( ( ( ( a6*r -a2 )*r +a3 )*r +a4 )*r +a7)*r + a8 + a9/r
   else
      gc = ( ( ( a1*r +a2)*r +a3 )*r -a4)*r*r + a5
   end if
 else
    gc = 0.0
 end if
end function gc

subroutine globeloc (aloc,lat,lon)
 implicit none
 real,intent(inout) :: aloc(:,:,:)
 real,intent(in)    :: lat(:),lon(:)
 real    pi
 integer i,j

 real,allocatable:: clat(:),clon(:),slat(:),slon(:)
 if(size(aloc,3)<3) then
   print *, 'globeloc error: check 2nd dim of aloc ', size(aloc,3)
 endif
 pi=4.0*atan(1.0)
 allocate(clat(size(lat)),clon(size(lon)))
 allocate(slat(size(lat)),slon(size(lon)))
 clat=cos(lat*pi/180.)
 slat=sin(lat*pi/180.)
 clon=cos(lon*pi/180.)
 slon=sin(lon*pi/180.)
 do j=1,size(aloc,2)
   do i=1,size(aloc,1)
      aloc(i,j,1) = clat(j)*clon(i)
      aloc(i,j,2) = clat(j)*slon(i)
      aloc(i,j,3) = slat(j)
   enddo
 enddo
 deallocate(slat,slon)
 deallocate(clat,clon)
end subroutine globeloc

subroutine globeadloc (fcstlen,aloc,lat,lon, u,v)
 implicit none
 real,intent(in)    :: fcstlen
 real,intent(inout) :: aloc(:,:,:)
 real,intent(in)    :: lat(:),lon(:)
 real,intent(in)    :: u(:,:),v(:,:)
 real    pi,pi2,halfpi
 real    adlat,adlon
 real    adtime
 integer i,j

 real,allocatable:: clat(:),clon(:),slat(:),slon(:)
 if(size(aloc,3)<3) then
   print *, 'globeloc error: check 2nd dim of aloc ', size(aloc,3)
 endif
 pi=4.0*atan(1.0)
 pi2=2.0*pi
 halfpi=0.5*pi
 adtime=adrate*fcstlen*3600./radius_earth
 allocate(clat(size(lat)),clon(size(lon)))
 allocate(slat(size(lat)),slon(size(lon)))
 clat=lat*pi/180.
 slat=lat*pi/180.
 clon=lon*pi/180.
 slon=lon*pi/180.
 do j=1,size(aloc,2)
   do i=1,size(aloc,1)
      adlon = clon(i) - u(i,j) * cos(clat(j)) * adtime
      adlat = clat(j) - v(i,j) * adtime
      if(adlat > halfpi) then
         adlat = pi - adlat
         adlon = adlon + pi
      else if(adlat < -halfpi) then
         adlat = -pi - adlat
         adlon = adlon + pi
      end if
      if(adlon > pi2) then
         adlon = mod(adlon,pi2)
      else if(adlon < 0.0) then
         adlon = mod(adlon,pi2) + pi2
      end if

      aloc(i,j,1) = cos(adlat)*cos(adlon)
      aloc(i,j,2) = cos(adlat)*sin(adlon)
      aloc(i,j,3) = sin(adlat)
   enddo
 enddo
 deallocate(slat,slon)
 deallocate(clat,clon)
end subroutine globeadloc

subroutine lon_shift(field,im,jm)
   Implicit NONE

   integer, intent(in) :: im
   integer, intent(in) :: jm

   real(4), intent(inout) :: field(im,jm)
   integer i, j
   real(4) tmp

   do j = 1, jm
      do i = 1, im/2
         tmp = field(i,j)
         field(i,j) = field(i+im/2,j)
         field(i+im/2,j) = tmp
      enddo
   enddo

end subroutine lon_shift

end program dyn_blob
