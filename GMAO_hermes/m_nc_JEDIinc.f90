module m_nc_JEDIinc
use netcdf
implicit none
private

public :: nc_JEDIinc_vars_init
public :: nc_JEDIinc_vars_final
public :: nc_JEDIinc_vars_comp
public :: nc_JEDIinc_vars_copy
public :: nc_JEDIinc_vars
public :: nc_JEDIinc_dims
public :: nc_JEDIinc_read
public :: nc_JEDIinc_write
public :: nc_JEDIinc_summary
public :: nc_JEDIinc_getpointer
public :: nc_JEDIinc_geos2jedi
public :: nc_JEDIinc_jedi2geos

type nc_JEDIinc_vars
   logical :: initialized=.false.
   integer :: nlon,nlat,nsig
   logical :: gsiset=.false.
   real(4),pointer,dimension(:):: ak=>NULL(),bk=>NULL()
   real(4),pointer,dimension(:,:,:):: dp=>NULL()
   real(4),pointer,dimension(:,:,:):: tv=>NULL()
   real(4),pointer,dimension(:,:,:):: t=>NULL()
   real(4),pointer,dimension(:,:,:):: u=>NULL(),v=>NULL()
   real(4),pointer,dimension(:,:,:):: qv=>NULL()
   real(4),pointer,dimension(:,:,:):: qi=>NULL(),ql=>NULL(),qr=>NULL(),qs=>NULL()
   real(4),pointer,dimension(:,:,:):: oz=>NULL()
   real(4),pointer,dimension(:,:)  :: ps=>NULL(),ts=>NULL()
!
   real(4),pointer,dimension(:)    :: v1d=>NULL()
   real(4),pointer,dimension(:,:)  :: v2d=>NULL()
   real(4),pointer,dimension(:,:,:):: v3d=>NULL()
end type nc_JEDIinc_vars

character(len=*), parameter :: myname = 'm_nc_JEDIinc'
real, parameter:: PPMV2GpG = 1.6571E-6 ! from ppmv to g/g
real, parameter:: mbar_per_Pa = 0.01   ! mb to Pa

integer, parameter :: nv2d = 2
character(len=4),parameter :: cvars2d(nv2d) = (/ 'ps  ', 'ts  ' /)

!integer, parameter :: nv3d = 10
!character(len=5),parameter :: cvars3d(nv3d) = (/ &
!                                              'tv   ', 'u    ', 'v    ', &
!                                              'sphu ', 'qitot', 'qltot', &
!                                              'qrtot', 'qstot', 'ozone', &
!                                              'delp '&
!                                              /)
!integer, parameter :: nv2d = 1
!character(len=4),parameter :: cvars2d(nv2d) = (/ 'ps  '/)

integer, parameter :: nv3d = 5
character(len=6),parameter :: cvars3d(nv3d) = (/ &
                                              't     ', 'ua    ', 'va    ', &
                                              'q     ', 'o3ppmv' &
                                              /)

interface nc_JEDIinc_dims; module procedure    &
  read_dims_ ; end interface
!interface nc_JEDIinc_vars; module procedure    &
!  set_vars_ ; end interface
interface nc_JEDIinc_read; module procedure    &
  read_JEDIinc_ ; end interface
interface nc_JEDIinc_write; module procedure    &
  write_JEDIinc_ ; end interface
interface nc_JEDIinc_vars_init; module procedure    &
  init_JEDIinc_vars_ ; end interface
interface nc_JEDIinc_vars_final; module procedure    &
  final_JEDIinc_vars_ ; end interface
interface nc_JEDIinc_vars_comp; module procedure    &
  comp_JEDIinc_vars_ ; end interface
interface nc_JEDIinc_vars_copy; module procedure    &
  copy_ ; end interface
interface nc_JEDIinc_summary; module procedure    &
  summary_ ; end interface
interface nc_JEDIinc_geos2jedi; module procedure    &
  geos2jedi_ ; end interface
interface nc_JEDIinc_jedi2geos; module procedure    &
  jedi2geos_ ; end interface
interface nc_JEDIinc_getpointer
  module procedure get_pointer_2d_
  module procedure get_pointer_3d_
end interface

! internal only
interface stddev_
  module procedure stddev2_
  module procedure stddev3_
end interface
contains

subroutine read_dims_ (fname,nlat,nlon,nlev,rc, myid,root)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  integer, intent(out) :: rc
  integer, intent(out) :: nlat,nlon,nlev
  integer, intent(in), optional :: myid, root

! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid, ier
  integer :: mype_,root_

! Local variables
  character(len=*), parameter :: myname_ = myname//"::dims_"
  logical :: verbose

! Return code (status)
  rc=0; mype_=0; root_=0
  if(present(myid) .and. present(root) ) then
     mype_ = myid
     root_ = root
  endif

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.

  call check_( nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Read global attributes
  call check_( nf90_inq_dimid(ncid, "lon", varid), rc, mype_, root_)
  call check_( nf90_inquire_dimension(ncid, varid, len=nlon), rc, mype_, root_ )
  call check_( nf90_inq_dimid(ncid, "lat", varid), rc, mype_, root_ )
  call check_( nf90_inquire_dimension(ncid, varid, len=nlat), rc, mype_, root_ )
  call check_( nf90_inq_dimid(ncid, "lev", varid), rc, mype_, root_ )
  call check_( nf90_inquire_dimension(ncid, varid, len=nlev), rc, mype_, root_ )

! Close the file, freeing all resources.
  call check_( nf90_close(ncid), rc, mype_, root_ )

  return

end subroutine read_dims_

subroutine read_JEDIinc_ (fname,bvars,rc, myid,root, gsiset)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  type(nc_JEDIinc_vars),intent(inout) :: bvars ! background error variables
  integer, intent(out) :: rc
  integer, intent(in), optional :: myid,root ! accommodate MPI calling programs
  logical, intent(in), optional :: gsiset

! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

! Local variables
  character(len=*), parameter :: myname_ = myname//"::read_"
  character(len=4) :: cindx
  integer :: kk,nv,nl,nlat,nlon,nlev
  integer :: ndims_, nvars_, ngatts_, unlimdimid_
  integer :: nlat_,nlon_,nlev_
  integer :: mype_,root_
  real(4), allocatable :: data_in(:,:,:)
  logical :: verbose
  logical :: init_
  logical :: gsi_

! Return code (status)
  rc=0; mype_=0; root_=0
  verbose=.true.
  init_=.false.
  if(present(myid).and.present(root) )then
    if(myid/=root) verbose=.false.
    mype_ = myid
    root_ = root
  endif

  gsi_=.false.
  if(present(gsiset)) then
     gsi_=gsiset
  endif

! Get dimensions
  call read_dims_ (fname,nlat_,nlon_,nlev_,rc, mype_,root_)

  init_ = bvars%initialized
  if ( init_ ) then
!   Set dims
    nlat=bvars%nlat
    nlon=bvars%nlon
    nlev=bvars%nsig

!   Consistency check
    if (nlon_ /= nlon .or. nlat_ /=nlat .or. nlev_/=nlev ) then
       rc=1
       if(myid==root) then
         print *, 'nlat(file) = ', nlat_, 'nlat(required) = ', nlat
         print *, 'nlon(file) = ', nlon_, 'nlon(required) = ', nlon
         print *, 'nlev(file) = ', nlev_, 'nlev(required) = ', nlev
         print *, myname_,  'Inconsistent dimensions, aborting ... '
       endif
       return
    endif
  else
!   Set dims
    nlat=nlat_
    nlon=nlon_
    nlev=nlev_
    call init_JEDIinc_vars_(bvars,nlon,nlat,nlev,gsi=gsi_)
  endif

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.

  call check_( nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Read global attributes
! call check_( nf90_inquire(ncid, ndims_, nvars_, ngatts_, unlimdimid_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lon", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlon_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lat", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlat_), rc, mype_, root_ )
! call check_( nf90_inq_dimid(ncid, "lev", varid), rc, mype_, root_ )
! call check_( nf90_inquire_dimension(ncid, varid, len=nlev_), rc, mype_, root_ )

! Read data to file
  allocate(data_in(nlon,nlat,1))
  do nv = 1, nv2d
     call check_( nf90_inq_varid(ncid, trim(cvars2d(nv)), varid), rc, mype_, root_ )
     call check_( nf90_get_var(ncid, varid, data_in(:,:,1)), rc, mype_, root_ )
     if (bvars%gsiset) then
       if(trim(cvars2d(nv))=="ps" ) bvars%ps = transpose(data_in(:,:,1))
       if(trim(cvars2d(nv))=="ts" ) bvars%ts = transpose(data_in(:,:,1))
     else
       if(trim(cvars2d(nv))=="ps" ) bvars%ps = data_in(:,:,1)
       if(trim(cvars2d(nv))=="ts" ) bvars%ts = data_in(:,:,1)
     endif
  enddo
  deallocate(data_in)
!
  allocate(data_in(nlon,nlat,nlev))
  do nv = 1, nv3d
     call check_( nf90_inq_varid(ncid, trim(cvars3d(nv)), varid), rc, mype_, root_  )
     call check_( nf90_get_var(ncid, varid, data_in(:,:,:)), rc, mype_, root_ )

     if (bvars%gsiset) then
        if(trim(cvars3d(nv))=="delp") then
           do kk=1,bvars%nsig
              bvars%dp(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="tv"  ) then
           do kk=1,bvars%nsig
              bvars%tv(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="t"  ) then
           do kk=1,bvars%nsig
              bvars%t(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="u" .or. trim(cvars3d(nv))=="ua") then
           do kk=1,bvars%nsig
              bvars%u(:,:,kk)  = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="v" .or. trim(cvars3d(nv))=="va") then
           do kk=1,bvars%nsig
              bvars%v(:,:,kk)  = transpose(data_in(:,:,kk))
           enddo
        endif
!
        if(trim(cvars3d(nv))=="sphu" .or. trim(cvars3d(nv))=="q") then
           do kk=1,bvars%nsig
              bvars%qv(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="qitot" .or. trim(cvars3d(nv))=="qi") then
           do kk=1,bvars%nsig
              bvars%qi(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="qltot" .or. trim(cvars3d(nv))=="ql") then
           do kk=1,bvars%nsig
              bvars%ql(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="qrtot" .or. trim(cvars3d(nv))=="qr") then
           do kk=1,bvars%nsig
              bvars%qr(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
        if(trim(cvars3d(nv))=="qstot" .or. trim(cvars3d(nv))=="qs") then
           do kk=1,bvars%nsig
              bvars%qs(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
!
        if(trim(cvars3d(nv))=="ozone" .or. trim(cvars3d(nv))=="o3ppmv") then
           do kk=1,bvars%nsig
              bvars%oz(:,:,kk) = transpose(data_in(:,:,kk))
           enddo
        endif
     else
        if(trim(cvars3d(nv))=="delp") bvars%dp = data_in(:,:,:)
        if(trim(cvars3d(nv))=="tv"  ) bvars%tv = data_in(:,:,:)
        if(trim(cvars3d(nv))=="t"   ) bvars%t  = data_in(:,:,:)
        if(trim(cvars3d(nv))=="u" .or. trim(cvars3d(nv))=="ua" ) bvars%u  = data_in(:,:,:)
        if(trim(cvars3d(nv))=="v" .or. trim(cvars3d(nv))=="va" ) bvars%v  = data_in(:,:,:)
!
        if(trim(cvars3d(nv))=="sphu"  .or. trim(cvars3d(nv))=="q" ) bvars%qv = data_in(:,:,:)
        if(trim(cvars3d(nv))=="qitot" .or. trim(cvars3d(nv))=="qi") bvars%qi = data_in(:,:,:)
        if(trim(cvars3d(nv))=="qltot" .or. trim(cvars3d(nv))=="ql") bvars%ql = data_in(:,:,:)
        if(trim(cvars3d(nv))=="qrtot" .or. trim(cvars3d(nv))=="qr") bvars%qr = data_in(:,:,:)
        if(trim(cvars3d(nv))=="qstot" .or. trim(cvars3d(nv))=="qs") bvars%qs = data_in(:,:,:)
!
        if(trim(cvars3d(nv))=="ozone" .or. trim(cvars3d(nv))=="o3ppmv") bvars%oz = data_in(:,:,:)
     endif
!
  enddo
  deallocate(data_in)

! Close the file, freeing all resources.
  call check_( nf90_close(ncid), rc, mype_, root_ )

  if(verbose) print *,"*** Finish reading file: ", trim(fname)

! Convert to GEOS units and orientation
  call jedi2geos_(bvars)

  return

end subroutine read_JEDIinc_

subroutine write_JEDIinc_ (fname,bvars,lats,lons,rc, myid,root,plevs)
  implicit none
  character(len=*), intent(in)    :: fname ! input filename
  type(nc_JEDIinc_vars),intent(in)    :: bvars ! background error variables
  real(4), intent(in) :: lats(:)           ! latitudes per GSI: increase index from South to North Pole
  real(4), intent(in) :: lons(:)           ! longitude per GSI: increase index from East to West
  integer, intent(out) :: rc
  real(4), intent(in), optional :: plevs(:)
  integer, intent(in), optional :: myid,root        ! accommodate MPI calling programs

  character(len=*), parameter :: myname_ = myname//"::read_"
  integer, parameter :: NDIMS = 3

! When we create netCDF files, variables and dimensions, we get back
! an ID for each one.
  character(len=4) :: cindx
  integer :: ncid, dimids(NDIMS)
  integer :: x_dimid, y_dimid, z_dimid
  integer :: lon_varid, lat_varid, lev_varid
  integer :: ii,jj,nl,nv,nn,nlat,nlon,nlev
  integer :: mype_,root_
  integer, allocatable :: varid2d(:), varid3d(:)
  logical :: verbose

! This is the data array we will write. It will just be filled with
! a progression of integers for this example.
  real(4), allocatable :: data_out(:,:,:)
  real(4), allocatable :: idlevs(:)

! Consistency check
  if (bvars%gsiset) then
     print *,myname,'write must be in GEOS orientation'
     rc=99
     return
  endif

! Convert to JEDI units and orientation
  call geos2jedi_(bvars)

! Return code (status)
  rc=0; mype_=0; root_=0
  verbose=.true.
  if(present(myid).and.present(root) )then
    if(myid/=root) verbose=.false.
    mype_ = myid
    root_ = root
  endif

! Set dims
  nlat=bvars%nlat
  nlon=bvars%nlon
  nlev=bvars%nsig

! Always check the return code of every netCDF function call. In
! this example program, wrapping netCDF calls with "call check()"
! makes sure that any return which is not equal to nf90_noerr (0)
! will print a netCDF error message and exit.

! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
  call check_( nf90_create(fname, NF90_CLOBBER, ncid), rc, mype_, root_ )
  if(rc/=0) return

! Define the dimensions. NetCDF will hand back an ID for each.
  call check_( nf90_def_dim(ncid, "lon", nlon, x_dimid), rc, mype_, root_ )
  call check_( nf90_def_dim(ncid, "lat", nlat, y_dimid), rc, mype_, root_ )
  call check_( nf90_def_dim(ncid, "lev", nlev, z_dimid), rc, mype_, root_ )

  call check_( nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, lon_varid), rc, mype_, root_ )
  call check_( nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, lat_varid), rc, mype_, root_ )
  call check_( nf90_def_var(ncid, "lev", NF90_REAL, z_dimid, lev_varid), rc, mype_, root_ )

  call check_( nf90_put_att(ncid, lon_varid, "units", "degress"), rc, mype_, root_ )
  call check_( nf90_put_att(ncid, lat_varid, "units", "degress"), rc, mype_, root_ )
  call check_( nf90_put_att(ncid, lev_varid, "units", "Pa"), rc, mype_, root_ )

! The dimids array is used to pass the IDs of the dimensions of
! the variables. Note that in fortran arrays are stored in
! column-major format.
  dimids =  (/ x_dimid, y_dimid, z_dimid /)

! Define variables.
  allocate(varid2d(nv2d))
  do nv = 1, nv2d
     call check_( nf90_def_var(ncid, trim(cvars2d(nv)), NF90_REAL, (/ x_dimid, y_dimid /), varid2d(nv)), rc, mype_, root_ )
  enddo
  allocate(varid3d(nv3d))
  do nv = 1, nv3d
     call check_( nf90_def_var(ncid, trim(cvars3d(nv)), NF90_REAL, (/ x_dimid, y_dimid, z_dimid /), varid3d(nv)), rc, mype_, root_ )
  enddo

! End define mode. This tells netCDF we are done defining metadata.
  call check_( nf90_enddef(ncid), rc, mype_, root_ )

! Write coordinate variables data
  call check_( nf90_put_var(ncid, lon_varid, lons ), rc, mype_, root_ )
  call check_( nf90_put_var(ncid, lat_varid, lats ), rc, mype_, root_ )
  if(present(plevs)) then
    call check_( nf90_put_var(ncid, lev_varid, plevs), rc, mype_, root_ )
  else
    allocate(idlevs(nlev))
    do ii = 1,nlev
       idlevs(ii) = ii
    enddo
    call check_( nf90_put_var(ncid, lev_varid, idlevs), rc, mype_, root_ )
    deallocate(idlevs)
  endif

! Write data to file
  allocate(data_out(nlon,nlat,1))
  do nv = 1, nv2d
     if(trim(cvars2d(nv))=="ps" ) data_out(:,:,1) = bvars%ps
     if(trim(cvars2d(nv))=="ts" ) data_out(:,:,1) = bvars%ts
     call check_( nf90_put_var(ncid, varid2d(nv), data_out(:,:,1)), rc, mype_, root_)
  enddo
  deallocate(data_out)
  allocate(data_out(nlon,nlat,nlev))
  do nv = 1, nv3d
     if(trim(cvars3d(nv))=="delp") data_out(:,:,:) = bvars%dp
     if(trim(cvars3d(nv))=="tv"  ) data_out(:,:,:) = bvars%tv
     if(trim(cvars3d(nv))=="u"   ) data_out(:,:,:) = bvars%u
     if(trim(cvars3d(nv))=="v"   ) data_out(:,:,:) = bvars%v
!
     if(trim(cvars2d(nv))=="sphu" ) data_out(:,:,:) = bvars%qv
     if(trim(cvars2d(nv))=="qitot") data_out(:,:,:) = bvars%qi
     if(trim(cvars2d(nv))=="qltot") data_out(:,:,:) = bvars%ql
     if(trim(cvars2d(nv))=="qrtot") data_out(:,:,:) = bvars%qr
     if(trim(cvars2d(nv))=="qstot") data_out(:,:,:) = bvars%qs
!
     if(trim(cvars2d(nv))=="ozone") data_out(:,:,:) = bvars%oz
!
     call check_( nf90_put_var(ncid, varid3d(nv), data_out(:,:,:)), rc, mype_, root_ )
  enddo
  deallocate(data_out)

! Close file
  call check_( nf90_close(ncid), rc, mype_, root_ )

  deallocate(varid3d)
  deallocate(varid2d)

  print *, "*** Finish writing file ", fname

  return

end subroutine write_JEDIinc_

subroutine init_JEDIinc_vars_(vr,nlon,nlat,nsig,gsi)

  integer,intent(in) :: nlon,nlat,nsig
  type(nc_JEDIinc_vars) vr
  logical,intent(in),optional :: gsi

  if(vr%initialized) return

  if(present(gsi)) then
    vr%gsiset = gsi
  endif

  vr%nlon=nlon
  vr%nlat=nlat
  vr%nsig=nsig

! allocate single precision arrays
  if (vr%gsiset) then
     print *, 'GSI-like 2d-array'
     allocate(vr%tv(nlat,nlon,nsig),vr%u (nlat,nlon,nsig),vr%v (nlat,nlon,nsig),vr%qv(nlat,nlon,nsig),&
              vr%qi(nlat,nlon,nsig),vr%ql(nlat,nlon,nsig),vr%qr(nlat,nlon,nsig),vr%qs(nlat,nlon,nsig),&
              vr%oz(nlat,nlon,nsig),vr%dp(nlat,nlon,nsig),vr%t(nlon,nlat,nsig) )
     allocate(vr%ps(nlat,nlon),vr%ts(nlat,nlon))
  else
     print *, 'GEOS-like 2d-array'
     allocate(vr%tv(nlon,nlat,nsig),vr%u (nlon,nlat,nsig),vr%v (nlon,nlat,nsig),vr%qv(nlon,nlat,nsig),&
              vr%qi(nlon,nlat,nsig),vr%ql(nlon,nlat,nsig),vr%qr(nlon,nlat,nsig),vr%qs(nlon,nlat,nsig),&
              vr%oz(nlon,nlat,nsig),vr%dp(nlon,nlat,nsig),vr%t(nlon,nlat,nsig) )
     allocate(vr%ps(nlon,nlat),vr%ts(nlon,nlat))
  endif
  vr%initialized=.true.
  end subroutine init_JEDIinc_vars_

  subroutine final_JEDIinc_vars_(vr)
  type(nc_JEDIinc_vars) vr
! deallocate arrays
  if(.not. vr%initialized) return
  deallocate(vr%tv,vr%u ,vr%v ,vr%qv,  &
             vr%qi,vr%ql,vr%qr,vr%qs,&
             vr%oz,vr%dp,vr%t)
  deallocate(vr%ps,vr%ts)
  vr%initialized=.false.
end subroutine final_JEDIinc_vars_

subroutine comp_JEDIinc_vars_(va,vb,rc, myid,root)
  type(nc_JEDIinc_vars) va
  type(nc_JEDIinc_vars) vb
  integer, intent(out) :: rc
  integer, intent(in), optional :: myid,root        ! accommodate MPI calling programs
  character(len=*), parameter :: myname_ = myname//"::comp_JEDIinc_vars_"
  integer :: ii,jj
  logical :: verbose, failed
  real :: tolerance = 10.e-10
  integer, allocatable :: ier(:)
!
  rc=0
  verbose=.true.
  if(present(myid).and.present(root) )then
    if(myid/=root) verbose=.false.
  endif
! Consistency check
  if (va%nlon/=vb%nlon .or. va%nlat/=vb%nlat .or. va%nsig/=vb%nsig ) then
     rc=1
     if(myid==root) then
       print *, 'nlat(va) = ', va%nlat, 'nlat(vb) = ', vb%nlat
       print *, 'nlon(va) = ', va%nlon, 'nlon(vb) = ', vb%nlon
       print *, 'nlev(va) = ', va%nsig, 'nlev(vb) = ', vb%nsig
       print *, myname_,  'Inconsistent dimensions, aborting ... '
     endif
     return
  endif

  allocate(ier(nv2d+nv3d))
  ii=0;ier=0
  ii=ii+1; if(abs(sum(va%dp - vb%dp)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%tv - vb%tv)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%t  - vb%t )) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%u  - vb%u )) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%v  - vb%v )) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%qv - vb%qv)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%qi - vb%qi)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%ql - vb%ql)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%qr - vb%qr)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%qs - vb%qs)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%oz - vb%oz)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%ps - vb%ps)) >tolerance) ier(ii)=ii
  ii=ii+1; if(abs(sum(va%ts - vb%ts)) >tolerance) ier(ii)=ii
  failed=.false.
  do jj=1,ii
     if(ier(jj)/=0.and.verbose) then
       print *, 'Found field ', jj, ' not to match'
       failed=.true.
     endif
  enddo
  deallocate(ier)
  if (.not.failed) then
       if(verbose) print *, 'Comp finds all fields to match'
  endif
end subroutine comp_JEDIinc_vars_

subroutine copy_(ivars,ovars,rc)
  type(nc_JEDIinc_vars) ivars
  type(nc_JEDIinc_vars) ovars
  integer, intent(out) :: rc
  integer :: kk

  rc=0
  if (ovars%nlon/=ivars%nlon .or. &
      ovars%nlat/=ivars%nlat .or. &
      ovars%nsig/=ivars%nsig ) then
      print*, 'copy_JEDIinc_vars_: Trying to copy inconsistent vectors, aborting ...'
      rc=99
      return
  endif

  if (ivars%gsiset .neqv. ovars%gsiset ) then
     do kk=1,ovars%nsig
        ovars%dp(:,:,kk) = transpose(ivars%dp(:,:,kk))
        ovars%tv(:,:,kk) = transpose(ivars%tv(:,:,kk))
        ovars%u (:,:,kk) = transpose(ivars%u (:,:,kk))
        ovars%v (:,:,kk) = transpose(ivars%v (:,:,kk))
        ovars%qv(:,:,kk) = transpose(ivars%qv(:,:,kk))
        ovars%qi(:,:,kk) = transpose(ivars%qi(:,:,kk))
        ovars%ql(:,:,kk) = transpose(ivars%ql(:,:,kk))
        ovars%qr(:,:,kk) = transpose(ivars%qr(:,:,kk))
        ovars%qs(:,:,kk) = transpose(ivars%qs(:,:,kk))
        ovars%oz(:,:,kk) = transpose(ivars%oz(:,:,kk))
     enddo
     ovars%ps = transpose(ivars%ps)
     ovars%ts = transpose(ivars%ts)
  else
     ovars%dp = ivars%dp
     ovars%tv = ivars%tv
     ovars%u  = ivars%u
     ovars%v  = ivars%v
     ovars%qv = ivars%qv
     ovars%qi = ivars%qi
     ovars%ql = ivars%ql
     ovars%qr = ivars%qr
     ovars%qs = ivars%qs
     ovars%oz = ivars%oz

     ovars%ps = ivars%ps
     ovars%ts = ivars%ts
  endif

end subroutine copy_

subroutine get_pointer_2d_ (vname, bvars, ptr, rc )
implicit none
character(len=*), intent(in) :: vname
type(nc_JEDIinc_vars) bvars
real(4),pointer,intent(inout) :: ptr(:,:)
integer,intent(out) :: rc
rc=-1
if(trim(vname)=='ps') then
  ptr => bvars%ps
  rc=0
endif
if(trim(vname)=='ts') then
  ptr => bvars%ts
  rc=0
endif
end subroutine get_pointer_2d_

subroutine get_pointer_3d_ (vname, bvars, ptr, rc )
implicit none
character(len=*), intent(in) :: vname
type(nc_JEDIinc_vars) bvars
real(4),pointer,intent(inout) :: ptr(:,:,:)
integer,intent(out) :: rc
character(len=5) :: var
rc=-1
!
var='delp'
if(trim(vname)==trim(var)) then
  ptr => bvars%dp
  rc=0
  return
endif
!
var='tv'
if(trim(vname)==trim(var)) then
  ptr => bvars%tv
  rc=0
  return
endif
!
if(trim(vname)=='u' .or. trim(vname)=='ua') then
  ptr => bvars%u
  rc=0
  return
endif
!
if(trim(vname)=='v' .or. trim(vname)=='va') then
  ptr => bvars%v
  rc=0
  return
endif
!
var='sphu'
if(trim(vname)=='sphu' .or. trim(vname)=='q') then
  ptr => bvars%qv
  rc=0
  return
endif
!
var='qitot'
if(trim(vname)==trim(var)) then
  ptr => bvars%qi
  rc=0
  return
endif
!
var='qltot'
if(trim(vname)==trim(var)) then
  ptr => bvars%ql
  rc=0
  return
endif
!
var='qrtot'
if(trim(vname)==trim(var)) then
  ptr => bvars%qr
  rc=0
  return
endif
!
var='qstot'
if(trim(vname)==trim(var)) then
  ptr => bvars%qs
  rc=0
  return
endif
!
var='ozone'
if(trim(vname)==trim(var)) then
  ptr => bvars%oz
  rc=0
  return
endif
end subroutine get_pointer_3d_

subroutine check_(status,rc, myid, root)
    integer, intent ( in) :: status
    integer, intent (out) :: rc
    integer, intent ( in) :: myid, root
    rc=0
    if(status /= nf90_noerr) then
      if(myid==root) print *, trim(nf90_strerror(status))
      rc=999
    end if
end subroutine check_

subroutine geos2jedi_ (x)
  implicit none
  type(nc_JEDIinc_vars) x

  !==> input ps in Pa - convert to hPa(mb)
! x%grid%ak  = x%grid%ak * mbar_per_Pa
! x%ps       = x%ps * mbar_per_Pa
! x%dp       = x%dp * mbar_per_Pa
! x%oz       = x%oz * PPMV2GpG
  ! need flip so localization function applies equaly to EnKF and Hybrid-GSI
  call flip_(x)
  ! still need a transpose
end subroutine geos2jedi_

subroutine jedi2geos_ (x)
  implicit none
  type(nc_JEDIinc_vars) x

  !==> input ps in mbar - convert to Pa
! x%grid%ak   = x%grid%ak    / mbar_per_Pa
! x%ps = x%ps / mbar_per_Pa
! x%dp = x%dp / mbar_per_Pa
! x%oz = x%oz / PPMV2GpG
  ! need flip: horiz. orientation in JEDI is [0,360]
  call flip_(x)
end subroutine jedi2geos_

subroutine flip_(x)
  implicit none
  type(nc_JEDIinc_vars) x
  integer im,jm,km
  im=x%nlon
  jm=x%nlat
  km=x%nsig
!
  call hflip2_(x%ps,im,jm,x%gsiset)
  call hflip2_(x%ts,im,jm,x%gsiset)
!
  call hflip3_(x%dp,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%dp,im,jm,km)

  call hflip3_(x%tv,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%tv,im,jm,km)

  call hflip3_(x%t,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%t,im,jm,km)

  call hflip3_(x%u ,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%u ,im,jm,km)

  call hflip3_(x%v ,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%v ,im,jm,km)

  call hflip3_(x%qv,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%qv,im,jm,km)

  call hflip3_(x%qi,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%qi,im,jm,km)

  call hflip3_(x%ql,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%ql,im,jm,km)

  call hflip3_(x%qr,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%qr,im,jm,km)

  call hflip3_(x%qs,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%qs,im,jm,km)

  call hflip3_(x%oz,im,jm,km,x%gsiset)
  if(x%gsiset) call vflip_ (x%oz,im,jm,km)

end subroutine flip_

subroutine hflip3_ ( q,im,jm,km, gsi )
    implicit none
    integer  im,jm,km,i,j,k
    logical  gsi
    real(4), intent(inout) :: q(:,:,:)
    real(4), allocatable   :: dum(:)
    allocate ( dum(im) )
    if (gsi) then
       do k=1,km
       do j=1,jm
       do i=1,im/2
          dum(i) = q(j,i+im/2,k)
          dum(i+im/2) = q(j,i,k)
       enddo
          q(j,:,k) = dum(:)
       enddo
       enddo
    else
       do k=1,km
       do j=1,jm
       do i=1,im/2
          dum(i) = q(i+im/2,j,k)
          dum(i+im/2) = q(i,j,k)
       enddo
          q(:,j,k) = dum(:)
       enddo
       enddo
    endif
    deallocate ( dum )
end subroutine hflip3_

subroutine hflip2_ ( q,im,jm, gsi )
    implicit none
    integer  im,jm,i,j
    logical  gsi
    real(4), intent(inout) :: q(:,:)
    real(4), allocatable   :: dum(:)
    allocate ( dum(im) )
    if (gsi) then
       do j=1,jm
       do i=1,im/2
          dum(i) = q(j,i+im/2)
          dum(i+im/2) = q(j,i)
       enddo
          q(j,:) = dum(:)
       enddo
    else
       do j=1,jm
       do i=1,im/2
          dum(i) = q(i+im/2,j)
          dum(i+im/2) = q(i,j)
       enddo
          q(:,j) = dum(:)
       enddo
    endif
    deallocate ( dum )
end subroutine hflip2_

subroutine vflip_(q,im,jm,km)
   implicit none
   integer,intent(in) :: im,jm,km
   real(4),intent(inout) :: q(im,jm,km)
   real(4), allocatable  :: dum(:)
   integer i,j
   allocate(dum(km))
   do j=1,jm
      do i=1,im
         dum      = q(i,j,:)
         q(i,j,:) = dum(km:1:-1)
     end do
  end do
  deallocate(dum)
end subroutine vflip_

subroutine summary_(x,myid)
implicit none
  type(nc_JEDIinc_vars) x
  integer, intent(in), optional :: myid
  integer lu,nxy,nxyz
  nxy  = x%nlat*x%nlon
  nxyz = nxy*x%nsig
  lu=6
  if(present(myid)) lu=myid
  write(lu,'(a)') "================================================"
  write(lu,'(a)') "var     min        max       mean      stddev   "
  write(lu,'(a)') "================================================"
  write(lu,'(a4,1p,4(e10.3,1x))') 'ps', minval(x%ps), maxval(x%ps), sum(x%ps)/nxy,  stddev_(x%ps)
  write(lu,'(a4,1p,4(e10.3,1x))') 'ts', minval(x%ts), maxval(x%ts), sum(x%ts)/nxy,  stddev_(x%ts)
  write(lu,'(a4,1p,4(e10.3,1x))') 'dp', minval(x%dp), maxval(x%dp), sum(x%dp)/nxyz, stddev_(x%dp)
  write(lu,'(a4,1p,4(e10.3,1x))') 'tv', minval(x%tv), maxval(x%tv), sum(x%tv)/nxyz, stddev_(x%tv)
  write(lu,'(a4,1p,4(e10.3,1x))') 'qv', minval(x%qv), maxval(x%qv), sum(x%qv)/nxyz, stddev_(x%qv)
  write(lu,'(a4,1p,4(e10.3,1x))') 'qi', minval(x%qi), maxval(x%qi), sum(x%qi)/nxyz, stddev_(x%qi)
  write(lu,'(a4,1p,4(e10.3,1x))') 'ql', minval(x%ql), maxval(x%ql), sum(x%ql)/nxyz, stddev_(x%ql)
  write(lu,'(a4,1p,4(e10.3,1x))') 'qr', minval(x%qr), maxval(x%qr), sum(x%qr)/nxyz, stddev_(x%qr)
  write(lu,'(a4,1p,4(e10.3,1x))') 'qs', minval(x%qs), maxval(x%qs), sum(x%qs)/nxyz, stddev_(x%qs)
  write(lu,'(a4,1p,4(e10.3,1x))') 'oz', minval(x%oz), maxval(x%oz), sum(x%oz)/nxyz, stddev_(x%oz)
  write(lu,'(a)') "================================================"
end subroutine summary_

real function stddev2_(x)
 implicit none
 real(4) :: x(:,:)
 integer im,jm
 real mean
 im = size(x,1)
 jm = size(x,2)
 mean = sum(x)/(im*jm)
 stddev2_= sqrt(sum((x-mean)*(x-mean)))/(im*jm-1)
end function stddev2_

real function stddev3_(x)
 implicit none
 real(4) :: x(:,:,:)
 integer im,jm,km
 real mean
 im = size(x,1)
 jm = size(x,2)
 km = size(x,3)
 mean = sum(x)/(im*jm*km)
 stddev3_= sqrt(sum((x-mean)*(x-mean)))/(im*jm*km-1)
end function stddev3_

end module m_nc_JEDIinc
