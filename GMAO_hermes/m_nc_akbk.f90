module m_nc_akbk
use netcdf
implicit none
private
public :: write_nc_akbk
interface write_nc_akbk
  module procedure write_nc_akbk_
end interface
integer :: mype_ = 0 
integer :: root_ = 0 
contains
subroutine write_nc_akbk_ (fname,ak, bk)

  character(len=*), intent(in) :: fname
  double precision, intent(in) :: ak(:), bk(:)

  ! NetCDF variables
  integer :: ncid, dimid_edge
  integer :: varid_edge, varid_ak, varid_bk
  integer :: rc
  integer :: i, km
  double precision, allocatable :: edge(:)
  
  km = size(ak)
  allocate(edge(km))
  do i=1,km
    edge(i) = dble(i)
  enddo
 
  ! Create the NetCDF file
  call check_( nf90_create(fname, NF90_CLOBBER, ncid), rc, mype_, root_ )
  if(rc/=0) return

  ! Define dimensions
  call check_ ( nf90_def_dim(ncid, "edge", km, dimid_edge), rc, mype_, root_ )
  if(rc/=0) return

  ! Define variables
  call check_ ( nf90_def_var(ncid, "edge", NF90_DOUBLE, (/dimid_edge/), varid_edge), rc, mype_, root_ )
  call check_ ( nf90_def_var(ncid, "ak",   NF90_DOUBLE, (/dimid_edge/), varid_ak), rc, mype_, root_ )
  call check_ ( nf90_def_var(ncid, "bk",   NF90_DOUBLE, (/dimid_edge/), varid_bk), rc, mype_, root_ )

  ! Add attributes to edge
  call check_ ( nf90_put_att(ncid, varid_edge, "units", "level"), rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_edge, "long_name", "sigma at layer edges"), rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_edge, "standard_name", "atmosphere_hybrid_sigma_pressure_coordinate"),&
                             rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_edge, "coordinate", "eta"), rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_edge, "positive", "down"), rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_edge, "formulaTerms", "ap: ak b: bk ps: ps p0: p00"),&
                             rc, mype_, root_ )

  ! Add attributes to ak
  call check_ ( nf90_put_att(ncid, varid_ak, "long_name", "hybrid_sigma_pressure_a"),&
                rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_ak, "units", "Pa"), rc, mype_, root_ )

  ! Add attributes to bk
  call check_ ( nf90_put_att(ncid, varid_bk, "long_name", "hybrid_sigma_pressure_b"), &
                rc, mype_, root_ )
  call check_ ( nf90_put_att(ncid, varid_bk, "units", "1"), rc, mype_, root_ )

  ! Global attributes
  call check_ ( nf90_put_att(ncid, NF90_GLOBAL, "NASA/GMAO", &
       "Homepage = http://gmao.gfsc.nasa.gov/"), rc, mype_, root_ )

  ! End define mode
  call check_ ( nf90_enddef(ncid), rc, mype_, root_ )

  ! Write data
  call check_ ( nf90_put_var(ncid, varid_edge, edge), rc, mype_, root_ )
  call check_ ( nf90_put_var(ncid, varid_ak, ak), rc, mype_, root_ )
  call check_ ( nf90_put_var(ncid, varid_bk, bk), rc, mype_, root_ )

  ! Close file
  call check_ ( nf90_close(ncid), rc, mype_, root_ )

  deallocate(edge)

  print *, "NetCDF file ",trim(fname)," written successfully."

end subroutine write_nc_akbk_
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
end module m_nc_akbk

