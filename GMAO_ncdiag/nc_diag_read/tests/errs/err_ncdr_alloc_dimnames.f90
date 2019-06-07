program err_ncdr_alloc_dimnames
    use kinds, only: i_long
    use nc_diag_read, only: nc_diag_read_init, &
        nc_diag_read_get_dim_names
    
    implicit none
    
    integer(i_long) :: ndims, ndims_len
    character(len=:), dimension(:), allocatable :: dim_names
    
    call nc_diag_read_init("test.nc")
    
    allocate(character(len=10) :: dim_names(10))
    call nc_diag_read_get_dim_names(ndims, ndims_len, dim_names)
end program err_ncdr_alloc_dimnames
