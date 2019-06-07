program err_ncdr_alloc_varnames
    use kinds, only: i_long
    use nc_diag_read, only: nc_diag_read_init, &
        nc_diag_read_get_dim_names, nc_diag_read_get_var_names
    
    implicit none
    
    integer(i_long) :: nvars, nvars_len
    character(len=:), dimension(:), allocatable :: var_names
    
    integer(i_long) :: ndims, ndims_len
    character(len=:), dimension(:), allocatable :: dim_names
    
    integer(i_long) :: i, var_type, var_ndims
    
    call nc_diag_read_init("test.nc")
    
    ! Using an invalid NCDR ID should fail:
    !print *, nc_diag_read_get_dim(1234567, "asdf")
    
    call nc_diag_read_get_dim_names(ndims, ndims_len, dim_names)
    write (*, "(A, I0, A, I0)") " ** Number of dimensions: ", ndims, &
        " | Maximum length of dimension names: ", ndims_len
    print *, "** All dimensions: **"
    print *, dim_names
    
    allocate(character(len=10) :: var_names(10))
    
    call nc_diag_read_get_var_names(nvars, nvars_len, var_names)
end program err_ncdr_alloc_varnames
