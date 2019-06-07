program test_ncdr_get
    use kinds, only: i_long
    use nc_diag_read, only: nc_diag_read_init, nc_diag_read_close, &
        nc_diag_read_lookup_dim, nc_diag_read_assert_dim, &
        nc_diag_read_check_dim, nc_diag_read_get_dim, &
        nc_diag_read_check_dim_unlim, nc_diag_read_id_init, &
        nc_diag_read_get_dim_names, ncdr_error
    
    implicit none
    
    integer(i_long) :: ndims, ndims_len
    character(len=:), dimension(:), allocatable :: dim_names
    
    integer(i_long) :: i, tmp_ncdr_id, tmp_ncdr_id_2, ind
    character(len=:), allocatable :: dim_name
    
    !------------------------------------------------------------------
    ! Make sure if we close with ncdr_id via caching, we actually clear
    ! the cache!
    !------------------------------------------------------------------
    call nc_diag_read_init("test.nc", tmp_ncdr_id)
    
    write (*, "(A)") " ** File: test.nc (using cached ncdr_id)"
    
    call nc_diag_read_get_dim_names(ndims, ndims_len, dim_names)
    write (*, "(A, I0, A, I0)") " ** Number of dimensions in test.nc: ", ndims, &
        " | Maximum length of dimension names: ", ndims_len
    print *, "** All dimensions: **"
    print *, dim_names
    
    print *, "** Dimension details: **"
    
    do i = 1, ndims
        dim_name = trim(dim_names(i))
        ind = nc_diag_read_lookup_dim(dim_name)
        ind = nc_diag_read_assert_dim(dim_name)
        if (nc_diag_read_check_dim(dim_name) == .FALSE.) &
            call ncdr_error("Can't find dim with check(), even when it's listed!")
        write (*, "(A, I0, A, L)") "    -> Dimension: " // dim_name // " | Size : ", &
            nc_diag_read_get_dim(dim_name), " | Unlimited? ", &
            nc_diag_read_check_dim_unlim(dim_name)
    end do
    
    if (nc_diag_read_lookup_dim("INVALID_DIM_INVALID") /= -1) &
        call ncdr_error("Invalid dimension lookup result check failed.")
    
    if (nc_diag_read_check_dim("INVALID_DIM_INVALID")) &
        call ncdr_error("Invalid dimension check result check = TRUE failed.")
    
    ! These will result in an error:
    !i = nc_diag_read_assert_dim("INVALID_DIM_INVALID")
    !i = nc_diag_read_get_dim("INVALID_DIM_INVALID")
    !print *, nc_diag_read_check_dim_unlim("INVALID_DIM_INVALID")
    !call nc_diag_read_init("invalid file name.nc/\/\/\")
    
    tmp_ncdr_id_2 = nc_diag_read_id_init("test_fixed.nc")
    
    deallocate(dim_names)
    
    write (*, "(A)") " ** File: test_fixed.nc (using ncdr_id)"
    
    call nc_diag_read_get_dim_names(tmp_ncdr_id_2, ndims, ndims_len, dim_names)
    write (*, "(A, I0, A, I0)") " ** Number of dimensions in test_fixed.nc: ", ndims, &
        " | Maximum length of dimension names: ", ndims_len
    print *, "** All dimensions: **"
    print *, dim_names
    
    print *, "** Dimension details: **"
    
    do i = 1, ndims
        dim_name = trim(dim_names(i))
        ind = nc_diag_read_lookup_dim(tmp_ncdr_id_2, dim_name)
        ind = nc_diag_read_assert_dim(tmp_ncdr_id_2, dim_name)
        if (nc_diag_read_check_dim(tmp_ncdr_id_2, dim_name) == .FALSE.) &
            call ncdr_error("Can't find dim with check(), even when it's listed!")
        write (*, "(A, I0, A, L)") "    -> Dimension: " // dim_name // " | Size : ", &
            nc_diag_read_get_dim(tmp_ncdr_id_2, dim_name), " | Unlimited? ", &
            nc_diag_read_check_dim_unlim(tmp_ncdr_id_2, dim_name)
    end do
    
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id)
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id_2)
end program test_ncdr_get
