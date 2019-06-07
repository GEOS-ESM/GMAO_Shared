program test_ncdr_vars
    use kinds, only: i_long
    
    use nc_diag_read, only: nc_diag_read_init, nc_diag_read_close, &
        nc_diag_read_lookup_var, nc_diag_read_assert_var, &
        nc_diag_read_check_var, nc_diag_read_get_var_ndims, &
        nc_diag_read_ret_var_dims, nc_diag_read_get_type_str, &
        nc_diag_read_get_var_type, nc_diag_read_id_init, &
        nc_diag_read_get_var_names, ncdr_error
    
    implicit none
    
    integer(i_long) :: nvars, nvars_len
    character(len=:), dimension(:), allocatable :: var_names
    
    integer(i_long) :: i, j, tmp_ncdr_id, tmp_ncdr_id_2, ind
    character(len=:), allocatable :: var_name
    
    integer(i_long), dimension(:), allocatable  :: var_dims
    
    call nc_diag_read_init("test.nc", tmp_ncdr_id)
    
    write (*, "(A)") " ** File: test.nc (using cached ncdr_id)"
    
    call nc_diag_read_get_var_names(nvars, nvars_len, var_names)
    write (*, "(A, I0, A, I0)") " ** Number of variables in test.nc: ", nvars, &
        " | Maximum length of variable names: ", nvars_len
    print *, "** All variables: **"
    print *, var_names
    
    print *, "** Variable details: **"
    
    do i = 1, nvars
        var_name = trim(var_names(i))
        ind = nc_diag_read_lookup_var(var_name)
        ind = nc_diag_read_assert_var(var_name)
        if (nc_diag_read_check_var(var_name) == .FALSE.) &
            call ncdr_error("Can't find var with check(), even when it's listed!")
        write (*, "(A, I0, A)", advance = "no") &
            "    -> Variable: " // var_name // " | Number of dimensions : ", &
            nc_diag_read_get_var_ndims(var_name), " | Type " // &
            nc_diag_read_get_type_str(nc_diag_read_get_var_type(var_name)) // &
            " | Dimensions: ("
        
        allocate(var_dims(nc_diag_read_get_var_ndims(var_name)))
        var_dims = nc_diag_read_ret_var_dims(var_name)
        
        do j = 1, nc_diag_read_get_var_ndims(var_name)
            if (j > 1) write (*, "(A)", advance = "no") ", "
            write (*, "(I0)", advance = "no") var_dims(j)
        end do
        
        write (*, "(A)") ")"
        deallocate(var_dims)
    end do
    
    if (nc_diag_read_lookup_var("INVALID_VAR_INVALID") /= -1) &
        call ncdr_error("Invalid variable lookup result check failed.")
    
    if (nc_diag_read_check_var("INVALID_VAR_INVALID")) &
        call ncdr_error("Invalid variable check result check = TRUE failed.")
    
    ! These will result in an error:
    !i = nc_diag_read_assert_var("INVALID_VAR_INVALID")
    !i = nc_diag_read_get_var_ndims("INVALID_VAR_INVALID")
    !i = nc_diag_read_get_var_type("INVALID_VAR_INVALID")
    !call nc_diag_read_get_var_dims("INVALID_VAR_INVALID")
    !print *, nc_diag_read_ret_var_dims("INVALID_VAR_INVALID")
    !call nc_diag_read_init("invalid file name.nc/\/\/\")
    
    tmp_ncdr_id_2 = nc_diag_read_id_init("test_fixed.nc")
    
    deallocate(var_names)
    
    write (*, "(A)") " ** File: test_fixed.nc (using ncdr_id)"
    
    call nc_diag_read_get_var_names(tmp_ncdr_id_2, nvars, nvars_len, var_names)
    write (*, "(A, I0, A, I0)") " ** Number of variables in test_fixed.nc: ", nvars, &
        " | Maximum length of variable names: ", nvars_len
    print *, "** All variables: **"
    print *, var_names
    
    print *, "** Variable details: **"
    
    do i = 1, nvars
        var_name = trim(var_names(i))
        ind = nc_diag_read_lookup_var(tmp_ncdr_id_2, var_name)
        ind = nc_diag_read_assert_var(tmp_ncdr_id_2, var_name)
        if (nc_diag_read_check_var(tmp_ncdr_id_2, var_name) == .FALSE.) &
            call ncdr_error("Can't find var with check(), even when it's listed!")
        write (*, "(A, I0, A)", advance = "no") &
            "    -> Variable: " // var_name // " | Number of dimensions : ", &
            nc_diag_read_get_var_ndims(tmp_ncdr_id_2, var_name), " | Type " // &
            nc_diag_read_get_type_str(nc_diag_read_get_var_type(tmp_ncdr_id_2, var_name)) // &
            " | Dimensions: ("
        
        allocate(var_dims(nc_diag_read_get_var_ndims(tmp_ncdr_id_2, var_name)))
        var_dims = nc_diag_read_ret_var_dims(tmp_ncdr_id_2, var_name)
        
        do j = 1, nc_diag_read_get_var_ndims(tmp_ncdr_id_2, var_name)
            if (j > 1) write (*, "(A)", advance = "no") ", "
            write (*, "(I0)", advance = "no") var_dims(j)
        end do
        
        write (*, "(A)") ")"
        deallocate(var_dims)
    end do
    
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id)
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id_2)
end program test_ncdr_vars
