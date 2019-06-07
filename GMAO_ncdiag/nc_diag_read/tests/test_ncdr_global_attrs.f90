program test_ncdr_attrs
    use kinds, only: i_byte, i_short, i_long, r_single, r_double
    
    use nc_diag_read, only: nc_diag_read_init, nc_diag_read_close, &
        nc_diag_read_check_global_attr, &
        nc_diag_read_get_global_attr_type, nc_diag_read_id_init, &
        nc_diag_read_get_global_attr_len, &
        nc_diag_read_ret_global_attr_len, &
        nc_diag_read_get_global_attr_names, &
        nc_diag_read_assert_global_attr, &
        nc_diag_read_get_global_attr, &
        nc_diag_read_noid_get_global_attr_1d_string, &
        nc_diag_read_id_get_global_attr_1d_string, &
        ncdr_error
    use netcdf, only: NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, &
        NF90_DOUBLE, NF90_CHAR, NF90_FILL_BYTE, NF90_FILL_SHORT, &
        NF90_FILL_INT, NF90_FILL_FLOAT, NF90_FILL_DOUBLE, &
        NF90_FILL_CHAR
    
    implicit none
    
    integer(i_long) :: nattrs, nattrs_len
    character(len=:), dimension(:), allocatable :: attr_names
    
    integer(i_long) :: i, tmp_ncdr_id, tmp_ncdr_id_2, ind, attr_len, attr_type
    character(len=:), allocatable :: attr_name
    
    call nc_diag_read_init("test.nc", tmp_ncdr_id)
    
    write (*, "(A)") " ** File: test.nc (using cached ncdr_id)"
    
    call nc_diag_read_get_global_attr_names(nattrs, nattrs_len, attr_names)
    write (*, "(A, I0, A, I0)") " ** Number of attributes in test.nc: ", nattrs, &
        " | Maximum length of attribute names: ", nattrs_len
    print *, "** All Attributes: **"
    print *, attr_names
    
    print *, "** Attribute details: **"
    
    do i = 1, nattrs
        attr_name = trim(attr_names(i))
        call nc_diag_read_assert_global_attr(attr_name, attr_type, attr_len)
        if (nc_diag_read_check_global_attr(attr_name) == .FALSE.) &
            call ncdr_error("Can't find attr with check(), even when it's listed!")
        call nc_diag_read_get_global_attr_len(attr_name, attr_len)
        write (*, "(A, I0)") "    -> Attribute: " // attr_name // " | Length: ", &
            attr_len
        
        if (attr_len /= nc_diag_read_ret_global_attr_len(attr_name)) &
            call ncdr_error("nc_diag_read_get_global_attr_len != nc_diag_read_ret_global_attr_len!")
        
        if (attr_type /= nc_diag_read_get_global_attr_type(attr_name)) &
            call ncdr_error("nc_diag_read_assert_global_attr != nc_diag_read_get_global_attr_type!")
        
        if (attr_type == NF90_BYTE) call display_1d_attr_byte(attr_name)
        if (attr_type == NF90_SHORT) call display_1d_attr_short(attr_name)
        if (attr_type == NF90_INT) call display_1d_attr_long(attr_name)
        if (attr_type == NF90_FLOAT) call display_1d_attr_float(attr_name)
        if (attr_type == NF90_DOUBLE) call display_1d_attr_double(attr_name)
        if (attr_type == NF90_CHAR) call display_1d_attr_string(attr_name)
    end do
    
    if (nc_diag_read_check_global_attr("INVALID_attr_INVALID")) &
        call ncdr_error("Invalid attribute check result check = TRUE failed.")
    
    ! These will result in an error:
    !i = nc_diag_read_assert_attr("INVALID_attr_INVALID")
    !i = nc_diag_read_get_global_attr("INVALID_attr_INVALID")
    !print *, nc_diag_read_check_attr_unlim("INVALID_attr_INVALID")
    !call nc_diag_read_init("invalid file name.nc/\/\/\")
    
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id)
    
    tmp_ncdr_id = nc_diag_read_id_init("test_fixed.nc")
    
    deallocate(attr_names)
    
    write (*, "(A)") " ** File: test_fixed.nc (using ncdr_id)"
    
    call nc_diag_read_get_global_attr_names(tmp_ncdr_id, nattrs, nattrs_len, attr_names)
    write (*, "(A, I0, A, I0)") " ** Number of attributes in test_fixed.nc: ", nattrs, &
        " | Maximum length of attribute names: ", nattrs_len
    print *, "** All Attributes: **"
    print *, attr_names
    
    print *, "** Attribute details: **"
    
    do i = 1, nattrs
        attr_name = trim(attr_names(i))
        call nc_diag_read_assert_global_attr(tmp_ncdr_id, attr_name, attr_type, attr_len)
        if (nc_diag_read_check_global_attr(tmp_ncdr_id, attr_name) == .FALSE.) &
            call ncdr_error("Can't find attr with check(), even when it's listed!")
        call nc_diag_read_get_global_attr_len(tmp_ncdr_id, attr_name, attr_len)
        
        if (attr_len /= nc_diag_read_ret_global_attr_len(tmp_ncdr_id, attr_name)) &
            call ncdr_error("nc_diag_read_get_global_attr_len != nc_diag_read_ret_global_attr_len!")
        
        write (*, "(A, I0, A, L)") "    -> Attribute: " // attr_name // " | Length: ", &
            attr_len
    end do
    
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id)
    
    tmp_ncdr_id = nc_diag_read_id_init("test_allcov.nc")
    
    deallocate(attr_names)
    
    write (*, "(A)") " ** File: test_allcov.nc (using ncdr_id)"
    
    call nc_diag_read_get_global_attr_names(tmp_ncdr_id, nattrs, nattrs_len, attr_names)
    write (*, "(A, I0, A, I0)") " ** Number of attributes in test_allcov.nc: ", nattrs, &
        " | Maximum length of attribute names: ", nattrs_len
    print *, "** All Attributes: **"
    print *, attr_names
    
    print *, "** Attribute details: **"
    
    do i = 1, nattrs
        attr_name = trim(attr_names(i))
        call nc_diag_read_assert_global_attr(tmp_ncdr_id, attr_name, attr_type, attr_len)
        if (nc_diag_read_check_global_attr(tmp_ncdr_id, attr_name) == .FALSE.) &
            call ncdr_error("Can't find attr with check(), even when it's listed!")
        call nc_diag_read_get_global_attr_len(tmp_ncdr_id, attr_name, attr_len)
        
        if (attr_len /= nc_diag_read_ret_global_attr_len(tmp_ncdr_id, attr_name)) &
            call ncdr_error("nc_diag_read_get_global_attr_len != nc_diag_read_ret_global_attr_len!")
        
        write (*, "(A, I0, A, L)") "    -> Attribute: " // attr_name // " | Length: ", &
            attr_len
    end do
    
    call nc_diag_read_close(file_ncdr_id = tmp_ncdr_id)
    
    contains
        subroutine display_1d_attr_byte(attr_name)
            character(len=*)                           :: attr_name
            integer(i_byte), dimension(:), allocatable :: attr_stor
            integer(i_byte), dimension(:), allocatable :: attr_stor2
            integer(i_byte)                            :: attr_stor_single
            integer(i_byte)                            :: attr_stor_single2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            
            if (any(attr_stor /= attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (size(attr_stor) == 1) then
                call nc_diag_read_get_global_attr(attr_name, attr_stor_single)
                call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor_single2)
                
                if ((attr_stor(1) /= attr_stor_single) .AND. (attr_stor(1) /= attr_stor_single2)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end if
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", size(attr_stor), ")"
            
            do i = 1, size(attr_stor)
                if (attr_stor(i) == NF90_FILL_BYTE) then
                    write (*, "(A4)") "(em)"
                else
                    write (*, "(I4)") attr_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
        end subroutine display_1d_attr_byte
        
        subroutine display_1d_attr_short(attr_name)
            character(len=*)                           :: attr_name
            integer(i_short),dimension(:), allocatable :: attr_stor
            integer(i_short),dimension(:), allocatable :: attr_stor2
            integer(i_short)                           :: attr_stor_single
            integer(i_short)                           :: attr_stor_single2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            
            if (any(attr_stor /= attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (size(attr_stor) == 1) then
                call nc_diag_read_get_global_attr(attr_name, attr_stor_single)
                call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor_single2)
                
                if ((attr_stor(1) /= attr_stor_single) .AND. (attr_stor(1) /= attr_stor_single2)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end if
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", size(attr_stor), ")"
            
            do i = 1, size(attr_stor)
                if (attr_stor(i) == NF90_FILL_SHORT) then
                    write (*, "(A6)") "(emp)"
                else
                    write (*, "(I6)") attr_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
        end subroutine display_1d_attr_short
        
        subroutine display_1d_attr_long(attr_name)
            character(len=*)                           :: attr_name
            integer(i_long), dimension(:), allocatable :: attr_stor
            integer(i_long), dimension(:), allocatable :: attr_stor2
            integer(i_long)                            :: attr_stor_single
            integer(i_long)                            :: attr_stor_single2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            
            if (any(attr_stor /= attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (size(attr_stor) == 1) then
                call nc_diag_read_get_global_attr(attr_name, attr_stor_single)
                call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor_single2)
                
                if ((attr_stor(1) /= attr_stor_single) .AND. (attr_stor(1) /= attr_stor_single2)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end if
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", size(attr_stor), ")"
            
            do i = 1, size(attr_stor)
                if (attr_stor(i) == NF90_FILL_INT) then
                    write (*, "(A12)") "(empty)"
                else
                    write (*, "(I12)") attr_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
            
        end subroutine display_1d_attr_long
        
        subroutine display_1d_attr_float(attr_name)
            character(len=*)                           :: attr_name
            real(r_single), dimension(:), allocatable  :: attr_stor
            real(r_single), dimension(:), allocatable  :: attr_stor2
            real(r_single)                             :: attr_stor_single
            real(r_single)                             :: attr_stor_single2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            if (any(attr_stor /= attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (size(attr_stor) == 1) then
                call nc_diag_read_get_global_attr(attr_name, attr_stor_single)
                call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor_single2)
                
                if ((attr_stor(1) /= attr_stor_single) .AND. (attr_stor(1) /= attr_stor_single2)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end if
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", size(attr_stor), ")"
            
            do i = 1, size(attr_stor)
                if (attr_stor(i) == NF90_FILL_FLOAT) then
                    write (*, "(A18)") "(empty)"
                else
                    write (*, "(F18.10)") attr_stor(i)
                end if
            end do
        end subroutine display_1d_attr_float
        
        subroutine display_1d_attr_double(attr_name)
            character(len=*)                           :: attr_name
            real(r_double), dimension(:), allocatable  :: attr_stor
            real(r_double), dimension(:), allocatable  :: attr_stor2
            real(r_double)                             :: attr_stor_single
            real(r_double)                             :: attr_stor_single2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            
            if (any(attr_stor /= attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (size(attr_stor) == 1) then
                call nc_diag_read_get_global_attr(attr_name, attr_stor_single)
                call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor_single2)
                
                if ((attr_stor(1) /= attr_stor_single) .AND. (attr_stor(1) /= attr_stor_single2)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end if
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", size(attr_stor), ")"
            
            do i = 1, size(attr_stor)
                if (attr_stor(i) == NF90_FILL_DOUBLE) then
                    write (*, "(A16)") "(empty)"
                else
                    write (*, "(F16.13)") attr_stor(i)
                end if
            end do
        end subroutine display_1d_attr_double
        
        subroutine display_1d_attr_string(attr_name)
            character(len=*)                           :: attr_name
            character(len=10000)                       :: attr_stor
            character(len=10000)                       :: attr_stor2
            character(len=:), allocatable              :: attr_stor_alloc
            character(len=:), allocatable              :: attr_stor_alloc2
            
            integer(i_long) :: i
            
            call nc_diag_read_get_global_attr(attr_name, attr_stor)
            call nc_diag_read_get_global_attr(tmp_ncdr_id, attr_name, attr_stor2)
            
            call nc_diag_read_noid_get_global_attr_1d_string(attr_name, attr_stor_alloc)
            call nc_diag_read_id_get_global_attr_1d_string(tmp_ncdr_id, attr_name, attr_stor_alloc2)
            
            if (len_trim(attr_stor) /= len_trim(attr_stor2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (len_trim(attr_stor) /= len_trim(attr_stor_alloc)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            if (len_trim(attr_stor) /= len_trim(attr_stor_alloc2)) &
                call ncdr_error("Storage with and without NCDR ID don't match!")
            
            do i = 1, len_trim(attr_stor)
                if (attr_stor(i:i) /= attr_stor2(i:i)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
                if (attr_stor(i:i) /= attr_stor_alloc(i:i)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
                if (attr_stor(i:i) /= attr_stor_alloc2(i:i)) &
                    call ncdr_error("Storage with and without NCDR ID don't match!")
            end do
            
            write (*, "(A, I0, A)") " ** Attribute (1D): " // attr_name // " (Elements: ", len_trim(attr_stor), ")"
            
            if ((attr_stor(1:1) == NF90_FILL_CHAR) .OR. (len(attr_stor) == 0)) then
                write (*, "(A)") "(empty)"
            else
                write (*, "(A)") '"' // trim(attr_stor) // '"'
            end if
            
            write (*, "(A)") ""
            
        end subroutine display_1d_attr_string
end program test_ncdr_attrs
