program test_ncdr_get
    use kinds, only: i_byte, i_short, i_long, r_single, r_double
    use nc_diag_read, only: nc_diag_read_init, nc_diag_read_close, &
        nc_diag_read_get_dim_names, &
        nc_diag_read_get_var, &
        nc_diag_read_get_var_type, &
        nc_diag_read_get_var_ndims, &
        nc_diag_read_get_var_names, &
        ncdr_error
    use netcdf, only: NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, &
        NF90_DOUBLE, NF90_CHAR, NF90_FILL_BYTE, NF90_FILL_SHORT, &
        NF90_FILL_INT, NF90_FILL_FLOAT, NF90_FILL_DOUBLE, NF90_FILL_CHAR
    
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
    
    call nc_diag_read_get_var_names(nvars, nvars_len, var_names)
    write (*, "(A, I0, A, I0)") " ** Number of variables: ", nvars, &
        " | Maximum length of variable names: ", nvars_len
    print *, "** All variables: **"
    print *, var_names
    
    !call nc_diag_read_init("test_fixed.nc")
    call display_1d_var_long("chaninfosimple1")
    call display_1d_var_long("chaninfosimple2")
    
    ! This won't work due to mismatched types:
    !call display_1d_var_byte("chaninfosimple2")
    
    call display_1d_var_float("chaninfosimple4_float")
    
    ! Mismatched types:
    !call display_1d_var_long("chaninfosimple4_float")
    !call display_1d_var_double("chaninfosimple4_float")
    
    call display_1d_var_double("chaninfosimple5_double")
    
    call display_1d_var_string("chaninfosimple6_str")
    call display_1d_var_string("chaninfosimple7_str")
    
    call display_1d_var_long("chaninfosimple3_notcomplete")
    
    call display_1d_var_string("chaninfosimple8_str")
    
    call display_1d_var_long("chaninfosimple9_buf")
    call display_1d_var_long("chaninfosimple10_notcomplete")
    
    call display_1d_var_long("metadatasimple1")
    
    call display_1d_var_string("metadatasimple6_str")
    call display_1d_var_string("metadatasimple8_str")
    
    call display_1d_var_long("metadata_notcomplete")
    
    call display_1d_var_string("metadata_str_notcomplete")
    
    call display_2d_var_long("data2dsimple1")
    call display_2d_var_long("data2dsimple2")
    
    call display_2d_var_float("data2dsimple4_float")
    call display_2d_var_float("data2dsimple4_float2")
    
    call display_2d_var_double("data2dsimple5_double")
    
    call display_2d_var_long("data2dsimple99")
    
    call display_2d_var_string("data2dsimple6_str")
    
    call display_2d_var_long("data2dsimple7")
    call display_2d_var_long("data2d_notcomplete")
    
    call nc_diag_read_close
    
    deallocate(dim_names)
    deallocate(var_names)
    
    call nc_diag_read_init("test_allcov.nc")
    
    call nc_diag_read_get_dim_names(ndims, ndims_len, dim_names)
    write (*, "(A, I0, A, I0)") " ** Number of dimensions: ", ndims, &
        " | Maximum length of dimension names: ", ndims_len
    print *, "** All dimensions: **"
    print *, dim_names
    
    call nc_diag_read_get_var_names(nvars, nvars_len, var_names)
    write (*, "(A, I0, A, I0)") " ** Number of variables: ", nvars, &
        " | Maximum length of variable names: ", nvars_len
    print *, "** All variables: **"
    print *, var_names
    
    do i = 1, nvars
        var_type = nc_diag_read_get_var_type(var_names(i))
        var_ndims = nc_diag_read_get_var_ndims(var_names(i))
        
        if ((var_ndims == 1) .OR. &
            ((var_ndims == 2) .AND. (var_type == NF90_CHAR))) then
            if (var_type == NF90_BYTE) then
                call display_1d_var_byte(var_names(i))
            else if (var_type == NF90_SHORT) then
                call display_1d_var_short(var_names(i))
            else if (var_type == NF90_INT) then
                call display_1d_var_long(var_names(i))
            else if (var_type == NF90_FLOAT) then
                call display_1d_var_float(var_names(i))
            else if (var_type == NF90_DOUBLE) then
                call display_1d_var_double(var_names(i))
            else if (var_type == NF90_CHAR) then
                call display_1d_var_string(var_names(i))
            else
                call ncdr_error("Invalid type!")
            end if
        else if (((var_ndims == 2) .AND. (var_type /= NF90_CHAR)) &
            .OR. ((var_ndims == 3) .AND. (var_type == NF90_CHAR))) then
            if (var_type == NF90_BYTE) then
                call display_2d_var_byte(var_names(i))
            else if (var_type == NF90_SHORT) then
                call display_2d_var_short(var_names(i))
            else if (var_type == NF90_INT) then
                call display_2d_var_long(var_names(i))
            else if (var_type == NF90_FLOAT) then
                call display_2d_var_float(var_names(i))
            else if (var_type == NF90_DOUBLE) then
                call display_2d_var_double(var_names(i))
            else if (var_type == NF90_CHAR) then
                call display_2d_var_string(var_names(i))
            else
                call ncdr_error("Invalid type!")
            end if
        else
            call ncdr_error("Invalid ndims!")
        end if
    end do
    
    contains
        subroutine display_1d_var_byte(var_name)
            character(len=*)                           :: var_name
            integer(i_byte), dimension(:), allocatable :: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if (var_stor(i) == NF90_FILL_INT) then
                    write (*, "(A4)") "(em)"
                else
                    write (*, "(I4)") var_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
        end subroutine display_1d_var_byte
        
        subroutine display_1d_var_short(var_name)
            character(len=*)                           :: var_name
            integer(i_short), dimension(:), allocatable :: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if (var_stor(i) == NF90_FILL_INT) then
                    write (*, "(A6)") "(emp)"
                else
                    write (*, "(I6)") var_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
        end subroutine display_1d_var_short
        
        subroutine display_1d_var_long(var_name)
            character(len=*)                           :: var_name
            integer(i_long), dimension(:), allocatable :: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if (var_stor(i) == NF90_FILL_INT) then
                    write (*, "(A12)") "(empty)"
                else
                    write (*, "(I12)") var_stor(i)
                end if
            end do
            
            write (*, "(A)") ""
            
        end subroutine display_1d_var_long
        
        subroutine display_1d_var_float(var_name)
            character(len=*)                           :: var_name
            real(r_single), dimension(:), allocatable  :: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if (var_stor(i) == NF90_FILL_FLOAT) then
                    write (*, "(A18)") "(empty)"
                else
                    write (*, "(F18.10)") var_stor(i)
                end if
            end do
        end subroutine display_1d_var_float
        
        subroutine display_1d_var_double(var_name)
            character(len=*)                           :: var_name
            real(r_double), dimension(:), allocatable  :: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if (var_stor(i) == NF90_FILL_DOUBLE) then
                    write (*, "(A16)") "(empty)"
                else
                    write (*, "(F16.13)") var_stor(i)
                end if
            end do
        end subroutine display_1d_var_double
        
        subroutine display_1d_var_string(var_name)
            character(len=*)                           :: var_name
            character(len=:), dimension(:), allocatable:: var_stor
            
            integer(i_long) :: i
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (1D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor)
                if ((var_stor(i)(1:1) == NF90_FILL_CHAR) .OR. (len(var_stor(i)) == 0)) then
                    write (*, "(A20)") "(empty)"
                else
                    write (*, "(A20)") '"' // var_stor(i) // '"'
                end if
            end do
            
            write (*, "(A)") ""
            
        end subroutine display_1d_var_string
        
        subroutine display_2d_var_byte(var_name)
            character(len=*)                           :: var_name
            integer(i_byte),dimension(:,:),allocatable :: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if (var_stor(j, i) == NF90_FILL_BYTE) then
                        write (*, "(A5)", advance = "no") "(e) "
                    else
                        write (*, "(I4, A)", advance = "no") var_stor(j, i), " "
                    end if
                end do
                write (*, "(A)") ""
            end do
        end subroutine display_2d_var_byte
        
        subroutine display_2d_var_short(var_name)
            character(len=*)                           :: var_name
            integer(i_short),dimension(:,:),allocatable :: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if (var_stor(j, i) == NF90_FILL_SHORT) then
                        write (*, "(A7)", advance = "no") "(emp) "
                    else
                        write (*, "(I6, A)", advance = "no") var_stor(j, i), " "
                    end if
                end do
                write (*, "(A)") ""
            end do
        end subroutine display_2d_var_short
        
        subroutine display_2d_var_long(var_name)
            character(len=*)                           :: var_name
            integer(i_long),dimension(:,:),allocatable :: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if (var_stor(j, i) == NF90_FILL_INT) then
                        write (*, "(A13)", advance = "no") "(empty) "
                    else
                        write (*, "(I12, A)", advance = "no") var_stor(j, i), " "
                    end if
                end do
                write (*, "(A)") ""
            end do
        end subroutine display_2d_var_long
        
        subroutine display_2d_var_float(var_name)
            character(len=*)                           :: var_name
            real(r_single), dimension(:,:), allocatable:: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if (var_stor(j, i) == NF90_FILL_FLOAT) then
                        write (*, "(A19)", advance = "no") "(empty) "
                    else
                        write (*, "(F18.10, A)", advance = "no") var_stor(j, i), " "
                    end if
                end do
                write (*, "(A)") ""
            end do
        end subroutine display_2d_var_float
        
        subroutine display_2d_var_double(var_name)
            character(len=*)                           :: var_name
            real(r_double), dimension(:,:), allocatable:: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if (var_stor(j, i) == NF90_FILL_DOUBLE) then
                        write (*, "(A17)", advance = "no") "(empty) "
                    else
                        write (*, "(F16.13, A)", advance = "no") var_stor(j, i), " "
                    end if
                end do
                write (*, "(A)") ""
            end do
            
        end subroutine display_2d_var_double
        
        ! NOTE - dimensions have to be flipped
        subroutine display_2d_var_string(var_name)
            character(len=*)                           :: var_name
            character(len=:),dimension(:,:),allocatable:: var_stor
            
            integer(i_long) :: i, j
            
            call nc_diag_read_get_var(var_name, var_stor)
            
            write (*, "(A, I0, A)") " ** Variable (2D): " // var_name // " (Elements: ", size(var_stor), ")"
            
            print *, shape(var_stor)
            
            do i = 1, size(var_stor, 2)
                do j = 1, size(var_stor, 1)
                    if ((j > 1) .AND. (mod(j - 1, 5) == 0)) write (*, "(A)") "..."
                    if ((var_stor(j, i)(1:1) == NF90_FILL_CHAR) .OR. (len(var_stor(j, i)) == 0)) then
                        write (*, "(A20)", advance = "no") "(empty) "
                    else
                        write (*, "(A20)", advance = "no") '"' // var_stor(j, i) // '" '
                    end if
                end do
                write (*, "(A)") ""
            end do
            
            write (*, "(A)") ""
            
        end subroutine display_2d_var_string
end program test_ncdr_get
