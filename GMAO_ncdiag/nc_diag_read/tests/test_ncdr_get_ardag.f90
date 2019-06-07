! ardag:
!   -> ALLOCATE dimensions by fetching dimensions and allocating
!   -> RETURN DIM: dimension fetching method by returning dimensions
!   -> ALLOCATED GET: allocated dimensions manually with ndims before
!                     fetch, with additional checks from get()
program test_ncdr_get_ardag
    use kinds, only: i_byte, i_short, i_long, r_single, r_double
    use nc_diag_read, only: nc_diag_read_init, nc_diag_read_close, &
        nc_diag_read_get_var, &
        nc_diag_read_get_var_ndims, &
        nc_diag_read_get_var_dims
    use netcdf, only: NF90_FILL_BYTE, NF90_FILL_SHORT, NF90_FILL_INT, &
        NF90_FILL_FLOAT, NF90_FILL_DOUBLE, NF90_FILL_CHAR
    
    implicit none
    
    call nc_diag_read_init("test.nc")
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
    
    contains
        subroutine display_1d_var_byte(var_name)
            character(len=*)                           :: var_name
            integer(i_byte), dimension(:), allocatable :: var_stor
            
            integer(i_long) :: i, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1)))
            
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
            
            integer(i_long) :: i, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1)))
            
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
            
            integer(i_long) :: i, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1)))
            
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
            
            integer(i_long) :: i, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1)))
            
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
            
            integer(i_long) :: i, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1)))
            
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
            
            integer(i_long) :: i, ndims, dim_1, dim_2
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            dim_1 = var_dims(1)
            dim_2 = var_dims(2)
            allocate(character(len=dim_1) :: var_stor(dim_2))
            
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
            
            integer(i_long) :: i, j, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1), var_dims(2)))
            
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
            
            integer(i_long) :: i, j, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1), var_dims(2)))
            
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
            
            integer(i_long) :: i, j, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1), var_dims(2)))
            
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
            
            integer(i_long) :: i, j, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1), var_dims(2)))
            
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
            
            integer(i_long) :: i, j, ndims
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            allocate(var_stor(var_dims(1), var_dims(2)))
            
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
            
            integer(i_long) :: i, j, ndims, dim_1, dim_2, dim_3
            integer(i_long), dimension(:), allocatable :: var_dims
            
            ndims = nc_diag_read_get_var_ndims(var_name)
            allocate(var_dims(ndims))
            call nc_diag_read_get_var_dims(var_name, var_dims = var_dims)
            
            dim_1 = var_dims(1)
            dim_2 = var_dims(2)
            dim_3 = var_dims(3)
            
            allocate(character(len=dim_1) :: var_stor(dim_2, dim_3))
            
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
end program test_ncdr_get_ardag
