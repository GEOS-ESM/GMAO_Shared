program test_nclayer_stress
    use kinds
    use netcdf
    use nc_diag_write_mod
    implicit none
    
    integer(i_long)   :: i, j
    
    integer(i_byte)   :: int_byte
    integer(i_short)  :: int_short
    integer(i_long)   :: int_long
    real(r_single)    :: real_single
    real(r_double)    :: real_double
    character(len=30) :: string
    
    integer(i_byte)  , dimension(5) :: int_byte_arr
    integer(i_short) , dimension(5) :: int_short_arr
    integer(i_long)  , dimension(5) :: int_long_arr
    real(r_single)   , dimension(5) :: real_single_arr
    real(r_double)   , dimension(5) :: real_double_arr
    character(len=30), dimension(5) :: string_arr
    
    ! Enable info messages
    call nc_set_action_display(.TRUE.)
    call nc_set_info_display(.TRUE.)
    
    call nc_diag_set_trim(.FALSE.)
    
    ! Fluid mode
    call nc_diag_set_trim(.TRUE.)
    
    call nc_diag_init("test_stress.nc")
    
    ! Uncomment below line to enable strict mode, aka strict
    ! variable bounds checking. When strict mode is enabled, if any
    ! variables have different lengths from each other, or if there are
    ! any differences between the variable rows (e.g. an uneven array),
    ! an error will occur and the program will halt.
    
    ! call nc_diag_set_strict(.TRUE.)
    
    print *, "===================="
    print *, "Single:"
    print *, "===================="
    !call nc_diag_header("test1", (/ 123, 234, 345, 456, 567, 678, 789 /))
    !call nc_diag_header("test2", (/ 123, 234, 345, 456, 567, 678, 789 /))
    ! 100,000,000
    call nc_diag_chaninfo_dim_set(10000)
    
    do i = 1, 4999
        int_byte    = i
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_chaninfo("chaninfosimple_ibyte",   int_byte)
        call nc_diag_chaninfo("chaninfosimple_ishort",  int_short)
        call nc_diag_chaninfo("chaninfosimple_ilong",   int_long)
        call nc_diag_chaninfo("chaninfosimple_rsingle", real_single)
        call nc_diag_chaninfo("chaninfosimple_rdouble", real_double)
        call nc_diag_chaninfo("chaninfosimple_string",  string)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        ! Header
        call nc_diag_header("headersimple_ibyte",   int_byte)
        call nc_diag_header("headersimple_ishort",  int_short)
        call nc_diag_header("headersimple_ilong",   int_long)
        call nc_diag_header("headersimple_rsingle", real_single)
        call nc_diag_header("headersimple_rdouble", real_double)
        call nc_diag_header("headersimple_string",  string)
        call nc_diag_header("headerarr_ibyte",   int_byte_arr)
        call nc_diag_header("headerarr_ishort",  int_short_arr)
        call nc_diag_header("headerarr_ilong",   int_long_arr)
        call nc_diag_header("headerarr_rsingle", real_single_arr)
        call nc_diag_header("headerarr_rdouble", real_double_arr)
    end do
    
    !------------------------------------------------------------------
    ! Variable attribute test! (With definition locking on the side!)
    !------------------------------------------------------------------
    
    ! In order for variable attributes to work, we MUST call
    ! nc_diag_lock_def! This is due to the fact that we need the NetCDF
    ! variable IDs in order for attribute defining to work, and
    ! the variable IDs aren't created until the variables definitions
    ! have been created (and locked)!
    call nc_diag_lock_def
    
    ! Now we can add variable attributes!
    do i = 5000, 6999
        int_byte    = i
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_chaninfo("chaninfosimple_ibyte",   int_byte)
        call nc_diag_chaninfo("chaninfosimple_ishort",  int_short)
        call nc_diag_chaninfo("chaninfosimple_ilong",   int_long)
        call nc_diag_chaninfo("chaninfosimple_rsingle", real_single)
        call nc_diag_chaninfo("chaninfosimple_rdouble", real_double)
        call nc_diag_chaninfo("chaninfosimple_string",  string)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        call nc_diag_varattr("chaninfosimple_ibyte", "int_byte_simple",    int_byte)
        call nc_diag_varattr("chaninfosimple_ibyte", "int_short_simple",   int_short)
        call nc_diag_varattr("chaninfosimple_ibyte", "int_long_simple",    int_long)
        call nc_diag_varattr("chaninfosimple_ibyte", "real_single_simple", real_single)
        call nc_diag_varattr("chaninfosimple_ibyte", "real_double_simple", real_double)
        call nc_diag_varattr("chaninfosimple_ibyte", "string_simple",      string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        call nc_diag_varattr("chaninfosimple_ibyte", "int_byte_arr",    int_byte_arr)
        call nc_diag_varattr("chaninfosimple_ibyte", "int_short_arr",   int_short_arr)
        call nc_diag_varattr("chaninfosimple_ibyte", "int_long_arr",    int_long_arr)
        call nc_diag_varattr("chaninfosimple_ibyte", "real_single_arr", real_single_arr)
        call nc_diag_varattr("chaninfosimple_ibyte", "real_double_arr", real_double_arr)
    end do
    
    !------------------------------------------------------------------
    ! Buffered writing test!
    !------------------------------------------------------------------
    ! NOTE: For now, data2d does NOT have buffered writing enabled.
    !       This will be fixed in a future release.
    
    do i = 7000, 7999
        int_byte    = 12
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_chaninfo("chaninfosimple_ibyte",   int_byte)
        call nc_diag_chaninfo("chaninfosimple_ishort",  int_short)
        call nc_diag_chaninfo("chaninfosimple_ilong",   int_long)
        call nc_diag_chaninfo("chaninfosimple_rsingle", real_single)
        call nc_diag_chaninfo("chaninfosimple_rdouble", real_double)
        call nc_diag_chaninfo("chaninfosimple_string",  string)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        ! Header
        call nc_diag_header("headersimple_ibyte",   int_byte)
        call nc_diag_header("headersimple_ishort",  int_short)
        call nc_diag_header("headersimple_ilong",   int_long)
        call nc_diag_header("headersimple_rsingle", real_single)
        call nc_diag_header("headersimple_rdouble", real_double)
        call nc_diag_header("headersimple_string",  string)
        call nc_diag_header("headerarr_ibyte",   int_byte_arr)
        call nc_diag_header("headerarr_ishort",  int_short_arr)
        call nc_diag_header("headerarr_ilong",   int_long_arr)
        call nc_diag_header("headerarr_rsingle", real_single_arr)
        call nc_diag_header("headerarr_rdouble", real_double_arr)
    end do
    
    print *, "Attempting to flush buf 1:"
    call nc_diag_flush_buffer
    
    do i = 8000, 8999
        int_byte    = 123
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_chaninfo("chaninfosimple_ibyte",   int_byte)
        call nc_diag_chaninfo("chaninfosimple_ishort",  int_short)
        call nc_diag_chaninfo("chaninfosimple_ilong",   int_long)
        call nc_diag_chaninfo("chaninfosimple_rsingle", real_single)
        call nc_diag_chaninfo("chaninfosimple_rdouble", real_double)
        call nc_diag_chaninfo("chaninfosimple_string",  string)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        ! Header
        call nc_diag_header("headersimple_ibyte",   int_byte)
        call nc_diag_header("headersimple_ishort",  int_short)
        call nc_diag_header("headersimple_ilong",   int_long)
        call nc_diag_header("headersimple_rsingle", real_single)
        call nc_diag_header("headersimple_rdouble", real_double)
        call nc_diag_header("headersimple_string",  string)
        call nc_diag_header("headerarr_ibyte",   int_byte_arr)
        call nc_diag_header("headerarr_ishort",  int_short_arr)
        call nc_diag_header("headerarr_ilong",   int_long_arr)
        call nc_diag_header("headerarr_rsingle", real_single_arr)
        call nc_diag_header("headerarr_rdouble", real_double_arr)
    end do
    
    print *, "Attempting to flush buf 2:"
    call nc_diag_flush_buffer
    
    do i = 9000, 10000
        int_byte    = 23
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_chaninfo("chaninfosimple_ibyte",   int_byte)
        call nc_diag_chaninfo("chaninfosimple_ishort",  int_short)
        call nc_diag_chaninfo("chaninfosimple_ilong",   int_long)
        call nc_diag_chaninfo("chaninfosimple_rsingle", real_single)
        call nc_diag_chaninfo("chaninfosimple_rdouble", real_double)
        call nc_diag_chaninfo("chaninfosimple_string",  string)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        ! Header
        call nc_diag_header("headersimple_ibyte",   int_byte)
        call nc_diag_header("headersimple_ishort",  int_short)
        call nc_diag_header("headersimple_ilong",   int_long)
        call nc_diag_header("headersimple_rsingle", real_single)
        call nc_diag_header("headersimple_rdouble", real_double)
        call nc_diag_header("headersimple_string",  string)
        call nc_diag_header("headerarr_ibyte",   int_byte_arr)
        call nc_diag_header("headerarr_ishort",  int_short_arr)
        call nc_diag_header("headerarr_ilong",   int_long_arr)
        call nc_diag_header("headerarr_rsingle", real_single_arr)
        call nc_diag_header("headerarr_rdouble", real_double_arr)
    end do
    
    print *, "==============================="
    print *, "Writing resulting NetCDF file:"
    print *, "==============================="
    
    call nc_diag_write
    
    ! Appending - we can reopen the same file in append mode!
    print *, "==============================="
    print *, "Appending to NetCDF file:"
    print *, "==============================="
    call nc_diag_init("test_stress.nc", .TRUE.)
    
    do i = 10001, 11000
        int_byte    = 45
        int_short   = i * 2
        int_long    = i * 3
        real_single = i * 4 + 0.1234
        real_double = i * 5 + 0.12345678
        
        write (string, "(A, I0)") "mystr_", (i * 6)
        
        call nc_diag_metadata("metadatasimple_ibyte",   int_byte)
        call nc_diag_metadata("metadatasimple_ishort",  int_short)
        call nc_diag_metadata("metadatasimple_ilong",   int_long)
        call nc_diag_metadata("metadatasimple_rsingle", real_single)
        call nc_diag_metadata("metadatasimple_rdouble", real_double)
        call nc_diag_metadata("metadatasimple_string",  string)
        
        ! 2D arrays
        int_byte_arr    = (/ i, i + 1, i + 2, i + 3, i + 4 /)
        int_short_arr   = (/ 2 * i, (2 * i) + 1, (2 * i) + 2, (2 * i) + 3, (2 * i) + 4 /)
        int_long_arr    = (/ 3 * i, (3 * i) + 1, (3 * i) + 2, (3 * i) + 3, (3 * i) + 4 /)
        real_single_arr = 4 * i
        real_double_arr = 5 * i
        
        do j = 1, 5
            real_single_arr(j) = real_single_arr(j) + 0.2345 + j
            real_double_arr(j) = real_double_arr(j) + 0.23456789 + j
        end do
        
        do j = 1, 5
            write (string_arr(j), "(A, I0, A, I0)") "mystr_", j, "_", (i * 6)
        end do
        
        call nc_diag_data2d("data2darr_ibyte",   int_byte_arr)
        call nc_diag_data2d("data2darr_ishort",  int_short_arr)
        call nc_diag_data2d("data2darr_ilong",   int_long_arr)
        call nc_diag_data2d("data2darr_rsingle", real_single_arr)
        call nc_diag_data2d("data2darr_rdouble", real_double_arr)
        call nc_diag_data2d("data2darr_string",  string_arr)
        
        ! Header
        call nc_diag_header("headersimple_ibyte",   int_byte)
        call nc_diag_header("headersimple_ishort",  int_short)
        call nc_diag_header("headersimple_ilong",   int_long)
        call nc_diag_header("headersimple_rsingle", real_single)
        call nc_diag_header("headersimple_rdouble", real_double)
        call nc_diag_header("headersimple_string",  string)
        call nc_diag_header("headerarr_ibyte",   int_byte_arr)
        call nc_diag_header("headerarr_ishort",  int_short_arr)
        call nc_diag_header("headerarr_ilong",   int_long_arr)
        call nc_diag_header("headerarr_rsingle", real_single_arr)
        call nc_diag_header("headerarr_rdouble", real_double_arr)
    end do
    
    call nc_diag_write
end program test_nclayer_stress
