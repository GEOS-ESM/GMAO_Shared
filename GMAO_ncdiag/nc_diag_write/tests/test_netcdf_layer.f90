program test_netcdf_layer
    use kinds
    use netcdf
    use nc_diag_write_mod
    implicit none
    
    integer :: i
    real(r_single) :: f
    real(r_double) :: d
    
    character(len=100) :: str_header
    character(len=100) :: str_chaninfo
    character(len=100) :: str_metadata
    character(len=100) :: str_data2d
    
    character(len=10)  :: str_chaninfo_fixed
    character(len=10)  :: str_metadata_fixed
    character(len=10)  :: str_data2d_fixed(3)
    
    character(len=11)  :: str_chaninfo_fixed_bad
    character(len=11)  :: str_metadata_fixed_bad
    character(len=11)  :: str_data2d_fixed_bad(3)
    
    f = 1.234
    d = 2.34567890
    
    ! Enable info messages
    call nc_set_action_display(.TRUE.)
    call nc_set_info_display(.TRUE.)
    
    call nc_diag_set_trim(.FALSE.)
    
    !-----------------------------------------------------------------
    ! Fixed checks
    !-----------------------------------------------------------------
    call nc_diag_init("test_fixed.nc")
    
    call nc_diag_chaninfo_dim_set(5)
    
    str_chaninfo_fixed = "one"
    str_metadata_fixed = "two"
    str_data2d_fixed   = (/ "three", "four", "five" /)
    
    str_chaninfo_fixed_bad = "one"
    str_metadata_fixed_bad = "two"
    str_data2d_fixed_bad   = (/ "three", "four", "five" /)
    
    call nc_diag_set_trim(.FALSE.)
    call nc_diag_chaninfo("chaninfo_strfix", str_chaninfo_fixed)
    call nc_diag_chaninfo("chaninfo_strfix1", str_chaninfo_fixed)
    str_chaninfo_fixed = "three"
    call nc_diag_chaninfo("chaninfo_strfix", str_chaninfo_fixed)
    
    call nc_diag_metadata("metadata_strfix", str_metadata_fixed)
    call nc_diag_metadata("metadata_strfix1", str_metadata_fixed)
    str_metadata_fixed = "four"
    call nc_diag_metadata("metadata_strfix", str_metadata_fixed)
    
    call nc_diag_data2d("data2d_strfix", str_data2d_fixed)
    call nc_diag_data2d("data2d_strfix1", str_data2d_fixed)
    str_data2d_fixed   = (/ "six", "seven", "eight" /)
    call nc_diag_data2d("data2d_strfix", str_data2d_fixed)
    
    ! This would cause errors, due to the change in un-trimmed string
    ! length! (Uncomment to trigger the error!)
    !call nc_diag_chaninfo("chaninfo_strfix", str_chaninfo_fixed_bad)
    !call nc_diag_metadata("metadata_strfix", str_metadata_fixed_bad)
    !call nc_diag_data2d("data2d_strfix", str_data2d_fixed_bad)
    
    call nc_diag_write
    
    ! Fluid mode
    call nc_diag_set_trim(.TRUE.)
    
    call nc_diag_init("test.nc")
    
    ! Test init checking + corresponding error:
    !call nc_diag_init("test2.nc")
    
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
    call nc_diag_chaninfo_dim_set(10)
    
    do i = 1, 10
        call nc_diag_chaninfo("chaninfosimple1", i)
        call nc_diag_chaninfo("chaninfosimple2", i*2)
        call nc_diag_chaninfo("chaninfosimple4_float", f + 1.00)
        call nc_diag_chaninfo("chaninfosimple5_double", d + 1.00)
        
        call nc_diag_metadata("metadatasimple1", i)
        !call nc_diag_metadata("metadatasimple2", i*2)
        !call nc_diag_metadata("metadatasimple4_float", f + 1.00 + i)
        !call nc_diag_metadata("metadatasimple4_float2", f + 2.00 + i)
        !call nc_diag_metadata("metadatasimple5_double", d + 1.00 + i)
        
        call nc_diag_data2d("data2dsimple1", (/ i, i+1, i+2 /))
        call nc_diag_data2d("data2dsimple2", (/ i*2, i*3, i*4 /))
        call nc_diag_data2d("data2dsimple4_float", (/ f + 1.00 + i, f + 2.00 + i, f + 3.00 + i, f + 4.00 + i /))
        call nc_diag_data2d("data2dsimple4_float2", (/ f + 2.00 + i, f + 4.00 + i /))
        call nc_diag_data2d("data2dsimple5_double", (/ d + 1.00 + i /))
        call nc_diag_data2d("data2dsimple99", (/ i /))
        
        write(str_chaninfo, "(A, I0)") "ci6_", i
        call nc_diag_chaninfo("chaninfosimple6_str", str_chaninfo)
        
        write(str_metadata, "(A, I0)") "hellometa_", i
        call nc_diag_metadata("metadatasimple6_str", str_metadata)
    end do
    
    do i = 1, 9
        write(str_chaninfo, "(A, I0)") "ci_strings_", i
        call nc_diag_chaninfo("chaninfosimple7_str", str_chaninfo)
    end do
    
    !print *, "str_chaninfo:"
    !print *, str_chaninfo
    
    do i = 1, 9
        call nc_diag_chaninfo("chaninfosimple3_notcomplete", i*3)
    end do
    
    !do i = 1, 10000000
    do i = 1, 10000!000
        call nc_diag_header("headertestsimple", 123)
        
        call nc_diag_header("headertestsimple2_float", f)
        call nc_diag_header("headertestsimple3_double", d)
        
        write(str_header, "(A, I0)") "header_", i
        call nc_diag_header("headertestsimple4_str", str_header)
        
        !call nc_diag_metadata("metadatasimple7_big", i*2)
    end do
    
    do i = 1, 10!0000
        write(str_metadata, "(A, I0)") "morehellometa_", i
        call nc_diag_metadata("metadatasimple8_str", str_metadata)
        
        write(str_data2d, "(A, I0)") "data2d_", i
        call nc_diag_data2d("data2dsimple6_str", (/ str_data2d, "fill1", "fill2" /))
        
        ! This is broken... but it's an interesting testcase, as it breaks
        ! a LOT of stuff!
        ! index_llong = i needs to be commented out
        !call nc_diag_data2d("data2dsimple7", index_llong, (/ i, i+1, i+2 /))
        call nc_diag_data2d("data2dsimple7", (/ i, i+1, i+2 /))
    end do
    
    ! Add one entry... so we can test out valid/invalid data adding
    ! below!
    call nc_diag_chaninfo("chaninfosimple8_str", "test1234")
    
    ! ...and another one, for fun with buffered writing!
    call nc_diag_chaninfo("chaninfosimple9_buf", 3)
    
    ! Appending data variables
    call nc_diag_chaninfo("chaninfosimple10_notcomplete", 5678)
    
    call nc_diag_metadata("metadata_notcomplete", 1234)
    call nc_diag_metadata("metadata_notcomplete", 2234)
    call nc_diag_metadata("metadata_notcomplete", 3234)
    call nc_diag_metadata("metadata_notcomplete", 4234)
    call nc_diag_metadata("metadata_notcomplete", 5234)
    call nc_diag_metadata("metadata_notcomplete", 6234)
    call nc_diag_metadata("metadata_notcomplete", 7234)
    call nc_diag_metadata("metadata_notcomplete", 8234)
    call nc_diag_metadata("metadata_notcomplete", 9234)
    call nc_diag_metadata("metadata_notcomplete", 1234)
    call nc_diag_metadata("metadata_notcomplete", 2234)
    call nc_diag_metadata("metadata_notcomplete", 3234)
    call nc_diag_metadata("metadata_notcomplete", 4234)
    call nc_diag_metadata("metadata_notcomplete", 5234)
    call nc_diag_metadata("metadata_notcomplete", 6234)
    call nc_diag_metadata("metadata_notcomplete", 7234)
    call nc_diag_metadata("metadata_notcomplete", 8234)
    
    call nc_diag_metadata("metadata_str_notcomplete", "abcd")
    call nc_diag_metadata("metadata_str_notcomplete", "abc2")
    call nc_diag_metadata("metadata_str_notcomplete", "abc3")
    call nc_diag_metadata("metadata_str_notcomplete", "abc4")
    call nc_diag_metadata("metadata_str_notcomplete", "abc5")
    call nc_diag_metadata("metadata_str_notcomplete", "abc6")
    call nc_diag_data2d("data2d_notcomplete", (/ 1, 2, 3 /))
    
    ! Invalid buffered write test - we can't do any buffered write
    ! until we lock definitions:
    !call nc_diag_flush_buffer
    
    ! This, combined with nc_diag_set_strict(.TRUE.), will result in
    ! an error. We can also add 4 elements instead of 2 to achieve the
    ! same error. Note that after definition locking, adding more than 3
    ! elements (assuming we didn't add 4 elements here) won't work,
    ! since it violates the length check. (That means we won't even
    ! see the strict bounds checking error!)
    call nc_diag_data2d("data2dsimple1", (/ 2000, 4000 /))
    
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
    call nc_diag_varattr("data2dsimple7", "data2dsimple7_testattr1", "hi")
    call nc_diag_varattr("data2dsimple7", "data2dsimple7_testattr2", (/ 1, 2, 3 /))
    
    ! We can still add more data, but now we must adhere to the maximum
    ! variable length (the array input length).
    
    ! This is fine:
    call nc_diag_data2d("data2dsimple6_str", (/ "data2d_11", "fill1", "fill2" /))
    call nc_diag_data2d("data2dsimple7", (/ -1, -2, -3 /))
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_11")
    call nc_diag_chaninfo("chaninfosimple8_str", "test5678")
    
    ! This, however, is not. (Note that the array/string is longer than
    ! the others above.) (Uncomment the below lines to see what will
    ! happen!)
    !call nc_diag_data2d("data2dsimple6_str", int8(12), (/ "data2d_122", "fill1", "fill2" /))
    !call nc_diag_data2d("data2dsimple7", int8(12), (/ -4, -5, -6, -7 /))
    !call nc_diag_metadata("metadatasimple8_str", "morehellometa_111")
    !call nc_diag_chaninfo("chaninfosimple8_str", "test9101112")
    
    !------------------------------------------------------------------
    ! Buffered writing test!
    !------------------------------------------------------------------
    ! NOTE: For now, data2d does NOT have buffered writing enabled.
    !       This will be fixed in a future release.
    
    call nc_diag_chaninfo("chaninfosimple9_buf", 6)
    call nc_diag_chaninfo("chaninfosimple9_buf", 9)
    
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b1")
    call nc_diag_metadata("metadatasimple6_str", "meta_b1")
    call nc_diag_metadata("metadatasimple1", 100)
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b2")
    call nc_diag_metadata("metadatasimple6_str", "meta_b2")
    
    call nc_diag_data2d("data2dsimple1", (/ 1000, 2000, 3000 /))
    call nc_diag_data2d("data2dsimple1", (/ 2000, 4000, 6000 /))
    call nc_diag_data2d("data2dsimple2", (/ 1111, 2222, 3333 /))
    call nc_diag_data2d("data2dsimple2", (/ 2222, 4444, 6666 /))
    call nc_diag_data2d("data2dsimple6_str", (/ "mwahahaha", "arrrrrgh", "grrrrowwl" /))
    call nc_diag_data2d("data2dsimple6_str", (/ "boink", "kabam", "peekaboo" /))
    call nc_diag_data2d("data2dsimple7", (/ 20, 40, 60 /))
    call nc_diag_data2d("data2dsimple7", (/ 40, 80, 120 /))
    
    print *, "Attempting to flush buf 1:"
    call nc_diag_flush_buffer
    
    call nc_diag_chaninfo("chaninfosimple9_buf", 12)
    call nc_diag_chaninfo("chaninfosimple9_buf", 15)
    call nc_diag_chaninfo("chaninfosimple9_buf", 18)
    call nc_diag_chaninfo("chaninfosimple9_buf", 21)
    
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b3")
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b4")
    call nc_diag_metadata("metadatasimple6_str", "meta_b3")
    call nc_diag_metadata("metadatasimple6_str", "meta_b4")
    call nc_diag_metadata("metadatasimple1", 200)
    
    ! We can add something in the future!
    call nc_diag_data2d("data2dsimple1", (/ -1000, -2000, -3000 /))
    call nc_diag_data2d("data2dsimple6_str", (/ "aaaaaaaaa", "bbbbbbbb", "ccccccccc" /))
    call nc_diag_data2d("data2dsimple7", (/ 4000, 8000, 12000 /))
    
    print *, "Attempting to flush buf 2:"
    call nc_diag_flush_buffer
    
    call nc_diag_chaninfo("chaninfosimple9_buf", 24)
    call nc_diag_chaninfo("chaninfosimple9_buf", 27)
    call nc_diag_chaninfo("chaninfosimple9_buf", 30)
    
    call nc_diag_metadata("metadatasimple1", 300)
    call nc_diag_metadata("metadatasimple6_str", "meta_b5")
    call nc_diag_metadata("metadatasimple6_str", "meta_b6")
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b5")
    call nc_diag_metadata("metadatasimple8_str", "morehellometa_b6")
    
    ! We can still change an old value at the end!
    call nc_diag_data2d("data2dsimple1", (/ 2000, 4000, 6000 /))
    call nc_diag_data2d("data2dsimple2", (/ 1111, 2222, 3333 /))
    
    call nc_diag_data2d("data2dsimple1", (/ 4000, 6000, 8000 /))
    call nc_diag_data2d("data2dsimple2", (/ 2222, 4444, 6666 /))
    
    call nc_diag_data2d("data2dsimple1", (/ 6000, 8000, 10000 /))
    call nc_diag_data2d("data2dsimple2", (/ 3333, 6666, 9999 /))
    
    ! Out of order is fine too!
    call nc_diag_data2d("data2dsimple6_str", (/ "mwahahaha", "arrrrrgh", "grrrrowwl" /))
    call nc_diag_data2d("data2dsimple7", (/ 20, 40, 60 /))
    
    call nc_diag_data2d("data2dsimple7", (/ 200, 400, 600 /))
    call nc_diag_data2d("data2dsimple6_str", (/ "asdfghjk", "zxcvbnm", "qwerty" /))
    
    call nc_diag_data2d("data2dsimple6_str", (/ "boink", "kabam", "peekaboo" /))
    call nc_diag_data2d("data2dsimple7", (/ 40, 80, 120 /))
    
    ! Even with buffering, you still can't overwrite nchans...
    ! (The following line, if uncommented, should result in an error!)
    !call nc_diag_chaninfo("chaninfosimple9_buf", 33)
    
    ! Back to header stuff...
    call nc_diag_header("headertestsimple5_str", "hello world")
    
    print *, "str_header:"
    print *, str_header
    
    print *, "===================="
    print *, "Vector:"
    print *, "===================="
    
    do i = 1, 1000
        call nc_diag_header("headertestarr1", (/ 123, 234, 345, 456, 567, 678, 789 /))
    end do
    
    call nc_diag_header("headertestarr2", (/ 222, 234, 345, 456, 567, 678, 789 /))
    call nc_diag_header("headertestarr3", (/ 333, 234, 345, 456, 567, 678, 789 /))
    call nc_diag_header("headertestarr4", (/ 444, 234, 345, 456, 567, 678, 789 /))
    call nc_diag_header("headertestarr5", (/ 111, 222, 333, 444, 555, 666, 777, 888, 999 /))
    call nc_diag_header("headertestarr6", (/ 999, 777, 555, 333, 111 /))
    call nc_diag_header("headertestsimple2", 123)
    call nc_diag_header("headertestsimple3", 321)
    call nc_diag_header("headertestarr7", (/ 111, 222, 333, 444, 555, 666, 777, 888, 999 /))
    call nc_diag_header("headertestarr7", (/ 222, 444, 666, 888 /))
    call nc_diag_header("headertestarr7", 999)
    
    call nc_diag_header("headertestarr8", (/ f, f*2, f*3, f*4 /))
    call nc_diag_header("headertestarr9", (/ d, d*2, d*3, d*4 /))
    
    ! nc_diag_header does not support arrays of strings... because
    ! NetCDF4 doesn't support it, either!
    ! (At least, that's what I remember from the docs... double check
    ! to make sure!)
    
    print *, "==============================="
    print *, "Writing resulting NetCDF file:"
    print *, "==============================="
    
    call nc_diag_write
    
    ! Appending - we can reopen the same file in append mode!
    print *, "==============================="
    print *, "Appending to NetCDF file:"
    print *, "==============================="
    call nc_diag_init("test.nc", .TRUE.)
    
    call nc_diag_chaninfo("chaninfosimple8_str", "test1")
    call nc_diag_chaninfo("chaninfosimple8_str", "test2")
    call nc_diag_chaninfo("chaninfosimple10_notcomplete", 1)
    call nc_diag_chaninfo("chaninfosimple10_notcomplete", 2)
    
    call nc_diag_metadata("metadata_notcomplete", 5678)
    call nc_diag_metadata("metadata_notcomplete", 6789)
    
    call nc_diag_metadata("metadata_str_notcomplete", "efgh")
    
    ! This will fail due to length constraints...
    !call nc_diag_metadata("metadata_str_notcomplete", "efghh")
    
    call nc_diag_data2d("data2d_notcomplete", (/ 2, 3, 4 /))
    
    ! This will also fail due to longer length
    !call nc_diag_data2d("data2d_notcomplete", (/ 2, 3, 4, 5 /))
    
    ! Flushing works, too!
    call nc_diag_flush_buffer
    
    call nc_diag_chaninfo("chaninfosimple8_str", "test3")
    call nc_diag_chaninfo("chaninfosimple8_str", "test4")
    call nc_diag_chaninfo("chaninfosimple10_notcomplete", 3)
    call nc_diag_chaninfo("chaninfosimple10_notcomplete", 4)
    
    call nc_diag_metadata("metadata_notcomplete", 7898)
    call nc_diag_metadata("metadata_notcomplete", 8989)
    
    call nc_diag_metadata("metadata_str_notcomplete", "ijkl")
    
    call nc_diag_write
end program test_netcdf_layer
