program gen_diff_size_netcdf
    use kinds
    use netcdf
    use nc_diag_write_mod
    implicit none
    
    integer :: i
    
    character(len=100) :: str_chaninfo
    character(len=100) :: str_metadata
    character(len=100) :: str_data2d
    
    ! Enable info messages
    call nc_set_action_display(.TRUE.)
    call nc_set_info_display(.TRUE.)
    
    call nc_diag_init("diff1.nc")
    
    print *, "===================="
    print *, "NetCDF file #1:"
    print *, "===================="
    
    call nc_diag_chaninfo_dim_set(10)
    
    do i = 1, 10
        write(str_chaninfo, "(A, I0)") "ci6_", i
        call nc_diag_chaninfo("chaninfosimple6_str", str_chaninfo)
        
        write(str_metadata, "(A, I0)") "hellometa_", i
        call nc_diag_metadata("metadatasimple6_str", str_metadata)
    end do
    
    do i = 1, 10!0000
        write(str_metadata, "(A, I0)") "morehellometa_", i
        call nc_diag_metadata("metadatasimple8_str", str_metadata)
        
        write(str_data2d, "(A, I0)") "data2d_", i
        call nc_diag_data2d("data2dsimple6_str", (/ str_data2d, "fill1", "fill2" /))
    end do
    
    print *, "==============================="
    print *, "Writing resulting NetCDF file:"
    print *, "==============================="
    
    call nc_diag_write
    
    ! Appending - we can reopen the same file in append mode!
    print *, "==============================="
    print *, "NetCDF file #2:"
    print *, "==============================="
    call nc_diag_init("diff2.nc")
    
    call nc_diag_chaninfo_dim_set(10)
    
    do i = 1, 10
        write(str_chaninfo, "(A, I0)") "ci66_", i
        call nc_diag_chaninfo("chaninfosimple6_str", str_chaninfo)
        
        write(str_metadata, "(A, I0)") "hellometameta_", i
        call nc_diag_metadata("metadatasimple6_str", str_metadata)
        
        write(str_metadata, "(A, I0)") "hellometameta_", i
        call nc_diag_metadata("metadatasimpleX_str", str_metadata)
    end do
    
    do i = 1, 10
        write(str_metadata, "(A, I0)") "morehellometameta_", i
        call nc_diag_metadata("metadatasimple8_str", str_metadata)
        
        write(str_data2d, "(A, I0)") "data2d_more_", i
        call nc_diag_data2d("data2dsimple6_str", (/ str_data2d, "fill1", "fill2" /))
    end do
    
    call nc_diag_write
end program gen_diff_size_netcdf
