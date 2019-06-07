program test_ncdres_fclose_crash_1
    use nc_diag_res, only: nc_diag_load_resource_file, &
        nc_diag_close_resource_file
    
    ! Nothing is open, so this should fail!
    call nc_diag_close_resource_file
end program test_ncdres_fclose_crash_1
