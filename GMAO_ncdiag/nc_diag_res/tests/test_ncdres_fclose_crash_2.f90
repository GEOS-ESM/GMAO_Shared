program test_ncdres_fclose_crash_2
    use nc_diag_res, only: nc_diag_load_resource_file, &
        nc_diag_close_resource_file
    
    call nc_diag_load_resource_file("test_invalid.json")
    
    call nc_diag_close_resource_file
    
    ! This should fail, since the file is already closed.
    call nc_diag_close_resource_file
end program test_ncdres_fclose_crash_2
