program test_ncdres_fopen_crash
    use nc_diag_res, only: nc_diag_load_resource_file
    
    call nc_diag_load_resource_file("test_invalid.json")
    
    ! This should fail
    call nc_diag_load_resource_file("test_invalid.json")
end program test_ncdres_fopen_crash
