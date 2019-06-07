program err_ncdr_invalid_ncdr_id_dim
    use kinds, only: i_long
    use nc_diag_read, only: nc_diag_read_get_dim
    
    implicit none
    
    integer(i_long) :: i
    
    i = nc_diag_read_get_dim(1234567, "asdf")
end program err_ncdr_invalid_ncdr_id_dim
