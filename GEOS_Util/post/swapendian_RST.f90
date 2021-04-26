program endian_convert
  implicit none
  
  real, pointer :: var(:)

  integer :: i, bpos, epos, status
  integer :: rsize
  character(128) :: str, f_in, f_out
  integer*4 :: ftell
  external  :: ftell

! Begin

  if (command_argument_count() /= 2) then
     call get_command_argument(0,str)
     write(*,*) "Usage:",trim(str)," <big_endian_in> <native_out>"
     error stop 2
  end if

  call get_command_argument(1,f_in)
  call get_command_argument(2,f_out)

  open(unit=10, file=trim(f_in),  form='unformatted', convert="big_endian")
  open(unit=20, file=trim(f_out), form='unformatted', convert="native")

  print *,'Converting File: ',trim(f_in)

  bpos=0
  do 
     read(10, end=100, err=200) ! skip to next record
     epos = ftell(10)          ! ending position of file pointer
     backspace(10)
     
     rsize = (epos-bpos)/4-2   ! record size (in 4 byte words; 
                               ! 2 is the number of fortran control words)
     bpos=epos
     allocate(var(rsize), stat=status)
     if (status  /= 0) then
        print *, 'Error: allocation ', rsize, ' failed!'
        error stop 11
     end if

     read (10) var
     write(20) var
     deallocate(var)
  end do

100 continue
  close(10)
  close(20)
  stop

200 print *,'Error reading file ',trim(f_in)
    error stop 11

end
