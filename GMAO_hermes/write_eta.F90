program write_eta
  use, intrinsic :: iso_fortran_env, only: REAL32, REAL64
  use m_set_eta
  implicit none
  character(len=128) :: arg
  integer :: nxt
  integer :: levels, use_sigma, use_ncep, stat
  character :: opt
  integer :: ks, k, ks_tmp, status, unit
  real(kind=REAL64), allocatable :: ak(:), bk(:)
  real(kind=REAL64) :: ak_tmp, bk_tmp
  real(kind=REAL64) :: p0, ptop, pint
  character(len=100) :: line
  character(len=100) :: filename
  !
  ! -l levels
  ! -s 1 use sigma
  ! -s 0 do not use sigma
  ! -c 1 use ncep levels
  ! -c 0 do not use ncep levels
  ! -p reference pressure
  ! -f output filename
  character(len=*), parameter :: Usage="Usage: ./write_eta.x -l n  -s 1 (or 0) -e 1 (or 0) -p ref_pressure -f filename \n &
                                        required: -l \n &
                                         default: -s 0 -e 0 -f eta.rc "

  p0 = 98400.d0 
  use_sigma = 0
  use_ncep  = 0
  levels    = 0
  nxt = 1
  filename = 'eta.rc'

  call getarg(nxt,arg)

  do while(arg(1:1)=='-')

     opt=arg(2:2)
     if(len(trim(arg))==2) then
        nxt = nxt + 1
        call getarg(nxt,arg)
     else
        arg = arg(3:)
     end if

     select case (opt)
     case ('l')
        read(arg,*, iostat=stat) levels
     case ('s')
        read(arg,*, iostat=stat) use_sigma
     case ('e')
        read(arg,*, iostat=stat) use_ncep
     case ('p')
        read(arg,*, iostat=stat) p0
     case ('f')
        read(arg,*, iostat=stat) filename
     case default
        stop Usage
     end select
     nxt = nxt + 1
     call getarg(nxt,arg)
  end do

  if (use_sigma /=0 ) then
     call set_sigma()
  endif

  if (use_ncep /=0 ) then
     call set_ncep72()
  endif

  if (levels == 0 ) stop "-l is a must"

  allocate(ak(levels+1),bk(levels+1))

  call set_eta(levels, ks, ptop, pint, ak, bk)

  open(newunit=unit, file=trim(filename), action="write",status='new', FORM="FORMATTED")

  write(unit,'(A)') "# The label bk-ak is for data for ak(k), bk(k)"
  write(unit,'(A)') "# the label TRANSIT_TO_P is ks for pint"
  write(unit,*)
  write(unit,'(A, I5)') "NUM_LEVELS:  ", levels
  write(unit,*)

  write(unit,'(A)') "ak-bk:  "
  do k = 1,levels+1
     write(unit,'(2ES23.15)'), ak(k), bk(k)
  enddo
 
  write(unit,'(A, ES23.15)') "REF_PRESSURE: ",p0
  close(unit) 

  ! read back date from the file and verify the data

  open(newunit=unit, file=trim(filename), action="read", FORM="FORMATTED")
  do 
    read(unit,'(A)', IOSTAT = status) line
    if (status < 0) exit
    if(index(line,"ak-bk") == 0)cycle

    ks_tmp = -1
    do k = 1, levels + 1
      read(unit,'(A)', IOSTAT = status) line
      read(line,*) ak_tmp, bk_tmp

      if(abs(ak_tmp-ak(k)) > 0.0d0) then
         print*, "WARNING!! double check ak(",k,") which is inconsistent"
         print*, ak(k), "is changed to ", ak_tmp
      endif
      if(abs(bk_tmp-bk(k)) > 0.0d0) then
         print*, "WARNING!! double check bk(",k,") which is inconsistent"
         print*, bk(k), "is changed to ", bk_tmp
      endif
      if (bk_tmp > 0.0d0 .and. ks_tmp == -1 ) then
           ks_tmp = k-2
           if (ks_tmp /= ks) then
              print*, "WARNING!! double check ks = ", ks
              print*, "ks is changed to ", ks_tmp
           endif
       endif
     enddo
     exit
   enddo

   close(unit)

end program

