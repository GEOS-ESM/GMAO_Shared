Program binarytile

  implicit none

  integer, parameter :: unitR=7
  integer, parameter :: unitW=8
  integer, parameter :: NumGlobalVars=4
  integer, parameter :: NumGridVars=3
  integer            :: N
  integer            :: NT, NPFAF, line_size
  integer            :: IM
  integer            :: JM
  integer            :: N_GRIDS,pfaf_number
  integer            :: grid_info(2),status
  real, allocatable  :: AVR(:,:)
  real               :: DUMMY
  character(len=128) :: NAME 
  character(len=128) :: filenameIN 
  character(len=128) :: filenameOUT
  integer, parameter :: max_rec=2 

  call getarg(1,filenameIN)
  if (filenameIN == "") filenameIN = 'input'
  call getarg(2,filenameOUT)
  if (filenameOUT == "") filenameOUT = 'output'

  open(unit=unitR, file=filenameIN, form='FORMATTED')
  open(unit=unitW, file=filenameOUT,form='UNFORMATTED')
  !READ (unitR, *) NT, NPFAF
  !WRITE(unitW   ) NT, NPFAF

! Number of grids that can be attached
!-------------------------------------

  grid_info(2)=-1
  do n=1,max_rec
     rewind(unitR)
     read(unitR,*,iostat=status)grid_info(1:n)
     if (status<0) exit
  enddo
  nt=grid_info(1)
  !pfaf_number=grid_info(2)
  WRITE (unitW) grid_info

! Number of grids that can be attached
!-------------------------------------

  READ (unitR, *) N_GRIDS
  WRITE(unitW   ) N_GRIDS

! The names and sizes of the grids to be tiled
!---------------------------------------------

  do N=1,N_GRIDS
     READ (unitR, *) NAME
     WRITE(unitW   ) NAME
     READ (unitR, *) IM
     WRITE(unitW   ) IM
     READ (unitR, *) JM
     WRITE(unitW   ) JM
  enddo


! Read location stream file into AVR
!---------------------------------------

  allocate(AVR(NT,NumGlobalVars+NumGridVars*N_GRIDS))

  do N=1, NT
     READ(unitR, *) AVR(N,1),AVR(N,10),AVR(N,2:6),DUMMY,AVR(N,7:9)
  end do

  close(unitR)

  do N=1,size(AVR,2)
     write(unitW) AVR(:,N)
  end do

!  write(unitW) AVR(:,1:3)
!  write(unitW) AVR(:,4:6)
!  write(unitW) AVR(:,7:9)
  close(unitW)

end Program binarytile
