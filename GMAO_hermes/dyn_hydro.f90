program dyn_hydro
use m_dyn
use m_die, only: die
implicit none

character(len=*),parameter :: myname = 'dyn_hydro'
integer,parameter :: mfiles=1
character(len=256) fnames(mfiles), odyn_file

integer nymd, nhms
integer freq, vectype, nstep, prec
integer ilm,ier
type(dyn_vect)  odyn   ! input non-compliant dyn-vector
type(dyn_vect)  idyn   ! input dyn-vector (typically at same date/time as above)

! initialization
  call init_(ier)
  if(ier<0) call usage_

! read file containing dyn-vector at same date/time as file above
  call dyn_get ( trim(fnames(1)), nymd, nhms, idyn, ier, timidx=1, freq=freq, nstep=nstep, vectype=vectype )

! Initialize dimension of output (interpolated) vector
! ----------------------------------------------------
  ilm=size(idyn%q,4)
  call dyn_init ( idyn, odyn, ier, copy=.true., vectype=vectype, lm=ilm+2 )
  if ( ier/=0 ) then
     call die (myname, ': Error initializing dyn vector(odyn)')
  endif

! write out
  call dyn_put ( odyn_file, nymd, nhms, prec, odyn, ier, freq=freq, nstep=nstep, vectype=vectype )
  if ( ier .ne. 0 ) then
       call dyn_clean ( odyn )
       call die(trim(myname), ': cannot write out ETA file, aborting ...')
  endif
  call dyn_clean ( idyn )

CONTAINS

  subroutine init_ (rc)

  integer, intent(out) :: rc

  character(len=*),parameter :: myname_ = 'init_'
  character(len=256) :: argv
  integer argc,i,iarg,iargc,nfiles
  rc=0
  argc =  iargc()
  if ( argc .lt. 1 ) then
     rc=-1
     return
  endif

! defaults
  vectype = 5 
  prec    = 0
  odyn_file = 'hydro.eta.nc4'

  nfiles=0
  do i = 1, 32767
     iarg = iarg + 1
     if ( iarg .gt. argc ) exit
     call GetArg ( iarg, argv )
     select case (argv)
       case ("-o")
         if ( iarg+1 .gt. argc ) then
             rc=-2
             return
         endif
         iarg = iarg + 1
         call GetArg ( iarg, odyn_file )
       case default
         nfiles = nfiles + 1
         if ( nfiles .gt. mfiles ) call die(myname_,': too many eta files')
         fnames(nfiles) = argv
     end select
  enddo
  if(nfiles/=mfiles) rc=-3

  end subroutine init_

  subroutine usage_()
  print *
  print *, '  -----------------------------------------------------------'
  print *, '  dyn_hydro.x - adds hydrometers to input file (no overwrite)'
  print *, '  -----------------------------------------------------------'
  print *
  print *,'Usage: '
  print *,'  dyn_hydro.x [-o odyn_file ] dyn_file'
  print *
  print *,' where '
  print *,'     dyn_file  original file'
  print *,' -o odyn_file  original w/ added hydrometeors'
  print *

  stop
  end subroutine usage_

end program  dyn_hydro
