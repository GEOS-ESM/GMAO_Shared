program dyn_jediupd

   use m_dyn, only: dyn_init
   use m_dyn, only: dyn_get
   use m_dyn, only: dyn_put
   use m_dyn, only: dyn_vect
   use m_dyn, only: dyn_clean

   use m_maph_pert, only: h_map_pert
   use m_set_eta, only: set_eta

   use m_const, only: zvir
   use m_die, only: die

   use m_nc_JEDIinc, only: nc_JEDIinc_vars
   use m_nc_JEDIinc, only: nc_JEDIinc_dims
   use m_nc_JEDIinc, only: nc_JEDIinc_read

   implicit none

   integer nymd, nhms, freq, dyntype
   integer nlat,nlon,nlev,lm
   integer k,ks,rc
   type(dyn_vect) x_i
   type(dyn_vect) x_x
   type(nc_JEDIinc_vars) y_i
   real(8) :: ptop8, pint8
   real(8),allocatable :: ak8(:),bk8(:)
   logical :: verbose

   character(len=*), parameter :: myname="dyn_jediupd"
   character(len=255) :: files(1)
   character(len=255) :: outinc

   call init_()

!  Get JEDI increment
!  ------------------
   call nc_JEDIinc_dims (files(1),nlat,nlon,nlev,rc)
   call nc_JEDIinc_read (files(1),y_i,rc, gsiset=.false. )

   allocate(ak8(nlev+1),bk8(nlev+1))
   call set_eta ( nlev, ks, ptop8, pint8, ak8, bk8 )
   print *, ak8, bk8

   if(verbose) print *, 'inc: ', nlon,nlat,nlev
   lm = 6 ! wire for now
   call dyn_init ( nlon, nlat, nlev, lm, x_x, rc, &
                   ptop8, ks, ak8, bk8, vectype=dyntype )
   print *, 'sum(ts) = ', sum(y_i%ts)
   x_x%ts = y_i%ts
   x_x%ps = y_i%ps
   x_x%pt = y_i%t
   x_x%u  = y_i%u
   x_x%v  = y_i%v
   x_x%q(:,:,:,1) = y_i%qv
   x_x%q(:,:,:,2) = y_i%oz

   do k=1,x_x%grid%km
      x_x%delp(:,:,k)= (x_x%grid%bk(k+1) - x_x%grid%bk(k))*x_x%ps(:,:)
   enddo

!  call h_map_pert ( x_x, x_i, 'tlm', rc )
!  print *, 'after interp sum(ts) = ', sum(x_i%ts)
!  call dyn_clean  ( x_x )

   call dyn_put ( trim(outinc), nymd, nhms, 0, x_x, rc, vectype=dyntype )
   
!  call dyn_clean (x_i)
   call dyn_clean (x_x)

CONTAINS

subroutine init_()
 implicit none

 integer i, iarg, argc, iargc
 integer ncount,nc
 character(len=255) :: argv

 verbose = .false.
 dyntype = 5
 outinc = 'inc.nc4'

 argc =  iargc()
 if ( argc .lt. 1 ) call usage_()

 iarg=0
 ncount=0; nc=0
 do i = 1, 32767
    iarg = iarg + 1
    if ( iarg .gt. argc ) exit
    call GetArg ( iarg, argv )
    select case (argv)

      case ('-verb')
         verbose=.true.

      case ("-o")
         if ( iarg+1 .gt. argc ) call usage_()
         iarg = iarg + 1
         call GetArg ( iarg, outinc )
      case default
        ncount = ncount + 1
        if (ncount==1) then
           call GetArg ( iarg, argv )
           read(argv,*) nymd
        else if (ncount==2 ) then
           call GetArg ( iarg, argv )
           read(argv,*) nhms
        else 
           nc=nc+1
           if(nc .gt. 1) then
!             call die(myname,'too many cases')
              print *, myname,'too many cases'
              call exit(1)
           endif
           files(nc)  = trim(argv)
        endif

    end select 

 enddo
 if (verbose) then
    write(*,'(a,i8.8,2x,i6.6)') 'Current time: ', nymd, nhms
 endif
 print *, 'Increment  file: ', trim(files(1))

end subroutine init_
subroutine usage_
  print *
  print *, 'Usage: ', trim(myname), '.x [args] nymd nhms inc ' 
  print *
  print *, ' -o FILE  output filename'
  print *
  stop 1
end subroutine usage_
end program dyn_jediupd

