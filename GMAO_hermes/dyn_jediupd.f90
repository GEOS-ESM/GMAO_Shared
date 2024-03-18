program dyn_jediupd

   use m_dyn, only: dyn_init
   use m_dyn, only: dyn_get
   use m_dyn, only: dyn_put
   use m_dyn, only: dyn_vect
   use m_dyn, only: dyn_clean

   use m_maph_pert, only: h_map_pert

   use m_const, only: zvir
   use m_die, only: die

   use m_nc_JEDIinc, only: nc_JEDIinc_vars
   use m_nc_JEDIinc, only: nc_JEDIinc_dims
   use m_nc_JEDIinc, only: nc_JEDIinc_read

   implicit none

   integer nymd, nhms, freq, dyntype
   integer nlat,nlon,nlev
   integer k,rc
   type(dyn_vect) x_a
   type(dyn_vect) x_i
   type(dyn_vect) x_x
   type(nc_JEDIinc_vars) y_i
   real, allocatable :: t(:,:,:)
   logical :: verbose

   character(len=*), parameter :: myname="dyn_jediupd"
   character(len=255) :: files(3)
   character(len=255) :: outinc

   call init_()

!  Read backgroud fields
!  ---------------------
   call dyn_get ( files(1), nymd, nhms, x_a, rc, timidx=1, freq=freq, vectype=dyntype )
   if ( rc .ne. 0 ) then
!     call die(myname,'cannot read background file')
      print *, myname,'cannot read background file'
      call exit(1)
   else
      print *, trim(myname), ': read background ', trim(files(1))
   end if

!  Get JEDI increment
!  ------------------
   call nc_JEDIinc_dims (files(2),nlat,nlon,nlev,rc)
   call nc_JEDIinc_read (files(2),y_i,rc, gsiset=.false. )

   if(x_a%grid%im/=nlon .or. x_a%grid%jm/=nlat .or. x_a%grid%km/=nlev) then

     if(verbose) print *, 'bkg: ', x_a%grid%im,x_a%grid%jm,x_a%grid%km
     call dyn_init ( x_a%grid%im, x_a%grid%jm, x_a%grid%km, x_a%grid%lm, x_i, rc, &
                     x_a%grid%ptop, x_a%grid%ks, x_a%grid%ak, &
                     x_a%grid%bk, vectype=dyntype )

     if(verbose) print *, 'inc: ', nlon,nlat,nlev
     call dyn_init ( nlon, nlat, x_a%grid%km, x_a%grid%lm, x_x, rc, &
                     x_a%grid%ptop, x_a%grid%ks, x_a%grid%ak, &
                     x_a%grid%bk, vectype=dyntype )
     if ( size(y_i%t,1)/=nlon .or. size(y_i%t,2)/=nlat ) then
!       call die(myname,'field transpose in unexpected way')
        print *, myname,'field transpose in unexpected way'
        call exit(1)
     endif
     print *, 'sum(ts) = ', sum(y_i%ts)
     x_x%ts = y_i%ts
     x_x%ps = y_i%ps
     x_x%pt = y_i%t
     x_x%u  = y_i%u
     x_x%v  = y_i%v
     x_x%q(:,:,:,1) = y_i%qv
     x_x%q(:,:,:,2) = 0.0*y_i%oz

     do k=1,x_x%grid%km
        x_x%delp(:,:,k)= (x_x%grid%bk(k+1) - x_x%grid%bk(k))*x_x%ps(:,:)
     enddo

     call h_map_pert ( x_x, x_i, 'tlm', rc )
     print *, 'after interp sum(ts) = ', sum(x_i%ts)
     call dyn_clean  ( x_x )
     if (trim(outinc) /= 'NULL') then
        call dyn_put ( trim(outinc), nymd, nhms, 0, x_i, rc, freq=freq, vectype=dyntype )
     endif
   else
!    call die(myname, 'not expecting same resolution bkg/inc')
     print *, myname, 'not expecting same resolution bkg/inc'
     call exit(1)
   endif
   
!  Get background temperature
!  --------------------------
   allocate(t(size(x_a%pt,1),size(x_a%pt,2),size(x_a%pt,3)))
   t = x_a%pt/(1.0+zvir*x_a%q(:,:,:,1))

!  Add increment to background
!  ---------------------------
   x_a%ts = x_a%ts + x_i%ts
   x_a%ps = x_a%ps + x_i%ps
   do k=1,x_a%grid%km
      x_a%delp(:,:,k)= (x_a%grid%ak(k+1) - x_a%grid%ak(k)) + &
                       (x_a%grid%bk(k+1) - x_a%grid%bk(k))*x_a%ps(:,:)
   enddo
   if (any(x_a%delp<0.0)) then
!    call die(myname, 'unacceptable delp')
     print *, myname, 'unacceptable delp'
     call exit(1)
   endif
   t = t + x_i%pt
   x_a%u  = x_a%u + x_i%u
   x_a%v  = x_a%v + x_i%v
   x_a%q(:,:,:,1) = x_a%q(:,:,:,1) + x_i%q(:,:,:,1)
   x_a%q(:,:,:,2) = x_a%q(:,:,:,2) + x_i%q(:,:,:,2)

   x_a%pt = t*(1.0+zvir*x_a%q(:,:,:,1))
   deallocate(t)

!  Write out analysis
!  ------------------
   call dyn_put ( trim(files(3)), nymd, nhms, 0, x_a, rc, freq=freq, vectype=dyntype )
   if ( rc .ne. 0 ) then
!     call die(myname,'cannot write analysis file')
      print *, myname,'cannot write analysis file'
      call exit(1)
   else
      print *, trim(myname), ': wrote analysis ', trim(files(3))
   end if

   call dyn_clean (x_i)
   call dyn_clean (x_a)

CONTAINS

subroutine init_()
 implicit none

 integer i, iarg, argc, iargc
 integer ncount,nc
 character(len=255) :: argv

 verbose = .false.
 dyntype = 5
 outinc = 'NULL'

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
           if(nc .gt. 3) then
!             call die(myname,'too many cases')
              print *, myname,'too many cases'
              call exit(1)
           endif
           files(nc)  = trim(argv)
        endif

    end select 

 enddo
 if (nc/=3) then
    
 endif
 if (verbose) then
    write(*,'(a,i8.8,2x,i6.6)') 'Current time: ', nymd, nhms
 endif
 print *, 'Background file: ', trim(files(1))
 print *, 'Increment  file: ', trim(files(2))
 print *, 'Analysis   file: ', trim(files(3))

end subroutine init_
subroutine usage_
  print *
  print *, 'Usage: ', trim(myname), '.x nymd nhms bkg inc ana ' 
  print *
  stop 1
end subroutine usage_
end program dyn_jediupd

