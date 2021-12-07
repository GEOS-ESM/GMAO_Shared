   program dyn_vordiv
!  must run as mpirun -np 1
   use m_ggGradientSP,only : ggGradientSP
   use m_ggGradientSP,only : ggGradientSP_init,clean
   use m_ggGradientSP,only : ggDivo

   use m_compact_diff, only: init_compact_diffs
   use m_compact_diff, only: create_cdiff_coefs
   use m_compact_diff, only: inisph
   use m_compact_diff, only: uv2vordiv

   use m_swapij, only: swapij

   use m_specgrid, only: vordiv2sfvp

   use m_dyn, only: dyn_vect
   use m_dyn, only: dyn_get
   use m_dyn, only: dyn_put
   use m_dyn, only: dyn_stat
   use m_dyn, only: dyn_flip

   use m_die, only: die

   implicit none

   character(len=*), parameter :: myname = 'dyn_vordiv'
   real, parameter :: rearth = 6.3712e+6 

   type(dyn_vect) :: dyn
   type(ggGradientSP) :: gr
   real,dimension(:,:,:), pointer :: vor
   real,dimension(:,:,:), pointer :: div
   real,dimension(:,:,:), pointer :: uw
   real,dimension(:,:,:), pointer :: vw
   real,dimension(:,:,:), pointer :: vp
   real,dimension(:,:,:), pointer :: sf
   real,dimension(:,:,:), pointer :: vx
   real,dimension(:,:,:), pointer :: sx

   real,allocatable,dimension(:)   :: glons
   real,allocatable,dimension(:)   :: glats
   real,allocatable,dimension(:,:,:) :: us,vs

   character(len=255) :: ifile,ofile

   character(len=256)  argv
   character(len=3),pointer::xtrnames(:)
   integer, parameter :: mfiles=2
   integer iargc,nfiles
   integer iarg,argc,nymd,nhms,freq
   integer im,jm,km,jcap
   integer i,ierr

   logical sfvp
   logical overwrite
 
   argc = iargc()
   if ( argc < 1 ) then
        call usage_()
   end if
   
   nfiles=0
   overwrite = .true.
   sfvp=.false.
   do i = 1, 32767
      iarg = iarg + 1
      if ( iarg .gt. argc ) exit
      call GetArg ( iarg, argv )
      select case (argv)
          case ("-noverw")
            overwrite = .false.
      case ("-jcap")
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iarg, argv )
          read(argv,*) jcap
          sfvp=.true.
      case default
          nfiles = nfiles + 1
          if ( nfiles .gt. mfiles ) call die(myname,'too many eta files')
          if (nfiles == 1) ifile = trim(argv)
          if (nfiles == 2) ofile = trim(argv)
      end select


   enddo

   print*, "Reading file: ", trim(ifile)
   allocate(xtrnames(2))
   xtrnames=(/'rh ','mrh'/)
   call dyn_get ( trim(ifile), nymd, nhms, dyn, ierr, timidx=1, freq=freq, vectype=5, xtrnames=xtrnames )
   deallocate(xtrnames)

   im = dyn%grid%im
   jm = dyn%grid%jm
   km = dyn%grid%km
   allocate(vor(im,jm,km),div(im,jm,km))
   allocate( sf(im,jm,km), vp(im,jm,km))
   vor=0.d0; div=0.d0
    sf=0.d0;  vp=0.d0
   print*, " BE AWARE: What you see in output is NOT what you get "
   if (sfvp) then
      print *, "Using compact diff to calculate div/vor"
      call dyn_flip(dyn,dover=.false.)
      allocate(glats(jm))
      allocate(glons(im))
      allocate(us(jm,im,km),vs(jm,im,km))
      allocate(sx(jm,im,km),vx(jm,im,km))
      call get_latlon_()
      call init_compact_diffs(jm,im)
      call create_cdiff_coefs(im,glons)
      call inisph(rearth,glats(2),im,jm-2)
      call swapij(dyn%u,us)
      call swapij(dyn%v,vs)
      call uv2vordiv(us,vs)
      call swapij(us,vor)
      call swapij(vs,div)
      call vordiv2sfvp(sx,vx,us,vs,jcap,.true.)
      call swapij(sx,sf)
      call swapij(vx,vp)
      deallocate(sx,vx)
      deallocate(us,vs)
      deallocate(glats)
      deallocate(glons)

      if(overwrite) then
        dyn%u=vp
        dyn%v=sf
        print*, " u/v fields now contain vp/sf respectively"
      else
        dyn%q(:,:,:,1)=div
        dyn%q(:,:,:,2)=vor
        dyn%q(:,:,:,3)=vp
        dyn%q(:,:,:,4)=sf
        print*, " q1/q2 fields now contain div/vor respectively"
        print*, " q3/q4 fields now contain vp/sf   respectively"
      endif
      call dyn_flip(dyn,dover=.false.)
   else
      print *, "Using direct diff to calculate div/vor"
      call ggGradientSP_init(gr,dyn%grid%im,dyn%grid%lat)
      call ggDivo(gr, &
           dyn%u,     &       ! in: u
           dyn%v,     &       ! in: v
           div,       &       ! div(u,v)
           vor)               ! vor(u,v)
      call clean(gr)
      if(overwrite) then
        dyn%u=div
        dyn%v=vor
        print*, " u/v fields now contain div/vor respectively"
      else
        dyn%q(:,:,:,1)=div
        dyn%q(:,:,:,2)=vor
        print*, " q1/q2 fields now contain div/vor respectively"
      endif
   endif


! write out vor/div and whatever else
! -----------------------------------
   print*, "Wrting file: ", trim(ofile)
   call dyn_put ( trim(ofile), nymd, nhms, 0, dyn, ierr, freq=freq, vectype=5, skip_setvec=.true. )

   deallocate( sf, vp, stat=ierr)
   deallocate(div,vor, stat=ierr)

contains
   subroutine get_latlon_
      implicit none
      real,parameter :: pi = 3.141593e+0 ! as in GSI
      real dlon,pih,deg2rad
      real,allocatable :: dlats(:)
      real,allocatable :: rlats(:)
      integer i
      pih=pi/2.0 
      deg2rad=pi/180.0
      allocate(rlats(jm))
      dlon=2*pi/float(im)
      do i=1,im
        glons(i)=float(i-1)*dlon
      end do
      do i=1, jm
        rlats(i)=-pih + (i-1)*pi / float(jm-1)
      end do
      allocate(dlats(jm))
      call gauss_lat_nmc(dlats(2:jm-1),jm-2)
      glats(2:jm-1) = dlats(2:jm-1) * deg2rad
      glats(1)=-pih
      glats(jm)=pih
      deallocate(dlats)
   end subroutine get_latlon_

   subroutine usage_
        print *
        print *, "Purpose: calculate vorticity and divergence from dyn-type file."
        print *, "         Ouput file is dyn-like file with: "
        print *, "            - u slot replaced with divergence"
        print *, "            - v slot replaced with  vorticity"
        print *, "         or                                  "
        print *, "            - u slot replaced with  vel. potential"
        print *, "            - v slot replaced with stream function"
        print *, "         When jcap present: "
        print *, "            - q(1) slot replaced with divergence"
        print *, "            - q(2) slot replaced with vorticity"
        print *, "            - q(3) slot replaced with velocity potential"
        print *, "            - q(4) slot replaced with stream function"
        print *, "         Note: poles are known not to be treated correctly "
        print *, "               in this case."

        print *
        print *, "Usage: dyn_divo.x [opts] ifile ofile"
        print *
        print *, "  ifile   - input filename "
        print *, "  ofile   - output filename with div/vor"
        print *
        print *, " Options: "
        print *, "  -jcap  JCAP  - set truncation (triggers calc vp/sf)"
        print *, "  -noverw      - do not overwrite u/v w/ vp/sf"
        print *
        stop
   end subroutine usage_

   end program dyn_vordiv

