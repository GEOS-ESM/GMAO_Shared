   program dyn_vordiv
!  must run as mpirun -np 1
   use m_ggGradientSP,only : ggGradientSP
   use m_ggGradientSP,only : ggGradientSP_init,clean
   use m_ggGradientSP,only : ggDivo

   use m_dyn, only: dyn_vect
   use m_dyn, only: dyn_get
   use m_dyn, only: dyn_put
   use m_dyn, only: dyn_stat

   use m_die, only: die

   implicit none

   character(len=*), parameter :: myname = 'dyn_vordiv'

   type(dyn_vect) :: dyn
   type(ggGradientSP) :: gr
   real,dimension(:,:,:), pointer :: vor
   real,dimension(:,:,:), pointer :: div
   real,dimension(:,:,:), pointer :: uw
   real,dimension(:,:,:), pointer :: vw
   real,dimension(:,:,:), pointer :: vp
   real,dimension(:,:,:), pointer :: sf

   character(len=255) :: ifile,ofile

   character(len=256)  argv
   integer iargc
   integer iarg,argc,nymd,nhms,freq
   integer im,jm,km,jcap
   integer ierr
   integer   IROMB
   integer   MAXWV
   integer   IDRTI
   integer   IMAXI
   integer   JMAXI
   integer   IDRTO
   integer   IMAXO
   integer   JMAXO
   integer   KMAX

   logical sfvp
 
   argc = iargc()
   if ( argc < 1 ) then
        print *
        print *, "Purpose: calculate vorticity and divergence from dyn-type file."
        print *, "         Ouput file is dyn-like file with: "
        print *, "            - u slot replaced with divergence"
        print *, "            - v slot replaced with  vorticity"
        print *, "         When jcap present: "
        print *, "            - u slot replaced with jcap-truncated u"
        print *, "            - v slot replaced with jcap-truncated v"
        print *, "            - q(1) slot replaced with divergence"
        print *, "            - q(2) slot replaced with vorticity"
        print *, "            - q(3) slot replaced with velocity potential"
        print *, "            - q(4) slot replaced with stream function"
        print *, "         Note: poles are known not to be treated correctly "
        print *, "               in this case."

        print *
        print *, "Usage: dyn_divo.x ifile ofile [jcap]"
        print *
        print *, "  ifile - input filename "
        print *, "  ofile - output filename with div/vor"
        print *, "  jcap  - optional, when present calc vp/sf"
        print *
        stop
   end if
   
   sfvp=.false.
   iarg = 1
   call GetArg ( iarg, ifile )
   iarg = iarg + 1
   call GetArg ( iarg, ofile )
   if ( argc > 2 ) then
       iarg = iarg + 1
       call GetArg ( iarg, argv )
       read(argv, * ) jcap
       print *, 'Using JCAP = ', jcap
       sfvp=.true.
   endif
      
   print*, "Reading file: ", trim(ifile)
   call dyn_get ( trim(ifile), nymd, nhms, dyn, ierr, timidx=1, freq=freq, vectype=5 )

   im = dyn%grid%im
   jm = dyn%grid%jm
   km = dyn%grid%km
   allocate(vor (im,jm,km), &
            div (im,jm,km), &
            stat=ierr)
   vor = 0.d0
   div = 0.d0

   call ggGradientSP_init(gr,dyn%grid%im,dyn%grid%lat)

   call ggDivo(gr, &
        dyn%u,     &       ! in: u
        dyn%v,     &       ! in: v
        div,       &       ! div(u,v)
        vor)               ! vor(u,v)

   call clean(gr)

   if (sfvp) then
      IROMB=0
      MAXWV=jcap
      IDRTI=0
      IMAXI=im
      JMAXI=jm
      IDRTO=0
      IMAXO=im
      JMAXO=jm
      KMAX= km
      allocate(uw(im,jm,km),vw(im,jm,km))
      allocate(sf(im,jm,km),vp(im,jm,km))
      call SPTRUNV (IROMB,MAXWV,IDRTI,IMAXI,JMAXI, &
                    IDRTO,IMAXO,JMAXO,KMAX, &
                    0,0,0,0, &
                    0,0,0,0,dyn%u,dyn%v, &
                    .true.,uw,vw,.true.,div,vor, &
                    .true.,vp,sf)

      dyn%u=uw
      dyn%v=vw
      dyn%q(:,:,:,1)=div
      dyn%q(:,:,:,2)=vor
      dyn%q(:,:,:,3)=vp
      dyn%q(:,:,:,4)=sf
      deallocate(sf,vp)
      deallocate(uw,vw)
   else
      dyn%u = div
      dyn%v = vor
   endif

! write out vor/div and whatever else
! -----------------------------------
   print*, "Wrting file: ", trim(ofile)
   call dyn_put ( trim(ofile), nymd, nhms, 0, dyn, ierr, freq=freq, vectype=5 )
   print*, " BE AWARE: What you see is NOT what you get "
   if(sfvp)then
      print*, " u/v   fields now contain JCAP-truncated u/v respectively"
      print*, " q1/q2 fields now contain div/vor respectively"
      print*, " q3/q4 fields now contain vp/sf   respectively"
   else
      print*, " u/v fields now contain div/vor respectively"
   endif

   deallocate(div,vor, stat=ierr)

   end program dyn_vordiv

