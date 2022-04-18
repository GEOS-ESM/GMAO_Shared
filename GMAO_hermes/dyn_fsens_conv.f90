! read reference state and fsens and convert fsens vars
! to JEDI analysis variables

    program dyn_fsens_conv 
    use m_dyn
    use m_ioutil, only : luavail
    use m_die, only: die
    use m_const, only : zvir 
    implicit none 
    character(len=*), parameter :: myname='dyn_fsen_conv'
    integer,parameter :: dyntype=5
    integer,parameter :: nfiles=2
    integer nymd, nhms, lu, n, freq, vectype, prec, ier, nstep
    integer i, nf, iarg, argc, ndim2, ndim3, intarg, iargc
    character(len=255) :: ofile
    character(len=255) :: dynfile(nfiles)
    character(len=255) argv
    type(dyn_vect) w_1
    type(dyn_vect) w_2
    real,allocatable :: t(:,:,:), t_tl(:,:,:), t_ad(:,:,:)
    logical use_ps
    real :: total, dotp(5)
    integer im,jm,km,lm,k,rc
    integer im1,jm1,km1,lm1
    integer im2,jm2,km2,lm2
    logical verbose

    ofile="NONE"
    prec = 0
    vectype = 5
    verbose = .false.

    iarg=0
    argc = iargc()
    if ( argc < 1 ) call usage_()
    nf=0
    do i = 1, 32767
       iarg = iarg + 1
       if ( iarg .gt. argc ) exit
       call GetArg ( iarg, argv )
       select case (trim(argv))
          case ("-verbose")
             verbose = .true.
          case ("-o")
             if ( iarg+1 .gt. argc ) call usage_()
             iarg = iarg + 1
             call GetArg ( iarg, ofile )
          case default
             nf = nf + 1
             if ( nf .gt. nfiles ) call die(myname,'too many eta files:',nf)
             dynfile(nf) = trim(argv)
       end select
    enddo

    if (verbose) then
       do nf=1,nfiles
          write(6,'(a,i3,2a)') "Dyn Vector ", nf , ": ", trim(dynfile(nf))
       enddo
    endif

    if (nf==0) call die(myname,'missing input files')

    call dyn_getdim ( trim(dynfile(1)), im1, jm1, km1, lm1, rc )
    call dyn_getdim ( trim(dynfile(2)), im2, jm2, km2, lm2, rc )
    if(km1/=km2) then ! ignore diff in lm for now
       print *, trim(myname), ': km/lm file 1 = ', km1,lm1
       print *, trim(myname), ': km/lm file 2 = ', km2,lm2
       call die(myname,'inconsistent km/lm')
    endif
    im=min(im1,im2)
    jm=min(jm1,jm2)
    km=min(km1,km2)
    lm=min(lm1,lm2)
    ndim2 = im*jm

!   check dimensions and take proper action
    n=1
    if ( im1==im2 .and. jm1==jm2 ) then
       call dyn_get ( trim(dynfile(1)), nymd, nhms, w_1, ier, timidx=n, &
                                      freq=freq, nstep=nstep, vectype=vectype )
       call dyn_get ( trim(dynfile(2)), nymd, nhms, w_2, ier, timidx=n, &
                                      freq=freq, nstep=nstep, vectype=vectype )
    else
       call die(myname,'cannot handle inconsistent dims')
    endif
    ndim3 = ndim2 * km

    allocate (t_ad(im,jm,km))
    call t2tv_ad (w_1%pt,w_1%q(:,:,:,1),w_2%pt,t_ad)
    w_2%pt = t_ad  ! overwrite field
    deallocate (t_ad)

    if (trim(ofile) == "NONE" ) then
       print *, "Overwriting input fsens with converted fields ..."
       ofile = trim(dynfile(2))
    else
       print *, "Writing converted fsens file ..."
    endif

!   change sign of sensitivity to adjust what for JEDI
!   temporary
!   --------------------------------------------------
!   w_2%delp = -w_2%delp 
!   w_2%ps = -w_2%ps 
!   w_2%pt = -w_2%pt 
!   w_2%u  = -w_2%u 
!   w_2%v  = -w_2%v 
!   w_2%q  = -w_2%q 
    call dyn_put ( trim(ofile), nymd, nhms, 0, w_2, rc, freq=freq, vectype=dyntype )

    contains

    subroutine usage_
       print *, "Purpose: Convert variables in fsens to require JEDI fields"
       print *, "   "
       print *, "Usage: dyn_fsens_conv.x [options] reference_state fsens "
       print *, "   "
       print *, " Optional arguments:  "
       print *, "   "
       print *, " -o         Output converted file (DEFAULT: NONE: overwrite input)"
       print *, " -verbose   Echo actions"
       print *, "   "
       stop
    end subroutine usage_

    subroutine tv2t (tv,qv,t)

    real,intent(in) :: tv(:,:,:), qv(:,:,:)
    real,intent(out) :: t(:,:,:)

    t=tv/(1.d0+zvir*qv)

    end subroutine tv2t

    subroutine tv2t_tl (qv,t, tv_tl,qv_tl, t_tl)

    real,intent(in) :: t    (:,:,:), qv  (:,:,:)
    real,intent(in) :: tv_tl(:,:,:), qv_tl(:,:,:)
    real,intent(out) :: t_tl(:,:,:)

    t_tl = (tv_tl - zvir*t*qv_tl)/(1.d0+zvir*qv)

    end subroutine tv2t_tl

!   subroutine tv2t_ad (qv,t, t_tl, tv_ad, qv_ad)

!   real,intent(in) :: t    (:,:,:), qv  (:,:,:)
!   real,intent(in) :: tv_ad(:,:,:), qv_tl(:,:,:)
!   real,intent(out) :: t_ad(:,:,:)

!   tv_ad = t_tl*(1.d0+zvir*qv)
!   qv_ad = t_tl*(1.d0+zvir*qv)/(zvir*t)

!   end subroutine tv2t_ad

    subroutine t2tv (t,qv,tv)

    real,intent(in) :: t(:,:,:), qv(:,:,:)
    real,intent(out) :: tv(:,:,:)

    tv=t*(1.d0+zvir*qv)

    end subroutine t2tv

    subroutine t2tv_tl (t,qv,t_tl,qv_tl,tv_tl)

    real,intent(in) :: t   (:,:,:), qv   (:,:,:)
    real,intent(in) :: t_tl(:,:,:), qv_tl(:,:,:)
    real,intent(out) :: tv_tl(:,:,:)

    tv_tl=t_tl*(1.d0+zvir*qv)+zvir*t*qv_tl

    end subroutine t2tv_tl

    subroutine t2tv_ad (t,qv,tv_ad,t_ad)

    real,intent(in) :: t    (:,:,:), qv(:,:,:)
    real,intent(in) :: tv_ad(:,:,:)
    real,intent(out) :: t_ad(:,:,:)

    t_ad=tv_ad/(1.d0+zvir*qv)
!   q_ad=tv_tl/(zvir*t)

    end subroutine t2tv_ad

    end program dyn_fsens_conv
