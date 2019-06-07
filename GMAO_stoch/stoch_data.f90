module stoch_data

! set up and initialize stochastic random patterns.

 use mod_param, only : nlevs=>levs
 use mod_param, only : r_kind => kind_io4, kind_evod
 use mod_param, only : ens_nam,lonr,latr,jcap
 use mod_param, only : sppt,sppt_tau,sppt_lscale,iseed_sppt,sppt_sigtop1, &
                       sppt_sigtop2,sppt_sfclimit,sppt_logit,stochini,    & 
                       skeb,skeb_tau,skeb_lscale,iseed_skeb,skeb_sigtop1, & 
                       skeb_sigtop2,sppt_sfclimit,sppt_logit,skeb_vfilt,  & 
                       skeb_varspect_opt,skeb_diss_smooth,opt_jff
 use mod_param, only : me, me_l_0,ls_dim, len_trie_ls,len_trio_ls
 use patterngenerator,only: random_pattern, patterngenerator_init,      &
                            getnoise, patterngenerator_advance,ndimspec,&
                            read_pattern,write_pattern,patterngenerator_advance_skeb
 implicit none 
 private
 public :: init_stochdata,destroy_stochdata,init_sppt,init_skeb
 public :: dump_patterns,restore_patterns
 public :: get_vfact
 real(kind_evod),    allocatable, public, dimension(:,:,:)   ::spec_sppt_e,spec_sppt_o
 real(kind_evod),    allocatable, public, dimension(:,:,:,:) ::spec_skeb_e,spec_skeb_o
 real(r_kind),       allocatable, public, dimension(:)       ::vfact_sppt,vfact_skeb
 type(random_pattern),public,save,allocatable,dimension(:)   ::rpattern_sppt,rpattern_skeb
 integer, public :: nsppt=0
 integer, public :: nskeb=0
 logical :: iexist

     namelist/spptparm/sppt,sppt_tau,sppt_lscale,iseed_sppt,sppt_sigtop1,sppt_sigtop2, & 
                       sppt_sfclimit,sppt_logit,stochini
     namelist/skebparm/skeb,skeb_tau,skeb_lscale,iseed_skeb,skeb_sigtop1,skeb_sigtop2, & 
                       skeb_diss_smooth,skeb_vfilt,stochini,skeb_varspect_opt,opt_jff
 contains
 
 subroutine init_sppt

    integer ::  n

    inquire(file="stoch.rc", exist=iexist)

    if (iexist ) then 
       print*,'init_sppt: Reading from namelist'
       open(11, file='stoch.rc')
       read(11,spptparm)
       ! determine number of random patterns to be used for each scheme.
       nsppt = 0
       do n=1,size(sppt)
         if (sppt(n) > 0) then
           nsppt=nsppt+1
         else
           exit
         endif
       enddo
       print*, 'init_sppt : Using ', nsppt ,'patterns'
       close(11)
    else  ! using default 
      nsppt          = 1
      sppt(:)        = -1
      sppt_tau(:)    = -1
      sppt_lscale(:) = -1
      sppt(1)        = 0.7 
      sppt_tau(1)    = 6*3600
      sppt_lscale(1) = 500*1000.
      iseed_sppt     = 4357
      ! parameters to control vertical tapering of stochastic physics with height 
      sppt_sigtop1   = 0.1        ! parameters to control vertical tapering of stochastic physics with height
      sppt_sigtop2   = 0.025      ! //          //           //       
      sppt_sfclimit  = .true.     ! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_logit     = .true.     ! logit transform for sppt to bounded interval [-1,+1]
      stochini       = .false.    ! true= read in pattern, false=initialize from seed
    endif 

 end subroutine init_sppt
 
 subroutine init_skeb


    integer ::  n

    inquire(file="stoch.rc", exist=iexist)

    if (iexist ) then 
       print*,'init_skeb: Reading from namelist'
       open(11, file='stoch.rc')
       read(11,skebparm)
       ! determine number of random patterns to be used for each scheme.
       nskeb = 0
       do n=1,size(skeb)
         if (skeb(n) > 0) then
           nskeb=nskeb+1
         else
           exit
         endif
       enddo
       print*, 'init_skeb : Using ', nskeb ,'patterns'
       close(11)
    else  ! using default 
      nskeb          = 1
      skeb(:)        = -1
      skeb_tau(:)    = -1
      skeb_lscale(:) = -1
      skeb(1)        = 0.28
      skeb_tau(1)    = 6*3600
      iseed_skeb     = 4357
      skeb_sigtop1   = 0.1
      skeb_sigtop2   = 0.025

      stochini          = .false. ! true= read in pattern, false=initialize from seed
      skeb_diss_smooth  = 12.     ! wavenumber defining smoothing of skeb dissipation estimate.
      skeb_varspect_opt = 1       ! gaussian or power law variance spectrum for skeb (0: gaussian, 1:
                                  ! power law). if power law, skeb_lscale interpreted as a power not a
                                  ! length scale.
      skeb_lscale(1)    = -1.27
      skeb_vfilt        = 1       ! number of passes of 1-2-1 vertical filter for skeb random patterns.
      opt_jff           = .false.
    endif 

 end subroutine init_skeb
 
 subroutine init_stochdata(delt,ls_node)

   ! initialize random patterns.  a spinup period of spinup_efolds times the
   ! temporal time scale is run for each pattern.

   real(r_kind) sigtopl,sigtop,sigbot,si(nlevs+1),sl
   real(kind_evod) delt
   integer, intent(in) :: ls_node(ls_dim,3)
   integer nn,nspinup,k,nm,spinup_efolds,stochlun,ierr,iret,n
   stochlun=99

   if (nsppt > 0) allocate(rpattern_sppt(nsppt))
   if (nskeb > 0) allocate(rpattern_skeb(nskeb))

!  if stochini is true, then read in pattern from a file
   if (me .eq. me_l_0) then
      if (stochini) then
         print*,'opening stoch_ini'//trim(ens_nam)
         open(stochlun,file='stoch_ini'//trim(ens_nam),form='unformatted',iostat=ierr,status='old')
         if (ierr .ne. 0) then
            print*,'error opening stoch_ini, error=',ierr
            call glob_abort(ierr,'stoch_date: error in opening file',1)
         endif
      endif
   endif

   ! no spinup needed if initial patterns are defined correctly.
   spinup_efolds = 0
   if (nsppt > 0) then
      allocate(spec_sppt_e(len_trie_ls,2,nsppt))
      allocate(spec_sppt_o(len_trio_ls,2,nsppt))
      spec_sppt_e=0.0
      spec_sppt_o=0.0       
      if (me .eq. me_l_0) print *, 'initialize random pattern for sppt'
      call patterngenerator_init(sppt_lscale,delt,sppt_tau,sppt,iseed_sppt,rpattern_sppt, &
                                 lonr,latr,jcap,ls_node,nsppt,0)
      do n=1,nsppt
         nspinup = spinup_efolds*sppt_tau(n)/delt
         if (stochini) then
            call read_pattern(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),stochlun)
         else
            call getnoise(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n))
            do nn=1,len_trie_ls
               nm = rpattern_sppt(n)%idx_e(nn)
               if (nm .eq. 0) cycle
               spec_sppt_e(nn,1,n) = rpattern_sppt(n)%stdev*spec_sppt_e(nn,1,n)*rpattern_sppt(n)%varspectrum(nm)
               spec_sppt_e(nn,2,n) = rpattern_sppt(n)%stdev*spec_sppt_e(nn,2,n)*rpattern_sppt(n)%varspectrum(nm)
            enddo
            do nn=1,len_trio_ls
               nm = rpattern_sppt(n)%idx_o(nn)
               if (nm .eq. 0) cycle
               spec_sppt_o(nn,1,n) = rpattern_sppt(n)%stdev*spec_sppt_o(nn,1,n)*rpattern_sppt(n)%varspectrum(nm)
               spec_sppt_o(nn,2,n) = rpattern_sppt(n)%stdev*spec_sppt_o(nn,2,n)*rpattern_sppt(n)%varspectrum(nm)
            enddo
            !do nn=1,nspinup
            do nn=1,100
            !  print*,'going through first spinup'
               call patterngenerator_advance(spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),rpattern_sppt(n))
            enddo
         endif
      enddo
   endif

   if (nskeb > 0) then
       ! backscatter noise.
       allocate(spec_skeb_e(len_trie_ls,2,nlevs,nskeb))
       allocate(spec_skeb_o(len_trio_ls,2,nlevs,nskeb))
       spec_skeb_e=0.0
       spec_skeb_o=0.0
       if (me .eq. me_l_0) print *, 'initialize random pattern for skeb'
       call patterngenerator_init(skeb_lscale,delt,skeb_tau,skeb,iseed_skeb,rpattern_skeb, &
           lonr,latr,jcap,ls_node,nskeb,skeb_varspect_opt)
       do n=1,nskeb
       do k=1,nlevs
          nspinup = spinup_efolds*skeb_tau(n)/delt
          if (stochini) then
             call read_pattern(rpattern_skeb(n),spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),stochlun)
          else
             call getnoise(rpattern_skeb(n),spec_skeb_e(1,1,k,n),spec_skeb_o(1,1,k,n))
             do nn=1,len_trie_ls
                nm = rpattern_skeb(n)%idx_e(nn)
                if (nm .eq. 0) cycle
                spec_skeb_e(nn,1,k,n) = rpattern_skeb(n)%stdev*spec_skeb_e(nn,1,k,n)*rpattern_skeb(n)%varspectrum(nm)
                spec_skeb_e(nn,2,k,n) = rpattern_skeb(n)%stdev*spec_skeb_e(nn,2,k,n)*rpattern_skeb(n)%varspectrum(nm)
             enddo
             do nn=1,len_trio_ls
                nm = rpattern_skeb(n)%idx_o(nn)
                if (nm .eq. 0) cycle
                spec_skeb_o(nn,1,k,n) = rpattern_skeb(n)%stdev*spec_skeb_o(nn,1,k,n)*rpattern_skeb(n)%varspectrum(nm)
                spec_skeb_o(nn,2,k,n) = rpattern_skeb(n)%stdev*spec_skeb_o(nn,2,k,n)*rpattern_skeb(n)%varspectrum(nm)
             enddo
             if( opt_jff .eqv. .true.) then 
               do nn=1,nspinup
                  call patterngenerator_advance(spec_skeb_e(1,1,k,n),spec_skeb_o(1,1,k,n),rpattern_skeb(n))
               enddo
             else 
               do nn=1,nspinup
                  call patterngenerator_advance_skeb(spec_skeb_e(1,1,k,n),spec_skeb_o(1,1,k,n),rpattern_skeb(n))
               enddo 
             endif
          endif
       enddo
       enddo
   endif ! skeb > 0

   if (me .eq. me_l_0 .and. stochini) close(stochlun)
 end subroutine init_stochdata

 subroutine get_vfact(PREF)
  
  implicit none
  real*4, intent(in)    :: PREF(nlevs+1)
  real*8  :: P00
  real(r_kind) sigtopl,sigtop,sigbot,si(nlevs+1),sl
  integer :: k,l 
  integer :: lun,ierr

  P00 = 100000.
   do k=1,nlevs+1
       si(k)=PREF(k)/P00
   enddo
   if( nsppt > 0 ) then 
     allocate(vfact_sppt(nlevs))
     do k=1,nlevs
        sl = 0.5*(si(k)+si(k+1))
        if (sl .lt. sppt_sigtop1 .and. sl .gt. sppt_sigtop2) then
           vfact_sppt(k) = (sl-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
        else if (sl .lt. sppt_sigtop2) then
           vfact_sppt(k) = 0.0
        else
           vfact_sppt(k) = 1.0
        endif
     enddo
     if (sppt_sfclimit) then
        ! reduce sppt amplitude in lowest two levels
        ! to prevent instability 
        vfact_sppt(nlevs-1) = vfact_sppt(nlevs-1)*0.5
        vfact_sppt(nlevs)  =  vfact_sppt(nlevs-1)*0.5
!       vfact_sppt(nlevs)  =  0.0
     endif
     print *, 'nlevs check ', nlevs
     prinT*, vfact_sppt
   endif
   if( nskeb > 0 ) then 
     allocate(vfact_skeb(nlevs))
     do k=1,nlevs
        sl = 0.5*(si(k)+si(k+1))
        if (sl .lt. skeb_sigtop1 .and. sl .gt. skeb_sigtop2) then
           vfact_skeb(k) = (sl-skeb_sigtop2)/(skeb_sigtop1-skeb_sigtop2)
        else if (sl .lt. skeb_sigtop2) then
           vfact_skeb(k) = 0.0
        else
           vfact_skeb(k) = 1.0
        endif
     enddo
   endif


 end subroutine get_vfact
 
 subroutine dump_patterns(sfile)
    character*9 :: sfile
    integer :: stochlun,k,n
    stochlun=99
    if (me .eq. me_l_0) then
       if (nsppt > 0 .or. nskeb > 0) then
          open(stochlun,file=sfile//trim(ens_nam),form='unformatted')
          print*,'open ',sfile,' for output'
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt 
       call write_pattern(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),stochlun)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
       do k=1,nlevs
          call write_pattern(rpattern_skeb(n),spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),stochlun)
       enddo
       enddo
    endif

    close(stochlun)
 end subroutine dump_patterns

 subroutine restore_patterns(sfile)
    character*9 :: sfile
    integer :: stochlun,k,n
    stochlun=99
    if (me .eq. me_l_0) then
       if (nsppt > 0 .or. nskeb > 0) then
          open(stochlun,file=sfile//trim(ens_nam),form='unformatted',status='old')
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt
       call read_pattern(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),stochlun)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
       do k=1,nlevs
          call read_pattern(rpattern_skeb(n),spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),stochlun)
       enddo
       enddo
    endif
    close(stochlun)
 end subroutine restore_patterns

 subroutine destroy_stochdata()
    ! deallocate arrays.
    if (allocated(spec_sppt_e))   deallocate(spec_sppt_e,spec_sppt_o)
    if (allocated(rpattern_sppt)) deallocate(rpattern_sppt)
    if (allocated(spec_skeb_e))   deallocate(spec_skeb_e,spec_skeb_o,vfact_skeb)
    if (allocated(rpattern_skeb)) deallocate(rpattern_skeb)
 end subroutine destroy_stochdata

end module stoch_data
