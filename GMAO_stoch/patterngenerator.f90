module patterngenerator

 ! generate random patterns with specified temporal and spatial auto-correlation
 ! in spherical harmonic space.
 use mersenne_twister, only: random_setseed,random_gauss,random_stat
 use mod_param
 use mod_param, only: r_kind => kind_io4, kind_evod
 use mod_param, only:  pi => con_pi, rerth => con_rerth

 implicit none
 private

 public :: computevarspec, computevargrid, setvarspect,&
           patterngenerator_init, patterngenerator_destroy, getnoise, &
           patterngenerator_advance, getvarspectrum, random_pattern, ndimspec, &
           read_pattern,write_pattern,patterngenerator_advance_skeb

 type random_pattern
    real(kind_evod), public :: lengthscale
    real(kind_evod), public :: tau
    real(kind_evod), public :: dt
    real(kind_evod), public :: phi
    real(kind_evod), public :: stdev
    real(4), allocatable, dimension(:), public :: varspectrum, varspectrum1d, lap
    integer, allocatable, dimension(:), public      ::    degree,order,idx_e,idx_o
    integer, allocatable, dimension(:,:), public    :: idx
    integer, public :: seed
    type(random_stat), public :: rstate
 end type random_pattern

 integer :: nlons,nlats,ntrunc,ndimspec
  
 contains

 subroutine patterngenerator_init(lscale, delt, tscale, stdev, iseed, rpattern,&
                                  nlon, nlat, jcap, ls_node, npatterns,&
                                  varspect_opt)
   real(kind_evod), intent(in),dimension(npatterns) :: lscale,tscale,stdev
   real(kind_evod), intent(in) :: delt
   integer, intent(in) :: nlon,nlat,jcap,npatterns,varspect_opt
   integer, intent(in) :: ls_node(ls_dim,3)
   type(random_pattern), intent(out), dimension(npatterns) :: rpattern
   integer(8), intent(in) :: iseed(npatterns)
   integer m,j,l,n,nm,nn,np,indev1,indev2,indod1,indod2
   integer(8) count, count_rate, count_max, count_trunc
   integer(8) :: iscale = 10000000000
   integer count4, member_id, ierr
   integer indlsod,indlsev,jbasev,jbasod
   include 'function_indlsod.h'
   include 'function_indlsev.h'

   nlons = nlon
   nlats = nlat
   ntrunc = jcap

   ndimspec = (ntrunc+1)*(ntrunc+2)/2   
   do np=1,npatterns
      allocate(rpattern(np)%idx(0:ntrunc,0:ntrunc))
      allocate(rpattern(np)%idx_e(len_trie_ls))
      allocate(rpattern(np)%idx_o(len_trio_ls))
      rpattern(np)%idx_e = 0; rpattern(np)%idx_o = 0; rpattern(np)%idx = 0
      nm = 0
      do m=0,ntrunc
         do n=m,ntrunc
            nm = nm + 1
            rpattern(np)%idx(m,n) = nm
         enddo
      enddo
      do j = 1, ls_max_node
         l=ls_node(j,1) ! zonal wavenumber
         jbasev=ls_node(j,2)
         jbasod=ls_node(j,3)
         indev1 = indlsev(l,l)
         indod1 = indlsod(l+1,l)
         if (mod(l,2) .eq. mod(ntrunc+1,2)) then
            indev2 = indlsev(ntrunc+1,l)
            indod2 = indlsod(ntrunc  ,l)
         else
            indev2 = indlsev(ntrunc  ,l)
            indod2 = indlsod(ntrunc+1,l)
         endif
         n = l ! degree
         do nn=indev1,indev2
            if (n <= ntrunc .and. l <= ntrunc) then
              nm = rpattern(np)%idx(l,n)
              rpattern(np)%idx_e(nn) = nm
            endif
            n = n + 2
         enddo
         n = l+1
         do nn=indod1,indod2
            if (n <= ntrunc .and. l <= ntrunc) then
              nm = rpattern(np)%idx(l,n)
              rpattern(np)%idx_o(nn) = nm
            endif
            n = n + 2
         enddo
      enddo
      allocate(rpattern(np)%degree(ndimspec),rpattern(np)%order(ndimspec),rpattern(np)%lap(ndimspec))
      rpattern(np)%degree = (/((n,n=m,ntrunc),m=0,ntrunc)/)
      rpattern(np)%order = (/((m,n=m,ntrunc),m=0,ntrunc)/)
      rpattern(np)%lap = -rpattern(np)%degree*(rpattern(np)%degree+1.0)
      rpattern(np)%tau = tscale(np)
      rpattern(np)%lengthscale = lscale(np)
      rpattern(np)%dt = delt
      rpattern(np)%phi = exp(-delt/tscale(np))
      rpattern(np)%stdev = stdev(np)
      allocate(rpattern(np)%varspectrum(ndimspec))
      allocate(rpattern(np)%varspectrum1d(0:ntrunc))
      ! seed computed on root, then bcast to all tasks and set.
      if (me == me_l_0) then
         member_id=ens_mem
         print*,' member_id is ',member_id
         if (iseed(np) == 0) then
           ! generate a random seed from system clock and ens member number
           call system_clock(count, count_rate, count_max)
           ! iseed is elapsed time since unix epoch began (secs)
           ! truncate to 4 byte integer
           count_trunc = iscale*(count/iscale)
           count4 = count - count_trunc + member_id
           print *,'using seed',count4
         else
           !count4 = iseed(np) + member_id
           ! don't rely on compiler to truncate integer(8) to integer(4) on
           ! overflow, do wrap around explicitly.
           count4 = mod(iseed(np) + member_id + 2147483648, 4294967296) - 2147483648
           print *,'using seed iseed member_id ',count4,iseed(np),member_id
         endif
      endif
      ! broadcast seed to all tasks.
      ! call mpi_bcast(count4,1,mpi_integer,me_l_0,mc_comp,ierr)
      rpattern(np)%seed = count4
      ! set seed (to be the same) on all tasks. save random state.
      call random_setseed(rpattern(np)%seed,rpattern(np)%rstate)
   
      if (varspect_opt .ne. 0 .and. varspect_opt .ne. 1) then
         if (me == me_l_0) then
            print *,'warning: illegal value for varspect_opt (should be 0 or 1), using 0 (gaussian spectrum)...'
         endif
         call setvarspect(rpattern(np),0)
      else if ( (varspect_opt == 1) .and. (opt_jff .eqv. .false.) ) then 
         call setvarspect_skeb(rpattern(np),varspect_opt)
      else
         call setvarspect(rpattern(np),varspect_opt)
      endif
   enddo ! n=1,npatterns
 end subroutine patterngenerator_init

 subroutine patterngenerator_destroy(rpattern,npatterns)
   type(random_pattern), intent(inout) :: rpattern(npatterns)
   integer, intent(in) :: npatterns
   integer n
   do n=1,npatterns
   deallocate(rpattern(n)%varspectrum,rpattern(n)%varspectrum1d)
   deallocate(rpattern(n)%degree,rpattern(n)%order,rpattern(n)%lap)
   deallocate(rpattern(n)%idx,rpattern(n)%idx_e,rpattern(n)%idx_o)
   enddo
 end subroutine patterngenerator_destroy

 subroutine computevarspec(rpattern,dataspec,var)
    ! compute globally integrated variance from spectral coefficients
    complex(r_kind), intent(in) :: dataspec(ndimspec)
    real(r_kind), intent(out) ::  var
    type(random_pattern), intent(in) :: rpattern
    integer n
    var = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
           var = var + dataspec(n)*conjg(dataspec(n))
       else
           var = var + 0.5*dataspec(n)*conjg(dataspec(n))
       endif
    enddo
 end subroutine computevarspec

 subroutine computevargrid(datagrid,var,areawts)
    ! compute globally integrated variance from data on gaussian grid
    real(r_kind), intent(in) :: datagrid(nlons,nlats),areawts(nlons,nlats)
    real(r_kind), intent(out) ::  var
    var = sum(areawts*datagrid**2)
 end subroutine computevargrid

 subroutine getvarspectrum(rpattern,dataspec,varspect)
    type(random_pattern), intent(in) :: rpattern
    complex(r_kind), intent(in) :: dataspec(ndimspec)
    real(r_kind), intent(out) :: varspect(0:ntrunc)
    integer n
    varspect = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
          varspect(rpattern%degree(n)) = varspect(rpattern%degree(n)) + &
          dataspec(n)*conjg(dataspec(n))
       else
          varspect(rpattern%degree(n)) = varspect(rpattern%degree(n)) + &
          0.5*dataspec(n)*conjg(dataspec(n))
       endif
    enddo
 end subroutine getvarspectrum

 subroutine getnoise(rpattern,noise_e,noise_o)
   real(kind_evod), intent(inout) :: noise_e(len_trie_ls,2)
   real(kind_evod), intent(inout) :: noise_o(len_trio_ls,2)
   ! generate white noise with unit variance in spectral space
   type(random_pattern), intent(inout) :: rpattern
   real :: noise(2*ndimspec)
   integer nm,nn
   call random_gauss(noise,rpattern%rstate)
   noise(1) = 0.; noise(ndimspec+1) = 0.
   noise = noise*sqrt(1./ntrunc)
   noise_e = 0.; noise_o = 0.
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      noise_e(nn,1) = noise(nm)/sqrt(2.*rpattern%degree(nm)+1)
      noise_e(nn,2) = noise(ndimspec+nm)/sqrt(2.*rpattern%degree(nm)+1)
      if (rpattern%order(nm) .eq. 0) then
        noise_e(nn,1) = sqrt(2.)*noise_e(nn,1)
        noise_e(nn,2) = 0.
      endif
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      noise_o(nn,1) = noise(nm)/sqrt(2.*rpattern%degree(nm)+1)
      noise_o(nn,2) = noise(ndimspec+nm)/sqrt(2.*rpattern%degree(nm)+1)
      if (rpattern%order(nm) .eq. 0) then
        noise_o(nn,1) = sqrt(2.)*noise_o(nn,1)
        noise_o(nn,2) = 0.
      endif
   enddo
 end subroutine getnoise
 
subroutine read_pattern(rpattern,pattern2d_e,pattern2d_o,lunptn)
   real(kind_evod), intent(inout) :: pattern2d_e(len_trie_ls,2)
   real(kind_evod), intent(inout) :: pattern2d_o(len_trio_ls,2)
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn
   real(kind_evod),allocatable  :: pattern2d(:)
   integer nm,nn,ierr

   allocate(pattern2d(2*ndimspec))

   ! read only on root process, and send to all tasks
   if (me .eq. me_l_0) then
      read(lunptn) pattern2d
      print*,'reading in random pattern (min/max/size)',minval(pattern2d),maxval(pattern2d),size(pattern2d)
   endif
   call mpi_bcast(pattern2d,2*ndimspec,mpi_r_io_r,me_l_0,mc_comp,ierr)
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d_e(nn,1) = pattern2d(nm)
      pattern2d_e(nn,2) = pattern2d(ndimspec+nm)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d_o(nn,1) = pattern2d(nm)
      pattern2d_o(nn,2) = pattern2d(ndimspec+nm)
   enddo
   !print*,'after scatter...',me,maxval(pattern2d_e),maxval(pattern2d_o) &
   ! ,minval(pattern2d_e),minval(pattern2d_o)
   deallocate(pattern2d) 
 end subroutine read_pattern
 
 subroutine write_pattern(rpattern,pattern2d_e,pattern2d_o,lunptn)
   real(kind_evod), intent(in) :: pattern2d_e(len_trie_ls,2)
   real(kind_evod), intent(in) :: pattern2d_o(len_trio_ls,2)
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn
   real(kind_evod), allocatable  :: pattern2d(:),pattern2d_out(:)
   integer nm,nn,ierr

   allocate(pattern2d(2*ndimspec),pattern2d_out(2*ndimspec))
   pattern2d=0.0
   pattern2d_out=0.0
   ! fill in apprpriate pieces of array
   !print*,'before collection...',me,maxval(pattern2d_e),maxval(pattern2d_o) &
   ! ,minval(pattern2d_e),minval(pattern2d_o)
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = pattern2d_e(nn,1)
      pattern2d(ndimspec+nm) = pattern2d_e(nn,2)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = pattern2d_o(nn,1)
      pattern2d(ndimspec+nm) = pattern2d_o(nn,2)
   enddo
   call mpi_reduce(pattern2d,pattern2d_out,2*ndimspec,&
                   mpi_r_io_r,mpi_sum,me_l_0,mc_comp,ierr)
  !  write only on root process
   if (me .eq. me_l_0) then
      print*,'writing out random pattern (min/max/size)',&
      minval(pattern2d_out),maxval(pattern2d_out),size(pattern2d_out)
      write(lunptn) pattern2d_out
   endif
   deallocate(pattern2d,pattern2d_out)
 end subroutine write_pattern
 
 subroutine patterngenerator_advance(dataspec_e,dataspec_o,rpattern)
    ! advance 1st-order autoregressive process with
    ! specified autocorrelation (phi) and variance spectrum (spectrum)
    real(kind_evod), intent(inout) :: dataspec_e(len_trie_ls,2)
    real(kind_evod), intent(inout) :: dataspec_o(len_trio_ls,2)
    real(kind_evod) :: noise_e(len_trie_ls,2)
    real(kind_evod) :: noise_o(len_trio_ls,2)
    type(random_pattern), intent(inout) :: rpattern
    integer j,l,n,nn,nm
    call getnoise(rpattern,noise_e,noise_o)
    do nn=1,len_trie_ls
       nm = rpattern%idx_e(nn)
       if (nm == 0) cycle
       dataspec_e(nn,1) =  rpattern%phi*dataspec_e(nn,1) + &
                           rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_e(nn,1)
       dataspec_e(nn,2) =  rpattern%phi*dataspec_e(nn,2) + &
                           rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_e(nn,2)
    enddo
    do nn=1,len_trio_ls
       nm = rpattern%idx_o(nn)
       if (nm == 0) cycle
       dataspec_o(nn,1) =  rpattern%phi*dataspec_o(nn,1) + &
                           rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_o(nn,1)
       dataspec_o(nn,2) =  rpattern%phi*dataspec_o(nn,2) + &
                           rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_o(nn,2)
    enddo
 end subroutine patterngenerator_advance

   subroutine setvarspect(rpattern,varspect_opt)

      use mod_param, only: rerth => con_rerth
     ! define variance spectrum (isotropic covariance)
     ! normalized to unit global variance
     type(random_pattern), intent(inout) :: rpattern
     integer, intent(in) :: varspect_opt
     complex(r_kind) noise(ndimspec)
     real(r_kind) var
     integer n


     ! 1d variance spectrum (as a function of total wavenumber)
     if (varspect_opt == 0) then ! gaussian
       ! rpattern%lengthscale is interpreted as an efolding length
       ! scale, in meters.
       do n=0,ntrunc
          rpattern%varspectrum1d(n) = exp(-rpattern%lengthscale**2*(float(n)*(float(n)+1.))/(4.*rerth**2))
       enddo
       ! scaling factors for spectral coeffs of white noise pattern with unit variance
       rpattern%varspectrum = sqrt(ntrunc*exp(rpattern%lengthscale**2*rpattern%lap/(4.*rerth**2)))
     else if (varspect_opt == 1) then ! power law
       ! rpattern%lengthscale is interpreted as a power, not a length.
       do n=0,ntrunc
         !   rpattern%varspectrum1d(n) = float(n)**(rpattern%lengthscale)
         rpattern%varspectrum1d(n) = float(n+1)**(rpattern%lengthscale)
       enddo
       ! scaling factors for spectral coeffs of white noise pattern with unit variance
       !rpattern%varspectrum = sqrt(ntrunc*(rpattern%degree)**(rpattern%lengthscale))
       rpattern%varspectrum = sqrt(ntrunc*((rpattern%degree+1)**(rpattern%lengthscale)))
     endif
 
     noise = 0.
     do n=1,ndimspec
       if (rpattern%order(n) .ne. 0.) then
         noise(n) = cmplx(1.,1.)/sqrt(2.*rpattern%degree(n)+1)
       else
         noise(n) = sqrt(2.)/sqrt(2.*rpattern%degree(n)+1.)
       endif
     enddo
     noise(1) = 0 ! no global mean.
     ! make sure global mean variance is 1.
     noise = noise*sqrt(1./ntrunc)
     noise = rpattern%varspectrum*noise
     call computevarspec(rpattern,noise,var)
     rpattern%varspectrum = rpattern%varspectrum/sqrt(var)
     rpattern%varspectrum1d = rpattern%varspectrum1d/var

   end subroutine setvarspect

   !=================================================================================================== 
   subroutine patterngenerator_advance_skeb(dataspec_e,dataspec_o,rpattern)
   
     ! advance 1st-order autoregressive process with
     ! specified autocorrelation (phi) and variance spectrum (spectrum)
     real(kind_evod), intent(inout) :: dataspec_e(len_trie_ls,2)
     real(kind_evod), intent(inout) :: dataspec_o(len_trio_ls,2)
     real(kind_evod) :: noise_e(len_trie_ls,2)
     real(kind_evod) :: noise_o(len_trio_ls,2)
     real(r_kind) :: alpha
     type(random_pattern), intent(inout) :: rpattern
     integer j,l,n,nn,nm
  
     alpha = rpattern%dt/rpattern%tau
     call getnoise(rpattern,noise_e,noise_o)
     do nn=1,len_trie_ls
        nm = rpattern%idx_e(nn)
        if (nm == 0) cycle
        dataspec_e(nn,1) =  (1-alpha)*dataspec_e(nn,1) + sqrt(alpha)*rpattern%varspectrum(nm)*rpattern%stdev*noise_e(nn,1)
        dataspec_e(nn,2) =  (1-alpha)*dataspec_e(nn,2) + sqrt(alpha)*rpattern%varspectrum(nm)*rpattern%stdev*noise_e(nn,2)
     enddo
     do nn=1,len_trio_ls
        nm = rpattern%idx_o(nn)
        if (nm == 0) cycle
        dataspec_o(nn,1) =  (1-alpha)*dataspec_o(nn,1) + sqrt(alpha)*rpattern%varspectrum(nm)*rpattern%stdev*noise_o(nn,1)
        dataspec_o(nn,2) =  (1-alpha)*dataspec_o(nn,2) + sqrt(alpha)*rpattern%varspectrum(nm)*rpattern%stdev*noise_o(nn,2)
     enddo

   end subroutine    ! patterngenerator_advance_skeb

   !=================================================================================================== 

   subroutine setvarspect_skeb(rpattern,varspect_opt)

   ! define variance spectrum (isotropic covariance)
   ! normalized to unit global variance

     type(random_pattern), intent(inout) :: rpattern
     integer, intent(in) :: varspect_opt
     complex(r_kind) noise(ndimspec)
     real(r_kind) var,rgamma, F0, delE,delEp, TOT_Backscat, pi,alpha
     integer n

     pi = atan(1.0)*4.0 
     alpha = rpattern%dt/rpattern%tau
     delE = 1.E-5    ! m2 s-3 (Berner et al 2009) 
                     ! This may depend on the model resolution!
     delEp = delE*rpattern%dt

     rpattern%varspectrum = (rpattern%degree+1)**(rpattern%lengthscale)
     rgamma = 0. 
     do n=1,ndimspec
        rgamma = rgamma +  rpattern%degree(n)*(2.*rpattern%degree(n)+1)*      & 
                           (rpattern%degree(n)+1)**(2*rpattern%lengthscale+1)
     enddo

     F0 = sqrt(4*pi*rerth**2*alpha*delEp/rgamma)/rpattern%stdev

     rpattern%varspectrum   = F0*rpattern%varspectrum

   end subroutine setvarspect_skeb

end module patterngenerator
