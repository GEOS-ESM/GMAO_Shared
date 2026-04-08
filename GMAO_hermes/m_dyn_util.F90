      module m_dyn_util

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_dyn_util: utilities for dyn-vect
!
! !USES:
!
      use m_dyn

      implicit NONE
      private

! !DESCRIPTION: Collection of utilities related to dyn-vect
!
! !REVISION HISTORY:
!
!  30Aug2017   Todling    Initial wrap from collected utils
!  10Jun2020   Todling    Add two versions of qsat
!
!-------------------------------------------------------------------------
!EOP

   public :: dyn_util_tv2t
   public :: Dyn_Scale_by_TotEne
   public :: Dyn_Qsat
   public :: Dyn_Get_Energy

   interface Dyn_Util_Tv2T
      module procedure tv2t_
   end interface
   interface Dyn_Scale_by_TotEne
      module procedure scale_by_totene_
   end interface
   interface Dyn_TotEne_Weights
      module procedure TotEne_Weights_
   end interface
   interface Dyn_TotEne_Dotp
      module procedure dotp_
   end interface
   interface Dyn_Get_Energy
      module procedure get_ene_field_
   end interface
   interface Dyn_Qsat
      module procedure becov_qsat_
      module procedure gsi_qsat_
   end interface

   integer, parameter :: nkxe = 1   ! index for kinetic energy contribution
   integer, parameter :: nape = 2   ! index for avail. pot. energy contribution
   integer, parameter :: ngpe = 3   ! index for geopotential. energy contribution
   integer, parameter :: nqxe = 4   ! index for moist energy contribution
   integer, parameter :: ntxe = 5   ! index for total dry energy
   integer, parameter :: ntwe = 6   ! index for total wet energy
   integer, parameter :: nall = 6   ! max index from those above

   character(len=*), dimension(nall), parameter :: vnames = &
                      (/'kxe','ape','pse','qxe','txe','twe'/)
   character(len=*), dimension(nall), parameter :: vunits = &
                      (/'J/Kg','J/Kg','J/Kg','J/Kg','J/Kg','J/Kg'/)
   character(len=*), dimension(nall), parameter :: vtitles = &
                       (/'Kinetic Energy            ', &
                         'Available Potential Energy', &
                         'Potential Energy          ', &
                         'Latent Heat Energy        ', &
                         'Total Dry Energy          ', &
                         'Total Wet Energy          '/)


CONTAINS

!.................................................................
      subroutine tv2t_ (t,q)
      use m_const, only : zvir
      implicit none
      real, intent(inout) :: t(:,:,:)
      real, intent(in)    :: q(:,:,:)
      integer :: i,j,k
      integer :: im,jm,km
      im=size(t,1)
      jm=size(t,2)
      km=size(t,3)
      do k=1,km
        do j=1,jm
          do i=1,im
            t(i,j,k)=t(i,j,k)/(1.d0+zvir*q(i,j,k))
          enddo
        enddo
      enddo
      end subroutine tv2t_ 
!.................................................................
      subroutine scale_by_totene_(x,eps_eer,anorm,jnorm,projlon,projlat,projlev,&
                                  nymd,nhms,ntype,optene,normlz,vnorm,ps,delp)
      use m_const, only: cp=>cpm
      use m_const, only: rd=>rgas
      use m_const, only: tref=>tstd
      use m_const, only: pstd
      use m_const, only: alhl

!
! This code has been checked and shown to agree (very closely) with the
! energy result from initadj. Note that initadj embeds the following 
! definition of J and dJ/dx:
!           J = e^T E e
!       dJ/dx = 2 E e
! with E = 1/2 sum(energy_density) dA dp
!
      implicit none

      type(dyn_vect),intent(inout) :: x
      real, intent(in) :: eps_eer
      character(len=*), intent(in) :: anorm
      character(len=*), intent(in) :: jnorm
      real,    intent(in) :: projlat(2), projlon(2)
      integer, intent(in) :: projlev(2)
      integer, intent(in) :: nymd, nhms
      character(len=*), intent(in), optional :: ntype ! norm type: L2 or Ene (default)
      integer, intent(in), optional :: optene ! -2=E^-1; -1=E^-1/2; 1=E^1/2; 2=E
      integer, intent(in), optional :: vnorm  ! 0=mass-weight; 1=height-weight
      logical, intent(in), optional :: normlz ! vector will come out normalized to 1
      real, intent(in), optional :: ps(:,:)
      real, intent(in), optional :: delp(:,:,:)

      real, allocatable :: w2d(:,:)
      real, allocatable :: w3d(:,:,:)

      real, parameter :: pref=100.*pstd
      real fac,fact,ufac,tfac,qfac,pfac
      real dot(6)
      integer :: i,j,k,npt
      integer :: optene_
      logical :: normlz_
 
      normlz_ = .false.
      if(present(normlz)) then
        normlz_ = normlz
      endif
      optene_ = 1
      if(present(optene)) then
        optene_ = optene
      endif
      allocate(w2d(x%grid%im,x%grid%jm))
      allocate(w3d(x%grid%im,x%grid%jm,x%grid%km))

      call TotEne_Weights_(x,projlon,projlat,projlev,w2d,w3d,optene=optene_,&
                           vnorm=vnorm,ps=ps,delp=delp)

      print *, 'norm used: ', trim(anorm)
      print *, 'ene-scale: eps_eer: ', eps_eer

      if ( trim(ntype) == 'L2' .or. trim(ntype) == 'l2' ) then
         pfac = 1.0
         ufac = 1.0
         tfac = 1.0
         qfac = 1.0
      else
         fact=1.0/sqrt(2.) ! 1/2 factor in front of energy definition
         pfac = fact*sqrt(rd*tref)/pref
         ufac = fact
         tfac = fact*sqrt(cp/tref)
         qfac = fact*alhl*sqrt(eps_eer/(cp*tref))
      endif

      if(optene_==-1) then
        pfac = 1.0/pfac
        ufac = 1.0/ufac
        tfac = 1.0/tfac
        if(eps_eer>1.e-10) then
           qfac = 1.0/qfac
        else
           qfac = 0.0
        endif
      endif
      if(optene_==-2) then
        pfac = 1.0/(pfac*pfac)
        ufac = 1.0/(ufac*ufac)
        tfac = 1.0/(tfac*tfac)
        if(eps_eer>1.e-10) then
           qfac = 1.0/(qfac*qfac)
        else
           qfac = 0.0
        endif
      endif
      if(optene_==2) then
        pfac = pfac*pfac
        ufac = ufac*ufac
        tfac = tfac*tfac
        qfac = qfac*qfac
      endif

!     LPO: to be done properly ...
      do k=1,x%grid%km
         x%u (:,:,k)   = x%u (:,:,k)   * ufac * w3d(:,:,k)
         x%v (:,:,k)   = x%v (:,:,k)   * ufac * w3d(:,:,k)
         x%pt(:,:,k)   = x%pt(:,:,k)   * tfac * w3d(:,:,k)
      enddo
      if (trim(anorm)=='twe') then
         do k=1,x%grid%km
            x%q (:,:,k,1) = x%q (:,:,k,1) * qfac * w3d(:,:,k)
         enddo
      else
         x%q (:,:,:,1) = 0.0
      endif
      x%q(:,:,:,2:) = 0.0
      x%ts = 0.0
      x%delp = 0.0
      do k=1,x%grid%km
         x%delp(:,:,k) = x%ps * pfac * w3d(:,:,k)
      enddo
      x%ps = x%ps * pfac * w2d
      if (optene_>0.or.normlz_) then
         call dotp_(x,dot,anorm,jnorm,nymd,nhms)
      endif
      if ( normlz_ ) then
         fac = 1.0/sqrt(dot(1))
         x%u = fac * x%u
         x%v = fac * x%v
         x%pt= fac * x%pt
         x%q = fac * x%q
         x%ps= fac * x%ps
         x%delp = fac * x%delp
         ! univariate normalization ...
!        fac = 1.0
!        x%u = fac * x%u /sqrt(dot(2))
!        x%v = fac * x%v /sqrt(dot(3))
!        x%pt= fac * x%pt/sqrt(dot(4))
!        x%ps= fac * x%ps/sqrt(dot(5))
!        x%q = fac * x%q /sqrt(dot(6))
!        x%delp = fac * x%delp /sqrt(dot(5))
         call dotp_(x,dot,anorm,jnorm,nymd,nhms)
      endif

      deallocate(w3d)
      deallocate(w2d)
      return
      end subroutine scale_by_totene_

!.................................................................

      subroutine TotEne_Weights_(x,projlon,projlat,projlev,w2d,w3d,optene,vnorm,ps,delp)
      use m_const, only: pstd
      implicit none
      type(dyn_vect) x
      real,    intent(in) :: projlat(2), projlon(2)
      integer, intent(in) :: projlev(2)
      integer, intent(in), optional :: optene
      integer, intent(in), optional :: vnorm
      real,intent(in),optional :: ps(:,:)
      real,intent(in),optional :: delp(:,:,:)
      real :: w2d(:,:), w3d(:,:,:)
      
      real, parameter :: pref=100.*pstd
      real, allocatable :: ple (:,:,:)
      real, allocatable :: dsig(:,:,:)
      real, allocatable :: hlpo(:,:)
      real, allocatable :: vlpo(:)
      real, allocatable :: rlat(:),rlon(:)
      real, allocatable :: jweights(:,:),glats(:,:)
      real pi,crlat,sumcl
      real clat,clon
      real rlon1,rlon2,rlat1,rlat2
      real rlon0
      integer im,jm,km
      integer i,j,k
      integer optene_,vnorm_

      optene_=1
      if(present(optene)) then
         optene_=optene
      endif
      vnorm_=0
      if(present(vnorm)) then
         vnorm_=vnorm
         if(.not.present(delp)) then
           print *, 'need delp for height-weights'
           call exit (99)
         endif
      endif
      im=x%grid%im
      jm=x%grid%jm
      km=x%grid%km
      pi=4.0*atan(1.0)
      if (x%grid%lon(1)<0.0) then
         rlon1 = pi * (projlon(1)-180.)/180.
         rlon2 = pi * (projlon(2)-180.)/180.
      else
         rlon1 = pi * projlon(1)/180.
         rlon2 = pi * projlon(2)/180.
      endif
      rlat1 = pi * projlat(1)/180.
      rlat2 = pi * projlat(2)/180.

      allocate(hlpo(im,jm))
      allocate(vlpo(km))
      allocate(rlon(im),rlat(jm))

      vlpo=0.0
      do k=1,km
         if(k >= projlev(1) .and. k <= min(km,projlev(2)) ) vlpo(k) = 1.0
      enddo

      rlon0 = pi*x%grid%lon(1)/180.
      rlat  = pi*x%grid%lat/180.
      rlon  = pi*x%grid%lon/180.
      hlpo=0.0
      if (rlon1 < rlon2) then
         do j = 1, jm
            clat = rlat(j)
            if ( clat >= rlat1 .and. clat <= rlat2 ) then
              do i = 1, im
                 if ( rlon(i) >= rlon1 .and. rlon(i) <= rlon2 ) hlpo(i,j) = 1.0
              end do
            end if
         end do
      else
         do j = 1, jm
            clat = rlat(j)
            if ( clat >= rlat1 .and. clat <= rlat2 ) then
              do i = 1, im
                 if ( rlon(i) >= rlon1 .or. rlon(i) <= rlon2 ) hlpo(i,j) = 1.0
              end do
            end if
         end do
      endif ! (rlon1 < rlon2 for t)

      allocate (jweights(jm,2),glats(jm,2))
      call horiz_grid_ (jm, jweights, glats)
!     2d weights
      do j=1,jm
         crlat=cos(pi*x%grid%lat(j)/180.)
         do i=1,im
           w2d(i,j) = jweights(j,2)
         enddo
      enddo
      deallocate (jweights,glats)
      sumcl = 1.0 / im ! this is so it match Ron''s weights more closely
      if (optene_==2) then
         do j=1,jm
            w2d(:,j)=hlpo(:,j)*(w2d(:,j)*sumcl)
         enddo
      endif
      if (optene_==-2) then
         do j=1,jm
            w2d(:,j)=hlpo(:,j)/(w2d(:,j)*sumcl)
         enddo
      endif
      if (optene_==-1) then
         do j=1,jm
            w2d(:,j)=hlpo(:,j)/sqrt(w2d(:,j)*sumcl)
         enddo
      endif
      if (optene_==1) then
         w2d = sqrt(w2d*sumcl)
         w2d = w2d*hlpo
      endif

!     3d weights
      allocate(ple (im,jm,km+1))
      allocate(dsig(im,jm,km  ))
      if (present(ps)) then
         if ( present(delp) ) then
            if (vnorm_>0) then
               print *, 'using ref delp in height-weights'
               ple(:,:,1) = x%grid%ak(1)
               do k=2,km+1
                  ple(:,:,k) = ple(:,:,k-1) + delp(:,:,k-1)
               enddo
               do k=1,km
                  dsig(:,:,k) = log(ple(:,:,k+1)/ple(:,:,k))/log(ple(:,:,km+1)/ple(:,:,1))
               enddo
            else
               print *, 'using ref delp and ps in mass-weights'
               do k=1,km
                  dsig(:,:,k)=delp(:,:,k)/ps(:,:)
               enddo
            endif 
         else
            do k=1,km+1
               ple(:,:,k)=x%grid%ak(k)+x%grid%bk(k)*ps     !state dependent as in initadj
            enddo
            do k=1,km
               dsig(:,:,k)=(ple(:,:,k+1)-ple(:,:,k))/ps
            enddo
         endif
      else
         do k=1,km+1
            ple(:,:,k)=x%grid%ak(k)+x%grid%bk(k)*pref      ! set to ref pressure
         enddo
         do k=1,km
            dsig(:,:,k)=(ple(:,:,k+1)-ple(:,:,k))/pref
         enddo
      endif
      if (optene_==1) then
         dsig=sqrt(dsig)
      endif
      if (optene_==-1) then
         dsig=1.0/sqrt(dsig)
      endif
      if (optene_==-2) then
         dsig=1.0/dsig
      endif
      if (optene_==2) then
         dsig=dsig
      endif
      do k=1,km
         w3d(:,:,k)=vlpo(k)*dsig(:,:,k)*w2d
      enddo

      deallocate(rlon,rlat)
      deallocate(vlpo)
      deallocate(hlpo)
      deallocate(dsig)
      deallocate(ple) 
      end subroutine TotEne_Weights_

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: horiz_grid_ --- Determine some D-grid information required
!
!
! !INTERFACE:
!
      subroutine horiz_grid_ (jn1, jweights, glats)

!USES:

      implicit none

! !INPUT PARAMETERS:

      integer,  intent(in)  :: jn1    

! !OUTPUT PARAMETERS:

      real, intent(out) :: jweights(jn1,2) ! area weights (1) u, (2) T   
      real, intent(out) :: glats(jn1,2)    ! degrees lats (1) u, (2) T   

! !INPUT/OUTPUT PARAMETERS:
 
! !DESCRIPTION:
!
!  Determine some D-grid information required
!
!  i=1 on v field corresponds to long=0
!  i=1 on T,u fields corresponds to long=0+360/(2*im1)
!  v defined on same lats as T, excluding poles
!  u defined between T lats, but on T lons
!

! !SEE ALSO:
!
  
! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!
!EOP
!-------------------------------------------------------------------------
!
!  local variables

      integer   :: i,j
      real  :: slats(jn1,2)  ! sines of latitude for (1)u and (2)T
      real  :: pi, pi180
      real  :: rlat, rlat_half
      real  :: tlat, ulat
      
      pi=4.d0*datan(1.d0)
      pi180=180.d0/pi
      rlat=pi/dble(jn1-1)
      rlat_half=0.5d0*rlat

      tlat=-0.5d0*pi   ! latitude of south pole
      glats(1,1)=pi180*rlat  ! a special value since otherwise not used
      glats(1,2)=pi180*tlat  
      slats(1,1)=0.d0  ! a value not used
      slats(1,2)=-1.d0
      do j=2,jn1
        ulat=tlat+rlat_half
        tlat=tlat+rlat
        glats(j,1)=pi180*ulat
        glats(j,2)=pi180*tlat
        slats(j,1)=dsin(real(ulat,8))
        slats(j,2)=dsin(real(tlat,8))
      enddo
!
       jweights(1,1)=0.d0  ! not used
       jweights(1,2)=0.5d0*(1.d0+slats(2,1))
       do j=2,jn1-1
         jweights(j,1)=0.5d0*(slats(j,2)-slats(j-1,2))
         jweights(j,2)=0.5d0*(slats(j+1,1)-slats(j,1))
       enddo
       jweights(jn1,1)=0.5d0*(slats(jn1,2)-slats(jn1-1,2))
       jweights(jn1,2)=0.5d0*(1.d0-slats(jn1,1))
    
       end subroutine horiz_grid_

!.................................................................

      subroutine dotp_(x,dot,anorm,jnorm,nymd,nhms)
      use m_ioutil, only : luavail
      implicit none
      character(len=*), intent(in) :: anorm, jnorm
      integer nymd, nhms
      real dot(:)
      type(dyn_vect) x
      integer lu
      dot(2) = sum(x%u *x%u )
      dot(3) = sum(x%v *x%v )
      dot(4) = sum(x%pt*x%pt)
      dot(5) = sum(x%ps*x%ps)
      dot(6) = sum(x%q(:,:,:,1:1)*x%q(:,:,:,1:1))
      dot(1) = sum(dot(2:6))
      write(*,'(a)') 'sum, u, v, t, ps, q'
      write(*,'(6(1x,f15.12))') dot
      if (trim(jnorm) /= "NULL" ) then
         lu = luavail()
         open  (lu, file=trim(jnorm), form='formatted')
         write (lu, '(a,i8.8,a,i6.6,2a)') 'Date: ', nymd, ' Time: ', nhms, ' Norm type: ', trim(anorm)
         write (lu, '(a)') 'sum, u, v, t, ps, q'
         write (lu, '(6f20.12)') dot
         close (lu)
      endif
      call get_ene_field_(x,nymd,nhms)
      end subroutine dotp_

   subroutine get_ene_field_(x,nymd,nhms)

      use m_GFIO_PutFld,only : GFIO_PutFld
      implicit none

      integer, intent(in) :: nymd,nhms
      type(dyn_vect) x

      real, allocatable :: energy_field(:,:,:,:)

      integer k,im,jm,km,ier

      im = x%grid%im
      jm = x%grid%jm
      km = x%grid%km
      allocate(energy_field(im,jm,km,nall))

      energy_field = 0.0

      energy_field(:,:,:,nkxe) = x%u*x%u + x%v*x%v
      energy_field(:,:,:,nape) = x%pt*x%pt
      energy_field(:,:,:,nqxe) = x%q(:,:,:,1)*x%q(:,:,:,1)
      energy_field(:,:,:,ngpe) = x%delp*x%delp

      energy_field(:,:,:,nape) = energy_field(:,:,:,nape) + energy_field(:,:,:,ngpe)
      energy_field(:,:,:,ntxe) = energy_field(:,:,:,nape) + energy_field(:,:,:,nkxe)
      energy_field(:,:,:,ntwe) = energy_field(:,:,:,nqxe) + energy_field(:,:,:,ntxe)

      call GFIO_PutFld ( 'gmao_pert_energy',  'GMAO',  'enepert.nc4',  &
                          nymd, nhms, 240000, &
                          im,jm,km,x%grid%ptop,x%grid%ks,x%grid%ak,x%grid%bk, &
                          nall, vnames, vtitles, vunits, &
                          energy_field, ier, untag=.true. )

      print *, 'Energy partition: '
      print *, 'Kinetic          Energy: ', sum(energy_field(:,:,:,nkxe))
      print *, 'Moist Static     Energy: ', sum(energy_field(:,:,:,nqxe))
      print *, 'Potential        Energy: ', sum(energy_field(:,:,:,ngpe))
      print *, 'Avail. Potential Energy: ', sum(energy_field(:,:,:,nape))
      print *, 'Total dry        Energy: ', sum(energy_field(:,:,:,ntxe))
      print *, 'Total wet        Energy: ', sum(energy_field(:,:,:,ntwe))

      deallocate(energy_field)

  end subroutine get_ene_field_

!.................................................................

  subroutine becov_qsat_(qsat,temp,pmk,lat2,lon2,icesat)
!
!   input argument list:
!     temp     - virtual temperature (K)
!     pmk      - mean-layer pressure (Pa)
!     lat2     - number of in the sub-domain array
!     lon2     - number of longitudes in the sub-domain array
!     icesat   - logical flag:  T=include ice and ice-water effects,
!                depending on t, in qsat and esat calcuations.
!                otherwise, compute esat and qsat with respect to water surface
!                depending on t, in qsat and esat calcuations.
!                otherwise, compute esat and qsat with respect to water surface
!
!   output argument list:
!     qsat     - specific humidity (input), saturation specific humidity (output)
!  10Jun2020 Todling  pulled it from NCEP_bkgecov into here
!                     (a) 3d into 2d
!                     (b) note that routine is invariant to lat/lon and lon/lat 
!                     (c) note that routine is invariant wrt vertical grid
!                     orientation
!
  use m_realkinds, only : fp_kind => kind_r4 
! use m_const, only: fv=>zvir 
! use m_const, only: rv=>rvap
! use m_const, only: rd=>rgas
! use m_const, only: cvap=>cpv
! use m_const, only: hvap=>alhl
  implicit none

  logical,intent(in):: icesat
  integer,intent(in):: lat2,lon2
  real(fp_kind),intent(inout),dimension(:,:):: qsat
  real(fp_kind),intent(in),dimension(:,:):: temp
  real(fp_kind),intent(in),dimension(:,:):: pmk

  real,parameter :: rd     = 2.8705e+2_fp_kind              !  gas constant of dry air (J/kg/K)
  real,parameter :: rv     = 4.6150e+2_fp_kind              !  gas constant of h2o vapor (J/kg/K)
  real,parameter :: fv     = rv/rd-1._fp_kind               !  used in virtual temp.  equation (adim)
  real,parameter :: cvap   = 1.8460e+3_fp_kind              !  specific heat of h2o vapor      (J/kg/K)
  real,parameter :: hvap   = 2.5000e+6_fp_kind              !  latent heat of h2o condensation (J/kg)

  real,parameter :: cliq   = 4.1855e+3_fp_kind              !  specific heat of liquid h2o     (J/kg/K)
  real,parameter :: psat   = 6.1078e+2_fp_kind              !  pressure at h2o triple point    (Pa)
  real,parameter :: csol   = 2.1060e+3_fp_kind              !  specific heat of solid h2o (ice)(J/kg/K)
  real,parameter :: ttp    = 2.7316e+2_fp_kind              !  temperature at h2o triple point (K)
  real,parameter :: hfus   = 3.3358e+5_fp_kind              !  latent heat of h2o fusion       (J/kg)

  real,parameter :: zero=0.0_fp_kind
  real,parameter :: one=1.0_fp_kind

  real,parameter :: dldt=cvap-cliq
  real,parameter :: xa=-dldt/rv
  real,parameter :: xb=xa+hvap/(rv*ttp)
  real,parameter :: hsub=hvap+hfus
  real,parameter :: dldti=cvap-csol 
  real,parameter :: xai=-dldti/rv
  real,parameter :: xbi=xai+hsub/(rv*ttp)
  real,parameter :: tmix=ttp-20.0_fp_kind
  real,parameter :: eps=rd/rv
  real,parameter :: omeps = 1._fp_kind-eps
  integer k,j,i
  real(fp_kind) pw,q,tdry,tr,es,qs,esi,esw
  real(fp_kind) w,pscl,esmax
!
  pscl = 1.00_fp_kind ! input pressure in Pa as required

  if (icesat) then
!    do k = 1,nsig
        do j = 1,lon2
           do i = 1,lat2

              pw  = pmk(i,j)
              pw  = pscl*pw
! maximum vapor pressure 5% of atmospheric pressure
              esmax=0.05*pw

              q  = qsat(i,j)
              if (q.lt.zero) q=zero

              tdry = temp(i,j)/(one+fv*q)
              tr = ttp/tdry
              if (tdry >= ttp) then
                 es = psat * (tr**xa) * exp(xb*(one-tr))
              elseif (tdry < tmix) then
                 es = psat * (tr**xai) * exp(xbi*(one-tr))
              else
                 w  = (tdry - tmix) / (ttp - tmix)
                 es =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                      + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))
              endif

              es = min(es,esmax)
              qs = eps * es / (pw - omeps * es)

              if (qs.lt.qsat(i,j)) then
                 tdry = temp(i,j)/(one+fv*qs)
                 tr = ttp/tdry
                 if (tdry >= ttp) then
                    es = psat * (tr**xa) * exp(xb*(one-tr))
                 elseif (tdry < tmix) then
                    es = psat * (tr**xai) * exp(xbi*(one-tr))
                 else
                    w  = (tdry - tmix) / (ttp - tmix)
                    es =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                         + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))
                 endif

                 es = min(es,esmax)
                 qs = eps * es / (pw - omeps * es)

                 tdry = temp(i,j)/(one+fv*qs)
                 tr = ttp/tdry
                 if (tdry >= ttp) then
                    es = psat * (tr**xa) * exp(xb*(one-tr))
                 elseif (tdry < tmix) then
                    es = psat * (tr**xai) * exp(xbi*(one-tr))
                 else
                    w  = (tdry - tmix) / (ttp - tmix)
                    es =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                         + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))
                 endif
                 es = min(es,esmax)
                 qs = eps * es / (pw - omeps * es)

              end if

              qsat(i,j) = qs

           end do
        end do
!    end do

!
!     Compute saturation values with respect to water surface
  else
!    do k = 1,nsig
        do j = 1,lon2
           do i = 1,lat2

              pw = pmk(i,j)
              pw  = pscl*pw
! maximum vapor pressure 5% of atmospheric pressure
              esmax=0.05*pw

              q  = qsat(i,j)
              if (q.lt.zero) q=zero

              tdry = temp(i,j)/(one+fv*q)
              tr = ttp/tdry
              es = psat * (tr**xa) * exp(xb*(one-tr))
              es = min(es,esmax)
              qs = eps * es / (pw - omeps * es)

              if (qs.lt.qsat(i,j)) then
                 tdry = temp(i,j)/(one+fv*qs)
                 tr = ttp/tdry
                 es = psat * (tr**xa) * exp(xb*(one-tr))
                 es = min(es,esmax)
                 qs = eps * es / (pw - omeps * es)
                 tdry = temp(i,j)/(one+fv*qs)
                 tr = ttp/tdry
                 es = psat * (tr**xa) * exp(xb*(one-tr))
                 es = min(es,esmax)
                 qs = eps * es / (pw - omeps * es)

              end if
              qsat(i,j) = qs


           end do
        end do
!    end do

  endif

  return
  end subroutine becov_qsat_

!.................................................................

  subroutine gsi_qsat_(qsat,tsen,prsl,lat2,lon2,nsig,ice,ntop)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    genqsat
!   prgmmr: derber           org: np23                date: 1998-01-14
!
! abstract: obtain saturation specific humidity for given temperature.
!
! program history log:
!   1998-01-14  derber
!   1998-04-05  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   1903-10-07  Wei Gu, bug fixes,if qs<0,then set qs=0; merge w/ GSI by R Todling
!   2003-12-23  kleist, use guess pressure, adapt module framework
!   2004-05-13  kleist, documentation
!   2004-06-03  treadon, replace ggrid_g3 array with ges_* arrays
!   2005-02-23  wu, output dlnesdtv
!   2005-11-21  kleist, derber  add dmax array to decouple moisture from temp and
!               pressure for questionable qsat
!   2006-02-02  treadon - rename prsl as ges_prsl
!   2006-09-18  derber - modify to limit saturated values near top
!   2006-11-22  derber - correct bug:  es<esmax should be es<=esmax
!   2008-06-04  safford - rm unused vars
!   2010-03-23  derber - simplify and optimize
!   2010-03-24  derber - generalize so that can be used for any lat,lon,nsig and any tsen and prsl (for hybrid)
!   2010-12-17  pagowski - add cmaq
!   2011-08-15  gu/todling - add pseudo-q2 options
!   2014-12-03  derber - add additional threading
!   2018-02-15  wu - add code for fv3_regional option
!   2020-05-14  todling - opt arg ntop to accommodate flipped pressure levels wrt to GSI
!                       - update top mid-pressure to 0 Pa
!
!   input argument list:
!     tsen      - input sensibile temperature field (lat2,lon2,nsig)
!     prsl      - input layer mean pressure field (lat2,lon2,nsig)
!     lat2      - number of latitudes                              
!     lon2      - number of longitudes                             
!     nsig      - number of levels                              
!     ice       - logical flag:  T=include ice and ice-water effects,
!                 depending on t, in qsat calcuations.
!                 otherwise, compute qsat with respect to water surface
!     ntop      - index of top level (optional)
!
!   output argument list:
!     qsat      - saturation specific humidity (output)
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use m_realkinds, only : r_kind => kind_r4
  use m_const, only: fv=>zvir
  use m_const, only: rv=>rvap
  use m_const, only: rd=>rgas
  use m_const, only: cvap=>cpv
  use m_const, only: hvap=>alhl
  use m_const, only: hfus=>alhf
  use m_const, only: cliq=>capwtr
  use m_const, only: csol=>capice
! use constants, only: xai,tmix,xb,omeps,eps,xbi,one,xa,psat,ttp
  implicit none

  logical                      ,intent(in   ) :: ice
  integer                      ,intent(in   ) :: lat2,lon2,nsig
  real(r_kind),dimension(lat2,lon2,nsig),intent(  out) :: qsat
  real(r_kind),dimension(lat2,lon2,nsig),intent(in   ) :: tsen,prsl
  integer,optional             ,intent(in   ) :: ntop

! real,parameter :: rd     = 2.8705e+2_r_kind              !  gas constant of dry air (J/kg/K)
! real,parameter :: rv     = 4.6150e+2_r_kind              !  gas constant of h2o vapor (J/kg/K)
! real,parameter :: fv     = rv/rd-1._r_kind               !  used in virtual temp.  equation (adim)
! real,parameter :: cvap   = 1.8460e+3_r_kind              !  specific heat of h2o vapor      (J/kg/K)
! real,parameter :: hvap   = 2.5000e+6_r_kind              !  latent heat of h2o condensation (J/kg)

! real,parameter :: cliq   = 4.1855e+3_r_kind              !  specific heat of liquid h2o     (J/kg/K)
  real,parameter :: psat   = 6.1078e+2_r_kind              !  pressure at h2o triple point    (Pa)
! real,parameter :: csol   = 2.1060e+3_r_kind              !  specific heat of solid h2o (ice)(J/kg/K)
  real,parameter :: ttp    = 2.7316e+2_r_kind              !  temperature at h2o triple point (K)
! real,parameter :: hfus   = 3.3358e+5_r_kind              !  latent heat of h2o fusion       (J/kg)

  real,parameter :: zero=0.0_r_kind
  real,parameter :: one=1.0_r_kind

  real,parameter :: dldt=cvap-cliq
  real,parameter :: xa=-dldt/rv
  real,parameter :: xb=xa+hvap/(rv*ttp)
  real,parameter :: hsub=hvap+hfus
  real,parameter :: dldti=cvap-csol
  real,parameter :: xai=-dldti/rv
  real,parameter :: xbi=xai+hsub/(rv*ttp)
  real,parameter :: tmix=ttp-20.0_r_kind
  real,parameter :: eps=rd/rv
  real,parameter :: omeps = 1._r_kind-eps

! Declare local parameters
  integer         k,j,i,nbot,ntop_
  real(r_kind) pw,tdry,tr,es,es2
  real(r_kind) w,onep3,esmax
  real(r_kind) esi,esw
  real(r_kind),dimension(lat2):: mint,estmax
  integer,dimension(lat2):: lmint

  onep3 = 1.0_r_kind
  ntop_ = nsig ! GSI top level index
  if(present(ntop)) then
    ntop_ = ntop
  endif
  nbot = nsig-ntop_+1

  do j=1,lon2
     do i=1,lat2
        mint(i)=340._r_kind
        lmint(i)=ntop_
     end do
     do k=1,nsig
        do i=1,lat2
           if((prsl(i,j,k) < 30._r_kind  .and.  &
               prsl(i,j,k) >  0._r_kind) .and.  &
               tsen(i,j,k) < mint(i))then
              lmint(i)=k
              mint(i)=tsen(i,j,k)
           end if
        end do
     end do
     do i=1,lat2
        tdry = mint(i)
        if( abs(tdry) < 1.0e-8_r_kind ) tdry = 1.0e-8_r_kind
        tr = ttp/tdry
        if (tdry >= ttp .or. .not. ice) then
           estmax(i) = psat * (tr**xa) * exp(xb*(one-tr))
        elseif (tdry < tmix) then
           estmax(i) = psat * (tr**xai) * exp(xbi*(one-tr))
        else
           w  = (tdry - tmix) / (ttp - tmix)
           estmax(i) =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                   + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))
        endif
     end do

     do k = 1,nsig
        do i = 1,lat2
           tdry = tsen(i,j,k)
           if( abs(tdry) < 1.0e-8_r_kind ) tdry = 1.0e-8_r_kind
           tr = ttp/tdry
           if (tdry >= ttp .or. .not. ice) then
              es = psat * (tr**xa) * exp(xb*(one-tr))
           elseif (tdry < tmix) then
              es = psat * (tr**xai) * exp(xbi*(one-tr))
           else
              esw = psat * (tr**xa) * exp(xb*(one-tr)) 
              esi = psat * (tr**xai) * exp(xbi*(one-tr)) 
              w  = (tdry - tmix) / (ttp - tmix)
!             es =  w * esw + (one-w) * esi
              es =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                       + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))

           endif

           pw = onep3*prsl(i,j,k)
           esmax = es
           if (nbot==1) then
              if(lmint(i) < k)then
                 esmax=0.1_r_kind*pw
                 esmax=min(esmax,estmax(i))
              end if
           else
              if(lmint(i) > k)then
                 esmax=0.1_r_kind*pw
                 esmax=min(esmax,estmax(i))
              end if
           end if
           es2=min(es,esmax)
           qsat(i,j,k) = eps * es2 / (pw - omeps * es2)

        end do
     end do
  end do
  return
  end subroutine gsi_qsat_

!.................................................................

  end module m_dyn_util
