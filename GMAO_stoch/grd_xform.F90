   module grd_xform

! !USES:
   use mod_param
   use m_const, only : undef

   implicit none

   contains

!-------------------------------------------------------------------------
!
   subroutine  gs2gd ( gfld,glon,glat,nsig,rfld,rlon,rlat)
!
! !INPUT PARAMETERS:
  integer,              intent(in)  :: glon,glat,rlon,rlat,nsig
  real(kind=kind_evod), intent(in)  :: gfld(glon,glat,nsig)
  real(kind=kind_evod), intent(out) :: rfld(rlon,rlat,nsig)
!-------------------------------------------------------------------------
  real, pointer :: dlam(:) => null()
  real, pointer :: dphi(:) => null()
  real, pointer :: lons(:) => null()
  real, pointer :: lats(:) => null()
  real, pointer :: dumr(:,:,:) => null()
  real, pointer :: dumg(:,:,:) => null()
  integer, parameter :: ONE=1  
  integer, parameter :: SFL=1      ! flag for scalar 
  integer, parameter :: VFL=-1     ! flag for vector 
  integer, parameter :: ORD=3      ! order of interpolation
  real, parameter    :: ROT = 0.0  ! Rotation parameter for lam
  real, parameter    :: TLT = 90.0 ! Rotation parameter for phi 
  real, parameter    :: PRE = 0.0  ! Rotation parameter lam_0
  logical, parameter :: FLG = .true.

   allocate(dumg(glon,glat,nsig)) 
   allocate(dumr(rlon,rlat,nsig)) 
   allocate(dlam(glon))
   allocate(dphi(glat))
   allocate(lons(rlon*rlat))
   allocate(lats(rlon*rlat))

   call comput_gaus(glon,glat,rlon,rlat,dphi,dlam,lons,lats,.false.)

   dumg = gfld
   call interp_h ( dumg,glon,glat,nsig,dlam,dphi,ROT,TLT,PRE,dumr, &
                   rlon*rlat,lons,lats,VFL,ORD,FLG,undef)
   rfld = dumr

   deallocate(lats)
   deallocate(lons)
   deallocate(dphi)
   deallocate(dlam)
   deallocate(dumr)
   deallocate(dumg)
    
   end subroutine gs2gd

!
   subroutine  gd2gs ( rfld,rlon,rlat,nsig,gfld,glon,glat)
!
! !INPUT PARAMETERS:
!
   integer, intent(in)               :: glon, glat 
   integer, intent(in)               :: rlon, rlat 
   real(kind=kind_evod), intent(in)  :: rfld(rlon,rlat,nsig)
   real(kind=kind_evod), intent(out) :: gfld(glon,glat,nsig)
   integer, intent(in)               :: nsig      
!
! !OUTPUT PARAMETERS:
!
!
!-------------------------------------------------------------------------
   integer ::  gxg
   real, pointer :: dlam(:) => null()
   real, pointer :: dphi(:) => null()
   real, pointer :: lons(:,:) => null()
   real, pointer :: lats(:,:) => null()
   real, pointer :: dumr(:,:,:) => null()
   real, pointer :: dumg(:,:,:) => null()
   integer, parameter :: ONE=1
   integer, parameter :: SFL=1      ! flag for scalar
   integer, parameter :: VFL=-1     ! flag for vector
   integer, parameter :: ORD=3      ! order of interpolation
   real, parameter    :: ROT = 0.0 ! Rotation parameter for lam 
!   real, parameter    :: ROT = -180.0 ! Rotation parameter for lam 
   real, parameter    :: TLT = 90.0 ! Rotation parameter for phi
   real, parameter    :: PRE = 0.0  ! Rotation parameter lam_0
   logical, parameter :: FLG = .true.

   gxg = glon*glat
 
   allocate(dumr(rlon,rlat,nsig)) 
   allocate(dumg(glon,glat,nsig)) 
   allocate(dlam(rlon))
   allocate(dphi(rlat))
   allocate(lons(glon,glat))
   allocate(lats(glon,glat))

   call comput_gaus(rlon,rlat,glon,glat,dphi,dlam,lons,lats,.false.)

! regrid from latlon grid to a gaussian grid

   dumr= rfld
   call interp_h ( dumr,rlon,rlat,nsig,dlam,dphi,ROT,TLT,PRE,dumg, &
                   gxg,lons,lats,SFL,ORD,FLG,undef)
   gfld = dumg

   deallocate(lats)
   deallocate(lons)
   deallocate(dphi)
   deallocate(dlam)
   deallocate(dumg)
   deallocate(dumr)
   
   end subroutine gd2gs

   subroutine  poleuv ( ufld,vfld,rlon,rlat,nsig )

   integer, intent(in)                 :: rlon, rlat,nsig 
   real(kind=kind_evod), intent(inout) :: ufld(rlon,rlat,nsig),vfld(rlon,rlat,nsig)
   real(kind=kind_evod), allocatable   :: sinl(:),cosl(:)
   real(kind=kind_evod)                :: DL,LON,PI,UP,VP
   integer                             :: i,m,J,N,L

   PI = atan(1.0)*4.0         !constant
   DL = 2*PI/rlon
   allocate(sinl(rlon))
   allocate(cosl(rlon))
   do i=1,rlon
          LON = -PI + (i-1)*DL 
      cosl(i) = cos(LON)
      sinl(i) = sin(LON)
   enddo
   do L=1,nsig
      do m = 1,2
         N = (-1)**m
         if(m.eq.1) J = 1
         if(m.eq.2) J = rlat
         UP = 0.0
         VP = 0.0
         do i = 1,rlon
            UP = UP -   ufld(i,J-N,L)*sinl(i) - N*vfld(i,J-N,L)*cosl(i)
            VP = VP + N*ufld(i,J-N,L)*cosl(i) -   vfld(i,J-N,L)*sinl(i) 
         enddo 
         UP = UP / rlon
         VP = VP / rlon
         do i = 1,rlon
            ufld(i,J,L) = -  UP*sinl(i) + N*VP*cosl(i)
            vfld(i,J,L) = -N*UP*cosl(i) -   VP*sinl(i)
         enddo 
      enddo 
   enddo 

   deallocate(sinl,cosl)

   end subroutine poleuv
 end module grd_xform
   
