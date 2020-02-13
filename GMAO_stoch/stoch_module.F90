#include "MAPL_Generic.h"

module stoch_module

  use ESMF
  use MAPL
  use mod_param
  use stoch_data

  implicit none 

  public setup_pattern, clear_pattern
  public skeb_pattern
  public writout8,writout4

  type (LatlonGridFactory) :: ll_factory
  class(AbstractRegridder), pointer :: L2C => null()
  type(ESMF_Grid)          :: LLgrid
  integer                  :: isp,isk,ishutoff
  logical                  :: cubed,isallset

  CONTAINS
  !-----------------------------------------------------------------------------------
  subroutine setup_pattern(CF,GRID)

    use initialize

    implicit none
    !
    type(ESMF_Config),intent(inout)  :: CF
    type(ESMF_Grid)  ,intent(inout)  :: GRID

    real*8  :: DT 
    integer :: NX,NY,MEMID,ISPPT,ISKEB,ISHUT
    integer :: STATUS,ierr,rc
    integer :: im_world,jm_world,lm_world, lm
    integer :: im_gg,jm_gg
      
    call ESMF_ConfigGetAttribute(CF, NX      , Label="NX:"      , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, NY      , Label="NY:"      , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, im_world, Label="AGCM_IM:" , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, jm_world, Label="AGCM_JM:" , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, lm_world, Label="AGCM_LM:" , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, DT      , Label="RUN_DT:"  , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, MEMID   , Label="MEMID:"   , DEFAULT=999, RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, ISPPT   , Label="SPPT:"    , DEFAULT=0  , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, ISKEB   , Label="SKEB:"    , DEFAULT=0  , RC=STATUS)
    call ESMF_ConfigGetAttribute(CF, ISHUT   , Label="STOCH_shutoff:" , DEFAULT=1  , RC=STATUS)
  
    lm = lm_world 
    cubed=.false.
    if (jm_world==6*im_world) cubed=.true.

    call get_pert_grid(im_world,jm_world,im_ll,jm_ll,jcap,im_gg,jm_gg)

    if (MAPL_am_I_root()) then
      print*, 'Will be running with the options SPPT=', ISPPT ,'and SKEB=',ISKEB 
      call sub_init(im_gg,jm_gg,im_ll,jm_ll,lm,MEMID,DT)
      call sub_create()
      call sub_setup()                    
      if(ISPPT /= 0) call init_sppt
      if(ISKEB /= 0) call init_skeb
      call init_stochdata(DT,ls_node)
    endif 

    if (cubed) then
      ll_factory = LatLonGridFactory(grid_name='XYLLgridDC', &
                          Nx = Nx, Ny = Ny,  &
                          IM_World = IM_ll,  &
                          JM_World = JM_ll,  &
                          LM = LM, pole='PC', dateline='DC')
       LLgrid = grid_manager%make_grid(ll_factory)
       L2C => regridder_manager%make_regridder(LLGrid, Grid, REGRID_METHOD_BILINEAR)
    endif

    isallset = .true.
    isp = 0 ; isk = 0
    ishutoff = 999                 ! beyond 5days
    if (ISHUT /= 0) then 
       ishutoff = int(2*21600/DT)  ! shutoff after 12h model integration
       if (MAPL_am_I_root()) print*,' stochastic perts shutoff after ',ishutoff ,'steps'
    endif

  end subroutine setup_pattern
!-----------------------------------------------------------------------------------
  subroutine clear_pattern 

    use initialize
    implicit none

    if (MAPL_am_I_root()) then 
      call sub_clear
      call destroy_stochdata 
    endif 

  end subroutine clear_pattern
!-----------------------------------------------------------------------------------
  subroutine sppt_pattern(CF, GRID, RNDPTR, PREF, IM,JM,LM,DT)

    use initialize
    use m_die, only: die
    use stoch_data

    implicit none
    !
    character(len=ESMF_MAXSTR) :: COMP_NAME
    character(len=*),parameter :: fnout = 'rndpttn_latlon.grd'
    character(len=*),parameter :: fnoutc = 'rndpttn.grd'
    character(len=40)          :: fname
    character(len=8)           :: indx
    integer,parameter          :: lu=20
    !
    type(ESMF_Config),intent(inout)  :: CF
    type(ESMF_Grid)  ,intent(inout)  :: GRID
    !
    real*4, intent(inout) :: RNDPTR(IM,JM,LM)
    real*4, intent(in)    :: PREF(LM+1)
    real*4, intent(in)    :: DT 
    real*8, allocatable   :: RNDPTR_WORLD(:,:)
    real*4, allocatable   :: work2d(:,:)
    real*4, allocatable   :: rndptr_local(:,:,:)
    !
    integer, intent(in)   :: IM,JM,LM
    integer ::  DIMS(ESMF_MAXGRIDDIM)
    integer :: NX,NY,MEMID
    !
    integer :: STATUS,ierr,rc
    integer :: IM_world,JM_world,LM_world
    integer :: im_gg,jm_gg
    integer :: im_ll_local, jm_ll_local
    integer :: L 

    if (.not. isallset) call setup_pattern (CF,GRID)

    if (isp .ge. ishutoff -1) return 
    isp = isp + 1

    if (MAPL_am_I_root()) then
      allocate(RNDPTR_WORLD(IM_ll,JM_ll))
      call get_sppt(RNDPTR_WORLD,PREF)

     ! write (indx, "(I3.3)") isp
     ! fname='rndpttn_latlon'//trim(indx)//'.grd'
     ! print*, 'output fname ', trim(fname)
     ! open (lu,file=trim(fname),form='unformatted',access='sequential',convert='little_endian')
     !   call writout8(lu,RNDPTR_WORLD,im_ll,jm_ll,1)
     ! close(lu)
    endif 
      
    if (cubed) then
      call MAPL_GridGet(LLgrid, localCellCountPerDim=DIMS, RC=STATUS)
      im_ll_local=dims(1)
      jm_ll_local=dims(2)
      allocate(rndptr_local(IM_ll_local,JM_ll_local,lm))
      rndptr_local = 0.
    else
      im_ll_local = im
      jm_ll_local = jm
    endif

    allocate(work2d(im_ll,jm_ll))
    work2d =0.

    do L=1,lm !
      if (MAPL_am_I_root()) work2d=vfact_sppt(L)*RNDPTR_World(:,:)
      if(cubed) then
         call ArrayScatter(RNDPTR_local(:,:,L), work2d, LLgrid, rc=status)
      else
         call ArrayScatter(RNDPTR(:,:,L),       work2d,   grid, rc=status)
      endif
      if(status/=0) then
        ierr = 99
        call die('sppt_pattern','error',ierr)     
      endif
    enddo
    deallocate(work2d)

    if (cubed) then
      call L2C%regrid(RNDPTR_local, RNDPTR)
      deallocate(RNDPTR_local)
    endif 

    if (MAPL_am_I_root()) then 
      deallocate(RNDPTR_WORLD)
    endif 

    !open (lu,file=trim(fnoutc),form='unformatted',access='sequential',convert='little_endian')
    !  call writout4(lu,RNDPTR,Grid,JM,IM,LM)
    !close(lu)

  end subroutine sppt_pattern
!-----------------------------------------------------------------------------------
  subroutine skeb_pattern(CF,GRID,SKEBU_WT,SKEBV_WT,PREF,IM,JM,LM,DT)

    use initialize
    use m_die, only: die
    use stoch_data

    implicit none

    character(len=ESMF_MAXSTR) :: COMP_NAME
!
    type(ESMF_Config),intent(inout)  :: CF
    type(ESMF_Grid)  ,intent(inout)  :: GRID
    real*4,dimension(IM,JM,LM), intent(inout) :: SKEBU_WT, SKEBV_WT
    real*4, intent(in)    :: PREF(LM+1)
    real*4, intent(in)    :: DT 
    real*8, allocatable   :: u_world(:,:,:),v_world(:,:,:)
    real*4, allocatable   :: u_local(:,:,:),v_local(:,:,:)
    real*4, allocatable   :: worku2d(:,:),workv2d(:,:)

    integer, intent(in)   :: IM,JM,LM
    integer :: DIMS(ESMF_MAXGRIDDIM)
    integer :: NX,NY,MEMID

    integer :: STATUS,ierr,rc
    integer :: IM_world,JM_world,LM_world
    integer :: im_gg,jm_gg
    integer :: im_ll_local, jm_ll_local
    integer :: L,myid

    if (.not. isallset) call setup_pattern (CF,GRID)
    if (isk .ge. ishutoff -1) return 
    isk = isk + 1

    if (MAPL_am_I_root()) then
      allocate(u_world(im_ll,jm_ll,lm))
      allocate(v_world(im_ll,jm_ll,lm))
      u_world=0.
      v_world=0.

      call get_skeb(u_world,v_world,PREF)
    endif 

    ! Convert back to original grid
    if (cubed) then 
      call MAPL_GridGet(LLgrid, localCellCountPerDim=DIMS, RC=STATUS)
      im_ll_local=dims(1)
      jm_ll_local=dims(2)
      allocate(u_local(im_ll_local,jm_ll_local,lm))
      allocate(v_local(im_ll_local,jm_ll_local,lm))
    else
      allocate(u_local(IM,JM,LM))
      allocate(v_local(IM,JM,LM))
    endif

    allocate(worku2d(im_ll,jm_ll))
    allocate(workv2d(im_ll,jm_ll))

    worku2d = 0. ; workv2d = 0.   
    u_local = 0. ; v_local = 0. 

    do L=1,lm
      if (MAPL_am_I_root()) then
         worku2d=u_world(:,:,L)
         workv2d=v_world(:,:,L)
      endif 
      if(cubed) then
         call ArrayScatter(u_local(:,:,L), worku2d, LLgrid, rc=status)
         call ArrayScatter(v_local(:,:,L), workv2d, LLgrid, rc=status)
      else
         call ArrayScatter(u_local(:,:,L), worku2d,   GRID, rc=status)
         call ArrayScatter(v_local(:,:,L), workv2d,   GRID, rc=status)
      endif
      if(status/=0) then
         ierr = 99
         call die('skeb_pattern','error',ierr)     
      endif
    enddo
    deallocate(worku2d,workv2d)

    if (cubed) then
      call L2C%regrid(u_local,v_local,SKEBU_WT,SKEBV_WT,rotate=.false.)
    else
      do L=1,LM
        SKEBU_WT(:,:,L) = u_local(:,:,L)
        SKEBV_WT(:,:,L) = v_local(:,:,L)
      enddo
    endif
    deallocate(u_local,v_local)

    if (MAPL_am_I_root())  deallocate(u_world,v_world)
      
  end subroutine skeb_pattern
!-----------------------------------------------------------------------------------
  subroutine get_sppt(sppt_wt,PREF)

    use mod_param
    use stoch_data
    use get_pattern
    use initialize, only: ls_node,ls_nodes,max_ls_nodes,lats_nodes_r, & 
                          global_lats_r,lonsperlar, plnev_r,plnod_r,  & 
                          lonr,lats_node_r,colrad_r,latr2,nlon,nlat,levs

    implicit none
 
    real*4, intent(in)         :: PREF(levs+1)
    real(kind=kind_evod),allocatable,dimension(:,:)   :: sppt_2d
    real(kind=kind_evod)       :: sppt_wt(nlon,nlat)
    character(len=*),parameter :: fnout = 'sppt_ptrn.grd'
    integer,parameter          :: lu=20

    allocate(sppt_2d(lonr,lats_node_r))
    sppt_2d = 0.
    call get_pattern_sppt(sppt_2d,PREF,ls_node,ls_nodes,max_ls_nodes, &
                          lats_nodes_r,global_lats_r,lonsperlar,      &
                          plnev_r,plnod_r)                            

    ! Convert from Gaussian to regular grid
    sppt_wt  = 0.
    call gauss2grid(sppt_2d,lonr,lats_node_r,colrad_r,latr2,sppt_wt,nlon,nlat)

 !  open (lu,file=trim(fnout),form='unformatted',access='sequential',convert='little_endian')
 !    call writout8(lu,sppt_wt,nlon,nlat,1)
 !  close(lu)

    deallocate(sppt_2d)

  end subroutine get_sppt
!-----------------------------------------------------------------------------------
  subroutine get_skeb(skebu_3d,skebv_3d,PREF)
 
    use mod_param
    use stoch_data
    use grd_xform
    use get_pattern
    use initialize, only: ls_node,ls_nodes,max_ls_nodes,lats_nodes_r,     & 
                          global_lats_r,lonsperlar, plnev_r,plnod_r,      & 
                          lonr,lats_node_r,colrad_r,latr2,nlon,nlat,levs, &
                          plnew_r,plnow_r,plnev_r,plnod_r,lonll,latll,    &
                          epsedn,epsodn,snnp1ev,snnp1od

    use m_ggGradientSP,only : ggGradientSP
    use m_ggGradientSP,only : ggGradientSP_init,clean
    use m_ggGradientSP,only : ggDivo

    implicit none
 
    integer              iprint,locl,node,i,k,ie,io
    integer              n,l,j

    real(kind=kind_evod),allocatable,dimension(:,:,:) :: pertu_3d,pertv_3d
    real(kind=kind_evod),dimension(nlon,nlat,levs)    :: skebu_3d,skebv_3d
    real(kind=kind_evod),dimension(lonr,latr,levs)    :: uur,vvr
    real(kind=kind_evod),dimension(lonll,latll,levs)  :: ull,vll
    real*4, intent(in)    :: PREF(levs+1)

    type(ggGradientSP) :: gr
    real(kind=kind_evod),allocatable,dimension(:,:,:)::div,vor

    character(len=*),parameter :: fnout = 'uv_ll_out.grd'
    character(len=*),parameter :: fnke  = 'ke.grd'
    character(len=40) :: fname
    character(len=8) :: indx
    integer,parameter          :: lu=20

    !-------------------------------------------------
    uur =0. ; vvr =0.
    call get_pattern_skeb(uur,vvr,PREF,                          & 
                          ls_node,ls_nodes,max_ls_nodes,         &
                          lats_nodes_r,global_lats_r,lonsperlar, & 
                          epsedn,epsodn,snnp1ev,snnp1od,         & 
                          plnew_r,plnow_r,plnev_r,plnod_r)

    skebu_3d  = 0.
    skebv_3d  = 0.
    call gs2gd(uur,lonr,latr,levs,skebu_3d,lonll,latll)
    call gs2gd(vvr,lonr,latr,levs,skebv_3d,lonll,latll)
    call poleuv(skebu_3d,skebv_3d,lonll,latll,levs)

  end subroutine get_skeb
!-----------------------------------------------------------------------------------
  subroutine writout8 (lu,fld,im,jm,km)
    implicit none

    integer, intent(in) :: lu,im,jm,km
    real(kind=kind_evod),    intent(in) :: fld(jm,im,km)
    integer  k
    real(4),allocatable,dimension(:,:):: fout

    allocate(fout(jm,im))
    do k=1,km
       fout(:,:)=fld(:,:,k) 
       write(lu) fout
    enddo
    deallocate(fout)

  end subroutine writout8
!-----------------------------------------------------------------------------------
  subroutine writout4 (lu,fld,Grid,im,jm,km)
    use ESMF

    implicit none
    integer, intent(in) :: lu,im,jm,km
    real(kind=4),    intent(in) :: fld(jm,im,km)
    type(ESMF_GRID) :: Grid
    integer  k,status,dims(3),imglo,jmglo
    real(4),allocatable,dimension(:,:):: fout

    call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
    jmglo=dims(1)
    imglo=dims(2)
    allocate(fout(jmglo,imglo))
    do k=1,km
       call ArrayGather(fld(:,:,k), fout, Grid, rc=status)
       if (MAPL_am_I_root()) write(lu) fout
    enddo
    deallocate(fout)

  end subroutine writout4
!-----------------------------------------------------------------------------------
end module stoch_module

