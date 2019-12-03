module initialize

  use mod_param

  IMPLICIT NONE
  
  integer,allocatable,dimension(:,:) :: ls_node,ls_nodes
  integer,allocatable,dimension(:)   :: lonsperlar,global_lats_r,lats_nodes_r
  integer,allocatable,dimension(:)   :: max_ls_nodes

! ! Spherical harmonics
  real(kind=kind_evod),allocatable ::  epse(:),epso(:),epsedn(:),epsodn(:)
  real(kind=kind_evod),allocatable ::  plnev_r(:,:),plnod_r(:,:),plnew_r(:,:),plnow_r(:,:)
  real(kind=kind_evod),allocatable ::  pddev_r(:,:),pddod_r(:,:) 
  real(kind=kind_evod),allocatable ::  snnp1ev(:),snnp1od(:)
! ! Gaussin Grid gg_def --------------------------------------------!
  real(kind=kind_evod), allocatable, dimension(:) :: colrad_r,sinlat_r,coslat_r
  real(kind=kind_evod), allocatable, dimension(:) :: wgt_r,wgtcs_r,rcs2_r
  real(kind=kind_evod), allocatable, dimension(:) :: colrad_a,sinlat_a,wgt_a,wgtcs_a

  public sub_init, sub_create ,sub_setup, sub_clear

  CONTAINS

  !==================================================================!
  subroutine sub_init(im_gg,jm_gg,im_ll,jm_ll,lm,memid,dt)

    implicit none

    real*8,  intent(in) :: dt
    integer, intent(in) :: im_gg,jm_gg,im_ll,jm_ll,lm,memid
    integer, parameter   :: lun=222
    integer ierr
    logical exist

    ! Initialize grid -----------------------------------!
    delt  = dt
    nodes = 1
    levs  = lm
    nlat  = jm_ll ; nlon  = im_ll
    latll = jm_ll ; lonll = im_ll
    latr  = jm_gg ; lonr  = im_gg
    latg  = latr
    latr2 = latr/2 ; latg2 = latg/2
    lonrx = lonr+2
 
    jcap1 = jcap+1

    len_trie_ls = (jcap+1)*(jcap+2)/2
    len_trio_ls = (jcap+1)*(jcap+2)/2
 
    ls_dim = (jcap1-1)/nodes+1
    ls_max_node = ls_dim
    ipt_lats_node_r = 1

    ! Ensemble ID ---this needs to change-----------------------------!  
    ens_mem = memid
   
  end subroutine sub_init
  !==================================================================!
  subroutine sub_create

    implicit none
    ! Allocate arrays 

    allocate(epse(len_trie_ls),epso(len_trio_ls),epsedn(len_trie_ls),epsodn(len_trio_ls))
    allocate(snnp1ev(len_trie_ls),snnp1od(len_trio_ls))
    allocate(plnev_r(len_trie_ls,latr2),plnod_r(len_trio_ls,latr2))
    allocate(pddev_r(len_trie_ls,latr2),pddod_r(len_trio_ls,latr2))
    allocate(plnew_r(len_trie_ls,latr2),plnow_r(len_trio_ls,latr2))

    allocate(lonsperlar(latr),global_lats_r(latr),lats_nodes_r(nodes))
    allocate(ls_node(ls_dim,3),ls_nodes(ls_dim,nodes))
    allocate(max_ls_nodes(nodes))
    allocate(colrad_r(latr), sinlat_r(latr), coslat_r(latr))
    allocate(wgt_r(latr2), wgtcs_r(latr2), rcs2_r(latr2))

  end subroutine sub_create
  !==================================================================!
  subroutine sub_setup( )

    implicit none

    integer :: iprint,locl,node,i,k,ie,io
    integer :: lat,n,l,j
    integer :: indev,indod
    integer :: indlsev,jbasev
    integer :: indlsod,jbasod

    real(kind=kind_evod), parameter :: cons0 = 0.d0
    real(kind=kind_evod)            ::  colat1

    include 'function_indlsev.h'
    include 'function_indlsod.h'

  ! Fill in coordinate quantities ----------------------!
    iprint = 0
    epse = 0. ; epso = 0. ; epsedn = 0. ; epsodn=0.
  !------------------------------------------------------------------
    lonsperlar=lonr    
    allocate(lat1s_r(0:jcap))
    do l=0,jcap
      do lat = 1, latr2
        if ( l .le. min(jcap,lonsperlar(lat)/2) ) then
          lat1s_r(l) = lat
          go to 220
        endif
      end do
220   continue
      end do

    do node=1,nodes
       call get_ls_node( node-1, ls_nodes(1,node),max_ls_nodes(node), iprint )
    enddo

    ie=0; io=0
    ls_node(:,:) = 0
    do i =1,jcap+1
      ls_node(i,1) = ls_nodes(i,1)
      ls_node(i,2) = (i-1)*(JCAP+2)/2
      ls_node(i,3) = (i-1)*(JCAP+2)/2+1
    enddo
    call setlats_r(lats_nodes_r,global_lats_r,iprint,lonsperlar)
    lats_dim_r = 0
    do node=1,nodes
      lats_dim_r = max(lats_dim_r,lats_nodes_r(node))
    enddo
    lats_node_r = lats_nodes_r(1)

    call glats(latr2,colrad_r,wgt_r,wgtcs_r,rcs2_r,iprint)
    colat1 = colrad_r(1)
    do i=latr2+1,latr
      colrad_r(i) = colrad_r(latr+1-i)
    enddo

    do j=1,latr
      if (j.le.latr2) then
        sinlat_r(j) = cos(colrad_r(j))
      else
        sinlat_r(j) = -cos(colrad_r(j))
      endif
      coslat_r(j) = sqrt(1. -sinlat_r(j)*sinlat_r(j))
    enddo
      
    !------------------------------------------------------------------
    plnev_r = 0.;plnod_r =0.
    call epslon(epse,epso,epsedn,epsodn,ls_node)
    call pln2eo_r(plnev_r,plnod_r,epse,epso,colrad_r,ls_node,latr2)
    call gozrineo_r(plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r, &
                    epse,epso,rcs2_r,wgt_r,ls_node,latr2)
    do locl=1,ls_max_node
      l = ls_node(locl,1)
      jbasev = ls_node(locl,2)
      indev  = indlsev(l,l)
      do n = l, jcap, 2
        snnp1ev(indev) = n*(n+1)
        indev        = indev+1
      end do
    end do
!
    do locl=1,ls_max_node
      l = ls_node(locl,1)
      jbasod = ls_node(locl,3)
      if ( l .le. jcap-1 ) then
        indod = indlsod(l+1,l)
        do n = l+1, jcap, 2
          snnp1od(indod) = n*(n+1)
          indod        = indod+1
        end do
      end if
    end do
    do locl=1,ls_max_node
      l = ls_node(locl,1)
      jbasev = ls_node(locl,2)
      jbasod = ls_node(locl,3)
      if (mod(l,2).eq.mod(jcap+1,2)) then
      ! set the even (n-l) terms of the top row to zero
        snnp1ev(indlsev(jcap+1,l)) = cons0     !constant
      else
      ! set the  odd (n-l) terms of the top row to zero
        snnp1od(indlsod(jcap+1,l)) = cons0     !constant
      endif
    enddo

  end subroutine sub_setup
!==================================================================!
  subroutine get_pert_grid(im_world,jm_world,im_ll,jm_ll,jcap,im_gg,jm_gg)

    use ESMF
    use MAPL
    use m_die, only: die

    implicit none
    integer,intent(in)  :: im_world,jm_world
    integer,intent(out) :: im_ll,jm_ll,im_gg,jm_gg
    integer,intent(out) :: jcap
    integer             :: ierr
    logical             :: cubed

    cubed=.false.
    if (im_world*6==jm_world) cubed=.true.

    if (cubed) then
      im_ll=-1
      jm_ll=-1
      if (im_world==48) then
        im_ll=144
        jm_ll=91
      elseif (im_world==90) then
        im_ll=288
        jm_ll=181
      elseif (im_world==180) then
        im_ll=576
        jm_ll=361
      elseif (im_world==360) then
        im_ll=1152
        jm_ll=721
      elseif (im_world==720) then
        im_ll=2304
        jm_ll=1441
      endif
      if(im_ll<0.or.jm_ll<0) then
        print*,' NOT a known grid'
        ierr = 10
        call die('get_pert_grid','error',ierr)
      endif
    else   ! if not cubed
      im_ll=im_world
      jm_ll=jm_world
    endif  ! End if cubed

    if(im_ll==2304.and.jm_ll==1441) then
      jcap =878
      im_gg=2304
      jm_gg=1152
    endif
    if(im_ll==1152.and.jm_ll==721) then
      jcap =574
      im_gg=1760
      jm_gg=880
    endif
    if(im_ll==576.and.jm_ll==361) then
      jcap =254
      im_gg=768
      jm_gg=384
    endif
    if(im_ll==288.and.jm_ll==181) then
      jcap =126
      im_gg=384
      jm_gg=190
    endif
    if(im_ll==144.and.jm_ll==91) then
      jcap =62
      im_gg=192
      jm_gg=94
    endif

  end subroutine get_pert_grid
!==================================================================!
  subroutine sub_clear

    implicit none

    deallocate(epse,epso,epsedn,epsodn)
    deallocate(snnp1ev,snnp1od)
    deallocate(plnev_r,plnod_r)
    deallocate(pddev_r,pddod_r)
    deallocate(plnew_r,plnow_r)
    deallocate(colrad_r,sinlat_r,coslat_r,lat1s_r)
    deallocate(wgt_r,wgtcs_r,rcs2_r)
    deallocate(lonsperlar,global_lats_r,lats_nodes_r)
    deallocate(ls_node,ls_nodes,max_ls_nodes)

  end subroutine sub_clear

end module initialize
