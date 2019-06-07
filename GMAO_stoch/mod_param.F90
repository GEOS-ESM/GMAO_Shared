  module mod_param

    use ESMF
    implicit none
    save

    include 'mpif.h'

!   ! Machine   ------------------------------------------------------! 
    integer, parameter :: kind_evod = 8, kind_dbl_prec = 8, kind_io4 = 4
    integer, parameter :: kind_real = 8, kind_integer = 4                     
    integer, parameter :: kind_phys = selected_real_kind(13,60) 
    integer, parameter :: kind_qdt_prec = selected_real_kind(30,90)

!   ! Constants ------------------------------------------------------! 
    real(kind=kind_phys),parameter:: con_pi     =3.1415926535897931 ! pi
    real(kind=kind_phys),parameter:: con_rerth  =6.3712e+6 ! radius of earth (m)

!   ! Resolution parameters resol_def --------------------------------!
    integer :: jcap,jcap1,jcap2
    integer :: latg,latg2
    integer :: latr,latr2,latll
    integer :: lonr,lonrx,lonll
    integer :: nlon,nlat,levs

!   ! Global dimensions
    integer ::  im_ll,jm_ll

!   ! Previously in layout1 ------------------------------------------!
    integer ,allocatable :: lat1s_r(:)
    integer :: nodes, me, me_l_0, ls_dim,ls_max_node
    integer :: lats_dim_r, lats_node_r, ipt_lats_node_r
    integer :: len_trie_ls, len_trio_ls, lon_dim_r

!
    ! Previously in Internal -----------------------------------------!
    real(kind=kind_evod) :: delt

!   ! stochastic parameters ------------------------------------------!
    integer      :: ens_mem
    character*20 :: ens_nam
    real(kind=kind_evod), dimension(5) :: sppt,sppt_lscale,sppt_tau, &
                                          skeb,skeb_lscale,skeb_tau 
    real(kind=kind_evod)     :: sppt_sigtop1,sppt_sigtop2, & 
                                skeb_sigtop1,skeb_sigtop2
    real(kind=kind_evod)     :: skeb_diss_smooth
    integer,    dimension(5) :: skeb_vfilt
    integer(8), dimension(5) :: iseed_sppt,iseed_skeb
    integer                  :: skeb_varspect_opt
    logical                  :: sppt_sfclimit,sppt_logit,stochini
    logical                  :: opt_jff

!   ! Coordinate_def -------------------------------------------------!
    real(kind=kind_evod), allocatable :: ak5(:),bk5(:)

    ! MPI-related mpi_def --------------------------------------------!
    integer, parameter :: mpi_r_io =mpi_real4
    integer, parameter :: mpi_r_io_r=mpi_real8
    integer, parameter :: mpi_r_mpi=mpi_real4
    integer, parameter :: mc_comp=mpi_comm_world
    integer, parameter :: kind_mpi=4
    integer :: icolor
    logical :: liope

  end module mod_param
