!=======================================================================
!BOP
!
! !MODULE: ice_domain_size
!
! !DESCRIPTION:
!
! Defines the global domain size and number of categories and layers.
! Code originally based on domain_size.F in POP
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author Elizabeth C. Hunke, LANL
! 2004: Block structure and snow parameters added by William Lipscomb
!       Renamed (used to be ice_model_size)
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!       Removed hardwired sizes (NX...can now be set in compile scripts)
!
! !INTERFACE:
!
      module ice_domain_size
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!=======================================================================

      implicit none
      save
#ifdef GEOS
      integer (kind=int_kind) :: &
        nx_global , & ! i-axis size
        ny_global     ! j-axis size

      integer (kind=int_kind) :: &
        ncat ,  & ! number of categories
        nilyr,  & ! number of ice layers per category
        ntilyr, & ! number of ice layers in all categories
        nslyr,  & ! number of snow layers per category
        ntslyr, & ! number of snow layers in all categories
        ntrcr     ! number of tracers (defined in ice_state)
                  ! 1 = surface temperature

      integer (kind=int_kind)   :: &
        block_size_x, & ! size of block in first horiz dimension
        block_size_y    ! size of block in second horiz dimension
#else
      integer (kind=int_kind), parameter :: &
        nx_global = NXGLOB    , & ! i-axis size
        ny_global = NYGLOB        ! j-axis size

      integer (kind=int_kind), parameter :: &
        ncat      = n_ice_categories, & ! number of categories
        nilyr     = n_ice_layers, & ! number of ice layers per category
        ntilyr    = ncat*nilyr, & ! number of ice layers in all categories
        nslyr     =   1       , & ! number of snow layers per category
        ntslyr    = ncat*nslyr, & ! number of snow layers in all categories
        ntrcr     =   3           ! number of tracers (defined in ice_state)
                                  ! 1 = surface temperature

      integer (kind=int_kind), parameter :: &
        block_size_x = BLCKX  , & ! size of block in first horiz dimension
        block_size_y = BLCKY      ! size of block in second horiz dimension
#endif


   !*** The model will inform the user of the correct
   !*** values for the parameter below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max_blocks = (nx_global/block_size_x)*(ny_global/block_size_y)/
   !***               num_procs

#ifdef GEOS 
      integer (kind=int_kind) :: &
        max_blocks                ! max number of blocks per processor

      integer (int_kind) :: &
         nghost           ! number of ghost cells around each block

      integer (int_kind) :: &! size of block domain in
         nx_block, &         !  x,y dir including ghost
         ny_block            !  cells 

#else
      integer (kind=int_kind), parameter :: &
        max_blocks = MXBLCKS      ! max number of blocks per processor
#endif

#ifdef GEOS
!=======================================================================

      contains

!=======================================================================
      subroutine init_domain_size(nxg, nyg, blx, bly, nprocs)
         
        integer (kind=int_kind), intent(in) :: &
          nxg, nyg, blx, bly, nprocs

        nx_global = nxg        ! i-axis size
        ny_global = nyg        ! j-axis size
        block_size_x = (nxg/blx)    ! size of block in first horiz dimension
        block_size_y = (nyg/bly)    ! size of block in second horiz dimension

        max_blocks = max( (( (nx_global/block_size_x)*(ny_global/block_size_y) ) / nprocs), 1)
                               ! max number of blocks per processor

        nghost = 1             ! number of ghost cells around each block

        nx_block = block_size_x + 2*nghost     !  x,y dir including ghost
        ny_block = block_size_y + 2*nghost     !  cells 

        !write (*,*) "nx_global: ", nx_global
        !write (*,*) "ny_global: ", ny_global
        !write (*,*) "block_size_x: ", block_size_x
        !write (*,*) "block_size_y: ", block_size_y

        !write (*,*) "max_blocks: ", max_blocks

        !write (*,*) "nghost: ", nghost

        !write (*,*) "nx_block: ", nx_block
        !write (*,*) "ny_block: ", ny_block

      end subroutine init_domain_size

      subroutine init_column_physics(n_ice_cat, n_ice_lay)
        integer (kind=int_kind), intent(in) :: n_ice_cat, n_ice_lay

        ncat   = n_ice_cat     ! number of categories
        nilyr  = n_ice_lay     ! number of ice layers per category
        ntilyr = ncat*nilyr    ! number of ice layers in all categories
        nslyr  =   1           ! number of snow layers per category
        ntslyr = ncat*nslyr    ! number of snow layers in all categories
        ntrcr  =   3           ! number of tracers (defined in ice_state)
                               ! 1 = surface temperature
      end subroutine init_column_physics
      
#endif

!=======================================================================

      end module ice_domain_size

!=======================================================================
