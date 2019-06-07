!=======================================================================
!
!BOP
!
! !MODULE: ice_work - globally accessible, temporary work arrays
!
! !DESCRIPTION:
!
! Declare globally accessible, temporary work arrays to conserve memory.
! These arrays should be used only within a single subroutine!
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_work
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_domain_size
!
!EOP
!
      implicit none

      ! global

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, &
         work_g2, &
         work_g3

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      real (kind=real_kind), dimension(:,:,:), allocatable :: &
         work_gr3

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      ! all blocks
      real (kind=dbl_kind), dimension (:,:,:) , allocatable:: &
         work1, &
         work2

      ! local (single block)
      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         worka, &
         workb, &
         workc, &
         workd

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_work - initialize work arrays
!
! !DESCRIPTION:
!
! Initialize work arrays
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_work
!
! !USES:
!
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      work1(:,:,:) = c0
      work2(:,:,:) = c0

      worka(:,:) = c0
      workb(:,:) = c0
      workc(:,:) = c0
      workd(:,:) = c0

      end subroutine init_work

!=======================================================================
!BOP
!
! !ROUTINE: alloc_work - allocate work arrays
!
! !DESCRIPTION:
!
! allocate work arrays
!
! !REVISION HISTORY:
!
! author: Matthew A. Thompson
!
! !INTERFACE:
!
      subroutine alloc_work
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      allocate(work1(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)
      allocate(work2(nx_block,ny_block,max_blocks), source=0.0_dbl_kind)

      allocate(worka(nx_block,ny_block), source=0.0_dbl_kind)
      allocate(workb(nx_block,ny_block), source=0.0_dbl_kind)
      allocate(workc(nx_block,ny_block), source=0.0_dbl_kind)
      allocate(workd(nx_block,ny_block), source=0.0_dbl_kind)

      end subroutine alloc_work

!=======================================================================

!=======================================================================
!BOP
!
! !ROUTINE: dealloc_work - deallocate work arrays
!
! !DESCRIPTION:
!
! deallocate work arrays
!
! !REVISION HISTORY:
!
! author: Matthew A. Thompson
!
! !INTERFACE:
!
      subroutine dealloc_work
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      deallocate(work1)
      deallocate(work2)

      deallocate(worka)
      deallocate(workb)
      deallocate(workc)
      deallocate(workd)

      end subroutine dealloc_work

!=======================================================================

      end module ice_work

!=======================================================================
