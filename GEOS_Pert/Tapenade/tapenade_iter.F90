module tapenade_iter

! 20170413 | D.Holdaway | Tool for customised interfacing with Tapenade checkpointing in an iterative environment
! 20171010 | D.Holdaway | Update to use derived types and generalize for access by multiple modules


!Functionality this module (will) provide
!
! - Hold checkpoints in memory accross multiple iterations.
! - Compile time precision choice for real variables.
! - Profile a subroutine to see how much memory it uses (coming soon).
!
!-------------------------------------------------------------------------------------------------------------
!               | START             | FORWARD PHASE                  | BACKWARD PHASE                        |
!-------------------------------------------------------------------------------------------------------------
! Iteration  0: | N/A               | Do as Tapenade would do        | Do as Tapenade would do               |
! Iteration  1: | Alloc counter     | Count checkpoints and bytes    | Count down sanity check               |
! Iteration  2: | Alloc do cp array | Count checkpoints and bytes    | Count down and check if cp was needed |
! Iteration  3: | Alloc cp array    | Save cp to array if needed     | Read necessary cps from array         |
! Iteration >3: |                   | OFF                            | Read necessary cps from array         |
!-------------------------------------------------------------------------------------------------------------

#ifdef SPMD
 use MPI
#endif

 implicit none
 private

!Derived type for holding all the controls for the iterative use of checkpoitning
!--------------------------------------------------------------------------------
 type cp_iter_controls_type
   integer :: cp_nt           !Number of time steps
   integer :: cp_i            !Iteration number
   integer :: cp_t            !Current time step
   real    :: cp_gb           !User estimated number of gigabytes avaialable per processor
   integer :: cp_nm           !Number of modules using this routine
 end type

 !Global for all modules using custom checkpointing
 type(cp_iter_controls_type), target :: cp_iter_controls

 !Index in arrays for different types of checkpointed variables
 integer, parameter :: idx_control = 1
 integer, parameter :: idx_integer = 2
 integer, parameter :: idx_real_r4 = 3
 integer, parameter :: idx_real_r8 = 4
 integer, parameter :: total_types = 4 !Total variable types

 integer, parameter :: status_kind = 1

!Derived type for holding all the checkpoints and 
!information about checkpoints for a given module
!------------------------------------------------
 type cp_iter_type

   !User provided
   character(len=3) :: my_name            !User defined name for this module
   logical :: cp_test                     !Testing mode
   logical :: cp_rep                      !Write reports on checkpointing
   logical :: cp_nscp                     !Checkpoint all n-splits?
   logical :: cp_adm_pp                   !Make sure _adm calls do not use custom checkpointing
   logical :: check_st_control            !Check whether control push/pops are necessary
   logical :: check_st_integer            !Check whether control push/pops are necessary
   logical :: check_st_real_r4            !Check whether control push/pops are necessary
   logical :: check_st_real_r8            !Check whether control push/pops are necessary
   integer :: test_dim_st_control         !Dimension of arrays in test mode
   integer :: test_dim_st_integer
   integer :: test_dim_st_real_r4
   integer :: test_dim_st_real_r8
   integer :: test_dim_cp_control
   integer :: test_dim_cp_integer
   integer :: test_dim_cp_real_r4
   integer :: test_dim_cp_real_r8

   !Counts for both number of checkpoints and dimension of all checkpoints
   integer, allocatable :: count_psh(:,:)
   integer, allocatable :: count_pop(:,:)
   integer, allocatable :: index_psh(:,:)
   integer, allocatable :: index_pop(:,:)

   !Arrays for holding status of checkpoints
   integer(kind=status_kind), allocatable :: st_control(:,:)
   integer(kind=status_kind), allocatable :: st_integer(:,:)
   integer(kind=status_kind), allocatable :: st_real_r4(:,:)
   integer(kind=status_kind), allocatable :: st_real_r8(:,:)

   !Arrays for holding checkpoints
   integer, allocatable :: cp_control(:,:)
   integer, allocatable :: cp_integer(:,:)
   real(4), allocatable :: cp_real_r4(:,:)
   real(8), allocatable :: cp_real_r8(:,:)
 
   !Ability to turn recording on and off for stages of loops
   logical :: recording = .true.
   logical :: loop_last_step = .false.

   integer, allocatable :: count_pop_start(:,:)
   integer, allocatable :: index_pop_start(:,:)

 end type cp_iter_type

 !All the checkpoints are contained here
 logical, save :: cp_iter_initialized = .false.
 logical, save :: cp_iter_finalized = .false.
 type(cp_iter_type), allocatable, target :: cp_iter(:)

 !Pointer to active module checkpointing
 type(cp_iter_type), pointer :: am

 !Bytes to mb and gb, RAM so 1024^x convention
 real(8), parameter :: b2mb = 9.536743164062500d-7  !1024.0**-2
 real(8), parameter :: b2gb = 9.313225746154785d-10 !1024.0**-3

 !Counter checking that no push attempts in backwards phase
 integer(8) :: count_psh_mid(total_types)

 !Convenience pointers
 integer, pointer :: cp_nt, cp_i, cp_t, cp_nm
 real, pointer :: cp_gb

#ifdef SPMD
 integer :: MPIERR
#endif
 logical, save :: root_pe = .false.

!Publicize routines and variables
!--------------------------------

 !Checkpointing controls
 public cp_iter_controls, cp_iter

 !Subroutines for controls and testing
 public initialize_cp_iter, finalize_cp_iter, cp_mod_ini, cp_mod_mid, cp_mod_end

 !Checkpoint interfaces for fwd/bwd
 public pushcontrol, popcontrol, &
        pushinteger, popinteger, &
        pushrealarray, poprealarray

 !Checkpoint interfaces for adm
 public pushrealarray_adm, poprealarray_adm

!Interfaces for the fwd/bwd/adm checkpointing subroutines
!--------------------------------------------------------
 
 !Checkpointing used in the fwd/bwd routines
 interface pushinteger
   module procedure psh_integer_k0
   module procedure psh_integer_k1
 end interface

 interface popinteger
   module procedure pop_integer_k0
   module procedure pop_integer_k1
 end interface

 interface pushrealarray
   module procedure psh_real_r4_k0
   module procedure psh_real_r4_k1
   module procedure psh_real_r4_k2
   module procedure psh_real_r4_k3
   module procedure psh_real_r4_k4
   module procedure psh_real_r8_k0
   module procedure psh_real_r8_k1
   module procedure psh_real_r8_k2
   module procedure psh_real_r8_k3
   module procedure psh_real_r8_k4
 end interface

 interface poprealarray
   module procedure pop_real_r4_k0
   module procedure pop_real_r4_k1
   module procedure pop_real_r4_k2
   module procedure pop_real_r4_k3
   module procedure pop_real_r4_k4
   module procedure pop_real_r8_k0
   module procedure pop_real_r8_k1
   module procedure pop_real_r8_k2
   module procedure pop_real_r8_k3
   module procedure pop_real_r8_k4
 end interface

 !Checkpointing used in the adm routines, included only to provide compile time precision choice
 interface pushrealarray_adm
   module procedure psh_adm_real_r4_k0
   module procedure psh_adm_real_r4_k1
   module procedure psh_adm_real_r4_k2
   module procedure psh_adm_real_r4_k3
   module procedure psh_adm_real_r4_k4
   module procedure psh_adm_real_r8_k0
   module procedure psh_adm_real_r8_k1
   module procedure psh_adm_real_r8_k2
   module procedure psh_adm_real_r8_k3
   module procedure psh_adm_real_r8_k4
 end interface

 interface poprealarray_adm
   module procedure pop_adm_real_r4_k0
   module procedure pop_adm_real_r4_k1
   module procedure pop_adm_real_r4_k2
   module procedure pop_adm_real_r4_k3
   module procedure pop_adm_real_r4_k4
   module procedure pop_adm_real_r8_k0
   module procedure pop_adm_real_r8_k1
   module procedure pop_adm_real_r8_k2
   module procedure pop_adm_real_r8_k3
   module procedure pop_adm_real_r8_k4
 end interface

contains


! Global initialize, called once
! ------------------------------

 subroutine initialize_cp_iter

#ifdef SPMD
 logical :: IS_MPI_INIT
 logical :: MPI_INIT_HERE
 integer :: pe_id
#endif

  if (.not. cp_iter_initialized) then

     !Set convenience pointers
     cp_nt => cp_iter_controls%cp_nt
     cp_i  => cp_iter_controls%cp_i 
     cp_t  => cp_iter_controls%cp_t 
     cp_nm => cp_iter_controls%cp_nm
     cp_gb => cp_iter_controls%cp_gb
   
     if (cp_i == 0) return !Not using the tool
   
     allocate(cp_iter(cp_nm))
   
     !Intiailize MPI
#ifdef SPMD
     !Check if initialized
     call MPI_Initialized(IS_MPI_INIT,MPIERR)
    
     !If MPI not itnialized then do that here
     MPI_INIT_HERE = .false.
     if (.not. IS_MPI_INIT) then
        call MPI_init(MPIERR)
        IS_MPI_INIT = .true.
        MPI_INIT_HERE = .true.
     endif
   
     call MPI_Comm_rank(MPI_COMM_WORLD, pe_id, MPIERR)
     if (pe_id == 0) then
        root_pe = .true.
     endif
     if (root_pe .and. MPI_INIT_HERE) write(*,*) 'MPI initialized by iterative checkpointing tool'
     if (root_pe .and. .not.MPI_INIT_HERE) write(*,*) 'MPI already initialized'
#else
     root_pe = .true.
#endif
   
     cp_iter_initialized = .true.

  endif

 end subroutine initialize_cp_iter

 subroutine finalize_cp_iter

  if (.not. cp_iter_finalized) then

     nullify(cp_nt)
     nullify(cp_i )
     nullify(cp_t )
     nullify(cp_nm)
     nullify(cp_gb)

     if (allocated(cp_iter)) deallocate(cp_iter)

     cp_iter_finalized = .true.

  endif

 end subroutine finalize_cp_iter


! Initialize called by each module accessing this tool
! ----------------------------------------------------

 subroutine cp_mod_ini(cp_mod_index)

  implicit none
  integer, intent(in) :: cp_mod_index

  !Call intialize, usually already done
  !------------------------------------
  if (.not.cp_iter_initialized) call initialize_cp_iter


  !Set pointer to the tool wide active calling module
  !--------------------------------------------------
  am => cp_iter(cp_mod_index)


  !Write information at user request
  !---------------------------------
  if (cp_i .ne. 0 .and. am%cp_rep) then
     if ( root_pe ) then
        write(*,"(A)") '  '
        write(*,"(A)") 'Checkpointing information...'
        write(*,"(A,I4,A,I4)") 'Iteration number : ', cp_i
        write(*,"(A,I4,A,I4)") 'Time step in iter: ', cp_t, ' of, ', cp_nt
        write(*,"(A)") ' '
        write(*,"(A)") 'Active module information:'
        write(*,"(A)") 'Name: '//am%my_name
        write(*,"(A,L)") 'Running test mode for module: ', am%cp_test
        write(*,"(A)") '  '
     endif
  endif


  !Iteration specific initializations
  !----------------------------------
  if (cp_t == 1) then

     if (cp_i == 1) then
   
        !Allocate space for push counters
        allocate(am%count_psh(cp_nt,total_types)); am%count_psh = 0
        allocate(am%count_pop(cp_nt,total_types)); am%count_pop = 0
        allocate(am%index_psh(cp_nt,total_types)); am%index_psh = 0
        allocate(am%index_pop(cp_nt,total_types)); am%index_pop = 0
   
     elseif (cp_i == 2) then
   
        if (am%cp_test) then
           !Not allocated as running only this iteration
           allocate(am%count_psh(cp_nt,total_types)); am%count_psh = 0
           allocate(am%count_pop(cp_nt,total_types)); am%count_pop = 0
           allocate(am%index_psh(cp_nt,total_types)); am%index_psh = 0
           allocate(am%index_pop(cp_nt,total_types)); am%index_pop = 0
           !Maximum over all time steps and processors, from previous iter and user provided
           am%count_psh(:,idx_control) = am%test_dim_st_control
           am%count_psh(:,idx_integer) = am%test_dim_st_integer
           am%count_psh(:,idx_real_r4) = am%test_dim_st_real_r4
           am%count_psh(:,idx_real_r8) = am%test_dim_st_real_r8
        endif
      
        !Allocate space for checkpoint status
        if (am%check_st_control) then
           allocate(am%st_control(cp_nt,maxval(am%count_psh(:,idx_control))))
           am%st_control = 1
        endif
        if (am%check_st_integer) then
           allocate(am%st_integer(cp_nt,maxval(am%count_psh(:,idx_integer))))
           am%st_integer = 1
        endif
        if (am%check_st_real_r4) then
           allocate(am%st_real_r4(cp_nt,maxval(am%count_psh(:,idx_real_r4))))
           am%st_real_r4 = 1
        endif
        if (am%check_st_real_r8) then
           allocate(am%st_real_r8(cp_nt,maxval(am%count_psh(:,idx_real_r8))))
           am%st_real_r8 = 1
        endif
   
        !Once allocated reset counters
        am%count_psh = 0
        am%index_psh = 0

     elseif (cp_i == 3) then

        if (am%cp_test) then
           !Not allocated as running only this iteration
           allocate(am%count_psh(cp_nt,total_types)); am%count_psh = 0
           allocate(am%count_pop(cp_nt,total_types)); am%count_pop = 0
           allocate(am%index_psh(cp_nt,total_types)); am%index_psh = 0
           allocate(am%index_pop(cp_nt,total_types)); am%index_pop = 0
           !Maximum over all time steps and processors, from previous iter and user provided
           am%count_psh(:,idx_control) = am%test_dim_st_control
           am%count_psh(:,idx_integer) = am%test_dim_st_integer
           am%count_psh(:,idx_real_r4) = am%test_dim_st_real_r4
           am%count_psh(:,idx_real_r8) = am%test_dim_st_real_r8
           !Maximum over all time steps and processors, from previous iter and user provided
           am%index_psh(:,idx_control) = am%test_dim_cp_control
           am%index_psh(:,idx_integer) = am%test_dim_cp_integer
           am%index_psh(:,idx_real_r4) = am%test_dim_cp_real_r4
           am%index_psh(:,idx_real_r8) = am%test_dim_cp_real_r8
        endif
   
        allocate(am%count_pop_start(cp_nt,total_types)); am%count_pop_start = 0
        allocate(am%index_pop_start(cp_nt,total_types)); am%index_pop_start = 0

        !Allocate space for checkpoints
        allocate(am%cp_control(cp_nt,maxval(am%index_pop(:,idx_control)))); am%cp_control = 0
        allocate(am%cp_integer(cp_nt,maxval(am%index_pop(:,idx_integer)))); am%cp_integer = 0
        allocate(am%cp_real_r4(cp_nt,maxval(am%index_pop(:,idx_real_r4)))); am%cp_real_r4 = 0.0_4
        allocate(am%cp_real_r8(cp_nt,maxval(am%index_pop(:,idx_real_r8)))); am%cp_real_r8 = 0.0_8
   
        !Once allocated reset counters
        am%count_psh = 0
        am%index_psh = 0
      
     endif

  endif

 end subroutine cp_mod_ini


! Mid point after fwd and before bwd parts of the adjoint
! -------------------------------------------------------

 subroutine cp_mod_mid

  implicit none

#ifdef SPMD
  integer(8) :: count_tmp(total_types)
  integer :: types
#endif

 
  if (cp_i > 1) then

     !Reset countdowns
     am%count_pop(cp_t,:) = am%count_psh(cp_t,:)
     am%index_pop(cp_t,:) = am%index_psh(cp_t,:)

  elseif (cp_i == 1) then

     !Counting up for testing so set to zero
     am%count_pop(cp_t,:) = 0

     count_psh_mid(idx_control) = int(am%count_psh(cp_t,idx_control),8)
     count_psh_mid(idx_integer) = int(am%count_psh(cp_t,idx_integer),8)
     count_psh_mid(idx_real_r4) = int(am%count_psh(cp_t,idx_real_r4),8)
     count_psh_mid(idx_real_r8) = int(am%count_psh(cp_t,idx_real_r8),8)

#ifdef SPMD
     do types = 1,total_types
        call MPI_ALLREDUCE(count_psh_mid(types),count_tmp(types),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)     
        count_psh_mid(types) = count_tmp(types)
     enddo
#endif

  endif 

  if (cp_i == 3) then

     am%count_pop_start(cp_t,:) = am%count_psh(cp_t,:)
     am%index_pop_start(cp_t,:) = am%index_psh(cp_t,:)

  endif

  if (cp_i >= 3) then

    am%count_pop(cp_t,:) = am%count_pop_start(cp_t,:)
    am%index_pop(cp_t,:) = am%index_pop_start(cp_t,:)

  endif


 end subroutine cp_mod_mid


!End of module information output and clean up
!---------------------------------------------

 subroutine cp_mod_end

  implicit none

  integer(8) :: count_psh(total_types)
  integer(8) :: index_psh(total_types)
  integer(8) :: count_pop(total_types)

#ifdef SPMD
  integer(8) :: count_tmp(total_types)
  integer(8) :: index_tmp(total_types)
#endif
 
  real(8) :: memuse_pe_st(total_types+1,2)
  real(8) :: memuse_pe_cp(total_types+1,2)
  real(8) :: memuse_mn_st(total_types+1,2)
  real(8) :: memuse_mn_cp(total_types+1,2)
  real(8) :: memuse_mx_st(total_types+1,2)
  real(8) :: memuse_mx_cp(total_types+1,2)
  real(8) :: memuse_su_st(total_types+1,2)
  real(8) :: memuse_su_cp(total_types+1,2)
  real(8) :: memsav_pe(total_types+1)
  real(8) :: memsav_mn(total_types+1)
  real(8) :: memsav_mx(total_types+1)
  real(8) :: memsav_su(total_types+1)
  integer :: types 

  character(len=19) :: int2char

  if (cp_i == 1) then

     !Sometimes Tapenade will make use of push to stack during the backward
     !sweep of the code. This would not be compatible with this utility
     !since the indexing of the arrays that hold the checkpoints would be
     !off. It could be gotten around by creating another stack to hold these
     !checkpoints in, with intenal adjusting of the pointer to the current stack.
     !However it is often easiest to just modify the Tapenade code. If the below
     !checks fail the utility will revert to recompting the forward sweep and
     !will not make any attempt to write the reference state to static memory.

     count_psh(idx_control) = int(am%count_psh(cp_t,idx_control),8)
     count_psh(idx_integer) = int(am%count_psh(cp_t,idx_integer),8)
     count_psh(idx_real_r4) = int(am%count_psh(cp_t,idx_real_r4),8)
     count_psh(idx_real_r8) = int(am%count_psh(cp_t,idx_real_r8),8)

#ifdef SPMD
     do types = 1,total_types
        call MPI_ALLREDUCE(count_psh(types),count_tmp(types),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)     
        count_psh(types) = count_tmp(types)
     enddo
#endif

     if (maxval(count_psh - count_psh_mid) .ne. 0) then

        if (root_pe) then
           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' !!!!!!!!!!!! WARNING !!!!!!!!!!!'
           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' Attempts to push during backward'
           write(*,"(A)"  ) ' sweep detected, this utility    '
           write(*,"(A)"  ) ' will not work so reverting to   '
           write(*,"(A)"  ) ' doing recomputation.            '
           write(*,"(A)"  ) ' '
           write(*,"(A,I19)") ' control: ', count_psh(idx_control)-count_psh_mid(idx_control)
           write(*,"(A,I19)") ' integer: ', count_psh(idx_integer)-count_psh_mid(idx_integer)
           write(*,"(A,I19)") ' real_r4: ', count_psh(idx_real_r4)-count_psh_mid(idx_real_r4)
           write(*,"(A,I19)") ' real_r8: ', count_psh(idx_real_r8)-count_psh_mid(idx_real_r8)
           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        endif

        cp_i = 0
        deallocate(am%count_psh)
        deallocate(am%count_pop)
        deallocate(am%index_psh)
        deallocate(am%index_pop)

     endif

  endif


  if (cp_i == 1 .and. am%cp_rep) then

     count_psh(idx_control) = int(am%count_psh(cp_t,idx_control),8)
     count_psh(idx_integer) = int(am%count_psh(cp_t,idx_integer),8)
     count_psh(idx_real_r4) = int(am%count_psh(cp_t,idx_real_r4),8)
     count_psh(idx_real_r8) = int(am%count_psh(cp_t,idx_real_r8),8)

     count_pop(idx_control) = int(am%count_pop(cp_t,idx_control),8)
     count_pop(idx_integer) = int(am%count_pop(cp_t,idx_integer),8)
     count_pop(idx_real_r4) = int(am%count_pop(cp_t,idx_real_r4),8)
     count_pop(idx_real_r8) = int(am%count_pop(cp_t,idx_real_r8),8)

#ifdef SPMD
     call MPI_ALLREDUCE(count_psh(idx_control),count_tmp(idx_control),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_integer),count_tmp(idx_integer),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_real_r4),count_tmp(idx_real_r4),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_real_r8),count_tmp(idx_real_r8),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     count_psh(idx_control) = count_tmp(idx_control)
     count_psh(idx_integer) = count_tmp(idx_integer)
     count_psh(idx_real_r4) = count_tmp(idx_real_r4)
     count_psh(idx_real_r8) = count_tmp(idx_real_r8)
     call MPI_ALLREDUCE(count_pop(idx_control),count_tmp(idx_control),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_pop(idx_integer),count_tmp(idx_integer),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_pop(idx_real_r4),count_tmp(idx_real_r4),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_pop(idx_real_r8),count_tmp(idx_real_r8),1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     count_pop(idx_control) = count_tmp(idx_control)
     count_pop(idx_integer) = count_tmp(idx_integer)
     count_pop(idx_real_r4) = count_tmp(idx_real_r4)
     count_pop(idx_real_r8) = count_tmp(idx_real_r8)
#endif

     if (root_pe) then
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) 'COUNTS FORWARD VERSUS COUNTS BACKWARD, SHOULD BE EQUAL '
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) 'Forward:'
        write(*,"(A,I19)") 'control = ', count_psh(idx_control) 
        write(*,"(A,I19)") 'integer = ', count_psh(idx_integer)
        write(*,"(A,I19)") 'real_r4 = ', count_psh(idx_real_r4)
        write(*,"(A,I19)") 'real_r8 = ', count_psh(idx_real_r8)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) 'Backward:'
        write(*,"(A,I19)") 'control = ', count_pop(idx_control) 
        write(*,"(A,I19)") 'integer = ', count_pop(idx_integer)
        write(*,"(A,I19)") 'real_r4 = ', count_pop(idx_real_r4)
        write(*,"(A,I19)") 'real_r8 = ', count_pop(idx_real_r8)
        write(*,"(A)"  ) ' '
     endif

  endif

  if (cp_i == 1 .and. cp_t == cp_nt .and. (am%cp_test .or. am%cp_rep)) then

     count_psh(idx_control) = int(maxval(am%count_psh(:,idx_control)),8)
     count_psh(idx_integer) = int(maxval(am%count_psh(:,idx_integer)),8)
     count_psh(idx_real_r4) = int(maxval(am%count_psh(:,idx_real_r4)),8)
     count_psh(idx_real_r8) = int(maxval(am%count_psh(:,idx_real_r8)),8)

#ifdef SPMD
     call MPI_ALLREDUCE(count_psh(idx_control),count_tmp(idx_control),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_integer),count_tmp(idx_integer),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_real_r4),count_tmp(idx_real_r4),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(count_psh(idx_real_r8),count_tmp(idx_real_r8),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     count_psh(idx_control) = count_tmp(idx_control)
     count_psh(idx_integer) = count_tmp(idx_integer)
     count_psh(idx_real_r4) = count_tmp(idx_real_r4)
     count_psh(idx_real_r8) = count_tmp(idx_real_r8)
#endif

     if (root_pe) then
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) 'Array dimensions for next iteration'
        write(*,"(A)"  ) 'Module ID: '//am%my_name
        write(int2char,"(I19)") count_psh(idx_control)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_st_control: ', adjustl(trim(int2char))
        write(int2char,"(I19)") count_psh(idx_integer)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_st_integer: ', adjustl(trim(int2char))
        write(int2char,"(I19)") count_psh(idx_real_r4)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_st_real_r4: ', adjustl(trim(int2char))
        write(int2char,"(I19)") count_psh(idx_real_r8)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_st_real_r8: ', adjustl(trim(int2char))
        write(*,"(A)"  ) ' '
     endif

  elseif (cp_i == 2 .and. cp_t == cp_nt .and. (am%cp_test .or. am%cp_rep)) then

     index_psh(idx_control) = int(maxval(am%index_psh(:,idx_control)),8)
     index_psh(idx_integer) = int(maxval(am%index_psh(:,idx_integer)),8)
     index_psh(idx_real_r4) = int(maxval(am%index_psh(:,idx_real_r4)),8)
     index_psh(idx_real_r8) = int(maxval(am%index_psh(:,idx_real_r8)),8)

#ifdef SPMD
     call MPI_ALLREDUCE(index_psh(idx_control),index_tmp(idx_control),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(index_psh(idx_integer),index_tmp(idx_integer),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(index_psh(idx_real_r4),index_tmp(idx_real_r4),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     call MPI_ALLREDUCE(index_psh(idx_real_r8),index_tmp(idx_real_r8),1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD,MPIERR)
     index_psh(idx_control) = index_tmp(idx_control)
     index_psh(idx_integer) = index_tmp(idx_integer)
     index_psh(idx_real_r4) = index_tmp(idx_real_r4)
     index_psh(idx_real_r8) = index_tmp(idx_real_r8)
#endif

     !Memory use by each type
     != Positions*precision*bytestogb
     memuse_pe_st(idx_control,1) = real(maxval(am%count_psh(:,idx_control)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_integer,1) = real(maxval(am%count_psh(:,idx_integer)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_real_r4,1) = real(maxval(am%count_psh(:,idx_real_r4)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_real_r8,1) = real(maxval(am%count_psh(:,idx_real_r8)) * status_kind ,8) * b2gb
     memuse_pe_st(total_types+1,1) = sum(memuse_pe_st(1:total_types,1))

     memuse_pe_st(idx_control,2) = real(maxval(am%count_pop(:,idx_control)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_integer,2) = real(maxval(am%count_pop(:,idx_integer)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_real_r4,2) = real(maxval(am%count_pop(:,idx_real_r4)) * status_kind ,8) * b2gb
     memuse_pe_st(idx_real_r8,2) = real(maxval(am%count_pop(:,idx_real_r8)) * status_kind ,8) * b2gb
     memuse_pe_st(total_types+1,2) = sum(memuse_pe_st(1:total_types,1))

     memuse_pe_cp(idx_control,1) = real(maxval(am%index_psh(:,idx_control)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_integer,1) = real(maxval(am%index_psh(:,idx_integer)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_real_r4,1) = real(maxval(am%index_psh(:,idx_real_r4)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_real_r8,1) = real(maxval(am%index_psh(:,idx_real_r8)) * 8 ,8) * b2gb
     memuse_pe_cp(total_types+1,1) = sum(memuse_pe_cp(1:total_types,1))

     memuse_pe_cp(idx_control,2) = real(maxval(am%index_pop(:,idx_control)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_integer,2) = real(maxval(am%index_pop(:,idx_integer)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_real_r4,2) = real(maxval(am%index_pop(:,idx_real_r4)) * 4 ,8) * b2gb
     memuse_pe_cp(idx_real_r8,2) = real(maxval(am%index_pop(:,idx_real_r8)) * 8 ,8) * b2gb
     memuse_pe_cp(total_types+1,2) = sum(memuse_pe_cp(1:total_types,2))

     !Memory that could be saved by checking the status of a checkpoint
     memsav_pe(idx_control) = memuse_pe_cp(idx_control,1) - (memuse_pe_cp(idx_control,2) + memuse_pe_st(idx_control,1))
     memsav_pe(idx_integer) = memuse_pe_cp(idx_integer,1) - (memuse_pe_cp(idx_integer,2) + memuse_pe_st(idx_integer,1))
     memsav_pe(idx_real_r4) = memuse_pe_cp(idx_real_r4,1) - (memuse_pe_cp(idx_real_r4,2) + memuse_pe_st(idx_real_r4,1))
     memsav_pe(idx_real_r8) = memuse_pe_cp(idx_real_r8,1) - (memuse_pe_cp(idx_real_r8,2) + memuse_pe_st(idx_real_r8,1))
     memsav_pe(total_types+1) = sum(memsav_pe(1:total_types))

     if (.not. am%cp_test) then
        !Multiply for all time steps
        memuse_pe_st = memuse_pe_st * real(cp_nt,8)
        memuse_pe_cp = memuse_pe_cp * real(cp_nt,8)
        memsav_pe    = memsav_pe    * real(cp_nt,8)
     endif

#ifdef SPMD
     !Sum, max and min across processors
     do types = 1,total_types+1
        call MPI_ALLREDUCE(memuse_pe_st(types,1),memuse_mn_st(types,1),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,1),memuse_mn_cp(types,1),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_st(types,1),memuse_mx_st(types,1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,1),memuse_mx_cp(types,1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_st(types,1),memuse_su_st(types,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,1),memuse_su_cp(types,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_st(types,2),memuse_mn_st(types,2),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,2),memuse_mn_cp(types,2),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_st(types,2),memuse_mx_st(types,2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,2),memuse_mx_cp(types,2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_st(types,2),memuse_su_st(types,2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memuse_pe_cp(types,2),memuse_su_cp(types,2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memsav_pe(types),memsav_mn(types),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memsav_pe(types),memsav_mx(types),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
        call MPI_ALLREDUCE(memsav_pe(types),memsav_su(types),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIERR)
     enddo

     if (root_pe) then
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' MEMORY REQUIREMENTS FOR REFERENCE STATE PUSH (GB of RAM)'
        write(*,"(A)"  ) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
        write(*,"(A)"  ) ' Minumum across processors (control) '
        write(*,"(E15.7)"  )   memuse_mn_cp(idx_control,1)
        write(*,"(A)"  ) ' Maximum across processors (control) '
        write(*,"(E15.7)"  )   memuse_mx_cp(idx_control,1)
        write(*,"(A)"  ) ' Sum across processors     (control) '
        write(*,"(E15.7)"  )   memuse_su_cp(idx_control,1)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (integer) '
        write(*,"(E15.7)"  )   memuse_mn_cp(idx_integer,1)
        write(*,"(A)"  ) ' Maximum across processors (integer) '
        write(*,"(E15.7)"  )   memuse_mx_cp(idx_integer,1)
        write(*,"(A)"  ) ' Sum across processors     (integer) '
        write(*,"(E15.7)"  )   memuse_su_cp(idx_integer,1)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (real_r4) '
        write(*,"(E15.7)"  )   memuse_mn_cp(idx_real_r4,1)
        write(*,"(A)"  ) ' Maximum across processors (real_r4) '
        write(*,"(E15.7)"  )   memuse_mx_cp(idx_real_r4,1)
        write(*,"(A)"  ) ' Sum across processors     (real_r4) '
        write(*,"(E15.7)"  )   memuse_su_cp(idx_real_r4,1)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (real_r8) '
        write(*,"(E15.7)"  )   memuse_mn_cp(idx_real_r8,1)
        write(*,"(A)"  ) ' Maximum across processors (real_r8) '
        write(*,"(E15.7)"  )   memuse_mx_cp(idx_real_r8,1)
        write(*,"(A)"  ) ' Sum across processors     (real_r8) '
        write(*,"(E15.7)"  )   memuse_su_cp(idx_real_r8,1)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (Total) '
        write(*,"(E15.7)"  )   memuse_mn_cp(total_types+1,1)
        write(*,"(A)"  ) ' Maximum across processors (Total) '
        write(*,"(E15.7)"  )   memuse_mx_cp(total_types+1,1)
        write(*,"(A)"  ) ' Sum across processors     (Total) '
        write(*,"(E15.7)"  )   memuse_su_cp(total_types+1,1)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' MEMORY SAVED BY CHECKING THE STATUS OF CHECKPOINTS (GB of RAM)'
        write(*,"(A)"  ) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
        write(*,"(A)"  ) ' Minumum across processors (control) '
        write(*,"(E15.7)"  )   memsav_mn(idx_control)
        write(*,"(A)"  ) ' Maximum across processors (control) '
        write(*,"(E15.7)"  )   memsav_mx(idx_control)
        write(*,"(A)"  ) ' Sum across processors (control) '
        write(*,"(E15.7)"  )   memsav_su(idx_control)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (integer) '
        write(*,"(E15.7)"  )   memsav_mn(idx_integer)
        write(*,"(A)"  ) ' Maximum across processors (integer) '
        write(*,"(E15.7)"  )   memsav_mx(idx_integer)
        write(*,"(A)"  ) ' Sum across processors (integer) '
        write(*,"(E15.7)"  )   memsav_su(idx_integer)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (real_r4) '
        write(*,"(E15.7)"  )   memsav_mn(idx_real_r4)
        write(*,"(A)"  ) ' Maximum across processors (real_r4) '
        write(*,"(E15.7)"  )   memsav_mx(idx_real_r4)
        write(*,"(A)"  ) ' Sum across processors (real_r4) '
        write(*,"(E15.7)"  )   memsav_su(idx_real_r4)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (real_r8) '
        write(*,"(E15.7)"  )   memsav_mn(idx_real_r8)
        write(*,"(A)"  ) ' Maximum across processors (real_r8) '
        write(*,"(E15.7)"  )   memsav_mx(idx_real_r8)
        write(*,"(A)"  ) ' Sum across processors (real_r8) '
        write(*,"(E15.7)"  )   memsav_su(idx_real_r8)
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) ' Minumum across processors (total) '
        write(*,"(E15.7)"  )   memsav_mn(total_types+1)
        write(*,"(A)"  ) ' Maximum across processors (total) '
        write(*,"(E15.7)"  )   memsav_mx(total_types+1)
        write(*,"(A)"  ) ' Sum across processors (total) '
        write(*,"(E15.7)"  )   memsav_su(total_types+1)
        write(*,"(A)"  ) ' '
     endif

#else

     memuse_mn_st = memuse_pe_st
     memuse_mn_cp = memuse_pe_cp
     memuse_mx_st = memuse_pe_st
     memuse_mx_cp = memuse_pe_cp
     memuse_su_st = memuse_pe_st
     memuse_su_cp = memuse_pe_cp
     memsav_mn = memsav_pe
     memsav_mx = memsav_pe
     memsav_su = memsav_pe

     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' MEMORY REQUIREMENTS FOR REFERENCE STATE (Gb of RAM)'
     write(*,"(A)"  ) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
     write(*,"(A)"  ) ' Control '
     write(*,"(E15.7)"  )   memuse_su_cp(idx_control,1)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' Integer '
     write(*,"(E15.7)"  )   memuse_su_cp(idx_integer,1)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' real_r4 '
     write(*,"(E15.7)"  )   memuse_su_cp(idx_real_r4,1)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' real_r8 '
     write(*,"(E15.7)"  )   memuse_su_cp(idx_real_r8,1)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' Total '
     write(*,"(E15.7)"  )   memuse_su_cp(total_types+1,1)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' MEMORY SAVED BY CHECKING THE STATUS OF CHECKPOINTS (Gb of RAM)'
     write(*,"(A)"  ) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
     write(*,"(A)"  ) ' Control '
     write(*,"(E15.7)"  )   memsav_su(idx_control)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' Integer '
     write(*,"(E15.7)"  )   memsav_su(idx_integer)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' real_r4 '
     write(*,"(E15.7)"  )   memsav_su(idx_real_r4)
     write(*,"(A)"  ) ' '
     write(*,"(A)"  ) ' real_r8 '
     write(*,"(E15.7)"  )   memsav_su(idx_real_r8)
     write(*,"(A)"  ) ' '

#endif

     if (root_pe) then
        write(*,"(A)"  ) ' '
        write(*,"(A)"  ) 'Array dimensions for next iteration'
        write(*,"(A)"  ) 'Module ID: '//am%my_name
        write(int2char,"(I19)") index_psh(idx_control)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_cp_control: ', adjustl(trim(int2char))
        write(int2char,"(I19)") index_psh(idx_integer)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_cp_integer: ', adjustl(trim(int2char))
        write(int2char,"(I19)") index_psh(idx_real_r4)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_cp_real_r4: ', adjustl(trim(int2char))
        write(int2char,"(I19)") index_psh(idx_real_r8)
        write(*,"(A,A)") 'CP_'//am%my_name//'_dim_cp_real_r8: ', adjustl(trim(int2char))
        write(*,"(A)"  ) ' '
     endif

     if (cp_gb > 0) then
        if (memuse_mx_cp(total_types+1,1) > cp_gb .and. root_pe) then

           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' !!!! WARNING !!!!'
           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' Maximum expected memory use by this module is '
           write(*,"(A)"  ) ' greater than the user provided estimation for '
           write(*,"(A)"  ) ' the amount available per processor.           '
           write(*,"(A)"  ) ' '
           write(*,"(A)"  ) ' Possibility of crash at the next iteration.   '
           write(*,"(A)"  ) ' '

        endif
     endif

  endif

  !Nullify the active module pointer
  nullify(am)

 end subroutine cp_mod_end


! pushcontrol
! -----------

 subroutine pushcontrol(ctype,field)

  implicit none

  integer, intent(in) :: ctype,field

  integer(status_kind) :: docp

  if (cp_i == 0) then

     if (ctype == 1) then
        call PUSHCONTROL1B(field)
     elseif (ctype == 2) then
        call PUSHCONTROL2B(field)
     elseif (ctype == 3) then
        call PUSHCONTROL3B(field)
     elseif (ctype == 4) then
        call PUSHCONTROL4B(field)
     elseif (ctype == 5) then
        call PUSHCONTROL5B(field)
     endif

  elseif (cp_i == 1 .or. cp_i == 2) then

     am%count_psh(cp_t,idx_control) = am%count_psh(cp_t,idx_control) + 1
     am%index_psh(cp_t,idx_control) = am%index_psh(cp_t,idx_control) + 1

     if (ctype == 1) then
        call PUSHCONTROL1B(field)
     elseif (ctype == 2) then
        call PUSHCONTROL2B(field)
     elseif (ctype == 3) then
        call PUSHCONTROL3B(field)
     elseif (ctype == 4) then
        call PUSHCONTROL4B(field)
     elseif (ctype == 5) then
        call PUSHCONTROL5B(field)
     endif

  elseif (cp_i == 3) then

     am%count_psh(cp_t,idx_control) = am%count_psh(cp_t,idx_control) + 1

     docp = 1
     if (am%check_st_control) then
        docp = am%st_control(cp_t,am%count_psh(cp_t,idx_control)) 
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_control) = am%index_psh(cp_t,idx_control) + 1
        am%cp_control(cp_t,am%index_psh(cp_t,idx_control)) = field
     endif

  endif

 end subroutine pushcontrol


! popcontrol
! -----------

 subroutine popcontrol(ctype,field)

  implicit none

  integer, intent(in) :: ctype
  integer, intent(inout) :: field
  integer :: tmp
  integer(status_kind) :: docp

  if (cp_i == 0) then

     if (ctype == 1) then
        call POPCONTROL1B(field)
     elseif (ctype == 2) then
        call POPCONTROL2B(field)
     elseif (ctype == 3) then
        call POPCONTROL3B(field)
     elseif (ctype == 4) then
        call POPCONTROL4B(field)
     elseif (ctype == 5) then
        call POPCONTROL5B(field)
     endif

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_control) = am%count_pop(cp_t,idx_control) + 1

     if (am%check_st_control .and. cp_i == 2) tmp = field

     if (ctype == 1) then
        call POPCONTROL1B(field)
     elseif (ctype == 2) then
        call POPCONTROL2B(field)
     elseif (ctype == 3) then
        call POPCONTROL3B(field)
     elseif (ctype == 4) then
        call POPCONTROL4B(field)
     elseif (ctype == 5) then
        call POPCONTROL5B(field)
     endif

     if (am%check_st_control .and. cp_i == 2) then
        if (field - tmp == 0) then
           am%st_control(cp_t,am%count_pop(cp_t,idx_control)) = 0
           am%index_pop(cp_t,idx_control) = am%index_pop(cp_t,idx_control) - 1
        elseif (field == 0) then
           am%st_control(cp_t,am%count_pop(cp_t,idx_control)) = -1
           am%index_pop(cp_t,idx_control) = am%index_pop(cp_t,idx_control) - 1
        endif
        am%count_pop(cp_t,idx_control) = am%count_pop(cp_t,idx_control) - 1
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_control) then
        docp = am%st_control(cp_t,am%count_pop(cp_t,idx_control))
     endif
     am%count_pop(cp_t,idx_control) = am%count_pop(cp_t,idx_control) - 1

     if (docp == 1) then
        field = am%cp_control(cp_t,am%index_pop(cp_t,idx_control))
        am%index_pop(cp_t,idx_control) = am%index_pop(cp_t,idx_control) - 1
     elseif (docp == -1) then
        field = 0
     endif
 
  endif

 end subroutine popcontrol



! pushinteger
! -----------

 subroutine psh_integer_k0(field,skip)

  implicit none

  integer, intent(in) :: field
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHINTEGER4(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_integer) = am%count_psh(cp_t,idx_integer) + 1
        am%index_psh(cp_t,idx_integer) = am%index_psh(cp_t,idx_integer) + 1
     endif

     if(.not.skipcp) CALL PUSHINTEGER4(field)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_integer) = am%count_psh(cp_t,idx_integer) + 1

     docp = 1
     if (am%check_st_integer) then
        docp = am%st_integer(cp_t,am%count_psh(cp_t,idx_integer))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_integer) = am%index_psh(cp_t,idx_integer) + 1
        am%cp_integer(cp_t,am%index_psh(cp_t,idx_integer)) = field
     endif

  endif

 end subroutine psh_integer_k0

 subroutine psh_integer_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  integer, intent(in) :: field(dimen)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHINTEGER4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_integer) = am%count_psh(cp_t,idx_integer) + 1
        am%index_psh(cp_t,idx_integer) = am%index_psh(cp_t,idx_integer) + dimen
     endif

     if(.not.skipcp) CALL PUSHINTEGER4ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_integer) = am%count_psh(cp_t,idx_integer) + 1

     docp = 1
     if (am%check_st_integer) then
        docp = am%st_integer(cp_t,am%count_psh(cp_t,idx_integer))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_integer) = am%index_psh(cp_t,idx_integer) + dimen
        am%cp_integer(cp_t,am%index_psh(cp_t,idx_integer)-dimen+1:am%index_psh(cp_t,idx_integer)) = field
     endif

  endif

 end subroutine psh_integer_k1

! popinteger
! -----------

 subroutine pop_integer_k0(field,skip)

  implicit none

  integer, intent(inout) :: field
  logical, optional, intent(in) :: skip
  integer :: tmp
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPINTEGER4(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) + 1

     if (am%check_st_integer .and. cp_i == 2 .and. am%recording) tmp = field

     if(.not.skipcp) CALL POPINTEGER4(field)

     if (am%check_st_integer .and. cp_i == 2 .and. am%recording) then
        if ((field - tmp == 0)) then
           am%st_integer(cp_t,am%count_pop(cp_t,idx_integer)) = 0
           am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - 1
        elseif (field == 0) then
           am%st_integer(cp_t,am%count_pop(cp_t,idx_integer)) = -1
           am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - 1
        endif
        am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) - 1
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_integer) then
        docp = am%st_integer(cp_t,am%count_pop(cp_t,idx_integer))
     endif
     am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) - 1

     if (docp == 1) then
        field = am%cp_integer(cp_t,am%index_pop(cp_t,idx_integer))
        am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - 1
     elseif (docp == -1) then
        field = 0
     endif

  endif

 end subroutine pop_integer_k0

 subroutine pop_integer_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  integer, intent(inout) :: field(dimen)
  logical, optional, intent(in) :: skip
  integer :: tmp(dimen)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPINTEGER4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) + 1

     if (am%check_st_integer .and. cp_i == 2 .and. am%recording) tmp = field

     if(.not.skipcp) CALL POPINTEGER4ARRAY(field,dimen)

     if (am%check_st_integer .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(field - tmp)) == 0) then
           am%st_integer(cp_t,am%count_pop(cp_t,idx_integer)) = 0
           am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - dimen
        elseif (maxval(abs(field)) == 0) then
           am%st_integer(cp_t,am%count_pop(cp_t,idx_integer)) = -1
           am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - dimen
        endif
        am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) - 1
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_integer) then
        docp = am%st_integer(cp_t,am%count_pop(cp_t,idx_integer))
     endif
     am%count_pop(cp_t,idx_integer) = am%count_pop(cp_t,idx_integer) - 1

     if (docp == 1) then
        field = am%cp_integer(cp_t,am%index_pop(cp_t,idx_integer)-dimen+1:am%index_pop(cp_t,idx_integer))
        am%index_pop(cp_t,idx_integer) = am%index_pop(cp_t,idx_integer) - dimen
     elseif (docp == -1) then
        field = 0
     endif

  endif

 end subroutine pop_integer_k1

! pushrealarray - r4
! ------------------

 subroutine psh_real_r4_k0(field,skip)

  implicit none

  real(4), intent(in) :: field
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL4(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + 1
     endif

     if(.not.skipcp) CALL PUSHREAL4(field)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_psh(cp_t,idx_real_r4))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + 1
        am%cp_real_r4(cp_t,am%index_psh(cp_t,idx_real_r4)) = field
     endif

  endif

 end subroutine psh_real_r4_k0

 subroutine psh_real_r4_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(in) :: field(dimen)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_psh(cp_t,idx_real_r4))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
        am%cp_real_r4(cp_t,am%index_psh(cp_t,idx_real_r4)-dimen+1:am%index_psh(cp_t,idx_real_r4)) = field
     endif

  endif

 end subroutine psh_real_r4_k1

 subroutine psh_real_r4_k2(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(in) :: field(:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_psh(cp_t,idx_real_r4))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
        am%cp_real_r4(cp_t,am%index_psh(cp_t,idx_real_r4)-dimen+1:am%index_psh(cp_t,idx_real_r4)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r4_k2

 subroutine psh_real_r4_k3(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(in) :: field(:,:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_psh(cp_t,idx_real_r4))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
        am%cp_real_r4(cp_t,am%index_psh(cp_t,idx_real_r4)-dimen+1:am%index_psh(cp_t,idx_real_r4)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r4_k3

 subroutine psh_real_r4_k4(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(in) :: field(:,:,:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL4ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r4) = am%count_psh(cp_t,idx_real_r4) + 1

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_psh(cp_t,idx_real_r4))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r4) = am%index_psh(cp_t,idx_real_r4) + dimen
        am%cp_real_r4(cp_t,am%index_psh(cp_t,idx_real_r4)-dimen+1:am%index_psh(cp_t,idx_real_r4)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r4_k4

! poprealarray - r4
! -----------------

 subroutine pop_real_r4_k0(field,skip)

  implicit none

  real(4), intent(inout) :: field
  logical, optional, intent(in) :: skip
  real(4) :: tmp
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL4(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) + 1

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) tmp = field

     if(.not.skipcp) CALL POPREAL4(field)

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        if ((field - tmp == 0.0_4)) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = 0
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - 1
        elseif (field == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = -1
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - 1
        endif
        am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4))
     endif
     am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1

     if (docp == 1) then
        if(.not.skipcp) field = am%cp_real_r4(cp_t,am%index_pop(cp_t,idx_real_r4))
        am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - 1
     elseif (docp == -1) then
        if(.not.skipcp) field = 0.0_4
     endif

  endif

 end subroutine pop_real_r4_k0

 subroutine pop_real_r4_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(inout) :: field(dimen)
  logical, optional, intent(in) :: skip
  real(4), allocatable :: tmp(:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) + 1

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen))
        tmp = field
     endif

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(field - tmp)) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = 0
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        elseif (maxval(abs(field)) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = -1
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        endif
        am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4))
     endif
     am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1

     if (docp == 1) then
        field = am%cp_real_r4(cp_t,am%index_pop(cp_t,idx_real_r4)-dimen+1:am%index_pop(cp_t,idx_real_r4))
        am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
     elseif (docp == -1) then
        field = 0.0_4
     endif

  endif

 end subroutine pop_real_r4_k1

 subroutine pop_real_r4_k2(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(inout) :: field(:,:)
  logical, optional, intent(in) :: skip
  real(4), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) + 1

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = 0
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        elseif (maxval(abs(field)) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = -1
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        endif
        am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4))
     endif
     am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r4(cp_t,am%index_pop(cp_t,idx_real_r4)-dimen+1:am%index_pop(cp_t,idx_real_r4)),&
                       (/size(field,1),size(field,2)/))
        am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
     elseif (docp == -1) then
        field = 0.0_4
     endif

  endif

 end subroutine pop_real_r4_k2

 subroutine pop_real_r4_k3(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(inout) :: field(:,:,:)
  logical, optional, intent(in) :: skip
  real(4), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) + 1

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = 0
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        elseif (maxval(abs(field)) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = -1
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        endif
        am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4))
     endif
     am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r4(cp_t,am%index_pop(cp_t,idx_real_r4)-dimen+1:am%index_pop(cp_t,idx_real_r4)),&
                       (/size(field,1),size(field,2),size(field,3)/))
        am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
     elseif (docp == -1) then
        field = 0.0_4
     endif

  endif

 end subroutine pop_real_r4_k3

 subroutine pop_real_r4_k4(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(4), intent(inout) :: field(:,:,:,:)
  logical, optional, intent(in) :: skip
  real(4), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) + 1

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL4ARRAY(field,dimen)

     if (am%check_st_real_r4 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = 0
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        elseif (maxval(abs(field)) == 0.0_4) then
           am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4)) = -1
           am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
        endif
        am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r4) then
        docp = am%st_real_r4(cp_t,am%count_pop(cp_t,idx_real_r4))
     endif
     am%count_pop(cp_t,idx_real_r4) = am%count_pop(cp_t,idx_real_r4) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r4(cp_t,am%index_pop(cp_t,idx_real_r4)-dimen+1:am%index_pop(cp_t,idx_real_r4)),&
                       (/size(field,1),size(field,2),size(field,3),size(field,4)/))
        am%index_pop(cp_t,idx_real_r4) = am%index_pop(cp_t,idx_real_r4) - dimen
     elseif (docp == -1) then
        field = 0.0_4
     endif

  endif

 end subroutine pop_real_r4_k4

! pushrealarray - r8
! ------------------

 subroutine psh_real_r8_k0(field,skip)

  implicit none

  real(8), intent(in) :: field
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL8(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + 1
     endif

     if(.not.skipcp) CALL PUSHREAL8(field)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_psh(cp_t,idx_real_r8))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + 1
        am%cp_real_r8(cp_t,am%index_psh(cp_t,idx_real_r8)) = field
     endif

  endif

 end subroutine psh_real_r8_k0

 subroutine psh_real_r8_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(in) :: field(dimen)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_psh(cp_t,idx_real_r8))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
        am%cp_real_r8(cp_t,am%index_psh(cp_t,idx_real_r8)-dimen+1:am%index_psh(cp_t,idx_real_r8)) = field
     endif

  endif

 end subroutine psh_real_r8_k1

 subroutine psh_real_r8_k2(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(in) :: field(:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_psh(cp_t,idx_real_r8))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
        am%cp_real_r8(cp_t,am%index_psh(cp_t,idx_real_r8)-dimen+1:am%index_psh(cp_t,idx_real_r8)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r8_k2

 subroutine psh_real_r8_k3(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(in) :: field(:,:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_psh(cp_t,idx_real_r8))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
        am%cp_real_r8(cp_t,am%index_psh(cp_t,idx_real_r8)-dimen+1:am%index_psh(cp_t,idx_real_r8)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r8_k3

 subroutine psh_real_r8_k4(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(in) :: field(:,:,:,:)
  logical, optional, intent(in) :: skip
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (am%recording) then
        am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
     endif

     if(.not.skipcp) CALL PUSHREAL8ARRAY(field,dimen)

  elseif (cp_i == 3 .and. am%recording) then

     am%count_psh(cp_t,idx_real_r8) = am%count_psh(cp_t,idx_real_r8) + 1

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_psh(cp_t,idx_real_r8))
     endif

     if (docp == 1) then
        am%index_psh(cp_t,idx_real_r8) = am%index_psh(cp_t,idx_real_r8) + dimen
        am%cp_real_r8(cp_t,am%index_psh(cp_t,idx_real_r8)-dimen+1:am%index_psh(cp_t,idx_real_r8)) = reshape(field,(/dimen/))
     endif

  endif

 end subroutine psh_real_r8_k4

! poprealarray - r8
! -----------------

 subroutine pop_real_r8_k0(field,skip)

  implicit none

  real(8), intent(inout) :: field
  logical, optional, intent(in) :: skip
  real(8) :: tmp
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL8(field)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) + 1

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) tmp = field

     if(.not.skipcp) CALL POPREAL8(field)

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        if ((field - tmp == 0.0_8)) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = 0
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - 1
        elseif (field == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = -1
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - 1
        endif
        am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8))
     endif
     am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1

     if (docp == 1) then
        if(.not.skipcp) field = am%cp_real_r8(cp_t,am%index_pop(cp_t,idx_real_r8))
        am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - 1
     elseif (docp == -1) then
        if(.not.skipcp) field = 0.0_8
     endif

  endif

 end subroutine pop_real_r8_k0

 subroutine pop_real_r8_k1(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(inout) :: field(dimen)
  logical, optional, intent(in) :: skip
  real(8), allocatable :: tmp(:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) + 1

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen))
        tmp = field
     endif

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(field - tmp)) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = 0
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        elseif (maxval(abs(field)) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = -1
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        endif
        am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8))
     endif
     am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1

     if (docp == 1) then
        field = am%cp_real_r8(cp_t,am%index_pop(cp_t,idx_real_r8)-dimen+1:am%index_pop(cp_t,idx_real_r8))
        am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
     elseif (docp == -1) then
        field = 0.0_8
     endif

  endif

 end subroutine pop_real_r8_k1

 subroutine pop_real_r8_k2(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(inout) :: field(:,:)
  logical, optional, intent(in) :: skip
  real(8), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) + 1

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = 0
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        elseif (maxval(abs(field)) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = -1
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        endif
        am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8))
     endif
     am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r8(cp_t,am%index_pop(cp_t,idx_real_r8)-dimen+1:am%index_pop(cp_t,idx_real_r8)),&
                       (/size(field,1),size(field,2)/))
        am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
     elseif (docp == -1) then
        field = 0.0_8
     endif

  endif

 end subroutine pop_real_r8_k2

 subroutine pop_real_r8_k3(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(inout) :: field(:,:,:)
  logical, optional, intent(in) :: skip
  real(8), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) + 1

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = 0
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        elseif (maxval(abs(field)) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = -1
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        endif
        am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8))
     endif
     am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r8(cp_t,am%index_pop(cp_t,idx_real_r8)-dimen+1:am%index_pop(cp_t,idx_real_r8)),&
                       (/size(field,1),size(field,2),size(field,3)/))
        am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
     elseif (docp == -1) then
        field = 0.0_8
     endif

  endif

 end subroutine pop_real_r8_k3

 subroutine pop_real_r8_k4(field,dimen,skip)

  implicit none

  integer, intent(in) :: dimen
  real(8), intent(inout) :: field(:,:,:,:)
  logical, optional, intent(in) :: skip
  real(8), allocatable :: tmp(:,:)
  logical :: skipcp
  integer(status_kind) :: docp

  skipcp = .false.
  if (present(skip)) then
     skipcp = skip
  endif

  if (cp_i == 0) then

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

  elseif (cp_i == 1 .or. cp_i == 2) then

     if (cp_i == 1) am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) + 1

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        allocate(tmp(dimen,1))
        tmp = reshape(field,(/dimen, 1/))
     endif

     if(.not.skipcp) CALL POPREAL8ARRAY(field,dimen)

     if (am%check_st_real_r8 .and. cp_i == 2 .and. am%recording) then
        if (maxval(abs(tmp - reshape(field,(/dimen, 1/)))) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = 0
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        elseif (maxval(abs(field)) == 0.0_8) then
           am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8)) = -1
           am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
        endif
        am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1
        deallocate(tmp)
     endif

  elseif (cp_i >= 3) then

     docp = 1
     if (am%check_st_real_r8) then
        docp = am%st_real_r8(cp_t,am%count_pop(cp_t,idx_real_r8))
     endif
     am%count_pop(cp_t,idx_real_r8) = am%count_pop(cp_t,idx_real_r8) - 1

     if (docp == 1) then
        field = reshape(am%cp_real_r8(cp_t,am%index_pop(cp_t,idx_real_r8)-dimen+1:am%index_pop(cp_t,idx_real_r8)),&
                       (/size(field,1),size(field,2),size(field,3),size(field,4)/))
        am%index_pop(cp_t,idx_real_r8) = am%index_pop(cp_t,idx_real_r8) - dimen
     elseif (docp == -1) then
        field = 0.0_8
     endif

  endif

 end subroutine pop_real_r8_k4



! pushrealarray_adm
! -----------------

 subroutine psh_adm_real_r4_k0(field)

  implicit none
  real(4), intent(inout) :: field

  CALL PUSHREAL4(field)

 end subroutine psh_adm_real_r4_k0

 subroutine psh_adm_real_r4_k1(field,dimen)

  implicit none
  real(4), intent(inout) :: field(:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL4ARRAY(field(:),dimen)

 end subroutine psh_adm_real_r4_k1

 subroutine psh_adm_real_r4_k2(field,dimen)

  implicit none

  real(4), intent(inout) :: field(:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL4ARRAY(field(:,:),dimen)

 end subroutine psh_adm_real_r4_k2

 subroutine psh_adm_real_r4_k3(field,dimen)

  implicit none
  real(4), intent(inout) :: field(:,:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL4ARRAY(field(:,:,:),dimen)

 end subroutine psh_adm_real_r4_k3

 subroutine psh_adm_real_r4_k4(field,dimen)

  implicit none
  real(4), intent(inout) :: field(:,:,:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL4ARRAY(field(:,:,:,:),dimen)

 end subroutine psh_adm_real_r4_k4

 subroutine psh_adm_real_r8_k0(field)

  implicit none
  real(8), intent(inout) :: field

  CALL PUSHREAL8(field)

 end subroutine psh_adm_real_r8_k0

 subroutine psh_adm_real_r8_k1(field,dimen)

  implicit none
  real(8), intent(inout) :: field(:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL8ARRAY(field(:),dimen)

 end subroutine psh_adm_real_r8_k1

 subroutine psh_adm_real_r8_k2(field,dimen)

  implicit none
  real(8), intent(inout) :: field(:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL8ARRAY(field(:,:),dimen)

 end subroutine psh_adm_real_r8_k2

 subroutine psh_adm_real_r8_k3(field,dimen)

  implicit none
  real(8), intent(inout) :: field(:,:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL8ARRAY(field(:,:,:),dimen)

 end subroutine psh_adm_real_r8_k3

 subroutine psh_adm_real_r8_k4(field,dimen)

  implicit none
  real(8), intent(inout) :: field(:,:,:,:)
  integer, intent(in   ) :: dimen

  CALL PUSHREAL8ARRAY(field(:,:,:,:),dimen)

 end subroutine psh_adm_real_r8_k4



! poprealarray
! -------------

!These routines are included to provide an abstract layer for real 
!number precision when checkpointing the reference state. They allow
!for compile-time precision choice by an end user, rather than at 
!the time Tapenade is used to generate the adjont code.

 subroutine pop_adm_real_r4_k0(field)

  implicit none

  real(4), intent(inout) :: field

  CALL POPREAL4(field)

 end subroutine pop_adm_real_r4_k0

 subroutine pop_adm_real_r4_k1(field,dimen)

  implicit none

  real(4), intent(inout) :: field(:)
  integer, intent(in   ) :: dimen

  CALL POPREAL4ARRAY(field(:),dimen)

 end subroutine pop_adm_real_r4_k1


 subroutine pop_adm_real_r4_k2(field,dimen)

  implicit none

  real(4), intent(inout) :: field(:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL4ARRAY(field(:,:),dimen)

 end subroutine pop_adm_real_r4_k2


 subroutine pop_adm_real_r4_k3(field,dimen)

  implicit none

  real(4), intent(inout) :: field(:,:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL4ARRAY(field(:,:,:),dimen)

 end subroutine pop_adm_real_r4_k3


 subroutine pop_adm_real_r4_k4(field,dimen)

  implicit none

  real(4), intent(inout) :: field(:,:,:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL4ARRAY(field(:,:,:,:),dimen)

 end subroutine pop_adm_real_r4_k4

 subroutine pop_adm_real_r8_k0(field)

  implicit none

  real(8), intent(inout) :: field

  CALL POPREAL8(field)

 end subroutine pop_adm_real_r8_k0

 subroutine pop_adm_real_r8_k1(field,dimen)

  implicit none

  real(8), intent(inout) :: field(:)
  integer, intent(in   ) :: dimen

  CALL POPREAL8ARRAY(field(:),dimen)

 end subroutine pop_adm_real_r8_k1


 subroutine pop_adm_real_r8_k2(field,dimen)

  implicit none

  real(8), intent(inout) :: field(:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL8ARRAY(field(:,:),dimen)

 end subroutine pop_adm_real_r8_k2


 subroutine pop_adm_real_r8_k3(field,dimen)

  implicit none

  real(8), intent(inout) :: field(:,:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL8ARRAY(field(:,:,:),dimen)

 end subroutine pop_adm_real_r8_k3


 subroutine pop_adm_real_r8_k4(field,dimen)

  implicit none

  real(8), intent(inout) :: field(:,:,:,:)
  integer, intent(in   ) :: dimen

  CALL POPREAL8ARRAY(field(:,:,:,:),dimen)

 end subroutine pop_adm_real_r8_k4

endmodule tapenade_iter
