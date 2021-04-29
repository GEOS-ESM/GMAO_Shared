esma_set_this()

set (srcs
   m_SHMEM.F90
   mpi_abort.F90
   mpi_comm_dup.F90
   mpi_comm_free.F90
   mpi_comm_rank.F90
   mpi_comm_size.F90
   mpi_error_string.F90
   mpi_finalize.F90
   mpi_init.F90
   mpi_initialized.F90
   mpi_isend.F90
   mpi_send.F90
   mpi_sendrecv.F90
   mpi_recv.F90
   mpi_irecv.F90
   mpi_waitall.F90
   mpi_wtick.F90
   mpi_wtime.F90
   mpi_barrier.F90
   mpi_bcast.F90
   mpi_gather.F90
   mpi_gatherv.F90
   mpi_scatter.F90
   mpi_scatterv.F90
   mpi_allgather.F90
   mpi_allgatherv.F90
   mpi_alltoall.F90
   mpi_alltoallv.F90
   mpi_reduce.F90
   mpi_allreduce.F90
   mpi_reduce_scatter.F90
   mpi0_copy.F90
   mpi_zzz.F90
   mpi_wait.F90
   mpi_comm_split.F90
   m_mpi0.F
   Const_MPItype.H
   )

esma_add_library(${this} SRCS ${srcs})
target_include_directories (${this} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../jic>)
