 module num_parthd

 contains
 function num_parthds()
       integer :: num_parthds

!      use omp_lib
!!$omp parallel
!     num_parthds = omp_get_num_threads()
      num_parthds = 8
!     num_parthds = 6
!     num_parthds = 4
!!$omp end parallel
      return
 end function
 end module num_parthd


 
