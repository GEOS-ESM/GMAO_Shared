
      character*256,  parameter :: c_slp  = 'slp'
      character*256,  parameter :: c_ps   = 'ps'
      character*256,  parameter :: c_dp   = 'dp'
      character*256,  parameter :: c_pl   = 'pl'
      character*256,  parameter :: c_ple  = 'ple'
      character*256,  parameter :: c_u    = 'u'
      character*256,  parameter :: c_v    = 'v'
      character*256,  parameter :: c_t    = 't'
      character*256,  parameter :: c_q    = 'q'
      character*256,  parameter :: c_h    = 'h'
      character*256,  parameter :: c_gz   = 'gz'
      character*256,  parameter :: c_th   = 'th'
      character*256,  parameter :: c_tv   = 'tv'
      character*256,  parameter :: c_thv  = 'thv'
      character*256,  parameter :: c_phis = 'phis'

