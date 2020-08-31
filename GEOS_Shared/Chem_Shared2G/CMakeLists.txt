esma_set_this ()

set (srcs
  Chem_AeroGeneric.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

#include_directories (${esma_include}/${this})

